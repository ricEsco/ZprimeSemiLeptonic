#!/usr/bin/env python
"""
ttbargen_statusflags_miniaod.py

UHH2 GenParticles do not store reco::GenStatusFlags, so this script works on MiniAOD:

  1. reads prunedGenParticles (FWLite),
  2. emulates UHH2's prunedPrunedGenParticles (same GenParticlePruner rules and the same
     mother/daughter re-linking), i.e. reproduces the GenParticles content of the ntuple,
  3. runs the TTbarGen selection on it (line-by-line port of common/src/TTbarGen.cxx),
  4. reports, for every TTbarGen object, the Pythia status and the GenStatusFlags of the
     corresponding MiniAOD particle, and compares it with the *last copy* of that particle
     (followed through prunedGenParticles), i.e. what an "isLastCopy" selection would pick.

Run inside CMSSW (cmsenv, CMSSW_10_6_X), with a grid proxy for xrootd files:
  python ttbargen_statusflags_miniaod.py root://cms-xrd-global.cern.ch//store/mc/.../file.root \
      [--max-events 5000] [--print-events 5] [--uhh2-file analysis_output.root]

--uhh2-file restricts the MiniAOD events to the (run, lumi, event) present in an UHH2
AnalysisTree (optional; needs those events to be in the MiniAOD files you pass).
"""
from __future__ import print_function, division

import argparse
import math
import sys
from collections import Counter, defaultdict, namedtuple

FLAGS = ["isPrompt", "isDecayedLeptonHadron", "isTauDecayProduct", "isPromptTauDecayProduct",
         "isDirectTauDecayProduct", "isDirectPromptTauDecayProduct", "isDirectHadronDecayProduct",
         "isHardProcess", "fromHardProcess", "isHardProcessTauDecayProduct",
         "isDirectHardProcessTauDecayProduct", "fromHardProcessBeforeFSR",
         "isFirstCopy", "isLastCopy", "isLastCopyBeforeFSR"]

LABELS = ["Top", "Antitop", "WTop", "WAntitop", "bTop", "bAntitop",
          "Wdecay1", "Wdecay2", "WMinusdecay1", "WMinusdecay2"]

# Pruned-collection record used by the TTbarGen port (same fields as UHH2 GenParticle)
GP = namedtuple("GP", "index pdgId status mother1 mother2 daughter1 daughter2 orig")


# ---------------------------------------------------------------------------
# 1) read prunedGenParticles into plain python lists
# ---------------------------------------------------------------------------
def read_event(gps):
    n = gps.size()
    pdg, st, moth, daus = [0] * n, [0] * n, [None] * n, [None] * n
    for i in range(n):
        p = gps[i]
        pdg[i] = p.pdgId()
        st[i] = p.status()
        moth[i] = [p.motherRef(j).key() for j in range(p.numberOfMothers())]
        daus[i] = [p.daughterRef(j).key() for j in range(p.numberOfDaughters())]
    return pdg, st, moth, daus


# ---------------------------------------------------------------------------
# 2) emulate UHH2 prunedPrunedGenParticles (core/python/ntuple_generator.py)
# ---------------------------------------------------------------------------
def uhh2_keep(i, pdg, st, moth, daus):
    a = abs(pdg[i])
    keep = False                                       # 'drop *'
    if st[i] == 3:                                     # 'keep status == 3'
        keep = True
    if 20 <= st[i] <= 30:                              # 'keep 20 <= status <= 30'
        keep = True
    is_lep = 11 <= a <= 16
    one_mother = len(moth[i]) == 1
    mpdg = abs(pdg[moth[i][0]]) if one_mother else 0
    if is_lep and one_mother and 23 <= mpdg <= 25:
        keep = True
    if is_lep and one_mother and mpdg == 6:
        keep = True
    if is_lep and one_mother and mpdg == 42:
        keep = True
    if is_lep and one_mother and mpdg == 6:            # 'drop ... mother == top'
        keep = False
    if is_lep and one_mother and mpdg == 6:            # '... unless the top has no W daughter'
        md = daus[moth[i][0]]
        if len(md) > 2 and all(abs(pdg[md[k]]) != 24 for k in range(3)):
            keep = True
    return keep


def prune(pdg, st, moth, daus):
    """GenParticlePruner: keep selected particles, re-link mothers/daughters to the nearest
    kept ancestor/descendant (recursively through dropped particles)."""
    n = len(pdg)
    kept = [uhh2_keep(i, pdg, st, moth, daus) for i in range(n)]
    new_idx = {}
    for i in range(n):
        if kept[i]:
            new_idx[i] = len(new_idx)

    def kept_daughters(i, out, seen):
        for d in daus[i]:
            if d in seen:
                continue
            seen.add(d)
            if kept[d]:
                out.append(new_idx[d])
            else:
                kept_daughters(d, out, seen)
        return out

    def kept_mothers(i, out, seen):
        for m in moth[i]:
            if m in seen:
                continue
            seen.add(m)
            if kept[m]:
                out.append(new_idx[m])
            else:
                kept_mothers(m, out, seen)
        return out

    pruned = []
    for i in range(n):
        if not kept[i]:
            continue
        ms = kept_mothers(i, [], set())
        ds = kept_daughters(i, [], set())
        pruned.append(GP(new_idx[i], pdg[i], st[i],
                         ms[0] if len(ms) > 0 else -1, ms[1] if len(ms) > 1 else -1,
                         ds[0] if len(ds) > 0 else -1, ds[1] if len(ds) > 1 else -1, i))
    return pruned


# ---------------------------------------------------------------------------
# 3) TTbarGen port (common/src/TTbarGen.cxx, as in the UHH2 version you use)
# ---------------------------------------------------------------------------
class TTbarGenFailure(Exception):
    pass


def ttbargen(gps):
    by_index = {g.index: g for g in gps}
    dau = lambda g, k: by_index.get(g.daughter1 if k == 1 else g.daughter2)
    mom = lambda g, k: by_index.get(g.mother1 if k == 1 else g.mother2)

    def has_top_mother(g, top):
        m1, m2 = mom(g, 1), mom(g, 2)
        return (m1 is not None and m1.index == top.index) or (m2 is not None and m2.index == top.index)

    res, n_top, n_antitop, top, antitop = {}, 0, 0, None, None
    for genp in gps:
        if abs(genp.pdgId) != 6:
            continue
        if top is None and genp.pdgId == 6:
            top = genp
        if antitop is None and genp.pdgId == -6:
            antitop = genp
        w, b = dau(genp, 1), dau(genp, 2)
        if w is None or b is None:
            raise TTbarGenFailure("top has not ==2 daughters")
        if abs(w.pdgId) != 24:
            w, b = b, w
        if abs(w.pdgId) != 24:
            for g in gps:
                if has_top_mother(g, genp) and abs(g.pdgId) == 24:
                    w = g
                    break
        if abs(b.pdgId) not in (5, 3, 1):
            for g in gps:
                if has_top_mother(g, genp) and abs(g.pdgId) in (5, 3, 1):
                    b = g
                    break
        if abs(b.pdgId) not in (5, 3, 1):
            raise TTbarGenFailure("top has no b daughter")
        wd1, wd2 = dau(w, 1), dau(w, 2)
        while True:
            if wd1 is not None and wd2 is None:
                w = wd1
                wd1, wd2 = dau(w, 1), dau(w, 2)
            elif wd1 is not None and wd2 is not None:
                break
            else:
                raise TTbarGenFailure("W has no daughters")
        if genp.pdgId == 6:
            res.update(Top=top, WTop=w, bTop=b, Wdecay1=wd1, Wdecay2=wd2)
            n_top += 1
        else:
            res.update(Antitop=antitop, WAntitop=w, bAntitop=b, WMinusdecay1=wd1, WMinusdecay2=wd2)
            n_antitop += 1
    if n_top != 1 or n_antitop != 1:
        raise TTbarGenFailure("did not find exactly one top and one antitop")
    return res


# ---------------------------------------------------------------------------
# 4) helpers: last copy, kinematics
# ---------------------------------------------------------------------------
def last_copy(i, pdg, daus):
    """Follow daughters with the same pdgId until there are none (last copy of the chain)."""
    seen = set()
    while True:
        nxt = [d for d in daus[i] if pdg[d] == pdg[i]]
        if not nxt or nxt[0] in seen:
            return i
        seen.add(i)
        i = nxt[0]


def dR(a, b):
    dphi = abs(a.phi() - b.phi())
    if dphi > math.pi:
        dphi = 2 * math.pi - dphi
    return math.hypot(a.eta() - b.eta(), dphi)


class Stat(object):
    def __init__(self):
        self.n = 0; self.s = 0.0; self.s2 = 0.0
    def add(self, x):
        self.n += 1; self.s += x; self.s2 += x * x
    def mean(self):
        return self.s / self.n if self.n else float("nan")
    def rms(self):
        return math.sqrt(max(self.s2 / self.n - self.mean() ** 2, 0.0)) if self.n else float("nan")


# ---------------------------------------------------------------------------
# Four-momentum consistency checks (pure python, (E, px, py, pz) tuples)
# ---------------------------------------------------------------------------
TOL_M, TOL_RELPT, TOL_DR = 0.05, 1e-3, 1e-3      # "identical" = all three within tolerance


def p4_of(p):
    return (p.energy(), p.px(), p.py(), p.pz())


def p4_add(*vs):
    return tuple(sum(c) for c in zip(*vs))


def p4_m(v):
    m2 = v[0] ** 2 - v[1] ** 2 - v[2] ** 2 - v[3] ** 2
    return math.copysign(math.sqrt(abs(m2)), m2)


def p4_pt(v):
    return math.hypot(v[1], v[2])


def p4_dr(a, b):
    def etaphi(v):
        pt = p4_pt(v)
        eta = math.asinh(v[3] / pt) if pt > 0 else math.copysign(1e9, v[3])
        return eta, math.atan2(v[2], v[1])
    ea, pa = etaphi(a)
    eb, pb = etaphi(b)
    dphi = abs(pa - pb)
    if dphi > math.pi:
        dphi = 2 * math.pi - dphi
    return math.hypot(ea - eb, dphi)


class Compare(object):
    """Compares a reference four-vector a with b: mass, relative pT, dR."""
    def __init__(self):
        self.n = self.same_m = self.identical = 0
        self.dm, self.dpt, self.dr = Stat(), Stat(), Stat()

    def add(self, a, b):
        dm = p4_m(b) - p4_m(a)
        dpt = (p4_pt(b) - p4_pt(a)) / p4_pt(a) if p4_pt(a) > 0 else 0.0
        dr = p4_dr(a, b)
        self.n += 1
        self.dm.add(dm); self.dpt.add(dpt); self.dr.add(dr)
        if abs(dm) < TOL_M:
            self.same_m += 1
            if abs(dpt) < TOL_RELPT and dr < TOL_DR:
                self.identical += 1
        return dm, dpt, dr

    def line(self, name):
        if not self.n:
            return "  %-40s (no entries)" % name
        return "  %-40s %7d %8.1f%% %9.1f%%  %8.3f+-%-7.3f %8.4f+-%-7.4f %7.4f+-%-7.4f" % (
            name, self.n, 100.0 * self.same_m / self.n, 100.0 * self.identical / self.n,
            self.dm.mean(), self.dm.rms(), self.dpt.mean(), self.dpt.rms(), self.dr.mean(), self.dr.rms())


def print_compare_header():
    print("  %-40s %7s %9s %10s  %16s %16s %16s" % ("comparison (reference vs other)", "N", "same m", "identical",
                                                   "dm [GeV]", "dpT/pT", "dR"))
    print("  %-40s %7s %9s %10s  %16s %16s %16s" % ("", "", "", "", "mean +- rms", "mean +- rms", "mean +- rms"))


def object_role(label, part):
    """Group leptonic/hadronic W daughters for the summary."""
    if label.startswith("W") and "decay" in label:
        a = abs(part.pdgId)
        if a in (11, 13):
            return "W daughter: e/mu"
        if a == 15:
            return "W daughter: tau"
        if a in (12, 14, 16):
            return "W daughter: neutrino"
        return "W daughter: quark"
    return {"Top": "top", "Antitop": "top", "WTop": "W", "WAntitop": "W",
            "bTop": "b", "bAntitop": "b"}[label]


ROLES = ["top", "W", "b", "W daughter: quark", "W daughter: e/mu", "W daughter: tau", "W daughter: neutrino"]


# ---------------------------------------------------------------------------
# Event dump (for building a Feynman-style diagram of one event)
# ---------------------------------------------------------------------------
_NAMES = {1: "d", 2: "u", 3: "s", 4: "c", 5: "b", 6: "t", 21: "g", 22: "gamma", 23: "Z", 24: "W+", 25: "H",
          11: "e-", 12: "nu_e", 13: "mu-", 14: "nu_mu", 15: "tau-", 16: "nu_tau", 2212: "p", 2112: "n",
          211: "pi+", 111: "pi0", 321: "K+", 130: "K0L", 310: "K0S", 311: "K0"}


def pname(pid):
    a = abs(pid)
    if a in (24,):
        return "W+" if pid > 0 else "W-"
    if a in (11, 13, 15):
        n = _NAMES[a]
        return n if pid > 0 else n[:-1] + "+"
    if a in (211, 321):
        n = _NAMES[a]
        return n if pid > 0 else n[:-1] + "-"
    n = _NAMES.get(a, str(a))
    return n + ("~" if pid < 0 and a not in (21, 22, 23, 25, 111, 130, 310) else "")


def short_flags(p):
    sf = p.statusFlags()
    out = []
    for f, tag in (("isHardProcess", "HP"), ("fromHardProcess", "fromHP"), ("isFirstCopy", "first"),
                   ("isLastCopy", "last"), ("isPrompt", "prompt")):
        if getattr(sf, f)():
            out.append(tag)
    return ",".join(out)


def copy_chain(i, pdg, daus):
    """Same-pdgId copy chain starting at i (inclusive) down to the last copy."""
    chain, seen = [i], set([i])
    while True:
        nxt = [d for d in daus[chain[-1]] if pdg[d] == pdg[chain[-1]] and d not in seen]
        if not nxt:
            return chain
        chain.append(nxt[0]); seen.add(nxt[0])


def p4_rap(v):
    return 0.5 * math.log((v[0] + v[3]) / (v[0] - v[3])) if v[0] > abs(v[3]) else float("nan")


def dump_event(aux, gps, pdg, st, moth, daus, sel, max_tree_lines=400):
    """Print a full, diagram-oriented dump of one event and return it as a JSON-serialisable dict."""
    n = len(pdg)
    roles = defaultdict(list)
    for lab in LABELS:
        roles[sel[lab].orig].append(lab + " (TTbarGen)")
    chains = {}
    for tl in ("Top", "Antitop"):
        ch = copy_chain(sel[tl].orig, pdg, daus)
        chains[tl] = ch
        for k, i in enumerate(ch[1:], 1):
            roles[i].append("%s copy %d" % (tl, k))
        roles[ch[-1]].append("%s LAST COPY" % tl)
        roles[ch[0]].append("%s FIRST COPY" % tl)

    hdr = "run %d lumi %d event %d" % (aux.run(), aux.luminosityBlock(), aux.event())
    print("\n" + "#" * 100)
    print("# EVENT DUMP  %s   (prunedGenParticles, %d particles)" % (hdr, n))
    print("#" * 100)

    # (1) full table
    print("\n(1) Full particle table  [mothers/daughters are indices in this table; flags: HP=isHardProcess, "
          "fromHP, first/last=isFirst/LastCopy]")
    print("  %4s %8s %-7s %4s %-14s %-22s %8s %7s %7s %8s  %-28s %s" % (
        "idx", "pdgId", "name", "st", "mothers", "daughters", "pt", "eta", "phi", "mass", "flags", "role"))
    parts = []
    for i in range(n):
        p = gps[i]
        fl = short_flags(p)
        ms = ",".join(str(m) for m in moth[i]) or "-"
        ds = daus[i]
        dstr = ",".join(str(d) for d in ds[:6]) + (",..(%d)" % len(ds) if len(ds) > 6 else "") if ds else "-"
        print("  %4d %8d %-7s %4d %-14s %-22s %8.2f %7.2f %7.2f %8.2f  %-28s %s" % (
            i, pdg[i], pname(pdg[i]), st[i], ms, dstr, p.pt(), p.eta(), p.phi(), p.mass(), fl, "; ".join(roles.get(i, []))))
        parts.append({"idx": i, "pdgId": pdg[i], "name": pname(pdg[i]), "status": st[i],
                      "mothers": list(moth[i]), "daughters": list(ds),
                      "pt": p.pt(), "eta": p.eta(), "phi": p.phi(), "mass": p.mass(),
                      "px": p.px(), "py": p.py(), "pz": p.pz(), "E": p.energy(),
                      "flags": fl, "roles": roles.get(i, [])})

    # (2) decay/history tree from the beam protons (hadrons summarised, partons at status 71/72 not expanded)
    print("\n(2) History tree (from beam particles; repeated nodes shown as '-> see idx'; hadrons/photons summarised)")
    roots = [i for i in range(n) if not moth[i]]
    lines, seen = [], set()

    def is_hadronic_leaf(i):
        return abs(pdg[i]) >= 100 and pdg[i] != 2212 or (pdg[i] == 22 and st[i] == 1)

    def walk(i, depth):
        if len(lines) >= max_tree_lines:
            return
        p = gps[i]
        tag = ("   <== " + "; ".join(roles[i])) if i in roles else ""
        if i in seen:
            lines.append("  " + "   " * depth + "%s [%d] -> see idx %d" % (pname(pdg[i]), i, i))
            return
        seen.add(i)
        lines.append("  " + "   " * depth + "%s [idx %d, st %d] pt=%.1f eta=%.2f phi=%.2f m=%.1f%s" % (
            pname(pdg[i]), i, st[i], p.pt(), p.eta(), p.phi(), p.mass(), tag))
        if st[i] in (71, 72) or (is_hadronic_leaf(i) and i not in roles):
            if daus[i]:
                lines.append("  " + "   " * (depth + 1) + "... %d further descendants (hadronisation/decays) not expanded" % len(daus[i]))
            return
        leaves = [d for d in daus[i] if is_hadronic_leaf(d) and d not in roles and not daus[d]]
        for d in daus[i]:
            if d in leaves:
                continue
            walk(d, depth + 1)
        if leaves:
            cnt = Counter(pname(pdg[d]) for d in leaves)
            lines.append("  " + "   " * (depth + 1) + "+ %d final-state hadrons/photons: %s" % (
                len(leaves), ", ".join("%dx%s" % (c, k) for k, c in cnt.most_common(6))))
    for r in roots:
        walk(r, 0)
    if len(lines) >= max_tree_lines:
        lines.append("  ... (tree truncated at %d lines)" % max_tree_lines)
    print("\n".join(lines))

    # (3) top-quark histories: every copy, momentum change per step, extra (non-top) daughters = radiation
    print("\n(3) Top-quark copy chains (first copy = TTbarGen Top(), last copy = isLastCopy = W + b)")
    chain_info = {}
    for tl in ("Top", "Antitop"):
        ch = chains[tl]
        print("  %s:" % tl)
        prev = None
        steps = []
        for k, i in enumerate(ch):
            v = p4_of(gps[i])
            dv = tuple(a - b for a, b in zip(v, prev)) if prev else None
            extra = [d for d in daus[i] if pdg[d] != pdg[i]]
            extra_s = ", ".join("%s[%d,st%d,pt%.1f]" % (pname(pdg[d]), d, st[d], gps[d].pt()) for d in extra)
            print("    copy %d: idx %3d st %2d  pt=%7.2f eta=%6.3f phi=%6.3f m=%7.2f E=%8.2f%s%s" % (
                k, i, st[i], gps[i].pt(), gps[i].eta(), gps[i].phi(), gps[i].mass(), gps[i].energy(),
                ("   dp=(%+.2f,%+.2f,%+.2f) dE=%+.2f" % (dv[1], dv[2], dv[3], dv[0])) if dv else "",
                ("   other daughters: " + extra_s) if extra_s else ""))
            steps.append({"idx": i, "status": st[i], "extra_daughters": extra})
            prev = v
        mothers_first = moth[ch[0]]
        print("    first copy mothers: %s" % ", ".join("%s[%d,st%d]" % (pname(pdg[m]), m, st[m]) for m in mothers_first))
        chain_info[tl] = steps

    # (4) ttbar system before/after showering, and the hard-process recoil partner(s)
    fc = p4_add(p4_of(gps[chains["Top"][0]]), p4_of(gps[chains["Antitop"][0]]))
    lc = p4_add(p4_of(gps[chains["Top"][-1]]), p4_of(gps[chains["Antitop"][-1]]))
    print("\n(4) ttbar system:  FIRST copies: m=%.2f pT=%.2f y=%.3f   |   LAST copies: m=%.2f pT=%.2f y=%.3f   "
          "|   dm=%+.2f dpT=%+.2f" % (p4_m(fc), p4_pt(fc), p4_rap(fc), p4_m(lc), p4_pt(lc), p4_rap(lc),
                                     p4_m(lc) - p4_m(fc), p4_pt(lc) - p4_pt(fc)))
    hp_other = [i for i in range(n) if st[i] == 23 and moth[i] and all(st[m] == 21 for m in moth[i])]
    for i in hp_other:
        print("    other hard-process outgoing parton (POWHEG emission): %s[%d] pt=%.2f eta=%.2f phi=%.2f" % (
            pname(pdg[i]), i, gps[i].pt(), gps[i].eta(), gps[i].phi()))
    print("    (at LHE level pT(ttbar, first copies) is balanced by the POWHEG emission; the change to the last copies")
    print("     comes from the parton shower: ISR recoil and radiation off the tops / colour partners)")

    return {"run": aux.run(), "lumi": aux.luminosityBlock(), "event": aux.event(), "particles": parts,
            "top_chain": chains["Top"], "antitop_chain": chains["Antitop"],
            "ttbar_first": {"m": p4_m(fc), "pt": p4_pt(fc), "y": p4_rap(fc)},
            "ttbar_last": {"m": p4_m(lc), "pt": p4_pt(lc), "y": p4_rap(lc)},
            "hard_process_extra": hp_other,
            "ttbargen": dict((lab, sel[lab].orig) for lab in LABELS)}


def read_uhh2_ids(fname, tree):
    import ROOT
    for lib in ("libTreePlayer", "libProof", "libProofPlayer"):
        ROOT.gSystem.Load(lib)
    f = ROOT.TFile.Open(fname)
    t = f.Get(tree)
    t.SetBranchStatus("*", 0)
    for b in ("run", "luminosityBlock", "event"):
        t.SetBranchStatus(b, 1)
    ids = set()
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        ids.add((int(t.run), int(t.luminosityBlock), int(t.event)))
    f.Close()
    return ids


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="MiniAOD files (local or root://...)")
    ap.add_argument("--max-events", type=int, default=5000, help="-1 = all")
    ap.add_argument("--print-events", type=int, default=5)
    ap.add_argument("--progress", type=int, default=500)
    ap.add_argument("--uhh2-file", default=None, help="restrict to events present in this UHH2 AnalysisTree")
    ap.add_argument("--uhh2-tree", default="AnalysisTree")
    ap.add_argument("--dump-events", type=int, default=0,
                    help="write a full diagram-oriented dump (table, history tree, top copy chains) for N events")
    ap.add_argument("--dump-ids", default="",
                    help='dump these events: comma-separated run:lumi:event, e.g. "1:142148:142147009"')
    ap.add_argument("--dump-min-dr", type=float, default=0.0,
                    help="only dump events where dR(first, last copy) of the top or antitop is at least this (e.g. 0.3)")
    ap.add_argument("--dump-max-particles", type=int, default=400,
                    help="skip events with more prunedGenParticles than this when choosing dump events (readability)")
    ap.add_argument("--dump-json", default="event_dump.json", help="JSON file with the dumped events")
    ap.add_argument("--dump-only", action="store_true",
                    help="stop as soon as all requested dumps are written (skip the statistics summary)")
    args = ap.parse_args()

    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.FWLiteEnabler.enable()
    from DataFormats.FWLite import Events, Handle

    wanted = read_uhh2_ids(args.uhh2_file, args.uhh2_tree) if args.uhh2_file else None
    if wanted is not None:
        print("Restricting to %d (run, lumi, event) from %s" % (len(wanted), args.uhh2_file))

    handle = Handle("std::vector<reco::GenParticle>")
    events = Events(args.files)

    n_read = n_proc = n_ok = 0
    fails = Counter()
    status = defaultdict(Counter)                 # role -> status -> n
    flags = defaultdict(Counter)                  # role -> flag -> n(True)
    n_role = Counter()
    is_last = Counter()                           # role -> TTbarGen particle == end of copy chain
    end_flag_ok = Counter()                       # role -> chain end has isLastCopy flag
    dpt = defaultdict(Stat); ddr = defaultdict(Stat); dm = defaultdict(Stat)
    tt = {"mtt": Stat(), "pttt": Stat()}
    # four-momentum consistency checks
    CMP_KEYS = ["t_last vs sum(daughters of t_last)",
                "t_TTbarGen vs W+b (TTbarGen)",
                "t_last vs W+b (TTbarGen)",
                "W_TTbarGen vs dau1+dau2 [W->lnu]", "W_last vs dau1+dau2 [W->lnu]",
                "W_TTbarGen vs dau1+dau2 [W->qq']", "W_last vs dau1+dau2 [W->qq']"]
    cmp_ = dict((k, Compare()) for k in CMP_KEYS)
    tlast_daus = Counter()        # (|pdgId|, status) signature of the last-copy top's daughters
    wb_parent = Counter()         # TTbarGen W and b are direct daughters of: t_TTbarGen / t_last / ...
    wd_parent = Counter()         # TTbarGen W daughters are direct daughters of: W_TTbarGen / W_last / ...
    # event dump bookkeeping
    dumped = []
    dump_ids = set()
    for tok in args.dump_ids.split(","):
        if tok.strip():
            r, l, e = (int(x) for x in tok.strip().split(":"))
            dump_ids.add((r, l, e))

    for ev in events:
        n_read += 1
        if wanted is not None:
            aux = ev.eventAuxiliary()
            if (aux.run(), aux.luminosityBlock(), aux.event()) not in wanted:
                continue
        if 0 <= args.max_events <= n_proc:
            break
        n_proc += 1
        if n_proc % args.progress == 0:
            print("  ... processed %d events" % n_proc)
            sys.stdout.flush()

        ev.getByLabel("prunedGenParticles", handle)
        gps = handle.product()
        pdg, st, moth, daus = read_event(gps)
        pruned = prune(pdg, st, moth, daus)
        try:
            sel = ttbargen(pruned)
        except TTbarGenFailure as e:
            fails[str(e)] += 1
            continue
        n_ok += 1

        p4_tg, p4_last = {}, {}
        lines = []
        for lab in LABELS:
            part = sel[lab]
            role = object_role(lab, part)
            orig = gps[part.orig]
            sf = orig.statusFlags()
            n_role[role] += 1
            status[role][part.status] += 1
            fl = [f for f in FLAGS if getattr(sf, f)()]
            for f in fl:
                flags[role][f] += 1

            li = last_copy(part.orig, pdg, daus)
            last = gps[li]
            if li == part.orig:
                is_last[role] += 1
            if last.statusFlags().isLastCopy():
                end_flag_ok[role] += 1
            if orig.pt() > 0:
                dpt[role].add((last.pt() - orig.pt()) / orig.pt())
            ddr[role].add(dR(orig, last))
            dm[role].add(last.mass() - orig.mass())
            if lab in ("Top", "Antitop"):
                p4_tg[lab] = orig.p4(); p4_last[lab] = last.p4()
            lines.append("  %-13s pdgId=%6d status=%3d pt=%8.2f | last copy: status=%3d pt=%8.2f | %s" % (
                lab, part.pdgId, part.status, orig.pt(), last.status(), last.pt(), ",".join(fl)))

        a = p4_tg["Top"] + p4_tg["Antitop"]
        b = p4_last["Top"] + p4_last["Antitop"]
        tt["mtt"].add(b.M() - a.M())
        tt["pttt"].add(b.Pt() - a.Pt())

        # ---- four-momentum consistency checks ----
        checks = []
        for tl, wl, bl, d1, d2 in (("Top", "WTop", "bTop", "Wdecay1", "Wdecay2"),
                                   ("Antitop", "WAntitop", "bAntitop", "WMinusdecay1", "WMinusdecay2")):
            t1 = sel[tl].orig
            tL = last_copy(t1, pdg, daus)
            w1, b1 = sel[wl].orig, sel[bl].orig
            wL = last_copy(w1, pdg, daus)
            wd = (sel[d1].orig, sel[d2].orig)

            dsig = tuple(sorted((abs(pdg[d]), st[d]) for d in daus[tL]))
            tlast_daus[dsig] += 1
            if w1 in daus[tL] and b1 in daus[tL]:
                wb_parent["W and b are daughters of t_last"] += 1
            elif w1 in daus[t1] and b1 in daus[t1]:
                wb_parent["W and b are daughters of t_TTbarGen"] += 1
            else:
                wb_parent["other (e.g. via intermediate copies)"] += 1
            if wd[0] in daus[wL] and wd[1] in daus[wL]:
                wd_parent["W daughters are daughters of W_last"] += 1
            elif wd[0] in daus[w1] and wd[1] in daus[w1]:
                wd_parent["W daughters are daughters of W_TTbarGen"] += 1
            else:
                wd_parent["other (e.g. via intermediate copies)"] += 1

            t1_4, tL_4 = p4_of(gps[t1]), p4_of(gps[tL])
            wb_4 = p4_add(p4_of(gps[w1]), p4_of(gps[b1]))
            dsum_4 = p4_add(*[p4_of(gps[d]) for d in daus[tL]]) if daus[tL] else (0.0, 0.0, 0.0, 0.0)
            wd_4 = p4_add(p4_of(gps[wd[0]]), p4_of(gps[wd[1]]))
            tag = "[W->qq']" if abs(pdg[wd[0]]) <= 5 else "[W->lnu]"
            checks.append(("%s: t_last vs sum(dau)" % tl, cmp_["t_last vs sum(daughters of t_last)"].add(tL_4, dsum_4)))
            checks.append(("%s: t_TTbarGen vs W+b" % tl, cmp_["t_TTbarGen vs W+b (TTbarGen)"].add(t1_4, wb_4)))
            checks.append(("%s: t_last vs W+b" % tl, cmp_["t_last vs W+b (TTbarGen)"].add(tL_4, wb_4)))
            checks.append(("%s: W_TTbarGen vs dau sum" % wl,
                           cmp_["W_TTbarGen vs dau1+dau2 %s" % tag].add(p4_of(gps[w1]), wd_4)))
            checks.append(("%s: W_last vs dau sum" % wl, cmp_["W_last vs dau1+dau2 %s" % tag].add(p4_of(gps[wL]), wd_4)))
            lines.append("  %s last copy daughters (|pdgId|, status): %s" % (tl, list(dsig)))

        if n_ok <= args.print_events:
            aux = ev.eventAuxiliary()
            print("\n--- run %d lumi %d event %d ---" % (aux.run(), aux.luminosityBlock(), aux.event()))
            print("\n".join(lines))
            for name, (cdm, cdpt, cdr) in checks:
                print("  check %-30s dm=%8.3f GeV  dpT/pT=%8.4f  dR=%7.4f" % (name, cdm, cdpt, cdr))

        # ---- optional diagram-oriented event dump ----
        if dump_ids or len(dumped) < args.dump_events:
            aux = ev.eventAuxiliary()
            eid = (aux.run(), aux.luminosityBlock(), aux.event())
            want = eid in dump_ids
            if not want and len(dumped) < args.dump_events and len(pdg) <= args.dump_max_particles:
                drs = [dR(gps[sel[tl].orig], gps[last_copy(sel[tl].orig, pdg, daus)]) for tl in ("Top", "Antitop")]
                want = max(drs) >= args.dump_min_dr
            if want:
                dumped.append(dump_event(aux, gps, pdg, st, moth, daus, sel))
                dump_ids.discard(eid)
            if not dump_ids and len(dumped) >= args.dump_events and args.dump_only:
                break

    # ---------------- event dumps ----------------
    if args.dump_events or args.dump_ids:
        import json
        with open(args.dump_json, "w") as f:
            json.dump(dumped, f, indent=1)
        print("\nWrote %d event dump(s) to %s" % (len(dumped), args.dump_json))
        if dump_ids:
            print("WARNING: requested events not found (check the files / --max-events -1): %s" % sorted(dump_ids))
        if args.dump_only:
            return

    # ---------------- summary ----------------
    print("\n" + "=" * 78)
    print("MiniAOD events read %d, processed %d, TTbarGen (on emulated UHH2 collection) OK in %d" % (
        n_read, n_proc, n_ok))
    for m, n in fails.most_common():
        print("  failure: %6d  %s" % (n, m))

    print("\nPythia status of the TTbarGen objects (cross-check with the UHH2 ntuple results):")
    for r in ROLES:
        if n_role[r]:
            print("  %-22s %s" % (r, ", ".join("%d: %.1f%%" % (s, 100.0 * c / n_role[r])
                                               for s, c in sorted(status[r].items()))))

    print("\nGenStatusFlags of the TTbarGen objects (fraction of objects with flag = true):")
    hdr = "  %-36s" % "flag" + "".join("%11s" % r.replace("W daughter: ", "Wd:") for r in ROLES if n_role[r])
    print(hdr)
    for f in FLAGS:
        row = "  %-36s" % f
        for r in ROLES:
            if n_role[r]:
                row += "%10.1f%%" % (100.0 * flags[r][f] / n_role[r])
        print(row)
    print("  %-36s" % "N objects" + "".join("%11d" % n_role[r] for r in ROLES if n_role[r]))

    print("\nComparison with the last copy (end of the same-pdgId copy chain in prunedGenParticles):")
    print("  %-22s %10s %12s %16s %14s %14s" % ("object", "TTbarGen=", "chain end", "(pT_last-pT)/pT", "dR",
                                                "m_last - m"))
    print("  %-22s %10s %12s %16s %14s %14s" % ("", "last copy", "isLastCopy", "mean +- rms", "mean +- rms",
                                                "mean +- rms"))
    for r in ROLES:
        if n_role[r]:
            print("  %-22s %9.1f%% %11.1f%% %8.4f+-%-7.4f %6.4f+-%-7.4f %6.3f+-%-7.3f" % (
                r, 100.0 * is_last[r] / n_role[r], 100.0 * end_flag_ok[r] / n_role[r],
                dpt[r].mean(), dpt[r].rms(), ddr[r].mean(), ddr[r].rms(), dm[r].mean(), dm[r].rms()))
    print("\n  ttbar system (last-copy tops minus TTbarGen tops):")
    print("    m_tt  difference [GeV]: mean %.3f  rms %.3f" % (tt["mtt"].mean(), tt["mtt"].rms()))
    print("    pT_tt difference [GeV]: mean %.3f  rms %.3f" % (tt["pttt"].mean(), tt["pttt"].rms()))

    print("\nFour-momentum consistency checks (other minus reference;")
    print("  'identical' = |dm| < %.2f GeV, |dpT/pT| < %g, dR < %g):" % (TOL_M, TOL_RELPT, TOL_DR))
    print_compare_header()
    for k in CMP_KEYS:
        print(cmp_[k].line(k))
    ntops = sum(wb_parent.values())
    print("\n  Mother of TTbarGen's W and b in prunedGenParticles (per top, N = %d):" % ntops)
    for k, n in wb_parent.most_common():
        print("    %-45s %6.1f%%" % (k, 100.0 * n / ntops))
    print("  Mother of TTbarGen's W daughters in prunedGenParticles (per W, N = %d):" % sum(wd_parent.values()))
    for k, n in wd_parent.most_common():
        print("    %-45s %6.1f%%" % (k, 100.0 * n / sum(wd_parent.values())))
    print("  Daughters of the last-copy top, (|pdgId|, status) - most common:")
    for sig, n in tlast_daus.most_common(8):
        print("    %6.2f%%  %s" % (100.0 * n / ntops, list(sig)))
    print("\n  Hypothesis under test: TTbarGen's W (22) and b (23) are the decay products of the LAST-copy top")
    print("  (status 62). Then 't_last vs sum(dau)' and 't_last vs W+b' are identical in ~100% of tops,")
    print("  't_TTbarGen vs W+b' has the same mass but is identical in ~0%, and 'W_last vs dau sum' is")
    print("  identical in ~100% while 'W_TTbarGen vs dau sum' is identical only when W_TTbarGen is the last copy (~24%).")


if __name__ == "__main__":
    main()