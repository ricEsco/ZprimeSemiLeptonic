#!/usr/bin/env python
"""
ttbargen_status_check.py

Prints the Pythia status codes of the gen-particles that UHH2's TTbarGen selects for the
ttbar system, read from the GenParticles branch of UHH2 AnalysisTree output files.
The selection itself is done by your validated TTbarGenPy class (ttbargen_py.py), which
must sit in the same directory as this script (or on PYTHONPATH).

Run inside your UHH2 CMSSW environment (cmsenv / singularity):
  python ttbargen_status_check.py file1.root [file2.root ...] \
      [--tree AnalysisTree] [--branch GenParticles] [--max-events 2000] \
      [--print-events 10] [--dump 1]
"""
from __future__ import print_function

import argparse
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ttbargen_py import TTbarGenPy

CHANNEL_NAMES = {
    TTbarGenPy.e_had: "had", TTbarGenPy.e_ehad: "ehad", TTbarGenPy.e_muhad: "muhad",
    TTbarGenPy.e_tauhad: "tauhad", TTbarGenPy.e_ee: "ee", TTbarGenPy.e_mumu: "mumu",
    TTbarGenPy.e_tautau: "tautau", TTbarGenPy.e_emu: "emu", TTbarGenPy.e_etau: "etau",
    TTbarGenPy.e_mutau: "mutau", TTbarGenPy.e_notfound: "notfound",
}

# (label, TTbarGenPy member)
OBJECTS = [
    ("Top", "m_Top"), ("Antitop", "m_Antitop"),
    ("WTop", "m_WTop"), ("WAntitop", "m_WAntitop"),
    ("bTop", "m_bTop"), ("bAntitop", "m_bAntitop"),
    ("Wdecay1", "m_Wdecay1"), ("Wdecay2", "m_Wdecay2"),
    ("WMinusdecay1", "m_WMinusdecay1"), ("WMinusdecay2", "m_WMinusdecay2"),
]
LABELS = [lab for lab, _ in OBJECTS]


def setup_root():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    # libSFrameCore (auto-loaded while reading the file's streamer info) needs
    # TProofOutputFile from libProof. Load the PROOF/TreePlayer libs *first* so that
    # the auto-load does not fail and segfault.
    for lib in ("libTreePlayer", "libProof", "libProofPlayer"):
        ROOT.gSystem.Load(lib)
    return ROOT


def ensure_genparticle_dict(ROOT):
    """Called after the first file is open. Normally the dictionary for GenParticle is
    auto-loaded by then; only load libSUHH2core explicitly if it is still missing."""
    cl = ROOT.TClass.GetClass("GenParticle")
    if cl and cl.IsLoaded():
        return
    for lib in ("libSFrameCore", "libSUHH2core"):
        ROOT.gSystem.Load(lib)
    cl = ROOT.TClass.GetClass("GenParticle")
    if not (cl and cl.IsLoaded()):
        sys.exit("ERROR: no dictionary for GenParticle. Run inside your UHH2 environment (cmsenv).")


# ---------------------------------------------------------------------------
# Four-momentum consistency checks (pure python, (E, px, py, pz) tuples)
# ---------------------------------------------------------------------------
TOL_M, TOL_RELPT, TOL_DR = 0.05, 1e-3, 1e-3      # "identical" = all three within tolerance


def p4_of(gp):
    pt, eta, phi, e = gp.pt(), gp.eta(), gp.phi(), gp.energy()
    return (e, pt * math.cos(phi), pt * math.sin(phi), pt * math.sinh(eta))


def p4_add(*vs):
    return tuple(sum(c) for c in zip(*vs))


def p4_m(v):
    m2 = v[0] ** 2 - v[1] ** 2 - v[2] ** 2 - v[3] ** 2
    return math.copysign(math.sqrt(abs(m2)), m2)


def p4_pt(v):
    return math.hypot(v[1], v[2])


def p4_etaphi(v):
    pt = p4_pt(v)
    eta = math.asinh(v[3] / pt) if pt > 0 else math.copysign(1e9, v[3])
    return eta, math.atan2(v[2], v[1])


def p4_dr(a, b):
    ea, pa = p4_etaphi(a)
    eb, pb = p4_etaphi(b)
    dphi = abs(pa - pb)
    if dphi > math.pi:
        dphi = 2 * math.pi - dphi
    return math.hypot(ea - eb, dphi)


class Stat(object):
    def __init__(self):
        self.n = 0; self.s = 0.0; self.s2 = 0.0
    def add(self, x):
        self.n += 1; self.s += x; self.s2 += x * x
    def mean(self):
        return self.s / self.n if self.n else float("nan")
    def rms(self):
        return math.sqrt(max(self.s2 / self.n - self.mean() ** 2, 0.0)) if self.n else float("nan")


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
            return "  %-34s (no entries)" % name
        return "  %-34s %8d %9.1f%% %9.1f%%  %8.3f+-%-7.3f %8.4f+-%-7.4f %7.4f+-%-7.4f" % (
            name, self.n, 100.0 * self.same_m / self.n, 100.0 * self.identical / self.n,
            self.dm.mean(), self.dm.rms(), self.dpt.mean(), self.dpt.rms(), self.dr.mean(), self.dr.rms())


def print_compare_header():
    print("  %-34s %8s %10s %10s  %16s %16s %16s" % ("comparison (reference vs sum)", "N", "same m", "identical",
                                                    "dm [GeV]", "dpT/pT", "dR"))
    print("  %-34s %8s %10s %10s  %16s %16s %16s" % ("", "", "", "", "mean +- rms", "mean +- rms", "mean +- rms"))


def all_daughters(parent, genparticles):
    """All particles listing `parent` as mother1/mother2. Shows W daughters beyond the two
    TTbarGen uses. Compares the stored mother indices directly: GenParticle::mother() does a
    linear search of the whole list per call, which makes this O(N^2) per event and very slow
    when the GenParticles branch also holds GenJet constituents / stable particles."""
    out = []
    pidx = parent.index()
    for gp in genparticles:
        if gp.mother1() == pidx or gp.mother2() == pidx:
            out.append(gp)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--tree", default="AnalysisTree")
    ap.add_argument("--branch", default="GenParticles")
    ap.add_argument("--max-events", type=int, default=5000, help="total events to process (all files); -1 = all")
    ap.add_argument("--progress", type=int, default=500, help="print a progress line every N events")
    ap.add_argument("--no-extra", action="store_true", help="skip the (slower) scan for extra W daughters")
    ap.add_argument("--print-events", type=int, default=10, help="print selected particles for the first N good events")
    ap.add_argument("--dump", type=int, default=0, help="dump the full GenParticles list for the first N events")
    args = ap.parse_args()

    ROOT = setup_root()

    status_count = {lab: Counter() for lab in LABELS}
    combo_count = Counter()
    channel_count = Counter()
    extra_wdau = Counter()
    cmp_top = Compare()                                  # TTbarGen top  vs  W + b
    cmp_w = {"W -> l nu": Compare(), "W -> q q'": Compare()}   # TTbarGen W vs its two daughters
    n_proc = n_ok = 0
    dict_checked = False

    for fname in args.files:
        f = ROOT.TFile.Open(fname)
        if not f or f.IsZombie():
            print("WARNING: cannot open %s, skipping" % fname)
            continue
        if not dict_checked:
            ensure_genparticle_dict(ROOT)
            dict_checked = True
        t = f.Get(args.tree)
        if not t:
            print("WARNING: no tree '%s' in %s, skipping" % (args.tree, fname))
            f.Close()
            continue

        for ievt in range(t.GetEntries()):
            if 0 <= args.max_events <= n_proc:
                break
            t.GetEntry(ievt)
            n_proc += 1
            if n_proc % args.progress == 0:
                print("  ... processed %d events" % n_proc)
                sys.stdout.flush()
            genparticles = getattr(t, args.branch)

            if n_proc <= args.dump:
                print("\n=== %s event %d: full %s list (%d particles) ===" % (
                    os.path.basename(fname), ievt, args.branch, genparticles.size()))
                for j in range(genparticles.size()):
                    g = genparticles[j]
                    print("  idx=%3d pdgId=%6d status=%3d mothers=(%5d,%5d) daughters=(%5d,%5d) pt=%8.2f" % (
                        g.index(), g.pdgId(), g.status(), g.mother1(), g.mother2(),
                        g.daughter1(), g.daughter2(), g.pt()))

            ttgen = TTbarGenPy(genparticles, throw_on_failure=False)
            ch = ttgen.DecayChannel()
            channel_count[CHANNEL_NAMES[ch]] += 1
            if ch == TTbarGenPy.e_notfound:
                continue
            n_ok += 1

            parts = dict((lab, getattr(ttgen, attr)) for lab, attr in OBJECTS)
            statuses = tuple(parts[lab].status() for lab in LABELS)
            for lab, st in zip(LABELS, statuses):
                status_count[lab][st] += 1
            combo_count[statuses] += 1

            for wlab, d1, d2 in (() if args.no_extra else
                                 (("WTop", "Wdecay1", "Wdecay2"),
                                  ("WAntitop", "WMinusdecay1", "WMinusdecay2"))):
                used = (parts[d1].index(), parts[d2].index())
                for d in all_daughters(parts[wlab], genparticles):
                    if d.index() not in used:
                        extra_wdau[(d.pdgId(), d.status())] += 1

            # ---- four-momentum consistency checks ----
            checks = []
            for tl, wl, bl, d1, d2 in (("Top", "WTop", "bTop", "Wdecay1", "Wdecay2"),
                                       ("Antitop", "WAntitop", "bAntitop", "WMinusdecay1", "WMinusdecay2")):
                t4, w4, b4 = p4_of(parts[tl]), p4_of(parts[wl]), p4_of(parts[bl])
                res = cmp_top.add(t4, p4_add(w4, b4))
                checks.append(("%s vs %s+%s" % (tl, wl, bl), res))
                key = "W -> q q'" if abs(parts[d1].pdgId()) <= 5 else "W -> l nu"
                res = cmp_w[key].add(w4, p4_add(p4_of(parts[d1]), p4_of(parts[d2])))
                checks.append(("%s vs %s+%s" % (wl, d1, d2), res))

            if n_ok <= args.print_events:
                print("\n--- %s event %d  channel=%s ---" % (os.path.basename(fname), ievt, CHANNEL_NAMES[ch]))
                for lab in LABELS:
                    g = parts[lab]
                    print("  %-13s pdgId=%6d status=%3d idx=%3d pt=%8.2f eta=%6.2f m=%8.2f" % (
                        lab, g.pdgId(), g.status(), g.index(), g.pt(), g.eta(), p4_m(p4_of(g))))
                for name, (dm, dpt, dr) in checks:
                    print("  check %-36s dm=%8.3f GeV  dpT/pT=%8.4f  dR=%7.4f" % (name, dm, dpt, dr))

        f.Close()
        if 0 <= args.max_events <= n_proc:
            break

    # ---------------- summary ----------------
    print("\n" + "=" * 70)
    print("Processed %d events, TTbarGen found the ttbar system in %d" % (n_proc, n_ok))
    print("\nDecay channels:")
    for ch, n in channel_count.most_common():
        print("  %-9s %6d" % (ch, n))
    print("\nStatus codes per TTbarGen object (status: count):")
    for lab in LABELS:
        print("  %-13s %s" % (lab, ", ".join("%d: %d" % (s, n) for s, n in sorted(status_count[lab].items()))))
    print("\nMost common status combinations")
    print("  (%s):" % " ".join(LABELS))
    for combo, n in combo_count.most_common(10):
        print("  %6d  %s" % (n, combo))
    print("\nExtra W daughters ignored by TTbarGen (pdgId, status): count")
    if extra_wdau:
        for (pid, st), n in extra_wdau.most_common(20):
            print("  (%6d, %3d): %d" % (pid, st, n))
    else:
        print("  none")

    print("\nFour-momentum consistency of the TTbarGen objects (sum minus reference;")
    print("  'identical' = |dm| < %.2f GeV, |dpT/pT| < %g, dR < %g):" % (TOL_M, TOL_RELPT, TOL_DR))
    print_compare_header()
    print(cmp_top.line("top (TTbarGen) vs W + b"))
    for key in ("W -> l nu", "W -> q q'"):
        print(cmp_w[key].line("W (TTbarGen) vs daughters, %s" % key))
    print("\n  Interpretation: if TTbarGen's W and b are the decay products of the *last-copy* top while the")
    print("  top itself is the first (status-22) copy, then 'top vs W + b' has the same mass in ~100% of")
    print("  events but is identical in ~0% (W + b = last-copy top). Likewise 'W vs daughters' is identical")
    print("  only when the W was not copied in the shower (~24% of events in the MiniAOD flags study).")


if __name__ == "__main__":
    main()