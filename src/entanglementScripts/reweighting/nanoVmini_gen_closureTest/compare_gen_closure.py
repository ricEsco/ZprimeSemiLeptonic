# compare_gen_closure.py -- match the two gen dumps by (run,lumi,event) and report residuals.
#
#   python compare_gen_closure.py  gendump_uhh2.root  gendump_nano.root
#
# A correct closure: residuals at the NanoAOD float-precision level (~1e-3) with NO population
# of large-residual outliers, and no channel mis-classifications. A cluster of large residuals
# (or channel flips), especially in one decay channel, flags a genuine selection difference.
from __future__ import print_function
import sys
import ROOT
ROOT.gROOT.SetBatch(True)

# what counts as a "selection-level" discrepancy (well above NanoAOD rounding)
TOL_MTT_REL = 0.01     # 1% on m_tt
TOL_ANG     = 0.05     # absolute on cosTheta*, cHel, cHel_P3n

CH = {0: "had", 1: "ehad", 2: "muhad", 3: "tauhad", 4: "ee", 5: "mumu",
      6: "tautau", 7: "emu", 8: "etau", 9: "mutau", 10: "notfound"}


def load(fn):
    f = ROOT.TFile.Open(fn); t = f.Get("gendump")
    d = {}
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        d[(int(t.run), int(t.lumi), int(t.event))] = (
            int(t.gen_channel), int(t.gen_semilep), float(t.gen_mtt),
            float(t.gen_cosThetaStar), float(t.gen_cHel), float(t.gen_cHel_P3n), float(t.genweight))
    f.Close()
    return d


def pct(vals, q):
    if not vals:
        return 0.0
    s = sorted(vals)
    return s[min(len(s) - 1, int(q * len(s)))]


def main():
    uhh2 = load(sys.argv[1])            # dict keyed by (run,lumi,event)
    nano = load(sys.argv[2])
    print("loaded: uhh2=%d  nano=%d" % (len(uhh2), len(nano)))

    matched = both_semilep = chan_mismatch = semilep_mismatch = 0
    res = {"mtt": [], "cosTS": [], "cHel": [], "cHelP3n": [], "gw": []}
    outliers = {}          # channel -> count of selection-level discrepancies
    outlier_examples = []

    for key, u in uhh2.items():
        n = nano.get(key)
        if n is None:
            continue
        matched += 1
        uc, us, umtt, ucts, uch, ucp, ugw = u
        nc, ns, nmtt, ncts, nch, ncp, ngw = n
        if us != ns:                   # do the two codes agree it IS semileptonic? (key check)
            semilep_mismatch += 1
        if not (us and ns):            # residuals + channel check only where both are semileptonic
            continue
        both_semilep += 1
        if uc != nc:                   # do they agree on ehad/muhad/tauhad?
            chan_mismatch += 1
        d_mtt = abs(nmtt - umtt) / umtt if umtt else 0.0
        d_cts = abs(ncts - ucts)
        d_ch  = abs(nch - uch)
        d_cp  = abs(ncp - ucp)
        d_gw  = abs(ngw - ugw) / abs(ugw) if ugw else 0.0
        res["mtt"].append(d_mtt); res["cosTS"].append(d_cts)
        res["cHel"].append(d_ch); res["cHelP3n"].append(d_cp); res["gw"].append(d_gw)
        if d_mtt > TOL_MTT_REL or d_cts > TOL_ANG or d_ch > TOL_ANG or d_cp > TOL_ANG:
            outliers[uc] = outliers.get(uc, 0) + 1
            if len(outlier_examples) < 10:
                outlier_examples.append((key, CH.get(uc, uc), CH.get(nc, nc),
                                         umtt, nmtt, uch, nch, ucp, ncp))

    print("\n=== coverage ===")
    print("  matched events (in both)        : %d" % matched)
    print("  unmatched uhh2 (no nano partner): %d" % (len(uhh2) - matched))
    print("\n=== classification agreement (matched) ===")
    print("  gen_channel mismatches : %d  (%.3f%%)" % (chan_mismatch, 100.0 * chan_mismatch / matched if matched else 0))
    print("  gen_semilep mismatches : %d  (%.3f%%)" % (semilep_mismatch, 100.0 * semilep_mismatch / matched if matched else 0))
    print("\n=== residuals over %d both-semileptonic matches (|nano - uhh2|) ===" % both_semilep)
    print("  %-9s %10s %10s %10s" % ("quantity", "median", "99th pct", "max"))
    for k, lbl in (("mtt", "m_tt (rel)"), ("cosTS", "cosTheta*"), ("cHel", "cHel"),
                   ("cHelP3n", "cHel_P3n"), ("gw", "genweight (rel)")):
        v = res[k]
        print("  %-9s %10.2e %10.2e %10.2e" % (lbl, pct(v, 0.50), pct(v, 0.99), max(v) if v else 0.0))

    nout = sum(outliers.values())
    print("\n=== selection-level discrepancies (m_tt>%.0f%% or angle>%.2f) ===" % (100 * TOL_MTT_REL, TOL_ANG))
    print("  total: %d / %d  (%.3f%%)" % (nout, both_semilep, 100.0 * nout / both_semilep if both_semilep else 0))
    for ch, c in sorted(outliers.items(), key=lambda kv: -kv[1]):
        print("    channel %-8s : %d" % (CH.get(ch, ch), c))
    for ex in outlier_examples:
        print("    e.g. id=%s  chan u/n=%s/%s  mtt u/n=%.1f/%.1f  cHel u/n=%.3f/%.3f  cHelP3n u/n=%.3f/%.3f" % ex)

    verdict = "CLOSURE OK" if (chan_mismatch == 0 and nout == 0) else \
              "INVESTIGATE (discrepancies above -- likely tau/pruning; check the channel breakdown)"
    print("\n%s" % verdict)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: python compare_gen_closure.py gendump_uhh2.root gendump_nano.root"); sys.exit(1)
    main()
