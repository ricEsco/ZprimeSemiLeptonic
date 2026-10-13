"""
plotSpinAnalyzerBoostShift.py

How much does the DIRECTION of each spin analyser change when the boosts in gen_observables() (makeAnalysisSkims.py)
use the first-copy tops (TTbarGen Top()/Antitop(), gen_ branches) instead of the last-copy tops (W + b, genLC_ branches)?

The two boost steps of gen_observables() are reproduced exactly, once per top definition:
  step 1 (CM):  boost lepton and hadronic b (and both tops) by -(t + tbar).BoostVector()        -> ttbar CM frame
  step 2 (RF):  boost the CM-frame lepton by -t_lep(CM).BoostVector() and the CM-frame b_had by -t_had(CM).BoostVector()
The lepton and b_had lab four-vectors are identical in both definitions; only the tops (and therefore both boosts) differ.

For the SAME object it then compares the direction obtained with first-copy (FC) and last-copy (LC) boosts via the
3D opening angle between the two unit vectors (what enters cos(theta) = SA . axis), for
  * spin analyser 1 = decay product of the TOP     (lepton if the top decays leptonically, else b_had)
  * spin analyser 2 = decay product of the ANTITOP (as in the cosTheta1* / cosTheta2* projections)
plus, as a reference, the change of the helicity axis k = top direction in the CM frame.

Leptonic/hadronic top assignment: per event from the skim itself. The last-copy top is exactly W + b, so
m(t_LC - b_had) equals the W mass for the hadronic top only; the side closer to 80.4 GeV is taken as hadronic.
The agreement with the gen_Lepton_charge branch (either sign convention) is printed as a cross-check.

Inputs: SC_SkimTree from makeAnalysisSkims.py with DO_LASTCOPY = True (gen_/genLC_ Top, Antitop; gen_ Lepton, bHad).

Usage (CMSSW_10_6_28 / python 2.7 / ROOT 6.14, numpy):
  python plotSpinAnalyzerBoostShift.py TTToSemiLeptonicUL18_SCskim_*.root --out SpinAnalyzerBoostShift
  python plotSpinAnalyzerBoostShift.py skims*.root --cut "reco_is_toptag_reconstruction>0.5" --out SAshift_merged
  python plotSpinAnalyzerBoostShift.py skims*.root --max-entries 500000          # quick look
Outputs: <out>.pdf, <out>.png, <out>_summary.txt
"""
from __future__ import print_function, division

import argparse
import math
from array import array

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)

MW = 80.4
# C_CM, C_RF, C_K = "#2a78d6", "#eb6834", "#1baf7a"   # validated categorical slots 1-3; also distinguished by line style + legend
C_CM, C_RF, C_K = "#8000ff", "#FF7F00", "#830000"   # colors: purple, orange, dark red (distinguishable in B/W printouts)

PARTICLES = ["Top", "Antitop", "Lepton", "bHad"]


# ------------------------------------------------------------------ reading
def read_chunk(tree, exprs, sel, first, n):
    """Columns for entries [first, first+n) passing sel, via TTree::Draw (4 columns per call, same row order)."""
    tree.SetEstimate(n + 1)
    out = []
    for k in range(0, len(exprs), 4):
        chunk = exprs[k:k + 4]
        m = tree.Draw(":".join(chunk), sel, "goff", n, first)
        if m < 0:
            raise RuntimeError("TTree::Draw failed for %s (check branch names)" % chunk)
        for j in range(len(chunk)):
            if m == 0:
                out.append(np.zeros(0))
                continue
            buf = getattr(tree, "GetV%d" % (j + 1))()
            buf.SetSize(m)
            out.append(np.array(buf, dtype=np.float64, copy=True))
    return out


# ------------------------------------------------------------------ kinematics (same formulas as TLorentzVector)
def p4_from(pt, eta, phi, E):
    return np.stack([E, pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)], axis=1)


def boost(v, b):
    """TLorentzVector::Boost(b) for arrays: v (N,4) = (E,px,py,pz), b (N,3)."""
    b2 = np.sum(b * b, axis=1)
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(b * v[:, 1:], axis=1)
    g2 = np.where(b2 > 0, (g - 1.0) / np.where(b2 > 0, b2, 1.0), 0.0)
    return np.column_stack([g * (v[:, 0] + bp), v[:, 1:] + ((g2 * bp + g * v[:, 0])[:, None]) * b])


def beta(v):
    return v[:, 1:] / v[:, :1]


def mass(v):
    m2 = v[:, 0] ** 2 - np.sum(v[:, 1:] ** 2, axis=1)
    return np.sign(m2) * np.sqrt(np.abs(m2))


def angle3d(a, b):
    """opening angle between the 3-vectors of a and b (atan2 form, precise for tiny angles)."""
    a3, b3 = a[:, 1:], b[:, 1:]
    return np.arctan2(np.linalg.norm(np.cross(a3, b3), axis=1), np.sum(a3 * b3, axis=1))


def gen_observables_frames(top, atop, lep, bhad, top_is_lep):
    """Exactly the boost chain of gen_observables(): returns CM-frame and parent-rest-frame lepton / b_had and k."""
    tt = top + atop
    bcm = -beta(tt)
    top_cm, atop_cm = boost(top, bcm), boost(atop, bcm)
    lep_cm, bhad_cm = boost(lep, bcm), boost(bhad, bcm)
    tlep_cm = np.where(top_is_lep[:, None], top_cm, atop_cm)
    thad_cm = np.where(top_is_lep[:, None], atop_cm, top_cm)
    lep_rf = boost(lep_cm, -beta(tlep_cm))
    bhad_rf = boost(bhad_cm, -beta(thad_cm))
    return dict(lep_cm=lep_cm, bhad_cm=bhad_cm, lep_rf=lep_rf, bhad_rf=bhad_rf, k=top_cm)


# ------------------------------------------------------------------ histogram helpers
def color(h, c, ls):
    col = ROOT.TColor.GetColor(c)
    h.SetLineColor(col); h.SetLineWidth(2); h.SetLineStyle(ls); h.SetFillStyle(0)


def fill(h, x, w):
    if len(x):
        h.FillN(len(x), array("d", x.tolist()), array("d", w.tolist()))


def quantile(h, q):
    tot = h.Integral(0, h.GetNbinsX() + 1)
    if tot <= 0:
        return float("nan")
    acc = h.GetBinContent(0)
    for i in range(1, h.GetNbinsX() + 1):
        c = h.GetBinContent(i)
        if acc + c >= q * tot:
            lo, wd = h.GetXaxis().GetBinLowEdge(i), h.GetXaxis().GetBinWidth(i)
            return lo + (q * tot - acc) / c * wd if c > 0 else lo
        acc += c
    return h.GetXaxis().GetXmax()


def mean(h):
    tot = h.Integral(0, h.GetNbinsX() + 1)
    return sum(h.GetBinContent(i) * h.GetXaxis().GetBinCenter(i) for i in range(1, h.GetNbinsX() + 1)) / tot if tot else float("nan")


def mrad(x):
    return "%.1f mrad" % (1e3 * x)


# ------------------------------------------------------------------ main
# key, legend label, colour, line style
QUANTS = [
    ("sa1_cm", "in t#bar{t} CM frame",  C_CM, 1),
    ("sa1_rf", "in t rest frame",       C_RF, 1),
    ("sa2_cm", "in t#bar{t} CM frame",  C_CM, 1),
    ("sa2_rf", "in #bar{t} rest frame", C_RF, 1),
    ("k",      "k in t#bar{t} CM frame", C_K, 1),
]
QLABEL = dict((q[0], q) for q in QUANTS)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--tree", default="SC_SkimTree")
    ap.add_argument("--cut", default="1", help='extra selection, e.g. "reco_is_toptag_reconstruction>0.5"')
    ap.add_argument("--weight", default="genweight")
    ap.add_argument("--max-entries", type=int, default=-1, help="entries read per file; -1 = all")
    ap.add_argument("--chunk", type=int, default=1000000, help="entries per TTree::Draw chunk (memory)")
    ap.add_argument("--xmax", type=float, default=0.6, help="upper x range of the plots [rad]")
    ap.add_argument("--out", default="SpinAnalyzerBoostShift")
    args = ap.parse_args()

    sel = "(gen_semilep == 1) && (gen_Top_pt > 0) && (genLC_Top_pt > 0) && (%s)" % args.cut
    exprs = ["%s_%s_%s" % (pre, p, c) for pre, p in [("gen", x) for x in PARTICLES] + [("genLC", "Top"), ("genLC", "Antitop")]
             for c in ("pt", "eta", "phi", "E")] + [args.weight, "gen_Lepton_charge"]

    # display histograms (100 bins) + fine histograms (exact means / quantiles)
    H = {}
    for key, _, _, _ in QUANTS:
        H[key] = ROOT.TH1D("h_" + key, "", 100, 0.0, args.xmax)
        H[key + "_fine"] = ROOT.TH1D("f_" + key, "", 20000, 0.0, math.pi)
        H[key].SetDirectory(0); H[key + "_fine"].SetDirectory(0)

    n_sel, n_agree, sum_w, sum_w_toplep, printed = 0, 0, 0.0, 0.0, False
    for fn in args.files:
        f = ROOT.TFile.Open(fn)
        t = f.Get(args.tree)
        n_tot = t.GetEntries() if args.max_entries < 0 else min(t.GetEntries(), args.max_entries)
        for first in range(0, n_tot, args.chunk):
            n = min(args.chunk, n_tot - first)
            cols = read_chunk(t, exprs, sel, first, n)
            if not len(cols[-1]):
                continue
            v = [p4_from(*cols[4 * i:4 * i + 4]) for i in range(6)]
            top_fc, atop_fc, lep, bhad, top_lc, atop_lc = v
            w, qbranch = cols[-2], cols[-1]

            # leptonic/hadronic assignment from the W-mass test on the last-copy tops (exact: t_LC = W + b)
            m_top_minus_b = mass(top_lc - bhad)      # = m_W if the TOP is the hadronic one
            m_atop_minus_b = mass(atop_lc - bhad)    # = m_W if the ANTITOP is the hadronic one
            top_is_lep = np.abs(m_atop_minus_b - MW) < np.abs(m_top_minus_b - MW)
            n_agree += int(np.sum(top_is_lep == (qbranch < 0)))     # convention: Lepton_charge = -1 for l+
            sum_w += float(np.sum(w)); sum_w_toplep += float(np.sum(w[top_is_lep]))

            res = {}
            for tag, (tp, at) in (("FC", (top_fc, atop_fc)), ("LC", (top_lc, atop_lc))):
                fr = gen_observables_frames(tp, at, lep, bhad, top_is_lep)
                sel_l = top_is_lep[:, None]
                # spin analyser 1 = decay product of the TOP, 2 = of the ANTITOP (as in cosTheta1* / cosTheta2*)
                res[tag] = dict(sa1_cm=np.where(sel_l, fr["lep_cm"], fr["bhad_cm"]),
                                sa1_rf=np.where(sel_l, fr["lep_rf"], fr["bhad_rf"]),
                                sa2_cm=np.where(sel_l, fr["bhad_cm"], fr["lep_cm"]),
                                sa2_rf=np.where(sel_l, fr["bhad_rf"], fr["lep_rf"]),
                                k=fr["k"])
            for key, _, _, _ in QUANTS:
                ang = angle3d(res["FC"][key], res["LC"][key])
                fill(H[key], ang, w); fill(H[key + "_fine"], ang, w)

            if not printed:
                printed = True
                mw_right = np.where(top_is_lep, m_atop_minus_b, m_top_minus_b)
                print("Self-check on the first chunk of %s (%d events):" % (fn, len(w)))
                print("  m(t_had,LC - b_had): median %.2f GeV, %.1f%% within 80.4 +- 15 GeV (must be ~W line shape)"
                      % (np.median(mw_right), 100.0 * np.mean(np.abs(mw_right - MW) < 15)))
                for i in range(min(3, len(w))):
                    angs = tuple(angle3d(res["FC"][k][i:i + 1], res["LC"][k][i:i + 1])[0]
                                 for k in ("sa1_cm", "sa1_rf", "sa2_cm", "sa2_rf"))
                    print("  event %d: SA1 = %-6s  angle(FC,LC) [rad]: SA1 CM %.4f RF %.4f   SA2 CM %.4f RF %.4f"
                          % ((i, "lepton" if top_is_lep[i] else "b_had") + angs))
            n_sel += len(w)
        f.Close()
        print("  %s done (%d selected events so far)" % (fn, n_sel))

    agree = n_agree / float(max(n_sel, 1))
    f_toplep = sum_w_toplep / sum_w if sum_w else float("nan")
    print("Leptonic-top assignment (W-mass test) agrees with gen_Lepton_charge(-1 = l+) in %.2f%% of events" % (100 * agree))

    # ------------------------------------------------------------------ plot
    canv = ROOT.TCanvas("c", "c", 1800, 620)
    summary = ["Selection: %s   Weight: %s   Events: %d" % (sel, args.weight, n_sel),
               "Leptonic-top assignment from m(t_LC - b_had) = m_W; agreement with gen_Lepton_charge (-1 = l+): %.2f%%" % (100 * agree),
               "SA1 (from t) is the lepton in %.1f%% of events, SA2 (from tbar) in %.1f%% (weighted)" % (100 * f_toplep, 100 * (1 - f_toplep)),
               "", "3D angle between the FC- and LC-boosted directions of the same object:"]
    keep = []

    canv.Divide(3, 1)
    panels = [("Spin analyser from t",       ["sa1_cm", "sa1_rf"], "#Delta#theta(SA_{1}^{FC}, SA_{1}^{LC}) [rad]", " (lepton %.0f%%, b_{had} %.0f%%)" % (100 * f_toplep, 100 * (1 - f_toplep))),
              ("Spin analyser from #bar{t}", ["sa2_cm", "sa2_rf"], "#Delta#theta(SA_{2}^{FC}, SA_{2}^{LC}) [rad]", " (lepton %.0f%%, b_{had} %.0f%%)" % (100 * (1 - f_toplep), 100 * f_toplep)),
              ("Reference: helicity axis k", ["k"], "#Delta#theta(k^{FC}, k^{LC}) [rad]", "")]
    


    # # PLotting one transformation at a time
    # canv.Divide(2, 1)
    # canv.SetCanvasSize(1200, 620) # Width, Height

    # # Plotting CM frame only for clarity, as the RF frame is commented out in the QUANTS list
    # panels = [("Spin analyser from t",       ["sa1_cm"], "#Delta#theta(SA_{1}^{FC}, SA_{1}^{LC}) [rad]", " (lepton %.0f%%, b_{had} %.0f%%)" % (100 * f_toplep, 100 * (1 - f_toplep))),
    #           ("Spin analyser from #bar{t}", ["sa2_cm"], "#Delta#theta(SA_{2}^{FC}, SA_{2}^{LC}) [rad]", " (lepton %.0f%%, b_{had} %.0f%%)" % (100 * (1 - f_toplep), 100 * f_toplep))]

    # # Plotting RF frame only for clarity, as the CM frame is commented out in the QUANTS list
    # panels = [("Spin analyser from t",       ["sa1_rf"], "#Delta#theta(SA_{1}^{FC}, SA_{1}^{LC}) [rad]", " (lepton %.0f%%, b_{had} %.0f%%)" % (100 * f_toplep, 100 * (1 - f_toplep))),
    #           ("Spin analyser from #bar{t}", ["sa2_rf"], "#Delta#theta(SA_{2}^{FC}, SA_{2}^{LC}) [rad]", " (lepton %.0f%%, b_{had} %.0f%%)" % (100 * (1 - f_toplep), 100 * f_toplep))]
    




    for ip, (ptitle, keys, xtitle, mix) in enumerate(panels):
        pad = canv.cd(ip + 1)
        pad.SetLeftMargin(0.14); pad.SetRightMargin(0.03); pad.SetTopMargin(0.12); pad.SetBottomMargin(0.13)
        hs = []
        for key in keys:
            h = H[key].Clone("d_" + key)
            tot = H[key + "_fine"].Integral(0, H[key + "_fine"].GetNbinsX() + 1)
            if tot > 0:
                h.Scale(1.0 / tot)
            q = QLABEL[key]
            color(h, q[2], q[3])
            hs.append((key, h, q))
        ymax = 1.25 * max(h.GetMaximum() for _, h, _ in hs)
        for i, (key, h, q) in enumerate(hs):
            h.SetMaximum(ymax); h.SetMinimum(0)
            h.GetXaxis().SetTitle(xtitle); h.GetYaxis().SetTitle("norm. events")
            h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
            h.GetXaxis().SetLabelSize(0.040); h.GetYaxis().SetLabelSize(0.040)
            h.GetYaxis().SetTitleOffset(1.5); h.GetXaxis().SetTitleOffset(1.2)
            h.Draw("hist" if i == 0 else "hist same")
        # legend with mean, median, 90% quantile
        leg = ROOT.TLegend(0.18, 0.86 - 0.08 * len(hs), 0.97, 0.86)
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.036)
        for key, h, q in hs:
            hf = H[key + "_fine"]
            med, q90, mu = quantile(hf, 0.5), quantile(hf, 0.9), mean(hf)
            leg.AddEntry(h, q[1], "l")
            blank = ROOT.TH1D("blank_" + key, "", 1, 0, 1); blank.SetLineColor(0); keep.append(blank)
            leg.AddEntry(blank, "#mu = %.0f, med. = %.0f, q_{90} = %.0f mrad" % (1e3 * mu, 1e3 * med, 1e3 * q90), "")
            summary.append("  %-7s %-24s mean %-11s median %-11s 90%% quantile %s"
                           % (key, q[1].replace("#bar{t}", "tbar").replace("t#", "t"), mrad(mu), mrad(med), mrad(q90)))
        leg.Draw()
        # title
        tl = ROOT.TLatex(); tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.048)
        tl.DrawLatex(0.14, 0.91, ptitle+mix)
        keep.extend([h for _, h, _ in hs] + [leg, tl])
    canv.cd(0)
    hdr = ROOT.TLatex(); hdr.SetNDC(); hdr.SetTextFont(42); hdr.SetTextSize(0.026)
    hdr.DrawLatex(0.005, 0.965, "Spin Analyzers boosted with FC (TTbarGen Top()) vs LC (W + b) tops: %d events,  selection: %s" % (n_sel, args.cut))
    keep.append(hdr)

    # Save plots
    canv.Print(args.out + ".pdf")
    canv.Print(args.out + ".png")
    txt = "\n".join(summary)
    print("\n" + txt)
    # with open(args.out + "_summary.txt", "w") as fo:
    #     fo.write(txt + "\n")
    # print("Wrote %s.pdf, %s.png, %s_summary.txt" % (args.out, args.out, args.out))


if __name__ == "__main__":
    main()