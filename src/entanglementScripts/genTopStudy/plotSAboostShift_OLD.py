"""
plotSpinAnalyzerBoostShift.py

How much does the DIRECTION of each spin analyser change when the boosts in gen_observables() (makeAnalysisSkims.py)
use the first-copy tops (TTbarGen Top()/Antitop(), gen_ branches) instead of the last-copy tops (W + b, genLC_ branches)?

The two boost steps of gen_observables() are reproduced exactly, once per top definition:
  step 1 (CM):  boost lepton and hadronic b (and both tops) by -(t + tbar).BoostVector()        -> ttbar CM frame
  step 2 (RF):  boost the CM-frame lepton by -t_lep(CM).BoostVector() and the CM-frame b_had by -t_had(CM).BoostVector()
The lepton and b_had lab four-vectors are identical in both definitions; only the tops (and therefore both boosts) differ.

For the SAME particle it then compares the direction obtained with first-copy and last-copy boosts:
  * 3D opening angle between the two unit vectors (what enters cos(theta) = SA . axis)      <- headline quantity
  * DeltaR(eta, phi) in that frame (frame-dependent, shown for comparison with the lab-frame DeltaR plots)
plus, as a reference, the change of the helicity axis k = top direction in the CM frame.

Leptonic/hadronic top assignment: per event from the skim itself. The last-copy top is exactly W + b, so
m(t_LC - b_had) equals the W mass for the hadronic top only; the side closer to 80.4 GeV is taken as hadronic.
The agreement with the gen_Lepton_charge branch (either sign convention) is printed as a cross-check.

Inputs: SC_SkimTree from makeAnalysisSkims.py with DO_LASTCOPY = True (gen_/genLC_ Top, Antitop; gen_ Lepton, bHad).

Usage (CMSSW_10_6_28 / python 2.7 / ROOT 6.14, numpy):
  python plotSpinAnalyzerBoostShift.py TTToSemiLeptonicUL18_SCskim_*.root --out SpinAnalyzerBoostShift
  python plotSpinAnalyzerBoostShift.py skims*.root --cut "reco_is_toptag_reconstruction>0.5" --out SAshift_merged
  python plotSpinAnalyzerBoostShift.py skims*.root --max-entries 500000          # quick look
Outputs: <out>.pdf (page 1: 3D angle, page 2: DeltaR), <out>_angle.png, <out>_deltaR.png, <out>_summary.txt
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
C_CM, C_RF, C_K = "#2a78d6", "#eb6834", "#1baf7a"   # validated categorical slots 1-3; also distinguished by line style + legend

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


def delta_r(a, b):
    def etaphi(v):
        pt = np.hypot(v[:, 1], v[:, 2])
        return np.arcsinh(v[:, 3] / np.where(pt > 0, pt, 1e-12)), np.arctan2(v[:, 2], v[:, 1])
    ea, pa = etaphi(a)
    eb, pb = etaphi(b)
    dphi = np.abs(pa - pb)
    dphi = np.where(dphi > math.pi, 2 * math.pi - dphi, dphi)
    return np.hypot(ea - eb, dphi)


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
QUANTS = [  # key, label, colour, line style, row
    ("lep_cm",  "in t#bar{t} CM frame",          C_CM, 1),
    ("lep_rf",  "in leptonic-Top Rest frame",    C_RF, 2),
    ("bhad_cm", "in t#bar{t} CM frame",          C_CM, 1),
    ("bhad_rf", "in hadronic-Top Rest frame",    C_RF, 2),
    ("k",       "k = top direction in t#bar{t} CM frame",   C_K, 1),
]


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

    # display histograms (100 bins) + fine histograms for quantiles; both for the 3D angle and DeltaR
    H = {}
    for kind in ("ang", "dr"):
        for key, _, _, _ in QUANTS:
            H[(kind, key)] = ROOT.TH1D("h_%s_%s" % (kind, key), "", 100, 0.0, args.xmax)
            H[(kind, key, "fine")] = ROOT.TH1D("f_%s_%s" % (kind, key), "", 20000, 0.0, math.pi * (2 if kind == "dr" else 1))
            H[(kind, key)].SetDirectory(0); H[(kind, key, "fine")].SetDirectory(0)

    n_sel, n_agree_pos, printed = 0, 0, False
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
            n_agree_pos += int(np.sum(top_is_lep == (qbranch < 0)))     # convention: Lepton_charge = -1 for l+

            fc = gen_observables_frames(top_fc, atop_fc, lep, bhad, top_is_lep)
            lc = gen_observables_frames(top_lc, atop_lc, lep, bhad, top_is_lep)
            for key, _, _, _ in QUANTS:
                a_ang = angle3d(fc[key], lc[key])
                a_dr = delta_r(fc[key], lc[key])
                fill(H[("ang", key)], a_ang, w); fill(H[("ang", key, "fine")], a_ang, w)
                fill(H[("dr", key)], a_dr, w); fill(H[("dr", key, "fine")], a_dr, w)

            if not printed:
                printed = True
                mw_right = np.where(top_is_lep, m_atop_minus_b, m_top_minus_b)
                print("Self-check on the first chunk of %s (%d events):" % (fn, len(w)))
                print("  m(t_had,LC - b_had): median %.2f GeV, %.1f%% within 80.4 +- 15 GeV (must be ~W line shape)"
                      % (np.median(mw_right), 100.0 * np.mean(np.abs(mw_right - MW) < 15)))
                for i in range(min(3, len(w))):
                    print("  event %d: top_is_leptonic=%s  angle(lep) CM %.4f RF %.4f   angle(b_had) CM %.4f RF %.4f rad"
                          % (i, bool(top_is_lep[i]), angle3d(fc["lep_cm"][i:i+1], lc["lep_cm"][i:i+1])[0],
                             angle3d(fc["lep_rf"][i:i+1], lc["lep_rf"][i:i+1])[0],
                             angle3d(fc["bhad_cm"][i:i+1], lc["bhad_cm"][i:i+1])[0],
                             angle3d(fc["bhad_rf"][i:i+1], lc["bhad_rf"][i:i+1])[0]))
            n_sel += len(w)
        f.Close()
        print("  %s done (%d selected events so far)" % (fn, n_sel))

    agree = n_agree_pos / float(max(n_sel, 1))
    print("Leptonic-top assignment (W-mass test) agrees with gen_Lepton_charge(-1 = l+) in %.2f%% of events" % (100 * agree))

    # ------------------------------------------------------------------ plots
    canv = ROOT.TCanvas("c", "c", 1600, 620)
    pdf = args.out + ".pdf"
    canv.Print(pdf + "[")
    summary = ["Selection: %s   Weight: %s   Events: %d" % (sel, args.weight, n_sel),
               "Leptonic-top assignment from m(t_LC - b_had) = m_W; agreement with gen_Lepton_charge (-1 = l+): %.2f%%" % (100 * agree), ""]
    keep = []
    for kind, xt, tag in (("ang", "3D angle between FC and LC boosted directions [rad]", "angle"),
                          ("dr", "#DeltaR(#eta,#phi) between FC and LC boosted directions", "deltaR")):
        canv.Clear(); canv.Divide(3, 1)
        panels = [("Lepton (spin analyser 1)", ["lep_cm", "lep_rf"]),
                  ("Hadronic b quark (spin analyser 2)", ["bhad_cm", "bhad_rf"]),
                  ("Reference: helicity axis k", ["k"])]
        summary.append("== %s ==" % ("3D angle" if kind == "ang" else "DeltaR(eta,phi)"))
        for ip, (ptitle, keys) in enumerate(panels):
            pad = canv.cd(ip + 1)
            pad.SetLeftMargin(0.14); pad.SetRightMargin(0.03); pad.SetTopMargin(0.12); pad.SetBottomMargin(0.13)
            pad.SetGridy(1)
            hs = []
            for key in keys:
                h = H[(kind, key)].Clone("d_%s_%s" % (kind, key))
                tot = H[(kind, key, "fine")].Integral(0, H[(kind, key, "fine")].GetNbinsX() + 1)
                if tot > 0:
                    h.Scale(1.0 / tot)
                q = [k for k in QUANTS if k[0] == key][0]
                color(h, q[2], q[3])
                hs.append((key, h, q))
            ymax = 1.3 * max(h.GetMaximum() for _, h, _ in hs)
            for i, (key, h, q) in enumerate(hs):
                h.SetMaximum(ymax); h.SetMinimum(0)
                h.GetXaxis().SetTitle(xt); h.GetYaxis().SetTitle("fraction of events / bin")
                h.GetXaxis().SetTitleSize(0.040); h.GetYaxis().SetTitleSize(0.045)
                h.GetXaxis().SetLabelSize(0.040); h.GetYaxis().SetLabelSize(0.040); h.GetYaxis().SetTitleOffset(1.5)
                h.Draw("hist" if i == 0 else "hist same")
            leg = ROOT.TLegend(0.20, 0.60, 0.96, 0.86)
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.034)
            for key, h, q in hs:
                hf = H[(kind, key, "fine")]
                med, q90, mu = quantile(hf, 0.5), quantile(hf, 0.9), mean(hf)
                leg.AddEntry(h, q[1], "l")
                blank = ROOT.TH1D("blank_%s_%s" % (kind, key), "", 1, 0, 1); blank.SetLineColor(0); keep.append(blank)
                leg.AddEntry(blank, "   mean %s, median %s, 90%% < %s" % (mrad(mu), mrad(med), mrad(q90)), "")
                summary.append("  %-9s %-45s mean %-11s median %-11s 90%% quantile %s"
                               % (key, q[1].replace("#bar{t}", "tbar").replace("t#", "t"), mrad(mu), mrad(med), mrad(q90)))
            leg.Draw()
            tl = ROOT.TLatex(); tl.SetNDC(); tl.SetTextFont(62); tl.SetTextSize(0.045)
            tl.DrawLatex(0.14, 0.92, ptitle)
            keep.extend([h for _, h, _ in hs] + [leg, tl])
        canv.cd(0)
        hdr = ROOT.TLatex(); hdr.SetNDC(); hdr.SetTextFont(42); hdr.SetTextSize(0.026)
        hdr.DrawLatex(0.005, 0.965, "Same particle, boosted with first-copy (TTbarGen Top()) vs last-copy (W + b) tops, "
                      "following gen_observables():  %d events,  selection: %s" % (n_sel, args.cut))
        keep.append(hdr)
        canv.Print(pdf)
        canv.Print("%s_%s.png" % (args.out, tag))
        summary.append("")
    canv.Print(pdf + "]")
    txt = "\n".join(summary)
    print("\n" + txt)
    with open(args.out + "_summary.txt", "w") as fo:
        fo.write(txt + "\n")
    print("Wrote %s, %s_angle.png, %s_deltaR.png, %s_summary.txt" % (pdf, args.out, args.out, args.out))


if __name__ == "__main__":
    main()