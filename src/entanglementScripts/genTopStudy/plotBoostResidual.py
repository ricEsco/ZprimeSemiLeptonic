"""
plotBoostResidual.py

How good is the rest frame given by each gen-level top definition?
Boost two complementary decay products A and B of a top into the rest frame of that top and measure the residual
angle  delta = pi - angle(A, B)  between them. For a correct parent frame A and B are exactly back-to-back (delta = 0).

  FC : parent = gen_<parent>    (TTbarGen Top()/Antitop() = first copy, status 22)
  LC : parent = genLC_<parent>  (W + b = last copy, status 62 = statusFlags().isLastCopy())

Presets (--presets, comma-separated; one page each):
  Wb_lep : W_lep  vs b_lep             in the leptonic-top rest frame
  Wb_had : W_had  vs b_had             in the hadronic-top rest frame
  lep    : lepton vs (neutrino + b_lep) in the leptonic-top rest frame   <- lepton spin analyser
  bhad   : b_had  vs (q1 + q2)          in the hadronic-top rest frame   <- hadronic-b spin analyser
Custom: --parent TopLep --A Lepton --B Neutrino+bLep   (names = the gen_<name>_{pt,eta,phi,E} branches of the skim)

Notes
  * A and B must be the two sides of ONE parent's decay.
  * Wb_*: the last-copy top IS W + b, so LC gives delta = 0 by construction (up to float precision of the stored
    four-vectors); the FC curve shows how far the TTbarGen Top() frame is from the frame the W and b were decayed in.
  * lep / bhad: l + nu and q1 + q2 sum to the LAST copy of the W, which differs from TTbarGen's W (status 22) when the W
    absorbed recoil from radiation off the b (~70% of events). LC is then small but not exactly zero; this is a genuine,
    independent check of the frame.

Needs skims from makeAnalysisSkims.py with DO_LASTCOPY = True and the TopLep/TopHad/WLep/WHad/bLep/Neutrino/Q1/Q2 branches.

Usage (CMSSW_10_6_28 / python 2.7 / ROOT 6.14 or newer):
  python plotBoostResidual.py TTToSemiLeptonicUL18_SCskim_*.root --out BoostResidual
  python plotBoostResidual.py skims*.root --presets Wb_lep,Wb_had --cut "reco_is_toptag_reconstruction>0.5"
  python plotBoostResidual.py skims*.root --presets custom --parent TopLep --A Lepton --B Neutrino+bLep
Outputs: <out>.pdf (one page per preset), <out>_<preset>.png, <out>_summary.txt
"""
from __future__ import print_function, division

import argparse
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)

PRESETS = {
    "Wb_lep": dict(parent="TopLep", A="WLep", B="bLep",
                   title="W_{lep} vs b_{lep} in the leptonic-top rest frame"),
    "Wb_had": dict(parent="TopHad", A="WHad", B="bHad",
                   title="W_{had} vs b_{had} in the hadronic-top rest frame"),
    "lep":    dict(parent="TopLep", A="Lepton", B="Neutrino+bLep",
                   title="lepton vs (#nu + b_{lep}) in the leptonic-top rest frame"),
    "bhad":   dict(parent="TopHad", A="bHad", B="Q1+Q2",
                   title="b_{had} vs (q_{1} + q_{2}) in the hadronic-top rest frame"),
}

import numpy as np


def read_columns(tree, exprs, sel, nmax):
    """Values of the TTree expressions for the selected entries, via TTree::Draw (4 columns per call, same order)."""
    n_all = tree.GetEntries()
    tree.SetEstimate(n_all + 1)
    out = []
    for k in range(0, len(exprs), 4):
        chunk = exprs[k:k + 4]
        if nmax >= 0:
            n = tree.Draw(":".join(chunk), sel, "goff", nmax)
        else:
            n = tree.Draw(":".join(chunk), sel, "goff")
        if n < 0:
            raise RuntimeError("TTree::Draw failed for %s (check branch names)" % chunk)
        for j in range(len(chunk)):
            buf = getattr(tree, "GetV%d" % (j + 1))()
            if n == 0:
                out.append(np.zeros(0))
                continue
            buf.SetSize(n)
            out.append(np.array(buf, dtype=np.float64, copy=True))
    return out


def p4_from(pt, eta, phi, E):
    """(N, 4) array of (E, px, py, pz) from pt, eta, phi, E arrays."""
    return np.stack([E, pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)], axis=1)


def boost_to_rest(v, parent):
    """Boost four-vectors v (N,4) into the rest frame of parent (N,4); same formula as TLorentzVector::Boost(-beta)."""
    b = -parent[:, 1:] / parent[:, :1]
    b2 = np.sum(b * b, axis=1)
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(b * v[:, 1:], axis=1)
    g2 = np.where(b2 > 0, (g - 1.0) / np.where(b2 > 0, b2, 1.0), 0.0)
    E = g * (v[:, 0] + bp)
    p = v[:, 1:] + ((g2 * bp + g * v[:, 0])[:, None]) * b
    return np.column_stack([E, p])


def residual_angle(parent, a, b):
    """delta = angle between A and -B in the parent rest frame (0 = back-to-back); atan2 form, precise for tiny angles."""
    va = boost_to_rest(a, parent)[:, 1:]
    vc = -boost_to_rest(b, parent)[:, 1:]
    return np.arctan2(np.linalg.norm(np.cross(va, vc), axis=1), np.sum(va * vc, axis=1))


def mass(v):
    m2 = v[:, 0] ** 2 - np.sum(v[:, 1:] ** 2, axis=1)
    return np.sign(m2) * np.sqrt(np.abs(m2))


def fill_hist(h, x, w):
    from array import array
    if len(x):
        h.FillN(len(x), array("d", x.tolist()), array("d", w.tolist()))


C_FC, C_LC = "#eb6834", "#2a78d6"     # validated categorical pair (orange / blue); identity also via legend + line style


def log_edges(lo_exp=-7, hi=math.pi, per_decade=20):
    import array
    edges, k = [], 0
    while True:
        e = 10 ** (lo_exp + k / float(per_decade))
        if e >= hi:
            edges.append(hi + 1e-9)
            break
        edges.append(e)
        k += 1
    return array.array("d", edges)


def frac_below(h, x):
    """weighted fraction of entries (incl. underflow, excl. nothing) with value < x (x should be a bin edge)."""
    tot = h.Integral(0, h.GetNbinsX() + 1)
    b = h.GetXaxis().FindBin(x * (1 + 1e-9)) - 1
    return h.Integral(0, b) / tot if tot > 0 else float("nan")


def quantile(h, q):
    """quantile from the (fine, log-binned) histogram incl. underflow; returns 0 if it falls in the underflow."""
    tot = h.Integral(0, h.GetNbinsX() + 1)
    if tot <= 0:
        return float("nan")
    if h.GetBinContent(0) / tot >= q:
        return 0.0
    acc = h.GetBinContent(0)
    for i in range(1, h.GetNbinsX() + 1):
        c = h.GetBinContent(i)
        if acc + c >= q * tot:
            lo, hi = h.GetXaxis().GetBinLowEdge(i), h.GetXaxis().GetBinUpEdge(i)
            f = (q * tot - acc) / c if c > 0 else 0.0
            return math.exp(math.log(lo) + f * (math.log(hi) - math.log(lo)))
        acc += c
    return h.GetXaxis().GetXmax()


def style(h, color, lstyle):
    c = ROOT.TColor.GetColor(color)
    h.SetLineColor(c); h.SetMarkerColor(c); h.SetLineWidth(3); h.SetLineStyle(lstyle); h.SetFillStyle(0)


def cumulative(h, name):
    """fraction of events with residual < x (underflow included)."""
    c = h.Clone(name)
    tot = h.Integral(0, h.GetNbinsX() + 1)
    acc = h.GetBinContent(0)
    for i in range(1, h.GetNbinsX() + 1):
        acc += h.GetBinContent(i)
        c.SetBinContent(i, acc / tot if tot > 0 else 0.0)
        c.SetBinError(i, 0.0)
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--tree", default="SC_SkimTree")
    ap.add_argument("--presets", default="Wb_lep,Wb_had,lep,bhad", help="comma-separated: %s, custom" % ", ".join(PRESETS))
    ap.add_argument("--parent", default="TopLep", help="custom: parent name (TopLep, TopHad, Top, Antitop)")
    ap.add_argument("--A", default="Lepton", help="custom: decay product A, sum with '+', e.g. Lepton")
    ap.add_argument("--B", default="Neutrino+bLep", help="custom: decay product B, sum with '+', e.g. Neutrino+bLep")
    ap.add_argument("--cut", default="1", help='extra selection, e.g. "reco_is_toptag_reconstruction>0.5"')
    ap.add_argument("--weight", default="genweight")
    ap.add_argument("--max-entries", type=int, default=-1, help="entries read per file; -1 = all")
    ap.add_argument("--out", default="BoostResidual")
    args = ap.parse_args()

    edges = log_edges()
    cut = "(gen_semilep == 1) && (%s)" % args.cut

    configs = []
    for p in [x.strip() for x in args.presets.split(",") if x.strip()]:
        if p == "custom":
            configs.append(("custom", dict(parent=args.parent, A=args.A, B=args.B,
                                           title="%s vs %s in the %s rest frame" % (args.A, args.B, args.parent))))
        else:
            configs.append((p, PRESETS[p]))

    canv = ROOT.TCanvas("c", "c", 1500, 650)
    pdf = args.out + ".pdf"
    canv.Print(pdf + "[")
    summary = ["Selection: %s   Weight: %s   Files: %d" % (cut, args.weight, len(args.files)), ""]
    keep = []
    for key, cfg in configs:
        hFC = ROOT.TH1D("hFC_" + key, "", len(edges) - 1, edges)
        hLC = ROOT.TH1D("hLC_" + key, "", len(edges) - 1, edges)
        names_A = [n.strip() for n in cfg["A"].split("+")]
        names_B = [n.strip() for n in cfg["B"].split("+")]
        prefixes = (["gen_" + n for n in names_A] + ["gen_" + n for n in names_B]
                    + ["gen_" + cfg["parent"], "genLC_" + cfg["parent"]])
        exprs = ["%s_%s" % (pre, c) for pre in prefixes for c in ("pt", "eta", "phi", "E")] + [args.weight]
        nsel = 0
        for fn in args.files:
            f = ROOT.TFile.Open(fn)
            t = f.Get(args.tree)
            cols = read_columns(t, exprs, cut, args.max_entries)
            f.Close()
            vecs = [p4_from(*cols[4 * i:4 * i + 4]) for i in range(len(prefixes))]
            w = cols[-1]
            a = sum(vecs[:len(names_A)])
            b = sum(vecs[len(names_A):len(names_A) + len(names_B)])
            pFC, pLC = vecs[-2], vecs[-1]
            rFC, rLC = residual_angle(pFC, a, b), residual_angle(pLC, a, b)
            fill_hist(hFC, rFC, w)
            fill_hist(hLC, rLC, w)
            if nsel == 0 and len(w):
                # self-check printout: if these look wrong, the skim content (not the boost) is the problem
                print("[%s] first selected events of %s:" % (key, fn))
                for i in range(min(3, len(w))):
                    print("   A=(%.1f,%.1f,%.1f,%.1f) B=(%.1f,%.1f,%.1f,%.1f) m(A+B)=%.2f  m(FC)=%.2f  m(LC)=%.2f"
                          "   delta FC=%.3e  LC=%.3e  w=%.3g"
                          % (tuple(a[i]) + tuple(b[i]) + (mass((a + b)[i:i + 1])[0], mass(pFC[i:i + 1])[0],
                                                           mass(pLC[i:i + 1])[0], rFC[i], rLC[i], w[i])))
                print("   median |m(A+B) - m(parent)| over this file:  FC %.3g GeV,  LC %.3g GeV  "
                      "(for W+b presets LC must be ~0)" % (np.median(np.abs(mass(a + b) - mass(pFC))),
                                                           np.median(np.abs(mass(a + b) - mass(pLC)))))
            nsel += len(w)
        print("[%s] %d selected events" % (key, nsel))

        stats = {}
        for tag, h in (("FC", hFC), ("LC", hLC)):
            stats[tag] = dict(med=quantile(h, 0.5), q90=quantile(h, 0.9),
                              f1m=frac_below(h, 1e-3), f10m=frac_below(h, 1e-2), f100m=frac_below(h, 1e-1),
                              und=h.GetBinContent(0) / max(h.Integral(0, h.GetNbinsX() + 1), 1e-300))

        # ---------------- page ----------------
        canv.Clear()
        canv.Divide(2, 1)
        title = ROOT.TLatex()
        title.SetNDC(); title.SetTextFont(62); title.SetTextSize(0.034)
        for ipad, (kind) in enumerate(("dist", "cdf")):
            pad = canv.cd(ipad + 1)
            pad.SetLogx(1); pad.SetTopMargin(0.16); pad.SetLeftMargin(0.12); pad.SetRightMargin(0.04)
            pad.SetBottomMargin(0.13); pad.SetGridy(1)
            if kind == "dist":
                a = hFC.Clone("dFC_" + key); b = hLC.Clone("dLC_" + key)
                for h in (a, b):
                    tot = h.Integral(0, h.GetNbinsX() + 1)
                    if tot > 0:
                        h.Scale(1.0 / tot)
                ytitle = "fraction of events / bin"
            else:
                a = cumulative(hFC, "cFC_" + key); b = cumulative(hLC, "cLC_" + key)
                ytitle = "fraction of events with #delta < x"
            style(a, C_FC, 1); style(b, C_LC, 2)
            a.GetXaxis().SetTitle("residual angle #delta = #pi - #angle(A, B) in the parent rest frame [rad]")
            a.GetYaxis().SetTitle(ytitle)
            for ax in (a.GetXaxis(), a.GetYaxis()):
                ax.SetTitleSize(0.042); ax.SetLabelSize(0.038)
            a.GetYaxis().SetTitleOffset(1.35); a.GetXaxis().SetTitleOffset(1.25)
            a.GetXaxis().SetMoreLogLabels(False)
            if kind == "dist":
                a.SetMaximum(1.25 * max(a.GetMaximum(), b.GetMaximum())); a.SetMinimum(0)
            else:
                a.SetMaximum(1.05); a.SetMinimum(0)
            a.Draw("hist"); b.Draw("hist same")
            leg = ROOT.TLegend(0.14, 0.62, 0.62, 0.82) if kind == "dist" else ROOT.TLegend(0.14, 0.66, 0.60, 0.82)
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.031)
            blank = ROOT.TH1D("blank_%s_%s" % (key, kind), "", 1, 0, 1)   # text-only legend rows
            blank.SetLineColor(0); keep.append(blank)
            for tag, h, name in (("FC", a, "first copy (TTbarGen Top())"), ("LC", b, "last copy (W + b)")):
                s = stats[tag]
                if kind == "dist":
                    leg.AddEntry(h, "%s:  median %s" % (name, fmt_angle(s["med"])), "l")
                    leg.AddEntry(blank, "    #delta < 1 mrad: %.1f%%,  < 10 mrad: %.1f%%" % (100 * s["f1m"], 100 * s["f10m"]), "")
                else:
                    leg.AddEntry(h, "%s:  90%% below %s" % (name, fmt_angle(s["q90"])), "l")
            leg.Draw()
            keep.extend([a, b, leg])
        canv.cd(0)
        title.DrawLatex(0.01, 0.955, "Back-to-back test: %s" % cfg["title"])
        sub = ROOT.TLatex(); sub.SetNDC(); sub.SetTextFont(42); sub.SetTextSize(0.024); sub.SetTextColor(ROOT.kGray + 2)
        # sub.DrawLatex(0.01, 0.915, "A = %s,  B = %s,  parent = gen_%s (first copy) / genLC_%s (last copy).   "
        #               "A correct parent frame gives #delta = 0.   %d events,  selection: %s"
        #               % (cfg["A"], cfg["B"], cfg["parent"], cfg["parent"], nsel, args.cut))
        keep.extend([title, sub])
        canv.Print(pdf)
        canv.Print("%s_%s.png" % (args.out, key))

        summary.append("[%s] %s  (A = %s, B = %s, parent = %s, N = %d)" % (key, cfg["title"].replace("#", ""), cfg["A"], cfg["B"],
                                                                              cfg["parent"], nsel))
        for tag in ("FC", "LC"):
            s = stats[tag]
            summary.append("   %s: median %-12s 90%% quantile %-12s  <1 mrad %6.2f%%  <10 mrad %6.2f%%  <100 mrad %6.2f%%  "
                           "(<1e-7 rad: %.2f%%)" % (tag, fmt_angle(s["med"], root=False), fmt_angle(s["q90"], root=False),
                                                     100 * s["f1m"], 100 * s["f10m"], 100 * s["f100m"], 100 * s["und"]))
        summary.append("")
    canv.Print(pdf + "]")
    txt = "\n".join(summary)
    print("\n" + txt)
    with open(args.out + "_summary.txt", "w") as f:
        f.write(txt + "\n")
    print("Wrote %s, %s_<preset>.png, %s_summary.txt" % (pdf, args.out, args.out))


def fmt_angle(x, root=True):
    if x != x:
        return "n/a"
    if x == 0.0:
        return "< 10^{-7} rad" if root else "<1e-7 rad"
    if x < 1e-3:
        return ("%.1f #murad" if root else "%.1f urad") % (x * 1e6)
    if x < 1.0:
        return "%.1f mrad" % (x * 1e3)
    return "%.2f rad" % x


if __name__ == "__main__":
    main()