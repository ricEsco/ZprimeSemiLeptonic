"""
plotTopDefinitionComparison.py

Compares the gen-level spin-correlation observables computed with the two top definitions stored in the
SC_SkimTree produced by makeAnalysisSkims.py (DO_LASTCOPY = True):

  gen_<Obs>    : TTbarGen Top()/Antitop()  = FIRST-copy (status-22) tops   ("FC")
  genLC_<Obs>  : W + b / W- + bbar         = LAST-copy (status-62) tops    ("LC", == statusFlags().isLastCopy())

Lepton and hadronic b are the same in both; only the top four-vectors (ttbar frame, helicity basis, rest-frame boosts) differ.

Per observable it produces
  - normalised distributions FC vs LC with an LC/FC ratio panel (each definition with its own valid-event population,
    i.e. exactly what an analysis using that definition would see; matters for the windowed slices)
  - the per-event difference LC - FC (events valid in both)
  - the 2D migration LC vs FC (events valid in both)
and a summary table of weighted means (mean of the per-event projections/products is what enters the spin-correlation
coefficients), their shift LC - FC and its statistical significance.
Also plots the per-event first-vs-last-copy top diagnostics (dR, relative pT, mass difference).

Usage (CMSSW_10_6_28 / python 2.7 / ROOT 6.14, or any newer ROOT):
  python plotTopDefinitionComparison.py TTToSemiLeptonicUL18_SCskim_*.root \
      [--weight genweight] [--cut "reco_is_toptag_reconstruction>0.5"] [--out TopDefComparison] \
      [--vars cHel,Ckk,Crr,Cnn] [--max-entries 1000000]

Outputs: <out>.pdf (all plots), <out>.root (all histograms), <out>_summary.txt and <out>_summary.csv
"""
from __future__ import print_function, division

import argparse
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.TH1.StatOverflows(True)          # means/RMS include under/overflow -> exact, independent of the plotting range

SENT, NOTp = -10.0, 99.0
PI = math.pi

# Observables stored with both gen_ and genLC_ prefixes
OBSERVABLES = [
    "M_tt", "beta", "absDeltaY", "cosThetaStar", "pt_hadTop",
    "cosTheta1k_antiLep", "cosTheta1r_antiLep", "cosTheta1n_antiLep", "cosTheta1kStar_antiLep", "cosTheta1rStar_antiLep",
    "cosTheta2k_Lep", "cosTheta2r_Lep", "cosTheta2n_Lep", "cosTheta2kStar_Lep", "cosTheta2rStar_Lep",
    "cosTheta1k", "cosTheta1r", "cosTheta1n", "cosTheta1kStar", "cosTheta1rStar",
    "cosTheta2k", "cosTheta2r", "cosTheta2n", "cosTheta2kStar", "cosTheta2rStar",
    "Cnn", "Cnr", "Cnk", "Crn", "Crr", "Crk", "Ckn", "Ckr", "Ckk",
    "Crk_plus", "Crk_minus", "Cnr_plus", "Cnr_minus", "Cnk_plus", "Cnk_minus",
    "cHel", "cHel_Mtt300_400", "cHel_Mtt300_400_betaLT0p9",
    "cHel_P3n", "cHel_P3n_Mtt800_Inf", "cHel_P3n_Mtt800_Inf_cosThetaLT0p4",
    "Sigma_phi", "Sigma_phi_low", "Sigma_phi_high",
    "Delta_phi", "Delta_phi_low", "Delta_phi_high",
]

DIAGNOSTICS = [  # (branch, title, bins)
    ("gen_topFCvsLC_dR",        "#DeltaR(t_{FC}, t_{LC})",                 (100, 0.0, 1.0)),
    ("gen_antitopFCvsLC_dR",    "#DeltaR(#bar{t}_{FC}, #bar{t}_{LC})",     (100, 0.0, 1.0)),
    ("gen_topFCvsLC_relPt",     "(p_{T}^{LC} - p_{T}^{FC}) / p_{T}^{FC}  (top)",      (100, -1.0, 1.0)),
    ("gen_antitopFCvsLC_relPt", "(p_{T}^{LC} - p_{T}^{FC}) / p_{T}^{FC}  (antitop)",  (100, -1.0, 1.0)),
    ("gen_topFCvsLC_dM",        "m_{LC} - m_{FC} [GeV] top (expect 0)",     (100, -0.5, 0.5)),
    ("gen_antitopFCvsLC_dM",    "m_{LC} - m_{FC} [GeV] antitop (expect 0)", (100, -0.5, 0.5)),
]


def binning(var):
    """(nbins, lo, hi) for the distribution and for the per-event difference."""
    if var.startswith("Sigma_phi") or var.startswith("Delta_phi"):
        return (40, -PI, PI), (100, -2 * PI, 2 * PI)
    if var == "M_tt":
        return (60, 250.0, 1500.0), (100, -100.0, 100.0)
    if var == "beta":
        return (50, 0.0, 1.0), (100, -0.3, 0.3)
    if var == "absDeltaY":
        return (60, -3.0, 3.0), (100, -1.0, 1.0)
    if var == "pt_hadTop":
        return (60, 0.0, 900.0), (100, -200.0, 200.0)
    # cosines, products of cosines, cHel variants: all in [-1, 1]
    return (40, -1.0, 1.0), (100, -2.0, 2.0)


def valid(branch):
    return "(abs(%s - (%g)) > 1e-4 && abs(%s - (%g)) > 1e-4)" % (branch, SENT, branch, NOTp)


def draw(chain, expr, hist, weight, cut, nmax):
    ROOT.gROOT.cd()
    hist.SetDirectory(ROOT.gROOT)
    sel = "(%s)*(%s)" % (weight, cut)
    if nmax > 0:
        chain.Draw("%s>>%s" % (expr, hist.GetName()), sel, "goff", nmax)
    else:
        chain.Draw("%s>>%s" % (expr, hist.GetName()), sel, "goff")
    return hist


def neff(h):
    return h.GetEffectiveEntries()


def style(h, color, marker=20):
    h.SetLineColor(color); h.SetMarkerColor(color); h.SetMarkerStyle(marker); h.SetMarkerSize(0.6); h.SetLineWidth(1)


def latex(x, y, text, size=0.035):
    l = ROOT.TLatex()
    l.SetNDC(); l.SetTextSize(size); l.SetTextFont(42)
    l.DrawLatex(x, y, text)
    return l


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--tree", default="SC_SkimTree")
    ap.add_argument("--weight", default="genweight", help='event weight expression, e.g. "genweight", "1", "eventweight"')
    ap.add_argument("--cut", default="1", help='extra selection, e.g. "reco_is_toptag_reconstruction>0.5"')
    ap.add_argument("--out", default="TopDefComparison")
    ap.add_argument("--vars", default="", help="comma-separated subset of observables (default: all)")
    ap.add_argument("--max-entries", type=int, default=-1)
    args = ap.parse_args()

    chain = ROOT.TChain(args.tree)
    for f in args.files:
        chain.Add(f)
    print("Entries in chain: %d" % chain.GetEntries())

    branches = set(b.GetName() for b in chain.GetListOfBranches())
    obs = [v.strip() for v in args.vars.split(",") if v.strip()] if args.vars else OBSERVABLES
    obs = [v for v in obs if ("gen_" + v) in branches and ("genLC_" + v) in branches]
    if not obs:
        raise SystemExit("No gen_/genLC_ observable pairs found - were the skims produced with DO_LASTCOPY = True?")

    base = "(gen_semilep == 1) && (%s)" % args.cut
    outroot = ROOT.TFile(args.out + ".root", "RECREATE")
    pdf = args.out + ".pdf"
    canv = ROOT.TCanvas("c", "c", 1600, 800)
    canv.Print(pdf + "[")

    # ---------------- title page ----------------
    canv.Clear()
    latex(0.05, 0.85, "#bf{Gen-level observables: first_copy (FC) vs last_copy (LC) top definition}", 0.060)
    latex(0.05, 0.72, "FC Top = TTbarGen   Top()    = status 22 is from hard scatter (before ISR recoil / primordial kT)", 0.04)
    latex(0.05, 0.64, "LC Top = TTbarGen W() + b()  = status 62 is parent of TTbarGen's W&b and is the last copy before decay", 0.04)
    latex(0.05, 0.52, "Selection: %s" % base.replace("&&", " && "), 0.035)
    latex(0.05, 0.45, "Entries: %d" % chain.GetEntries(), 0.035)
    latex(0.05, 0.38, "Weight: %s " % args.weight, 0.035)
    canv.Print(pdf)

    # ---------------- diagnostics ----------------
    canv.Clear()
    canv.SetCanvasSize(600, 600)
    canv.Divide(2, 3)
    keep = []
    for i, (br, title, b) in enumerate(DIAGNOSTICS):
        if br not in branches:
            continue
        canv.cd(i + 1)
        ROOT.gPad.SetTopMargin(0.05); ROOT.gPad.SetBottomMargin(0.2)
        h = ROOT.TH1D("diag_" + br, ";%s;events (weighted)" % title, *b)
        draw(chain, br, h, args.weight, "%s && %s" % (base, valid(br)), args.max_entries)
        style(h, ROOT.kAzure + 2)
        h.GetXaxis().SetLabelSize(0.04); h.GetXaxis().SetTitleSize(0.07); h.GetXaxis().SetTitleOffset(1.0)
        h.GetYaxis().SetLabelSize(0.04); h.GetYaxis().SetTitleSize(0.05); h.GetYaxis().SetTitleOffset(0.8)
        h.Draw("hist")
        latex(0.25, 0.85, "mean %.4f   RMS %.4f" % (h.GetMean(), h.GetRMS()), 0.065)
        outroot.cd(); h.Write(); keep.append(h)
    canv.Print(pdf)

    # ---------------- per-observable comparison ----------------
    rows = []
    for var in obs:
        # apply binning and validity cuts on histograms
        (nb, lo, hi), (nbd, lod, hid) = binning(var)
        fc_b, lc_b = "gen_" + var, "genLC_" + var
        cut_fc = "%s && %s" % (base, valid(fc_b))
        cut_lc = "%s && %s" % (base, valid(lc_b))
        cut_both = "%s && %s && %s" % (base, valid(fc_b), valid(lc_b))
        # draw all histograms
        h_fc = draw(chain, fc_b, ROOT.TH1D("FC_" + var, ";%s;area normalised" % var, nb, lo, hi), args.weight, cut_fc, args.max_entries)
        h_lc = draw(chain, lc_b, ROOT.TH1D("LC_" + var, ";%s;area normalised" % var, nb, lo, hi), args.weight, cut_lc, args.max_entries)
        h_d = draw(chain, "%s - %s" % (lc_b, fc_b),
                   ROOT.TH1D("diff_" + var, ";%s_LC - %s_FC;events (weighted)" % (var, var), nbd, lod, hid),
                   args.weight, cut_both, args.max_entries)
        h_2d = draw(chain, "%s:%s" % (lc_b, fc_b),
                    ROOT.TH2D("mig_" + var, ";%s_FC;%s_LC" % (var, var), nb, lo, hi, nb, lo, hi),
                    args.weight, cut_both, args.max_entries)
        ### ---------------------------------------- ###
        ### Compute relevant statistical quantitites ###

        # normalized distribution: mean and mean error for FC and LC
        m_fc, e_fc = h_fc.GetMean(), h_fc.GetMeanError()
        m_lc, e_lc = h_lc.GetMean(), h_lc.GetMeanError()

        # per-event difference: mean, mean error, RMS
        m_d, e_d, r_d = h_d.GetMean(), h_d.GetMeanError(), h_d.GetRMS()
        # shift in mean between FC and LC normalized distributions
        shift = m_lc - m_fc
        # relative shift relative to FC mean
        rel = shift / m_fc if m_fc != 0 else float("nan")
        # statistical significance of the shift (mean difference divided by its error)
        sig = m_d / e_d if e_d > 0 else float("nan")

        # Correlation factor between FC and LC per-event values
        corr = h_2d.GetCorrelationFactor()

        # append to summary table
        rows.append((var, neff(h_fc), m_fc, e_fc, neff(h_lc), m_lc, e_lc, shift, rel, m_d, e_d, sig, r_d, corr))
        # write histograms to out root file
        for h in (h_fc, h_lc, h_d, h_2d):
            outroot.cd(); h.Write()

        # ---- page ----
        canv.Clear()
        canv.SetCanvasSize(2400, 600); canv.Divide(3, 1)

        # Normalized distributions FC vs LC with ratio panel, canvas 1 ----------------------------------------------
        p1 = canv.cd(1)
        up = ROOT.TPad("up_" + var, "", 0, 0.30, 1, 1); up.SetBottomMargin(0.02); up.SetTopMargin(0.10)   ; up.SetLeftMargin(0.15); up.Draw()
        dn = ROOT.TPad("dn_" + var, "", 0, 0, 1, 0.30); dn.SetTopMargin(0.02);    dn.SetBottomMargin(0.4) ; dn.SetLeftMargin(0.15); dn.Draw()

        up.cd()
        a, b_ = h_fc.Clone("FCn_" + var), h_lc.Clone("LCn_" + var)
        for h in (a, b_):
            if h.Integral(0, h.GetNbinsX() + 1) > 0:
                h.Scale(1.0 / h.Integral(0, h.GetNbinsX() + 1))
        style(a, ROOT.kBlack, 20); style(b_, ROOT.kRed + 1, 24)
        a.SetMaximum(1.35 * max(a.GetMaximum(), b_.GetMaximum())); a.SetMinimum(0)
        a.GetXaxis().SetLabelSize(0)
        a.GetYaxis().SetLabelSize(0.04); a.GetYaxis().SetTitleSize(0.05); a.GetYaxis().SetTitleOffset(1.1)
        a.Draw("hist"); b_.Draw("hist same")
        leg = ROOT.TLegend(0.15, 0.74, 0.88, 0.89); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04)
        leg.AddEntry(a,  "FC: #LT x #GT = %.4f #pm %.4f" % (m_fc, e_fc), "l")
        leg.AddEntry(b_, "LC: #LT x #GT = %.4f #pm %.4f" % (m_lc, e_lc), "l")
        leg.Draw()

        dn.cd()
        r = b_.Clone("ratio_" + var); r.Divide(a)
        r.SetTitle(";%s;LC / FC" % var); r.SetMinimum(0.8); r.SetMaximum(1.2)
        style(r, ROOT.kRed + 1, 20)
        r.GetXaxis().SetLabelSize(0.08); r.GetXaxis().SetTitleSize(0.13); r.GetXaxis().SetTitleOffset(0.8)
        r.GetYaxis().SetLabelSize(0.08); r.GetYaxis().SetTitleSize(0.05); r.GetYaxis().SetTitleOffset(0.9); r.GetYaxis().SetNdivisions(505)
        r.Draw("ep")
        line = ROOT.TLine(lo, 1.0, hi, 1.0); line.SetLineStyle(2); line.Draw()
        # Per-event difference LC - FC, canvas 2 ----------------------------------------------------------------
        canv.cd(2)
        ROOT.gPad.SetTopMargin(0.083); ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetLogy(1)
        style(h_d, ROOT.kAzure + 2); h_d.Draw("hist")
        h_d.GetXaxis().SetLabelSize(0.04); h_d.GetXaxis().SetTitleSize(0.06); h_d.GetXaxis().SetTitleOffset(0.7)
        h_d.GetYaxis().SetLabelSize(0.04); h_d.GetYaxis().SetTitleSize(0.05); h_d.GetYaxis().SetTitleOffset(1.1)
        latex(0.14, 0.87, "#LT LC - FC #GT = %.4f #pm %.4f  (%.1f#sigma)" % (m_d, e_d, sig), 0.048)
        # latex(0.14, 0.80, "RMS = %.4f" % r_d, 0.045)
        # 2D migration LC vs FC, canvas 3 --------------------------------------------------------------------------
        canv.cd(3)
        ROOT.gPad.SetTopMargin(0.083); ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetLeftMargin(0.10); ROOT.gPad.SetRightMargin(0.13)
        ROOT.gPad.SetLogz(1)
        h_2d.GetXaxis().SetLabelSize(0.04); h_2d.GetXaxis().SetTitleSize(0.06); h_2d.GetXaxis().SetTitleOffset(0.7)
        h_2d.GetYaxis().SetLabelSize(0.04); h_2d.GetYaxis().SetTitleSize(0.05); h_2d.GetYaxis().SetTitleOffset(0.9)
        h_2d.Draw("colz")
        latex(0.14, 0.93, "correlation = %.3f" % corr, 0.06)
        canv.Print(pdf)
        keep.extend([a, b_, r, up, dn, leg, line])

    canv.Print(pdf + "]")
    outroot.Close()

    # ---------------- summary ----------------
    hdr = ("%-34s %10s %18s %10s %18s %10s %9s %20s %7s %8s %6s" %                                                       ("Observable", "Neff(FC)", "<Obs>_FC", "Neff(LC)", "<Obs>_LC", "LC-FC", "rel."   , "<LC-FC>", "signif", "RMS(d)", "corr"))
    lines = [hdr, "-" * len(hdr)]
    for (var, n1, m1, e1, n2, m2, e2, sh, rel, md, ed, sig, rd, corr) in rows:
        lines.append("%-34s %10.0f %9.5f+-%-7.5f %10.0f %9.5f+-%-7.5f %10.5f %8.1f%% %10.5f+-%-8.5f %7.1f %8.4f %6.3f" % (var         , n1        , m1, e1    , n2        , m2, e2    , sh     , 100 * rel, md, ed    , sig      , rd     , corr))
    lines.append("")
    lines.append(" <Obs>: weighted mean")
    lines.append(" <LC-FC> per event uses valid events only")
    lines.append(" 'signif' = <LC-FC>/error.")
    txt = "\n".join(lines)
    print("\n" + txt)
    with open(args.out + "_summary.txt", "w") as f:
        f.write("Selection: %s   Weight: %s\n\n" % (base, args.weight) + txt + "\n")
    with open(args.out + "_summary.csv", "w") as f:
        f.write("observable,neff_FC,mean_FC,err_FC,neff_LC,mean_LC,err_LC,shift,rel_shift,mean_diff,err_diff,signif,rms_diff,corr\n")
        for row in rows:
            f.write(",".join([row[0]] + ["%.6g" % x for x in row[1:]]) + "\n")
    print("\nWrote %s, %s.root, %s_summary.txt, %s_summary.csv" % (pdf, args.out, args.out, args.out))


if __name__ == "__main__":
    main()