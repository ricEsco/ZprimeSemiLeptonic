# check_ntuple_inclusive.py -- is the UHH2 pre-selection ntuple inclusive at gen level, or lepton-skimmed?
# Compare the NORMALIZED gen shapes of the two closure dumps.
#
#   python check_ntuple_inclusive.py  gendump_uhh2.root  gendump_nano.root
#
# Decisive signal: if the ntuple is lepton-skimmed, its gen_mtt is SUPPRESSED near threshold
# (300-400 GeV) relative to NanoAOD's truly-inclusive falling spectrum, and gen cosTheta* is
# sculpted (domed) -- exactly the acceptance effect we already characterized. If inclusive, the
# normalized shapes overlay and the ratio sits at 1.
from __future__ import print_function
import sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)

CUT = "gen_semilep && gen_mtt > 0"          # semileptonic, sentinel-cleaned


def norm(fn, expr, nb, lo, hi, tag):
    f = ROOT.TFile.Open(fn); t = f.Get("gendump")
    t.Draw("%s>>htmp(%d,%g,%g)" % (expr, nb, lo, hi), CUT, "goff")
    h = t.GetHistogram().Clone("h_%s_%s" % (expr, tag)); h.SetDirectory(0)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    f.Close()
    return h


def ratio_canvas(base, hu, hn, xlabel):
    c = ROOT.TCanvas("c_" + base, base, 720, 720); ROOT.SetOwnership(c, False)
    p1 = ROOT.TPad("p1", "", 0, 0.32, 1, 1); p1.SetBottomMargin(0.02); p1.Draw()
    p2 = ROOT.TPad("p2", "", 0, 0.0, 1, 0.32); p2.SetTopMargin(0.02); p2.SetBottomMargin(0.32); p2.SetGridy(); p2.Draw()
    ROOT.SetOwnership(p1, False); ROOT.SetOwnership(p2, False)
    p1.cd()
    hn.SetLineColor(ROOT.kAzure + 1); hn.SetLineWidth(2)     # NanoAOD (inclusive reference)
    hu.SetLineColor(ROOT.kOrange + 7); hu.SetLineWidth(2)    # UHH2 ntuple
    hn.SetMaximum(1.35 * max(hn.GetMaximum(), hu.GetMaximum())); hn.SetMinimum(0)
    hn.SetTitle(base); hn.GetYaxis().SetTitle("a.u. (unit area)"); hn.GetXaxis().SetLabelSize(0)
    hn.Draw("hist"); hu.Draw("hist same")
    leg = ROOT.TLegend(0.66, 0.76, 0.89, 0.90); leg.SetBorderSize(0); ROOT.SetOwnership(leg, False)
    leg.AddEntry(hn, "NanoAOD (inclusive)", "l"); leg.AddEntry(hu, "UHH2 ntuple", "l"); leg.Draw()
    c.cd(); p2.cd()
    hr = hu.Clone("r_" + base); hr.SetDirectory(0); hr.Divide(hn)
    hr.SetLineColor(ROOT.kBlack); hr.SetMarkerStyle(20); hr.SetMarkerSize(0.5); hr.SetTitle("")
    hr.SetMinimum(0.9); hr.SetMaximum(1.1)
    ya = hr.GetYaxis(); ya.SetTitle("ntuple/nano"); ya.SetNdivisions(505); ya.SetTitleSize(0.11); ya.SetTitleOffset(0.42); ya.SetLabelSize(0.09)
    xa = hr.GetXaxis(); xa.SetTitle(xlabel); xa.SetTitleSize(0.12); xa.SetLabelSize(0.10)
    hr.Draw("ep")
    line = ROOT.TLine(xa.GetXmin(), 1.0, xa.GetXmax(), 1.0); line.SetLineStyle(2); line.Draw(); ROOT.SetOwnership(line, False)
    c._keep = (p1, p2, hn, hu, hr, leg, line)
    return c


def frac(h, a, b):
    lo, hi = h.FindBin(a), h.FindBin(b)
    return h.Integral(lo, hi)          # h is unit-normalized, so this is the fraction in [a,b]


def main(uhh2_f, nano_f):
    pdf = "ntuple_inclusiveness.pdf"
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    c0 = ROOT.TCanvas(); c0.Print(pdf + "[")

    # --- gen_mtt: threshold behaviour is the decisive test ---
    hu = norm(uhh2_f, "gen_mtt", 100, 300, 2000, "u")
    hn = norm(nano_f, "gen_mtt", 100, 300, 2000, "n")
    c = ratio_canvas("gen_mtt", hu, hn, "gen m_{t#bar{t}} [GeV]"); c.Print(pdf)
    fu, fn = frac(hu, 300, 400), frac(hn, 300, 400)
    print("\n=== gen_mtt (normalized) ===")
    print("  fraction in threshold [300,400] GeV :  ntuple=%.4f   nano=%.4f   ntuple/nano=%.2f" %
          (fu, fn, fu / fn if fn else 0))
    print("  mean m_tt                           :  ntuple=%.1f   nano=%.1f  GeV" % (hu.GetMean(), hn.GetMean()))

    # --- gen cosTheta*: sculpting cross-check ---
    hu2 = norm(uhh2_f, "gen_cosThetaStar", 40, -1, 1, "u")
    hn2 = norm(nano_f, "gen_cosThetaStar", 40, -1, 1, "n")
    c2 = ratio_canvas("gen_cosThetaStar", hu2, hn2, "gen cos#theta*"); c2.Print(pdf)
    edge = 0.5 * (frac(hu2, -1.0, -0.8) + frac(hu2, 0.8, 1.0)) / (0.5 * (frac(hn2, -1.0, -0.8) + frac(hn2, 0.8, 1.0)) or 1)
    print("\n=== gen cosTheta* (normalized) ===")
    print("  forward/backward edge ntuple/nano   :  %.2f   (much < 1 => ntuple sculpted by a lepton-pt skim)" % edge)

    c0.Print(pdf + "]")
    print("\nwrote %s" % pdf)

    verdict = "INCLUSIVE  -> Path 1 (build the map with TTbarGenPy on the ntuples) is clean." \
        if (fn and 0.85 <= fu / fn <= 1.15) else \
        "SKIMMED    -> Path 2 (keep the NanoAOD map, align nanoGen<->TTbarGen copies)."
    print("\nverdict (from gen_mtt threshold ratio): %s" % verdict)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: python check_ntuple_inclusive.py gendump_uhh2.root gendump_nano.root"); sys.exit(1)
    main(sys.argv[1], sys.argv[2])