from __future__ import print_function
import ROOT
ROOT.gROOT.SetBatch(True)

f = ROOT.TFile.Open("spincorr_skim_TTToSemiLeptonic_UL18.root")
t = f.Get("entSkim")
# gen distributions must use the GENERATOR weight (genweight)
GENW = "genweight" if t.GetBranch("genweight") else "1" 

# --- 1) a weighted histogram, boosted region, reco vs gen ---------------------
hr = ROOT.TH1F("hr", "boosted #tilde{D} witness;cHel_P3n;events", 24, -1, 1)
hg = ROOT.TH1F("hg", "gen;cHel_P3n;events", 24, -1, 1)
# selection strings act directly on the branches; weight = eventweight
t.Draw("reco_cHel_P3n>>hr", "eventweight*(pass_reco && reco_M_tt>800)", "goff")
t.Draw("gen_cHel_P3n>>hg",  "%s*(gen_semilep && gen_M_tt>800)" % GENW, "goff")
print("entries  reco(boosted)=%d  gen(boosted)=%d" % (hr.GetEntries(), hg.GetEntries()))

# --- 2) a coefficient (D-tilde = 5 * <sign(cHel_P3n)>) via a TProfile ---------
# truth-binned boosted yardstick, identical (matched) events -> reproduces 1.26
pr = ROOT.TProfile("pr", "reco", 1, 0, 1)
pg = ROOT.TProfile("pg", "gen",  1, 0, 1)
sel = "pass_reco && gen_semilep && gen_M_tt>800"
t.Draw("5*(reco_cHel_P3n>=0 ? 1 : -1):0.5>>pr", "eventweight*(%s)" % sel, "goff")
t.Draw("5*(gen_cHel_P3n>=0 ? 1 : -1):0.5>>pg",  "%s*(%s)" % (GENW, sel),  "goff")
rr, gg = pr.GetBinContent(1), pg.GetBinContent(1)
print("boosted Dtilde: reco=%.4f  gen=%.4f  reco/gen=%.2f" % (rr, gg, rr/gg if gg else 0))

# --- 3) per-component comparison is equally easy, e.g.:
#   t.Draw("reco_cosTheta1k", "eventweight*(pass_reco && gen_semilep)")
#   t.Draw("gen_cosTheta1k",  "eventweight*(gen_semilep)")

# --- 4) systematic variation: reweight by a shifted-over-nominal SF ratio ----
# The skim carries the nominal SF and its up/down siblings, so a varied event weight is eventweight * (weight_X_up / weight_X_nominal)
hnom = ROOT.TH1F("hnom", "PU up vs nominal;cHel_P3n;events", 24, -1, 1)
hpu  = ROOT.TH1F("hpu",  "pu up", 24, -1, 1)
sel  = "pass_reco && reco_M_tt>800"
t.Draw("reco_cHel_P3n>>hnom", "eventweight*(%s)" % sel, "goff")
t.Draw("reco_cHel_P3n>>hpu",  "eventweight*(weight_pu_up/weight_pu)*(%s)" % sel, "goff")
print("PU-up integral / nominal integral = %.4f" % (hpu.Integral()/hnom.Integral() if hnom.Integral() else 0))
f.Close()


# =============================================================================
# Dump every skim variable to its own TH1F, and for each variable that has both
# a reco_ and gen_ version, build a reco-vs-gen overlay with a ratio panel below.
#   outputs:  skim_hists.root        -- all TH1F + the overlay canvases
#             skim_reco_vs_gen.pdf   -- one page per overlaid variable
# Physics variables are eventweighted and sentinel-cleaned (drop -10 and 99);
# weights/flags are shown as their raw distributions.
# =============================================================================
import math as _math
 
 
def _binning(name):
    """(nbins, lo, hi) for a variable by name pattern, or None to auto-range."""
    base = name.replace("reco_", "").replace("gen_", "")
    if base == "chi2":            return None              # reco-only, range unknown -> auto
    if base == "M_tt":            return (100, 300.,  2500.)
    if base == "pt_hadTop":       return (100,   0.,  1200.)
    if base == "beta":            return ( 50,   0.,     1.)
    if base == "dyreco":          return ( 80,  -4.,     4.)
    if base.startswith("Sigma_phi") or base.startswith("Delta_phi"):
        return (64, -_math.pi, _math.pi)
    if base.endswith("_plus") or base.endswith("_minus"):  # C cross-combos live in [-2,2]
        return (80, -2., 2.)
    return (50, -1., 1.)          # cosTheta*, cHel(_P3n)(+slices), C{ij}, cosThetaStar
 
 
def _valid(name):
    """selection dropping the -10 / 99 sentinels and requiring the right flag."""
    flag = "pass_reco" if name.startswith("reco_") else "gen_semilep"
    return "(%s && %s > -9.9 && %s < 98.)" % (flag, name, name)
 
 
def _ratio_canvas(c, base, hr, hg):
    """Draw reco (blue) vs gen (orange) unit-area shapes on top and the reco/gen ratio below, onto the REUSED canvas c. 
    Ownership of every drawn object is handed to ROOT (SetOwnership False) / detached from the file (SetDirectory 0) so Python's GC never double-frees them. 
    Returns the sub-objects so the caller can keep refs alive."""
    c.Clear()
    p1 = ROOT.TPad("p1", "", 0, 0.30, 1, 1.0); p1.SetBottomMargin(0.02); p1.Draw()
    p2 = ROOT.TPad("p2", "", 0, 0.0, 1, 0.30)
    p2.SetTopMargin(0.02); p2.SetBottomMargin(0.32); p2.SetGridy(); p2.Draw()
    ROOT.SetOwnership(p1, False); ROOT.SetOwnership(p2, False)
 
    p1.cd()
    hr.SetLineColor(ROOT.kAzure + 1);  hr.SetLineWidth(2)
    hg.SetLineColor(ROOT.kOrange + 7); hg.SetLineWidth(2)
    ymax = max(hr.GetMaximum(), hg.GetMaximum())
    hr.SetMaximum(1.25 * ymax if ymax > 0 else 1.0); hr.SetMinimum(0.)
    hr.SetTitle(base); hr.GetYaxis().SetTitle("a.u. (unit area)"); hr.GetXaxis().SetLabelSize(0)
    hr.Draw("hist"); hg.Draw("hist same")
    leg = ROOT.TLegend(0.70, 0.78, 0.89, 0.90); leg.SetBorderSize(0)
    leg.AddEntry(hr, "reco", "l"); leg.AddEntry(hg, "gen", "l"); leg.Draw()
    ROOT.SetOwnership(leg, False)
 
    c.cd(); p2.cd()
    hra = hr.Clone("ratio_" + base); hra.SetDirectory(0); hra.Divide(hg)
    hra.SetLineColor(ROOT.kBlack); hra.SetMarkerStyle(20); hra.SetMarkerSize(0.5); hra.SetTitle("")
    hra.SetMinimum(0.0); hra.SetMaximum(2.0)
    yax = hra.GetYaxis(); yax.SetTitle("reco/gen"); yax.SetNdivisions(505)
    yax.SetTitleSize(0.11); yax.SetTitleOffset(0.42); yax.SetLabelSize(0.09)
    xax = hra.GetXaxis(); xax.SetTitle(base); xax.SetTitleSize(0.12); xax.SetLabelSize(0.10)
    hra.Draw("ep")
    line = ROOT.TLine(xax.GetXmin(), 1.0, xax.GetXmax(), 1.0); line.SetLineStyle(2); line.Draw()
    ROOT.SetOwnership(line, False)
    c.cd()
    return (p1, p2, hra, leg, line)
 
 
def dump_all_histograms(infile="spincorr_skim_TTToSemiLeptonic_UL18.root",
                        outroot="skim_hists.root", outpdf="skim_reco_vs_gen.pdf"):
    import gc
    fin = ROOT.TFile.Open(infile)
    tt = fin.Get("entSkim")
    fout = ROOT.TFile(outroot, "RECREATE")
    ROOT.gStyle.SetOptStat(0)
    _prev = ROOT.gErrorIgnoreLevel
    ROOT.gErrorIgnoreLevel = ROOT.kWarning     # quiet the per-page "Info in <TCanvas::Print>" lines
 
    allbr = [b.GetName() for b in tt.GetListOfBranches()]
    FLAGS = ("pass_reco", "gen_semilep", "gen_channel")   # gen_* but NOT physics -> treat as raw
    reco  = sorted(b for b in allbr if b.startswith("reco_"))
    gen   = sorted(b for b in allbr if b.startswith("gen_") and b not in FLAGS)
    other = [b for b in allbr if b not in reco and b not in gen]   # eventweight, genweight, flags, weights
    GENW = "genweight" if tt.GetBranch("genweight") else "1"        # gen weight (fallback: unweighted)
 
    # 1) one TH1F per branch -> outroot
    for b in reco + gen:                                       # physics: eventweighted, sentinel-cleaned
        bn = _binning(b)
        spec = "(%d,%g,%g)" % bn if bn else ""
        wexpr = GENW if b.startswith("gen_") else "eventweight"   # reco SFs must NOT touch gen
        tt.Draw("%s>>h_%s%s" % (b, b, spec), "%s*%s" % (wexpr, _valid(b)), "goff")
        h = ROOT.gDirectory.Get("h_%s" % b)
        if h: h.SetTitle(b); h.SetXTitle(b); fout.cd(); h.Write()
    for b in other:                                           # weights/flags: raw distribution, auto range
        tt.Draw("%s>>h_%s" % (b, b), "", "goff")
        h = ROOT.gDirectory.Get("h_%s" % b)
        if h: h.SetTitle(b); h.SetXTitle(b); fout.cd(); h.Write()
 
    # 2) reco-vs-gen overlay + ratio -> outpdf (one page each) and into outroot as c_<var>.
    #    Reuse ONE canvas, disable the cyclic GC, and bracket-open/close the PDF so the document is finalized even if a variable is skipped
    shared = sorted(set(b[5:] for b in reco) & set(b[4:] for b in gen))
    gc.disable()
    keep = []                                  # strong refs so nothing is collected mid-loop
    c = ROOT.TCanvas("ratio_canvas", "", 700, 700); ROOT.SetOwnership(c, False)
    c.Print(outpdf + "[")                       # open the multipage document (adds no page)
    npage = 0
    for base in shared:
        nb, lo, hi = _binning("reco_" + base)
        tt.Draw("reco_%s>>ovr_%s(%d,%g,%g)" % (base, base, nb, lo, hi),
                "eventweight*%s" % _valid("reco_" + base), "goff")
        hr_ = ROOT.gDirectory.Get("ovr_%s" % base)
        tt.Draw("gen_%s>>ovg_%s(%d,%g,%g)" % (base, base, nb, lo, hi),
                "%s*%s" % (GENW, _valid("gen_" + base)), "goff")
        hg_ = ROOT.gDirectory.Get("ovg_%s" % base)
        if (not hr_) or (not hg_):
            continue
        hr_.SetDirectory(0); hg_.SetDirectory(0)
        if hr_.Integral() > 0: hr_.Scale(1.0 / hr_.Integral())   # unit-area shape comparison
        if hg_.Integral() > 0: hg_.Scale(1.0 / hg_.Integral())
        extra = _ratio_canvas(c, base, hr_, hg_)
        fout.cd(); c.Write("c_" + base)
        c.Print(outpdf)                         # add one page
        keep += [hr_, hg_] + list(extra)
        npage += 1
    c.Print(outpdf + "]")                        # close (finalize) the multipage document
    gc.enable()
    ROOT.gErrorIgnoreLevel = _prev
 
    fout.Close(); fin.Close()
    print("wrote %s : %d individual TH1F (%d reco + %d gen + %d weights/flags)"
          % (outroot, len(reco) + len(gen) + len(other), len(reco), len(gen), len(other)))
    print("wrote %s : %d reco-vs-gen overlay pages" % (outpdf, npage))
 
 
if __name__ == "__main__":
    dump_all_histograms()