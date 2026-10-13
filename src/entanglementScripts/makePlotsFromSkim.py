from __future__ import print_function
import ROOT
ROOT.gROOT.SetBatch(True)
import math as _math

# =============================================================================
# Dump every skim variable to its own TH1F, save it in a root file and PDF
#   outputs:  skim_hists.root        -- all TH1F + the overlay canvases
#             skim_reco_vs_gen.pdf   -- one page per branch
# Physics variables are event weighted and sentinel-cleaned (drop -10 and 99)
# weights/flags are shown as their raw distributions.
# =============================================================================


# Custom binning for some variables, otherwise auto-range (ROOT default)
def _binning(name):
    """(nbins, lo, hi) for a variable by name pattern, or None to auto-range."""
    base = name.replace("reco_", "").replace("gen_", "")

    if base == "chi2":      return (100,   0.,   100.)
    if base == "M_tt":      return (100, 300.,  2500.)
    if base == "pt_hadTop": return (100,   0.,  1500.)
    if base == "beta":      return ( 50,   0.,     1.)
    if base == "absDeltaY": return ( 36,  -5.,     5.)
    # Kinematics
    if base.endswith("_pt"):     return (100,  0.,  1500.)
    if base.endswith("_eta"):    return ( 50, -5.,     5.)
    if base.endswith("_phi"):    return ( 12, -_math.pi, _math.pi)
    if base.endswith("_E"):      return (100,   0., 2000.)
    if base.endswith("_charge"): return (  3, -1.5,   1.5)
    # Spin correlation variables
    if base.startswith("cHel"): return (6,   -1.,   1.)
    if base.startswith("cos")       or base.startswith("C"):         return (24,       -1.,       1.)
    if base.startswith("Sigma_phi") or base.startswith("Delta_phi"): return (12, -_math.pi, _math.pi)

    # topology and gen ttbar system flags (0 or 1)
    if base == "is_toptag_reconstruction" or base == "semilep": return (2,  -0.5, 1.5)
    # gen-decay channel (0=hadronic, 1-3=semilep, 4-9=dilep)
    if base == "channel":   return (10, -0.5, 9.5)
    # event weights
    if base == "eventweight":  return (200,   -1.,   1.)
    if base == "genweight":    return (800, -400., 400.) # genweight is pm 303.35800
    # systematic weights
    if base.startswith("weight_") or "weight" in base: return (500, 0., 5.) 

    return None  # auto-range for others, mostly weights


# Selections require:
#     gen_semilep == 1, means the gen-level ttbar system was found/defined as semileptonic
#     reco_chi2 < 30, means the reco-level semileptonic system passed the chi2 cut
#     toplogyFlag == 1, means there was a top-tagged jet in the event (i.e. merged topology)
#     obs != -10, means the physics observable was defined in the event (not a sentinel)
#     obs != 99, means the lepton projection was defined in the event (not a sentinel)
def _selections(name):
    """Event selection plus removal of the sentinel values that apply to this branch."""
    sel = ["gen_semilep", "reco_chi2 < 30", "reco_is_toptag_reconstruction"]
    # SENT (-10): physics observables only, never weights
    if name.startswith("reco_") or name.startswith("gen_"):
        sel.append("%s != -10" % name)
    # NOTp (99): only the lepton projections, set when that lepton charge isn't in the event
    if name.endswith("_Lep") or name.endswith("_antiLep"):
        sel.append("%s != 99" % name)
    return "(" + " && ".join(sel) + ")"
 
    
def write_hists2rootfile(infile= "/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/genTopStudy/TTToSemiLeptonicUL18_SCskim_wDoLC.root", 
                        outroot="/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/skimPlots/bothLeptons/Skim_wLastCopy_hists.root"):
    print("")
    fin = ROOT.TFile.Open(infile)
    TTree = fin.Get("SC_SkimTree")
    fout = ROOT.TFile(outroot, "RECREATE")
    print("Opened input file: %s" % infile)
    print("Writing to output file: %s" % outroot)

    # Event weights
    genWeight = "genweight" if TTree.GetBranch("genweight") else "1" 
    if TTree.GetBranch("genweight"): print("Found genweight branch, will use it for gen histograms")
    else: print("No genweight branch found, gen histograms will be unweighted")

    recoWeight = "eventweight" if TTree.GetBranch("eventweight") else "1"
    if TTree.GetBranch("eventweight"): print("Found eventweight branch, will use it for reco histograms")
    else: print("No eventweight branch found, reco histograms will be unweighted")
 
    # Arrays sorting all, flags, reco, gen, and other branches
    allbr = [branch.GetName() for branch in TTree.GetListOfBranches()]
    FLAGS = ("gen_semilep", "gen_channel", "reco_is_toptag_reconstruction")
    reco  = sorted(branch for branch in allbr if branch.startswith("reco_")and branch not in FLAGS)
    gen   = sorted(branch for branch in allbr if branch.startswith("gen_") and branch not in FLAGS)
    other = [branch for branch in allbr if branch not in reco and branch not in gen]   # eventweight, genweight, flags, and syst weights
 
    # ==================================== #
    # Write one TH1F per branch to outroot #
    # ==================================== #
    # gen kinematics
    for branch in gen:
        binning = _binning(branch)
        spec = "(%d,%g,%g)" % binning if binning else ""
        TTree.Draw("%s>>h_%s%s" % (branch, branch, spec), "%s*%s" % (genWeight, _selections(branch)), "goff")
        hist = ROOT.gDirectory.Get("h_%s" % branch)
        if hist: hist.SetTitle(branch); hist.SetXTitle(branch); fout.cd(); hist.Write()
    # gen kinematics
    for branch in reco:
        binning = _binning(branch)
        spec = "(%d,%g,%g)" % binning if binning else ""
        TTree.Draw("%s>>h_%s%s" % (branch, branch, spec), "%s*%s" % (recoWeight, _selections(branch)), "goff")
        hist = ROOT.gDirectory.Get("h_%s" % branch)
        if hist: hist.SetTitle(branch); hist.SetXTitle(branch); fout.cd(); hist.Write()
    # weights and flags
    for branch in other:
        binning = _binning(branch)
        spec = "(%d,%g,%g)" % binning if binning else ""
        TTree.Draw("%s>>h_%s%s" % (branch, branch, spec), "%s" % _selections(branch), "goff")
        hist = ROOT.gDirectory.Get("h_%s" % branch)
        if hist: hist.SetTitle(branch); hist.SetXTitle(branch); fout.cd(); hist.Write()
    fout.Close()
    fin.Close()

    # Summary
    print("wrote %s : %d individual TH1F (%d reco + %d gen + %d weights/flags)" % (outroot, len(reco) + len(gen) + len(other), len(reco), len(gen), len(other)))
 

def make_histPlots(infile="/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/genTopStudy/TTToSemiLeptonicUL18_SCskim_wDoLC.root", 
                    outpdf="/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/skimPlots/bothLeptons/Skim_wLastCopy_hists.pdf"):
    import gc
    _prev = ROOT.gErrorIgnoreLevel         # restore at end of function
    ROOT.gErrorIgnoreLevel = ROOT.kWarning # quiet the per-page "Info in <TCanvas::Print>" lines

    print("")
    fin = ROOT.TFile.Open(infile)
    TTree = fin.Get("SC_SkimTree")
    ROOT.gStyle.SetOptStat("iemou")        # entries, integral, mean, rms, underflow, overflow
    print("Opened input file: %s" % infile)
    print("Writing hitograms to: %s" % outpdf)

    # Event weights
    genWeight = "genweight" if TTree.GetBranch("genweight") else "1" 
    if TTree.GetBranch("genweight"): print("Found genweight branch, will use it for gen histograms")
    else: print("No genweight branch found, gen histograms will be unweighted")

    recoWeight = "eventweight" if TTree.GetBranch("eventweight") else "1"
    if TTree.GetBranch("eventweight"): print("Found eventweight branch, will use it for reco histograms")
    else: print("No eventweight branch found, reco histograms will be unweighted")
 
    # Arrays sorting all, flags, reco, gen, and other branches
    allbr = [branch.GetName() for branch in TTree.GetListOfBranches()]
    FLAGS = ("gen_semilep", "gen_channel", "reco_is_toptag_reconstruction")
    reco  = sorted(branch for branch in allbr if branch.startswith("reco_")and branch not in FLAGS)
    gen   = sorted(branch for branch in allbr if branch.startswith("gen_") and branch not in FLAGS)
    other = [branch for branch in allbr if branch not in reco and branch not in gen]   # eventweight, genweight, flags, and syst weights

    # ============================================================================ #
    # Plot the hist and save in outpdf (one page each) and into outroot as c_<var> #
    # ============================================================================ #
    # Reuse ONE canvas, disable the cyclic GC, and bracket-open/close the PDF so the document is finalized even if a variable is skipped
    gc.disable()
    canv = ROOT.TCanvas("canvas", "", 600, 600); ROOT.SetOwnership(canv, False)
    canv.Print(outpdf + "[")  # open the multipage document (adds no page)
    keep = [] # strong refs so nothing is collected mid-loop
    npage = 0
    for branch in allbr:
        # print status
        print("Plotting %s branch" % branch)

        # binning
        binning = _binning(branch)
        spec = "(%d,%g,%g)" % binning if binning else ""

        # Draw gen hist (normalized to unit area)
        if branch in gen:
            TTree.Draw("%s>>h_%s%s" % (branch, branch, spec),"%s*%s" % (genWeight, _selections(branch)), "goff")
            hist = ROOT.gDirectory.Get("h_%s" % branch)
            if hist.Integral() > 0: hist.Scale(1.0 / hist.Integral())
        # Draw reco hist
        elif branch in reco:
            TTree.Draw("%s>>h_%s%s" % (branch, branch, spec),"%s*%s" % (recoWeight, _selections(branch)), "goff")
            hist = ROOT.gDirectory.Get("h_%s" % branch)
        # Draw weights/flags hist (normalized to unit area)
        else:
            TTree.Draw("%s>>h_%s%s" % (branch, branch, spec), "%s" % _selections(branch), "goff")
            hist = ROOT.gDirectory.Get("h_%s" % branch)
            if hist.Integral() > 0: hist.Scale(1.0 / hist.Integral())

        # check hist exists
        if not hist:
            continue
        hist.SetDirectory(0)

        canv.Clear()
        p1 = ROOT.TPad("p1", "", 0, 0, 1.0, 1.0); p1.Draw()
        ROOT.SetOwnership(p1, False)
        p1.cd()

        hist.SetTitle(branch)
        if branch in reco: hist.GetYaxis().SetTitle("Events")
        else: hist.GetYaxis().SetTitle("au (normalized)")
        hist.GetXaxis().SetTitle(branch)
        hist.SetLineWidth(2)
        if branch in gen: hist.SetLineColor(ROOT.kRed); hist.SetFillColor(ROOT.kRed); hist.SetFillStyle(0)
        elif branch in reco: hist.SetLineColor(ROOT.kRed); hist.SetFillColor(ROOT.kRed); hist.SetFillStyle(1001)
        elif branch in other: hist.SetLineColor(ROOT.kBlack); hist.SetFillColor(ROOT.kGray); hist.SetFillStyle(3254); ROOT.gStyle.SetHatchesSpacing(0.75)
        ymax = hist.GetMaximum() * 1.2; hist.SetMaximum(ymax); hist.SetMinimum(0.0)
        hist.Draw("HIST")
        canv.cd(); canv.Update()
        canv.Print(outpdf) # add one page
        keep += [hist, p1] 
        npage += 1
    canv.Print(outpdf + "]") # close (finalize) the multipage document
    gc.enable()
    ROOT.gErrorIgnoreLevel = _prev
    # Summary
    print("wrote %s : %d pages" % (outpdf, npage))
 
 
if __name__ == "__main__":
    make_histPlots()
    write_hists2rootfile()