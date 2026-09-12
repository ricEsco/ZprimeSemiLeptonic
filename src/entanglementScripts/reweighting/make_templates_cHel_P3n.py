"""
Description:
Build the reco-level templates for Combine using the reweighting method (docs/reweighting_method.md)
following the dilepton entanglement analysis arXiv:2406.03976

Procedure: 
for every event in the boosted-central RECO region:
  1) look up alpha(m_tt, cosTheta*) using the Ntuple alpha_maps (built from INCLUSIVE gen-kinematics)
  2) read its gen-level Observable (gen_<Obs>)
  3) construct the per-event reweight factor and fill the reco-level Observable (reco_<Obs>) with: recoWeight * w(f)
     f = +1  -> "NoSC up" (spin corr removed)
     f = 0   -> nominal
     f = -1  -> "NoSC down" (spin corr enhanced)

Implementation:  
    density  d  = 1 + alpha * Obs
    varied  d_f = 1 + alpha * (1-f) * Obs
    w(f) = d_f / d

Bare pyROOT (plain skim + plain TProfile2D map); py2.7/ROOT 6.14 (EL7 container) and py3 both fine.
==================================================================================================================================
* this script currently targets cHel_P3n, the boosted/central entanglement proxy, but can be generalized to other SC Observables *
* boosted + central -> (|cosTheta*| < 0.4 & m_tt > 800) -> pre-windowed branch                                                   *
* Binning: 6 bins on [-1, 1]                                                                                                     *
==================================================================================================================================
"""
from __future__ import print_function
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Inputs
skimFile = "TTToSemiLeptonicUL18_SCskim.root"
skimTree = "SC_SkimTree"
mapFile  = "strengthmaps_uhh2/NoSC_alphaStrengthMap_lb.root"
mapObj   = "TTToSemiLeptonic/alpha_lb_cHel_P3n" # map lookup syntax: TFile.Open(out_path).Get(<process> + '/alpha_' + <Observable_name>)

# Outputs
OUT_ROOT = "cHel_P3n_templates.root" # ROOT file with stored up/nominal/down templates
OUT_PDF  = "cHel_P3n_templates.pdf"  # PDF with the templates and their ratios to nominal


# reco-level boosted/central (m_tt>800 & |cosTheta*|<0.4) Observable is defined within the RECO-level phase space
RECO_Obs = "reco_cHel_P3n_Mtt800_Inf_cosThetaLT0p4"

# Inclusive gen-level inputs for reweighting
GEN_Obs   = "gen_cHel_P3n"
GEN_Mtt   = "gen_M_tt"
GEN_cosTS = "gen_cosThetaStar"

# Binning for the reco-level templates (6 bins on [-1, 1])
NBINS, XLO, XHI = 6, -1.0, 1.0

# Degree of "NoSC" variation (keep f inside [-1, +1] to avoid negative weights)
f_up   = +1.0 # remove spin correlations
f_down = -1.0 # enhance spin correlations

# Spin analyzer product for (lepton, hadronic-b) pair
kappaProduct   = 0.4  # kappa_a * kappa_b  for (lepton, hadronic-b)
alphaPrefactor = 2.0 # PrefactorMap = 2 for cHel_P3n (D-tilde)

# Precaution flags
SENTINEL   = -3.0           # branch value <= SENTINEL means "not available / outside window"
D_FLOOR    = 1e-6           # guard on |density|
MIN_MAPENT = 10             # flag map cells with fewer entries than this (noisy alpha)

# Precaution clamp function to keep (m_tt, cosTheta*) inside the map range for genuinely-boosted / edge events
def clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)


##### Main function that looks up alpha, reads skim, books and fills the reco-level templates #####
#############################################################################################
def main(skim=skimFile, mapf=mapFile, out_root=OUT_ROOT, out_pdf=OUT_PDF):

    # load alpha_map(gen m_tt, cosTheta*) and extract its axis ranges
    mapFile = ROOT.TFile.Open(mapf)
    if (not mapFile) or mapFile.IsZombie():
        raise RuntimeError("cannot open strength map: %s" % mapf)
    TProf2D = mapFile.Get(mapObj)
    if not TProf2D:
        raise RuntimeError("cannot find %s in %s" % (mapObj, mapf))
    TProf2D.SetDirectory(0)
    mapFile.Close()
    ax_m, ax_c = TProf2D.GetXaxis(), TProf2D.GetYaxis()
    m_lo, m_hi = ax_m.GetXmin(), ax_m.GetXmax()
    c_lo, c_hi = ax_c.GetXmin(), ax_c.GetXmax()
    # epsilon * bin width to keep (m_tt, cosTheta*) inside the map range
    eps_m = 1e-3 * (m_hi - m_lo) / TProf2D.GetNbinsX()
    eps_c = 1e-3 * (c_hi - c_lo) / TProf2D.GetNbinsY()
    # look up alpha
    def alpha_lookup(mtt, cts):
        """
        alpha (=kappa_a*kappa_b*D =2*A) at gen-level (m_tt, cosTheta*)
        clamped into range so extremely-boosted (m_tt>2000) / edge events use the nearest valid cell instead of under/overflow.
        """
        bin = TProf2D.FindBin(clamp(mtt, m_lo + eps_m, m_hi - eps_m),
                              clamp(cts, c_lo + eps_c, c_hi - eps_c))
        return TProf2D.GetBinContent(bin), TProf2D.GetBinEntries(bin)

    # Read skim file and enable only the branches we need for the reweighting
    skimFile = ROOT.TFile.Open(skim)
    if (not skimFile) or skimFile.IsZombie():
        raise RuntimeError("cannot open skim: %s" % skim)
    skimTTree = skimFile.Get(skimTree)
    skimTTree.SetBranchStatus("*", 0)
    for branch in ("eventweight", "reco_is_toptag_reconstruction", "gen_semilep", RECO_Obs, GEN_Obs, GEN_Mtt, GEN_cosTS):
        skimTTree.SetBranchStatus(branch, 1)
    setFilter = "reco_is_toptag_reconstruction > 0 && %s > %g" % (RECO_Obs, SENTINEL) # filter makes sure reco Obs is in Merged Topology and was well defined
    skimTTree.Draw(">>eventList", setFilter, "goff")                                  # list of events that pass filter
    eventList = ROOT.gDirectory.Get("eventList")
    npass = eventList.GetN()
    print("Applied '%s' filter -> looping over %d events passed" % (setFilter, npass))

    # Book the reco-level templates (nominal, NoSC up, NoSC down)
    def book(name, title):
        hist = ROOT.TH1F(name, title, NBINS, XLO, XHI)
        hist.Sumw2()
        hist.SetDirectory(0)
        return hist
    h_nom = book("cHel_P3n",          "nominal;cHel_{P3n}^{RECO};Events")
    h_up  = book("cHel_P3n_NoSCUp",   "NoSC up   (f=+1);cHel_{P3n}^{RECO};Events")
    h_dn  = book("cHel_P3n_NoSCDown", "NoSC down (f=-1);cHel_{P3n}^{RECO};Events")

    # Counters for reweighting diagnostics
    n_gen_used = n_guard = n_neg_up = n_neg_dn = n_lowstat = 0
    density_min = 1e30 # will ultimately save smallest density value
    
    # Store components of Asymmetry fraction: A = (N+ - N-) / (N+ + N-) = sum(weight*sign(Obs)) / sum(weight)
    sumWeights = {"nom": [0.0, 0.0], "up": [0.0, 0.0], "down": [0.0, 0.0]}  # [sum(weight*sign(Obs)), sum(weight)]
    # Accumulate Asymmetry components fo calculate entanglement witness
    def accumulate(tag, weight, Obs):
        sumWeights[tag][0] += weight * (1.0 if Obs >= 0.0 else -1.0) # SignedSum
        sumWeights[tag][1] += weight                                 # sumOfWeights

    # Loop over skimmed events to fill the reco-level templates
    for i in range(npass):
        skimTTree.GetEntry(eventList.GetEntry(i))
        recoWeight  = skimTTree.eventweight
        recoObs = getattr(skimTTree, RECO_Obs)
        upWeight = downWeight = 1.0
        # Only construct the reweighting factor if the gen-level Observable was defined and set
        if skimTTree.gen_semilep and getattr(skimTTree, GEN_Obs) > SENTINEL:
            genObs = getattr(skimTTree, GEN_Obs)
            alpha, ent = alpha_lookup(getattr(skimTTree, GEN_Mtt), getattr(skimTTree, GEN_cosTS))
            if ent < MIN_MAPENT: # count low-stat map cells -> noisy alpha values
                n_lowstat += 1
            # density of Observable
            density = 1.0 + alpha * genObs
            # density = 1.0 + alpha*kappaProduct * genObs
            # save smallest density value for diagnostics
            if density_min > density:
                density_min = density
            # save the number of events where |density| < 1e-6
            if abs(density) < D_FLOOR:
                n_guard += 1
            # construct up/down weights
            else:
                n_gen_used += 1
                upWeight   = (1.0 + alpha * genObs * (1.0 - f_up)   ) / density
                downWeight = (1.0 + alpha * genObs * (1.0 - f_down) ) / density
                # Count negative up/down weights
                if upWeight < 0.0: n_neg_up += 1
                if downWeight < 0.0: n_neg_dn += 1
        # Fill reco-level Observable and reweighted up/down templates
        h_nom.Fill(recoObs, recoWeight)
        h_up.Fill(recoObs, recoWeight * upWeight)
        h_dn.Fill(recoObs, recoWeight * downWeight)
        accumulate("nom",  recoWeight,              recoObs)
        accumulate("up",   recoWeight * upWeight,   recoObs)
        accumulate("down", recoWeight * downWeight, recoObs)
    skimFile.Close()

    # Write the templates to output ROOT file
    outFile = ROOT.TFile(out_root, "RECREATE")
    for hist in (h_nom, h_up, h_dn):
        hist.Write()
    outFile.Close()

    # Convert alpha to entanglement witness
    # alpha = prefactor * Asymmetry = prefactor * sum(weight*sign(Obs)) / sum(weight)
    # D = alpha / (kappa_a * kappa_b) = alpha / kappaProduct
    def coeff(tag):
        signedSum    = sumWeights[tag][0]
        sumOfWeights = sumWeights[tag][1]
        return (alphaPrefactor/kappaProduct) * (signedSum/sumOfWeights) if sumOfWeights else 0.0
    

    # Print summary of yields, D-tilde estimator, and reweighting diagnostics to the terminal
    print("\n=================== event yields  ===================")
    print("  nominal: %.1f   NoSC_up: %.1f   NoSC_down: %.1f" % (h_nom.Integral(), h_up.Integral(), h_dn.Integral()))
    coeff_nom, coeff_up, coeff_down = coeff("nom"), coeff("up"), coeff("down")
    print("\n=== reco D-tilde estimate  (prefactor/kappaProduct) * Asymmetry ===")
    print("  nominal      : %+.4f" % coeff_nom)
    print("  NoSC_up      : %+.4f  ->spin correlation removed" % coeff_up)
    print("  NoSC_down    : %+.4f  ->spin correlation enhanced" % coeff_down)
    print("  --> ordering should be: |up| < |nominal| < |down| : %s" %
          ("OK" if abs(coeff_up) < abs(coeff_nom) < abs(coeff_down) else "REVERSED -> flip density sign"))
    print("\n=== reweighting weights diagnostic (gen-semileptonic subset) ===")
    print("  events reweighted          : %d"       % n_gen_used)
    print("  min density (1+alpha*Obs)  : %+.4f"    % density_min)
    print("  |density| < %g counter     : %d"       % (D_FLOOR, n_guard))
    print("  negative up / down weights : %d / %d"  % (n_neg_up, n_neg_dn))
    print("  low-stat map cells (<%d)   : %d"       % (MIN_MAPENT, n_lowstat))
    print("\nwrote %s  (cHel_P3n, cHel_P3n_NoSCup, cHel_P3n_NoSCdown)" % out_root)

    # Plot the templates and their ratios to nominal, with D-tilde in the legend
    make_plot(h_nom, h_up, h_dn, out_pdf, coeff_nom, coeff_up, coeff_down)
    print("wrote %s" % out_pdf)

# plotting function
def make_plot(h_nom, h_up, h_dn, out_pdf, coeff_nom, coeff_up, coeff_down):
    # setup canvas and pads
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas("c_tmpl", "templates", 720, 780); ROOT.SetOwnership(canvas, False)
    pad1 = ROOT.TPad("pad1", "", 0, 0.30, 1, 1); pad1.SetBottomMargin(0.02); pad1.Draw()
    pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.30); pad2.SetTopMargin(0.02); pad2.SetBottomMargin(0.34); pad2.SetGridy(); pad2.Draw()
    ROOT.SetOwnership(pad1, False); ROOT.SetOwnership(pad2, False)
    
    pad1.cd() # top pad
    # Colors and line styles
    h_nom.SetLineColor(ROOT.kBlack);     h_nom.SetLineWidth(2); h_nom.SetMarkerStyle(20); h_nom.SetMarkerSize(0.7)
    h_up.SetLineColor(ROOT.kAzure + 1);  h_up.SetLineWidth(2)
    h_dn.SetLineColor(ROOT.kOrange + 7); h_dn.SetLineWidth(2)
    # y-axis
    ymax = max(h_nom.GetMaximum(), h_up.GetMaximum(), h_dn.GetMaximum())
    h_nom.SetMaximum(1.40 * ymax)
    h_nom.SetMinimum(0.0)
    # titles
    h_nom.SetTitle("Merged + boosted + central (Top-tagged, m_{t#bar{t}}>800, |cos#theta*|<0.4)")
    h_nom.GetXaxis().SetLabelSize(0)
    h_nom.GetYaxis().SetTitle("Events")
    # Draw templates
    h_nom.Draw("hist e"); h_up.Draw("hist same"); h_dn.Draw("hist same"); h_nom.Draw("hist e same")
    # legend
    leg = ROOT.TLegend(0.15, 0.70, 0.55, 0.90); leg.SetBorderSize(0); leg.SetFillStyle(0); ROOT.SetOwnership(leg, False)
    leg.AddEntry(h_up,  "NoSC up   (f=+1) #tilde{D}=%+.3f" % coeff_up, "le")
    leg.AddEntry(h_nom, "nominal   (f= 0) #tilde{D}=%+.3f" % coeff_nom, "le")
    leg.AddEntry(h_dn,  "NoSC down (f=-1) #tilde{D}=%+.3f" % coeff_down, "le")
    leg.Draw()
    
    pad2.cd() # bottom pad (ratio wrt nominal)
    ratio_up = h_up.Clone("ratio_up"); ratio_up.SetDirectory(0); ratio_up.Divide(h_nom)
    ratio_down = h_dn.Clone("ratio_down"); ratio_down.SetDirectory(0); ratio_down.Divide(h_nom)
    
    ratio_up.SetTitle("")
    # y-axis
    ratio_up.SetMinimum(0.8); ratio_up.SetMaximum(1.2)
    ratio_up.GetYaxis().SetTitle("NoSC / nom"); ratio_up.GetYaxis().SetNdivisions(505)
    ratio_up.GetYaxis().SetTitleSize(0.12); ratio_up.GetYaxis().SetTitleOffset(0.40); ratio_up.GetYaxis().SetLabelSize(0.10)
    # x-axis
    ratio_up.GetXaxis().SetTitle("cHel_{P3n}^{RECO}"); ratio_up.GetXaxis().SetTitleSize(0.13); ratio_up.GetXaxis().SetLabelSize(0.11)
    # Draw ratios and horizontal line at 1
    ratio_up.Draw("hist")
    ratio_down.Draw("hist same")
    line = ROOT.TLine(XLO, 1.0, XHI, 1.0); line.SetLineStyle(2); line.Draw(); ROOT.SetOwnership(line, False)

    # canvas._keep = (pad1, pad2, ratio_up, ratio_down, leg, line)
    canvas._keep = (pad1, pad2, ratio_up, ratio_down, leg, line)
    # the _keep attribute is a trick to prevent the objects from being garbage collected when the function exits, which can happen in PyROOT. 
    # By storing references to these objects in an attribute of the canvas, we ensure they remain alive for the lifetime of the canvas.
    canvas.Print(out_pdf)

if __name__ == "__main__":
    a = sys.argv[1:]
    kw = {}
    if len(a) >= 1: kw["skim"] = a[0]
    if len(a) >= 2: kw["mapf"] = a[1]
    if len(a) >= 3: kw["out_root"] = a[2]
    if len(a) >= 4: kw["out_pdf"] = a[3]
    main(**kw)
