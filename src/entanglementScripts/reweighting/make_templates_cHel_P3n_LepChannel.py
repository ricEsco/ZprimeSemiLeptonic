"""
Description:
Build the reco-level templates for Combine.
This is the lepton-flavor specific version.

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

# Choose a lepton flavor: electron or muon
LepFlavor = "muon"

# Inputs
if LepFlavor == "electron":
    # skimFile = "TTToSemiLeptonicUL18_SCskim_electron.root"
    skimFile = "../genTopStudy/TTToSemiLeptonicUL18_SCskim_wDoLC_electron.root"
elif LepFlavor == "muon":
    # skimFile = "TTToSemiLeptonicUL18_SCskim_muon.root"
    skimFile = "../genTopStudy/TTToSemiLeptonicUL18_SCskim_wDoLC_muon.root"
skimTree = "SC_SkimTree"

# Alpha maps are inclusive by construction -> no lepton-flavor splitting
# mapFile  = "strengthmaps_uhh2/NoSC_alphaStrengthMap_lb.root"
mapFile  = "strengthmaps_wLastCopy/NoSC_alphaStrengthMap_lb_wLastCopy.root"
mapObj   = "TTToSemiLeptonic/alpha_lb_cHel_P3n"

# Outputs
if LepFlavor == "electron":
    # OUT_ROOT = "cHel_P3n_templates_electron.root"
    # OUT_PDF  = "cHel_P3n_templates_electron.pdf"
    OUT_ROOT = "cHel_P3n_templates_wLastCopy_electron.root"
    OUT_PDF  = "cHel_P3n_templates_wLastCopy_electron.pdf"
elif LepFlavor == "muon":
    # OUT_ROOT = "cHel_P3n_templates_muon.root"
    # OUT_PDF  = "cHel_P3n_templates_muon.pdf"
    OUT_ROOT = "cHel_P3n_templates_wLastCopy_muon.root"
    OUT_PDF  = "cHel_P3n_templates_wLastCopy_muon.pdf"


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
    print("Reading alpha map from \n%s:%s" % (mapf, mapObj))
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
    print("Reading skim from \n%s:%s" % (skim, skimTree))
    skimFile = ROOT.TFile.Open(skim)
    if (not skimFile) or skimFile.IsZombie():
        raise RuntimeError("cannot open skim: %s" % skim)
    skimTTree = skimFile.Get(skimTree)
    skimTTree.SetBranchStatus("*", 0)
    for branch in ("eventweight", "reco_is_toptag_reconstruction", "reco_chi2", "gen_semilep", RECO_Obs, GEN_Obs, GEN_Mtt, GEN_cosTS):
        skimTTree.SetBranchStatus(branch, 1)
    # Event filter defines SR: Top-tagged && Chi2<30 && not SENT value
    setFilter = "reco_is_toptag_reconstruction > 0 && reco_chi2 < 30. && %s > %g" % (RECO_Obs, SENTINEL)
    skimTTree.Draw(">>eventList", setFilter, "goff") # list of events that pass filter
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
    # Manually book up/down weight distributions for finer binning
    h_up_wt = ROOT.TH1F("cHel_P3n_NoSCUp_wt",   "NoSCUp (f=+1) weights;NoSCUp_weight;Events",     100, -5.0, 5.0)
    h_dn_wt = ROOT.TH1F("cHel_P3n_NoSCDown_wt", "NoSCDown (f=-1) weights;NoSCDown_weight;Events", 100, -5.0, 5.0)
    h_up_wt.SetDirectory(0); h_dn_wt.SetDirectory(0)

    # Counters for reweighting diagnostics
    n_gen_used = n_guard = n_neg_up = n_neg_dn = n_lowstat = 0
    density_min = 1e30 # will ultimately save smallest density value
    
    # Store components of Asymmetry fraction: A = (N+ - N-) / (N+ + N-) = sum(weight*sign(Obs)) / sum(weight)
    sumWeights = {"nom": [0.0, 0.0], "up": [0.0, 0.0], "down": [0.0, 0.0]}  # [sum(weight*sign(Obs)), sum(weight)]
    # Accumulate Asymmetry components fo calculate entanglement witness
    def accumulate(tag, weight, Obs):
        sumWeights[tag][0] += weight * (1.0 if Obs >= 0.0 else -1.0) # SignedSum
        sumWeights[tag][1] += weight                                 # sumOfWeights

    # firstTen = 0
    # Loop over skimmed events to fill the reco-level templates
    for i in range(npass):
        # firstTen += 1
        skimTTree.GetEntry(eventList.GetEntry(i))
        # if firstTen <= 10:
        #     print("Currently processing event %d" % skimTTree.GetReadEntry())
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
                if upWeight < 0.0: 
                    n_neg_up += 1
                    print("Negative up weight is %.4f: density %.4f = 1 + alpha %.4f * genObs %.4f and eventweight is %.4f" % (upWeight, density, alpha, genObs, recoWeight))
                if downWeight < 0.0: 
                    n_neg_dn += 1
                    print("Negative down weight is %.4f: density %.4f = 1 + alpha %.4f * genObs %.4f and eventweight is %.4f" % (downWeight, density, alpha, genObs, recoWeight))
        # Fill reco-level Observable and reweighted up/down templates
        h_nom.Fill(recoObs, recoWeight)
        h_up.Fill(recoObs, recoWeight * upWeight)
        h_dn.Fill(recoObs, recoWeight * downWeight)
        accumulate("nom",  recoWeight,              recoObs)
        accumulate("up",   recoWeight * upWeight,   recoObs)
        accumulate("down", recoWeight * downWeight, recoObs)
        h_up_wt.Fill(upWeight); h_dn_wt.Fill(downWeight)
    skimFile.Close()

    # Write the templates to output ROOT file
    outFile = ROOT.TFile(out_root, "RECREATE")
    for hist in (h_nom, h_up, h_dn, h_up_wt, h_dn_wt):
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
    print("\nwrote %s  (cHel_P3n, cHel_P3n_NoSCup, cHel_P3n_NoSCdown, cHel_P3n_NoSCUp_wt, cHel_P3n_NoSCDown_wt)" % out_root)

    # Plot the templates and their ratios to nominal, with D-tilde in the legend
    make_plot(h_nom, h_up, h_dn, out_pdf, coeff_nom, coeff_up, coeff_down)
    print("wrote %s" % out_pdf)


"""
electron channel output:
Reading alpha map from
strengthmaps_uhh2/NoSC_alphaStrengthMap_lb.root:TTToSemiLeptonic/alpha_lb_cHel_P3n
Reading skim from
TTToSemiLeptonicUL18_SCskim_electron.root:SC_SkimTree
Applied 'reco_is_toptag_reconstruction > 0 && reco_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 > -3' filter -> looping over 67053 events passed

=================== event yields  ===================
  nominal: 1895.7   NoSC_up: 1892.2   NoSC_down: 1899.2

=== reco D-tilde estimate  (prefactor/kappaProduct) * Asymmetry ===
  nominal      : +0.4785
  NoSC_up      : +0.1418  ->spin correlation removed
  NoSC_down    : +0.8139  ->spin correlation enhanced
  --> ordering should be: |up| < |nominal| < |down| : OK

=== reweighting weights diagnostic (gen-semileptonic subset) ===
  events reweighted          : 67029
  min density (1+alpha*Obs)  : +0.3841
  |density| < 1e-06 counter     : 0
  negative up / down weights : 0 / 5
  low-stat map cells (<10)   : 0

wrote cHel_P3n_templates_electron.root  (cHel_P3n, cHel_P3n_NoSCup, cHel_P3n_NoSCdown)
wrote cHel_P3n_templates_electron.pdf
"""

"""
compare electron output to using isLastCopy top quarks:
Reading alpha map from
strengthmaps_wLastCopy/NoSC_alphaStrengthMap_lb_wLastCopy.root:TTToSemiLeptonic/alpha_lb_cHel_P3n
Reading skim from
../genTopStudy/TTToSemiLeptonicUL18_SCskim_wDoLC_electron.root:SC_SkimTree
Applied 'reco_is_toptag_reconstruction > 0 && reco_chi2 < 30. && reco_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 > -3' filter -> looping over 67053 events passed
Negative down weight is -0.6623: density 0.3756 = 1 + alpha 1.0189 * genObs -0.6128 and eventweight is 0.0282
Negative down weight is -0.0243: density 0.4940 = 1 + alpha 0.5137 * genObs -0.9851 and eventweight is 0.0117
Negative down weight is -0.2701: density 0.4405 = 1 + alpha 0.6392 * genObs -0.8753 and eventweight is 0.0162
Negative down weight is -1.0396: density 0.3290 = 1 + alpha 0.8070 * genObs -0.8315 and eventweight is 0.0274
Negative down weight is -0.1771: density 0.4593 = 1 + alpha 0.6763 * genObs -0.7995 and eventweight is -0.0066
Negative down weight is -0.0481: density 0.4882 = 1 + alpha 0.7634 * genObs -0.6704 and eventweight is 0.0266
Negative down weight is -0.8571: density 0.3500 = 1 + alpha 1.0189 * genObs -0.6380 and eventweight is 0.0126
Negative down weight is -0.2184: density 0.4508 = 1 + alpha 0.6122 * genObs -0.8971 and eventweight is 0.0360
Negative down weight is -0.7126: density 0.3686 = 1 + alpha 0.7097 * genObs -0.8896 and eventweight is 0.0129
Negative down weight is -2.7740: density 0.2095 = 1 + alpha 0.9020 * genObs -0.8765 and eventweight is 0.0177
Negative down weight is -0.0818: density 0.4804 = 1 + alpha 0.5487 * genObs -0.9471 and eventweight is 0.0263
Negative down weight is -0.1108: density 0.4738 = 1 + alpha 0.6457 * genObs -0.8150 and eventweight is 0.0166
Negative down weight is -0.4259: density 0.4122 = 1 + alpha 0.6119 * genObs -0.9606 and eventweight is 0.0128
Negative down weight is -0.3134: density 0.4323 = 1 + alpha 0.6457 * genObs -0.8793 and eventweight is -0.0141
Negative down weight is -0.2086: density 0.4528 = 1 + alpha 0.6667 * genObs -0.8208 and eventweight is -0.0195
Negative down weight is -0.0391: density 0.4904 = 1 + alpha 0.6392 * genObs -0.7973 and eventweight is 0.0305
Negative down weight is -1.8986: density 0.2565 = 1 + alpha 0.7634 * genObs -0.9740 and eventweight is 0.0224
Negative down weight is -0.1527: density 0.4645 = 1 + alpha 0.8070 * genObs -0.6635 and eventweight is 0.0268
Negative down weight is -0.4697: density 0.4049 = 1 + alpha 0.7097 * genObs -0.8385 and eventweight is 0.0184
Negative down weight is -0.5586: density 0.3908 = 1 + alpha 0.6122 * genObs -0.9950 and eventweight is -0.0288
Negative down weight is -0.0011: density 0.4997 = 1 + alpha 0.8070 * genObs -0.6199 and eventweight is 0.0142
Negative down weight is -0.0926: density 0.4779 = 1 + alpha 0.6122 * genObs -0.8528 and eventweight is 0.0211
Negative down weight is -1.0776: density 0.3249 = 1 + alpha 0.8070 * genObs -0.8365 and eventweight is 0.0309
Negative down weight is -1.9675: density 0.2520 = 1 + alpha 0.7634 * genObs -0.9798 and eventweight is 0.0345
Negative down weight is -0.0478: density 0.4883 = 1 + alpha 0.5246 * genObs -0.9754 and eventweight is -0.0201

=================== event yields  ===================
  nominal: 1895.7   NoSC_up: 1896.0   NoSC_down: 1895.5

=== reco D-tilde estimate  (prefactor/kappaProduct) * Asymmetry ===
  nominal      : +0.4785
  NoSC_up      : +0.1088  ->spin correlation removed
  NoSC_down    : +0.8483  ->spin correlation enhanced
  --> ordering should be: |up| < |nominal| < |down| : OK

=== reweighting weights diagnostic (gen-semileptonic subset) ===
  events reweighted          : 67029
  min density (1+alpha*Obs)  : +0.2095
  |density| < 1e-06 counter     : 0
  negative up / down weights : 0 / 25
  low-stat map cells (<10)   : 0

wrote cHel_P3n_templates_wLastCopy_electron.root  (cHel_P3n, cHel_P3n_NoSCup, cHel_P3n_NoSCdown, cHel_P3n_NoSCUp_wt, cHel_P3n_NoSCDown_wt)
wrote cHel_P3n_templates_wLastCopy_electron.pdf

real	0m9.793s
user	0m8.981s
sys	0m0.562s
"""

"""
---------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

"""
muon channel output:
Reading alpha map from
strengthmaps_uhh2/NoSC_alphaStrengthMap_lb.root:TTToSemiLeptonic/alpha_lb_cHel_P3n
Reading skim from
TTToSemiLeptonicUL18_SCskim_muon.root:SC_SkimTree
Applied 'reco_is_toptag_reconstruction > 0 && reco_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 > -3' filter -> looping over 89859 events passed

=================== event yields  ===================
  nominal: 2506.4   NoSC_up: 2502.9   NoSC_down: 2509.9

=== reco D-tilde estimate  (prefactor/kappaProduct) * Asymmetry ===
  nominal      : +0.4156
  NoSC_up      : +0.0770  ->spin correlation removed
  NoSC_down    : +0.7533  ->spin correlation enhanced
  --> ordering should be: |up| < |nominal| < |down| : OK

=== reweighting weights diagnostic (gen-semileptonic subset) ===
  events reweighted          : 89839
  min density (1+alpha*Obs)  : +0.3578
  |density| < 1e-06 counter     : 0
  negative up / down weights : 0 / 11
  low-stat map cells (<10)   : 0

wrote cHel_P3n_templates_muon.root  (cHel_P3n, cHel_P3n_NoSCup, cHel_P3n_NoSCdown)
wrote cHel_P3n_templates_muon.pdf
"""

"""
compare electron output to using isLastCopy top quarks:
Reading alpha map from
strengthmaps_wLastCopy/NoSC_alphaStrengthMap_lb_wLastCopy.root:TTToSemiLeptonic/alpha_lb_cHel_P3n
Reading skim from
../genTopStudy/TTToSemiLeptonicUL18_SCskim_wDoLC_muon.root:SC_SkimTree
Applied 'reco_is_toptag_reconstruction > 0 && reco_chi2 < 30. && reco_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 > -3' filter -> looping over 89859 events passed
Negative down weight is -0.2858: density 0.4375 = 1 + alpha 0.7097 * genObs -0.7926 and eventweight is 0.0195
Negative down weight is -2.9833: density 0.2007 = 1 + alpha 0.8493 * genObs -0.9411 and eventweight is 0.0178
Negative down weight is -9.4865: density 0.0871 = 1 + alpha 1.0189 * genObs -0.8960 and eventweight is -0.0194 <-- second largest negative down weight in the muon channel skim
Negative down weight is -0.0018: density 0.4995 = 1 + alpha 0.6122 * genObs -0.8174 and eventweight is 0.0212
Negative down weight is -0.0564: density 0.4863 = 1 + alpha 0.5455 * genObs -0.9418 and eventweight is 0.0171
Negative down weight is -0.5311: density 0.3951 = 1 + alpha 0.6122 * genObs -0.9880 and eventweight is 0.0128
Negative down weight is -0.0376: density 0.4908 = 1 + alpha 0.5314 * genObs -0.9583 and eventweight is 0.0215
Negative down weight is -0.8348: density 0.3528 = 1 + alpha 0.6667 * genObs -0.9709 and eventweight is 0.0220
Negative down weight is -1.3579: density 0.2978 = 1 + alpha 1.0189 * genObs -0.6892 and eventweight is 0.0250
Negative down weight is -0.0481: density 0.4882 = 1 + alpha 0.5344 * genObs -0.9577 and eventweight is 0.0204
Negative down weight is -0.4506: density 0.4081 = 1 + alpha 0.6119 * genObs -0.9674 and eventweight is -0.0203
Negative down weight is -0.0154: density 0.4962 = 1 + alpha 0.6122 * genObs -0.8229 and eventweight is -0.0181
Negative down weight is -0.3317: density 0.4289 = 1 + alpha 1.0189 * genObs -0.5606 and eventweight is 0.0162
Negative down weight is -0.6891: density 0.3719 = 1 + alpha 0.7634 * genObs -0.8229 and eventweight is -0.0249
Negative down weight is -0.2021: density 0.4541 = 1 + alpha 0.6457 * genObs -0.8455 and eventweight is -0.0266
Negative down weight is -7.3544: density 0.1069 = 1 + alpha 1.0189 * genObs -0.8766 and eventweight is -0.0237 <-- third largest negative down weight in the muon channel skim
Negative down weight is -0.0731: density 0.4824 = 1 + alpha 0.5833 * genObs -0.8874 and eventweight is 0.0178
Negative down weight is -2.6212: density 0.2164 = 1 + alpha 1.0189 * genObs -0.7691 and eventweight is 0.0136
Negative down weight is -0.1555: density 0.4639 = 1 + alpha 0.5763 * genObs -0.9302 and eventweight is 0.0192
Negative down weight is -0.2604: density 0.4424 = 1 + alpha 0.5841 * genObs -0.9547 and eventweight is 0.0176
Negative down weight is -0.0988: density 0.4765 = 1 + alpha 0.5424 * genObs -0.9653 and eventweight is 0.0340
Negative down weight is -1.6546: density 0.2736 = 1 + alpha 1.0189 * genObs -0.7129 and eventweight is 0.0126
Negative down weight is -0.9943: density 0.3340 = 1 + alpha 0.6763 * genObs -0.9849 and eventweight is 0.0289
Negative down weight is -0.4501: density 0.4081 = 1 + alpha 0.6077 * genObs -0.9739 and eventweight is 0.0203
Negative down weight is -0.4438: density 0.4092 = 1 + alpha 0.6119 * genObs -0.9656 and eventweight is 0.0217
Negative down weight is -1.2805: density 0.3048 = 1 + alpha 0.7634 * genObs -0.9107 and eventweight is 0.0206
Negative down weight is -20.1335: density 0.0452 = 1 + alpha 1.0189 * genObs -0.9371 and eventweight is 0.0166 <-- first largest negative down weight in the muon channel skim
Negative down weight is -0.0291: density 0.4928 = 1 + alpha 0.6119 * genObs -0.8289 and eventweight is 0.0215
Negative down weight is -0.2494: density 0.4446 = 1 + alpha 0.6392 * genObs -0.8690 and eventweight is 0.0131
Negative down weight is -0.1028: density 0.4756 = 1 + alpha 0.6122 * genObs -0.8566 and eventweight is 0.0244
Negative down weight is -0.0013: density 0.4997 = 1 + alpha 1.0189 * genObs -0.4911 and eventweight is 0.0276
Negative down weight is -2.6362: density 0.2157 = 1 + alpha 1.0189 * genObs -0.7698 and eventweight is 0.0124
Negative down weight is -0.1751: density 0.4598 = 1 + alpha 0.6122 * genObs -0.8824 and eventweight is 0.0225
Negative down weight is -0.0775: density 0.4814 = 1 + alpha 0.5254 * genObs -0.9871 and eventweight is 0.0161
Negative down weight is -5.1154: density 0.1405 = 1 + alpha 1.0189 * genObs -0.8435 and eventweight is 0.0281
Negative down weight is -0.0548: density 0.4867 = 1 + alpha 0.5254 * genObs -0.9770 and eventweight is 0.0194
Negative down weight is -0.2786: density 0.4389 = 1 + alpha 1.0189 * genObs -0.5507 and eventweight is -0.0180
Negative down weight is -0.5971: density 0.3851 = 1 + alpha 1.0189 * genObs -0.6036 and eventweight is -0.0232
Negative down weight is -1.5668: density 0.2804 = 1 + alpha 0.7797 * genObs -0.9230 and eventweight is -0.0326
Negative down weight is -0.0611: density 0.4852 = 1 + alpha 0.5246 * genObs -0.9814 and eventweight is 0.0160
Negative down weight is -0.0859: density 0.4794 = 1 + alpha 0.5841 * genObs -0.8913 and eventweight is 0.0210
Negative down weight is -0.7377: density 0.3653 = 1 + alpha 0.7097 * genObs -0.8944 and eventweight is 0.0195
Negative down weight is -0.1278: density 0.4700 = 1 + alpha 0.6763 * genObs -0.7838 and eventweight is 0.0271
Negative down weight is -1.1956: density 0.3129 = 1 + alpha 0.7634 * genObs -0.9001 and eventweight is -0.0202

=================== event yields  ===================
  nominal: 2506.4   NoSC_up: 2507.9   NoSC_down: 2504.9

=== reco D-tilde estimate  (prefactor/kappaProduct) * Asymmetry ===
  nominal      : +0.4156
  NoSC_up      : +0.0440  ->spin correlation removed
  NoSC_down    : +0.7877  ->spin correlation enhanced
  --> ordering should be: |up| < |nominal| < |down| : OK

=== reweighting weights diagnostic (gen-semileptonic subset) ===
  events reweighted          : 89839
  min density (1+alpha*Obs)  : +0.0452
  |density| < 1e-06 counter     : 0
  negative up / down weights : 0 / 44
  low-stat map cells (<10)   : 0

wrote cHel_P3n_templates_wLastCopy_muon.root  (cHel_P3n, cHel_P3n_NoSCup, cHel_P3n_NoSCdown, cHel_P3n_NoSCUp_wt, cHel_P3n_NoSCDown_wt)
wrote cHel_P3n_templates_wLastCopy_muon.pdf

real	0m12.283s
user	0m11.337s
sys	0m0.611s
"""




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
