"""
Brief description:
A modular script that will book, fill, and write TProfile2D strength maps ("alpha" maps) to a ROOT file.

Maps are created for the following Observables:
 - Correlation Matrix elements C_ij     (9)
 - Top/antitop polarization B1_i, B2_i  (6)
 - Entanglement witnesses D and D-tilde (2)
These maps are used in the template reweighting method (docs/reweighting_method.md, Section 7).

alpha(m_tt, cosTheta*) = prefactor * <gen_weight * sign_Obs> 
-> <gen_weight * sign_Obs> = per-cell MEAN  = frac{sum(gen_weight * sign_Obs)} {sum(gen_weight)} = frac{N(Obs>0) - N(Obs<0)}{N(Obs>0) + N(Obs<0)} = Asymmetry at gen-level
-> sigma_alpha                              = per-cell ERROR = frac{s}{sqrt[sum(gen_weight)]} : s = std. dev -> ref: https://root.cern/doc/v636/classTProfile2D.html
Prefactor scales the measured Asymmetry strength to the corresponding Observable's SC variable strength

Output ROOT file is plain TProfile2D -> readable by the EL7/ROOT 6.14 skim as a lookup table.
"""
from __future__ import print_function
import ROOT
ROOT.TH1.AddDirectory(False)   # keep the profiles out of any open TFile until we write them

# Map of Observable -> prefactor 
PrefactorMap = {
    "lb_cos_theta1k": 2.0, "lb_cos_theta1r": 2.0, "lb_cos_theta1n": 2.0,   # B1_i (top polarization)
    "lb_cos_theta2k": 2.0, "lb_cos_theta2r": 2.0, "lb_cos_theta2n": 2.0,   # B2_i (antitop polarization)
    "lb_Cnn": 4.0, "lb_Cnr": 4.0, "lb_Cnk": 4.0,                           # C_ij (9)
    "lb_Crn": 4.0, "lb_Crr": 4.0, "lb_Crk": 4.0,
    "lb_Ckn": 4.0, "lb_Ckr": 4.0, "lb_Ckk": 4.0,
    "lb_cHel":       2.0,                        # D       (threshold entanglement witness)
    "lb_cHel_P3n":   2.0,                        # D-tilde (boosted entanglement witness)
}

# Binning for the 2D maps: m_tt vs cosTheta*
Mtt_NBINS, Mtt_LO, Mtt_HI = 52, 300.0, 2000.0    # 52 bins of 32 GeV width from Mtt [300, 2000] GeV
CosTS_NBINS, CosTS_LO, CosTS_HI = 48, -1.0, 1.0  # 48 bins of 0.0416667 width from cosTheta* [-1, 1]

# Book one TProfile2D per observable (once on import-> accumulated across all input files)
STRENGTH_MAPS = {}
for Obs_name in PrefactorMap:
    # TProfile2D("name", "title;xaxis;y-axis;z-axis", nxbins, xlow, xup, nybins, ylow, yup)
    STRENGTH_MAPS[Obs_name] = ROOT.TProfile2D(
        "alpha_" + Obs_name,
        "#alpha " + Obs_name + ";m_{t#bar{t}} [GeV];cos#theta*;#alpha",
        Mtt_NBINS, Mtt_LO, Mtt_HI, CosTS_NBINS, CosTS_LO, CosTS_HI)

# Fill the maps for every gen event
def fill(mtt, cosTS, genweight, obs):
    """
    Fill every map for one gen event with obs = {Obs_name: gen value of Obs} 
    Fills prefactor*sign(O) weighted by genweight, so each (m_tt, cosTS) cell accumulates the per-cell mean = prefactor*<sign(Obs)> = alpha
    """
    for Obs_name, TProf2D in STRENGTH_MAPS.items():
        Obs_val = obs.get(Obs_name)
        # skip if Observable value is missing / NaN
        if Obs_val is None or Obs_val != Obs_val:
            continue 
        TProf2D.Fill(mtt, cosTS, PrefactorMap[Obs_name] * (1.0 if Obs_val >= 0.0 else -1.0), genweight)

# Write maps to output ROOT file AFTER loop over all input files
def write(out_path, process="TTToSemiLeptonic"):
    """
    Write all maps into out_path under a <process>/ directory
    Lookup syntax: TFile.Open(out_path).Get(process + '/alpha_' + Observable_name)
    """
    fout = ROOT.TFile(out_path, "RECREATE")
    dir = fout.mkdir(process)
    dir.cd()
    for TProf2D in STRENGTH_MAPS.values():
        TProf2D.Write()
    fout.Close()
    print("wrote %d strength maps to %s:%s/" % (len(STRENGTH_MAPS), out_path, process))