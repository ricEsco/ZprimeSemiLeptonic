# strength_map.py -- alpha(m_tt, cosTheta*) spin-correlation strength maps for the reweighting method (docs/reweighting_method.md), one TProfile2D per observable.
#
# alpha = prefactor * <sign(O)> per (m_tt, cosTheta*) cell: each map is a TProfile2D filled with prefactor*sign(O) weighted by the generator weight, 
# so the per-cell MEAN is alpha and the per-cell ERROR is its uncertainty (matches "value + variance" in the note).
#
# Prefactors are the lb (lepton + hadronic-b) values from extractCoeff.py coeff_multipliers
# -- these already fold in the b-quark spin-analyzing power for our semileptonic analyzers.
#
# Built from the INCLUSIVE NanoAOD gen (nanoGen with the fiducial gate DISABLED), so the map is defined over the full gen production phase space.
#
# Output ROOT file is plain TProfile2D -> readable by the EL7/ROOT 6.14 skim as a lookup table.
from __future__ import print_function
import ROOT
ROOT.TH1.AddDirectory(False)   # keep the profiles out of any open TFile until we write them

# observable (nanoGen lb_ variable name) -> prefactor (extractCoeff.py coeff_multipliers, lb pairing)
MAP_MULT = {
    "lb_cHel":       5.0,   # D (helicity opening angle)
    "lb_cHel_P3n":   5.0,   # D-tilde (boosted entanglement witness)
    "lb_Cnn": -10.0, "lb_Cnr": -10.0, "lb_Cnk": -10.0,   # C_ij (9)
    "lb_Crn": -10.0, "lb_Crr": -10.0, "lb_Crk": -10.0,
    "lb_Ckn": -10.0, "lb_Ckr": -10.0, "lb_Ckk": -10.0,
    "lb_cos_theta1k": 5.0, "lb_cos_theta1r": 5.0, "lb_cos_theta1n": 5.0,   # B1_i (top polarization)
    "lb_cos_theta2k": 5.0, "lb_cos_theta2r": 5.0, "lb_cos_theta2n": 5.0,   # B2_i (antitop polarization)
}

# production-axis binning (docs/reweighting_method.md, Section 7)
MTT_NBINS, MTT_LO, MTT_HI = 52, 300.0, 2000.0
CTS_NBINS, CTS_LO, CTS_HI = 48, -1.0, 1.0

# book one TProfile2D per observable, once on import; accumulated across all input files
STRENGTH_MAPS = {}
for _name in MAP_MULT:
    STRENGTH_MAPS[_name] = ROOT.TProfile2D(
        "alpha_" + _name,
        "#alpha " + _name + ";m_{t#bar{t}} [GeV];cos#theta*;#alpha",
        MTT_NBINS, MTT_LO, MTT_HI, CTS_NBINS, CTS_LO, CTS_HI)


def fill(mtt, cosTS, genweight, obs):
    """Fill every map for one gen event.
    obs = {observable_name: gen value of O}. Fills prefactor*sign(O) weighted by genweight,
    so each (m_tt, cosTheta*) cell accumulates the per-cell mean = prefactor*<sign(O)> = alpha."""
    for name, prof in STRENGTH_MAPS.items():
        val = obs.get(name)
        if val is None or val != val:            # skip missing / NaN
            continue
        prof.Fill(mtt, cosTS, MAP_MULT[name] * (1.0 if val >= 0.0 else -1.0), genweight)


def write(out_path, process="TTToSemiLeptonic"):
    """Write all maps into out_path under a <process>/ directory (mirrors the note's process/var layout).
    Lookup later: TFile.Open(out_path).Get(process + '/alpha_' + observable)."""
    fout = ROOT.TFile(out_path, "RECREATE")
    d = fout.mkdir(process)
    d.cd()
    for prof in STRENGTH_MAPS.values():
        prof.Write()
    fout.Close()
    print("wrote %d strength maps to %s:%s/" % (len(STRENGTH_MAPS), out_path, process))
