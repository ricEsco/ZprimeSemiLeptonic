"""
Description:
Skim the Analysis output files (TTToSemiLeptonic, UL18) into a light flat TTree carrying the reco and gen spin-correlation variables
needed for the entanglement template reweighting method (docs/reweighting_method.md).

What this script does, per event:
  1) Copies the reco spin-correlation branches straight through from the AnalysisTree (computed by ZprimeSemiLeptonicModules.cxx).
  2) Recomputes the SAME set of variables at gen level directly from GenParticles via TTbarGenPy (identical to the C++ code), 
     exact same procedure used to construct alpha the inclusive gen-level strength maps (strength_map.py) for the template reweighting method.
  3) Writes everything into ONE flat output TTree ("SC_SkimTree").

Scope / design choices:
  - every tree entry is written so all selection can be done offline when building templates, e.g. pass/fail flags are stored instead of pre-filtering
  - reco_<Observable> = the reco branch already computed by the analysis framework, copied verbatim
  - gen_<Observable>  = the same physical quantity, computed from each event's ttbar system gen-particles

Topology-cut included for the boosted signal region:
  The "Merged" signal region is defined by the presence of a top-tagged jet and identified by the boolean is_toptag_reconstruction during candidate building.
  Downstream template-building code will need this flag to reproduce the same signal region as analysis output.

Sentinel conventions (mirror ZprimeSemiLeptonicModules.cxx exactly):
    -10 : variable not available this event
            reco: not reconstructed/b-tagged
            gen: not semileptonic OR a masked slice whose window is not satisfied
     99 : a (anti)lepton-only projection where that (anti)lepton is NOT present, e.g. cosTheta2k_Lep on an l+ event

CMSSW_10_6_28 / Python 2.7 / ROOT 6.14  (run inside the EL7 container)
"""
from __future__ import print_function
import ROOT
import math
import sys
from array import array
from ttbargen_py import TTbarGenPy

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# py2/py3 shim: define xrange regardless of python version, so the event loop below works in both.
try:
    xrange  # type: ignore
except NameError:
    xrange = range


# ============================================================================================================================== #
#                                                          Configuration                                                         #
# ============================================================================================================================== #
# Absolute path to TTToSemiLeptonic analysis output base directory
BASE = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/reprocessedSignalWtopologycut/"

# All four files make up the full TTToSemiLeptonic sample and all map to the "TTbar" process
TAGS = ["electron", "electron2", "muon", "muon2"]
ALL_FILES = [BASE + "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root" % tag for tag in TAGS]

# Uses the condor $(Process) = {0,1,2,3} index to skim each file in it's own condor job
# if no argument is given, all four files are skimmed together into one output
if len(sys.argv) > 1:
    idx = int(sys.argv[1])
    allFILES = [ALL_FILES[idx]]
    outF  = "TTToSemiLeptonicUL18_SCskim_%s.root" % TAGS[idx]
else:
    allFILES = ALL_FILES
    outF  = "TTToSemiLeptonicUL18_SCskim.root"

NMAX = -1                  # cap on events processed per file; -1 = full pass (use e.g. 200000 for a quick test)
SENT = -10.0               # sentinel: indicates event Obs outside kinematic window (e.g. Mtt>800) or gen-event not semileptonic
NOTp = 99.0                # sentinel: (anti-)lepton specific observable where (anti-)lepton is not present
pT_HADTOP_Thresh = 150.0   # low/high pT split for the Baumgart variables, see ZprimeSemiLeptonicModules.cxx:1517


# ====================================================================================================================== #
#                                                      Branch lists                                                      #
# ====================================================================================================================== #
# Reconstruction level spin-correlation branches
RECO = [
    # Reconstruction quality and event topology flag
    "chi2", "is_toptag_reconstruction", # is_toptag_reconstruction converted from bool->float for consistency

    # ttbar-system kinematics and charge asymmetry
    "M_tt", "beta", "absDeltaY", "cos_PosTop_beam", "pt_hadTop",

    # lepton-only projections and their *Star (charge-asymmetry) versions
    "cosTheta1k_antiLep", "cosTheta1r_antiLep", "cosTheta1n_antiLep", "cosTheta1kStar_antiLep", "cosTheta1rStar_antiLep",
    "cosTheta2k_Lep",     "cosTheta2r_Lep",     "cosTheta2n_Lep",     "cosTheta2kStar_Lep",     "cosTheta2rStar_Lep",

    # spin-analyzer 1/2 projections and their *Star asymmetry versions
    "cosTheta1k", "cosTheta1r", "cosTheta1n", "cosTheta1kStar", "cosTheta1rStar",
    "cosTheta2k", "cosTheta2r", "cosTheta2n", "cosTheta2kStar", "cosTheta2rStar",

    # correlation matrix elements and their sum/difference combinations
    "Cnn", "Cnr", "Cnk", "Crn", "Crr", "Crk", "Ckn", "Ckr", "Ckk",
    "Crk_plus", "Crk_minus", "Cnr_plus", "Cnr_minus", "Cnk_plus", "Cnk_minus",

    # entanglement witnesses and their differential (windowed) slices
    "cHel",     "cHel_Mtt300_400",     "cHel_Mtt300_400_betaLT0p9",
    "cHel_P3n", "cHel_P3n_Mtt800_Inf", "cHel_P3n_Mtt800_Inf_cosThetaLT0p4",

    # Baumgart et al. angular variables
    "Sigma_phi", "Sigma_phi_low", "Sigma_phi_high",
    "Delta_phi", "Delta_phi_low", "Delta_phi_high",
]

# Generator level spin correlation branches
GEN = [
    # ttbar-system kinematics and charge asymmetry
    "M_tt", "beta", "absDeltaY", "cosThetaStar", "pt_hadTop", "Lepton_charge",

    # lepton-only projections and *Star (99 on the slot the lepton is not in, exactly like reco)
    "cosTheta1k_antiLep", "cosTheta1r_antiLep", "cosTheta1n_antiLep", "cosTheta1kStar_antiLep", "cosTheta1rStar_antiLep",
    "cosTheta2k_Lep",     "cosTheta2r_Lep",     "cosTheta2n_Lep",     "cosTheta2kStar_Lep",     "cosTheta2rStar_Lep",

    # spin-analyzer 1/2 projections and their *Star asymmetry versions
    "cosTheta1k", "cosTheta1r", "cosTheta1n", "cosTheta1kStar", "cosTheta1rStar",
    "cosTheta2k", "cosTheta2r", "cosTheta2n", "cosTheta2kStar", "cosTheta2rStar",

    # correlation matrix elements and their sum/difference combinations
    "Cnn", "Cnr", "Cnk", "Crn", "Crr", "Crk", "Ckn", "Ckr", "Ckk",
    "Crk_plus", "Crk_minus", "Cnr_plus", "Cnr_minus", "Cnk_plus", "Cnk_minus",

    # entanglement witnesses and their differential (windowed) slices
    "cHel",     "cHel_Mtt300_400",     "cHel_Mtt300_400_betaLT0p9",
    "cHel_P3n", "cHel_P3n_Mtt800_Inf", "cHel_P3n_Mtt800_Inf_cosThetaLT0p4",

    # Baumgart et al. angular variables
    "Sigma_phi", "Sigma_phi_low", "Sigma_phi_high",
    "Delta_phi", "Delta_phi_low", "Delta_phi_high",
]

# 4-vectors for gen-level ttbar-system particles in Lab Frame (pt/eta/phi/E)
GEN_PARTICLES = ["Top", "Antitop", "bHad", "Lepton"]

# Systematic-correction weights (also copied through verbatim (no prefix))
# NB: Since the analysis runs per lepton flavor we adopt the following protocol for the weights;
# e.g. since weight_sfelec_* branches only exist in the electron files and weight_sfmu_* only in the muon files,
#      a weight absent from a given input file is written as 1.0 (neutral) for every event in that file.
WEIGHTS_RAW = """
prefiringweight prefiringweightdown prefiringweightup
weight_btagdisc_central
weight_btagdisc_cferr1_down weight_btagdisc_cferr1_up
weight_btagdisc_cferr2_down weight_btagdisc_cferr2_up
weight_btagdisc_hf_down weight_btagdisc_hf_up
weight_btagdisc_hfstats1_down weight_btagdisc_hfstats1_up
weight_btagdisc_hfstats2_down weight_btagdisc_hfstats2_up
weight_btagdisc_lf_down weight_btagdisc_lf_up
weight_btagdisc_lfstats1_down weight_btagdisc_lfstats1_up
weight_btagdisc_lfstats2_down weight_btagdisc_lfstats2_up
weight_fsr_2_down weight_fsr_2_up weight_fsr_4_down weight_fsr_4_up
weight_fsr_g2gg_cns_down weight_fsr_g2gg_cns_up weight_fsr_g2gg_mur_down weight_fsr_g2gg_mur_up
weight_fsr_g2qq_cns_down weight_fsr_g2qq_cns_up weight_fsr_g2qq_mur_down weight_fsr_g2qq_mur_up
weight_fsr_q2qg_cns_down weight_fsr_q2qg_cns_up weight_fsr_q2qg_mur_down weight_fsr_q2qg_mur_up
weight_fsr_sqrt2_down weight_fsr_sqrt2_up
weight_fsr_x2xg_cns_down weight_fsr_x2xg_cns_up weight_fsr_x2xg_mur_down weight_fsr_x2xg_mur_up
weight_isr_2_down weight_isr_2_up weight_isr_4_down weight_isr_4_up
weight_isr_g2gg_cns_down weight_isr_g2gg_cns_up weight_isr_g2gg_mur_down weight_isr_g2gg_mur_up
weight_isr_g2qq_cns_down weight_isr_g2qq_cns_up weight_isr_g2qq_mur_down weight_isr_g2qq_mur_up
weight_isr_q2qg_cns_down weight_isr_q2qg_cns_up weight_isr_q2qg_mur_down weight_isr_q2qg_mur_up
weight_isr_sqrt2_down weight_isr_sqrt2_up
weight_isr_x2xg_cns_down weight_isr_x2xg_cns_up weight_isr_x2xg_mur_down weight_isr_x2xg_mur_up
weight_isrfsr_2_down weight_isrfsr_2_up weight_isrfsr_4_down weight_isrfsr_4_up
weight_isrfsr_sqrt2_down weight_isrfsr_sqrt2_up
weight_murmuf_downdown weight_murmuf_downnone
weight_murmuf_dyn1_downdown weight_murmuf_dyn1_downnone weight_murmuf_dyn1_nonedown weight_murmuf_dyn1_noneup weight_murmuf_dyn1_upnone weight_murmuf_dyn1_upup
weight_murmuf_dyn2_downdown weight_murmuf_dyn2_downnone weight_murmuf_dyn2_nonedown weight_murmuf_dyn2_noneup weight_murmuf_dyn2_upnone weight_murmuf_dyn2_upup
weight_murmuf_dyn3_downdown weight_murmuf_dyn3_downnone weight_murmuf_dyn3_nonedown weight_murmuf_dyn3_noneup weight_murmuf_dyn3_upnone weight_murmuf_dyn3_upup
weight_murmuf_dyn4_downdown weight_murmuf_dyn4_downnone weight_murmuf_dyn4_nonedown weight_murmuf_dyn4_noneup weight_murmuf_dyn4_upnone weight_murmuf_dyn4_upup
weight_murmuf_nonedown weight_murmuf_noneup weight_murmuf_upnone weight_murmuf_upup
weight_pu weight_pu_down weight_pu_up
weight_sfelec_id weight_sfelec_id_down weight_sfelec_id_up
weight_sfelec_reco weight_sfelec_reco_down weight_sfelec_reco_up
weight_sfelec_trigger weight_sfelec_trigger_down weight_sfelec_trigger_up
weight_sfmu_id_stat weight_sfmu_id_stat_down weight_sfmu_id_stat_up
weight_sfmu_id_syst weight_sfmu_id_syst_down weight_sfmu_id_syst_up
weight_sfmu_iso_stat weight_sfmu_iso_stat_down weight_sfmu_iso_stat_up
weight_sfmu_iso_syst weight_sfmu_iso_syst_down weight_sfmu_iso_syst_up
weight_sfmu_reco weight_sfmu_reco_down weight_sfmu_reco_up
weight_sfmu_trigger_stat weight_sfmu_trigger_stat_down weight_sfmu_trigger_stat_up
weight_sfmu_trigger_syst weight_sfmu_trigger_syst_down weight_sfmu_trigger_syst_up
weight_topmistagsf weight_topmistagsf_down weight_topmistagsf_up
weight_toppt_a_down weight_toppt_a_up weight_toppt_b_down weight_toppt_b_up weight_toppt_nominal
weight_toptagsf weight_toptagsf_corr_down weight_toptagsf_corr_up weight_toptagsf_uncorr_down weight_toptagsf_uncorr_up
"""
WEIGHTS = [w for w in WEIGHTS_RAW.split() if w]


# ============================================================================== #
#                                Helper functions                                #
# ============================================================================== #
def TLvec(genp):
    """Convert a TTbarGenPy generator-level particle to a ROOT TLorentzVector."""
    particle = genp.v4()
    vector = ROOT.TLorentzVector()
    vector.SetPtEtaPhiE(particle.pt(), particle.eta(), particle.phi(), particle.energy())
    return vector


def _wrap(angle):
    """Wrap an angle back into (-pi, pi], applying a single +/- 2*pi correction"""
    if angle > math.pi:
        angle -= 2.0 * math.pi
    if angle < -math.pi:
        angle += 2.0 * math.pi
    return angle


def gen_observables(ttg):
    """
    Gen-level spin correlation variables
    Returns a dictionary keyed by the names in GEN branches
    """
    lep = ttg.ChargedLepton()         # lepton from W decay
    qpos = (lep.pdgId() < 0)          # lepton charge -> True(False) = Negative(Positive) -> top(antitop) decays leptonically
    TopV = TLvec(ttg.Top())           # top quark TLorentzVector
    antiTopV = TLvec(ttg.Antitop())   # antitop quark TLorentzVector
    LepV = TLvec(lep)                 # lepton TLorentzVectors
    hadBquarkV = TLvec(ttg.BHad())    # hadronic b-quark TLorentzVector
    ttbar = TopV + antiTopV           # ttbar system TLorentzVector
    Mtt = ttbar.M()                   # ttbar invariant mass
    Ett = TopV.E() + antiTopV.E()     # ttbar system energy

    # Lab Frame quantities: beta, absDeltaY, pt_hadTop
    beta = abs(TopV.Pz() + antiTopV.Pz()) / Ett if Ett > 0 else SENT # relativistic-beta of the ttbar system
    absDeltaY = abs(TopV.Rapidity()) - abs(antiTopV.Rapidity())         # |y_t| - |y_tbar|, the charge asymmetry
    pt_hadTop = TLvec(ttg.TopHad()).Pt()                             # pT of the hadronically-decaying top

    # Boost all 4-vectors to the ttbar CM frame
    CMboostV = -ttbar.BoostVector()
    for v in (TopV, antiTopV, LepV, hadBquarkV): v.Boost(CMboostV)

    # Bernreuther (k, r, n) helicity basis: k = top-quark direction, r/n span the production plane forming RH basis
    beam = ROOT.TVector3(0., 0., 1.)                           # beam direction
    k = TopV.Vect().Unit()                                     # top quark direction in CM frame
    cosTS = k.Dot(beam)                                        # cosine of top quark scattering angle; reco counterpart is cos_PosTop_beam
    if (1.0 - cosTS * cosTS) <= 1e-12:                         # protect against numerical instability (top quark parallel to beam)
        return None
    r = (beam - k * cosTS).Unit()                              # transverse direction in the scattering plane
    n = beam.Cross(k).Unit()                                   # normal direction to the scattering plane
    sign_cosTS = 1.0 if cosTS > 0 else -1.0                    # Bose symmetry factor
    k_basis, r_basis, n_basis = k, r*sign_cosTS, n*sign_cosTS  # Bose symmetry corrected Bernreuther basis

    sign_absDeltaY = 1.0 if absDeltaY > 0 else -1.0                  # charge-asymmetry factor
    kStar, rStar = k*sign_absDeltaY, (r*sign_cosTS) * sign_absDeltaY # charge-asymmetry-corrected basis

    # prepare rest-frame vectors for the lepton and hadronic b-quark
    LepV_RF = ROOT.TLorentzVector(LepV)
    hadBquarkV_RF = ROOT.TLorentzVector(hadBquarkV)

    # Boost Spin Analyzers to their parent top/antitop rest frame based on lepton charge
    if qpos:
        LepV_RF.Boost(-TopV.BoostVector())
        hadBquarkV_RF.Boost(-antiTopV.BoostVector())
    else:
        LepV_RF.Boost(-antiTopV.BoostVector())
        hadBquarkV_RF.Boost(-TopV.BoostVector())

    # prepare Spin Analyzer unit vectors in top/antitop rest frame
    LepV_unit, hadBquarkV_unit = LepV_RF.Vect().Unit(), hadBquarkV_RF.Vect().Unit()
    SA_1, SA_2 = (LepV_unit, hadBquarkV_unit) if qpos else (hadBquarkV_unit, LepV_unit)

    # Spin Analyzer projections onto the Bernreuther basis
    SA1_k = SA_1.Dot(k_basis)
    SA1_r = SA_1.Dot(r_basis)
    SA1_n = SA_1.Dot(n_basis)
    SA1_kStar = SA_1.Dot(kStar)
    SA1_rStar = SA_1.Dot(rStar)

    SA2_k = SA_2.Dot(k_basis)
    SA2_r = SA_2.Dot(r_basis)
    SA2_n = SA_2.Dot(n_basis)
    SA2_kStar = SA_2.Dot(kStar)
    SA2_rStar = SA_2.Dot(rStar)

    # Lepton exclusive projections (observable set to 99 if the lepton is not present in this event)
    if qpos:
        # no neg. charge leptons in this decay
        Lep_k = Lep_r = Lep_n = Lep_kStar = Lep_rStar = NOTp
        antiLep_k = LepV_unit.Dot(k_basis)
        antiLep_r = LepV_unit.Dot(r_basis)
        antiLep_n = LepV_unit.Dot(n_basis)
        antiLep_kStar = LepV_unit.Dot(kStar)
        antiLep_rStar = LepV_unit.Dot(rStar)
    else:
        # no pos. charge leptons in this decay
        antiLep_k = antiLep_r = antiLep_n = antiLep_kStar = antiLep_rStar = NOTp
        Lep_k = LepV_unit.Dot(k_basis)
        Lep_r = LepV_unit.Dot(r_basis)
        Lep_n = LepV_unit.Dot(n_basis)
        Lep_kStar = LepV_unit.Dot(kStar)
        Lep_rStar = LepV_unit.Dot(rStar)

    # Opening angles
    cHel = LepV_unit.Dot(hadBquarkV_unit)
    cHel_P3n = SA1_k * SA2_k + SA1_r * SA2_r - SA1_n * SA2_n # flipped n-component for boosted entanglement witness

    # Baumgart et al. angular variables
    lep_phi = math.atan2(LepV_unit.Dot(r_basis), LepV_unit.Dot(n_basis))
    b_phi = math.atan2(hadBquarkV_unit.Dot(r_basis), hadBquarkV_unit.Dot(n_basis))
    sphi = _wrap(lep_phi + b_phi)
    dphi = _wrap(lep_phi - b_phi) if qpos else _wrap(b_phi - lep_phi)
    hi = pt_hadTop > pT_HADTOP_Thresh   # "high" pT-of-hadronic-top slice
    lo = pt_hadTop < pT_HADTOP_Thresh   # "low"  pT-of-hadronic-top slice

    return {
        "M_tt": Mtt, "beta": beta, "absDeltaY": absDeltaY, "cosThetaStar": cosTS, "pt_hadTop": pt_hadTop, "Lepton_charge": -1.0 if qpos else 1.0,

        "cosTheta1k_antiLep": antiLep_k, "cosTheta1r_antiLep": antiLep_r, "cosTheta1n_antiLep": antiLep_n, "cosTheta1kStar_antiLep": antiLep_kStar, "cosTheta1rStar_antiLep": antiLep_rStar,
        "cosTheta2k_Lep": Lep_k,         "cosTheta2r_Lep": Lep_r,         "cosTheta2n_Lep": Lep_n,         "cosTheta2kStar_Lep": Lep_kStar,         "cosTheta2rStar_Lep": Lep_rStar,

        "cosTheta1k": SA1_k, "cosTheta1r": SA1_r, "cosTheta1n": SA1_n, "cosTheta1kStar": SA1_kStar, "cosTheta1rStar": SA1_rStar,
        "cosTheta2k": SA2_k, "cosTheta2r": SA2_r, "cosTheta2n": SA2_n, "cosTheta2kStar": SA2_kStar, "cosTheta2rStar": SA2_rStar,

        "Cnn": SA1_n * SA2_n, "Cnr": SA1_n * SA2_r, "Cnk": SA1_n * SA2_k,
        "Crn": SA1_r * SA2_n, "Crr": SA1_r * SA2_r, "Crk": SA1_r * SA2_k,
        "Ckn": SA1_k * SA2_n, "Ckr": SA1_k * SA2_r, "Ckk": SA1_k * SA2_k,
        "Crk_plus": SA1_r * SA2_k + SA1_k * SA2_r, "Crk_minus": SA1_r * SA2_k - SA1_k * SA2_r,
        "Cnr_plus": SA1_n * SA2_r + SA1_r * SA2_n, "Cnr_minus": SA1_n * SA2_r - SA1_r * SA2_n,
        "Cnk_plus": SA1_n * SA2_k + SA1_k * SA2_n, "Cnk_minus": SA1_n * SA2_k - SA1_k * SA2_n,

        "cHel": cHel,
        "cHel_Mtt300_400": cHel if Mtt < 400.0 else SENT,
        "cHel_Mtt300_400_betaLT0p9": cHel if (Mtt < 400.0 and beta < 0.9) else SENT,
        "cHel_P3n": cHel_P3n,
        "cHel_P3n_Mtt800_Inf": cHel_P3n if Mtt > 800.0 else SENT,
        "cHel_P3n_Mtt800_Inf_cosThetaLT0p4": cHel_P3n if (Mtt > 800.0 and abs(cosTS) < 0.4) else SENT,

        "Sigma_phi": sphi, "Sigma_phi_low": (sphi if lo else SENT), "Sigma_phi_high": (sphi if hi else SENT),
        "Delta_phi": dphi, "Delta_phi_low": (dphi if lo else SENT), "Delta_phi_high": (dphi if hi else SENT),
    }


def gen_particle_vectors(ttg):
    """
    Lab-frame 4-vectors of the ttbar-system gen particles listed in GEN_PARTICLES.
    Only call when ttg.IsSemiLeptonicDecay() is True -- ChargedLepton()/BHad()/BLep() are only well defined for a semileptonic decay
    Returns a dict of ROOT.TLorentzVector keyed by the names in GEN_PARTICLES.
    """
    return {
        "Top":     TLvec(ttg.Top()),
        "Antitop": TLvec(ttg.Antitop()),
        "bHad":    TLvec(ttg.BHad()),
        "Lepton":  TLvec(ttg.ChargedLepton()),
    }


# ============================================================================================== #
#                                      Book the output tree                                      #
# ============================================================================================== #
outFile = ROOT.TFile(outF, "RECREATE")
outFile.SetCompressionLevel(5)
tree = ROOT.TTree("SC_SkimTree", "reco+gen Spin Correlation per-event Skim in UL18")

# Event-weight branches: the nominal reco/gen weights
recoWeight_buffer = array('f', [0.0])
tree.Branch("eventweight", recoWeight_buffer, "eventweight/F")
genWeight_buffer = array('f', [0.0])
tree.Branch("genweight", genWeight_buffer, "genweight/F")

# Integer flags: indicates properly defined Observables and  gen-level ttbar decay channel
flag_buffer = {n: array('i', [0]) for n in ("gen_semilep", "gen_channel")}
for n in flag_buffer:
    tree.Branch(n, flag_buffer[n], n + "/I")

# Reco branches: verbatim copies of the RECO observables above with a reco_ prefix
Reco_buffer = {n: array('f', [0.0]) for n in RECO}
for n in RECO:
    tree.Branch("reco_" + n, Reco_buffer[n], "reco_" + n + "/F")

# Gen branches: the recomputed GEN observables above, written with a gen_ prefix
Gen_buffer = {n: array('f', [0.0]) for n in GEN}
for n in GEN:
    tree.Branch("gen_" + n, Gen_buffer[n], "gen_" + n + "/F")

# Gen-particle 4-vector branches: pt/eta/phi/E for each name in GEN_PARTICLES
GenParticle_buffer = {}
for pname in GEN_PARTICLES:
    for comp in ("pt", "eta", "phi", "E"):
        key = pname + "_" + comp
        GenParticle_buffer[key] = array('f', [0.0])
        tree.Branch("gen_" + key, GenParticle_buffer[key], "gen_" + key + "/F")

# Systematic-correction weight branches, copied through verbatim under their original names
SystWeight_buffer = {n: array('f', [0.0]) for n in WEIGHTS}
for n in WEIGHTS:
    tree.Branch(n, SystWeight_buffer[n], n + "/F")


# ============================================================================================================ #
#                                Loop over input files and fill the output tree                                #
# ============================================================================================================ #
for inFile in allFILES:
    name = inFile.split("MC.")[-1].replace(".root", "")
    rootFile = ROOT.TFile.Open(inFile)
    # Check that the file opened successfully and is not a zombie (corrupted)
    if (not rootFile) or rootFile.IsZombie():
        print("SKIP (cannot open):", name)
        continue
    
    # Grab the AnalysisTree from the input file
    AnalysisTree = rootFile.Get("AnalysisTree")

    # Only enable relevant branches
    AnalysisTree.SetBranchStatus("*", 0)
    AnalysisTree.SetBranchStatus("GenParticles*", 1)      # Generator Particle collection
    AnalysisTree.SetBranchStatus("eventweight", 1)        # nominal reconstruction weight
    AnalysisTree.SetBranchStatus("*m_originalXWGTUP*", 1) # nominal generator weight
    for recoBranch in RECO:
        AnalysisTree.SetBranchStatus(recoBranch, 1)       # Reconstruction level spin-correlation branches

    # Systematic weights: enable only those actually present in THIS lepton flavor file
    allbr = set(b.GetName() for b in AnalysisTree.GetListOfBranches())
    present_w = set(w for w in WEIGHTS if w in allbr)
    for w in present_w: AnalysisTree.SetBranchStatus(w, 1)

    # report exactly which weights were missing
    missing_w = [w for w in WEIGHTS if w not in present_w]
    if missing_w:
        print("  NOTE %d/%d weights absent here (written as 1.0): %s%s"
              % (len(missing_w), len(WEIGHTS), ", ".join(missing_w[:5]), " ..." if len(missing_w) > 5 else ""))

    # Get nominal generator-weight 
    # genWeight = AnalysisTree.GetLeaf("genInfo.m_originalXWGTUP") # GetLeaf() does not resolve "branch.leaf" in ROOT 6.14
    genWeight = AnalysisTree.GetLeaf("m_originalXWGTUP")
    if not genWeight:
        print("  WARNING: m_originalXWGTUP not found -> genWeight set to 1.0")

    # Print number of entries in AnalysisTree
    N = AnalysisTree.GetEntries()
    N = N if NMAX < 0 else min(N, NMAX)
    print("Starting loop over %d entries in %s" % (N, name))

    # =================================================================================== #
    #                                   main event loop                                   #
    # =================================================================================== #
    for i in xrange(N):
        # Grab the event
        AnalysisTree.GetEntry(i)

        # Event Weights
        recoWeight = AnalysisTree.eventweight
        recoWeight_buffer[0] = recoWeight
        genWeight_buffer[0] = genWeight.GetValue() if genWeight else 1.0

        # Reco branches: fill directly from AnalysisTree observables
        for recoBranch in RECO:
            Reco_buffer[recoBranch][0] = getattr(AnalysisTree, recoBranch)

        # Systematic weights: fill for those present in this file, default the rest to 1.0
        for n in WEIGHTS:
            SystWeight_buffer[n][0] = getattr(AnalysisTree, n) if n in present_w else 1.0

        # Gen branches: rebuild the ttbar decay from GenParticles and fill with recomputed observables
        ttg = TTbarGenPy(AnalysisTree.GenParticles, False)
        flag_buffer["gen_channel"][0] = ttg.DecayChannel()
        genObs = gen_observables(ttg) if ttg.IsSemiLeptonicDecay() else None
        if genObs is not None:
            flag_buffer["gen_semilep"][0] = 1
            for genBranch in GEN:
                Gen_buffer[genBranch][0] = genObs[genBranch]
        else:
            flag_buffer["gen_semilep"][0] = 0
            for genBranch in GEN:
                Gen_buffer[genBranch][0] = SENT
            print("  WARNING: event %d/%d is not semileptonic or genObs was None -> gen branches set to SENT" % (i, N))

        # Gen-particle 4-vectors: filled whenever the decay is semileptonic
        if ttg.IsSemiLeptonicDecay():
            genParticles = gen_particle_vectors(ttg)
            for pname, vec in genParticles.items():
                GenParticle_buffer[pname + "_pt"][0]  = vec.Pt()
                GenParticle_buffer[pname + "_eta"][0] = vec.Eta()
                GenParticle_buffer[pname + "_phi"][0] = vec.Phi()
                GenParticle_buffer[pname + "_E"][0]   = vec.E()
        else:
            for pname in GEN_PARTICLES:
                GenParticle_buffer[pname + "_pt"][0]  = SENT
                GenParticle_buffer[pname + "_eta"][0] = SENT
                GenParticle_buffer[pname + "_phi"][0] = SENT
                GenParticle_buffer[pname + "_E"][0]   = SENT
            print("  WARNING: event %d/%d is not semileptonic, is %d -> gen particle 4-vectors set to SENT" % (i, N, ttg.DecayChannel()))
            # decay channel codes: hadronic = {0}, semileptonic = {1, 2, 3} dileptonic = {4, 5, 6, 7, 8, 9}

        # Fill output tree with all buffers
        tree.Fill()

        # Periodic progress printout
        if i % 100000 == 0 and i:
            print("  Processed entry %d/%d of %-12s" % (i, N, name))

    # Close the input file
    rootFile.Close()

outFile.cd()    # move to outFile
tree.Write()    # Write output tree to outFile
outFile.Close() # close outFile