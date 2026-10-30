"""
Shared configuration of the three plotting scripts: input locations, sample grouping and colours, labels, and the
lists of distributions. The datacard processes and the nuisance definitions are imported from
systematics/buildCombineInput.py, so the plots always use the same samples and nuisance names as the Combine input.
"""
import fnmatch
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "systematics"))
import buildCombineInput as combine  # noqa: E402  (PROCESSES, DATA_SAMPLES, EXPERIMENTAL, THEORY_GROUP, ...)


# =================================================================================================================== #
#                                                   Input / output                                                    #
# =================================================================================================================== #
ANALYSIS_DIR     = "/data/dust/group/cms/zprime-uhh/Analysis_SC/{year}/{lep}/workdir_Analysis_{year}_{lep}/"
PRESELECTION_DIR = "/data/dust/group/cms/zprime-uhh/Presel_SC/{year}/workdir_Preselection_{year}/"
OUTPUT_DIR       = "plots/{year}"   # relative to where the script is run, override with --outdir

CHANNELS = ["muon", "electron"]


# =================================================================================================================== #
#                                                       Labels                                                        #
# =================================================================================================================== #
# --year -> iYear argument of CMS_lumi (luminosity values are the ones in CMS_lumi.py)
LUMI_PERIOD = {"UL16": 16, "UL16preVFP": 16, "UL16postVFP": 16, "UL17": 17, "UL18": 18, "Run2": 2}

CHANNEL_LABEL = {"muon": "#mu+jets", "electron": "e+jets"}
REGION_LABEL  = {"Chi2cut_PASS": "#chi^{2} selection", "Merged": "Merged"}

EXTRA_TEXT_DATA = "Private Work"
EXTRA_TEXT_SIM  = "Simulation Private Work"


# =================================================================================================================== #
#                                       Stacked processes for the control plots                                       #
# =================================================================================================================== #
# (label, colour, sample-name patterns), drawn bottom -> top of the stack; legend in reverse order.
# Colours: Using custom colour palette.
STACK_GROUPS = [
    ("QCD",                          "#5CC9CA", ["QCD_*"]),            # kCyan+1 HEX equivalent
    ("Diboson",                      "#EA33F7", ["WW", "WZ", "ZZ"]),   # kMagenta HEX equivalent
    ("Drell-Yan",                    "#F7CE46", ["DYJetsToLL_*"]),     # kOrange HEX equivalent
    ("W+jets",                       "#75FB4C", ["WJetsToLNu_*"]),     # kGreen HEX equivalent
    ("Single t",                     "#0000F5", ["ST_*"]),             # kBlue HEX equivalent
    ("t#bar{t} (had.)",              "#ED706B", ["TTToHadronic"]),     # kRed-7 HEX equivalent
    ("t#bar{t} (2#font[12]{l})",     "#8C1A11", ["TTTo2L2Nu"]),        # kRed+2 HEX equivalent
    ("t#bar{t} (#font[12]{l}+jets)", "#EA3323", ["TTToSemiLeptonic"]), # kRed HEX equivalent
]


def stack_group(version):
    """index of the STACK_GROUPS entry of a sample"""
    for i, (_, _, patterns) in enumerate(STACK_GROUPS):
        if any(fnmatch.fnmatch(version, p) for p in patterns):
            return i
    raise RuntimeError("sample %s is not in any STACK_GROUPS entry" % version)


def cells():
    """
    Every MC sample belongs to one 'cell' = (stack group, datacard process); histograms are summed per cell, which is
    enough to build both the stacked groups and the per-process systematic variations.
    returns {cell_name: [sample versions]}, cell_name = "g<stack group index>__<process>" (also used in the cache file)
    """
    out = {}
    for proc, versions in combine.PROCESSES.items():
        for v in versions:
            out.setdefault("g%d__%s" % (stack_group(v), proc), []).append(v)
    return out


def cell_group(cell):
    """stack-group label of a cell"""
    return STACK_GROUPS[int(cell.split("__")[0][1:])][0]


def cell_process(cell):
    return cell.split("__")[1]


# =================================================================================================================== #
#                                                Distributions to plot                                                #
# =================================================================================================================== #
def binning(name, source):
    """
    (nbins, lo, hi) for a distribution, or None to keep the binning it was booked with.
    source: "preselection", "analysis" (the _General folders) or "sc" (the <region>_SC_* folders of the systematics plots).

    The histograms are already binned, so a binning here can only merge booked bins and restrict the range: every new bin edge has to coincide with a booked edge. 
    Otherwise the plot is skipped with a message that states the booked binning. 
    
    Find original binning for Preselection hists in ZprimeSemiLeptonicPreselectionHists.cxx
    Find original binning for Analysis Selection hists in ZprimeSemiLeptonicHists.cxx

    Tune the entries below while inspecting the plots, e.g.
        if name.startswith("pt_mu"):  return (30, 0., 900.)
    """
    if source == "preselection":
        # Common binning
        if name.startswith("pt_jet"):  return (17,    0,  510) # (50, 0, 1500) binwidth=30 same for jet & jet1
        if name.startswith("pt_ele"):  return ( 6,    0,  180) # (50, 0, 1500) binwidth=30 same for ele & ele1
        if name.startswith("eta_"): return (25, -2.5,  2.5) # (50, -2.5, 2.5)  binwidth=0.1 same for eta_* where *={jet, jet1, mu, mu1, ele, ele1, AK8Puppijet1}
        if name.startswith("phi_"): return ( 7, -3.5,  3.5) # (35, -3.5, 3.5)  binwidth=0.2 same for phi_* where *={jet, jet1, mu, mu1, ele, ele1, AK8Puppijet1}
        if name.startswith("m_jet"): return (10,    0,  100)   # (50, 0,  500) binwidth=10 same for jet & jet1
        if name.startswith("dRmin_"): return (20, 0,   1) # (60, 0, 3) binwidth=0.05 for dRmin_* where *={mu1_jet, ele1_jet}
        if name.startswith("ptrel_"): return (20, 0, 200) # (50, 0, 500) binwidth=10 for ptrel_* where *={mu1_jet, ele1_jet}

        # Custom binning
        if name == "N_jets":  return (21, -0.5, 20.5) # binwidth=1
        if name == "N_mu":    return ( 7, -0.5,  6.5) # (11, -0.5, 10.5) binwidth=1
        if name == "N_ele":   return ( 7, -0.5,  6.5) # (11, -0.5, 10.5) binwidth=1
        if name == "N_AK8Puppijets":          return ( 6, -0.5,  5.5) # binwidth=1
        if name == "N_daughters_AK8Puppijet": return (51, -0.5, 50.5) # binwidth=1
        if name == "N_subjets_AK8Puppijet1":  return ( 6, -0.5,  5.5) # binwidth=1 need to run with --no-overflow or extend range in UHH2

        if name == "pt_mu":           return (17, 0, 510) # (50, 0, 1500) binwidth=30
        if name == "pt_mu1":          return (50, 0, 180) # (50, 0, 1500) binwidth=30
        if name == "pt_AK8Puppijet1": return (24, 0, 720) # (50, 0, 1500) binwidth=30

        if name == "mSD_AK8Puppijet1": return (20,    0,  200) # (50, 0, 500) binwidth=10

        if name == "MET": return (10, 0, 1400) # (50, 0, 7000) binwidth=140 need to run with logy
        if name == "ST":  return (10, 0, 1400) # (50, 0, 7000) binwidth=140


    if source == "analysis":
        # Common binning
        if name.startswith("N_"): return None # default okay for all N_*
        if name.startswith("pt_jet"): return (45, 0,  900) # binWidth=20 same for jet & jet1
        if name.startswith("pt_mu"):  return (30, 0,  300) # (90, 0,  900) binWidth=10 same for mu & mu1
        if name.startswith("pt_ele"): return (30, 0,  300) # (90, 0,  900) binWidth=10 same for ele & ele1
        if name.startswith("pt_AK8Puppi"): return None # default (45, 0,  900) binWidth=20 okay for AK8Puppijet1 & AK8PuppiTaggedjet1
        if name.startswith("eta_"): return (25, -2.5,  2.5) # (50, -2.5, 2.5)  binwidth=0.1 same for eta_* where *={jet, jet1, mu, mu1, ele, ele1, AK8Puppijet1, AK8PuppiTaggedjet1}
        if name.startswith("phi_"): return ( 7, -3.5,  3.5) # (35, -3.5, 3.5)  binwidth=0.2 same for phi_* where *={jet, jet1, mu, mu1, ele, ele1, AK8Puppijet1, AK8PuppiTaggedjet1}
        if name.startswith("m_jet"):  return (10,   0,  100)   # (50, 0, 500) binwidth=10 same for jet & jet1
        if name.startswith("dRmin_"): return None # (60, 0, 3) binWidth=0.05
        if name.startswith("ptrel_"): return None # (50, 0, 500) binWidth=10
        if name.startswith("massdecordeepjet_"): return None # (20, 0, 1) binWidth=0.05
        if name.startswith("deepjetbscore_"):    return None # (20, 0, 1) binWidth=0.05
        if name.startswith("cos_"): return None # (24, -1, 1) binWidth=0.083
        if name.startswith("C"):    return None # (24, -1, 1) binWidth=0.083

        # Custom binning
        if name =="toplep_pt":  return (10, 0, 1000) # (70, 0, 7000) binWidth=100
        if name =="tophad_pt":  return (10, 0, 1000) # (70, 0, 7000) binWidth=100
        if name =="toplep_eta": return (30, -3.0, 3.0) # (60, -3.0, 3.0) binWidth=0.1
        if name =="tophad_eta": return (30, -3.0, 3.0) # (60, -3.0, 3.0) binWidth=0.1
        if name =="toplep_phi": return (35, -3.5, 3.5) # (70, -3.5, 3.5) binWidth=0.1
        if name =="tophad_phi": return (35, -3.5, 3.5) # (70, -3.5, 3.5) binWidth=0.1
        if name =="toplep_m":   return ( 5, 0, 500) # (70, 0, 7000) binWidth=100
        if name =="tophad_m":   return ( 5, 0, 500) # (70, 0, 7000) binWidth=100

        if name =="mSD_AK8Puppijet1":       return (20, 0, 200) # (50, 0, 500) binWidth=10
        if name =="mSD_AK8PuppiTaggedjet1": return (30, 0, 300) # (50, 0, 500) binWidth=10
        
        if name =="ditop_mass": return (20, 0, 2000) # (70, 0, 7000) binWidth=100
        if name =="M_Zprime":   return (80, 0, 2000) # (280, 0, 7000) binWidth=25

        if name =="MET": return (10, 0, 1400) # (50, 0, 7000) binWidth=140
        if name =="ST":  return (20, 0, 2800) # (50, 0, 7000) binWidth=140

        if name =="beta_ttbar": return (50, 0, 1) # binWidth=0.02

        if name =="deltaR_hadTop_bGen": return (60, 0, 3) # binWidth=0.05

        if name =="Sigma_phi": return None # (16, -3.2, 3.2) binWidth=0.4
        if name =="Delta_phi": return None # (16, -3.2, 3.2) binWidth=0.4

        

    if source == "sc":
        return None          # None -> rebinned by SC_REBIN like the Combine input
    return None


# Under/overflow is added to the first/last bin for open-ended quantities (pt, masses, MET, ST, multiplicities, ...).
# Bounded quantities are never folded: their under/overflow only holds sentinel values (-10 / 99) of events where the
# quantity is not defined. --no-overflow switches the folding off everywhere; OVERFLOW_OVERRIDE sets it per name.
BOUNDED_PATTERNS = ["cos*", "C??", "C??_*", "cHel*", "*eta*", "*phi*", "*Phi*", "Sigma_phi*", "Delta_phi*",
                    "*score*", "*Prob*", "*TvsQCD*", "beta_ttbar", "tau*"]
OVERFLOW_OVERRIDE = {
}


def fold_overflow(name):
    if name in OVERFLOW_OVERRIDE:
        return OVERFLOW_OVERRIDE[name]
    return not any(fnmatch.fnmatchcase(name, p) for p in BOUNDED_PATTERNS)


PRESELECTION_KINEMATICS = [
    "N_jets",
    "pt_jet", "eta_jet", "phi_jet", "m_jet",
    "pt_jet1", "eta_jet1", "phi_jet1", "m_jet1",
    "N_mu", "pt_mu", "eta_mu", "phi_mu", "pt_mu1", "eta_mu1", "phi_mu1",
    "N_ele", "pt_ele", "eta_ele", "phi_ele", "pt_ele1", "eta_ele1", "phi_ele1",
    "dRmin_mu1_jet", "dRmin_ele1_jet", "ptrel_mu1_jet", "ptrel_ele1_jet",
    "N_AK8Puppijets", "pt_AK8Puppijet1", "eta_AK8Puppijet1", "phi_AK8Puppijet1", "mSD_AK8Puppijet1",
    "N_subjets_AK8Puppijet1", "N_daughters_AK8Puppijet",
    "MET", "ST",
]

ANALYSIS_KINEMATICS = [
    "N_jets",
    "pt_jet", "eta_jet", "phi_jet", "m_jet",
    "pt_jet1", "eta_jet1", "phi_jet1", "m_jet1",
    "deepjetbscore_jet1", "N_bJetsDeepJet_med",
    "N_mu", "pt_mu", "eta_mu", "phi_mu", "pt_mu1", "eta_mu1", "phi_mu1",
    "N_ele", "pt_ele", "eta_ele", "phi_ele", "pt_ele1", "eta_ele1", "phi_ele1",
    "dRmin_mu1_jet", "dRmin_ele1_jet", "ptrel_mu1_jet", "ptrel_ele1_jet",
    "N_AK8Puppijets", "pt_AK8Puppijet1", "eta_AK8Puppijet1", "phi_AK8Puppijet1", "mSD_AK8Puppijet1",
    "N_subjets_AK8Puppijet1", "N_daughters_AK8Puppijet", "massdecordeepjet_topProb_jet1",
    "N_AK8PuppiTaggedjets", "pt_AK8PuppiTaggedjet1", "eta_AK8PuppiTaggedjet1", "phi_AK8PuppiTaggedjet1",
    "mSD_AK8PuppiTaggedjet1", "massdecordeepjet_TvsQCD_Taggedjet1",
    "MET", "ST",
    "toplep_pt", "toplep_eta", "toplep_phi", "toplep_m",
    "tophad_pt", "tophad_eta", "tophad_phi", "tophad_m",
    "ditop_mass", "beta_ttbar", "deltaR_hadTop_bGen", "cos_ThetaStar",
    "cos_theta1k_antiLep", "cos_theta1r_antiLep", "cos_theta1n_antiLep", "cos_theta1kStar_antiLep", "cos_theta1rStar_antiLep",
    "cos_theta2k_Lep", "cos_theta2r_Lep", "cos_theta2n_Lep", "cos_theta2kStar_Lep", "cos_theta2rStar_Lep",
    "cos_theta1k", "cos_theta1r", "cos_theta1n", "cos_theta1kStar", "cos_theta1rStar",
    "cos_theta2k", "cos_theta2r", "cos_theta2n", "cos_theta2kStar", "cos_theta2rStar",
    "Cnn", "Cnr", "Cnk", "Crn", "Crr", "Crk", "Ckn", "Ckr", "Ckk",
    "Crk_plus", "Crk_minus", "Cnr_plus", "Cnr_minus", "Cnk_plus", "Cnk_minus",
    "Sigma_phi", "Delta_phi", "M_Zprime",
]
# generator-level quantities: data is never drawn for these
GEN_ONLY = {"deltaR_hadTop_bGen"}

ANALYSIS_REGIONS = ["Chi2cut_PASS", "Merged"]   # regions with <region>_SC_* folders (systematics plots)

# fit observables: data is only drawn for them with --unblind (blinded analysis)
BLINDED = {"cHel_P3n_Mtt800_Inf_cosThetaLT0p4"}

SC_REBIN = 4  # same rebinning as the Combine input (24 -> 6 bins)


# =================================================================================================================== #
#                                                Systematic variations                                                #
# =================================================================================================================== #
def nuisances(channel, year):
    """
    Nuisances as Combine knows them. Each entry:
      name      Combine nuisance name
      up, down  folder suffixes (<region>_SC_<suffix>), or "pdf" for the Hessian PDF built from pdf0..pdf99
      procs     datacard processes it acts on (None = all MC)
    """
    out = []
    for stem in combine.EXPERIMENTAL["common"] + combine.EXPERIMENTAL[channel]:
        name = stem + "_" + year if stem in combine.YEAR_DECORRELATED else stem
        out.append({"name": name, "up": stem + "_up", "down": stem + "_down", "procs": None})
    groups = []
    for proc, grp in combine.THEORY_GROUP.items():
        if grp not in groups:
            groups.append(grp)
    for stem, (up, down) in combine.THEORY_UPDOWN.items():
        for grp in groups:
            procs = [p for p, g in combine.THEORY_GROUP.items() if g == grp]
            out.append({"name": "%s_%s" % (stem, grp), "up": up, "down": down, "procs": procs})
    out.append({"name": "pdf", "up": "pdf", "down": "pdf", "procs": None})
    out.append({"name": "NoSC_f", "up": "NoSC_up", "down": "NoSC_down", "procs": [combine.SIGNAL]})
    return out


N_PDF = 100
