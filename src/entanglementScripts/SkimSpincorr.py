# skim_spincorr.py  --  flat per-event reco+gen spin-correlation skim
# ---------------------------------------------------------------------------
# Reads the heavy UL18 TTToSemiLeptonic AnalysisTree files, copies the (already-saved) reco spin-correlation branches verbatim, recomputes the
# gen-level spin-correlation variables from GenParticles via TTbarGenPy, and writes ONE light flat TTree ("entSkim") with reco_* and gen_* branches.
#
# Goal: a small file (no jet/AK8/HOTVR/GenParticle collections) from which any spin-correlation histogram can be remade, instead of re-reading the 230 GB
# analysis output every time.
#
# Scope:
#   - every tree entry is written (sentinels + flags let you filter offline)
#   - TTToSemiLeptonic only (the gen ttbar recipe is exact here)
#   - reco side = the reco spin-correlation branch list, copied through verbatim
#   - gen side  = the SAME variables, recomputed at gen level (same definitions, same sentinel conventions) so reco_X and gen_X are directly comparable
#
# Sentinel conventions (mirror ZprimeSemiLeptonicModules.cxx exactly):
#   -10  : variable not available this event (reco: not reconstructed/b-tagged ;
#          gen: not semileptonic) OR a masked slice whose window is not satisfied
#    99  : a lepton-only projection slot that the lepton is NOT in this event
#          (e.g. cosTheta2k_Lep on an l+ event). Cut both out with abs(x) <= 1.
#
# CMSSW_10_6_28 / Python 2.7 / ROOT 6.14  (run inside the EL7 container)
# ---------------------------------------------------------------------------
from __future__ import print_function
import ROOT
import math
import sys
from array import array
from ttbargen_py import TTbarGenPy

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
# py2/py3 shim: Define xrange regardless of python version.
try:
    xrange
except NameError:
    xrange = range

# ----------------------------- config --------------------------------------
base  = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/"
TAGS  = ["electron", "electron2", "muon", "muon2"]
ALL_FILES = [base + "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root" % s for s in TAGS]
 
# Per-file parallelism: pass a file index 0-3 (condor $(Process)) to skim ONE file into its
# own output; with NO argument, process all four into the single merged-name file.
if len(sys.argv) > 1:
    idx   = int(sys.argv[1])
    files = [ALL_FILES[idx]]
    out   = "spincorr_skim_TTToSemiLeptonic_UL18_%s.root" % TAGS[idx]
else:
    files = ALL_FILES
    out   = "spincorr_skim_TTToSemiLeptonic_UL18.root"
 
NMAX  = -1          # 200000 for a quick test; -1 = full pass
MULT  = 5.0         # coefficient multiplier for the self-validation (cHel / cHel_P3n)
SENT  = -10.0       # sentinel: a side/window that does not exist this event
NOTF  = 99.0        # sentinel: lepton-only slot the lepton is not in (matches C++ 99.)
PT_HADTOP_THRESH = 150.0   # ZprimeSemiLeptonicModules.cxx:1517 (Baumgart low/high split)
# ---------------------------------------------------------------------------

# Reco spin-correlation branches (ZprimeSemiLeptonicModules.cxx declare_event_output,
# lines 1339-1416). Copied through verbatim with a reco_ prefix.
RECO = [# variables used to filter events with
        "chi2", "M_tt", "beta", "dyreco", # need to save 'cos_PosTop_beam' and 'pt_hadTop' next time I process analysis
        # lep-only projections and *Star asymmetry versions
        "cosTheta1k_antiLep", "cosTheta1r_antiLep", "cosTheta1n_antiLep", "cosTheta1kStar_antiLep", "cosTheta1rStar_antiLep",
        "cosTheta2k_Lep", "cosTheta2r_Lep", "cosTheta2n_Lep", "cosTheta2kStar_Lep", "cosTheta2rStar_Lep",
        # spin-analyzer 1/2 projections and *Star asymmetry versions
        "cosTheta1k", "cosTheta1r", "cosTheta1n", "cosTheta1kStar", "cosTheta1rStar",
        "cosTheta2k", "cosTheta2r", "cosTheta2n", "cosTheta2kStar", "cosTheta2rStar",
        # Correlatin matrix elements and plus/minus combinations
        "Cnn", "Cnr", "Cnk", "Crn", "Crr", "Crk", "Ckn", "Ckr", "Ckk",
        "Crk_plus", "Crk_minus", "Cnr_plus", "Cnr_minus", "Cnk_plus", "Cnk_minus",
        # Entanglement witnesses and differential slices
        "cHel", "cHel_Mtt300_400", "cHel_Mtt300_400_betaLT0p9",
        "cHel_P3n", "cHel_P3n_Mtt800_Inf", "cHel_P3n_Mtt800_Inf_cosThetaLT0p4",
        # Baumgart angular variables
        "Sigma_phi", "Sigma_phi_low", "Sigma_phi_high",
        "Delta_phi", "Delta_phi_low", "Delta_phi_high",
        # need to save ttbar system kinematics from BestChi2Candidate next time I process analysis
        ]

# gen counterparts, recomputed with the SAME definitions (validated by component_diag.py).
# Parallel to RECO except: chi2 has no gen analog (dropped); cosThetaStar & pt_hadTop are
# gen-available now (reco versions await the next reprocessing).
GEN = [# ttbar-system kinematics + charge asymmetry (gen analogs of the reco "filter" block)
       "M_tt", "beta", "dyreco", "cosThetaStar", "pt_hadTop",
       # lep-only projections and *Star (99 on the slot the lepton is not in, like reco)
       "cosTheta1k_antiLep", "cosTheta1r_antiLep", "cosTheta1n_antiLep", "cosTheta1kStar_antiLep", "cosTheta1rStar_antiLep",
       "cosTheta2k_Lep", "cosTheta2r_Lep", "cosTheta2n_Lep", "cosTheta2kStar_Lep", "cosTheta2rStar_Lep",
       # spin-analyzer 1/2 projections and *Star asymmetry versions
       "cosTheta1k", "cosTheta1r", "cosTheta1n", "cosTheta1kStar", "cosTheta1rStar",
       "cosTheta2k", "cosTheta2r", "cosTheta2n", "cosTheta2kStar", "cosTheta2rStar",
       # Correlation matrix elements and plus/minus combinations
       "Cnn", "Cnr", "Cnk", "Crn", "Crr", "Crk", "Ckn", "Ckr", "Ckk",
       "Crk_plus", "Crk_minus", "Cnr_plus", "Cnr_minus", "Cnk_plus", "Cnk_minus",
       # Entanglement witnesses and differential slices
       "cHel", "cHel_Mtt300_400", "cHel_Mtt300_400_betaLT0p9",
       "cHel_P3n", "cHel_P3n_Mtt800_Inf", "cHel_P3n_Mtt800_Inf_cosThetaLT0p4",
       # Baumgart angular variables
       "Sigma_phi", "Sigma_phi_low", "Sigma_phi_high",
       "Delta_phi", "Delta_phi_low", "Delta_phi_high",
       ]

# Systematic-correction weights (per-event scalars in the AnalysisTree), copied through
# verbatim (no prefix) so any reco_ or gen_ histogram can be reweighted for systematics.
# NOTE: the analysis runs per lepton flavor, so weight_sfelec_* and weight_sfmu_* may each
# be present in only one flavor's files. A weight absent from a given file is written as 1.0
# (neutral) and reported per file -- check that print to confirm what each file carries.
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
WEIGHTS = [w for w in WEIGHTS_RAW.split() if w]   # split() drops all whitespace/newlines


def tlv(genp):
    q = genp.v4(); v = ROOT.TLorentzVector()
    v.SetPtEtaPhiE(q.pt(), q.eta(), q.phi(), q.energy()); return v


def _wrap(x):
    # map an angle back into (-pi, pi], one correction (mirrors C++ 1846-1849)
    if x > math.pi:  x -= 2.0 * math.pi
    if x < -math.pi: x += 2.0 * math.pi
    return x


def gen_observables(ttg):
    """Gen-level spin-correlation variables, recomputed with the SAME definitions and
    sentinel conventions as ZprimeSemiLeptonicModules.cxx (1628-1893). Returns a dict
    keyed exactly by GEN, or None if the production geometry is degenerate (top || beam)."""
    lep = ttg.ChargedLepton()
    qpos = (lep.pdgId() < 0)                       # anti-lepton (l+, charge>0) <=> lepton from the top
    Top  = tlv(ttg.Top());  Anti = tlv(ttg.Antitop())
    lepv = tlv(lep);        bv   = tlv(ttg.BHad())
    ttbar = Top + Anti
    Mtt = ttbar.M()

    # --- lab-frame quantities (computed BEFORE the CoM boost) ---
    eS = Top.E() + Anti.E()
    beta = abs(Top.Pz() + Anti.Pz()) / eS if eS > 0 else SENT     # cxx:1629
    dyreco = abs(Top.Rapidity()) - abs(Anti.Rapidity())          # |y_t| - |y_tbar|  (cxx:1632-1636)
    pt_hadTop = tlv(ttg.TopHad()).Pt()                            # cxx:1518 (lab pt of hadronic top)

    # --- boost tops + analyzers into the ttbar CoM ---
    bcm = -ttbar.BoostVector()
    for v in (Top, Anti, lepv, bv): v.Boost(bcm)
    beam = ROOT.TVector3(0., 0., 1.)
    k = Top.Vect().Unit(); cosTS = k.Dot(beam)                   # cos_PosTop_beam (cxx:1656)
    if (1.0 - cosTS * cosTS) <= 1e-12:
        return None
    r = (beam - k * cosTS).Unit(); n = beam.Cross(k).Unit()
    sgn_cos = 1.0 if cosTS > 0 else -1.0                          # Bose factor      (cxx:1660)
    sgn_rap = 1.0 if dyreco > 0 else -1.0                         # charge-asym factor (cxx:1661)
    kbase, rbase, nbase = k, r * sgn_cos, n * sgn_cos            # cxx:1664-1666
    kStar, rStar = k * sgn_rap, (r * sgn_cos) * sgn_rap          # cxx:1668-1669 (no nStar)

    # --- analyzers into their parent-top rest frames ---
    lr = ROOT.TLorentzVector(lepv); br = ROOT.TLorentzVector(bv)
    if qpos:                                                      # lepton from Top, b from Antitop
        lr.Boost(-Top.BoostVector());  br.Boost(-Anti.BoostVector())
    else:                                                        # lepton from Antitop, b from Top
        lr.Boost(-Anti.BoostVector()); br.Boost(-Top.BoostVector())
    lu = lr.Vect().Unit(); bu = br.Vect().Unit()

    # --- charge-ordered spin analyzers: slot 1 = top-side, slot 2 = antitop-side (cxx:1748-1775) ---
    a1, a2 = (lu, bu) if qpos else (bu, lu)
    c1k, c1r, c1n = a1.Dot(kbase), a1.Dot(rbase), a1.Dot(nbase)
    c2k, c2r, c2n = a2.Dot(kbase), a2.Dot(rbase), a2.Dot(nbase)
    c1kS, c1rS = a1.Dot(kStar), a1.Dot(rStar)
    c2kS, c2rS = a2.Dot(kStar), a2.Dot(rStar)

    # --- lepton-only projections: lepton -> slot 1 (antiLep) if l+, slot 2 (Lep) if l-;
    #     the other slot is 99, exactly like the C++ (cxx:1717-1743) ---
    if qpos:
        L1k, L1r, L1n, L1kS, L1rS = lu.Dot(kbase), lu.Dot(rbase), lu.Dot(nbase), lu.Dot(kStar), lu.Dot(rStar)
        L2k = L2r = L2n = L2kS = L2rS = NOTF
    else:
        L1k = L1r = L1n = L1kS = L1rS = NOTF
        L2k, L2r, L2n, L2kS, L2rS = lu.Dot(kbase), lu.Dot(rbase), lu.Dot(nbase), lu.Dot(kStar), lu.Dot(rStar)

    # --- entanglement witnesses (cxx:1812-1813) ---
    cHel = lu.Dot(bu)
    cHel_P3n = c1k * c2k + c1r * c2r - c1n * c2n

    # --- Baumgart phi (cxx:1832-1849) ---
    lep_phi = math.atan2(lu.Dot(rbase), lu.Dot(nbase))
    b_phi   = math.atan2(bu.Dot(rbase), bu.Dot(nbase))
    sphi = _wrap(lep_phi + b_phi)
    dphi = _wrap(lep_phi - b_phi) if qpos else _wrap(b_phi - lep_phi)
    hi = pt_hadTop > PT_HADTOP_THRESH
    lo = pt_hadTop < PT_HADTOP_THRESH

    return {
        "M_tt": Mtt, "beta": beta, "dyreco": dyreco, "cosThetaStar": cosTS, "pt_hadTop": pt_hadTop,
        "cosTheta1k_antiLep": L1k, "cosTheta1r_antiLep": L1r, "cosTheta1n_antiLep": L1n,
        "cosTheta1kStar_antiLep": L1kS, "cosTheta1rStar_antiLep": L1rS,
        "cosTheta2k_Lep": L2k, "cosTheta2r_Lep": L2r, "cosTheta2n_Lep": L2n,
        "cosTheta2kStar_Lep": L2kS, "cosTheta2rStar_Lep": L2rS,
        "cosTheta1k": c1k, "cosTheta1r": c1r, "cosTheta1n": c1n, "cosTheta1kStar": c1kS, "cosTheta1rStar": c1rS,
        "cosTheta2k": c2k, "cosTheta2r": c2r, "cosTheta2n": c2n, "cosTheta2kStar": c2kS, "cosTheta2rStar": c2rS,
        "Cnn": c1n * c2n, "Cnr": c1n * c2r, "Cnk": c1n * c2k,
        "Crn": c1r * c2n, "Crr": c1r * c2r, "Crk": c1r * c2k,
        "Ckn": c1k * c2n, "Ckr": c1k * c2r, "Ckk": c1k * c2k,
        "Crk_plus": c1r * c2k + c1k * c2r, "Crk_minus": c1r * c2k - c1k * c2r,
        "Cnr_plus": c1n * c2r + c1r * c2n, "Cnr_minus": c1n * c2r - c1r * c2n,
        "Cnk_plus": c1n * c2k + c1k * c2n, "Cnk_minus": c1n * c2k - c1k * c2n,
        "cHel": cHel,
        "cHel_Mtt300_400":          cHel if Mtt < 400.0 else SENT,
        "cHel_Mtt300_400_betaLT0p9": cHel if (Mtt < 400.0 and beta < 0.9) else SENT,
        "cHel_P3n": cHel_P3n,
        "cHel_P3n_Mtt800_Inf":               cHel_P3n if Mtt > 800.0 else SENT,
        "cHel_P3n_Mtt800_Inf_cosThetaLT0p4": cHel_P3n if (Mtt > 800.0 and abs(cosTS) < 0.4) else SENT,
        "Sigma_phi": sphi, "Sigma_phi_low": (sphi if lo else SENT), "Sigma_phi_high": (sphi if hi else SENT),
        "Delta_phi": dphi, "Delta_phi_low": (dphi if lo else SENT), "Delta_phi_high": (dphi if hi else SENT),
    }


# ----------------------------- output tree ---------------------------------
fout = ROOT.TFile(out, "RECREATE")
fout.SetCompressionLevel(5)
tree = ROOT.TTree("entSkim", "reco+gen spin-correlation per-event skim (TTToSemiLeptonic UL18)")

wbuf = array('f', [0.0]); tree.Branch("eventweight", wbuf, "eventweight/F")
gwbuf = array('f', [0.0]); tree.Branch("genweight",  gwbuf, "genweight/F")   # nominal generator weight -> use for gen, NEVER eventweight
fbuf = {n: array('i', [0]) for n in ("pass_reco", "gen_semilep", "gen_channel")}

for n in fbuf:
    tree.Branch(n, fbuf[n], n + "/I")
rbuf = {n: array('f', [0.0]) for n in RECO}

for n in RECO:
    tree.Branch("reco_" + n, rbuf[n], "reco_" + n + "/F")
gbuf = {n: array('f', [0.0]) for n in GEN}

for n in GEN:
    tree.Branch("gen_" + n, gbuf[n], "gen_" + n + "/F")
wtbuf = {n: array('f', [0.0]) for n in WEIGHTS}      # systematic-correction weights (verbatim names)

for n in WEIGHTS:
    tree.Branch(n, wtbuf[n], n + "/F")

# ----------------------------- validation accumulators ---------------------
iW = irS = 0.0          # inclusive reco D-tilde  (expect 0.3519)
r8W = r8rS = 0.0        # reco D-tilde, reco M_tt>800  (expect ~0.1537)
g8rW = g8rS = 0.0       # boosted reco (eventweight), gen M_tt>800
g8gW = g8gS = 0.0       # boosted gen  (genweight),   gen M_tt>800
n_written = 0

# ----------------------------- loop ----------------------------------------
for fn in files:
    name = fn.split("MC.")[-1].replace(".root", "")
    f = ROOT.TFile.Open(fn)
    if (not f) or f.IsZombie():
        print("SKIP (cannot open):", name); continue
    t = f.Get("AnalysisTree")
    t.SetBranchStatus("*", 0)
    t.SetBranchStatus("GenParticles*", 1)
    t.SetBranchStatus("eventweight", 1)
    t.SetBranchStatus("*m_originalXWGTUP*", 1)   # nominal generator weight (cheap leaf; NEVER the 16 GB m_systweights)
    for n in RECO:
        t.SetBranchStatus(n, 1)
    # systematic weights: enable only those present in THIS file (flavor sets differ)
    allbr = set(b.GetName() for b in t.GetListOfBranches())
    present_w = set(w for w in WEIGHTS if w in allbr)
    for w in present_w:
        t.SetBranchStatus(w, 1)
    missing_w = [w for w in WEIGHTS if w not in present_w]
    if missing_w:
        print("  NOTE %d/%d weights absent here (written as 1.0): %s%s"
              % (len(missing_w), len(WEIGHTS), ", ".join(missing_w[:5]),
                 " ..." if len(missing_w) > 5 else ""))
    gw_leaf = t.GetLeaf("genInfo.m_originalXWGTUP")
    if not gw_leaf:
        gw_leaf = t.GetLeaf("m_originalXWGTUP")
    if not gw_leaf:
        print("  WARNING: genInfo.m_originalXWGTUP not found -> genweight set to 1.0")
    N = t.GetEntries(); N = N if NMAX < 0 else min(N, NMAX)
    print("%-30s entries=%d" % (name, N))
    for i in xrange(N):
        t.GetEntry(i)
        w = t.eventweight
        wbuf[0] = w
        gwbuf[0] = gw_leaf.GetValue() if gw_leaf else 1.0

        # ---- reco side: verbatim copy ----
        for n in RECO:
            rbuf[n][0] = getattr(t, n)
        rcp = rbuf["cHel_P3n"][0]
        pass_reco = 1 if rcp > -1.5 else 0
        fbuf["pass_reco"][0] = pass_reco

        # ---- systematic weights: copy present ones, default absent to 1.0 ----
        for n in WEIGHTS:
            wtbuf[n][0] = getattr(t, n) if n in present_w else 1.0

        # ---- gen side: rebuild from GenParticles ----
        ttg = TTbarGenPy(t.GenParticles, False)
        fbuf["gen_channel"][0] = ttg.DecayChannel()
        g = gen_observables(ttg) if ttg.IsSemiLeptonicDecay() else None
        if g is not None:
            fbuf["gen_semilep"][0] = 1
            for n in GEN:
                gbuf[n][0] = g[n]
        else:
            fbuf["gen_semilep"][0] = 0
            for n in GEN:
                gbuf[n][0] = SENT

        tree.Fill(); n_written += 1

        # ---- running validation (mirrors dilution_map.py) ----
        if pass_reco:
            rs = 1.0 if rcp >= 0 else -1.0
            iW += w; irS += w * rs
            if rbuf["M_tt"][0] > 800.0:
                r8W += w; r8rS += w * rs
            if g is not None and g["M_tt"] > 800.0:
                gw = gwbuf[0]
                gs = 1.0 if g["cHel_P3n"] >= 0 else -1.0
                g8rW += w;  g8rS += w * rs       # reco: eventweight
                g8gW += gw; g8gS += gw * gs      # gen:  genweight

        if i % 500000 == 0 and i:
            print("  %-12s %d/%d" % (name, i, N))
    f.Close()

fout.cd(); tree.Write(); fout.Close()

# ----------------------------- report --------------------------------------
D = lambda s, wv: MULT * s / wv if wv else 0.0
print("\nwrote %s : %d events" % (out, n_written))
print("\n--- self-validation (should match the locked baseline) ---")
print("  reco Dtilde inclusive    = %.4f   (analysis cHel_P3n        = 0.3519)" % D(irS, iW))
print("  reco Dtilde reco_Mtt>800 = %.4f   (analysis cHel_P3n_Mtt800 = 0.1539)" % D(r8rS, r8W))
rr, gg = D(g8rS, g8rW), D(g8gS, g8gW)
print("  boosted (gen_Mtt>800): reco=%.4f (eventwt)  gen=%.4f (genwt)  reco/gen=%.2f"
      % (rr, gg, (rr / gg if gg else 0.0)))
print("    [gen now uses genweight; gen value & ratio shift slightly vs the old eventweighted 0.2195/1.26]")
