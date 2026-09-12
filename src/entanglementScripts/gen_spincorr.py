from __future__ import print_function
import ROOT
from array import array
from ttbargen_py import TTbarGenPy

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

'''
Input Files: TTToSemiLeptonic (electron and muon channels) MC analysis output
This script computes the threshold & boosted entanglement proxies (D & Dtilde) for semi-leptonic ttbar decays at generator level
1D histograms of spin analyzer's opening angle (cos(phi_lb)) and reflected-nbase opening angle (cos(phi_(P3n)lb)) are filled and saved to output root file
2D profiles of the entanglement proxies vs top quark scattering angle (cosTS) and ttbar invariant mass (Mtt) are filled and also saved to output root file
'''


# ---------------------------------------------- config ---------------------------------------------- #
base  = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/"
inFiles = [base + "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root" % s
         for s in ("electron", "electron2", "muon", "muon2")]
outFile = "gen_spincorr_TTToSemiLeptonic_all.root"
NMAX  = 50000 # TEST first with a small number; set to -1 for the full run
MULT  = 5.0   # based on "forward-backward" asymmetry integral calculation using l&b spin analyzers
# ---------------------------------------------------------------------------------------------------- #

# Extract the TLorentzVector from a GenParticle
def tlv(genp):
    q = genp.v4(); v = ROOT.TLorentzVector()
    v.SetPtEtaPhiE(q.pt(), q.eta(), q.phi(), q.energy()); return v

# Compute the spin correlation projections at generator level for a semi-leptonic ttbar decay, using a single Bernreuther basis for top quark
# return spin analyzer's opening angle, reflected-nbase opening angle, top quark scattering angle, and ttbar invariant mass
def spincorr(ttg):
    # lab frame
    lep = ttg.ChargedLepton()                             # lepton
    qpos = (lep.pdgId() < 0)                              # lepton charge boolean
    PosTop = tlv(ttg.Top()); NegTop = tlv(ttg.Antitop())  # top and antitop 4vectors
    lepv, bv = tlv(lep), tlv(ttg.BHad())                  # lepton and b-jet 4vectors
    ttbar = PosTop + NegTop; Mtt = ttbar.M()              # ttbar 4vector and invariant mass
    bcm = -ttbar.BoostVector()                            # boost vector = -ttbar
    # ttbar rest frame
    for v in (PosTop, NegTop, lepv, bv): v.Boost(bcm)     # boost all 4vectors to ttbar rest frame

    # build Bernreuther basis
    p_beam = ROOT.TVector3(0., 0., 1.)                    # p_beam direction unit-3vector = p
    k_axis = PosTop.Vect().Unit()                         # direction of top quark in ttbar rest frame unit-3vector = k
    cosTS = k_axis.Dot(p_beam)                            # top quark scattering angle = y
    if (1.0 - cosTS*cosTS) <= 1e-12: return None          # veto if cosTS is too close to 1 or -1 (numerical issues)
    r_axis = (p_beam - k_axis*cosTS).Unit()               # Bernreuther r-axis unit-3vector: (1/|r|)*(p - yk)
    n_axis = (p_beam.Cross(k_axis)).Unit()                # Bernreuther n-axis unit-3vector: (1/|r|)*(p x k)
    sgn = 1.0 if cosTS > 0 else -1.0                      # sign of top scattering angle
    kbase, rbase, nbase = k_axis, r_axis*sgn, n_axis*sgn  # reference axes used in projections (only top quark reference axes are used for both spin analyzers)

    # build SC variables
    lep_rest, b_rest = ROOT.TLorentzVector(lepv), ROOT.TLorentzVector(bv)
    # boost spin analyzers to their parent top quark rest frame (based on lepton charge)
    if qpos:
        lep_rest.Boost(-PosTop.BoostVector()); b_rest.Boost(-NegTop.BoostVector())
    else:
        lep_rest.Boost(-NegTop.BoostVector()); b_rest.Boost(-PosTop.BoostVector())
    # define directions of spin analyzers (unit-3vectors) in their parent top quark rest frame
    lu, bu = lep_rest.Vect().Unit(), b_rest.Vect().Unit()
    # compute projections of spin analyzers onto the top quark reference axes, based on lepton charge
    if qpos:
        c1 = (lu.Dot(kbase), lu.Dot(rbase), lu.Dot(nbase))
        c2 = (bu.Dot(kbase), bu.Dot(rbase), bu.Dot(nbase))
    else:
        c1 = (bu.Dot(kbase), bu.Dot(rbase), bu.Dot(nbase))
        c2 = (lu.Dot(kbase), lu.Dot(rbase), lu.Dot(nbase))

    return lu.Dot(bu), c1[0]*c2[0]+c1[1]*c2[1]-c1[2]*c2[2], cosTS, Mtt

# -------------------------------------------------------------- book histograms --------------------------------------------------------------
# 1D histograms
h_cHel    = ROOT.TH1F("cHel_gen",     "gen cos(#phi_{lb});cHel;events", 24, -1, 1)          # threshold entanglement proxy
h_cHelP3n = ROOT.TH1F("cHel_P3n_gen", "gen cos(#phi_{(P3n)lb});cHel_P3n;events", 24, -1, 1) # boosted entanglement proxy
# cosTS vs Mtt 2D profiles
cosb = array('d', [-1.,-.8,-.6,-.4,-.2,0.,.2,.4,.6,.8,1.])
mb   = array('d', [200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,1000.,1100.,1200.,1300.,1400.,1500.,2000.])
p_cHel    = ROOT.TProfile2D("cHel_coeff_Mtt_vs_cosThetaStar_gen",     "gen D;cos#theta*;M_{tt} [GeV];coeff",      10, cosb, 22, mb)
p_cHelP3n = ROOT.TProfile2D("cHel_P3n_coeff_Mtt_vs_cosThetaStar_gen", "gen Dtilde;cos#theta*;M_{tt} [GeV];coeff", 10, cosb, 22, mb)


# total counters over all inFiles
tot_semiLep = 0  # SemiLeptonic events
tot_weight = 0.0 # sum of event weights
tot_weightXthreshSign = 0.0  # weightedEvents*ThresholdProxySign
tot_weightXboostSign = 0.0   # weightedEvents*BoostedProxySign

# ---------------- loop over inFiles ----------------
for fn in inFiles:
    name = fn.split("MC.")[-1].replace(".root", "")
    rootFile = ROOT.TFile.Open(fn)
    if (not rootFile) or rootFile.IsZombie():
        print("SKIP (cannot open):", name); continue
    AnalysisTree = rootFile.Get("AnalysisTree")
    # disable all branches except the ones we need (to speed up reading)
    AnalysisTree.SetBranchStatus("*", 0)
    AnalysisTree.SetBranchStatus("GenParticles*", 1)
    AnalysisTree.SetBranchStatus("eventweight", 1)
    N = AnalysisTree.GetEntries(); N = N if NMAX < 0 else min(N, NMAX)

    # per file counters:
    f_semiLep = 0  
    f_weight = f_weightXthreshSign = f_weightXboostSign = 0.0   

    # begin loop over events
    for i in range(N):
        AnalysisTree.GetEntry(i)
        eventWeight = AnalysisTree.eventweight
        ttg = TTbarGenPy(AnalysisTree.GenParticles, False)
        if not ttg.IsSemiLeptonicDecay(): continue         # only semi-leptonic ttbar decays
        SCvars = spincorr(ttg)
        if SCvars is None: continue
        CHel, CHelP3n, cosTS, Mtt = SCvars
        CHel_sign  = 1.0 if CHel    >= 0 else -1.0   # threshold proxy sign
        CHelP3n_sign = 1.0 if CHelP3n >= 0 else -1.0 # boosted proxy sign
        # fill histograms with proxy values
        h_cHel.Fill(CHel, eventWeight); h_cHelP3n.Fill(CHelP3n, eventWeight)
        # fill profiles with proxy_sign*MULTIPLIER
        p_cHel.Fill(   cosTS, Mtt, MULT*CHel_sign,    eventWeight)
        p_cHelP3n.Fill(cosTS, Mtt, MULT*CHelP3n_sign, eventWeight)
        # update counters
        f_semiLep += 1; f_weight += eventWeight; f_weightXthreshSign += eventWeight*CHel_sign; f_weightXboostSign += eventWeight*CHelP3n_sign
        if i % 500000 == 0: print("  %-12s %d/%d" % (name, i, N))
    rootFile.Close()
    # compute inclusive entanglement proxies for this file
    threshProxy = MULT*f_weightXthreshSign/f_weight  if f_weight else 0.0
    boostProxy =  MULT*f_weightXboostSign/f_weight if f_weight else 0.0
    print("%-30s SL=%d  D=%.4f  Dtilde=%.4f" % (name, f_semiLep, threshProxy, boostProxy))
    # update total counters over all inFiles
    tot_semiLep += f_semiLep; tot_weight += f_weight; tot_weightXthreshSign += f_weightXthreshSign; tot_weightXboostSign += f_weightXboostSign

print("\nCOMBINED  SL=%d" % tot_semiLep)
print("  gen D      = %.4f   (muon-only reco ref: 0.1271)" % (MULT*tot_weightXthreshSign /tot_weight if tot_weight else 0))
print("  gen Dtilde = %.4f   (muon-only reco ref: 0.3645)" % (MULT*tot_weightXboostSign/tot_weight if tot_weight else 0))

# create output file and write histograms
fout = ROOT.TFile(outFile, "RECREATE")
for hobj in (h_cHel, h_cHelP3n, p_cHel, p_cHelP3n): hobj.Write()
fout.Close()
print("wrote", outFile)