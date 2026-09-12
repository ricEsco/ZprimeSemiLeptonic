"""
Description:
Build inclusive alpha strength map with TTbarGenPy module using inclusive gen-colleciton of the UHH2 TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8 Ntuples
Use the SAME gen-particle ttbar system reconstruction & SC variable recipe as the Analysis output SC Skim -> map-Obs equals skim-Obs by construction
alpha is constructed with all gen-level quantities and this is itself a purely gen-level map

alpha(m_tt, cosTheta*) = prefactor * <gen_weight * sign_Obs> 

Run in the EL7 container (CMSSW_10_6_28 / py2.7 / ROOT 6.14)
Script produces one condor job per ntuple, hadd the per-job maps afterwards (TProfile2D merges uncertainties correctly under the hadd)

usage: python build_map_uhh2.py  0000/Ntuple_123.root
"""
from __future__ import print_function
import sys
import os
import ROOT
from ttbargen_py import TTbarGenPy
import strength_map
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# PNFS path to the UHH2 ntuples (TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunII_106X_v2_UL18 sample)
BASE   = "/pnfs/desy.de/cms/tier2/store/group/uhh/uhh2ntuples/RunII_106X_v2/UL18/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_TTToSemiLeptonic_CP5_powheg-pythia8_Summer20UL18_v2/211116_134348/"
# output directory
OUTDIR = "/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/reweighting/strengthmaps_uhh2/"
NMAX   = -1  # per-file cap for a quick test and -1 = all

# convert a TTbarGenPy particle to a ROOT TLorentzVector
def TLvec(genp):
    particle = genp.v4()
    vector = ROOT.TLorentzVector()
    vector.SetPtEtaPhiE(particle.pt(), particle.eta(), particle.phi(), particle.energy())
    return vector

# Define SC observables for the (lepton & hadronic-b) spin analyzers
def gen_lb_observables(ttg):
    # gen-particles from ttbar decay system
    lep        = ttg.ChargedLepton()  # lepton from W decay (e, mu, or tau)
    qpos       = (lep.pdgId() < 0)    # lepton charge -> True(False)=Negative(Positive) -> top(antitop) decays leptonically
    TopV       = TLvec(ttg.Top())     # top quark TLorentzVector
    antiTopV   = TLvec(ttg.Antitop()) # antitop quark TLorentzVector
    LepV       = TLvec(lep)           # lepton TLorentzVectors
    hadBquarkV = TLvec(ttg.BHad())    # hadronic b-quark TLorentzVector
    ttbar    = TopV + antiTopV        # ttbar system TLorentzVector
    Mtt      = ttbar.M()              # ttbar invariant mass
    CMboostV = -ttbar.BoostVector()   # -ttbar boost vector

    # Boost all 4-vectors to the ttbar CM frame
    for v in (TopV, antiTopV, LepV, hadBquarkV): v.Boost(CMboostV)
    
    # Bernreuther basis
    beam = ROOT.TVector3(0., 0., 1.) # beam direction
    k = TopV.Vect().Unit()           # top quark direction in CM frame
    cosTS = k.Dot(beam)              # cosine of top quark scattering angle
    if (1.0 - cosTS * cosTS) <= 1e-12: # protect against numerical instability (top quark parallel to beam)
        return None
    r = (beam - k * cosTS).Unit() # transverse direction in the scattering plane
    n = beam.Cross(k).Unit()      # normal direction to the scattering plane
    sign_cTS = 1.0 if cosTS > 0 else -1.0 # Bose symmetry factor
    k_basis, r_basis, n_basis = k, r * sign_cTS, n * sign_cTS # Bose corrected Bernreuther basis

    # prepare rest-frame vectors for the lepton and hadronic b-quark
    LepV_RF = ROOT.TLorentzVector(LepV)
    hadBquarkV_RF = ROOT.TLorentzVector(hadBquarkV)

    # Boost spin analyzers to their parent top/antitop rest frame based on lepton charge
    if qpos: 
        LepV_RF.Boost(-TopV.BoostVector())
        hadBquarkV_RF.Boost(-antiTopV.BoostVector())
    else:    
        LepV_RF.Boost(-antiTopV.BoostVector())
        hadBquarkV_RF.Boost(-TopV.BoostVector())

    # prepare spin-analyzer unit vectors in top/antitop rest frame
    lu, bu = LepV_RF.Vect().Unit(), hadBquarkV_RF.Vect().Unit()
    SA_1, SA_2 = (lu, bu) if qpos else (bu, lu)

    # Compute Spin Analyzer projections onto the Bernreuther basis
    c1k, c1r, c1n = SA_1.Dot(k_basis), SA_1.Dot(r_basis), SA_1.Dot(n_basis)
    c2k, c2r, c2n = SA_2.Dot(k_basis), SA_2.Dot(r_basis), SA_2.Dot(n_basis)
    Obs = {
        "lb_cos_theta1k": c1k, "lb_cos_theta1r": c1r, "lb_cos_theta1n": c1n,   # B1_i
        "lb_cos_theta2k": c2k, "lb_cos_theta2r": c2r, "lb_cos_theta2n": c2n,   # B2_i
        "lb_Cnn": c1n * c2n, "lb_Cnr": c1n * c2r, "lb_Cnk": c1n * c2k,         # C_ij
        "lb_Crn": c1r * c2n, "lb_Crr": c1r * c2r, "lb_Crk": c1r * c2k,
        "lb_Ckn": c1k * c2n, "lb_Ckr": c1k * c2r, "lb_Ckk": c1k * c2k,
        "lb_cHel": lu.Dot(bu),                           # D (threshold entanglement witness)
        "lb_cHel_P3n": c1k * c2k + c1r * c2r - c1n * c2n # D_tilde (boosted entanglement witness)
    }
    return (Mtt, cosTS, Obs)

# Input file -> command line arguments
args = sys.argv[1:] or ["0000/Ntuple_1.root"]
ntot = nfill = 0 # entry/fill counters


# Main loop over input file events
for Ntuple in args:
    print("Processing:", Ntuple)
    inFile = BASE + Ntuple
    file = ROOT.TFile.Open(inFile)
    if (not file) or file.IsZombie():
        print("SKIP (cannot open):", inFile); continue
    TTree = file.Get("AnalysisTree")
    TTree.SetBranchStatus("*", 0)
    # only enable necessary branches: gen record & nominal gen weights
    for Branch in ("GenParticles*", "*m_originalXWGTUP*"):
        TTree.SetBranchStatus(Branch, 1)
    genW_leaf = TTree.GetLeaf("genInfo.m_originalXWGTUP") or TTree.GetLeaf("m_originalXWGTUP")
    N = TTree.GetEntries(); N = N if NMAX < 0 else min(N, NMAX)
    print("%s : %d entries" % (inFile, N))

    # begin loop over events
    for i in range(N):
        TTree.GetEntry(i)
        ntot += 1
        # Use TTbarGenPy to reconstruct the ttbar system from the gen-particles
        ttbar_gen = TTbarGenPy(TTree.GenParticles, False)
        if not ttbar_gen.IsSemiLeptonicDecay():
            print("SKIP Entry", i, "is not semi-leptonic decay; decay channel is", ttbar_gen.DecayChannel())
            continue
        # Compute SC Observables
        SCobs = gen_lb_observables(ttbar_gen)
        if SCobs is None:
            print("SKIP Entry", i, "failed to compute SC observables since cosTS is too close to +/-1 (top quark parallel to beam)")
            continue
        mtt, cosTS, Obs = SCobs
        genW = genW_leaf.GetValue() if genW_leaf else 1.0
        # fill strength maps
        strength_map.fill(mtt, cosTS, genW, Obs)
        nfill += 1
    file.Close()

# Write maps to output ROOT file
try:
    os.makedirs(OUTDIR)
except OSError:
    pass
tag = args[0].replace("/", "_").replace(".root", "")
outpath = os.path.join(OUTDIR, "NoSC_alphaStrengthMap_lb_" + tag + ".root")
strength_map.write(outpath)

# report N-events looped, N-events filled into the maps, and the output path
print("processed %d events, filled %d semileptonic -> %s" % (ntot, nfill, outpath))