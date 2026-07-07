# build_map_uhh2.py -- build the inclusive alpha strength map with TTbarGenPy (the skim's gen code) on the UHH2 pre-selection ntuples.
# Because it uses the SAME selection + spin recipe as the skim, on the SAME-lineage GenParticles, the map-O equals the skim-O by construction --
# this resolves the nanoGen<->TTbarGen closure discrepancy. The ntuples are inclusive at gen level (verified: gen_mtt threshold ratio 0.99, cosTheta* edge 1.00), so alpha is unbiased.
#
# Run in the EL7 container (CMSSW_10_6_28 / py2.7 / ROOT 6.14), like the skim. One condor job per ntuple; 
# hadd the per-job maps afterwards (TProfile2D merges correctly under hadd):
#   python build_map_uhh2.py  0000/Ntuple_123.root
from __future__ import print_function
import sys
import os
import ROOT
from ttbargen_py import TTbarGenPy
import strength_map
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# ONE production only -- same-lineage GenParticles as the skim. Do NOT mix the three timestamp
# dirs (211116_134348 / 230626_211359 / 240306_122806): they re-ntuplize the SAME MiniAOD, so
# combining them would multiply-count events.
BASE   = "/pnfs/desy.de/cms/tier2/store/group/uhh/uhh2ntuples/RunII_106X_v2/UL18/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_TTToSemiLeptonic_CP5_powheg-pythia8_Summer20UL18_v2/211116_134348/"
OUTDIR = "/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/reweighting/strengthmaps_uhh2/"   # per-job maps land here
NMAX   = -1  # per-file cap for a quick test; -1 = all


def tlv(genp):
    q = genp.v4(); v = ROOT.TLorentzVector()
    v.SetPtEtaPhiE(q.pt(), q.eta(), q.phi(), q.energy()); return v


def gen_lb_observables(ttg):
    """(mtt, cosThetaStar, obs) with the 17 lb_ observables strength_map expects -- IDENTICAL
    recipe/pairing (lepton + hadronic-b, Bernreuther basis) to the skim. None if degenerate."""
    lep = ttg.ChargedLepton(); qpos = (lep.pdgId() < 0)
    Top = tlv(ttg.Top()); Anti = tlv(ttg.Antitop())
    lepv, bv = tlv(lep), tlv(ttg.BHad())
    ttbar = Top + Anti; Mtt = ttbar.M()
    bcm = -ttbar.BoostVector()
    for v in (Top, Anti, lepv, bv): v.Boost(bcm)
    beam = ROOT.TVector3(0., 0., 1.)
    k = Top.Vect().Unit(); cosTS = k.Dot(beam)
    if (1.0 - cosTS * cosTS) <= 1e-12:
        return None
    r = (beam - k * cosTS).Unit(); n = beam.Cross(k).Unit()
    s = 1.0 if cosTS > 0 else -1.0
    kb, rb, nb = k, r * s, n * s
    lr = ROOT.TLorentzVector(lepv); br = ROOT.TLorentzVector(bv)
    if qpos:
        lr.Boost(-Top.BoostVector());  br.Boost(-Anti.BoostVector())
    else:
        lr.Boost(-Anti.BoostVector()); br.Boost(-Top.BoostVector())
    lu, bu = lr.Vect().Unit(), br.Vect().Unit()
    a1, a2 = (lu, bu) if qpos else (bu, lu)
    c1k, c1r, c1n = a1.Dot(kb), a1.Dot(rb), a1.Dot(nb)
    c2k, c2r, c2n = a2.Dot(kb), a2.Dot(rb), a2.Dot(nb)
    obs = {
        "lb_cHel":        lu.Dot(bu),
        "lb_cHel_P3n":    c1k * c2k + c1r * c2r - c1n * c2n,
        "lb_cos_theta1k": c1k, "lb_cos_theta1r": c1r, "lb_cos_theta1n": c1n,   # B1_i
        "lb_cos_theta2k": c2k, "lb_cos_theta2r": c2r, "lb_cos_theta2n": c2n,   # B2_i
        "lb_Cnn": c1n * c2n, "lb_Cnr": c1n * c2r, "lb_Cnk": c1n * c2k,
        "lb_Crn": c1r * c2n, "lb_Crr": c1r * c2r, "lb_Crk": c1r * c2k,
        "lb_Ckn": c1k * c2n, "lb_Ckr": c1k * c2r, "lb_Ckk": c1k * c2k,
    }
    return (Mtt, cosTS, obs)


args = sys.argv[1:] or ["0000/Ntuple_1.root"]
ntot = nfill = 0
for rel in args:
    print("processing:", rel)
    fn = BASE + rel
    f = ROOT.TFile.Open(fn)
    if (not f) or f.IsZombie():
        print("SKIP (cannot open):", fn); continue
    t = f.Get("AnalysisTree")
    t.SetBranchStatus("*", 0)
    for br in ("GenParticles*", "*m_originalXWGTUP*"):    # gen record + nominal gen weight only (skip the 16 GB systweights)
        t.SetBranchStatus(br, 1)
    gw_leaf = t.GetLeaf("genInfo.m_originalXWGTUP") or t.GetLeaf("m_originalXWGTUP")
    N = t.GetEntries(); N = N if NMAX < 0 else min(N, NMAX)
    print("%s : %d entries" % (fn, N))
    for i in range(N):
        t.GetEntry(i)
        ntot += 1
        ttg = TTbarGenPy(t.GenParticles, False)
        if not ttg.IsSemiLeptonicDecay():
            continue
        res = gen_lb_observables(ttg)
        if res is None:
            continue
        mtt, cosTS, obs = res
        gw = gw_leaf.GetValue() if gw_leaf else 1.0
        strength_map.fill(mtt, cosTS, gw, obs)
        nfill += 1
    f.Close()

try:
    os.makedirs(OUTDIR)
except OSError:
    pass
tag = args[0].replace("/", "_").replace(".root", "")
outpath = os.path.join(OUTDIR, "NoSC_strengthmap_lb_" + tag + ".root")
strength_map.write(outpath)
print("processed %d events, filled %d semileptonic -> %s" % (ntot, nfill, outpath))
