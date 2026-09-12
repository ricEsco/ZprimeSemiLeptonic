# dump_ntuples_gen.py -- dump inclusive-gen spin-correlation variables from the UHH2 ntuples with
# TTbarGenPy, into a flat tree that gen_closure_cHel_P3n.py can run on for a TRUE reweighting closure.
#
# Why: the strength maps were built (build_map_uhh2.py) from the INCLUSIVE UHH2 ntuples. The per-event
# skim, by contrast, lives at the reco preDNNselection stage, so its gen collection is a reco-BIASED
# subset -- running the closure on it leaves a small residual asymmetry (A_up ~ 0.005). To close the
# reweighting to A_up = 0 we must run over the SAME inclusive gen sample the map came from. This dumps
# it with the IDENTICAL gen recipe (gen_lb_observables, copied verbatim from build_map_uhh2.py), so
# dump-O == map-O by construction.
#
# Output tree "gendump" carries: gen_semilep, gen_M_tt, gen_cosThetaStar, genweight, and gen_<obs>
# for the full 17 lb observables (gen_cHel, gen_cHel_P3n, gen_Cij, gen_cos_theta{1,2}{k,r,n}) -- the
# branch names gen_closure_cHel_P3n.py expects ("gen_%s" % obs). One condor job per CHUNK; hadd after.
#
# EL7 container (CMSSW_10_6_28 / py2.7 / ROOT 6.14), same env as the map build:
#   python dump_ntuples_gen.py  0000/Ntuple_1.root 0000/Ntuple_10.root ...
from __future__ import print_function
import sys
import os
import array
import ROOT
from ttbargen_py import TTbarGenPy
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# ONE production only -- same-lineage GenParticles as the map/skim. Do NOT mix the three timestamp dirs.
BASE = "/pnfs/desy.de/cms/tier2/store/group/uhh/uhh2ntuples/RunII_106X_v2/UL18/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_TTToSemiLeptonic_CP5_powheg-pythia8_Summer20UL18_v2/211116_134348/"
OUTDIR = "gendump_uhh2ntuples"
NMAX   = -1     # per-file cap for a quick test; -1 = all

# 17 lb observables, in the SAME order/definition as build_map_uhh2.py -> strength_map MAP_MULT
OBSKEYS = ["lb_cHel", "lb_cHel_P3n",
           "lb_cos_theta1k", "lb_cos_theta1r", "lb_cos_theta1n",
           "lb_cos_theta2k", "lb_cos_theta2r", "lb_cos_theta2n",
           "lb_Cnn", "lb_Cnr", "lb_Cnk",
           "lb_Crn", "lb_Crr", "lb_Crk",
           "lb_Ckn", "lb_Ckr", "lb_Ckk"]


def tlv(genp):
    q = genp.v4(); v = ROOT.TLorentzVector()
    v.SetPtEtaPhiE(q.pt(), q.eta(), q.phi(), q.energy()); return v


def gen_lb_observables(ttg):
    """(mtt, cosThetaStar, obs) -- IDENTICAL recipe to build_map_uhh2.py / the skim. None if degenerate."""
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
        "lb_cos_theta1k": c1k, "lb_cos_theta1r": c1r, "lb_cos_theta1n": c1n,
        "lb_cos_theta2k": c2k, "lb_cos_theta2r": c2r, "lb_cos_theta2n": c2n,
        "lb_Cnn": c1n * c2n, "lb_Cnr": c1n * c2r, "lb_Cnk": c1n * c2k,
        "lb_Crn": c1r * c2n, "lb_Crr": c1r * c2r, "lb_Crk": c1r * c2k,
        "lb_Ckn": c1k * c2n, "lb_Ckr": c1k * c2r, "lb_Ckk": c1k * c2k,
    }
    return (Mtt, cosTS, obs)


args = sys.argv[1:] or ["0000/Ntuple_1.root"]

try:
    os.makedirs(OUTDIR)
except OSError:
    pass
tag = args[0].replace("/", "_").replace(".root", "")
outpath = os.path.join(OUTDIR, "gendump_" + tag + ".root")

fout = ROOT.TFile(outpath, "RECREATE")
tr = ROOT.TTree("gendump", "inclusive gen spin-correlation dump (TTbarGenPy)")
a_sl  = array.array('i', [0]);  tr.Branch("gen_semilep",      a_sl,  "gen_semilep/I")
a_mtt = array.array('f', [0.]); tr.Branch("gen_M_tt",         a_mtt, "gen_M_tt/F")
a_cts = array.array('f', [0.]); tr.Branch("gen_cosThetaStar", a_cts, "gen_cosThetaStar/F")
a_gw  = array.array('f', [0.]); tr.Branch("genweight",        a_gw,  "genweight/F")
aobs = {}
for k in OBSKEYS:
    bn = "gen_" + k[3:]                       # strip "lb_" -> gen_cHel, gen_cHel_P3n, gen_Cnn, ...
    aobs[k] = array.array('f', [0.]); tr.Branch(bn, aobs[k], bn + "/F")

ntot = nwrite = 0
for rel in args:
    fn = BASE + rel
    print("processing:", rel)
    f = ROOT.TFile.Open(fn)
    if (not f) or f.IsZombie():
        print("SKIP (cannot open):", fn); continue
    t = f.Get("AnalysisTree")
    t.SetBranchStatus("*", 0)
    for br in ("GenParticles*", "*m_originalXWGTUP*"):
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
        a_sl[0] = 1
        a_mtt[0] = mtt
        a_cts[0] = cosTS
        a_gw[0] = gw_leaf.GetValue() if gw_leaf else 1.0
        for k in OBSKEYS:
            aobs[k][0] = obs[k]
        tr.Fill()
        nwrite += 1
    f.Close()

fout.cd()
tr.Write()
fout.Close()
print("processed %d events, wrote %d semileptonic -> %s" % (ntot, nwrite, outpath))