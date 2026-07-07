# dump_uhh2_gen.py -- UHH2/TTbarGenPy side of the gen closure test.
# Reads a UHH2 pre-selection ntuple (AnalysisTree) and, per event, writes:
#   run, luminosityBlock, event, gen_channel, gen_semilep, gen_mtt, gen_cosThetaStar,
#   gen_cHel, gen_cHel_P3n, genweight
# using the SAME selection (TTbarGenPy) and SAME spin recipe as the skim, so this is the
# actual skim gen chain. Match these to the NanoAOD dump (dump_nano) by (run,lumi,event).
#
# Run in the EL7 container (CMSSW_10_6_28 / py2.7 / ROOT 6.14), same as the skim:
#   python dump_uhh2_gen.py  Ntuple_1.root [Ntuple_2.root ...]
from __future__ import print_function
import sys
import ROOT
from array import array
from ttbargen_py import TTbarGenPy
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

NMAX = -1          # per-file cap for a quick test; -1 = all
base_infile = "/pnfs/desy.de/cms/tier2/store/group/uhh/uhh2ntuples/RunII_106X_v2/UL18/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_TTToSemiLeptonic_CP5_powheg-pythia8_Summer20UL18_v2/211116_134348/0000/"
OUT  = "gendump_uhh2.root"


def tlv(genp):
    q = genp.v4(); v = ROOT.TLorentzVector()
    v.SetPtEtaPhiE(q.pt(), q.eta(), q.phi(), q.energy()); return v


def spincorr(ttg):
    """(mtt, cosThetaStar, cHel, cHel_P3n) -- identical recipe to the skim's gen_observables
    (charge-ordered lb analyzers, Bernreuther basis). Returns None if geometry is degenerate."""
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
    return (Mtt, cosTS, lu.Dot(bu), c1k * c2k + c1r * c2r - c1n * c2n)


# ---- output tree ----
fout = ROOT.TFile(OUT, "RECREATE")
t = ROOT.TTree("gendump", "UHH2/TTbarGenPy gen closure dump")
b = {"run": array('i', [0]), "lumi": array('i', [0]), "event": array('l', [0]),
     "gen_channel": array('i', [0]), "gen_semilep": array('i', [0]),
     "gen_mtt": array('f', [0.]), "gen_cosThetaStar": array('f', [0.]),
     "gen_cHel": array('f', [0.]), "gen_cHel_P3n": array('f', [0.]), "genweight": array('f', [0.])}
spec = {"event": "L"}
for n in ("run", "lumi", "event", "gen_channel", "gen_semilep"):
    t.Branch(n, b[n], n + "/" + spec.get(n, "I"))
for n in ("gen_mtt", "gen_cosThetaStar", "gen_cHel", "gen_cHel_P3n", "genweight"):
    t.Branch(n, b[n], n + "/F")

SENT = -10.0
nfiles = sys.argv[1:] or ["Ntuple_1.root"]
nfiles = [base_infile + f for f in nfiles]
ntot = 0
for fn in nfiles:
    f = ROOT.TFile.Open(fn)
    if (not f) or f.IsZombie():
        print("SKIP (cannot open):", fn); continue
    tin = f.Get("AnalysisTree")
    tin.SetBranchStatus("*", 0)
    for br in ("GenParticles*", "*m_originalXWGTUP*", "run", "luminosityBlock", "event"):
        tin.SetBranchStatus(br, 1)
    gw_leaf = tin.GetLeaf("genInfo.m_originalXWGTUP") or tin.GetLeaf("m_originalXWGTUP")
    N = tin.GetEntries(); N = N if NMAX < 0 else min(N, NMAX)
    print("%s : %d entries" % (fn, N))
    for i in range(N):
        tin.GetEntry(i)
        b["run"][0]   = int(tin.run)
        b["lumi"][0]  = int(tin.luminosityBlock)
        b["event"][0] = int(tin.event)
        b["genweight"][0] = gw_leaf.GetValue() if gw_leaf else 1.0
        ttg = TTbarGenPy(tin.GenParticles, False)
        b["gen_channel"][0] = ttg.DecayChannel()
        res = spincorr(ttg) if ttg.IsSemiLeptonicDecay() else None
        if res is not None:
            b["gen_semilep"][0] = 1
            b["gen_mtt"][0], b["gen_cosThetaStar"][0], b["gen_cHel"][0], b["gen_cHel_P3n"][0] = res
        else:
            b["gen_semilep"][0] = 0
            for n in ("gen_mtt", "gen_cosThetaStar", "gen_cHel", "gen_cHel_P3n"):
                b[n][0] = SENT
        t.Fill(); ntot += 1
    f.Close()

fout.cd(); t.Write(); fout.Close()
print("wrote %s : %d events" % (OUT, ntot))
