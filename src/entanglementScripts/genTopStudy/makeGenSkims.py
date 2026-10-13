"""
makeGenSkims.py  --  minimal gen-only skim of the TTToSemiLeptonic (UL18) analysis output

Writes ONLY the gen-level ttbar system built by TTbarGenPy (lab-frame pt/eta/phi/E), in the tree / branch names used by
plotBoostResidual.py (tree "SC_SkimTree"):

  gen_<P>_{pt,eta,phi,E}    P = Top, Antitop, TopLep, TopHad      first copy  (TTbarGen Top()/Antitop(), status 22)
                                WLep, WHad, bLep, bHad,             TTbarGen W (status 22) and b (status 23)
                                Lepton, Neutrino, Q1, Q2            W decay products (Q1, Q2 = hadronic W daughters)
  genLC_<P>_{pt,eta,phi,E}  P = Top, Antitop, TopLep, TopHad      last copy   (W + b = statusFlags().isLastCopy() top)
  gen_semilep (I), gen_channel (I), genweight (F)

Non-semileptonic events are kept with gen_semilep = 0 and all four-vector components set to -10.

Usage (CMSSW_10_6_28 / python 2.7, same environment as makeAnalysisSkims.py; ttbargen_py.py in the same directory):
  python makeGenSkims.py              # all four files -> TTToSemiLeptonicUL18_GenSkim.root
  python makeGenSkims.py 2            # condor $(Process) = 0..3 -> TTToSemiLeptonicUL18_GenSkim_<tag>.root
  python makeGenSkims.py 2 1000000    # ... and stop after 1M events per file (quick turnaround)
"""
from __future__ import print_function
import sys
from array import array
import ROOT
from ttbargen_py import TTbarGenPy

ROOT.gROOT.SetBatch(True)
try:
    xrange  # type: ignore
except NameError:
    xrange = range

# ------------------------------------------------------------------ configuration (same inputs as makeAnalysisSkims.py)
BASE = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/noDNN/UL18/bothLeptons_mergedFiles/"
TAGS = ["electron1", "electron2", "muon1", "muon2"]
ALL_FILES = [BASE + "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root" % tag for tag in TAGS]

if len(sys.argv) > 1:
    idx = int(sys.argv[1])
    FILES = [ALL_FILES[idx]]
    OUTF = "TTToSemiLeptonicUL18_GenSkim_%s.root" % TAGS[idx]
else:
    FILES = ALL_FILES
    OUTF = "TTToSemiLeptonicUL18_GenSkim.root"
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else -1     # events per file; -1 = all
SENT = -10.0

GEN_FC = ["Top", "Antitop", "TopLep", "TopHad", "WLep", "WHad", "bLep", "bHad", "Lepton", "Neutrino", "Q1", "Q2"]
GEN_LC = ["Top", "Antitop", "TopLep", "TopHad"]
COMP = ("pt", "eta", "phi", "E")


# ------------------------------------------------------------------ helpers
def TLvec(genp):
    """TTbarGenPy gen particle -> ROOT.TLorentzVector (identical to makeAnalysisSkims.py)."""
    p = genp.v4()
    v = ROOT.TLorentzVector()
    v.SetPtEtaPhiE(p.pt(), p.eta(), p.phi(), p.energy())
    return v


def ttbar_vectors(ttg):
    """All lab-frame four-vectors of the semileptonic ttbar system: returns (first-copy dict, last-copy dict)."""
    qpos = (ttg.ChargedLepton().pdgId() < 0)          # l+ -> the top decays leptonically (same convention as gen_observables)
    WTop, WAnti = TLvec(ttg.m_WTop), TLvec(ttg.m_WAntitop)
    bTop, bAnti = TLvec(ttg.bTop()), TLvec(ttg.bAntitop())
    t_fc, tb_fc = TLvec(ttg.Top()), TLvec(ttg.Antitop())
    t_lc, tb_lc = WTop + bTop, WAnti + bAnti           # last-copy tops (exact, validated on MiniAOD)
    q1, q2 = (ttg.m_WMinusdecay1, ttg.m_WMinusdecay2) if qpos else (ttg.m_Wdecay1, ttg.m_Wdecay2)
    fc = {
        "Top": t_fc, "Antitop": tb_fc,
        "TopLep": t_fc if qpos else tb_fc, "TopHad": tb_fc if qpos else t_fc,
        "WLep": WTop if qpos else WAnti, "WHad": WAnti if qpos else WTop,
        "bLep": bTop if qpos else bAnti, "bHad": bAnti if qpos else bTop,
        "Lepton": TLvec(ttg.ChargedLepton()), "Neutrino": TLvec(ttg.Neutrino()),
        "Q1": TLvec(q1), "Q2": TLvec(q2),
    }
    lc = {"Top": t_lc, "Antitop": tb_lc, "TopLep": t_lc if qpos else tb_lc, "TopHad": tb_lc if qpos else t_lc}
    return fc, lc


def book(tree, prefix, names, store):
    for n in names:
        for c in COMP:
            key = "%s_%s_%s" % (prefix, n, c)
            store[key] = array("f", [SENT])
            tree.Branch(key, store[key], key + "/F")


def fill_vectors(store, prefix, vecs):
    for n, v in vecs.items():
        store["%s_%s_pt" % (prefix, n)][0] = v.Pt()
        store["%s_%s_eta" % (prefix, n)][0] = v.Eta()
        store["%s_%s_phi" % (prefix, n)][0] = v.Phi()
        store["%s_%s_E" % (prefix, n)][0] = v.E()


# ------------------------------------------------------------------ output tree
def main():
    outFile = ROOT.TFile(OUTF, "RECREATE")
    tree = ROOT.TTree("SC_SkimTree", "gen-level ttbar system (first and last copy tops), UL18")
    buf = {}
    book(tree, "gen", GEN_FC, buf)
    book(tree, "genLC", GEN_LC, buf)
    semilep, channel, gw = array("i", [0]), array("i", [0]), array("f", [0.0])
    tree.Branch("gen_semilep", semilep, "gen_semilep/I")
    tree.Branch("gen_channel", channel, "gen_channel/I")
    tree.Branch("genweight", gw, "genweight/F")

    for inFile in FILES:
        f = ROOT.TFile.Open(inFile)
        if (not f) or f.IsZombie():
            print("SKIP (cannot open):", inFile)
            continue
        t = f.Get("AnalysisTree")
        t.SetBranchStatus("*", 0)
        t.SetBranchStatus("GenParticles*", 1)
        t.SetBranchStatus("*m_originalXWGTUP*", 1)
        wleaf = t.GetLeaf("m_originalXWGTUP")
        N = t.GetEntries() if NMAX < 0 else min(t.GetEntries(), NMAX)
        print("Skimming %d entries of %s" % (N, inFile))
        nbad = 0
        for i in xrange(N):
            t.GetEntry(i)
            gw[0] = wleaf.GetValue() if wleaf else 1.0
            ttg = TTbarGenPy(t.GenParticles, False)
            channel[0] = ttg.DecayChannel()
            if ttg.IsSemiLeptonicDecay():
                semilep[0] = 1
                fc, lc = ttbar_vectors(ttg)
                fill_vectors(buf, "gen", fc)
                fill_vectors(buf, "genLC", lc)
            else:
                semilep[0] = 0
                nbad += 1
                for k in buf:
                    buf[k][0] = SENT
            tree.Fill()
            if i and i % 200000 == 0:
                print("  %d / %d" % (i, N))
        print("  done: %d events, %d not semileptonic (kept with gen_semilep = 0)" % (N, nbad))
        f.Close()

    outFile.cd()
    tree.Write()
    outFile.Close()
    print("Wrote", OUTF)


if __name__ == "__main__":
    main()