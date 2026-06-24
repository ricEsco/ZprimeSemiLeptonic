from __future__ import print_function
import ROOT
from ttbargen_py import TTbarGenPy

names = {0:"had",1:"ehad",2:"muhad",3:"tauhad",4:"ee",5:"mumu",6:"tautau",7:"emu",8:"etau",9:"mutau",10:"notfound"}
fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root"
f = ROOT.TFile.Open(fn); t = f.Get("AnalysisTree")

from collections import Counter
c = Counter(); nshow = 0
N = 3000
for i in range(N):
    t.GetEntry(i)
    ttg = TTbarGenPy(t.GenParticles, False)
    c[ttg.DecayChannel()] += 1
    if ttg.IsSemiLeptonicDecay() and nshow < 5:
        nshow += 1
        l, b = ttg.ChargedLepton(), ttg.BHad()
        tl, th = ttg.TopLep(), ttg.TopHad()
        print("evt %d: lep pdg=%d pt=%.1f | BHad pt=%.1f eta=%.2f | topLep pt=%.1f topHad pt=%.1f"
              % (i, l.pdgId(), l.v4().pt(), b.v4().pt(), b.v4().eta(), tl.v4().pt(), th.v4().pt()))
print("\n--- decay-channel breakdown over %d events ---" % N)
for k in sorted(c):
    print("  %-9s %6d  (%.1f%%)" % (names[k], c[k], 100.0*c[k]/N))