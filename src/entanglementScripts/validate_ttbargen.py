# validation script for TTbarGenPy class, to be run in a CMSSW environment with ROOT and UHH2 installed
# prints out the decay channel breakdown of ttbar events in a given ROOT file, along with some details of semi-leptonic decays

from __future__ import print_function
import ROOT
from ttbargen_py import TTbarGenPy

names = {0:"had",1:"ehad",2:"muhad",3:"tauhad",4:"ee",5:"mumu",6:"tautau",7:"emu",8:"etau",9:"mutau",10:"notfound"}
fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_electron.root"
# show which file we are probing
print("Opening file: %s" % fn.split("/")[-1])
f = ROOT.TFile.Open(fn); t = f.Get("AnalysisTree")

from collections import Counter
c = Counter(); nshow = 0
N = 100000
for i in range(N):
    t.GetEntry(i)
    ttg = TTbarGenPy(t.GenParticles, False)
    c[ttg.DecayChannel()] += 1
    if ttg.IsSemiLeptonicDecay() and nshow < 10:
        nshow += 1
        l, b = ttg.ChargedLepton(), ttg.BHad()
        tl, th = ttg.TopLep(), ttg.TopHad()
        print("  evt %d: lep pdg=%d pt=%.1f | BHad pt=%.1f eta=%.2f | topLep pt=%.1f topHad pt=%.1f"
              % (i, l.pdgId(), l.v4().pt(), b.v4().pt(), b.v4().eta(), tl.v4().pt(), th.v4().pt()))
print("\n--- decay-channel breakdown over %d events ---" % N)
for k in sorted(c):
    print("  %-9s %6d  (%.1f%%)" % (names[k], c[k], 100.0*c[k]/N))


# results for muon file of analysis output:
'''
Opening file: uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root
  evt 0: lep pdg=13 pt=109.1 | BHad pt=42.8 eta=-0.12 | topLep pt=179.3 topHad pt=654.9
  evt 1: lep pdg=13 pt=76.0 | BHad pt=72.1 eta=-0.65 | topLep pt=190.5 topHad pt=156.6
  evt 2: lep pdg=-13 pt=83.5 | BHad pt=47.8 eta=0.29 | topLep pt=126.6 topHad pt=24.3
  evt 3: lep pdg=-13 pt=99.9 | BHad pt=148.9 eta=-0.04 | topLep pt=116.6 topHad pt=156.2
  evt 4: lep pdg=13 pt=109.6 | BHad pt=40.8 eta=1.05 | topLep pt=299.6 topHad pt=40.9
  evt 5: lep pdg=13 pt=183.9 | BHad pt=149.7 eta=0.06 | topLep pt=291.2 topHad pt=280.0
  evt 6: lep pdg=13 pt=70.3 | BHad pt=47.6 eta=1.60 | topLep pt=91.7 topHad pt=141.9
  evt 7: lep pdg=13 pt=66.3 | BHad pt=152.0 eta=1.30 | topLep pt=140.6 topHad pt=190.7
  evt 8: lep pdg=-13 pt=94.6 | BHad pt=33.5 eta=-0.35 | topLep pt=258.3 topHad pt=207.8
  evt 9: lep pdg=13 pt=37.6 | BHad pt=27.5 eta=-0.41 | topLep pt=151.2 topHad pt=90.9

--- decay-channel breakdown over 100000 events ---
  ehad          22  (0.0%)
  muhad      92099  (92.1%)
  tauhad      7850  (7.8%)
  ee            29  (0.0%)
'''


# results for electron file of analysis output:
'''
Opening file: uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_electron.root
  evt 0: lep pdg=15 pt=193.8 | BHad pt=180.2 eta=-0.34 | topLep pt=318.1 topHad pt=325.7
  evt 1: lep pdg=11 pt=172.4 | BHad pt=75.1 eta=-1.93 | topLep pt=289.0 topHad pt=159.8
  evt 2: lep pdg=11 pt=84.7 | BHad pt=54.5 eta=-0.16 | topLep pt=157.6 topHad pt=110.6
  evt 3: lep pdg=11 pt=54.5 | BHad pt=36.0 eta=-1.19 | topLep pt=56.4 topHad pt=81.5
  evt 4: lep pdg=-11 pt=43.5 | BHad pt=81.2 eta=1.23 | topLep pt=143.8 topHad pt=183.8
  evt 5: lep pdg=-11 pt=103.8 | BHad pt=160.0 eta=0.69 | topLep pt=273.9 topHad pt=288.8
  evt 6: lep pdg=-11 pt=48.4 | BHad pt=77.6 eta=-0.69 | topLep pt=51.4 topHad pt=149.3
  evt 7: lep pdg=11 pt=40.1 | BHad pt=103.2 eta=0.50 | topLep pt=177.1 topHad pt=150.0
  evt 8: lep pdg=11 pt=66.0 | BHad pt=98.5 eta=-0.48 | topLep pt=177.3 topHad pt=178.6
  evt 9: lep pdg=11 pt=94.6 | BHad pt=102.5 eta=-1.33 | topLep pt=92.8 topHad pt=178.6

--- decay-channel breakdown over 100000 events ---
  had           18  (0.0%)
  ehad       93548  (93.5%)
  muhad          6  (0.0%)
  tauhad      6415  (6.4%)
  ee            13  (0.0%)
'''