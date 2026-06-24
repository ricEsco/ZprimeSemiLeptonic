from __future__ import print_function
import ROOT

# TTbarGen lives in UHH2/common; load it if it didn't auto-load
if not hasattr(ROOT, "TTbarGen"):
    ROOT.gSystem.Load("libUHH2common.so")     # adjust name if your build differs
print("TTbarGen available:", hasattr(ROOT, "TTbarGen"))

fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root"
f = ROOT.TFile.Open(fn)
t = f.Get("AnalysisTree")

n_sl = 0
for i in range(50):
    t.GetEntry(i)
    ttg = ROOT.TTbarGen(t.GenParticles, False)   # False = throw_on_failure off
    if ttg.IsSemiLeptonicDecay():
        n_sl += 1
        b = ttg.BHad().v4()
        l = ttg.ChargedLepton()
        print("evt %d SL: BHad pt=%.1f eta=%.2f phi=%.2f | lep pdg=%d pt=%.1f" % (
              i, b.pt(), b.eta(), b.phi(), l.pdgId(), l.v4().pt()))
print("semileptonic events in first 50:", n_sl)