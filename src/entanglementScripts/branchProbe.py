from __future__ import print_function
import ROOT
fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root"
f = ROOT.TFile.Open(fn); t = f.Get("AnalysisTree")
have = set(b.GetName() for b in t.GetListOfBranches())
for w in ["cHel","cHel_P3n","cHel_P3n_Mtt800_Inf","M_tt","beta","dyreco","eventweight"]:
    print("%-24s %s" % (w, "FOUND" if w in have else "MISSING"))
for i in range(200):                      # show one reconstructed event
    t.GetEntry(i)
    if t.cHel > -1.5:                      # -10 = not reconstructed
        print("evt %d: reco cHel=%.4f cHel_P3n=%.4f M_tt=%.1f beta=%.3f eventweight=%.4g"
              % (i, t.cHel, t.cHel_P3n, t.M_tt, t.beta, t.eventweight)); break