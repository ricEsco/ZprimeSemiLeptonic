from __future__ import print_function
import ROOT

# fn = "/data/dust/group/cms/zprime-uhh/Presel_UL18/workdir_Preselection_UL18/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_1.root"                                              # my Preselection Output
# fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root"                          # my Analysis Output
fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/muon/mergedFiles/reprocessedSignalWtopologycut/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_2.root" # my reprocessed Analysis Output

# fn = "/data/dust/group/cms/zprime-uhh/Presel_UL18_templatemethod/workdir_Preselection_UL18_templatemethod_ttbar_newsystematics/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_3.root"         # Beren's Preselection output
# fn = "/data/dust/group/cms/zprime-uhh/Analysis_UL18_templatemethod/muon/workdir_Analysis_UL18_muon_templatemethod_ttbar_newsystematics/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_3.root" # Beren's Analysis output

f = ROOT.TFile.Open(fn); t = f.Get("AnalysisTree")
have = set(b.GetName() for b in t.GetListOfBranches())
for w in ["cHel","cHel_P3n","cHel_P3n_Mtt800_Inf","M_tt","beta","dyreco","eventweight"]:
    print("%-24s %s" % (w, "FOUND" if w in have else "MISSING"))
for i in range(200):                      # show one reconstructed event
    t.GetEntry(i)
    if t.cHel > -1.5:                      # -10 = not reconstructed
        print("evt %d: reco cHel=%.4f cHel_P3n=%.4f M_tt=%.1f beta=%.3f eventweight=%.4g"
              % (i, t.cHel, t.cHel_P3n, t.M_tt, t.beta, t.eventweight)); break