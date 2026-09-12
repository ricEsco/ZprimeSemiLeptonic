from __future__ import print_function
import ROOT

# fn = "/data/dust/group/cms/zprime-uhh/Presel_UL18/workdir_Preselection_UL18/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_1.root"                                                 # my Preselection Output
# fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root"                             # my Analysis Output
fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/reprocessedSignalWtopologycut/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root" # my reprocessed Analysis Output

# fn = "/data/dust/group/cms/zprime-uhh/Presel_UL18_templatemethod/workdir_Preselection_UL18_templatemethod_ttbar_newsystematics/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_3.root"         # Beren's Preselection output
# fn = "/data/dust/group/cms/zprime-uhh/Analysis_UL18_templatemethod/muon/workdir_Analysis_UL18_muon_templatemethod_ttbar_newsystematics/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_3.root" # Beren's Analysis output

f = ROOT.TFile.Open(fn)
t = f.Get("AnalysisTree")
print("entries:", t.GetEntries())

print("\n=== gen-related branches ===")
for b in t.GetListOfBranches():
    name = b.GetName()
    if "gen" in name.lower():
        print("  {0:40s} class='{1}'".format(name, b.GetClassName()))

print("\n=== GenParticles structure ===")
br = t.GetBranch("GenParticles")     # adjust if the gen-particle branch is named differently above
if not br:
    print("  no branch named 'GenParticles' - use a name from the list above")
else:
    print("  class:", br.GetClassName())
    for sb in br.GetListOfBranches():                 # split storage -> sub-branches
        print("  subbranch: {0:45s} {1}".format(sb.GetName(), sb.GetClassName()))
    for lf in br.GetListOfLeaves():                   # the actual stored fields
        print("  leaf:      {0:45s} {1}".format(lf.GetName(), lf.GetTypeName()))

# optional: object access (your CMSSW/UHH2 env should auto-load the GenParticle dictionary)
try:
    t.GetEntry(0)
    gps = t.GenParticles
    print("\n=== event 0: %d gen particles; tops found: ===" % gps.size())
    for p in gps:
        if abs(p.pdgId()) == 6:
            print("  pdgId %d  pt %.1f" % (p.pdgId(), p.v4().pt()))
except Exception as e:
    print("\nobject access not available ->", e)
    print("(no problem - the leaf/subbranch list above is what we need)")


