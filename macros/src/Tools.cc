#include <TString.h>
#include <TFile.h>

#include "../include/Tools.h"

AnalysisTool::AnalysisTool(bool do_puppi_) : do_puppi(do_puppi_)
{
  tag = "withoutDNN";
  path = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/muon/mergedFiles/uhh2.AnalysisModuleRunner.MC.";

  CMcuts = {"1_lepJet", "3_hadJets", "DeltaR_lepB"};
  
  for(unsigned int i=0; i<CMcuts.size(); i++){
    TString t = "";
    t += CMcuts[i];
    CMcuts_str.emplace_back(t);
  }
}
