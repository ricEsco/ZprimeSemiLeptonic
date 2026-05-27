#include <iostream>
#include <memory>
#include <fstream>
#include <cmath>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include "UHH2/common/include/PrintingModules.h"

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include "UHH2/common/include/LuminosityHists.h"
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/TopPtReweight.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/LeptonScaleFactors.h>
#include <UHH2/common/include/PSWeights.h>
#include "TH1.h"
#include "TFile.h"

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/TTbarLJHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSystematicsHists.h>  // not in analysis module
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicPDFHists.h>          // not in analysis module
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicMulticlassNNHists.h> // not in analysis module
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicCHSMatchHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicMistagHists.h>       // not in analysis module
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>
#include <UHH2/ZprimeSemiLeptonic/include/ElecTriggerSF.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopTagScaleFactor.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopMistagScaleFactor.h>

#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <UHH2/HOTVR/include/HadronicTop.h>
#include <UHH2/HOTVR/include/HOTVRScaleFactor.h>
#include <UHH2/HOTVR/include/HOTVRIds.h>

#include "UHH2/common/include/NeuralNetworkBase.hpp"

using namespace std;
using namespace uhh2;

/*
██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/


class NeuralNetworkModule: public NeuralNetworkBase {
public:
  explicit NeuralNetworkModule(uhh2::Context&, const std::string & ModelName, const std::string& ConfigName);
  virtual void CreateInputs(uhh2::Event & event) override;

protected:
  uhh2::Event::Handle<float> h_Ak4_j1_E;
  uhh2::Event::Handle<float> h_Ak4_j1_eta;
  uhh2::Event::Handle<float> h_Ak4_j1_m;
  uhh2::Event::Handle<float> h_Ak4_j1_phi;
  uhh2::Event::Handle<float> h_Ak4_j1_pt;
  uhh2::Event::Handle<float> h_Ak4_j1_deepjetbscore;

  uhh2::Event::Handle<float> h_Ak4_j2_E;
  uhh2::Event::Handle<float> h_Ak4_j2_eta;
  uhh2::Event::Handle<float> h_Ak4_j2_m;
  uhh2::Event::Handle<float> h_Ak4_j2_phi;
  uhh2::Event::Handle<float> h_Ak4_j2_pt;
  uhh2::Event::Handle<float> h_Ak4_j2_deepjetbscore;

  uhh2::Event::Handle<float> h_Ak4_j3_E;
  uhh2::Event::Handle<float> h_Ak4_j3_eta;
  uhh2::Event::Handle<float> h_Ak4_j3_m;
  uhh2::Event::Handle<float> h_Ak4_j3_phi;
  uhh2::Event::Handle<float> h_Ak4_j3_pt;
  uhh2::Event::Handle<float> h_Ak4_j3_deepjetbscore;

  uhh2::Event::Handle<float> h_Ak4_j4_E;
  uhh2::Event::Handle<float> h_Ak4_j4_eta;
  uhh2::Event::Handle<float> h_Ak4_j4_m;
  uhh2::Event::Handle<float> h_Ak4_j4_phi;
  uhh2::Event::Handle<float> h_Ak4_j4_pt;
  uhh2::Event::Handle<float> h_Ak4_j4_deepjetbscore;

  uhh2::Event::Handle<float> h_Ak4_j5_E;
  uhh2::Event::Handle<float> h_Ak4_j5_eta;
  uhh2::Event::Handle<float> h_Ak4_j5_m;
  uhh2::Event::Handle<float> h_Ak4_j5_phi;
  uhh2::Event::Handle<float> h_Ak4_j5_pt;
  uhh2::Event::Handle<float> h_Ak4_j5_deepjetbscore;

  uhh2::Event::Handle<float> h_Ele_E;
  uhh2::Event::Handle<float> h_Ele_eta;
  uhh2::Event::Handle<float> h_Ele_phi;
  uhh2::Event::Handle<float> h_Ele_pt;

  uhh2::Event::Handle<float> h_MET_phi;
  uhh2::Event::Handle<float> h_MET_pt;

  uhh2::Event::Handle<float> h_Mu_E;
  uhh2::Event::Handle<float> h_Mu_eta;
  uhh2::Event::Handle<float> h_Mu_phi;
  uhh2::Event::Handle<float> h_Mu_pt;

  uhh2::Event::Handle<float> h_N_Ak4;

  uhh2::Event::Handle<float> h_Ak8_j1_E;
  uhh2::Event::Handle<float> h_Ak8_j1_eta;
  uhh2::Event::Handle<float> h_Ak8_j1_mSD;
  uhh2::Event::Handle<float> h_Ak8_j1_phi;
  uhh2::Event::Handle<float> h_Ak8_j1_pt;
  uhh2::Event::Handle<float> h_Ak8_j1_tau21;
  uhh2::Event::Handle<float> h_Ak8_j1_tau32;

  uhh2::Event::Handle<float> h_Ak8_j2_E;
  uhh2::Event::Handle<float> h_Ak8_j2_eta;
  uhh2::Event::Handle<float> h_Ak8_j2_mSD;
  uhh2::Event::Handle<float> h_Ak8_j2_phi;
  uhh2::Event::Handle<float> h_Ak8_j2_pt;
  uhh2::Event::Handle<float> h_Ak8_j2_tau21;
  uhh2::Event::Handle<float> h_Ak8_j2_tau32;

  uhh2::Event::Handle<float> h_Ak8_j3_E;
  uhh2::Event::Handle<float> h_Ak8_j3_eta;
  uhh2::Event::Handle<float> h_Ak8_j3_mSD;
  uhh2::Event::Handle<float> h_Ak8_j3_phi;
  uhh2::Event::Handle<float> h_Ak8_j3_pt;
  uhh2::Event::Handle<float> h_Ak8_j3_tau21;
  uhh2::Event::Handle<float> h_Ak8_j3_tau32;

  uhh2::Event::Handle<float> h_N_Ak8;
};


NeuralNetworkModule::NeuralNetworkModule(Context& ctx, const std::string & ModelName, const std::string& ConfigName): NeuralNetworkBase(ctx, ModelName, ConfigName){
  h_Ak4_j1_E   = ctx.get_handle<float>("Ak4_j1_E");
  h_Ak4_j1_eta = ctx.get_handle<float>("Ak4_j1_eta");
  h_Ak4_j1_m   = ctx.get_handle<float>("Ak4_j1_m");
  h_Ak4_j1_phi = ctx.get_handle<float>("Ak4_j1_phi");
  h_Ak4_j1_pt  = ctx.get_handle<float>("Ak4_j1_pt");
  h_Ak4_j1_deepjetbscore  = ctx.get_handle<float>("Ak4_j1_deepjetbscore");

  h_Ak4_j2_E   = ctx.get_handle<float>("Ak4_j2_E");
  h_Ak4_j2_eta = ctx.get_handle<float>("Ak4_j2_eta");
  h_Ak4_j2_m   = ctx.get_handle<float>("Ak4_j2_m");
  h_Ak4_j2_phi = ctx.get_handle<float>("Ak4_j2_phi");
  h_Ak4_j2_pt  = ctx.get_handle<float>("Ak4_j2_pt");
  h_Ak4_j2_deepjetbscore  = ctx.get_handle<float>("Ak4_j2_deepjetbscore");

  h_Ak4_j3_E   = ctx.get_handle<float>("Ak4_j3_E");
  h_Ak4_j3_eta = ctx.get_handle<float>("Ak4_j3_eta");
  h_Ak4_j3_m   = ctx.get_handle<float>("Ak4_j3_m");
  h_Ak4_j3_phi = ctx.get_handle<float>("Ak4_j3_phi");
  h_Ak4_j3_pt  = ctx.get_handle<float>("Ak4_j3_pt");
  h_Ak4_j3_deepjetbscore  = ctx.get_handle<float>("Ak4_j3_deepjetbscore");

  h_Ak4_j4_E   = ctx.get_handle<float>("Ak4_j4_E");
  h_Ak4_j4_eta = ctx.get_handle<float>("Ak4_j4_eta");
  h_Ak4_j4_m   = ctx.get_handle<float>("Ak4_j4_m");
  h_Ak4_j4_phi = ctx.get_handle<float>("Ak4_j4_phi");
  h_Ak4_j4_pt  = ctx.get_handle<float>("Ak4_j4_pt");
  h_Ak4_j4_deepjetbscore  = ctx.get_handle<float>("Ak4_j4_deepjetbscore");

  h_Ak4_j5_E   = ctx.get_handle<float>("Ak4_j5_E");
  h_Ak4_j5_eta = ctx.get_handle<float>("Ak4_j5_eta");
  h_Ak4_j5_m   = ctx.get_handle<float>("Ak4_j5_m");
  h_Ak4_j5_phi = ctx.get_handle<float>("Ak4_j5_phi");
  h_Ak4_j5_pt  = ctx.get_handle<float>("Ak4_j5_pt");
  h_Ak4_j5_deepjetbscore  = ctx.get_handle<float>("Ak4_j5_deepjetbscore");

  h_Ele_E    = ctx.get_handle<float>("Ele_E");
  h_Ele_eta  = ctx.get_handle<float>("Ele_eta");
  h_Ele_phi  = ctx.get_handle<float>("Ele_phi");
  h_Ele_pt   = ctx.get_handle<float>("Ele_pt");

  h_MET_phi = ctx.get_handle<float>("MET_phi");
  h_MET_pt = ctx.get_handle<float>("MET_pt");

  h_Mu_E    = ctx.get_handle<float>("Mu_E");
  h_Mu_eta  = ctx.get_handle<float>("Mu_eta");
  h_Mu_phi  = ctx.get_handle<float>("Mu_phi");
  h_Mu_pt   = ctx.get_handle<float>("Mu_pt");

  h_N_Ak4 = ctx.get_handle<float>("N_Ak4");

  h_Ak8_j1_E     = ctx.get_handle<float>("Ak8_j1_E");
  h_Ak8_j1_eta   = ctx.get_handle<float>("Ak8_j1_eta");
  h_Ak8_j1_mSD   = ctx.get_handle<float>("Ak8_j1_mSD");
  h_Ak8_j1_phi   = ctx.get_handle<float>("Ak8_j1_phi");
  h_Ak8_j1_pt    = ctx.get_handle<float>("Ak8_j1_pt");
  h_Ak8_j1_tau21 = ctx.get_handle<float>("Ak8_j1_tau21");
  h_Ak8_j1_tau32 = ctx.get_handle<float>("Ak8_j1_tau32");

  h_Ak8_j2_E     = ctx.get_handle<float>("Ak8_j2_E");
  h_Ak8_j2_eta   = ctx.get_handle<float>("Ak8_j2_eta");
  h_Ak8_j2_mSD   = ctx.get_handle<float>("Ak8_j2_mSD");
  h_Ak8_j2_phi   = ctx.get_handle<float>("Ak8_j2_phi");
  h_Ak8_j2_pt    = ctx.get_handle<float>("Ak8_j2_pt");
  h_Ak8_j2_tau21 = ctx.get_handle<float>("Ak8_j2_tau21");
  h_Ak8_j2_tau32 = ctx.get_handle<float>("Ak8_j2_tau32");

  h_Ak8_j3_E     = ctx.get_handle<float>("Ak8_j3_E");
  h_Ak8_j3_eta   = ctx.get_handle<float>("Ak8_j3_eta");
  h_Ak8_j3_mSD   = ctx.get_handle<float>("Ak8_j3_mSD");
  h_Ak8_j3_phi   = ctx.get_handle<float>("Ak8_j3_phi");
  h_Ak8_j3_pt    = ctx.get_handle<float>("Ak8_j3_pt");
  h_Ak8_j3_tau21 = ctx.get_handle<float>("Ak8_j3_tau21");
  h_Ak8_j3_tau32 = ctx.get_handle<float>("Ak8_j3_tau32");

  h_N_Ak8 = ctx.get_handle<float>("N_Ak8");
}

void NeuralNetworkModule::CreateInputs(Event & event){
  NNInputs.clear();
  NNoutputs.clear();

  string varname[59];
  string scal[59];
  string mean[59];
  string std[59];
  double mean_val[59];
  double std_val[59];

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// LEPTON-SPECIFIC NN SETTINGS /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// DON'T FORGET TO CHANGE! ///////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////// ONLY UNCOMMENT LINE CORRESPONDING TO LEPTON FLAVOR YOU WANT TO RUN ////////////////////////////////////////////////
  //Muon
  ifstream normfile ("/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_muon/NormInfo.txt", ios::in);
  //Electron
  // ifstream normfile ("/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_ele/NormInfo.txt", ios::in);
  
  if(!normfile.good()) throw runtime_error("NeuralNetworkModule: The specified norm file does not exist.");
  if (normfile.is_open()){
    for(int i = 0; i < 59; ++i)
    {
      // cout<<varname<<endl;
      normfile >> varname[i] >> scal[i] >> mean[i] >> std[i];
      mean_val[i] = std::stod(mean[i]);
      std_val[i] = std::stod(std[i]);
    }
    normfile.close();
  }
  NNInputs.push_back( tensorflow::Tensor(tensorflow::DT_FLOAT, {1, 59}));

  //////////////////////////////////////////////// ONLY UNCOMMENT LINE CORRESPONDING TO LEPTON FLAVOR YOU WANT TO RUN ////////////////////////////////////////////////
  //Muon
  vector<uhh2::Event::Handle<float>> inputs = {
    h_Ak4_j1_E, h_Ak4_j1_deepjetbscore, h_Ak4_j1_eta, h_Ak4_j1_m, h_Ak4_j1_phi, h_Ak4_j1_pt, 
    h_Ak4_j2_E, h_Ak4_j2_deepjetbscore, h_Ak4_j2_eta, h_Ak4_j2_m, h_Ak4_j2_phi, h_Ak4_j2_pt, 
    h_Ak4_j3_E, h_Ak4_j3_deepjetbscore, h_Ak4_j3_eta, h_Ak4_j3_m, h_Ak4_j3_phi, h_Ak4_j3_pt, 
    h_Ak4_j4_E, h_Ak4_j4_deepjetbscore, h_Ak4_j4_eta, h_Ak4_j4_m, h_Ak4_j4_phi, h_Ak4_j4_pt, 
    h_Ak4_j5_E, h_Ak4_j5_deepjetbscore, h_Ak4_j5_eta, h_Ak4_j5_m, h_Ak4_j5_phi, h_Ak4_j5_pt, 
    h_Ak8_j1_E, h_Ak8_j1_eta, h_Ak8_j1_mSD, h_Ak8_j1_phi, h_Ak8_j1_pt, h_Ak8_j1_tau21, h_Ak8_j1_tau32, 
    h_Ak8_j2_E, h_Ak8_j2_eta, h_Ak8_j2_mSD, h_Ak8_j2_phi, h_Ak8_j2_pt, h_Ak8_j2_tau21, h_Ak8_j2_tau32, 
    h_Ak8_j3_E, h_Ak8_j3_eta, h_Ak8_j3_mSD, h_Ak8_j3_phi, h_Ak8_j3_pt, h_Ak8_j3_tau21, h_Ak8_j3_tau32, 
    h_MET_phi, h_MET_pt, h_Mu_E, h_Mu_eta, h_Mu_phi, h_Mu_pt, h_N_Ak4, h_N_Ak8}; // in alphabetical order to match NormInfo.txt

  // //Electron
  // vector<uhh2::Event::Handle<float>> inputs = {
  //   h_Ak4_j1_E, h_Ak4_j1_deepjetbscore, h_Ak4_j1_eta, h_Ak4_j1_m, h_Ak4_j1_phi, h_Ak4_j1_pt, 
  //   h_Ak4_j2_E, h_Ak4_j2_deepjetbscore, h_Ak4_j2_eta, h_Ak4_j2_m, h_Ak4_j2_phi, h_Ak4_j2_pt, 
  //   h_Ak4_j3_E, h_Ak4_j3_deepjetbscore, h_Ak4_j3_eta, h_Ak4_j3_m, h_Ak4_j3_phi, h_Ak4_j3_pt, 
  //   h_Ak4_j4_E, h_Ak4_j4_deepjetbscore, h_Ak4_j4_eta, h_Ak4_j4_m, h_Ak4_j4_phi, h_Ak4_j4_pt, 
  //   h_Ak4_j5_E, h_Ak4_j5_deepjetbscore, h_Ak4_j5_eta, h_Ak4_j5_m, h_Ak4_j5_phi, h_Ak4_j5_pt, 
  //   h_Ak8_j1_E, h_Ak8_j1_eta, h_Ak8_j1_mSD, h_Ak8_j1_phi, h_Ak8_j1_pt, h_Ak8_j1_tau21, h_Ak8_j1_tau32, 
  //   h_Ak8_j2_E, h_Ak8_j2_eta, h_Ak8_j2_mSD, h_Ak8_j2_phi, h_Ak8_j2_pt, h_Ak8_j2_tau21, h_Ak8_j2_tau32, 
  //   h_Ak8_j3_E, h_Ak8_j3_eta, h_Ak8_j3_mSD, h_Ak8_j3_phi, h_Ak8_j3_pt, h_Ak8_j3_tau21, h_Ak8_j3_tau32, 
  //   h_Ele_E, h_Ele_eta, h_Ele_phi, h_Ele_pt, h_MET_phi, h_MET_pt, h_N_Ak4, h_N_Ak8}; // in alphabetical order to match NormInfo.txt

  ///////////////////////////////////////////////////////////// LEPTON-SPECIFIC NN SETTINGS /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// DON'T FORGET TO CHANGE! ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(int i = 0; i < 59; ++i){
    // cout<<"looping over NN inputs "<< i <<endl;
    NNInputs.at(0).tensor<float, 2>()(0,i)  = (event.get(inputs.at(i))   - mean_val[i]) / (std_val[i]);
  }
  // cout <<"NNinputs size : "<< NNInputs.size()<< " Layer: "<<LayerInputs.size()<<endl;
  if (NNInputs.size()!=LayerInputs.size()) throw logic_error("NeuralNetworkModule.cxx: Create a number of inputs diffetent wrt. LayerInputs.size()="+to_string(LayerInputs.size()));
}



class ZprimeAnalysisModule_applyNN : public ModuleBASE {

public:
  explicit ZprimeAnalysisModule_applyNN(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:

  bool debug;
  
  // Cleaners
  std::unique_ptr<MuonCleaner> muon_cleaner_low, muon_cleaner_high;
  std::unique_ptr<ElectronCleaner> electron_cleaner_low, electron_cleaner_high;

  // scale factors
  unique_ptr<AnalysisModule> sf_muon_iso_stat_low, sf_muon_id_stat_low, sf_muon_id_stat_high, sf_muon_trigger_stat_low, sf_muon_trigger_stat_high;
  unique_ptr<AnalysisModule> sf_muon_iso_syst_low, sf_muon_id_syst_low, sf_muon_id_syst_high, sf_muon_trigger_syst_low, sf_muon_trigger_syst_high;

  unique_ptr<AnalysisModule> sf_muon_iso_stat_low_dummy, sf_muon_id_stat_dummy, sf_muon_trigger_stat_dummy;
  unique_ptr<AnalysisModule> sf_muon_iso_syst_low_dummy, sf_muon_id_syst_dummy, sf_muon_trigger_syst_dummy;

  unique_ptr<AnalysisModule> sf_ele_id_low, sf_ele_id_high, sf_ele_reco;
  unique_ptr<AnalysisModule> sf_ele_id_dummy, sf_ele_reco_dummy;
  unique_ptr<MuonRecoSF> sf_muon_reco;
  unique_ptr<AnalysisModule> sf_ele_trigger;
  unique_ptr<AnalysisModule> sf_btagging;


  // AnalysisModules
  unique_ptr<AnalysisModule> LumiWeight_module, PUWeight_module, TopPtReweight_module, MCScale_module;
  unique_ptr<AnalysisModule> NLOCorrections_module;
  unique_ptr<PSWeights> ps_weights;
  
  // Structure Constants Calculator for EFT
  unique_ptr<StructureConstantsCalculator> structure_constants_calculator;
  uhh2::Event::Handle<std::vector<float>> h_structure_constants;

  // Top tagging
  unique_ptr<HOTVRTopTagger> TopTaggerHOTVR;
  unique_ptr<AnalysisModule> hadronic_top;
  unique_ptr<AnalysisModule> sf_toptag;
  unique_ptr<AnalysisModule> sf_topmistag;
  unique_ptr<DeepAK8TopTagger> TopTaggerDeepAK8;

  // TopTags veto
  unique_ptr<Selection> TopTagVetoSelection;

  // Mass reconstruction
  unique_ptr<ZprimeCandidateBuilder> CandidateBuilder;

  // Chi2 discriminator
  unique_ptr<ZprimeChi2Discriminator> Chi2DiscriminatorZprime;
  unique_ptr<ZprimeCorrectMatchDiscriminator> CorrectMatchDiscriminatorZprime;
  std::unique_ptr<Hists> h_CHSMatchHists;

  // Selections
  unique_ptr<Selection> Chi2_selection, TTbarMatchable_selection, Chi2CandidateMatched_selection, ZprimeTopTag_selection;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  unique_ptr<Selection> TwoDCut_selection_low1;
  std::unique_ptr<Selection> sel_1btag, sel_2btag;
  std::unique_ptr<Selection> HEM_selection, DeltaEta_selection;
  unique_ptr<Selection> ThetaStar_selection_bin1, ThetaStar_selection_bin2, ThetaStar_selection_bin3, ThetaStar_selection_bin4, ThetaStar_selection_bin5, ThetaStar_selection_bin6;

  // NN variables handles
  unique_ptr<Variables_NN> Variables_module;
  // unique_ptr<Variables_EFT_SR> VariablesEFTSR_module;
  // unique_ptr<Variables_EFT_CR1> VariablesEFTCR1_module;
  // unique_ptr<Variables_EFT_CR2> VariablesEFTCR2_module;

  // Gen-level variables
  Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<TTbarGenProducer> ttgenprod;

  //Handles
  Event::Handle<float> h_weight;
  Event::Handle<float> h_weight_pu, h_weight_pu_up, h_weight_pu_down;
  Event::Handle<float> h_eventweight_SR;
  Event::Handle<float> h_dyreco_SR, h_dyreco_1_SR, h_dyreco_2_SR;  
  Event::Handle<float> h_Sigma_phi_1_SR, h_Sigma_phi_2_SR, h_Sigma_phi_SR; 
  Event::Handle<float> h_Delta_phi_1_SR, h_Delta_phi_2_SR, h_Delta_phi_SR; 
  Event::Handle<float> h_eventweight_CR1;
  Event::Handle<float> h_dyreco_CR1, h_dyreco_1_CR1, h_dyreco_2_CR1;  
  Event::Handle<float> h_Sigma_phi_1_CR1, h_Sigma_phi_2_CR1, h_Sigma_phi_CR1; 
  Event::Handle<float> h_Delta_phi_1_CR1, h_Delta_phi_2_CR1, h_Delta_phi_CR1; 
  Event::Handle<float> h_eventweight_CR2;
  Event::Handle<float> h_dyreco_CR2, h_dyreco_1_CR2, h_dyreco_2_CR2;  
  Event::Handle<float> h_Sigma_phi_1_CR2, h_Sigma_phi_2_CR2, h_Sigma_phi_CR2; 
  Event::Handle<float> h_Delta_phi_1_CR2, h_Delta_phi_2_CR2, h_Delta_phi_CR2; 

  // uhh2::Event::Handle<float> h_xi_gen;
  // uhh2::Event::Handle<float> h_mtt_gen;
  // uhh2::Event::Handle<float> h_DeltaY_gen;

  // ttbar reconstruction
  Event::Handle<float> h_chi2;
  Event::Handle<bool> h_is_zprime_reconstructed_chi2;
  Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;

  Event::Handle<bool> h_is_zprime_reconstructed_correctmatch;
  Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateCorrectMatch;

  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

  // DNN multiclass output hist
  std::unique_ptr<Hists> h_MulticlassNN_output;

  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_Inclusive_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_Inclusive_SR;


  // ================ SR ==================================================================================================================================================================================================================
  //muon and ele systematics
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_0_500_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_0_350_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_350_500_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_500_750_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_750_1000_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_1000_1500_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_1500Inf_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_0_700_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_700_900_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_900Inf_SR;

  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_0_500_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_0_350_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_350_500_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_500_750_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_750_1000_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_1000_1500_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_1500Inf_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_0_700_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_700_900_SR;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_900Inf_SR;



  // ================ CR1 ==================================================================================================================================================================================================================
  //muon and ele systematics
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_0_500_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_0_350_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_350_500_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_500_750_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_750_1000_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_1000_1500_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_1500Inf_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_0_700_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_700_900_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_900Inf_CR1;

  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_0_500_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_0_350_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_350_500_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_500_750_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_750_1000_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_1000_1500_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_1500Inf_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_0_700_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_700_900_CR1;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_900Inf_CR1;

  

  // ================ CR2 ==================================================================================================================================================================================================================
  //muon and electron systematics
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_0_500_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_0_350_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_350_500_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_500_750_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_750_1000_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_1000_1500_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_1500Inf_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_0_700_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_700_900_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_SystVariations_900Inf_CR2;

  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_0_500_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_0_350_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_350_500_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_500_750_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_750_1000_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_1000_1500_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_1500Inf_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_0_700_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_700_900_CR2;
  std::unique_ptr<Hists> h_DeltaY_reco_PDFVariations_900Inf_CR2;

  // ================ CR2 ends ==================================================================================================================================================================================================================


  // Configuration
  bool isMC, ishotvr, isdeepAK8;
  //bool isEleTriggerMeasurement;
  string Sys_PU, Prefiring_direction, Sys_TopPt_a, Sys_TopPt_b;
  TString sample;
  int runnr_oldtriggers = 299368;

  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;
  bool isMuon, isElectron, isEFT;
  bool isPhoton;
  TString year;

  TH2F *ratio_hist_muon;
  TH2F *ratio_hist_ele;


  Event::Handle<float> h_Ak4_j1_E;
  Event::Handle<float> h_Ak4_j1_eta;
  Event::Handle<float> h_Ak4_j1_m;
  Event::Handle<float> h_Ak4_j1_phi;
  Event::Handle<float> h_Ak4_j1_pt;
  Event::Handle<float> h_Ak4_j1_deepjetbscore;

  Event::Handle<float> h_Ak4_j2_E;
  Event::Handle<float> h_Ak4_j2_eta;
  Event::Handle<float> h_Ak4_j2_m;
  Event::Handle<float> h_Ak4_j2_phi;
  Event::Handle<float> h_Ak4_j2_pt;
  Event::Handle<float> h_Ak4_j2_deepjetbscore;

  Event::Handle<float> h_Ak4_j3_E;
  Event::Handle<float> h_Ak4_j3_eta;
  Event::Handle<float> h_Ak4_j3_m;
  Event::Handle<float> h_Ak4_j3_phi;
  Event::Handle<float> h_Ak4_j3_pt;
  Event::Handle<float> h_Ak4_j3_deepjetbscore;

  Event::Handle<float> h_Ak4_j4_E;
  Event::Handle<float> h_Ak4_j4_eta;
  Event::Handle<float> h_Ak4_j4_m;
  Event::Handle<float> h_Ak4_j4_phi;
  Event::Handle<float> h_Ak4_j4_pt;
  Event::Handle<float> h_Ak4_j4_deepjetbscore;

  Event::Handle<float> h_Ak4_j5_E;
  Event::Handle<float> h_Ak4_j5_eta;
  Event::Handle<float> h_Ak4_j5_m;
  Event::Handle<float> h_Ak4_j5_phi;
  Event::Handle<float> h_Ak4_j5_pt;
  Event::Handle<float> h_Ak4_j5_deepjetbscore;

  Event::Handle<float> h_E;
  Event::Handle<float> h_eta;
  Event::Handle<float> h_phi;
  Event::Handle<float> h_pt;

  Event::Handle<float> h_MET_phi;
  Event::Handle<float> h_MET_pt;

  Event::Handle<float> h_Mu_E;
  Event::Handle<float> h_Mu_eta;
  Event::Handle<float> h_Mu_phi;
  Event::Handle<float> h_Mu_pt;

  Event::Handle<float> h_N_Ak4;

  Event::Handle<float> h_Ak8_j1_E;
  Event::Handle<float> h_Ak8_j1_eta;
  Event::Handle<float> h_Ak8_j1_mSD;
  Event::Handle<float> h_Ak8_j1_phi;
  Event::Handle<float> h_Ak8_j1_pt;
  Event::Handle<float> h_Ak8_j1_tau21;
  Event::Handle<float> h_Ak8_j1_tau32;

  Event::Handle<float> h_Ak8_j2_E;
  Event::Handle<float> h_Ak8_j2_eta;
  Event::Handle<float> h_Ak8_j2_mSD;
  Event::Handle<float> h_Ak8_j2_phi;
  Event::Handle<float> h_Ak8_j2_pt;
  Event::Handle<float> h_Ak8_j2_tau21;
  Event::Handle<float> h_Ak8_j2_tau32;

  Event::Handle<float> h_Ak8_j3_E;
  Event::Handle<float> h_Ak8_j3_eta;
  Event::Handle<float> h_Ak8_j3_mSD;
  Event::Handle<float> h_Ak8_j3_phi;
  Event::Handle<float> h_Ak8_j3_pt;
  Event::Handle<float> h_Ak8_j3_tau21;
  Event::Handle<float> h_Ak8_j3_tau32;

  Event::Handle<float> h_N_Ak8;


  Event::Handle<std::vector<tensorflow::Tensor> > h_NNoutput;
  Event::Handle<double> h_NNoutput0;
  Event::Handle<double> h_NNoutput1;
  Event::Handle<double> h_NNoutput2;

  std::unique_ptr<NeuralNetworkModule> NNModule;
};

void ZprimeAnalysisModule_applyNN::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag + "_Skimming";
    mytag = tag+"_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicHists(ctx,mytag));
  }
}

void ZprimeAnalysisModule_applyNN::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag + "_Skimming";
  mytag = tag+"_General";
  HFolder(mytag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

ZprimeAnalysisModule_applyNN::ZprimeAnalysisModule_applyNN(uhh2::Context& ctx){
  // Print out for running locally
  // debug = true;
  debug = false;

  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }

  // Configuration
  isMC =      (ctx.get("dataset_type") == "MC");

  ishotvr =   (ctx.get("is_hotvr") == "true");
  isdeepAK8 = (ctx.get("is_deepAK8") == "true");
  TString mode = "hotvr";
  if(isdeepAK8) mode = "deepAK8";
  string tmp = ctx.get("dataset_version");
  sample = tmp;
  isUL16preVFP  = (ctx.get("dataset_version").find("UL16preVFP")  != std::string::npos);
  isUL16postVFP = (ctx.get("dataset_version").find("UL16postVFP") != std::string::npos);
  isUL17        = (ctx.get("dataset_version").find("UL17")        != std::string::npos);
  isUL18        = (ctx.get("dataset_version").find("UL18")        != std::string::npos);
  if(isUL16preVFP)  year = "UL16preVFP";
  if(isUL16postVFP) year = "UL16postVFP";
  if(isUL17) year = "UL17";
  if(isUL18) year = "UL18";
  isPhoton = (ctx.get("dataset_version").find("SinglePhoton") != std::string::npos);
  isMuon = false; isElectron = false; isEFT=false;
  if(ctx.get("channel") == "muon")     isMuon = true;
  if(ctx.get("channel") == "electron") isElectron = true;
  if(ctx.get("sample") == "eft")       isEFT = true;
  
  
  // Access to gen-level particles
  if(isMC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen"));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");



  // Lepton IDs
  ElectronId eleID_low  = ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp80);
  ElectronId eleID_high = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp80);
  MuonId     muID_low   = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  MuonId     muID_high  = MuonID(Muon::CutBasedIdGlobalHighPt);

  // Pt thresholds for low/high lepton categories
  double electron_pt_low;
  if(isUL17){ // UL17 ele trigger threshold is higher than other eras
    electron_pt_low = 38.; 
  }
  else{
    electron_pt_low = 35.;
  }
  double muon_pt_low(30.);
  double electron_pt_high(120.);
  double muon_pt_high(55.);

  // Lepton cleaners and id for low and high categories
  const MuonId muonID_low(AndId<Muon>(PtEtaCut(muon_pt_low, 2.4), muID_low));
  const ElectronId electronID_low(AndId<Electron>(PtEtaSCCut(electron_pt_low, 2.5), eleID_low));
  const MuonId muonID_high(AndId<Muon>(PtEtaCut(muon_pt_high, 2.4), muID_high));
  const ElectronId electronID_high(AndId<Electron>(PtEtaSCCut(electron_pt_high, 2.5), eleID_high));

  muon_cleaner_low.reset(new MuonCleaner(muonID_low));
  electron_cleaner_low.reset(new ElectronCleaner(electronID_low));
  muon_cleaner_high.reset(new MuonCleaner(muonID_high));
  electron_cleaner_high.reset(new ElectronCleaner(electronID_high));

  // TTbar reconstruction and chi2 discriminator
  double chi2_max(30.);

  // Lepton triggers
  string trigger_mu_A, trigger_mu_B, trigger_mu_C, trigger_mu_D, trigger_mu_E, trigger_mu_F;
  string trigger_A, trigger_B;
  string trigger_ph_A;

  if(isMuon){//muon channel
    if(isUL17){
      trigger_mu_A = "HLT_IsoMu27_v*";}
    else{
      trigger_mu_A = "HLT_IsoMu24_v*";}
    trigger_mu_B = "HLT_IsoTkMu24_v*";
    trigger_mu_C = "HLT_Mu50_v*";
    trigger_mu_D = "HLT_TkMu50_v*";
    trigger_mu_E = "HLT_OldMu100_v*";
    trigger_mu_F = "HLT_TkMu100_v*";
  }
  if(isElectron){//electron channel
    trigger_B = "HLT115_CaloIdVT_GsfTrkIdT_v*";
    if(isUL16preVFP || isUL16postVFP){
      trigger_A = "HLT27_WPTight_Gsf_v*";}
    if(isUL17){
      trigger_A = "HLT35_WPTight_Gsf_v*";}
    if(isUL18){
      trigger_A = "HLT32_WPTight_Gsf_v*";}
    if(isUL16preVFP || isUL16postVFP){
      trigger_ph_A = "HLT_Photon175_v*";}
    else{
      trigger_ph_A = "HLT_Photon200_v*";}
  }

  // Top tagged jet ID: "and" combination of the HOTVR top tagger (WP 0.8, min_mSD, max_mSD, min_pT) and a tau32 cut on the groomed subjets (WP 0.56)
  const TopJetId toptagID = AndId<TopJet>(HOTVRTopTag(0.8, 140.0, 220.0, 50.0), Tau32Groomed(0.56));

  // Systematics
  Sys_PU = ctx.get("Sys_PU");
  Prefiring_direction = ctx.get("Sys_prefiring");
  Sys_TopPt_a = ctx.get("Systematic_TopPt_a");
  Sys_TopPt_b = ctx.get("Systematic_TopPt_b");

  // b-tagged jet ID parameters
  BTag::algo btag_algo = BTag::DEEPJET; // algorithm used, e.g. DeepCSV, DeepJet, CSVv2, etc.
  BTag::wp btag_wp = BTag::WP_MEDIUM;   // algorithm working point, e.g. loose, medium, tight
  JetId id_btag = BTag(btag_algo, btag_wp);

  // TopPt reweighting parameters for functional form w = exp(a - (b pt_gen)) 
  double a_toppt = 0.0615;  // parameter a
  double b_toppt = -0.0005; // parameter b

  // Modules
  LumiWeight_module.reset(new MCLumiWeight(ctx));
  PUWeight_module.reset(new MCPileupReweight(ctx, Sys_PU));
  TopPtReweight_module.reset(new TopPtReweighting(ctx, a_toppt, b_toppt, Sys_TopPt_a, Sys_TopPt_b, ""));
  
  MCScale_module.reset(new MCScaleVariation(ctx));
  hadronic_top.reset(new HadronicTop(ctx));
  // sf_toptag.reset(new HOTVRScaleFactor(ctx, toptagID, ctx.get("Sys_TopTag", "nominal"), "HadronicTop", "TopTagSF", "HOTVRTopTagSFs"));
  sf_toptag.reset(new TopTagScaleFactor(ctx));
  sf_topmistag.reset(new TopMistagScaleFactor(ctx));
  NLOCorrections_module.reset(new NLOCorrections(ctx));
  ps_weights.reset(new PSWeights(ctx));

  // b-tagging SFs
  sf_btagging.reset(new MCBTagDiscriminantReweighting(ctx, BTag::algo::DEEPJET, "CHS_matched"));

  // set lepton scale factors: see UHH2/common/include/LeptonScaleFactors.h
  sf_muon_iso_stat_low.reset(new uhh2::MuonIsoScaleFactors_stat(ctx, Muon::Selector::PFIsoTight, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_stat_low.reset(new uhh2::MuonIdScaleFactors_stat(ctx, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_stat_high.reset(new uhh2::MuonIdScaleFactors_stat(ctx, Muon::Selector::CutBasedIdGlobalHighPt, true));
  sf_muon_trigger_stat_low.reset(new uhh2::MuonTriggerScaleFactors_stat(ctx, false, true));
  sf_muon_trigger_stat_high.reset(new uhh2::MuonTriggerScaleFactors_stat(ctx, true, false));
  sf_muon_iso_syst_low.reset(new uhh2::MuonIsoScaleFactors_syst(ctx, Muon::Selector::PFIsoTight, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_syst_low.reset(new uhh2::MuonIdScaleFactors_syst(ctx, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_syst_high.reset(new uhh2::MuonIdScaleFactors_syst(ctx, Muon::Selector::CutBasedIdGlobalHighPt, true));
  sf_muon_trigger_syst_low.reset(new uhh2::MuonTriggerScaleFactors_syst(ctx, false, true));
  sf_muon_trigger_syst_high.reset(new uhh2::MuonTriggerScaleFactors_syst(ctx, true, false));

  sf_muon_reco.reset(new MuonRecoSF(ctx));
  sf_ele_id_low.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_iso_V2_wp80, true));
  sf_ele_id_high.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_noIso_V2_wp80, true));
  sf_ele_reco.reset(new uhh2::ElectronRecoScaleFactors(ctx, false, true));

  sf_ele_trigger.reset( new uhh2::ElecTriggerSF(ctx, "central", "eta_ptbins", year) );

  // dummies (needed to aviod set value errors)
  sf_muon_iso_stat_low_dummy.reset(new uhh2::MuonIsoScaleFactors_stat(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_id_stat_dummy.reset(new uhh2::MuonIdScaleFactors_stat(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_trigger_stat_dummy.reset(new uhh2::MuonTriggerScaleFactors_stat(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_iso_syst_low_dummy.reset(new uhh2::MuonIsoScaleFactors_syst(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_id_syst_dummy.reset(new uhh2::MuonIdScaleFactors_syst(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_trigger_syst_dummy.reset(new uhh2::MuonTriggerScaleFactors_syst(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));

  sf_ele_id_dummy.reset(new uhh2::ElectronIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_ele_reco_dummy.reset(new uhh2::ElectronRecoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));

  // Selection modules
  Chi2_selection.reset(new Chi2Cut(ctx, 0., chi2_max));
  TwoDCut_selection_low1.reset(new TwoDCut(0.3, 10.));
  TTbarMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
  Chi2CandidateMatched_selection.reset(new Chi2CandidateMatchedSelection(ctx));
  ZprimeTopTag_selection.reset(new ZprimeTopTagSelection(ctx));
  HEM_selection.reset(new HEMSelection(ctx)); // HEM issue in 2018, veto on leptons and jets

  



  // Top Taggers
  TopTaggerHOTVR.reset(new HOTVRTopTagger(ctx));
  TopTaggerDeepAK8.reset(new DeepAK8TopTagger(ctx));




  


  // ttbar reconstruction candidate builder
  CandidateBuilder.reset(new ZprimeCandidateBuilder(ctx, mode));

  // ttbar reconstruction discriminators
  Chi2DiscriminatorZprime.reset(new ZprimeChi2Discriminator(ctx)); // chi2 discriminator
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  CorrectMatchDiscriminatorZprime.reset(new ZprimeCorrectMatchDiscriminator(ctx)); // correctMatch discriminator
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
  h_BestZprimeCandidateCorrectMatch = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestCorrectMatch");




  // DNN Categorization modules
  Variables_module.reset(new Variables_NN(ctx, mode)); // variables for NN
  // VariablesEFTSR_module.reset(new Variables_EFT_SR(ctx, mode)); // variables for EFT SR
  // VariablesEFTCR1_module.reset(new Variables_EFT_CR1(ctx, mode)); // variables for EFT CR1
  // VariablesEFTCR2_module.reset(new Variables_EFT_CR2(ctx, mode)); // variables for EFT CR2

  // Variations for DeltaY variables?
  h_DeltaY_reco_SystVariations_Inclusive_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_Inclusive_SR"));
  h_DeltaY_reco_PDFVariations_Inclusive_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_Inclusive_SR"));
  if (debug) cout << "[ZprimeAnalysisModule_applyNN - DEBUG] DeltaY_reco_PDFVariations_Inclusive_SR created successfully!" << endl;



  // Event output handles
  h_weight = ctx.declare_event_output<float> ("weight");
  h_weight_pu = ctx.get_handle<float>("weight_pu");
  h_weight_pu_up = ctx.get_handle<float>("weight_pu_up");
  h_weight_pu_down = ctx.get_handle<float>("weight_pu_down");

  // Reconstructed ttbar system chi2 handle
  h_chi2 = ctx.declare_event_output<float> ("rec_chi2");

  // SR variable handles
  h_eventweight_SR = ctx.declare_event_output<float> ("eventweight");

  h_dyreco_SR = ctx.declare_event_output<float>("dyreco_SR");
  h_Sigma_phi_SR = ctx.declare_event_output<float>("Sigma_phi_SR");
  h_Delta_phi_SR = ctx.declare_event_output<float>("Delta_phi_SR");

  h_Delta_phi_1_SR = ctx.declare_event_output<float>("Delta_phi_1_SR");
  h_Delta_phi_2_SR = ctx.declare_event_output<float>("Delta_phi_2_SR");
  
  h_dyreco_1_SR = ctx.declare_event_output<float>("dyreco_1_SR");
  h_dyreco_2_SR = ctx.declare_event_output<float>("dyreco_2_SR");

  h_Sigma_phi_1_SR=ctx.declare_event_output<float>("Sigma_phi_1_SR");
  h_Sigma_phi_2_SR=ctx.declare_event_output<float>("Sigma_phi_2_SR");


  // CR1 variable handles
  h_eventweight_CR1 = ctx.declare_event_output<float> ("eventweight");

  h_dyreco_CR1 = ctx.declare_event_output<float>("dyreco_CR1");
  h_Sigma_phi_CR1 = ctx.declare_event_output<float>("Sigma_phi_CR1");
  h_Delta_phi_CR1 = ctx.declare_event_output<float>("Delta_phi_CR1");

  h_Delta_phi_1_CR1 = ctx.declare_event_output<float>("Delta_phi_1_CR1");
  h_Delta_phi_2_CR1 = ctx.declare_event_output<float>("Delta_phi_2_CR1");
  
  h_dyreco_1_CR1 = ctx.declare_event_output<float>("dyreco_1_CR1");
  h_dyreco_2_CR1 = ctx.declare_event_output<float>("dyreco_2_CR1");

  h_Sigma_phi_1_CR1=ctx.declare_event_output<float>("Sigma_phi_1_CR1");
  h_Sigma_phi_2_CR1=ctx.declare_event_output<float>("Sigma_phi_2_CR1");

  
  // CR2 variable handles
  h_eventweight_CR2 = ctx.declare_event_output<float> ("eventweight");

  h_dyreco_CR2 = ctx.declare_event_output<float>("dyreco_CR2");
  h_Sigma_phi_CR2 = ctx.declare_event_output<float>("Sigma_phi_CR2");
  h_Delta_phi_CR2 = ctx.declare_event_output<float>("Delta_phi_CR2");

  h_Delta_phi_1_CR2 = ctx.declare_event_output<float>("Delta_phi_1_CR2");
  h_Delta_phi_2_CR2 = ctx.declare_event_output<float>("Delta_phi_2_CR2");
  
  h_dyreco_1_CR2 = ctx.declare_event_output<float>("dyreco_1_CR2");
  h_dyreco_2_CR2 = ctx.declare_event_output<float>("dyreco_2_CR2");

  h_Sigma_phi_1_CR2=ctx.declare_event_output<float>("Sigma_phi_1_CR2");
  h_Sigma_phi_2_CR2=ctx.declare_event_output<float>("Sigma_phi_2_CR2");


  
  // Strings that define a histogram folder
  vector<string> histogram_tags = {
  "Weights_Init", 
  "Weights_HEM", 
  "Weights_PU", "Weights_Lumi", "Weights_TopPt", "Weights_MCScale", "Weights_Prefiring", "Weights_PS", 
  "Weights_TopTag_SF", "Weights_TopMistag_SF",

  "IdEle_SF", 
  "IsoMuon_SF", 
  "IdMuon_SF", 
  "RecoEle_SF", 
  "RecoMuon_SF", 
  
  "TriggerMuon_SF", 

  "BeforeBtagSF", "AfterBtagSF", "AfterCustomBtagSF",
  "NLOCorrections",
  "TriggerEle_SF", 
  "AfterBaseline",

  "TTbarCandidate", 
  "CorrectMatchDiscriminator", "Chi2Discriminator",
  "PassChi2", "PassChi2_Merged", "PassChi2_Resolved", "FailChi2",
  "SR", "SR_PassesChi2", "SR_PassesChi2_Merged", "SR_PassesChi2_Resolved", "SR_FailsChi2",
  "CR1","CR1_PassesChi2",
  "CR2","CR2_PassesChi2",
  };
  book_histograms(ctx, histogram_tags);


  // Multiclass NN output histograms
  h_MulticlassNN_output.reset(new ZprimeSemiLeptonicMulticlassNNHists(ctx, "MulticlassNN"));
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] MulticlassNNHists created successfully!" << endl;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Handling dataset sample name, includig EFT samples //////////////////////////////////////////////////////////////
  if(isMC){
    TString sample_name = "";
    vector<TString> names = {"MC_EFT_Mttbar_0-700_UL17", "MC_EFT_Mttbar_700-900_UL17", "MC_EFT_Mttbar_900-Inf_UL17","ST", "WJets", "DY", "QCD"};

    for(unsigned int i=0; i<names.size(); i++){
      if( ctx.get("dataset_version").find(names.at(i)) != std::string::npos ) sample_name = names.at(i);
    }
    if( (ctx.get("dataset_version").find("TTToHadronic") != std::string::npos)
      || (ctx.get("dataset_version").find("TTToSemiLeptonic") != std::string::npos)
      || (ctx.get("dataset_version").find("TTTo2L2Nu") != std::string::npos) ) {
      sample_name = "TTbar";
    }
    if( (ctx.get("dataset_version").find("MC_EFT_Mttbar_0-700_UL17") != std::string::npos)
      || (ctx.get("dataset_version").find("MC_EFT_Mttbar_700-900_UL17") != std::string::npos)
      || (ctx.get("dataset_version").find("MC_EFT_Mttbar_900-Inf_UL17") != std::string::npos) ) {
      sample_name = "TTbar_EFT";
    }
    if( (ctx.get("dataset_version").find("WW") != std::string::npos)
      || (ctx.get("dataset_version").find("ZZ") != std::string::npos)
      || (ctx.get("dataset_version").find("WZ") != std::string::npos) ) {
      sample_name = "Diboson";
    }  
    // *** CHANGED ***: set isEFT if sample_name == "TTbar_EFT"
    if(sample_name == "TTbar_EFT") {
      isEFT = true;
    } else {
      isEFT = false;
    }

    // 2D b-tag SF reading with the new logic (EFT or others):
    if(isMuon){
      TFile* f_btag2Dsf_muon = new TFile("/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_muon_"+year+".root");
      if(isEFT){
        ratio_hist_muon = (TH2F*)f_btag2Dsf_muon->Get("N_Jets_vs_HT_TTbar");
      }
      else{
        ratio_hist_muon = (TH2F*)f_btag2Dsf_muon->Get("N_Jets_vs_HT_" + sample_name);
      } 
      ratio_hist_muon->SetDirectory(0);
    }
    else if(!isMuon){
      TFile* f_btag2Dsf_ele = new TFile("/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_electron_"+year+".root");
      if(isEFT){
        ratio_hist_ele = (TH2F*)f_btag2Dsf_ele->Get("N_Jets_vs_HT_TTbar");
      }
      else{
        ratio_hist_ele = (TH2F*)f_btag2Dsf_ele->Get("N_Jets_vs_HT_" + sample_name);
      }
      ratio_hist_ele->SetDirectory(0);
    }
  }

  if(debug) cout << "[ZprimeAnalysisModule_applyNN - DEBUG] About to get handles..." << endl;
  h_Ak4_j1_E   = ctx.get_handle<float>("Ak4_j1_E");
  h_Ak4_j1_eta = ctx.get_handle<float>("Ak4_j1_eta");
  h_Ak4_j1_m   = ctx.get_handle<float>("Ak4_j1_m");
  h_Ak4_j1_phi = ctx.get_handle<float>("Ak4_j1_phi");
  h_Ak4_j1_pt  = ctx.get_handle<float>("Ak4_j1_pt");
  h_Ak4_j1_deepjetbscore  = ctx.get_handle<float>("Ak4_j1_deepjetbscore");

  h_Ak4_j2_E   = ctx.get_handle<float>("Ak4_j2_E");
  h_Ak4_j2_eta = ctx.get_handle<float>("Ak4_j2_eta");
  h_Ak4_j2_m   = ctx.get_handle<float>("Ak4_j2_m");
  h_Ak4_j2_phi = ctx.get_handle<float>("Ak4_j2_phi");
  h_Ak4_j2_pt  = ctx.get_handle<float>("Ak4_j2_pt");
  h_Ak4_j2_deepjetbscore  = ctx.get_handle<float>("Ak4_j2_deepjetbscore");

  h_Ak4_j3_E   = ctx.get_handle<float>("Ak4_j3_E");
  h_Ak4_j3_eta = ctx.get_handle<float>("Ak4_j3_eta");
  h_Ak4_j3_m   = ctx.get_handle<float>("Ak4_j3_m");
  h_Ak4_j3_phi = ctx.get_handle<float>("Ak4_j3_phi");
  h_Ak4_j3_pt  = ctx.get_handle<float>("Ak4_j3_pt");
  h_Ak4_j3_deepjetbscore  = ctx.get_handle<float>("Ak4_j3_deepjetbscore");

  h_Ak4_j4_E   = ctx.get_handle<float>("Ak4_j4_E");
  h_Ak4_j4_eta = ctx.get_handle<float>("Ak4_j4_eta");
  h_Ak4_j4_m   = ctx.get_handle<float>("Ak4_j4_m");
  h_Ak4_j4_phi = ctx.get_handle<float>("Ak4_j4_phi");
  h_Ak4_j4_pt  = ctx.get_handle<float>("Ak4_j4_pt");
  h_Ak4_j4_deepjetbscore  = ctx.get_handle<float>("Ak4_j4_deepjetbscore");

  h_Ak4_j5_E   = ctx.get_handle<float>("Ak4_j5_E");
  h_Ak4_j5_eta = ctx.get_handle<float>("Ak4_j5_eta");
  h_Ak4_j5_m   = ctx.get_handle<float>("Ak4_j5_m");
  h_Ak4_j5_phi = ctx.get_handle<float>("Ak4_j5_phi");
  h_Ak4_j5_pt  = ctx.get_handle<float>("Ak4_j5_pt");
  h_Ak4_j5_deepjetbscore  = ctx.get_handle<float>("Ak4_j5_deepjetbscore");

  h_E    = ctx.get_handle<float>("Ele_E");
  h_eta  = ctx.get_handle<float>("Ele_eta");
  h_phi  = ctx.get_handle<float>("Ele_phi");
  h_pt   = ctx.get_handle<float>("Ele_pt");

  h_MET_phi = ctx.get_handle<float>("MET_phi");
  h_MET_pt = ctx.get_handle<float>("MET_pt");

  h_Mu_E    = ctx.get_handle<float>("Mu_E");
  h_Mu_eta  = ctx.get_handle<float>("Mu_eta");
  h_Mu_phi  = ctx.get_handle<float>("Mu_phi");
  h_Mu_pt   = ctx.get_handle<float>("Mu_pt");

  h_N_Ak4 = ctx.get_handle<float>("N_Ak4");

  h_Ak8_j1_E     = ctx.get_handle<float>("Ak8_j1_E");
  h_Ak8_j1_eta   = ctx.get_handle<float>("Ak8_j1_eta");
  h_Ak8_j1_mSD   = ctx.get_handle<float>("Ak8_j1_mSD");
  h_Ak8_j1_phi   = ctx.get_handle<float>("Ak8_j1_phi");
  h_Ak8_j1_pt    = ctx.get_handle<float>("Ak8_j1_pt");
  h_Ak8_j1_tau21 = ctx.get_handle<float>("Ak8_j1_tau21");
  h_Ak8_j1_tau32 = ctx.get_handle<float>("Ak8_j1_tau32");

  h_Ak8_j2_E     = ctx.get_handle<float>("Ak8_j2_E");
  h_Ak8_j2_eta   = ctx.get_handle<float>("Ak8_j2_eta");
  h_Ak8_j2_mSD   = ctx.get_handle<float>("Ak8_j2_mSD");
  h_Ak8_j2_phi   = ctx.get_handle<float>("Ak8_j2_phi");
  h_Ak8_j2_pt    = ctx.get_handle<float>("Ak8_j2_pt");
  h_Ak8_j2_tau21 = ctx.get_handle<float>("Ak8_j2_tau21");
  h_Ak8_j2_tau32 = ctx.get_handle<float>("Ak8_j2_tau32");

  h_Ak8_j3_E     = ctx.get_handle<float>("Ak8_j3_E");
  h_Ak8_j3_eta   = ctx.get_handle<float>("Ak8_j3_eta");
  h_Ak8_j3_mSD   = ctx.get_handle<float>("Ak8_j3_mSD");
  h_Ak8_j3_phi   = ctx.get_handle<float>("Ak8_j3_phi");
  h_Ak8_j3_pt    = ctx.get_handle<float>("Ak8_j3_pt");
  h_Ak8_j3_tau21 = ctx.get_handle<float>("Ak8_j3_tau21");
  h_Ak8_j3_tau32 = ctx.get_handle<float>("Ak8_j3_tau32");

  h_N_Ak8 = ctx.get_handle<float>("N_Ak8");

  h_NNoutput = ctx.get_handle<std::vector<tensorflow::Tensor>>("NNoutput");
  h_NNoutput0 = ctx.declare_event_output<double>("NNoutput0");
  h_NNoutput1 = ctx.declare_event_output<double>("NNoutput1");
  h_NNoutput2 = ctx.declare_event_output<double>("NNoutput2");
  // cout <<"about to get models" << endl;



  ///////////////////////////////////////////////////////////// LEPTON-SPECIFIC NN SETTINGS /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// DON'T FORGET TO CHANGE! ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //muon
  NNModule.reset( new NeuralNetworkModule(ctx, "/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_muon/model.pb", "/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_muon/model.config.pbtxt"));
  
  //electron
  // NNModule.reset( new NeuralNetworkModule(ctx, "/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_ele/model.pb", "/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_ele/model.config.pbtxt"));



  // Structure Constants Calculator for EFT samples
  if(isEFT){structure_constants_calculator.reset(new StructureConstantsCalculator(ctx));}
  else     {structure_constants_calculator.reset(nullptr);}
  h_structure_constants = ctx.get_handle<std::vector<float>>("structure_constants");

  // // declare GEN inputs only for TTbar samples (including EFT) to avoid missing-branch errors on backgrounds/data
  // std::string dataset_version = ctx.get("dataset_version");
  // const bool is_ttbar_or_eft = (isMC && ((dataset_version.find("TTTo") != std::string::npos) || (dataset_version.find("EFT") != std::string::npos)));
  // const bool use_noac = (ctx.get("noac_apply_event_weight", "false") == "true");
  // const bool declare_gen_branches = is_ttbar_or_eft && use_noac;
  
  // if(declare_gen_branches) {
  //   // Declare GEN branches needed for template method NoAC weight calculations
  //   // Note: EFT samples must have these branches for NoAC weights to work
  //   h_xi_gen     = ctx.declare_event_input<float>("xi_gen");
  //   h_DeltaY_gen = ctx.declare_event_input<float>("DeltaY_gen");
  //   h_mtt_gen    = ctx.declare_event_input<float>("mtt_gen");
  // }

  if(debug) cout << "[ZprimeAnalysisModule_applyNN - DEBUG] Handles created successfully!" << endl;
}

/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool ZprimeAnalysisModule_applyNN::process(uhh2::Event& event){
  static int event_counter = 0;
  event_counter++;

 if(debug)cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
 if(debug) cout << " run.event: " << event.run << ". " << event.event << endl;

 // Process TTbarGen first
 if(isMC && ttgenprod) ttgenprod->process(event);
 if(debug) cout << "[ZprimeAnalysisModule_applyNN] TTbarGen processed: ok" << endl;

  // Initialize reco flags with false
  event.set(h_is_zprime_reconstructed_chi2, false);
  event.set(h_is_zprime_reconstructed_correctmatch, false);
  event.set(h_chi2,-100);
  event.set(h_weight,-100);

  //EFT vars SR
  event.set(h_dyreco_SR,-10);
  event.set(h_Sigma_phi_SR,-10);
  event.set(h_Delta_phi_SR,-10);
  event.set(h_Sigma_phi_1_SR,-10);
  event.set(h_Sigma_phi_2_SR,-10);
  event.set(h_dyreco_1_SR,-10);
  event.set(h_dyreco_2_SR,-10);
  event.set(h_Delta_phi_1_SR,-10);
  event.set(h_Delta_phi_2_SR,-10);
  //EFT vars CR1
  event.set(h_dyreco_CR1,-10); 
  event.set(h_Sigma_phi_CR1,-10);
  event.set(h_Delta_phi_CR1,-10);
  event.set(h_Sigma_phi_1_CR1,-10);
  event.set(h_Sigma_phi_2_CR1,-10);
  event.set(h_dyreco_1_CR1,-10);
  event.set(h_dyreco_2_CR1,-10);
  event.set(h_Delta_phi_1_CR1,-10);
  event.set(h_Delta_phi_2_CR1,-10);
  //EFT CR2
  event.set(h_dyreco_CR2,-10);
  event.set(h_Sigma_phi_CR2,-10);
  event.set(h_Delta_phi_CR2,-10);
  event.set(h_Sigma_phi_1_CR2,-10);
  event.set(h_Sigma_phi_2_CR2,-10);
  event.set(h_dyreco_1_CR2,-10);
  event.set(h_dyreco_2_CR2,-10);
  event.set(h_Delta_phi_1_CR2,-10);
  event.set(h_Delta_phi_2_CR2,-10);
  
  // DNN outputs
  event.set(h_NNoutput0, 0);
  event.set(h_NNoutput1, 0);
  event.set(h_NNoutput2, 0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// Application of various weights to event.weight ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Begin HOTVR module-- don't use, might delete //
  if(ishotvr){
    TopTaggerHOTVR->process(event);
    hadronic_top->process(event);
  }
  else if(isdeepAK8){ 
    TopTaggerDeepAK8->process(event); // Creates DeepAK8TopTag Jet collection based on pT, SDmass, and discriminant scores
    hadronic_top->process(event); // only needed for HOTVR jets-- will likely delete
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] TopTaggerDeepAK8: ok" << endl;
  }

  // Initial weight (gen weight * pileup weight from preselection)
  fill_histograms(event, "Weights_Init");

  // HEM veto for 2018 data and MC on eta and phi of jets and leptons in specific runs
  if(!HEM_selection->passes(event)){
    if(!isMC) return false;
    else event.weight = event.weight*(1-0.64774715284); // see https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis for calculation
  }
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] HEM Selection: passed" << endl;
  fill_histograms(event, "Weights_HEM");

  // pileup weight
  PUWeight_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] PUWeight: ok" << endl;
  fill_histograms(event, "Weights_PU");
  
  // lumi weight
  LumiWeight_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] LumiWeight: ok" << endl;
  fill_histograms(event, "Weights_Lumi");

  // top pt reweighting
  TopPtReweight_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] TopPtReweight: ok" << endl;
  fill_histograms(event, "Weights_TopPt");

  // MC scale
  MCScale_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] MCScale: ok" << endl;
  fill_histograms(event, "Weights_MCScale");

  // Prefiring weights
  if(isMC){
    if (Prefiring_direction == "nominal") event.weight *= event.prefiringWeight;
    else if (Prefiring_direction == "up") event.weight *= event.prefiringWeightUp;
    else if (Prefiring_direction == "down") event.weight *= event.prefiringWeightDown;
  }
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] Prefiring: ok" << endl;
  fill_histograms(event, "Weights_Prefiring");

  // Write PSWeights from genInfo to own branch in output tree
  ps_weights->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] Weights_PS: ok" << endl;
  fill_histograms(event, "Weights_PS");

  // DeepAK8 TopTag SFs
  if(isdeepAK8) sf_toptag->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] TopTag_SF: ok" << endl;
  fill_histograms(event, "Weights_TopTag_SF");

  // DeepAK8 TopMistag SFs
  if(isdeepAK8) sf_topmistag->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] TopMistag_SF: ok" << endl;
  fill_histograms(event, "Weights_TopMistag_SF");


  /////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// Lepton SF's //////////////////////////////////////////

  // categorize muons by pT
  double muon_pt_high(55.);  // muon category pT threshold
  bool muon_is_low = false;  // boolean for low-pT muon category
  bool muon_is_high = false; // boolean for high-pT muon category

  // Set muon booleans based on muon pT threshold
  if(isMuon){
    vector<Muon>* muons = event.muons;
    for(unsigned int i=0; i<muons->size(); i++){
      if(event.muons->at(i).pt() < muon_pt_high){
        muon_is_low = true;}
      else{muon_is_high = true;}
    }
  }
  sort_by_pt<Muon>(*event.muons);


  // categorize electrons by pT
  // (includes some isMuon statements for eleTrigger SF measurement in muon channel)
  double electron_pt_high(120.); // electron category pT threshold
  bool ele_is_low = false;       // boolean for low-pT electron category
  bool ele_is_high = false;      // boolean for high-pT electron category

  // Set electron category booleans based on pT
  // but first remove ECAL-gap electrons
  if(isElectron){
    vector<Electron>* electrons = event.electrons;
    for(unsigned int i=0; i<electrons->size(); i++){
      if(event.electrons->at(i).pt()<=electron_pt_high){
        ele_is_low = true;}
      else{ele_is_high = true;}
    }
    if(debug && event_counter <= 5) 
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] Finished looping over electrons" << endl;
  }
  sort_by_pt<Electron>(*event.electrons);
  

  // electronID SF
  if(isMuon){ sf_ele_id_dummy->process(event);}
  if(isElectron){
    if     (ele_is_low){ sf_ele_id_low->process(event);}  // low-pT eleID SF
    else if(ele_is_high){sf_ele_id_high->process(event);} // high-pT eleID SF
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] electronID SF: ok" << endl;
    fill_histograms(event, "IdEle_SF");
  }

  // muonISO SF, both stat&syst but low-pT ONLY
  if(isMuon){
    if(muon_is_low){
      sf_muon_iso_stat_low->process(event);
      sf_muon_iso_syst_low->process(event);}
    else if(muon_is_high){
      sf_muon_iso_stat_low_dummy->process(event);  // dummy stat SFs for high-pT muons
      sf_muon_iso_syst_low_dummy->process(event);} // dummy syst SFs for high-pT muons
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] muonISO SF: ok" << endl;
    fill_histograms(event, "IsoMuon_SF");
  }
  if(isElectron){ // dummy muonISO SFs for electron channel
    sf_muon_iso_stat_low_dummy->process(event);
    sf_muon_iso_syst_low_dummy->process(event);
  }

  // muonID SF, both stat&syst and high&low-pT
  if(isMuon){
    if(muon_is_low){
      sf_muon_id_stat_low->process(event);
      sf_muon_id_syst_low->process(event);}
    else if(muon_is_high){
      sf_muon_id_stat_high->process(event);
      sf_muon_id_syst_high->process(event);}
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] muonID SF: ok" << endl;
    fill_histograms(event, "IdMuon_SF");
  }
  if(isElectron){ // dummy muonID SFs for electron channel
    sf_muon_id_stat_dummy->process(event);
    sf_muon_id_syst_dummy->process(event);
  }

  // electronRECO SF
  if(isMuon){sf_ele_reco_dummy->process(event);}
  if(isElectron){sf_ele_reco->process(event);
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] electronRECO SF: ok" << endl;
    fill_histograms(event, "RecoEle_SF");
  }

  // muonRECO SF
  sf_muon_reco->process(event);
  fill_histograms(event, "RecoMuon_SF");


  ////////////////////// muonTrigger SF //////////////////////
  if(isMuon){
    if(muon_is_low){ // low-pT
      sf_muon_trigger_stat_low->process(event); // stat SF
      sf_muon_trigger_syst_low->process(event); // syst SF
    }
    if(muon_is_high){ // high-pT
      sf_muon_trigger_stat_high->process(event); // stat SF
      sf_muon_trigger_syst_high->process(event); // syst SF
    }
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] muonTrigger SF: ok" << endl;
    fill_histograms(event, "TriggerMuon_SF");
  }
  if(isElectron){ // dummy muonTrigger SFs for electron channel
    sf_muon_trigger_stat_dummy->process(event);
    sf_muon_trigger_syst_dummy->process(event);
  }


  // Fill histograms before b-tagging SF - used to extract Custom BTag SF in (NJets,HT)
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] before b-tagging" << endl;
  fill_histograms(event, "BeforeBtagSF");
  sf_btagging->process(event); // Apply b-tag Shape SF 
  fill_histograms(event, "AfterBtagSF");
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] b-tag SF: ok" << endl;

  // Apply custom SF to correct for b-tag SF shape effects on NJets/HT (HT := sum of jet pTs)
  if(isMC && isMuon){
    float custom_sf;
    vector<Jet>* jets = event.jets;
    int Njets = jets->size();
    double st_jets = 0.;
    for(const auto & jet : *jets) st_jets += jet.pt();
    custom_sf = ratio_hist_muon->GetBinContent( ratio_hist_muon->GetXaxis()->FindBin(Njets), ratio_hist_muon->GetYaxis()->FindBin(st_jets) );
    event.weight *= custom_sf;
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] custom SF in muon channel: ok" << endl;
  }
  if(isMC && !isMuon){
    float custom_sf;
    vector<Jet>* jets = event.jets;
    int Njets = jets->size();
    double st_jets = 0.;
    for(const auto & jet : *jets) st_jets += jet.pt();
    custom_sf = ratio_hist_ele->GetBinContent( ratio_hist_ele->GetXaxis()->FindBin(Njets), ratio_hist_ele->GetYaxis()->FindBin(st_jets) );
    event.weight *= custom_sf;
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] custom SF in electron channel: ok" << endl;
  }
  fill_histograms(event, "AfterCustomBtagSF");
  
  // Higher Order Corrections: EWK & QCD NLO
  NLOCorrections_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] NLO Correction: ok" << endl;
  fill_histograms(event, "NLOCorrections");


  // electronTrigger SF
  sf_ele_trigger->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] electronTrigger SF: ok" << endl;
  fill_histograms(event, "TriggerEle_SF");


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////// Done applying SF's to baseline selection /////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  fill_histograms(event, "AfterBaseline");

  // build all possible candidates in event
  CandidateBuilder->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] CandidateBuilder: ok" << endl;
  fill_histograms(event, "TTbarCandidate");
  

  // Extract CorrectMatch discriminator from ttbar candidates
  if(isMC){ // exclude EFT samples until we make sure they have the expected gen-collection
    CorrectMatchDiscriminatorZprime->process(event);
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] CorrectMatchDiscriminator: ok" << endl;
    fill_histograms(event, "CorrectMatchDiscriminator");
  }

  // Extract chi2 discriminator from ttbar candidates
  Chi2DiscriminatorZprime->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] Chi2Discriminator: ok" << endl;
  fill_histograms(event, "Chi2Discriminator");



  // Select events whose chi2 candidates have chi2 < 30
  if(Chi2_selection->passes(event)){
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] Chi2 Selection: passed" << endl;
    fill_histograms(event, "PassChi2");
    
    // Merged topology
    if(ZprimeTopTag_selection->passes(event)){ 
      if(debug) cout << "[ZprimeAnalysisModule_applyNN] TopTag Selection: passed" << endl;
      fill_histograms(event, "PassChi2_Merged");   
    }
    // Resolved topology
    else{
      if(debug) cout << "[ZprimeAnalysisModule_applyNN] TopTag Selection: failed" << endl;
      fill_histograms(event, "PassChi2_Resolved"); 
    } 
  }
  else{fill_histograms(event, "FailChi2");}



  // Variables for NN
  sort_by_pt<Jet>(*event.jets);
  Variables_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] Variables_module: ok" << endl;

  // NN module
  NNModule->process(event);
  std::vector<tensorflow::Tensor> NNoutputs = NNModule->GetOutputs();

  if(debug) cout << "[ZprimeAnalysisModule_applyNN] Starting DNN" << endl;
  event.set(h_NNoutput0, (double)(NNoutputs[0].tensor<float, 2>()(0,0)));
  event.set(h_NNoutput1, (double)(NNoutputs[0].tensor<float, 2>()(0,1)));
  event.set(h_NNoutput2, (double)(NNoutputs[0].tensor<float, 2>()(0,2)));
  event.set(h_NNoutput, NNoutputs);
  double out0 = (double)(NNoutputs[0].tensor<float, 2>()(0,0));
  double out1 = (double)(NNoutputs[0].tensor<float, 2>()(0,1));
  double out2 = (double)(NNoutputs[0].tensor<float, 2>()(0,2));
  vector<double> out_event = {out0, out1, out2};

  // Assign each event a DNN-category based on highest score of the three output nodes
  // h_MulticlassNN_output->fill(event);
  double max_score = 0.0;
  for (int i = 0; i < 3; i++){
    if (out_event[i] > max_score){
      max_score = out_event[i];
    }
  }
  if (debug) cout <<"[ZprimeAnalysisModule_applyNN] Done setting scores" <<endl;


  //////////////////////////////////////////////////////////////////
  ///////////////////////// DNN categories /////////////////////////
  ////////// out0=TTbar[SR], out1=ST[CR1], out2=WJets[CR2] /////////
  //////////////////////////////////////////////////////////////////

  // SR
  if(out0 == max_score){
    // SR events, inclusive in event topology
    fill_histograms(event, "SR");
    // VariablesEFTSR_module->process(event);
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] Processed EFT SR Variables module" << endl;

    // cut on chi2 < 30
    if(Chi2_selection->passes(event)){
      fill_histograms(event, "SR_PassesChi2");

      // Cut on event topology 
      if(ZprimeTopTag_selection->passes(event)){
           fill_histograms(event, "SR_PassesChi2_Merged");}   // Merged topology
      else{fill_histograms(event, "SR_PassesChi2_Resolved");} // Resolved topology
    }
    else{fill_histograms(event, "SR_FailsChi2");}
  }

  // CR1
  if( out1 == max_score ){
    fill_histograms(event, "CR1");
    // VariablesEFTCR1_module->process(event);
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] Done EFT CR1" << endl;

    if(Chi2_selection->passes(event)){ 
      fill_histograms(event,"CR1_PassesChi2");
    }
  }
 
  // CR2
  if( out2 == max_score ){
    fill_histograms(event, "CR2");
    // VariablesEFTCR2_module->process(event);
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] Done EFT CR2" << endl;

    if(Chi2_selection->passes(event)){ 
      fill_histograms(event,"CR2_PassesChi2");
    }
  }

  if(debug) cout << "[ZprimeAnalysisModule_applyNN] Done with DNN regions" << endl;


  // Calculate structure constants for EFT weights
  // This accesses EFT weights starting at index 202 in event.genInfo->systweights()
  // and calculates structure constants that can be used to compute weights for any WC values
  // calculates the structure constants for each event.
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] isEFT: " << isEFT << endl;
  if(isEFT){structure_constants_calculator->process(event);}


  // Shows the number of structure constants stored in the event
  // Displays the first few structure constants
  // Shows the constant term (SM point) and a few linear terms
  // Prints for the first 5 EFT events
  
  // Debug output for structure constants (only for first few events)
  if (debug && isEFT) {static int event_counter = 0;
    if(debug) cout << "[ZprimeAnalysisModule_applyNN] Checking structure constants" << endl;
    if (event_counter < 5) {
      // Get the structure constants from the event
      if (event.is_valid(h_structure_constants)) {
        std::vector<float> structure_constants = event.get(h_structure_constants);
        
        std::cout << "[ZprimeAnalysisModule_applyNN] ===== Structure Constants Debug (Event " << event_counter << ") =====" << std::endl;
        std::cout << "[ZprimeAnalysisModule_applyNN] Number of structure constants: " << structure_constants.size() << std::endl;
        
        if (!structure_constants.empty()) {
          // Print first few constants
          std::cout << "[ZprimeAnalysisModule_applyNN] First few constants: ";
          for (size_t i = 0; i < std::min(size_t(10), structure_constants.size()); ++i) {
            std::cout << structure_constants[i] << ", ";
          }
          std::cout << std::endl;
        }
        
        // Increment counter after printing
        event_counter++;
      } else {
        std::cout << "[ZprimeAnalysisModule_applyNN] Structure constants not found in event!" << std::endl;
      }
    }
  }
  if(debug) cout << "[ZprimeAnalysisModule_applyNN] Moving on to next event" << endl;
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule_applyNN)