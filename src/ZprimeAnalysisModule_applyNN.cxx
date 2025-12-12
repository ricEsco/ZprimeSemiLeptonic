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
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSystematicsHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicPDFHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicMulticlassNNHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicCHSMatchHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicMistagHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>
#include <UHH2/ZprimeSemiLeptonic/include/ElecTriggerSF.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopTagScaleFactor.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopMistagScaleFactor.h>

#include <UHH2/common/include/TTbarGen.h>
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

  //NN - DON'T FORGET TO CHANGE!
  //Muon
  // ifstream normfile ("/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_muon/NormInfo.txt", ios::in);
  //Electron
  ifstream normfile ("/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_ele/NormInfo.txt", ios::in);
  
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

  //Only Ele or Mu variables!! DON'T FORGET TO CHANGE!

  // /Muon
  // vector<uhh2::Event::Handle<float>> inputs = {h_Ak4_j1_E, h_Ak4_j1_deepjetbscore, h_Ak4_j1_eta, h_Ak4_j1_m, h_Ak4_j1_phi, h_Ak4_j1_pt, h_Ak4_j2_E, h_Ak4_j2_deepjetbscore, h_Ak4_j2_eta, h_Ak4_j2_m, h_Ak4_j2_phi, h_Ak4_j2_pt, h_Ak4_j3_E, h_Ak4_j3_deepjetbscore, h_Ak4_j3_eta, h_Ak4_j3_m, h_Ak4_j3_phi, h_Ak4_j3_pt, h_Ak4_j4_E, h_Ak4_j4_deepjetbscore, h_Ak4_j4_eta, h_Ak4_j4_m, h_Ak4_j4_phi, h_Ak4_j4_pt, h_Ak4_j5_E, h_Ak4_j5_deepjetbscore, h_Ak4_j5_eta, h_Ak4_j5_m, h_Ak4_j5_phi, h_Ak4_j5_pt, h_Ak8_j1_E, h_Ak8_j1_eta, h_Ak8_j1_mSD, h_Ak8_j1_phi, h_Ak8_j1_pt, h_Ak8_j1_tau21, h_Ak8_j1_tau32, h_Ak8_j2_E, h_Ak8_j2_eta, h_Ak8_j2_mSD, h_Ak8_j2_phi, h_Ak8_j2_pt, h_Ak8_j2_tau21, h_Ak8_j2_tau32, h_Ak8_j3_E, h_Ak8_j3_eta, h_Ak8_j3_mSD, h_Ak8_j3_phi, h_Ak8_j3_pt, h_Ak8_j3_tau21, h_Ak8_j3_tau32, h_MET_phi, h_MET_pt, h_Mu_E, h_Mu_eta, h_Mu_phi, h_Mu_pt, h_N_Ak4, h_N_Ak8}; // in alphabetical order to match NormInfo.txt
  //Electron
  vector<uhh2::Event::Handle<float>> inputs = {h_Ak4_j1_E, h_Ak4_j1_deepjetbscore, h_Ak4_j1_eta, h_Ak4_j1_m, h_Ak4_j1_phi, h_Ak4_j1_pt, h_Ak4_j2_E, h_Ak4_j2_deepjetbscore, h_Ak4_j2_eta, h_Ak4_j2_m, h_Ak4_j2_phi, h_Ak4_j2_pt, h_Ak4_j3_E, h_Ak4_j3_deepjetbscore, h_Ak4_j3_eta, h_Ak4_j3_m, h_Ak4_j3_phi, h_Ak4_j3_pt, h_Ak4_j4_E, h_Ak4_j4_deepjetbscore, h_Ak4_j4_eta, h_Ak4_j4_m, h_Ak4_j4_phi, h_Ak4_j4_pt, h_Ak4_j5_E, h_Ak4_j5_deepjetbscore, h_Ak4_j5_eta, h_Ak4_j5_m, h_Ak4_j5_phi, h_Ak4_j5_pt, h_Ak8_j1_E, h_Ak8_j1_eta, h_Ak8_j1_mSD, h_Ak8_j1_phi, h_Ak8_j1_pt, h_Ak8_j1_tau21, h_Ak8_j1_tau32, h_Ak8_j2_E, h_Ak8_j2_eta, h_Ak8_j2_mSD, h_Ak8_j2_phi, h_Ak8_j2_pt, h_Ak8_j2_tau21, h_Ak8_j2_tau32, h_Ak8_j3_E, h_Ak8_j3_eta, h_Ak8_j3_mSD, h_Ak8_j3_phi, h_Ak8_j3_pt, h_Ak8_j3_tau21, h_Ak8_j3_tau32, h_Ele_E, h_Ele_eta, h_Ele_phi, h_Ele_pt, h_MET_phi, h_MET_pt, h_N_Ak4, h_N_Ak8}; // in alphabetical order to match NormInfo.txt
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
  // bool isEFT;
  
  // Cleaners
  std::unique_ptr<MuonCleaner>     muon_cleaner_low, muon_cleaner_high;
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
  unique_ptr<Variables_EFT_SR> VariablesEFTSR_module;
  unique_ptr<Variables_EFT_CR1> VariablesEFTCR1_module;
  unique_ptr<Variables_EFT_CR2> VariablesEFTCR2_module;


  //Handles
  Event::Handle<bool> h_is_zprime_reconstructed_chi2, h_is_zprime_reconstructed_correctmatch;
  Event::Handle<float> h_chi2;
  Event::Handle<float> h_weight;
  Event::Handle<float> h_weight_pu, h_weight_pu_up, h_weight_pu_down;
  Event::Handle<float> h_eventweight_SR;
  Event::Handle<float> h_dyreco_SR, h_dyreco_1_SR, h_dyreco_2_SR;  
  Event::Handle<float> h_dyreco_1_SR_0_500, h_dyreco_1_SR_0_350, h_dyreco_1_SR_350_500, h_dyreco_1_SR_500_750, h_dyreco_1_SR_750_1000, h_dyreco_1_SR_1000_1500, h_dyreco_1_SR_1500_Inf, h_dyreco_1_SR_0_700, h_dyreco_1_SR_700_900, h_dyreco_1_SR_900_Inf;
  Event::Handle<float> h_dyreco_2_SR_0_500, h_dyreco_2_SR_0_350, h_dyreco_2_SR_350_500, h_dyreco_2_SR_500_750, h_dyreco_2_SR_750_1000, h_dyreco_2_SR_1000_1500, h_dyreco_2_SR_1500_Inf, h_dyreco_2_SR_0_700, h_dyreco_2_SR_700_900, h_dyreco_2_SR_900_Inf;
  Event::Handle<float> h_Sigma_phi_1_SR, h_Sigma_phi_2_SR, h_Sigma_phi_SR; 
  Event::Handle<float> h_Sigma_phi_1_SR_0_500, h_Sigma_phi_1_SR_0_350, h_Sigma_phi_1_SR_350_500, h_Sigma_phi_1_SR_500_750, h_Sigma_phi_1_SR_750_1000, h_Sigma_phi_1_SR_1000_1500, h_Sigma_phi_1_SR_1500_Inf, h_Sigma_phi_1_SR_0_700, h_Sigma_phi_1_SR_700_900, h_Sigma_phi_1_SR_900_Inf;
  Event::Handle<float> h_Delta_phi_1_SR, h_Delta_phi_2_SR, h_Delta_phi_SR; 
  Event::Handle<float> h_Sigma_phi_2_SR_0_500, h_Sigma_phi_2_SR_0_350, h_Sigma_phi_2_SR_350_500, h_Sigma_phi_2_SR_500_750, h_Sigma_phi_2_SR_750_1000, h_Sigma_phi_2_SR_1000_1500, h_Sigma_phi_2_SR_1500_Inf,h_Sigma_phi_2_SR_0_700, h_Sigma_phi_2_SR_700_900, h_Sigma_phi_2_SR_900_Inf;
  Event::Handle<float> h_eventweight_CR1;
  Event::Handle<float> h_dyreco_CR1, h_dyreco_1_CR1, h_dyreco_2_CR1;  
  Event::Handle<float> h_dyreco_1_CR1_0_500, h_dyreco_1_CR1_0_350, h_dyreco_1_CR1_350_500, h_dyreco_1_CR1_500_750, h_dyreco_1_CR1_750_1000, h_dyreco_1_CR1_1000_1500, h_dyreco_1_CR1_1500_Inf, h_dyreco_1_CR1_0_700, h_dyreco_1_CR1_700_900, h_dyreco_1_CR1_900_Inf;
  Event::Handle<float> h_dyreco_2_CR1_0_500, h_dyreco_2_CR1_0_350, h_dyreco_2_CR1_350_500, h_dyreco_2_CR1_500_750, h_dyreco_2_CR1_750_1000, h_dyreco_2_CR1_1000_1500, h_dyreco_2_CR1_1500_Inf, h_dyreco_2_CR1_0_700, h_dyreco_2_CR1_700_900, h_dyreco_2_CR1_900_Inf;
  Event::Handle<float> h_Sigma_phi_1_CR1, h_Sigma_phi_2_CR1, h_Sigma_phi_CR1; 
  Event::Handle<float> h_Sigma_phi_1_CR1_0_500, h_Sigma_phi_1_CR1_0_350, h_Sigma_phi_1_CR1_350_500, h_Sigma_phi_1_CR1_500_750, h_Sigma_phi_1_CR1_750_1000, h_Sigma_phi_1_CR1_1000_1500, h_Sigma_phi_1_CR1_1500_Inf,h_Sigma_phi_1_CR1_0_700, h_Sigma_phi_1_CR1_700_900, h_Sigma_phi_1_CR1_900_Inf;
  Event::Handle<float> h_Delta_phi_1_CR1, h_Delta_phi_2_CR1, h_Delta_phi_CR1; 
  Event::Handle<float> h_Sigma_phi_2_CR1_0_500, h_Sigma_phi_2_CR1_0_350, h_Sigma_phi_2_CR1_350_500, h_Sigma_phi_2_CR1_500_750, h_Sigma_phi_2_CR1_750_1000, h_Sigma_phi_2_CR1_1000_1500, h_Sigma_phi_2_CR1_1500_Inf, h_Sigma_phi_2_CR1_0_700, h_Sigma_phi_2_CR1_700_900, h_Sigma_phi_2_CR1_900_Inf;
  Event::Handle<float> h_eventweight_CR2;
  Event::Handle<float> h_dyreco_CR2, h_dyreco_1_CR2, h_dyreco_2_CR2;  
  Event::Handle<float> h_dyreco_1_CR2_0_500,h_dyreco_1_CR2_0_350, h_dyreco_1_CR2_350_500, h_dyreco_1_CR2_500_750, h_dyreco_1_CR2_750_1000, h_dyreco_1_CR2_1000_1500, h_dyreco_1_CR2_1500_Inf, h_dyreco_1_CR2_0_700, h_dyreco_1_CR2_700_900, h_dyreco_1_CR2_900_Inf;
  Event::Handle<float> h_dyreco_2_CR2_0_500, h_dyreco_2_CR2_0_350, h_dyreco_2_CR2_350_500, h_dyreco_2_CR2_500_750, h_dyreco_2_CR2_750_1000, h_dyreco_2_CR2_1000_1500, h_dyreco_2_CR2_1500_Inf, h_dyreco_2_CR2_0_700, h_dyreco_2_CR2_700_900, h_dyreco_2_CR2_900_Inf;
  Event::Handle<float> h_Sigma_phi_1_CR2, h_Sigma_phi_2_CR2, h_Sigma_phi_CR2; 
  Event::Handle<float> h_Sigma_phi_1_CR2_0_500, h_Sigma_phi_1_CR2_0_350, h_Sigma_phi_1_CR2_350_500, h_Sigma_phi_1_CR2_500_750, h_Sigma_phi_1_CR2_750_1000, h_Sigma_phi_1_CR2_1000_1500, h_Sigma_phi_1_CR2_1500_Inf,h_Sigma_phi_1_CR2_0_700, h_Sigma_phi_1_CR2_700_900, h_Sigma_phi_1_CR2_900_Inf;
  Event::Handle<float> h_Delta_phi_1_CR2, h_Delta_phi_2_CR2, h_Delta_phi_CR2; 
  Event::Handle<float> h_Sigma_phi_2_CR2_0_500, h_Sigma_phi_2_CR2_0_350, h_Sigma_phi_2_CR2_350_500, h_Sigma_phi_2_CR2_500_750, h_Sigma_phi_2_CR2_750_1000, h_Sigma_phi_2_CR2_1000_1500, h_Sigma_phi_2_CR2_1500_Inf, h_Sigma_phi_2_CR2_0_700, h_Sigma_phi_2_CR2_700_900, h_Sigma_phi_2_CR2_900_Inf;


  uhh2::Event::Handle<float> h_xi_gen;
  uhh2::Event::Handle<float> h_mtt_gen;
  uhh2::Event::Handle<float> h_DeltaY_gen;

  
  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;

  // Lumi hists
  std::unique_ptr<Hists> lumihists_Weights_Init, lumihists_Weights_PU, lumihists_Weights_Lumi, lumihists_Weights_TopPt, lumihists_Weights_MCScale, lumihists_Weights_PS, lumihists_Chi2;

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

   //bool isEleTriggerMeasurement;



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
  // debug = true;
  debug = false;
  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }
  // Configuration
  isMC = (ctx.get("dataset_type") == "MC");
  ishotvr = (ctx.get("is_hotvr") == "true");
  isdeepAK8 = (ctx.get("is_deepAK8") == "true");
  TString mode = "hotvr";
  if(isdeepAK8) mode = "deepAK8";
  string tmp = ctx.get("dataset_version");
  sample = tmp;
  isUL16preVFP  = (ctx.get("dataset_version").find("UL16preVFP")  != std::string::npos);
  isUL16postVFP = (ctx.get("dataset_version").find("UL16postVFP") != std::string::npos);
  isUL17        = (ctx.get("dataset_version").find("UL17")        != std::string::npos);
  isUL18        = (ctx.get("dataset_version").find("UL18")        != std::string::npos);
  if(isUL16preVFP) year = "UL16preVFP";
  if(isUL16postVFP) year = "UL16postVFP";
  if(isUL17) year = "UL17";
  if(isUL18) year = "UL18";
  

  isPhoton = (ctx.get("dataset_version").find("SinglePhoton") != std::string::npos);
  // isEleTriggerMeasurement = (ctx.get("isTriggerMeasurement") == "true");

  // Lepton IDs
  ElectronId eleID_low  = ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp80);
  ElectronId eleID_high = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp80);
  MuonId     muID_low   = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  MuonId     muID_high  = MuonID(Muon::CutBasedIdGlobalHighPt);

  double electron_pt_low;
  if(isUL17){
    electron_pt_low = 38.; // UL17 ele trigger threshold is 35 (HLT35WPTight _Gsf) -> be above turn on
  }
  else{
    electron_pt_low = 35.;
  }
  double muon_pt_low(30.);
  double electron_pt_high(120.);
  double muon_pt_high(55.);

  const MuonId muonID_low(AndId<Muon>(PtEtaCut(muon_pt_low, 2.4), muID_low));
  const ElectronId electronID_low(AndId<Electron>(PtEtaSCCut(electron_pt_low, 2.5), eleID_low));
  const MuonId muonID_high(AndId<Muon>(PtEtaCut(muon_pt_high, 2.4), muID_high));
  const ElectronId electronID_high(AndId<Electron>(PtEtaSCCut(electron_pt_high, 2.5), eleID_high));

  muon_cleaner_low.reset(new MuonCleaner(muonID_low));
  electron_cleaner_low.reset(new ElectronCleaner(electronID_low));
  muon_cleaner_high.reset(new MuonCleaner(muonID_high));
  electron_cleaner_high.reset(new ElectronCleaner(electronID_high));

  // Important selection values
  double chi2_max(30.);
  string trigger_mu_A, trigger_mu_B, trigger_mu_C, trigger_mu_D, trigger_mu_E, trigger_mu_F;
  string trigger_A, trigger_B;
  string trigger_ph_A;
  isMuon = false; isElectron = false, isEFT=false;
  if(ctx.get("channel") == "muon") isMuon = true;
  if(ctx.get("sample") == "eft") isEFT = true;
  if(ctx.get("channel") == "electron") isElectron = true;

  if(isMuon){//semileptonic muon channel
    if(isUL17){
      trigger_mu_A = "HLT_IsoMu27_v*";
    }
    else{
      trigger_mu_A = "HLT_IsoMu24_v*";
    }
    trigger_mu_B = "HLT_IsoTkMu24_v*";
    trigger_mu_C = "HLT_Mu50_v*";
    trigger_mu_D = "HLT_TkMu50_v*";
    trigger_mu_E = "HLT_OldMu100_v*";
    trigger_mu_F = "HLT_TkMu100_v*";

  }
  if(isElectron){//semileptonic electron channel
    trigger_B = "HLT115_CaloIdVT_GsfTrkIdT_v*";
    if(isUL16preVFP || isUL16postVFP){
      trigger_A = "HLT27_WPTight_Gsf_v*";
    }
    if(isUL17){
      trigger_A = "HLT35_WPTight_Gsf_v*";
    }
    if(isUL18){
      trigger_A = "HLT32_WPTight_Gsf_v*";
    }
    if(isUL16preVFP || isUL16postVFP){
      trigger_ph_A = "HLT_Photon175_v*";
    }
    else{
      trigger_ph_A = "HLT_Photon200_v*";
    }
  }


  const TopJetId toptagID = AndId<TopJet>(HOTVRTopTag(0.8, 140.0, 220.0, 50.0), Tau32Groomed(0.56));

  Sys_PU = ctx.get("Sys_PU");
  Prefiring_direction = ctx.get("Sys_prefiring");
  Sys_TopPt_a = ctx.get("Systematic_TopPt_a");
  Sys_TopPt_b = ctx.get("Systematic_TopPt_b");

  BTag::algo btag_algo = BTag::DEEPJET;
  BTag::wp btag_wp = BTag::WP_MEDIUM;
  JetId id_btag = BTag(btag_algo, btag_wp);

  double a_toppt = 0.0615; // par a TopPt Reweighting
  double b_toppt = -0.0005; // par b TopPt Reweighting

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
  
  DeltaEta_selection.reset(new DeltaEtaSelection()); // Cut on DeltaEta(j1,j2)<3. to reduce QCD spikes

  Variables_module.reset(new Variables_NN(ctx, mode)); // variables for NN
  VariablesEFTSR_module.reset(new Variables_EFT_SR(ctx, mode)); // variables for EFT SR
  VariablesEFTCR1_module.reset(new Variables_EFT_CR1(ctx, mode)); // variables for EFT CR1
  VariablesEFTCR2_module.reset(new Variables_EFT_CR2(ctx, mode)); // variables for EFT CR2


 //  if(!isEleTriggerMeasurement) SystematicsModule.reset(new ZprimeSemiLeptonicSystematicsModule(ctx));

  // Top Taggers
  TopTaggerHOTVR.reset(new HOTVRTopTagger(ctx));
  TopTaggerDeepAK8.reset(new DeepAK8TopTagger(ctx));

  // TopTags veto
  TopTagVetoSelection.reset(new TopTag_VetoSelection(ctx, mode));

  // Zprime candidate builder
  CandidateBuilder.reset(new ZprimeCandidateBuilder(ctx, mode));

  // Zprime discriminators
  Chi2DiscriminatorZprime.reset(new ZprimeChi2Discriminator(ctx));
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  CorrectMatchDiscriminatorZprime.reset(new ZprimeCorrectMatchDiscriminator(ctx));
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");

  h_chi2 = ctx.declare_event_output<float> ("rec_chi2");
  h_weight = ctx.declare_event_output<float> ("weight");
  h_weight_pu = ctx.get_handle<float>("weight_pu");
  h_weight_pu_up = ctx.get_handle<float>("weight_pu_up");
  h_weight_pu_down = ctx.get_handle<float>("weight_pu_down");
  h_eventweight_SR = ctx.declare_event_output<float> ("eventweight");

  h_dyreco_SR = ctx.declare_event_output<float>("dyreco_SR");
  h_Sigma_phi_SR = ctx.declare_event_output<float>("Sigma_phi_SR");
  h_Delta_phi_SR = ctx.declare_event_output<float>("Delta_phi_SR");

  h_Delta_phi_1_SR = ctx.declare_event_output<float>("Delta_phi_1_SR");
  h_Delta_phi_2_SR = ctx.declare_event_output<float>("Delta_phi_2_SR");
  
  h_dyreco_1_SR = ctx.declare_event_output<float>("dyreco_1_SR");
  h_dyreco_1_SR_0_500 = ctx.declare_event_output<float>("dyreco_1_SR_0_500");
  h_dyreco_1_SR_0_350 = ctx.declare_event_output<float>("dyreco_1_SR_0_350");
  h_dyreco_1_SR_350_500 = ctx.declare_event_output<float>("dyreco_1_SR_350_500");
  h_dyreco_1_SR_500_750 = ctx.declare_event_output<float>("dyreco_1_SR_500_750");
  h_dyreco_1_SR_750_1000 = ctx.declare_event_output<float>("dyreco_1_SR_750_1000");
  h_dyreco_1_SR_1000_1500 = ctx.declare_event_output<float>("dyreco_1_SR_1000_1500");
  h_dyreco_1_SR_1500_Inf = ctx.declare_event_output<float>("dyreco_1_SR_1500_Inf");
  h_dyreco_1_SR_0_700 = ctx.declare_event_output<float>("dyreco_1_SR_0_700");
  h_dyreco_1_SR_700_900 = ctx.declare_event_output<float>("dyreco_1_SR_700_900");
  h_dyreco_1_SR_900_Inf = ctx.declare_event_output<float>("dyreco_1_SR_900_Inf");

  h_dyreco_2_SR = ctx.declare_event_output<float>("dyreco_2_SR");
  h_dyreco_2_SR_0_500 = ctx.declare_event_output<float>("dyreco_2_SR_0_500");
  h_dyreco_2_SR_0_350 = ctx.declare_event_output<float>("dyreco_2_SR_0_350");
  h_dyreco_2_SR_350_500 = ctx.declare_event_output<float>("dyreco_2_SR_350_500");
  h_dyreco_2_SR_500_750 = ctx.declare_event_output<float>("dyreco_2_SR_500_750");
  h_dyreco_2_SR_750_1000 = ctx.declare_event_output<float>("dyreco_2_SR_750_1000");
  h_dyreco_2_SR_1000_1500 = ctx.declare_event_output<float>("dyreco_2_SR_1000_1500");
  h_dyreco_2_SR_1500_Inf = ctx.declare_event_output<float>("dyreco_2_SR_1500_Inf");
  h_dyreco_2_SR_0_700 = ctx.declare_event_output<float>("dyreco_2_SR_0_700");
  h_dyreco_2_SR_700_900 = ctx.declare_event_output<float>("dyreco_2_SR_700_900");
  h_dyreco_2_SR_900_Inf = ctx.declare_event_output<float>("dyreco_2_SR_900_Inf");

  h_Sigma_phi_1_SR=ctx.declare_event_output<float>("Sigma_phi_1_SR");
  h_Sigma_phi_1_SR_0_500=ctx.declare_event_output<float>("Sigma_phi_1_SR_0_500");
  h_Sigma_phi_1_SR_0_350=ctx.declare_event_output<float>("Sigma_phi_1_SR_0_350");
  h_Sigma_phi_1_SR_350_500=ctx.declare_event_output<float>("Sigma_phi_1_SR_350_500");
  h_Sigma_phi_1_SR_500_750=ctx.declare_event_output<float>("Sigma_phi_1_SR_500_750");
  h_Sigma_phi_1_SR_750_1000=ctx.declare_event_output<float>("Sigma_phi_1_SR_750_1000");
  h_Sigma_phi_1_SR_1000_1500=ctx.declare_event_output<float>("Sigma_phi_1_SR_1000_1500");
  h_Sigma_phi_1_SR_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_1_SR_1500_Inf");
  h_Sigma_phi_1_SR_0_700=ctx.declare_event_output<float>("Sigma_phi_1_SR_0_700");
  h_Sigma_phi_1_SR_700_900=ctx.declare_event_output<float>("Sigma_phi_1_SR_700_900");
  h_Sigma_phi_1_SR_900_Inf=ctx.declare_event_output<float>("Sigma_phi_1_SR_900_Inf");

  h_Sigma_phi_2_SR=ctx.declare_event_output<float>("Sigma_phi_2_SR");
  h_Sigma_phi_2_SR_0_500=ctx.declare_event_output<float>("Sigma_phi_2_SR_0_500");
  h_Sigma_phi_2_SR_0_350=ctx.declare_event_output<float>("Sigma_phi_2_SR_0_350");
  h_Sigma_phi_2_SR_350_500=ctx.declare_event_output<float>("Sigma_phi_2_SR_350_500");
  h_Sigma_phi_2_SR_500_750=ctx.declare_event_output<float>("Sigma_phi_2_SR_500_750");
  h_Sigma_phi_2_SR_750_1000=ctx.declare_event_output<float>("Sigma_phi_2_SR_750_1000");
  h_Sigma_phi_2_SR_1000_1500=ctx.declare_event_output<float>("Sigma_phi_2_SR_1000_1500");
  h_Sigma_phi_2_SR_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_2_SR_1500_Inf");
  h_Sigma_phi_2_SR_0_700=ctx.declare_event_output<float>("Sigma_phi_2_SR_0_700");
  h_Sigma_phi_2_SR_700_900=ctx.declare_event_output<float>("Sigma_phi_2_SR_700_900");
  h_Sigma_phi_2_SR_900_Inf=ctx.declare_event_output<float>("Sigma_phi_2_SR_900_Inf");

  
  //CR1//
  
  h_eventweight_CR1 = ctx.declare_event_output<float> ("eventweight");

  h_dyreco_CR1 = ctx.declare_event_output<float>("dyreco_CR1");
  h_Sigma_phi_CR1 = ctx.declare_event_output<float>("Sigma_phi_CR1");
  h_Delta_phi_CR1 = ctx.declare_event_output<float>("Delta_phi_CR1");

  h_Delta_phi_1_CR1 = ctx.declare_event_output<float>("Delta_phi_1_CR1");
  h_Delta_phi_2_CR1 = ctx.declare_event_output<float>("Delta_phi_2_CR1");
  
  h_dyreco_1_CR1 = ctx.declare_event_output<float>("dyreco_1_CR1");
  h_dyreco_1_CR1_0_500 = ctx.declare_event_output<float>("dyreco_1_CR1_0_500");
  h_dyreco_1_CR1_0_350 = ctx.declare_event_output<float>("dyreco_1_CR1_0_350");
  h_dyreco_1_CR1_350_500 = ctx.declare_event_output<float>("dyreco_1_CR1_350_500");
  h_dyreco_1_CR1_500_750 = ctx.declare_event_output<float>("dyreco_1_CR1_500_750");
  h_dyreco_1_CR1_750_1000 = ctx.declare_event_output<float>("dyreco_1_CR1_750_1000");
  h_dyreco_1_CR1_1000_1500 = ctx.declare_event_output<float>("dyreco_1_CR1_1000_1500");
  h_dyreco_1_CR1_1500_Inf = ctx.declare_event_output<float>("dyreco_1_CR1_1500_Inf");
  h_dyreco_1_CR1_0_700 = ctx.declare_event_output<float>("dyreco_1_CR1_0_700");
  h_dyreco_1_CR1_700_900 = ctx.declare_event_output<float>("dyreco_1_CR1_700_900");
  h_dyreco_1_CR1_900_Inf = ctx.declare_event_output<float>("dyreco_1_CR1_900_Inf");

  h_dyreco_2_CR1 = ctx.declare_event_output<float>("dyreco_2_CR1");
  h_dyreco_2_CR1_0_500 = ctx.declare_event_output<float>("dyreco_2_CR1_0_500");
  h_dyreco_2_CR1_0_350 = ctx.declare_event_output<float>("dyreco_2_CR1_0_350");
  h_dyreco_2_CR1_350_500 = ctx.declare_event_output<float>("dyreco_2_CR1_350_500");
  h_dyreco_2_CR1_500_750 = ctx.declare_event_output<float>("dyreco_2_CR1_500_750");
  h_dyreco_2_CR1_750_1000 = ctx.declare_event_output<float>("dyreco_2_CR1_750_1000");
  h_dyreco_2_CR1_1000_1500 = ctx.declare_event_output<float>("dyreco_2_CR1_1000_1500");
  h_dyreco_2_CR1_1500_Inf = ctx.declare_event_output<float>("dyreco_2_CR1_1500_Inf");
  h_dyreco_2_CR1_0_700 = ctx.declare_event_output<float>("dyreco_2_CR1_0_700");
  h_dyreco_2_CR1_700_900 = ctx.declare_event_output<float>("dyreco_2_CR1_700_900");
  h_dyreco_2_CR1_900_Inf = ctx.declare_event_output<float>("dyreco_2_CR1_900_Inf");

  h_Sigma_phi_1_CR1=ctx.declare_event_output<float>("Sigma_phi_1_CR1");
  h_Sigma_phi_1_CR1_0_500=ctx.declare_event_output<float>("Sigma_phi_1_CR1_0_500");
  h_Sigma_phi_1_CR1_0_350=ctx.declare_event_output<float>("Sigma_phi_1_CR1_0_350");
  h_Sigma_phi_1_CR1_350_500=ctx.declare_event_output<float>("Sigma_phi_1_CR1_350_500");
  h_Sigma_phi_1_CR1_500_750=ctx.declare_event_output<float>("Sigma_phi_1_CR1_500_750");
  h_Sigma_phi_1_CR1_750_1000=ctx.declare_event_output<float>("Sigma_phi_1_CR1_750_1000");
  h_Sigma_phi_1_CR1_1000_1500=ctx.declare_event_output<float>("Sigma_phi_1_CR1_1000_1500");
  h_Sigma_phi_1_CR1_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_1_CR1_1500_Inf");
  h_Sigma_phi_1_CR1_0_700=ctx.declare_event_output<float>("Sigma_phi_1_CR1_0_700");
  h_Sigma_phi_1_CR1_700_900=ctx.declare_event_output<float>("Sigma_phi_1_CR1_700_900");
  h_Sigma_phi_1_CR1_900_Inf=ctx.declare_event_output<float>("Sigma_phi_1_CR1_900_Inf");

  h_Sigma_phi_2_CR1=ctx.declare_event_output<float>("Sigma_phi_2_CR1");
  h_Sigma_phi_2_CR1_0_500=ctx.declare_event_output<float>("Sigma_phi_2_CR1_0_500");
  h_Sigma_phi_2_CR1_0_350=ctx.declare_event_output<float>("Sigma_phi_2_CR1_0_350");
  h_Sigma_phi_2_CR1_350_500=ctx.declare_event_output<float>("Sigma_phi_2_CR1_350_500");
  h_Sigma_phi_2_CR1_500_750=ctx.declare_event_output<float>("Sigma_phi_2_CR1_500_750");
  h_Sigma_phi_2_CR1_750_1000=ctx.declare_event_output<float>("Sigma_phi_2_CR1_750_1000");
  h_Sigma_phi_2_CR1_1000_1500=ctx.declare_event_output<float>("Sigma_phi_2_CR1_1000_1500");
  h_Sigma_phi_2_CR1_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_2_CR1_1500_Inf");
  h_Sigma_phi_2_CR1_0_700=ctx.declare_event_output<float>("Sigma_phi_2_CR1_0_700");
  h_Sigma_phi_2_CR1_700_900=ctx.declare_event_output<float>("Sigma_phi_2_CR1_700_900");
  h_Sigma_phi_2_CR1_900_Inf=ctx.declare_event_output<float>("Sigma_phi_2_CR1_900_Inf");

  
  //CR2///
  
  h_eventweight_CR2 = ctx.declare_event_output<float> ("eventweight");

  h_dyreco_CR2 = ctx.declare_event_output<float>("dyreco_CR2");
  h_Sigma_phi_CR2 = ctx.declare_event_output<float>("Sigma_phi_CR2");
  h_Delta_phi_CR2 = ctx.declare_event_output<float>("Delta_phi_CR2");

  h_Delta_phi_1_CR2 = ctx.declare_event_output<float>("Delta_phi_1_CR2");
  h_Delta_phi_2_CR2 = ctx.declare_event_output<float>("Delta_phi_2_CR2");
  
  h_dyreco_1_CR2 = ctx.declare_event_output<float>("dyreco_1_CR2");
  h_dyreco_1_CR2_0_500 = ctx.declare_event_output<float>("dyreco_1_CR2_0_500");
  h_dyreco_1_CR2_0_350 = ctx.declare_event_output<float>("dyreco_1_CR2_0_350");
  h_dyreco_1_CR2_350_500 = ctx.declare_event_output<float>("dyreco_1_CR2_350_500");
  h_dyreco_1_CR2_500_750 = ctx.declare_event_output<float>("dyreco_1_CR2_500_750");
  h_dyreco_1_CR2_750_1000 = ctx.declare_event_output<float>("dyreco_1_CR2_750_1000");
  h_dyreco_1_CR2_1000_1500 = ctx.declare_event_output<float>("dyreco_1_CR2_1000_1500");
  h_dyreco_1_CR2_1500_Inf = ctx.declare_event_output<float>("dyreco_1_CR2_1500_Inf");
  h_dyreco_1_CR2_0_700 = ctx.declare_event_output<float>("dyreco_1_CR2_0_700");
  h_dyreco_1_CR2_700_900 = ctx.declare_event_output<float>("dyreco_1_CR2_700_900");
  h_dyreco_1_CR2_900_Inf = ctx.declare_event_output<float>("dyreco_1_CR2_900_Inf");
  h_dyreco_1_CR2_0_700= ctx.declare_event_output<float>("dyreco_1_CR2_0_700");
  h_dyreco_1_CR2_700_900= ctx.declare_event_output<float>("dyreco_1_CR2_700_900");
  h_dyreco_1_CR2_900_Inf= ctx.declare_event_output<float>("dyreco_1_CR2_900_Inf");

  h_dyreco_2_CR2 = ctx.declare_event_output<float>("dyreco_2_CR2");
  h_dyreco_2_CR2_0_500 = ctx.declare_event_output<float>("dyreco_2_CR2_0_500");
  h_dyreco_2_CR2_0_350 = ctx.declare_event_output<float>("dyreco_2_CR2_0_350");
  h_dyreco_2_CR2_350_500 = ctx.declare_event_output<float>("dyreco_2_CR2_350_500");
  h_dyreco_2_CR2_500_750 = ctx.declare_event_output<float>("dyreco_2_CR2_500_750");
  h_dyreco_2_CR2_750_1000 = ctx.declare_event_output<float>("dyreco_2_CR2_750_1000");
  h_dyreco_2_CR2_1000_1500 = ctx.declare_event_output<float>("dyreco_2_CR2_1000_1500");
  h_dyreco_2_CR2_1500_Inf = ctx.declare_event_output<float>("dyreco_2_CR2_1500_Inf");
  h_dyreco_2_CR2_0_700= ctx.declare_event_output<float>("dyreco_2_CR2_0_700");
  h_dyreco_2_CR2_700_900= ctx.declare_event_output<float>("dyreco_2_CR2_700_900");
  h_dyreco_2_CR2_900_Inf= ctx.declare_event_output<float>("dyreco_2_CR2_900_Inf");


  h_Sigma_phi_1_CR2=ctx.declare_event_output<float>("Sigma_phi_1_CR2");
  h_Sigma_phi_1_CR2_0_500=ctx.declare_event_output<float>("Sigma_phi_1_CR2_0_500");
  h_Sigma_phi_1_CR2_0_350=ctx.declare_event_output<float>("Sigma_phi_1_CR2_0_350");
  h_Sigma_phi_1_CR2_350_500=ctx.declare_event_output<float>("Sigma_phi_1_CR2_350_500");
  h_Sigma_phi_1_CR2_500_750=ctx.declare_event_output<float>("Sigma_phi_1_CR2_500_750");
  h_Sigma_phi_1_CR2_750_1000=ctx.declare_event_output<float>("Sigma_phi_1_CR2_750_1000");
  h_Sigma_phi_1_CR2_1000_1500=ctx.declare_event_output<float>("Sigma_phi_1_CR2_1000_1500");
  h_Sigma_phi_1_CR2_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_1_CR2_1500_Inf");
  h_Sigma_phi_1_CR2_0_700=ctx.declare_event_output<float>("Sigma_phi_1_CR2_0_700");
  h_Sigma_phi_1_CR2_700_900=ctx.declare_event_output<float>("Sigma_phi_1_CR2_700_900"); 
  h_Sigma_phi_1_CR2_900_Inf=ctx.declare_event_output<float>("Sigma_phi_1_CR2_900_Inf");

  h_Sigma_phi_2_CR2=ctx.declare_event_output<float>("Sigma_phi_2_CR2");
  h_Sigma_phi_2_CR2_0_500=ctx.declare_event_output<float>("Sigma_phi_2_CR2_0_500");
  h_Sigma_phi_2_CR2_0_350=ctx.declare_event_output<float>("Sigma_phi_2_CR2_0_350");
  h_Sigma_phi_2_CR2_350_500=ctx.declare_event_output<float>("Sigma_phi_2_CR2_350_500");
  h_Sigma_phi_2_CR2_500_750=ctx.declare_event_output<float>("Sigma_phi_2_CR2_500_750");
  h_Sigma_phi_2_CR2_750_1000=ctx.declare_event_output<float>("Sigma_phi_2_CR2_750_1000");
  h_Sigma_phi_2_CR2_1000_1500=ctx.declare_event_output<float>("Sigma_phi_2_CR2_1000_1500");
  h_Sigma_phi_2_CR2_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_2_CR2_1500_Inf");
  h_Sigma_phi_2_CR2_0_700=ctx.declare_event_output<float>("Sigma_phi_2_CR2_0_700");
  h_Sigma_phi_2_CR2_700_900=ctx.declare_event_output<float>("Sigma_phi_2_CR2_700_900");
  h_Sigma_phi_2_CR2_900_Inf=ctx.declare_event_output<float>("Sigma_phi_2_CR2_900_Inf");


  
  if(debug) cout << "[DEBUG] About to create CHSMatchHists..." << endl;
  h_CHSMatchHists.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch"));
  if(debug) cout << "[DEBUG] CHSMatchHists created successfully!" << endl;

  sel_1btag.reset(new NJetSelection(1, -1, id_btag));
  sel_2btag.reset(new NJetSelection(2,-1, id_btag));

  if (debug) cout << "[DEBUG] About to create DeltaY_reco_SystVariations_Inclusive_SR..." << endl;
  h_DeltaY_reco_SystVariations_Inclusive_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_Inclusive_SR"));
  if (debug) cout << "[DEBUG] DeltaY_reco_SystVariations_Inclusive_SR created successfully!" << endl;
  h_DeltaY_reco_PDFVariations_Inclusive_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_Inclusive_SR"));
  if (debug) cout << "[DEBUG] DeltaY_reco_PDFVariations_Inclusive_SR created successfully!" << endl;

  // ================ SR ==================================================================================================================================================================================================================
  
  // h_DeltaY_reco_SystVariations_0_500_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_0_500_SR"));
  // h_DeltaY_reco_SystVariations_0_350_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_0_350_SR"));
  // h_DeltaY_reco_SystVariations_350_500_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_350_500_SR"));
  // h_DeltaY_reco_SystVariations_500_750_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_500_750_SR"));
  // h_DeltaY_reco_SystVariations_750_1000_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_750_1000_SR"));
  // h_DeltaY_reco_SystVariations_1000_1500_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_1000_1500_SR"));
  // h_DeltaY_reco_SystVariations_1500Inf_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_1500Inf_SR"));
  // h_DeltaY_reco_SystVariations_0_700_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_0_700_SR"));
  // h_DeltaY_reco_SystVariations_700_900_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_700_900_SR"));
  // h_DeltaY_reco_SystVariations_900Inf_SR.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_900Inf_SR"));

  if(debug) cout << "[DEBUG] DeltaY_reco_SystVariations_0_500_SR created successfully!" << endl;

  // h_DeltaY_reco_PDFVariations_0_500_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_0_500_SR"));
  // h_DeltaY_reco_PDFVariations_0_350_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_0_350_SR"));
  // h_DeltaY_reco_PDFVariations_350_500_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_350_500_SR"));
  // h_DeltaY_reco_PDFVariations_500_750_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_500_750_SR"));
  // h_DeltaY_reco_PDFVariations_750_1000_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_750_1000_SR"));
  // h_DeltaY_reco_PDFVariations_1000_1500_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_1000_1500_SR"));
  // h_DeltaY_reco_PDFVariations_1500Inf_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_1500Inf_SR"));
  // h_DeltaY_reco_PDFVariations_0_700_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_0_700_SR"));
  // h_DeltaY_reco_PDFVariations_700_900_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_700_900_SR"));
  // h_DeltaY_reco_PDFVariations_900Inf_SR.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_900Inf_SR"));
  // ================ SR ends ==================================================================================================================================================================================================================
  
  if(debug) cout << "[DEBUG] DeltaY_reco_SystVariations_0_500_CR1 created successfully!" << endl;

  // ================ CR1 ==================================================================================================================================================================================================================
  
  // h_DeltaY_reco_SystVariations_0_500_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_0_500_CR1"));
  // h_DeltaY_reco_SystVariations_0_350_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_0_350_CR1"));
  // h_DeltaY_reco_SystVariations_350_500_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_350_500_CR1"));
  // h_DeltaY_reco_SystVariations_500_750_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_500_750_CR1"));
  // h_DeltaY_reco_SystVariations_750_1000_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_750_1000_CR1"));
  // h_DeltaY_reco_SystVariations_1000_1500_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_1000_1500_CR1"));
  // h_DeltaY_reco_SystVariations_1500Inf_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_1500Inf_CR1"));
  // h_DeltaY_reco_SystVariations_0_700_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_0_700_CR1"));
  // h_DeltaY_reco_SystVariations_700_900_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_700_900_CR1"));
  // h_DeltaY_reco_SystVariations_900Inf_CR1.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_900Inf_CR1"));
  
  
  // h_DeltaY_reco_PDFVariations_0_500_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_0_500_CR1"));
  // h_DeltaY_reco_PDFVariations_0_350_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_0_350_CR1"));
  // h_DeltaY_reco_PDFVariations_350_500_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_350_500_CR1"));
  // h_DeltaY_reco_PDFVariations_500_750_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_500_750_CR1"));
  // h_DeltaY_reco_PDFVariations_750_1000_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_750_1000_CR1"));
  // h_DeltaY_reco_PDFVariations_1000_1500_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_1000_1500_CR1"));
  // h_DeltaY_reco_PDFVariations_1500Inf_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_1500Inf_CR1"));
  // h_DeltaY_reco_PDFVariations_0_700_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_0_700_CR1"));
  // h_DeltaY_reco_PDFVariations_700_900_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_700_900_CR1"));
  // h_DeltaY_reco_PDFVariations_900Inf_CR1.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_900Inf_CR1"));


  // ================ CR1 ends ==================================================================================================================================================================================================================


  // ================ CR2 ==================================================================================================================================================================================================================
  // h_DeltaY_reco_SystVariations_0_500_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_0_500_CR2"));
  // h_DeltaY_reco_SystVariations_0_350_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_0_350_CR2"));
  // h_DeltaY_reco_SystVariations_350_500_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_350_500_CR2"));
  // h_DeltaY_reco_SystVariations_500_750_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_500_750_CR2"));
  // h_DeltaY_reco_SystVariations_750_1000_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_750_1000_CR2"));
  // h_DeltaY_reco_SystVariations_1000_1500_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_1000_1500_CR2"));
  // h_DeltaY_reco_SystVariations_1500Inf_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_1500Inf_CR2"));
  // h_DeltaY_reco_SystVariations_0_700_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_0_700_CR2"));
  // h_DeltaY_reco_SystVariations_700_900_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_700_900_CR2"));
  // h_DeltaY_reco_SystVariations_900Inf_CR2.reset(new ZprimeSemiLeptonicSystematicsHists(ctx, "DeltaY_reco_SystVariations_900Inf_CR2"));


  // h_DeltaY_reco_PDFVariations_0_500_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_0_500_CR2"));
  // h_DeltaY_reco_PDFVariations_0_350_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_0_350_CR2"));
  // h_DeltaY_reco_PDFVariations_350_500_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_350_500_CR2"));
  // h_DeltaY_reco_PDFVariations_500_750_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_500_750_CR2"));
  // h_DeltaY_reco_PDFVariations_750_1000_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_750_1000_CR2"));
  // h_DeltaY_reco_PDFVariations_1000_1500_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_1000_1500_CR2"));
  // h_DeltaY_reco_PDFVariations_1500Inf_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_1500Inf_CR2"));
  // h_DeltaY_reco_PDFVariations_0_700_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_0_700_CR2"));
  // h_DeltaY_reco_PDFVariations_700_900_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_700_900_CR2"));
  // h_DeltaY_reco_PDFVariations_900Inf_CR2.reset(new ZprimeSemiLeptonicPDFHists(ctx, "DeltaY_reco_PDFVariations_900Inf_CR2"));


  // ================ CR2 ends ==================================================================================================================================================================================================================

  
  
  vector<string> histogram_tags = {
  "Initial","Weights_Init", "Weights_HEM", "Weights_PU", "Weights_Lumi", "Weights_TopPt", "Weights_MCScale", "Weights_Prefiring", "Weights_TopTag_SF", "Weights_TopMistag_SF", "Weights_PS", 
  // "NLOCorrections","TwoDCut_low1","IdMuon_SF", "IdEle_SF", "IsoMuon_SF",
  // "RecoEle_SF", "MuonReco_SF", "TriggerMuon_SF", "BeforeBtagSF", "AfterBtagSF", "AfterCustomBtagSF", "TriggerEle_SF", "NNInputsBeforeReweight", "TopTagVeto", "DeltaEtaCut",
  "AfterChi2", 
  // "AfterBaseline", "Chi2_passes","Chi2_withTopTag", "Chi2_noTopTag","Chi2_inverse",
  "DNN_output0_nochi2",
  "DNN_output0",
  // "DNN_output1","DNN_output2","DNN_output1_chi2","DNN_output2_chi2","DNN_output0_TopTag", "DNN_output0_NoTopTag",
  // "DeltaY_reco_1500Inf_SR" ,"DeltaY_reco_1000_1500_SR" ,"DeltaY_reco_750_1000_SR" ,"DeltaY_reco_500_750_SR", "DeltaY_reco_0_350_SR", "DeltaY_reco_350_500_SR","DeltaY_reco_0_500_SR", "DeltaY_reco_0_700_SR", "DeltaY_reco_700_900_SR", "DeltaY_reco_900Inf_SR",
  // "DeltaY_reco_1500Inf_CR1" ,"DeltaY_reco_1000_1500_CR1" ,"DeltaY_reco_750_1000_CR1" ,"DeltaY_reco_500_750_CR1","DeltaY_reco_0_350_CR1", "DeltaY_reco_350_500_CR1", "DeltaY_reco_0_500_CR1","DeltaY_reco_0_700_CR1", "DeltaY_reco_700_900_CR1", "DeltaY_reco_900Inf_CR1", 
  // "DeltaY_reco_1500Inf_CR2" ,"DeltaY_reco_1000_1500_CR2" ,"DeltaY_reco_750_1000_CR2" ,"DeltaY_reco_500_750_CR2", "DeltaY_reco_0_350_CR2", "DeltaY_reco_350_500_CR2","DeltaY_reco_0_500_CR2", "DeltaY_reco_0_700_CR2", "DeltaY_reco_700_900_CR2", "DeltaY_reco_900Inf_CR2",
};

  if(debug) cout << "[DEBUG] About to book histograms..." << endl;
  book_histograms(ctx, histogram_tags);
  if (debug) cout << "[DEBUG] After book histograms" << endl;

  h_MulticlassNN_output.reset(new ZprimeSemiLeptonicMulticlassNNHists(ctx, "MulticlassNN"));
  if(debug) cout << "[DEBUG] MulticlassNNHists created successfully!" << endl;

  // lumihists_Weights_Init.reset(new LuminosityHists(ctx, "Lumi_Weights_Init"));
  // lumihists_Weights_PU.reset(new LuminosityHists(ctx, "Lumi_Weights_PU"));
  // lumihists_Weights_Lumi.reset(new LuminosityHists(ctx, "Lumi_Weights_Lumi"));
  // lumihists_Weights_TopPt.reset(new LuminosityHists(ctx, "Lumi_Weights_TopPt"));
  // lumihists_Weights_MCScale.reset(new LuminosityHists(ctx, "Lumi_Weights_MCScale"));
  // lumihists_Weights_PS.reset(new LuminosityHists(ctx, "Lumi_Weights_PS"));
  // lumihists_Chi2.reset(new LuminosityHists(ctx, "Lumi_Chi2"));
  
  // *** CHANGED ***
  // bool isEFT = false; // default false

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
    // if(sample_name == "TTbar_EFT") {
    //   isEFT = true;
    // } else {
    //   isEFT = false;
    // }
    if (debug)cout << "is it EFT? " << isEFT << endl;
  
    // 2D b-tag SF reading with the new logic (EFT or others):
    if(isMuon){
      TFile* f_btag2Dsf_muon = new TFile("/data/dust/user/deleokse/RunII_106_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_muon_"+year+".root");
      if(isEFT){
        ratio_hist_muon = (TH2F*)f_btag2Dsf_muon->Get("N_Jets_vs_HT_TTbar");
      }
      else{
        ratio_hist_muon = (TH2F*)f_btag2Dsf_muon->Get("N_Jets_vs_HT_" + sample_name);
      } 
      ratio_hist_muon->SetDirectory(0);
    }
    else if(!isMuon){
      TFile* f_btag2Dsf_ele = new TFile("/data/dust/user/deleokse/RunII_106_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_electron_"+year+".root");
      if(isEFT){
        ratio_hist_ele = (TH2F*)f_btag2Dsf_ele->Get("N_Jets_vs_HT_TTbar");
      }
      else{
        ratio_hist_ele = (TH2F*)f_btag2Dsf_ele->Get("N_Jets_vs_HT_" + sample_name);
      }
      ratio_hist_ele->SetDirectory(0);
    }
  }

  if(debug) cout << "[DEBUG] About to get handles..." << endl;
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

   //Only Ele or Mu variables!! DON'T FORGET TO CHANGE!
  //muon
  // NNModule.reset( new NeuralNetworkModule(ctx, "/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_muon/model.pb", "/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_muon/model.config.pbtxt"));
  
  //electron
  NNModule.reset( new NeuralNetworkModule(ctx, "/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_ele/model.pb", "/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/KerasNN/NN_DeepAK8_UL17_ele/model.config.pbtxt"));

  // Structure Constants Calculator for EFT
  // This calculates structure constants for each event using EFT weights
  // with the correct mapping from configurations to weight indices
  if(isEFT){
    structure_constants_calculator.reset(new StructureConstantsCalculator(ctx));
  }
  else{
    structure_constants_calculator.reset(nullptr);
  }
  // structure_constants_calculator.reset(new StructureConstantsCalculator(ctx));
  h_structure_constants = ctx.get_handle<std::vector<float>>("structure_constants");

  // declare GEN inputs only for TTbar samples (including EFT) to avoid missing-branch errors on backgrounds/data
  std::string dataset_version = ctx.get("dataset_version");
  const bool is_ttbar_or_eft = (isMC && ((dataset_version.find("TTTo") != std::string::npos) || (dataset_version.find("EFT") != std::string::npos)));
  const bool use_noac = (ctx.get("noac_apply_event_weight", "false") == "true");
  const bool declare_gen_branches = is_ttbar_or_eft && use_noac;
  
  if(declare_gen_branches) {
    // Declare GEN branches needed for template method NoAC weight calculations
    // Note: EFT samples must have these branches for NoAC weights to work
    h_xi_gen     = ctx.declare_event_input<float>("xi_gen");
    h_DeltaY_gen = ctx.declare_event_input<float>("DeltaY_gen");
    h_mtt_gen    = ctx.declare_event_input<float>("mtt_gen");
  }

  if(debug) cout << "[DEBUG] Handles created successfully!" << endl;
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
  // Initialize reco flags with false
  event.set(h_is_zprime_reconstructed_chi2, false);
  event.set(h_is_zprime_reconstructed_correctmatch, false);
  event.set(h_chi2,-100);
  event.set(h_weight,-100);
  //EFT vars SR
  event.set(h_Sigma_phi_SR,-10);
  event.set(h_Delta_phi_SR,-10);
  event.set(h_Sigma_phi_1_SR,-10);
  event.set(h_Sigma_phi_2_SR,-10);
  event.set(h_dyreco_1_SR,-10);
  event.set(h_dyreco_2_SR,-10);
  event.set(h_Delta_phi_1_SR,-10);
  event.set(h_Delta_phi_2_SR,-10);
  event.set(h_Sigma_phi_1_SR_0_500,-10);
  event.set(h_Sigma_phi_1_SR_0_350,-10);
  event.set(h_Sigma_phi_1_SR_350_500,-10);
  event.set(h_Sigma_phi_1_SR_500_750,-10);
  event.set(h_Sigma_phi_1_SR_750_1000,-10);
  event.set(h_Sigma_phi_1_SR_1000_1500,-10);
  event.set(h_Sigma_phi_1_SR_1500_Inf,-10);
  event.set(h_Sigma_phi_1_SR_0_700,-10);
  event.set(h_Sigma_phi_1_SR_700_900,-10);
  event.set(h_Sigma_phi_1_SR_900_Inf,-10);
  event.set(h_Sigma_phi_2_SR_0_500,-10);
  event.set(h_Sigma_phi_2_SR_0_350,-10);
  event.set(h_Sigma_phi_2_SR_350_500,-10);
  event.set(h_Sigma_phi_2_SR_500_750,-10);
  event.set(h_Sigma_phi_2_SR_500_750,-10);
  event.set(h_Sigma_phi_2_SR_750_1000,-10);
  event.set(h_Sigma_phi_2_SR_1000_1500,-10);
  event.set(h_Sigma_phi_2_SR_1500_Inf,-10);
  event.set(h_Sigma_phi_2_SR_0_700,-10);
  event.set(h_Sigma_phi_2_SR_700_900,-10);
  event.set(h_Sigma_phi_2_SR_900_Inf,-10);
  event.set(h_dyreco_1_SR_0_500,-10);
  event.set(h_dyreco_1_SR_0_350,-10);
  event.set(h_dyreco_1_SR_350_500,-10);
  event.set(h_dyreco_1_SR_500_750,-10);
  event.set(h_dyreco_1_SR_750_1000,-10);
  event.set(h_dyreco_1_SR_1000_1500,-10);
  event.set(h_dyreco_1_SR_1500_Inf,-10);
  event.set(h_dyreco_1_SR_0_700,-10);
  event.set(h_dyreco_1_SR_700_900,-10);
  event.set(h_dyreco_1_SR_900_Inf,-10);
  event.set(h_dyreco_2_SR_0_500,-10);
  event.set(h_dyreco_2_SR_0_350,-10);
  event.set(h_dyreco_2_SR_350_500,-10);
  event.set(h_dyreco_2_SR_500_750,-10);
  event.set(h_dyreco_2_SR_750_1000,-10);
  event.set(h_dyreco_2_SR_1000_1500,-10);
  event.set(h_dyreco_2_SR_1500_Inf,-10);
  event.set(h_dyreco_2_SR_0_700,-10);
  event.set(h_dyreco_2_SR_700_900,-10);
  event.set(h_dyreco_2_SR_900_Inf,-10);
  event.set(h_dyreco_SR,-10); 
  if(debug) cout << " set CR1 vars " << endl;
 
  //EFT vars CR1
  event.set(h_Sigma_phi_CR1,-10);
  event.set(h_Delta_phi_CR1,-10);
  event.set(h_Sigma_phi_1_CR1,-10);
  event.set(h_Sigma_phi_2_CR1,-10);
  event.set(h_dyreco_1_CR1,-10);
  event.set(h_dyreco_2_CR1,-10);
  event.set(h_Delta_phi_1_CR1,-10);
  event.set(h_Delta_phi_2_CR1,-10);
  event.set(h_Sigma_phi_1_CR1_0_500,-10);
  event.set(h_Sigma_phi_1_CR1_0_350,-10);
  event.set(h_Sigma_phi_1_CR1_350_500,-10);
  event.set(h_Sigma_phi_1_CR1_500_750,-10);
  event.set(h_Sigma_phi_1_CR1_750_1000,-10);
  event.set(h_Sigma_phi_1_CR1_1000_1500,-10);
  event.set(h_Sigma_phi_1_CR1_1500_Inf,-10);
  event.set(h_Sigma_phi_2_CR1_0_500,-10);
  event.set(h_Sigma_phi_2_CR1_0_350,-10);
  event.set(h_Sigma_phi_2_CR1_350_500,-10);
  event.set(h_Sigma_phi_2_CR1_500_750,-10);
  event.set(h_Sigma_phi_2_CR1_750_1000,-10);
  event.set(h_Sigma_phi_2_CR1_1000_1500,-10);
  event.set(h_Sigma_phi_2_CR1_1500_Inf,-10);
  event.set(h_dyreco_1_CR1_0_500,-10);
  event.set(h_dyreco_1_CR1_0_350,-10);
  event.set(h_dyreco_1_CR1_350_500,-10);
  event.set(h_dyreco_1_CR1_500_750,-10);
  event.set(h_dyreco_1_CR1_750_1000,-10);
  event.set(h_dyreco_1_CR1_1000_1500,-10);
  event.set(h_dyreco_1_CR1_1500_Inf,-10);
  event.set(h_dyreco_2_CR1_0_500,-10);
  event.set(h_dyreco_2_CR1_0_350,-10);
  event.set(h_dyreco_2_CR1_350_500,-10);
  event.set(h_dyreco_2_CR1_500_750,-10);
  event.set(h_dyreco_2_CR1_750_1000,-10);
  event.set(h_dyreco_2_CR1_1000_1500,-10);
  event.set(h_dyreco_2_CR1_1500_Inf,-10);
  event.set(h_dyreco_CR1,-10); 
  event.set(h_Sigma_phi_2_CR1_0_700,-10);
  event.set(h_Sigma_phi_2_CR1_700_900,-10);
  event.set(h_Sigma_phi_2_CR1_900_Inf,-10);
  event.set(h_dyreco_1_CR1_0_700,-10);
  event.set(h_dyreco_1_CR1_700_900,-10);
  event.set(h_dyreco_1_CR1_900_Inf,-10);
  event.set(h_dyreco_2_CR1_0_700,-10);
  event.set(h_dyreco_2_CR1_700_900,-10);
  event.set(h_dyreco_2_CR1_900_Inf,-10);
  event.set(h_Sigma_phi_1_CR1_0_700,-10);
  event.set(h_Sigma_phi_1_CR1_700_900,-10);
  event.set(h_Sigma_phi_1_CR1_900_Inf,-10);
  
  if(debug) cout << " set CR2 vars " << endl;

  //EFT CR2
  event.set(h_Sigma_phi_CR2,-10);
  if(debug) cout << " set Sigma " << endl;
  event.set(h_Delta_phi_CR2,-10);
  if(debug) cout << " set delta " << endl;
  event.set(h_Sigma_phi_1_CR2,-10);
  if(debug) cout << " set Sigma phi 1 " << endl;
  event.set(h_Sigma_phi_2_CR2,-10);
  if(debug) cout << " set Sigma  phi 2" << endl;
  event.set(h_dyreco_1_CR2,-10);
  if(debug) cout << " set dyreco 1 " << endl;
  event.set(h_dyreco_2_CR2,-10);
  if(debug) cout << " set dyreco 2 " << endl;
  event.set(h_Delta_phi_1_CR2,-10);
  if(debug) cout << " set delta phi 1 " << endl;
  event.set(h_Delta_phi_2_CR2,-10);
  if(debug) cout << " set delta phi 2 " << endl;
  event.set(h_Sigma_phi_1_CR2_0_500,-10);
  event.set(h_Sigma_phi_1_CR2_0_350,-10);
  event.set(h_Sigma_phi_1_CR2_350_500,-10);
  if(debug) cout << " set Sigma_phi_1_CR2_0_500 " << endl;
  event.set(h_Sigma_phi_1_CR2_500_750,-10);
  if(debug) cout << " set Sigma_phi_1_CR2_500_750 " << endl;
  event.set(h_Sigma_phi_1_CR2_750_1000,-10);
  if(debug) cout << " set Sigma_phi_1_CR2_750_1000 " << endl;
  event.set(h_Sigma_phi_1_CR2_1000_1500,-10);
  if(debug) cout << " set Sigma_phi_1_CR2_1000_1500 " << endl;
  event.set(h_Sigma_phi_1_CR2_1500_Inf,-10);
  if(debug) cout << " set Sigma_phi_1_CR2_1500_Inf " << endl;
  event.set(h_Sigma_phi_2_CR2_0_500,-10);
  event.set(h_Sigma_phi_2_CR2_0_350,-10);
  event.set(h_Sigma_phi_2_CR2_350_500,-10);
  if(debug) cout << " set Sigma_phi_2_CR2_0_500 " << endl;
  event.set(h_Sigma_phi_2_CR2_500_750,-10);
  if(debug) cout << " set Sigma_phi_2_CR2_500_750 " << endl;
  event.set(h_Sigma_phi_2_CR2_750_1000,-10);

  event.set(h_Sigma_phi_2_CR2_1000_1500,-10);
  event.set(h_Sigma_phi_2_CR2_1500_Inf,-10);
  event.set(h_dyreco_1_CR2_0_500,-10);
  event.set(h_dyreco_1_CR2_0_350,-10);
  event.set(h_dyreco_1_CR2_350_500,-10);
  event.set(h_dyreco_1_CR2_500_750,-10);
  event.set(h_dyreco_1_CR2_750_1000,-10);
  event.set(h_dyreco_1_CR2_1000_1500,-10);
  event.set(h_dyreco_1_CR2_1500_Inf,-10);
  if(debug) cout << " set dyreco 1 CR2 " << endl;
  event.set(h_dyreco_2_CR2_0_500,-10);
  event.set(h_dyreco_2_CR2_0_350,-10);
  event.set(h_dyreco_2_CR2_350_500,-10);
  event.set(h_dyreco_2_CR2_500_750,-10);
  event.set(h_dyreco_2_CR2_750_1000,-10);
  event.set(h_dyreco_2_CR2_1000_1500,-10);
  event.set(h_dyreco_2_CR2_1500_Inf,-10);
  event.set(h_dyreco_CR2,-10);  
  if(debug) cout << " done with dyreco CR2 " << endl;
  event.set(h_Sigma_phi_2_CR2_0_700,-10);
  event.set(h_Sigma_phi_2_CR2_700_900,-10);
  event.set(h_Sigma_phi_2_CR2_900_Inf,-10);
  if(debug) cout << " done with Sigma phi 2 CR2 " << endl;
  event.set(h_dyreco_1_CR2_0_700,-10);
  event.set(h_dyreco_1_CR2_700_900,-10);
  event.set(h_dyreco_1_CR2_900_Inf,-10);
  if(debug) cout << " done with dyreco 1 CR2 " << endl;
  event.set(h_dyreco_2_CR2_0_700,-10);
  event.set(h_dyreco_2_CR2_700_900,-10); 
  event.set(h_dyreco_2_CR2_900_Inf,-10);
  if(debug) cout << " done with dyreco 2 CR2 " << endl;
  event.set(h_Sigma_phi_1_CR2_0_700,-10);
  event.set(h_Sigma_phi_1_CR2_700_900,-10);
  event.set(h_Sigma_phi_1_CR2_900_Inf,-10);

  //////////////end EFT vars ///////////
  
  event.set(h_NNoutput0, 0);
  event.set(h_NNoutput1, 0);
  event.set(h_NNoutput2, 0);

  
  //Tagger
  if(ishotvr){
    TopTaggerHOTVR->process(event);
    hadronic_top->process(event);
  }else if(isdeepAK8){
    TopTaggerDeepAK8->process(event);
    hadronic_top->process(event);
  }
  if(debug) cout<<"Top Tagger ok"<<endl;

  // fill_histograms(event, "Weights_Init");
  //Setting low and high pt points
  // double muon_pt_high(55.);
  // bool muon_is_low = false;
  // bool muon_is_high = false;
  // if(debug)  cout<<"sorting leptons"<<endl;
  // if(isMuon){
  //   vector<Muon>* muons = event.muons;
  //   for(unsigned int i=0; i<muons->size(); i++){
  //     if(event.muons->at(i).pt()<=muon_pt_high){
  //       muon_is_low = true;
  //     }else{
  //       muon_is_high = true;
  //     }
  //   }
  // }
  // sort_by_pt<Muon>(*event.muons);

  // double electron_pt_high(120.);
  // bool ele_is_low = false;
  // bool ele_is_high = false;

  // if(isElectron){
  //   vector<Electron>* electrons = event.electrons;
  //   for(unsigned int i=0; i<electrons->size(); i++){
  //     if(event.electrons->at(i).pt()<=electron_pt_high){
  //       ele_is_low = true;
  //     }else{
  //       ele_is_high = true;
  //     }
  //   }
  // }
  // sort_by_pt<Electron>(*event.electrons);
  

  // if(debug)  cout<<"2D cut for low pt"<<endl;
  // //TwoD for low pt:
  // if(isMuon && muon_is_low){
  //   if (debug)cout <<"two d for muon"<<endl;
  //   if(!TwoDCut_selection_low1->passes(event)) return false;
  // }
  // // fill_histograms(event, "TwoDCut_Muon_low1");

  
  // if(isElectron && ele_is_low){
  //   if (debug)cout <<"two d for ele"<<endl;
  //   if(!TwoDCut_selection_low1->passes(event)) return false;
  // }
  
  
  // fill_histograms(event, "TwoDCut_low1");
  
    

  
  // if(debug)  cout<<"done 2D low cut"<<endl;
  // Run top-tagging
//  if(ishotvr){
//    TopTaggerHOTVR->process(event);
//    hadronic_top->process(event);
//  }else if(isdeepAK8){
//    TopTaggerDeepAK8->process(event);
//    hadronic_top->process(event);
//  }
//  if(debug) cout<<"Top Tagger ok"<<endl;

  //fill_histograms(event, "Weights_Init");
 // if(debug)  cout<<"Weights_Init"<<endl;
  // lumihists_Weights_Init->fill(event);
  
  

  if(!HEM_selection->passes(event)){
    if(!isMC) return false;
    else event.weight = event.weight*(1-0.64774715284); // calculated following instructions ar https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis
  }
  // fill_histograms(event, "Weights_HEM");

  // pileup weight
    PUWeight_module->process(event);

  if(debug)  cout<<"PUWeight ok"<<endl;
  // fill_histograms(event, "Weights_PU");
  // lumihists_Weights_PU->fill(event);

  // lumi weight
  LumiWeight_module->process(event);
  if(debug)  cout<<"LumiWeight ok"<<endl;
  // fill_histograms(event, "Weights_Lumi");
  // lumihists_Weights_Lumi->fill(event);

  // top pt reweighting
  TopPtReweight_module->process(event);
  // fill_histograms(event, "Weights_TopPt");
  // lumihists_Weights_TopPt->fill(event);


  // MC scale
  MCScale_module->process(event);
  // fill_histograms(event, "Weights_MCScale");
  // lumihists_Weights_MCScale->fill(event);

  // Prefiring weights
  if (isMC) {
    if (Prefiring_direction == "nominal") event.weight *= event.prefiringWeight;
    else if (Prefiring_direction == "up") event.weight *= event.prefiringWeightUp;
    else if (Prefiring_direction == "down") event.weight *= event.prefiringWeightDown;
  }
  // fill_histograms(event, "Weights_Prefiring");

  // Write PSWeights from genInfo to own branch in output tree
  ps_weights->process(event);
  // fill_histograms(event, "Weights_PS");
  // lumihists_Weights_PS->fill(event);

  // DeepAK8 TopTag SFs
  if(isdeepAK8) sf_toptag->process(event);
  if(debug) cout << "Weights_TopTag_SF: ok" << endl;
  // fill_histograms(event, "Weights_TopTag_SF");
  if(isdeepAK8) sf_topmistag->process(event);
  double muon_pt_high(55.);
  bool muon_is_low = false;
  bool muon_is_high = false;
  if(debug)  cout<<"sorting leptons"<<endl;
  if(isMuon){
    vector<Muon>* muons = event.muons;
    for(unsigned int i=0; i<muons->size(); i++){
      if(event.muons->at(i).pt()<=muon_pt_high){
        muon_is_low = true;
      }else{
        muon_is_high = true;
      }
    }
  }
  sort_by_pt<Muon>(*event.muons);

  double electron_pt_high(120.);
  bool ele_is_low = false;
  bool ele_is_high = false;

  if(isElectron){
    vector<Electron>* electrons = event.electrons;
    for(unsigned int i=0; i<electrons->size(); i++){
      if(event.electrons->at(i).pt()<=electron_pt_high){
        ele_is_low = true;
      }else{
        ele_is_high = true;
      }
    }
    if(debug && event_counter <= 5) cout << "[DEBUG] Finished looping over electrons" << endl;
  }
  sort_by_pt<Electron>(*event.electrons);
  

  if(debug)  cout<<"2D cut for low pt"<<endl;
  //TwoD for low pt:
  if(isMuon && muon_is_low){
    if (debug)cout <<"two d for muon"<<endl;
    if(!TwoDCut_selection_low1->passes(event)) return false;
  }
  // fill_histograms(event, "TwoDCut_Muon_low1");

  
  // if(isElectron && ele_is_low){
  //   if (debug)cout <<"two d for ele"<<endl;
  //   if(!TwoDCut_selection_low1->passes(event)) return false;
  // }
  
  
  // fill_histograms(event, "TwoDCut_low1");
  
    

  
  if(debug)  cout<<"done 2D low cut"<<endl;



  // apply electron id scale factors
  if(isMuon){
    sf_ele_id_dummy->process(event);
  }
  if(isElectron){
    if(ele_is_low){
      sf_ele_id_low->process(event);
    }
    else if(ele_is_high){
      sf_ele_id_high->process(event);
    }
    // fill_histograms(event, "IdEle_SF");
  }

  // apply muon isolation scale factors (low pT only)
  if(isMuon){
    if(muon_is_low){
      if(debug)  cout<<"doing muon iso low"<<endl;
      sf_muon_iso_stat_low->process(event);
      sf_muon_iso_syst_low->process(event);
    }
    else if(muon_is_high){
      if(debug)  cout<<"doing muon iso high"<<endl;
      sf_muon_iso_stat_low_dummy->process(event);
      sf_muon_iso_syst_low_dummy->process(event);

    }
    // fill_histograms(event, "IsoMuon_SF");
  }
  if(isElectron){
     if(debug)  cout<<"doing muon iso dummy"<<endl;
    sf_muon_iso_stat_low_dummy->process(event);
    sf_muon_iso_syst_low_dummy->process(event);
  }
  // apply muon id scale factors
  if(isMuon){
    if(muon_is_low){
       if(debug)  cout<<"doing muon id low"<<endl;
      sf_muon_id_stat_low->process(event);
      sf_muon_id_syst_low->process(event);

    }
    else if(muon_is_high){
       if(debug)  cout<<"doing muon id high"<<endl;
      sf_muon_id_stat_high->process(event);
      sf_muon_id_syst_low->process(event);
    }
    // fill_histograms(event, "IdMuon_SF");
  }
  if(isElectron){
     if(debug)  cout<<"doing muon id dummy"<<endl;
    sf_muon_id_stat_dummy->process(event);
    sf_muon_id_syst_dummy->process(event);

  }

  // apply electron reco scale factors
  if(isMuon){
    sf_ele_reco_dummy->process(event);
  }
  if(isElectron){
    sf_ele_reco->process(event);
    // fill_histograms(event, "RecoEle_SF");
  }

  // apply muon reco scale factors
  sf_muon_reco->process(event);
  // fill_histograms(event, "MuonReco_SF");
   

  // apply lepton trigger scale factors
  if(isMuon){
    if(muon_is_low){
       if(debug)  cout<<"doing muon trigger low"<<endl;
      sf_muon_trigger_stat_low->process(event);
      sf_muon_trigger_syst_low->process(event);

    }
    if(muon_is_high){
      if(debug)  cout<<"doing muon trigger high"<<endl;
      sf_muon_trigger_stat_high->process(event);
      sf_muon_trigger_syst_high->process(event);

    }
    // fill_histograms(event, "TriggerMuon_SF");
  }
  if(isElectron){
    if(debug)  cout<<"doing muon trigger dummy"<<endl;
    sf_muon_trigger_stat_dummy->process(event);
    sf_muon_trigger_syst_dummy->process(event);

  }
  if(debug) cout << "leptons: ok" << endl;
  //Fill histograms before BTagging SF - used to extract Custom BTag SF in (NJets,HT)
  // fill_histograms(event, "BeforeBtagSF");

  // btag shape sf (Ak4 chs jets)
  // new: using new modules, with PUPPI-CHS matching
  sf_btagging->process(event);
  // fill_histograms(event, "AfterBtagSF");

  // apply custom SF to correct for BTag SF shape effects on NJets/HT
  if(isMC && isMuon){
    float custom_sf;

    vector<Jet>* jets = event.jets;
    int Njets = jets->size();
    double st_jets = 0.;
    for(const auto & jet : *jets) st_jets += jet.pt();
    custom_sf = ratio_hist_muon->GetBinContent( ratio_hist_muon->GetXaxis()->FindBin(Njets), ratio_hist_muon->GetYaxis()->FindBin(st_jets) );

    event.weight *= custom_sf;
  }
  if(isMC && !isMuon){
    float custom_sf;

    vector<Jet>* jets = event.jets;
    int Njets = jets->size();
    double st_jets = 0.;
    for(const auto & jet : *jets) st_jets += jet.pt();
    custom_sf = ratio_hist_ele->GetBinContent( ratio_hist_ele->GetXaxis()->FindBin(Njets), ratio_hist_ele->GetYaxis()->FindBin(st_jets) );

    event.weight *= custom_sf;
  }
  // fill_histograms(event, "AfterCustomBtagSF");
  
  // Higher order corrections - EWK & QCD NLO
  NLOCorrections_module->process(event);
  // fill_histograms(event, "NLOCorrections");
  
  //apply ele trigger sf
  sf_ele_trigger->process(event);
  // fill_histograms(event, "TriggerEle_SF");
  // fill_histograms(event, "AfterBaseline");

  CandidateBuilder->process(event);
  if(debug) cout << "CandidateBuilder: ok" << endl;
  Chi2DiscriminatorZprime->process(event);
  if(debug) cout << "Chi2DiscriminatorZprime: ok" << endl;
  // CorrectMatchDiscriminatorZprime->process(event);
  // if(debug) cout << "CorrectMatchDiscriminatorZprime: ok" << endl;
  
  //check SR and CR without DNN
  if(Chi2_selection->passes(event)){
    if(debug) cout << "chi2: ok" << endl;
    // fill_histograms(event, "Chi2_passes");
    // if(ZprimeTopTag_selection->passes(event)){
    //   if(debug) cout << "top tag: ok" << endl;
  
    //   fill_histograms(event, "Chi2_withTopTag");
    // }
    // else{
    //   if(debug) cout << "not top tag" << endl;
  
    //   fill_histograms(event, "Chi2_noTopTag");
    // }
  }
  else{
    if(debug) cout << "noy chi2" << endl;
  
    // fill_histograms(event, "Chi2_inverse");
  }

  // Variables for NN
  if(debug) cout << "before var module" << endl;
  // 
  Variables_module->process(event);
  // fill_histograms(event, "NNInputsBeforeReweight");
  if(debug) cout << "Variables_module: ok" << endl;

  // NN module
  NNModule->process(event);
  std::vector<tensorflow::Tensor> NNoutputs = NNModule->GetOutputs();
  ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
  // bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2); 
  float Mass_tt = BestZprimeCandidate->Zprime_v4().M();

  if(debug) cout << "starting DNN" << endl;
  //float Mttbar_reco =inv_mass(BestZprimeCandidate->BestZprimeCandidate->top_leptonic_v4()+BestZprimeCandidate->BestZprimeCandidate->top_hadronic_v4());
  //cout << "what is Mttbar:" << Mttbar_reco<<endl;
  event.set(h_NNoutput0, (double)(NNoutputs[0].tensor<float, 2>()(0,0)));
  event.set(h_NNoutput1, (double)(NNoutputs[0].tensor<float, 2>()(0,1)));
  event.set(h_NNoutput2, (double)(NNoutputs[0].tensor<float, 2>()(0,2)));
  event.set(h_NNoutput, NNoutputs);

  double out0 = (double)(NNoutputs[0].tensor<float, 2>()(0,0));
  double out1 = (double)(NNoutputs[0].tensor<float, 2>()(0,1));
  double out2 = (double)(NNoutputs[0].tensor<float, 2>()(0,2));
  vector<double> out_event = {out0, out1, out2};

  // h_MulticlassNN_output->fill(event);

  //

  double max_score = 0.0;
  for ( int i = 0; i < 3; i++ ) {
    if ( out_event[i] > max_score) {
      max_score = out_event[i];
    }
  }
  if (debug) cout <<"done setting scores" <<endl;
  // Veto events with >= 2 TopTagged large-R jets
  if(!TopTagVetoSelection->passes(event)) return false;
  // fill_histograms(event, "TopTagVeto");

  if(!DeltaEta_selection->passes(event)) return false;
  // fill_histograms(event, "DeltaEtaCut");

  if(Chi2_selection->passes(event)){ 
    fill_histograms(event, "AfterChi2");
  }
  // VariablesEFTSR_module->process(event);
  // if(debug) cout << "done EFT SR" << endl;
  // VariablesEFTCR1_module->process(event);
  // if(debug) cout << "done EFT CR1" << endl;

  // VariablesEFTCR2_module->process(event);
  // if(debug) cout << "done EFT CR2" << endl;



  // out0=TTbar, out1=ST, out2=WJets
  if( out0 == max_score ){
    if(debug) cout << "inside signal node, about to process EFT vars" << endl;
    VariablesEFTSR_module->process(event);
    fill_histograms(event, "DNN_output0_nochi2");
    if(debug) cout << "signal DNN output0" << endl;
    if(Chi2_selection->passes(event)){  // cut on chi2<30 - only in SR == out0)
      // cout << "signal DNN output0 chi2" << endl;
      // h_CHSMatchHists->fill(event);
      fill_histograms(event, "DNN_output0");
      // cout << "signal DNN output0 chi2 fill syst vars" << endl;
      // h_DeltaY_reco_SystVariations_Inclusive_SR->fill(event);
      // cout << "signal DNN output0 chi2 fill pdf vars" << endl;
      // h_DeltaY_reco_PDFVariations_Inclusive_SR->fill(event);
      // cout << "signal DNN output0 chi2 fill mass bins" << endl;
      // if(Mass_tt>=0 && Mass_tt < 500){
      //   fill_histograms(event, "DeltaY_reco_0_500_SR");
      //   if(debug) cout << "signal DNN output0 chi2 0_500" << endl;
      //   // h_DeltaY_reco_PDFVariations_0_500_SR->fill(event);
      //   if(debug) cout << "signal PDF 0_500" << endl;
      //   //h_DeltaY_reco_SystVariations_0_500_SR->fill(event);
      //   if(debug) cout << "signal all syst vars" << endl;
      // }
      // if(Mass_tt>=0 && Mass_tt < 350){
      //   if(debug) cout << "signal DNN output0 chi2 0_350" << endl;
      //   fill_histograms(event, "DeltaY_reco_0_350_SR");
      //   //h_DeltaY_reco_SystVariations_0_350_SR->fill(event);
      //   // h_DeltaY_reco_PDFVariations_0_350_SR->fill(event);
      // }
      // if(Mass_tt>=350 && Mass_tt < 500){
      //   if(debug) cout << "signal DNN output0 chi2 350_500" << endl;
      //   fill_histograms(event, "DeltaY_reco_350_500_SR");
      //   //h_DeltaY_reco_SystVariations_350_500_SR->fill(event);
      //   // h_DeltaY_reco_PDFVariations_350_500_SR->fill(event);
      // }
      // if (Mass_tt>=0 && Mass_tt < 700){
      //   if(debug) cout << "signal DNN output0 chi2 0_700" << endl;
      //   fill_histograms(event, "DeltaY_reco_0_700_SR");
      //   //h_DeltaY_reco_SystVariations_0_700_SR->fill(event);
      //   // h_DeltaY_reco_PDFVariations_0_700_SR->fill(event);
      // }
      // if(Mass_tt>=700 && Mass_tt < 900){
      //   if(debug) cout << "signal DNN output0 chi2 700_900" << endl;
      //   fill_histograms(event, "DeltaY_reco_700_900_SR");
      //   //h_DeltaY_reco_SystVariations_700_900_SR->fill(event);
      //   // h_DeltaY_reco_PDFVariations_700_900_SR->fill(event);
      // }
      // if(Mass_tt>=900 ){
      //   if(debug) cout << "signal DNN output0 chi2 900_1500" << endl;
      //   fill_histograms(event, "DeltaY_reco_900Inf_SR");
      //   //h_DeltaY_reco_SystVariations_900Inf_SR->fill(event);
      //   // h_DeltaY_reco_PDFVariations_900Inf_SR->fill(event);
      // }
      // if(Mass_tt>=500 && Mass_tt < 750){
      //   if(debug) cout << "signal DNN output0 chi2 500_750" << endl;
      //   fill_histograms(event, "DeltaY_reco_500_750_SR");
      //   //h_DeltaY_reco_SystVariations_500_750_SR->fill(event);
      //   // h_DeltaY_reco_PDFVariations_500_750_SR->fill(event);
      // }
      // if(Mass_tt>=750 && Mass_tt < 1000){
      //   if(debug) cout << "signal DNN output0 chi2 750_1000" << endl;
      //   fill_histograms(event, "DeltaY_reco_750_1000_SR");
      //   //h_DeltaY_reco_SystVariations_750_1000_SR->fill(event);
      //   // h_DeltaY_reco_PDFVariations_750_1000_SR->fill(event);
      // }
      // if(Mass_tt>=1000 && Mass_tt < 1500){
      //   if(debug) cout << "signal DNN output0 chi2 1000_1500" << endl;
      //   fill_histograms(event, "DeltaY_reco_1000_1500_SR");
      //   //h_DeltaY_reco_SystVariations_1000_1500_SR->fill(event);
      //   // h_DeltaY_reco_PDFVariations_1000_1500_SR->fill(event);
      // }
      // if(Mass_tt>=1500){
      //   if(debug) cout << "signal DNN output0 chi2 1500_1000" << endl;
      //   fill_histograms(event, "DeltaY_reco_1500Inf_SR");
      //   //h_DeltaY_reco_SystVariations_1500Inf_SR->fill(event);
      //   // h_DeltaY_reco_PDFVariations_1500Inf_SR->fill(event);
      // }
     
      // if( ZprimeTopTag_selection->passes(event) ){
      //   fill_histograms(event, "DNN_output0_TopTag");
      // }
      // else{
      //   fill_histograms(event, "DNN_output0_NoTopTag");
      // }
      // if(debug) cout << "done with output0" << endl;

    }//Chi2
  }//out0
  // if( out1 == max_score ){
  //   if(debug) cout << "inside ST node, about to process EFT vars" << endl;
  //   VariablesEFTCR1_module->process(event);
  //   fill_histograms(event, "DNN_output1");
  //   if (debug)cout<<"processed CR1"<<endl;

  //   if(Mass_tt>=0 && Mass_tt < 500){
  //     fill_histograms(event, "DeltaY_reco_0_500_CR1");
  //     //h_DeltaY_reco_SystVariations_0_500_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_0_500_CR1->fill(event);  
  //   }
  //   if (Mass_tt>=0 && Mass_tt < 350){
  //     fill_histograms(event, "DeltaY_reco_0_350_CR1");
  //     //h_DeltaY_reco_SystVariations_0_350_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_0_350_CR1->fill(event);
  //   }
  //   if(Mass_tt>=350 && Mass_tt < 500){
  //     fill_histograms(event, "DeltaY_reco_350_500_CR1");
  //     //h_DeltaY_reco_SystVariations_350_500_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_350_500_CR1->fill(event);
  //   }
  //   if (Mass_tt>=0 && Mass_tt < 700){
  //     fill_histograms(event, "DeltaY_reco_0_700_CR1");
  //     //h_DeltaY_reco_SystVariations_0_700_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_0_700_CR1->fill(event);
  //   }
  //   if(Mass_tt>=700 && Mass_tt < 900){
  //     fill_histograms(event, "DeltaY_reco_700_900_CR1");
  //     //h_DeltaY_reco_SystVariations_700_900_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_700_900_CR1->fill(event);
  //   }
  //   if(Mass_tt>=900 ){
  //     fill_histograms(event, "DeltaY_reco_900Inf_CR1");
  //     //h_DeltaY_reco_SystVariations_900Inf_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_900Inf_CR1->fill(event);
  //   }
  //   if(Mass_tt>=500 && Mass_tt < 750){
  //     fill_histograms(event, "DeltaY_reco_500_750_CR1");
  //     //h_DeltaY_reco_SystVariations_500_750_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_500_750_CR1->fill(event);
  //   }
  //   if(Mass_tt>=750 && Mass_tt < 1000){
  //     fill_histograms(event, "DeltaY_reco_750_1000_CR1");
  //     //h_DeltaY_reco_SystVariations_750_1000_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_750_1000_CR1->fill(event);
  //   }
  //   if(Mass_tt>=1000 && Mass_tt < 1500){
  //     fill_histograms(event, "DeltaY_reco_1000_1500_CR1");
  //     //h_DeltaY_reco_SystVariations_1000_1500_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_1000_1500_CR1->fill(event);
  //   }
  //   if(Mass_tt>=1500){
  //     fill_histograms(event, "DeltaY_reco_1500Inf_CR1");
  //     //h_DeltaY_reco_SystVariations_1500Inf_CR1->fill(event);
  //     // h_DeltaY_reco_PDFVariations_1500Inf_CR1->fill(event);
  //   }
  //   if(Chi2_selection->passes(event)){ 
  //     fill_histograms(event,"DNN_output1_chi2");
  //   }
  // }//out1
 
  // if(debug) cout << "check for WJets node" << endl;

  // if( out2 == max_score ){
  //   if (debug)cout<<"inside WJets node about to process EFT vars"<<endl;
  //   VariablesEFTCR2_module->process(event);
  //   fill_histograms(event, "DNN_output2");
  //   if(Mass_tt>=0 && Mass_tt < 500){
  //     fill_histograms(event, "DeltaY_reco_0_500_CR2");
  //     //h_DeltaY_reco_SystVariations_0_500_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_0_500_CR2->fill(event);
  //   }
  //   if (Mass_tt>=0 && Mass_tt < 350){
  //     fill_histograms(event, "DeltaY_reco_0_350_CR2");
  //     //h_DeltaY_reco_SystVariations_0_350_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_0_350_CR2->fill(event);
  //   }
  //   if(Mass_tt>=350 && Mass_tt < 500){
  //     fill_histograms(event, "DeltaY_reco_350_500_CR2");
  //     //h_DeltaY_reco_SystVariations_350_500_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_350_500_CR2->fill(event);
  //   }
  //   if (Mass_tt>=0 && Mass_tt < 700){
  //     fill_histograms(event, "DeltaY_reco_0_700_CR2");
  //     //h_DeltaY_reco_SystVariations_0_700_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_0_700_CR2->fill(event);
  //   }
  //   if(Mass_tt>=700 && Mass_tt < 900){
  //     fill_histograms(event, "DeltaY_reco_700_900_CR2");
  //     //h_DeltaY_reco_SystVariations_700_900_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_700_900_CR2->fill(event);
  //   }
  //   if(Mass_tt>=900 ){
  //     fill_histograms(event, "DeltaY_reco_900Inf_CR2");
  //     //h_DeltaY_reco_SystVariations_900Inf_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_900Inf_CR2->fill(event);
  //   }
  //   if(Mass_tt>=500 && Mass_tt < 750){
  //     fill_histograms(event, "DeltaY_reco_500_750_CR2");
  //     //h_DeltaY_reco_SystVariations_500_750_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_500_750_CR2->fill(event);
  //   }
  //   if(Mass_tt>=750 && Mass_tt < 1000){
  //     fill_histograms(event, "DeltaY_reco_750_1000_CR2");
  //     //h_DeltaY_reco_SystVariations_750_1000_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_750_1000_CR2->fill(event);
  //   }
  //   if(Mass_tt>=1000 && Mass_tt < 1500){
  //     fill_histograms(event, "DeltaY_reco_1000_1500_CR2");
  //     //h_DeltaY_reco_SystVariations_1000_1500_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_1000_1500_CR2->fill(event);
  //   }
  //   if(Mass_tt>=1500){
  //     fill_histograms(event, "DeltaY_reco_1500Inf_CR2");
  //     // h_DeltaY_reco_SystVariations_1500Inf_CR2->fill(event);
  //     // h_DeltaY_reco_PDFVariations_1500Inf_CR2->fill(event);
  //   }
  //   if(Chi2_selection->passes(event)){ 
  //     fill_histograms(event,"DNN_output2_chi2");
  //   }
  // }//out2
  if(debug) cout << "done with DNNs" << endl;
  // if(debug) cout << "done" << endl;

  // Calculate structure constants for EFT weights
  // This accesses EFT weights starting at index 202 in event.genInfo->systweights()
  // and calculates structure constants that can be used to compute weights for any WC values

  // calculates the structure constants for each event.
  if(debug) cout << "before structure constants" << endl;
  if(debug) cout << "isEFT: " << isEFT << endl;
  if(isEFT){
    if(debug) cout<<" should not be in here if not EFT" << endl;
    structure_constants_calculator->process(event);
  }
  // Shows the number of structure constants stored in the event
  // Displays the first few structure constants
  // Shows the constant term (SM point) and a few linear terms
  // Prints for the first 5 EFT events
  
  // Debug output for structure constants (only for first few events)
  // if (debug && isEFT) {static int event_counter = 0;
  //   cout << "about to check structure constants debug " << endl;
  //   if (event_counter < 5) {
  //     // Get the structure constants from the event
  //     if (event.is_valid(h_structure_constants)) {
  //       std::vector<float> structure_constants = event.get(h_structure_constants);
        
  //       std::cout << "===== Structure Constants Debug (Event " << event_counter << ") =====" << std::endl;
  //       std::cout << "Number of structure constants: " << structure_constants.size() << std::endl;
        
  //       if (!structure_constants.empty()) {
  //         // Print first few constants
  //         std::cout << "First few constants: ";
  //         for (size_t i = 0; i < std::min(size_t(10), structure_constants.size()); ++i) {
  //           std::cout << structure_constants[i] << " ";
  //         }
  //         std::cout << std::endl;
  //       }
        
  //       // Increment counter after printing
  //       event_counter++;
  //     } else {
  //       std::cout << "Structure constants not found in event!" << std::endl;
  //     }
  //   }
  // }
  if(debug) cout << "moving on to next event" << endl;
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule_applyNN)