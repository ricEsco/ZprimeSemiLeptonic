#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include <UHH2/common/include/PrintingModules.h>

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
#include <UHH2/common/include/LuminosityHists.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/TopPtReweight.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/LeptonScaleFactors.h>
#include <UHH2/common/include/PSWeights.h>

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/TTbarLJHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicCHSMatchHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>
#include <UHH2/ZprimeSemiLeptonic/include/ElecTriggerSF.h>
#include <UHH2/ZprimeSemiLeptonic/include/AK4JetCorrections.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopPuppiJetCorrections.h>
//#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSystematicsModule.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopTagScaleFactor.h>

#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <UHH2/HOTVR/include/HadronicTop.h>
#include <UHH2/HOTVR/include/HOTVRScaleFactor.h>
#include <UHH2/HOTVR/include/HOTVRIds.h>

using namespace std;
using namespace uhh2;

/*
██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class ZprimeAnalysisModule_AzCorr_GenStudy : public ModuleBASE {

public:

  explicit ZprimeAnalysisModule_AzCorr_GenStudy(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:

  bool debug;

  // Cleaners
  std::unique_ptr<MuonCleaner> muon_cleaner_low, muon_cleaner_high;
  std::unique_ptr<ElectronCleaner> electron_cleaner_low, electron_cleaner_high;

  // scale factors
  unique_ptr<AnalysisModule> sf_muon_iso_low, sf_muon_id_low, sf_muon_id_high, sf_muon_trigger_low, sf_muon_trigger_high;
  unique_ptr<AnalysisModule> sf_muon_iso_low_dummy, sf_muon_id_dummy, sf_muon_trigger_dummy;
  unique_ptr<AnalysisModule> sf_ele_id_low, sf_ele_id_high, sf_ele_reco;
  unique_ptr<AnalysisModule> sf_ele_id_dummy, sf_ele_reco_dummy;
  unique_ptr<MuonRecoSF> sf_muon_reco;
  unique_ptr<AnalysisModule> sf_ele_trigger;
  unique_ptr<AnalysisModule> sf_btagging;

  // AnalysisModules
  unique_ptr<AnalysisModule> LumiWeight_module, PUWeight_module, TopPtReweight_module, MCScale_module;
  unique_ptr<AnalysisModule> NLOCorrections_module;
  unique_ptr<PSWeights> ps_weights;

  // Top tagging
  unique_ptr<HOTVRTopTagger> TopTaggerHOTVR;
  unique_ptr<AnalysisModule> hadronic_top;
  unique_ptr<AnalysisModule> sf_toptag;
  unique_ptr<DeepAK8TopTagger> TopTaggerDeepAK8;

  // Mass reconstruction
  unique_ptr<ZprimeCandidateBuilder> CandidateBuilder;

  // Chi2 discriminator
  unique_ptr<ZprimeChi2Discriminator> Chi2DiscriminatorZprime;
  unique_ptr<ZprimeCorrectMatchDiscriminator> CorrectMatchDiscriminatorZprime;

  // Selections
  unique_ptr<Selection> MuonVeto_selection, EleVeto_selection, NMuon1_selection, NEle1_selection;
  unique_ptr<Selection> Trigger_mu_A_selection, Trigger_mu_B_selection, Trigger_mu_C_selection, Trigger_mu_D_selection, Trigger_mu_E_selection, Trigger_mu_F_selection;
  unique_ptr<Selection> Trigger_ele_A_selection, Trigger_ele_B_selection, Trigger_ph_A_selection;
  unique_ptr<Selection> TwoDCut_selection, Jet1_selection, Jet2_selection, Met_selection, Chi2_selection, TTbarMatchable_selection, Chi2CandidateMatched_selection, ZprimeTopTag_selection, BlindData_selection;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  std::unique_ptr<Selection> sel_1btag, sel_2btag;
  std::unique_ptr<Selection> HEM_selection;
  unique_ptr<Selection> SignSplit;

  // NN variables handles
  unique_ptr<Variables_NN> Variables_module;

  // systematics handles
  //unique_ptr<ZprimeSemiLeptonicSystematicsModule> SystematicsModule;

  //Handles
  Event::Handle<bool> h_is_zprime_reconstructed_chi2, h_is_zprime_reconstructed_correctmatch;
  Event::Handle<float> h_chi2;   Event::Handle<float> h_weight;
  Event::Handle<float> h_MET;   Event::Handle<int> h_NPV;
  Event::Handle<float> h_lep1_pt; Event::Handle<float> h_lep1_eta;
  Event::Handle<float> h_ak4jet1_pt; Event::Handle<float> h_ak4jet1_eta;
  Event::Handle<float> h_ak8jet1_pt; Event::Handle<float> h_ak8jet1_eta;
  Event::Handle<float> h_Mttbar;

  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;

  // SpinCorr variables
  uhh2::Event::Handle<std::vector<Jet>> h_CHSjets_matched;  // Collection of CHS matched jets
  uhh2::Event::Handle<std::vector<TopJet>> h_DeepAK8TopTags;  // Collection of DeepAK8TopTagged jets
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<TTbarGenProducer> ttgenprod;
  
  Event::Handle<double> h_pt_hadTop;     // pt of hadronic top-jet(s)

  // Plotting mass of ttbar system
  Event::Handle<double> h_ttbar_mass_LabFrame;

  // Plotting longitudinal boost of ttbar system
  Event::Handle<double> h_ttbar_boost_LabFrame;

  // Phi of lepton from leptonic leg
  Event::Handle<double> h_phi_lep_LabFrame;
  Event::Handle<double> h_phi_lep_CoMFrame;
  Event::Handle<double> h_phi_lep_helicityFrame;
  Event::Handle<double> h_phi_lep;
  Event::Handle<double> h_phi_lep_high; 
  Event::Handle<double> h_phi_lep_low;
  Event::Handle<double> h_phi_lep_highMass; 
  Event::Handle<double> h_phi_lep_lowMass;

  // // Phi of b-jet from hadronic leg
  // Event::Handle<double> h_phi_b_LabFrame;
  // Event::Handle<double> h_phi_b_CoMFrame;
  // Event::Handle<double> h_phi_b_helicityFrame;
  // Event::Handle<double> h_phi_b;
  // Event::Handle<double> h_phi_b_high; 
  // Event::Handle<double> h_phi_b_low;
  // Event::Handle<double> h_phi_b_highMass; 
  // Event::Handle<double> h_phi_b_lowMass;

  // Phi of q1 W daughter from hadronic leg
  Event::Handle<double> h_phi_hadTop_q1_LabFrame;
  Event::Handle<double> h_phi_hadTop_q1_CoMFrame;
  Event::Handle<double> h_phi_hadTop_q1_helicityFrame;

  // Phi of q2 W daughter from hadronic leg
  Event::Handle<double> h_phi_hadTop_q2_LabFrame;
  Event::Handle<double> h_phi_hadTop_q2_CoMFrame;
  Event::Handle<double> h_phi_hadTop_q2_helicityFrame;

  // Phi of less-energetic W daughter from hadronic leg
  Event::Handle<double> h_phi_qlow;
  Event::Handle<double> h_phi_qlow_high; 
  Event::Handle<double> h_phi_qlow_low;
  Event::Handle<double> h_phi_qlow_highMass; 
  Event::Handle<double> h_phi_qlow_lowMass;

  // Sum of phi-coordinates
  Event::Handle<double> h_sphi;
  Event::Handle<double> h_sphi_low;
  Event::Handle<double> h_sphi_high;
  Event::Handle<double> h_sphi_Mass1;
  Event::Handle<double> h_sphi_Mass2;
  Event::Handle<double> h_sphi_Mass3;
  Event::Handle<double> h_sphi_Mass4;
  Event::Handle<double> h_sphi_boost1;
  Event::Handle<double> h_sphi_boost2;
  Event::Handle<double> h_sphi_boost3;
  Event::Handle<double> h_sphi_boost4;

  // Difference of phi-coordinates
  Event::Handle<double> h_dphi;
  Event::Handle<double> h_dphi_low;
  Event::Handle<double> h_dphi_high;
  Event::Handle<double> h_dphi_Mass1;
  Event::Handle<double> h_dphi_Mass2;
  Event::Handle<double> h_dphi_Mass3;
  Event::Handle<double> h_dphi_Mass4;
  Event::Handle<double> h_dphi_boost1;
  Event::Handle<double> h_dphi_boost2;
  Event::Handle<double> h_dphi_boost3;
  Event::Handle<double> h_dphi_boost4;

  // Above variables now separated by charge of lepton in system
  Event::Handle<double> h_phi_lepPlus;
  Event::Handle<double> h_phi_lepMinus;
  // sphi_plus  is defined as (lepTop + hadTop) for positive leptons
  Event::Handle<double> h_sphi_plus;
  Event::Handle<double> h_sphi_plus_low;
  Event::Handle<double> h_sphi_plus_high;
  // dphi_plus is defined as (lepTop - hadTop) for positive leptons
  Event::Handle<double> h_dphi_plus;
  Event::Handle<double> h_dphi_plus_low;
  Event::Handle<double> h_dphi_plus_high;
  // sphi_minus is defined as (hadTop + lepTop) for negative leptons
  Event::Handle<double> h_sphi_minus;
  Event::Handle<double> h_sphi_minus_low;
  Event::Handle<double> h_sphi_minus_high;
  // dphi_minus is defined as (hadTop - lepTop) for negative leptons
  Event::Handle<double> h_dphi_minus;
  Event::Handle<double> h_dphi_minus_low;
  Event::Handle<double> h_dphi_minus_high;

  // Lumi hists
  std::unique_ptr<Hists> lumihists_Weights_Init, lumihists_Weights_PU, lumihists_Weights_Lumi, lumihists_Weights_TopPt, lumihists_Weights_MCScale, lumihists_Weights_PS, lumihists_Muon1_LowPt, lumihists_Muon1_HighPt, lumihists_Ele1_LowPt, lumihists_Ele1_HighPt, lumihists_TriggerMuon, lumihists_TriggerEle, lumihists_TwoDCut_Muon, lumihists_TwoDCut_Ele, lumihists_Jet1, lumihists_Jet2, lumihists_MET, lumihists_HTlep, lumihists_Chi2;

  // PUPPI CHS match module
  std::unique_ptr<PuppiCHS_matching> AK4PuppiCHS_matching;
  // PUPPI CHS match - btagging
  std::unique_ptr<Selection> AK4PuppiCHS_BTagging;
  // Hists with matched CHS jets
  std::unique_ptr<Hists> h_CHSMatchHists;
  std::unique_ptr<Hists> h_CHSMatchHists_beforeBTagSF;
  std::unique_ptr<Hists> h_CHSMatchHists_afterBTagSF;
  std::unique_ptr<Hists> h_CHSMatchHists_after2DBTagSF;
  std::unique_ptr<Hists> h_CHSMatchHists_afterBTag;

  // Configuration
  bool isMC, ishotvr, isdeepAK8;
  string Sys_PU, Prefiring_direction, Sys_TopPt_a, Sys_TopPt_b;
  TString sample;
  int runnr_oldtriggers = 299368;

  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;
  bool isMuon, isElectron;
  bool isPhoton;
  bool isEleTriggerMeasurement;
  TString year, channel;

  TH2F *ratio_hist_muon;
  TH2F *ratio_hist_ele;

};

void ZprimeAnalysisModule_AzCorr_GenStudy::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag + "_Skimming";
    mytag = tag + "_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicHists(ctx,mytag));
  }
}

void ZprimeAnalysisModule_AzCorr_GenStudy::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag + "_Skimming";
  mytag = tag + "_General";
  HFolder(mytag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

ZprimeAnalysisModule_AzCorr_GenStudy::ZprimeAnalysisModule_AzCorr_GenStudy(uhh2::Context& ctx){

  debug = false; // false/true

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

  isEleTriggerMeasurement = (ctx.get("is_EleTriggerMeasurement") == "true");

  // Lepton IDs
  ElectronId eleID_low  = ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp80);
  ElectronId eleID_high = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp80);
  MuonId     muID_low   = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  MuonId     muID_high  = MuonID(Muon::CutBasedIdGlobalHighPt);

  double electron_pt_low;
  if(isUL17){
    electron_pt_low = 38.; // UL17 ele trigger threshold is 35 (HLT_Ele35WPTight _Gsf) -> be above turn on
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

  MuonVeto_selection.reset(new NMuonSelection(0, 0));
  EleVeto_selection.reset(new NElectronSelection(0, 0));
  NMuon1_selection.reset(new NMuonSelection(1, 1));
  NEle1_selection.reset(new NElectronSelection(1, 1));

  // Important selection values
  double jet1_pt(50.);
  double jet2_pt;
  double chi2_max(30.);
  // double mtt_blind(3000.);
  string trigger_mu_A, trigger_mu_B, trigger_mu_C, trigger_mu_D, trigger_mu_E, trigger_mu_F;
  string trigger_ele_A, trigger_ele_B;
  string trigger_ph_A;
  double MET_cut, HT_lep_cut;
  isMuon = false; isElectron = false;
  channel = ctx.get("channel");
  if(channel == "muon") isMuon = true;
  if(channel == "electron") isElectron = true;


  if(isMuon){ // semileptonic muon channel
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

    jet2_pt = 50;
    MET_cut = 70;
    HT_lep_cut = 0;
  }
  if(isElectron){ // semileptonic electron channel
    trigger_ele_B = "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*";
    if(isUL16preVFP || isUL16postVFP){
      trigger_ele_A = "HLT_Ele27_WPTight_Gsf_v*";
    }
    if(isUL17){
      trigger_ele_A = "HLT_Ele35_WPTight_Gsf_v*";
    }
    if(isUL18){
      trigger_ele_A = "HLT_Ele32_WPTight_Gsf_v*";
    }
    if(isUL16preVFP || isUL16postVFP){
      trigger_ph_A = "HLT_Photon175_v*";
    }
    else{
      trigger_ph_A = "HLT_Photon200_v*";
    }
    jet2_pt = 40;
    MET_cut = 60;
    HT_lep_cut = 0;
  }


  double TwoD_dr = 0.4, TwoD_ptrel = 25.;

  // const TopJetId toptagID = AndId<TopJet>(HOTVRTopTag(0.8, 140.0, 220.0, 50.0), Tau32Groomed(0.56));

  Sys_PU = ctx.get("Sys_PU");
  Prefiring_direction = ctx.get("Sys_prefiring");
  Sys_TopPt_a = ctx.get("Systematic_TopPt_a");
  Sys_TopPt_b = ctx.get("Systematic_TopPt_b");

  BTag::algo btag_algo = BTag::DEEPJET;
  BTag::wp btag_wp = BTag::WP_MEDIUM;
  JetId id_btag = BTag(btag_algo, btag_wp);

  // double a_toppt = 0.0615; // par a TopPt Reweighting
  // double b_toppt = -0.0005; // par b TopPt Reweighting

  // Modules
  LumiWeight_module.reset(new MCLumiWeight(ctx));
  PUWeight_module.reset(new MCPileupReweight(ctx, Sys_PU));
  //TopPtReweight_module.reset(new TopPtReweighting(ctx, a_toppt, b_toppt, Sys_TopPt_a, Sys_TopPt_b, ""));
  MCScale_module.reset(new MCScaleVariation(ctx));
  hadronic_top.reset(new HadronicTop(ctx));

  // sf_toptag.reset(new HOTVRScaleFactor(ctx, toptagID, ctx.get("Sys_TopTag", "nominal"), "HadronicTop", "TopTagSF", "HOTVRTopTagSFs"));
  sf_toptag.reset(new TopTagScaleFactor(ctx));
  NLOCorrections_module.reset(new NLOCorrections(ctx));
  ps_weights.reset(new PSWeights(ctx));

  // b-tagging SFs
  sf_btagging.reset(new MCBTagDiscriminantReweighting(ctx, BTag::algo::DEEPJET, "CHS_matched"));

  // set lepton scale factors: see UHH2/common/include/LeptonScaleFactors.h
  sf_muon_iso_low.reset(new uhh2::MuonIsoScaleFactors(ctx, Muon::Selector::PFIsoTight, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_low.reset(new uhh2::MuonIdScaleFactors(ctx, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_high.reset(new uhh2::MuonIdScaleFactors(ctx, Muon::Selector::CutBasedIdGlobalHighPt, true));
  sf_muon_trigger_low.reset(new uhh2::MuonTriggerScaleFactors(ctx, false, true));
  sf_muon_trigger_high.reset(new uhh2::MuonTriggerScaleFactors(ctx, true, false));
  sf_muon_reco.reset(new MuonRecoSF(ctx));
  sf_ele_id_low.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_iso_V2_wp80, true));
  sf_ele_id_high.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_noIso_V2_wp80, true));
  sf_ele_reco.reset(new uhh2::ElectronRecoScaleFactors(ctx, false, true));

  if(!isEleTriggerMeasurement) sf_ele_trigger.reset( new uhh2::ElecTriggerSF(ctx, "central", "eta_ptbins", year) );

  // dummies (needed to aviod set value errors)
  sf_muon_iso_low_dummy.reset(new uhh2::MuonIsoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_id_dummy.reset(new uhh2::MuonIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_trigger_dummy.reset(new uhh2::MuonTriggerScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_ele_id_dummy.reset(new uhh2::ElectronIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_ele_reco_dummy.reset(new uhh2::ElectronRecoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));

  // Selection modules
  Trigger_mu_A_selection.reset(new TriggerSelection(trigger_mu_A));
  Trigger_mu_B_selection.reset(new TriggerSelection(trigger_mu_B));
  Trigger_mu_C_selection.reset(new TriggerSelection(trigger_mu_C));
  Trigger_mu_D_selection.reset(new TriggerSelection(trigger_mu_D));
  Trigger_mu_E_selection.reset(new TriggerSelection(trigger_mu_E));
  Trigger_mu_F_selection.reset(new TriggerSelection(trigger_mu_F));

  Trigger_ele_A_selection.reset(new TriggerSelection(trigger_ele_A));
  Trigger_ele_B_selection.reset(new TriggerSelection(trigger_ele_B));
  Trigger_ph_A_selection.reset(new TriggerSelection(trigger_ph_A));


  TwoDCut_selection.reset(new TwoDCut(TwoD_dr, TwoD_ptrel));
  Jet1_selection.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.5))));
  Jet2_selection.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.5))));
  met_sel.reset(new METCut  (MET_cut   , uhh2::infinity));
  htlep_sel.reset(new HTlepCut(HT_lep_cut, uhh2::infinity));

  Chi2_selection.reset(new Chi2Cut(ctx, 0., chi2_max));
  TTbarMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
  Chi2CandidateMatched_selection.reset(new Chi2CandidateMatchedSelection(ctx));
  ZprimeTopTag_selection.reset(new ZprimeTopTagSelection(ctx));
  // BlindData_selection.reset(new BlindDataSelection(ctx, mtt_blind));

  HEM_selection.reset(new HEMSelection(ctx)); // HEM issue in 2018, veto on leptons and jets

  Variables_module.reset(new Variables_NN(ctx, mode)); // variables for NN

  //if(!isEleTriggerMeasurement) SystematicsModule.reset(new ZprimeSemiLeptonicSystematicsModule(ctx));


  // Split interference signal samples by sign
  if(ctx.get("dataset_version").find("_int") != std::string::npos){
    if     (ctx.get("dataset_version").find("_pos") != std::string::npos) SignSplit.reset(new SignSelection("pos"));
    else if(ctx.get("dataset_version").find("_neg") != std::string::npos) SignSplit.reset(new SignSelection("neg"));
    else SignSplit.reset(new uhh2::AndSelection(ctx));
  }
  else SignSplit.reset(new uhh2::AndSelection(ctx));

  // Top Taggers
  TopTaggerHOTVR.reset(new HOTVRTopTagger(ctx));
  TopTaggerDeepAK8.reset(new DeepAK8TopTagger(ctx));

  // Zprime candidate builder
  CandidateBuilder.reset(new ZprimeCandidateBuilder(ctx, mode));

  // Zprime discriminators
  if(isMC) ttgenprod.reset(new TTbarGenProducer(ctx));
  Chi2DiscriminatorZprime.reset(new ZprimeChi2Discriminator(ctx));
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  CorrectMatchDiscriminatorZprime.reset(new ZprimeCorrectMatchDiscriminator(ctx));
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_chi2 = ctx.declare_event_output<float> ("rec_chi2");
  h_MET = ctx.declare_event_output<float> ("met_pt");
  h_Mttbar = ctx.declare_event_output<float> ("Mttbar");
  h_lep1_pt = ctx.declare_event_output<float> ("lep1_pt");
  h_lep1_eta = ctx.declare_event_output<float> ("lep1_eta");
  h_ak4jet1_pt = ctx.declare_event_output<float> ("ak4jet1_pt");
  h_ak4jet1_eta = ctx.declare_event_output<float> ("ak4jet1_eta");
  h_ak8jet1_pt = ctx.declare_event_output<float> ("ak8jet1_pt");
  h_ak8jet1_eta = ctx.declare_event_output<float> ("ak8jet1_eta");

  h_NPV = ctx.declare_event_output<int> ("NPV");
  h_weight = ctx.declare_event_output<float> ("weight");

  // Spin Corr Variables
  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");       // Collection of CHS matched jets
  h_DeepAK8TopTags = ctx.get_handle<std::vector<TopJet>>("DeepAK8TopTags");  // Collection of DeepAK8TopTagged jets
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");                         // Access to gen-level particles
  h_pt_hadTop=ctx.declare_event_output<double> ("pt_hadTop");                 // pt of hadronic top-jet(s)
  
  // Plotting mass of ttbar system
  h_ttbar_mass_LabFrame=ctx.declare_event_output<double> ("ttbar_mass_LabFrame");

  // Plotting longitudinal boost of ttbar system
  h_ttbar_boost_LabFrame=ctx.declare_event_output<double> ("ttbar_boost_LabFrame");

  // Phi of lepton from leptonic leg
  h_phi_lep_LabFrame=ctx.declare_event_output<double> ("phi_lep_LabFrame");
  h_phi_lep_CoMFrame=ctx.declare_event_output<double> ("phi_lep_CoMFrame");
  h_phi_lep_helicityFrame=ctx.declare_event_output<double> ("phi_lep_helicityFrame");
  h_phi_lep=ctx.declare_event_output<double> ("phi_lep");
  h_phi_lep_high=ctx.declare_event_output<double> ("phi_lep_high");
  h_phi_lep_low=ctx.declare_event_output<double> ("phi_lep_low");
  h_phi_lep_highMass=ctx.declare_event_output<double> ("phi_lep_highMass");
  h_phi_lep_lowMass=ctx.declare_event_output<double> ("phi_lep_lowMass");

  // // Phi of b-jet from hadronic leg
  // h_phi_b_LabFrame=ctx.declare_event_output<double> ("phi_b_LabFrame");
  // h_phi_b_CoMFrame=ctx.declare_event_output<double> ("phi_b_CoMFrame");
  // h_phi_b_helicityFrame=ctx.declare_event_output<double> ("phi_b_helicityFrame");
  // h_phi_b=ctx.declare_event_output<double> ("phi_b");
  // h_phi_b_high=ctx.declare_event_output<double> ("phi_b_high");
  // h_phi_b_low=ctx.declare_event_output<double> ("phi_b_low");
  // h_phi_b_highMass=ctx.declare_event_output<double> ("phi_b_highMass");
  // h_phi_b_lowMass=ctx.declare_event_output<double> ("phi_b_lowMass");

  // Phi of q1 W daughter from hadronic leg
  h_phi_hadTop_q1_LabFrame=ctx.declare_event_output<double> ("phi_hadTop_q1_LabFrame");
  h_phi_hadTop_q1_CoMFrame=ctx.declare_event_output<double> ("phi_hadTop_q1_CoMFrame");
  h_phi_hadTop_q1_helicityFrame=ctx.declare_event_output<double> ("phi_hadTop_q1_helicityFrame");

  // Phi of q2 W daughter from hadronic leg
  h_phi_hadTop_q2_LabFrame=ctx.declare_event_output<double> ("phi_hadTop_q2_LabFrame");
  h_phi_hadTop_q2_CoMFrame=ctx.declare_event_output<double> ("phi_hadTop_q2_CoMFrame");
  h_phi_hadTop_q2_helicityFrame=ctx.declare_event_output<double> ("phi_hadTop_q2_helicityFrame");

  // Phi of less-energetic (in top's rest frame) W daughter from hadronic leg
  h_phi_qlow=ctx.declare_event_output<double> ("phi_qlow");
  h_phi_qlow_high=ctx.declare_event_output<double> ("phi_qlow_high");
  h_phi_qlow_low=ctx.declare_event_output<double> ("phi_qlow_low");
  h_phi_qlow_highMass=ctx.declare_event_output<double> ("phi_qlow_highMass");
  h_phi_qlow_lowMass=ctx.declare_event_output<double> ("phi_qlow_lowMass");

  // Sum of phi-coordinates
  h_sphi=ctx.declare_event_output<double> ("sphi");
  h_sphi_high=ctx.declare_event_output<double> ("sphi_high");
  h_sphi_low=ctx.declare_event_output<double> ("sphi_low");
  h_sphi_Mass1=ctx.declare_event_output<double> ("sphi_Mass1");
  h_sphi_Mass2=ctx.declare_event_output<double> ("sphi_Mass2");
  h_sphi_Mass3=ctx.declare_event_output<double> ("sphi_Mass3");
  h_sphi_Mass4=ctx.declare_event_output<double> ("sphi_Mass4");
  h_sphi_boost1=ctx.declare_event_output<double> ("sphi_boost1");
  h_sphi_boost2=ctx.declare_event_output<double> ("sphi_boost2");
  h_sphi_boost3=ctx.declare_event_output<double> ("sphi_boost3");
  h_sphi_boost4=ctx.declare_event_output<double> ("sphi_boost4");

  // Difference of phi-coordinates
  h_dphi=ctx.declare_event_output<double> ("dphi");
  h_dphi_high=ctx.declare_event_output<double> ("dphi_high");
  h_dphi_low=ctx.declare_event_output<double> ("dphi_low");
  h_dphi_Mass1=ctx.declare_event_output<double> ("dphi_Mass1");
  h_dphi_Mass2=ctx.declare_event_output<double> ("dphi_Mass2");
  h_dphi_Mass3=ctx.declare_event_output<double> ("dphi_Mass3");
  h_dphi_Mass4=ctx.declare_event_output<double> ("dphi_Mass4");
  h_dphi_boost1=ctx.declare_event_output<double> ("dphi_boost1");
  h_dphi_boost2=ctx.declare_event_output<double> ("dphi_boost2");
  h_dphi_boost3=ctx.declare_event_output<double> ("dphi_boost3");
  h_dphi_boost4=ctx.declare_event_output<double> ("dphi_boost4");

  // Above variables now separated by charge of lepton in system
  h_phi_lepPlus=ctx.declare_event_output<double> ("phi_lepPlus");
  h_phi_lepMinus=ctx.declare_event_output<double> ("phi_lepMinus");
  // sphi_plus  is defined as (lepTop + hadTop) for positive leptons
  h_sphi_plus=ctx.declare_event_output<double> ("sphi_plus");
  h_sphi_plus_low=ctx.declare_event_output<double> ("sphi_plus_low");
  h_sphi_plus_high=ctx.declare_event_output<double> ("sphi_plus_high");
  // dphi_plus is defined as (lepTop - hadTop) for positive leptons
  h_dphi_plus=ctx.declare_event_output<double> ("dphi_plus");
  h_dphi_plus_low=ctx.declare_event_output<double> ("dphi_plus_low");
  h_dphi_plus_high=ctx.declare_event_output<double> ("dphi_plus_high");
  // sphi_minus is defined as (hadTop + lepTop) for negative leptons
  h_sphi_minus=ctx.declare_event_output<double> ("sphi_minus");
  h_sphi_minus_low=ctx.declare_event_output<double> ("sphi_minus_low");
  h_sphi_minus_high=ctx.declare_event_output<double> ("sphi_minus_high");
  // dphi_minus is defined as (hadTop - lepTop) for negative leptons
  h_dphi_minus=ctx.declare_event_output<double> ("dphi_minus");
  h_dphi_minus_low=ctx.declare_event_output<double> ("dphi_minus_low");
  h_dphi_minus_high=ctx.declare_event_output<double> ("dphi_minus_high");


  // Btagging modules (ask for at least 1 or 2 btagged jets in event)
  sel_1btag.reset(new NJetSelection(1, -1, id_btag));
  sel_2btag.reset(new NJetSelection(2,-1, id_btag));

  // PUPPI CHS match modules & hists
  AK4PuppiCHS_matching.reset(new PuppiCHS_matching(ctx)); // match AK4 PUPPI jets to AK$ CHS jets for b-tagging
  AK4PuppiCHS_BTagging.reset(new PuppiCHS_BTagging(ctx)); // b-tagging on matched CHS jets
  h_CHSMatchHists.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch"));
  h_CHSMatchHists_beforeBTagSF.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_beforeBTagSF"));
  h_CHSMatchHists_afterBTagSF.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_afterBTagSF"));
  h_CHSMatchHists_after2DBTagSF.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_after2DBTagSF"));
  h_CHSMatchHists_afterBTag.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_afterBTag"));

  // Book histograms
  vector<string> histogram_tags = {"Weights_Init", "Weights_HEM", "Weights_PU", "Weights_Lumi", "Weights_TopPt", "Weights_MCScale", "Weights_Prefiring", "Weights_PS", "Weights_TopTag_SF", "Muon1_LowPt", "Muon1_HighPt", "Muon1_Tot", "Ele1_LowPt", "Ele1_HighPt", "Ele1_Tot", "1Mu1Ele_LowPt", "1Mu1Ele_HighPt", "1Mu1Ele_Tot", "IdMuon_SF", "IdEle_SF", "IsoMuon_SF", "RecoEle_SF", "MuonReco_SF", "TriggerMuon", "TriggerEle", "TriggerMuon_SF", "TwoDCut_Muon", "TwoDCut_Ele", "Jet1", "Jet2", "MET", "HTlep", "BeforeBtagSF", "AfterBtagSF", "AfterCustomBtagSF", "Btags1", "NLOCorrections", "TriggerEle_SF", "TTbarCandidate", "CorrectMatchDiscriminator", "Chi2Discriminator", "NNInputsBeforeReweight"};
  book_histograms(ctx, histogram_tags);

  lumihists_Weights_Init.reset(new LuminosityHists(ctx, "Lumi_Weights_Init"));
  lumihists_Weights_PU.reset(new LuminosityHists(ctx, "Lumi_Weights_PU"));
  lumihists_Weights_Lumi.reset(new LuminosityHists(ctx, "Lumi_Weights_Lumi"));
  lumihists_Weights_TopPt.reset(new LuminosityHists(ctx, "Lumi_Weights_TopPt"));
  lumihists_Weights_MCScale.reset(new LuminosityHists(ctx, "Lumi_Weights_MCScale"));
  lumihists_Weights_PS.reset(new LuminosityHists(ctx, "Lumi_Weights_PS"));
  lumihists_Muon1_LowPt.reset(new LuminosityHists(ctx, "Lumi_Muon1_LowPt"));
  lumihists_Muon1_HighPt.reset(new LuminosityHists(ctx, "Lumi_Muon1_HighPt"));
  lumihists_Ele1_LowPt.reset(new LuminosityHists(ctx, "Lumi_Ele1_LowPt"));
  lumihists_Ele1_HighPt.reset(new LuminosityHists(ctx, "Lumi_Ele1_HighPt"));
  lumihists_TriggerMuon.reset(new LuminosityHists(ctx, "Lumi_TriggerMuon"));
  lumihists_TriggerEle.reset(new LuminosityHists(ctx, "Lumi_TriggerEle"));
  lumihists_TwoDCut_Muon.reset(new LuminosityHists(ctx, "Lumi_TwoDCut_Muon"));
  lumihists_TwoDCut_Ele.reset(new LuminosityHists(ctx, "Lumi_TwoDCut_Ele"));
  lumihists_Jet1.reset(new LuminosityHists(ctx, "Lumi_Jet1"));
  lumihists_Jet2.reset(new LuminosityHists(ctx, "Lumi_Jet2"));
  lumihists_MET.reset(new LuminosityHists(ctx, "Lumi_MET"));
  lumihists_HTlep.reset(new LuminosityHists(ctx, "Lumi_HTlep"));
  lumihists_Chi2.reset(new LuminosityHists(ctx, "Lumi_Chi2"));

  if(isMC){
    TString sample_name = "";
    vector<TString> names = {"ST", "WJets", "DY", "QCD", "ALP_ttbar_signal", "ALP_ttbar_interference", "HscalarToTTTo1L1Nu2J_m365_w36p5_res", "HscalarToTTTo1L1Nu2J_m400_w40p0_res", "HscalarToTTTo1L1Nu2J_m500_w50p0_res", "HscalarToTTTo1L1Nu2J_m600_w60p0_res", "HscalarToTTTo1L1Nu2J_m800_w80p0_res", "HscalarToTTTo1L1Nu2J_m1000_w100p0_res", "HscalarToTTTo1L1Nu2J_m365_w36p5_int_pos", "HscalarToTTTo1L1Nu2J_m400_w40p0_int_pos", "HscalarToTTTo1L1Nu2J_m500_w50p0_int_pos", "HscalarToTTTo1L1Nu2J_m600_w60p0_int_pos", "HscalarToTTTo1L1Nu2J_m800_w80p0_int_pos", "HscalarToTTTo1L1Nu2J_m1000_w100p0_int_pos", "HscalarToTTTo1L1Nu2J_m365_w36p5_int_neg", "HscalarToTTTo1L1Nu2J_m400_w40p0_int_neg", "HscalarToTTTo1L1Nu2J_m500_w50p0_int_neg", "HscalarToTTTo1L1Nu2J_m600_w60p0_int_neg", "HscalarToTTTo1L1Nu2J_m800_w80p0_int_neg", "HscalarToTTTo1L1Nu2J_m1000_w100p0_int_neg", "HpseudoToTTTo1L1Nu2J_m365_w36p5_res", "HpseudoToTTTo1L1Nu2J_m400_w40p0_res", "HpseudoToTTTo1L1Nu2J_m500_w50p0_res", "HpseudoToTTTo1L1Nu2J_m600_w60p0_res", "HpseudoToTTTo1L1Nu2J_m800_w80p0_res", "HpseudoToTTTo1L1Nu2J_m1000_w100p0_res", "HpseudoToTTTo1L1Nu2J_m365_w36p5_int_pos", "HpseudoToTTTo1L1Nu2J_m400_w40p0_int_pos", "HpseudoToTTTo1L1Nu2J_m500_w50p0_int_pos", "HpseudoToTTTo1L1Nu2J_m600_w60p0_int_pos", "HpseudoToTTTo1L1Nu2J_m800_w80p0_int_pos", "HpseudoToTTTo1L1Nu2J_m1000_w100p0_int_pos", "HpseudoToTTTo1L1Nu2J_m365_w36p5_int_neg", "HpseudoToTTTo1L1Nu2J_m400_w40p0_int_neg", "HpseudoToTTTo1L1Nu2J_m500_w50p0_int_neg", "HpseudoToTTTo1L1Nu2J_m600_w60p0_int_neg", "HpseudoToTTTo1L1Nu2J_m800_w80p0_int_neg", "HpseudoToTTTo1L1Nu2J_m1000_w100p0_int_neg", "HscalarToTTTo1L1Nu2J_m365_w91p25_res", "HscalarToTTTo1L1Nu2J_m400_w100p0_res", "HscalarToTTTo1L1Nu2J_m500_w125p0_res", "HscalarToTTTo1L1Nu2J_m600_w150p0_res", "HscalarToTTTo1L1Nu2J_m800_w200p0_res", "HscalarToTTTo1L1Nu2J_m1000_w250p0_res", "HscalarToTTTo1L1Nu2J_m365_w91p25_int_pos", "HscalarToTTTo1L1Nu2J_m400_w100p0_int_pos", "HscalarToTTTo1L1Nu2J_m500_w125p0_int_pos", "HscalarToTTTo1L1Nu2J_m600_w150p0_int_pos", "HscalarToTTTo1L1Nu2J_m800_w200p0_int_pos", "HscalarToTTTo1L1Nu2J_m1000_w250p0_int_pos", "HscalarToTTTo1L1Nu2J_m365_w91p25_int_neg", "HscalarToTTTo1L1Nu2J_m400_w100p0_int_neg", "HscalarToTTTo1L1Nu2J_m500_w125p0_int_neg", "HscalarToTTTo1L1Nu2J_m600_w150p0_int_neg", "HscalarToTTTo1L1Nu2J_m800_w200p0_int_neg", "HscalarToTTTo1L1Nu2J_m1000_w250p0_int_neg", "HpseudoToTTTo1L1Nu2J_m365_w91p25_res", "HpseudoToTTTo1L1Nu2J_m400_w100p0_res", "HpseudoToTTTo1L1Nu2J_m500_w125p0_res", "HpseudoToTTTo1L1Nu2J_m600_w150p0_res", "HpseudoToTTTo1L1Nu2J_m800_w200p0_res", "HpseudoToTTTo1L1Nu2J_m1000_w250p0_res", "HpseudoToTTTo1L1Nu2J_m365_w91p25_int_pos", "HpseudoToTTTo1L1Nu2J_m400_w100p0_int_pos", "HpseudoToTTTo1L1Nu2J_m500_w125p0_int_pos", "HpseudoToTTTo1L1Nu2J_m600_w150p0_int_pos", "HpseudoToTTTo1L1Nu2J_m800_w200p0_int_pos", "HpseudoToTTTo1L1Nu2J_m1000_w250p0_int_pos", "HpseudoToTTTo1L1Nu2J_m365_w91p25_int_neg", "HpseudoToTTTo1L1Nu2J_m400_w100p0_int_neg", "HpseudoToTTTo1L1Nu2J_m500_w125p0_int_neg", "HpseudoToTTTo1L1Nu2J_m600_w150p0_int_neg", "HpseudoToTTTo1L1Nu2J_m800_w200p0_int_neg", "HpseudoToTTTo1L1Nu2J_m1000_w250p0_int_neg", "HscalarToTTTo1L1Nu2J_m365_w9p125_res", "HscalarToTTTo1L1Nu2J_m400_w10p0_res", "HscalarToTTTo1L1Nu2J_m500_w12p5_res", "HscalarToTTTo1L1Nu2J_m600_w15p0_res", "HscalarToTTTo1L1Nu2J_m800_w20p0_res", "HscalarToTTTo1L1Nu2J_m1000_w25p0_res", "HscalarToTTTo1L1Nu2J_m365_w9p125_int_pos", "HscalarToTTTo1L1Nu2J_m400_w10p0_int_pos", "HscalarToTTTo1L1Nu2J_m500_w12p5_int_pos", "HscalarToTTTo1L1Nu2J_m600_w15p0_int_pos", "HscalarToTTTo1L1Nu2J_m800_w20p0_int_pos", "HscalarToTTTo1L1Nu2J_m1000_w25p0_int_pos", "HscalarToTTTo1L1Nu2J_m365_w9p125_int_neg", "HscalarToTTTo1L1Nu2J_m400_w10p0_int_neg", "HscalarToTTTo1L1Nu2J_m500_w12p5_int_neg", "HscalarToTTTo1L1Nu2J_m600_w15p0_int_neg", "HscalarToTTTo1L1Nu2J_m800_w20p0_int_neg", "HscalarToTTTo1L1Nu2J_m1000_w25p0_int_neg", "HpseudoToTTTo1L1Nu2J_m365_w9p125_res", "HpseudoToTTTo1L1Nu2J_m400_w10p0_res", "HpseudoToTTTo1L1Nu2J_m500_w12p5_res", "HpseudoToTTTo1L1Nu2J_m600_w15p0_res", "HpseudoToTTTo1L1Nu2J_m800_w20p0_res", "HpseudoToTTTo1L1Nu2J_m1000_w25p0_res", "HpseudoToTTTo1L1Nu2J_m365_w9p125_int_pos", "HpseudoToTTTo1L1Nu2J_m400_w10p0_int_pos", "HpseudoToTTTo1L1Nu2J_m500_w12p5_int_pos", "HpseudoToTTTo1L1Nu2J_m600_w15p0_int_pos", "HpseudoToTTTo1L1Nu2J_m800_w20p0_int_pos", "HpseudoToTTTo1L1Nu2J_m1000_w25p0_int_pos", "HpseudoToTTTo1L1Nu2J_m365_w9p125_int_neg", "HpseudoToTTTo1L1Nu2J_m400_w10p0_int_neg", "HpseudoToTTTo1L1Nu2J_m500_w12p5_int_neg", "HpseudoToTTTo1L1Nu2J_m600_w15p0_int_neg", "HpseudoToTTTo1L1Nu2J_m800_w20p0_int_neg", "HpseudoToTTTo1L1Nu2J_m1000_w25p0_int_neg", "RSGluonToTT_M-500", "RSGluonToTT_M-1000", "RSGluonToTT_M-1500", "RSGluonToTT_M-2000", "RSGluonToTT_M-2500", "RSGluonToTT_M-3000", "RSGluonToTT_M-3500", "RSGluonToTT_M-4000", "RSGluonToTT_M-4500", "RSGluonToTT_M-5000", "RSGluonToTT_M-5500", "RSGluonToTT_M-6000", "ZPrimeToTT_M400_W40", "ZPrimeToTT_M500_W50", "ZPrimeToTT_M600_W60", "ZPrimeToTT_M700_W70", "ZPrimeToTT_M800_W80", "ZPrimeToTT_M900_W90", "ZPrimeToTT_M1000_W100", "ZPrimeToTT_M1200_W120", "ZPrimeToTT_M1400_W140", "ZPrimeToTT_M1600_W160", "ZPrimeToTT_M1800_W180", "ZPrimeToTT_M2000_W200", "ZPrimeToTT_M2500_W250", "ZPrimeToTT_M3000_W300", "ZPrimeToTT_M3500_W350", "ZPrimeToTT_M4000_W400", "ZPrimeToTT_M4500_W450", "ZPrimeToTT_M5000_W500", "ZPrimeToTT_M6000_W600",  "ZPrimeToTT_M7000_W700", "ZPrimeToTT_M8000_W800", "ZPrimeToTT_M9000_W900", "ZPrimeToTT_M400_W120", "ZPrimeToTT_M500_W150", "ZPrimeToTT_M600_W180", "ZPrimeToTT_M700_W210", "ZPrimeToTT_M800_W240", "ZPrimeToTT_M900_W270", "ZPrimeToTT_M1000_W300", "ZPrimeToTT_M1200_W360", "ZPrimeToTT_M1400_W420", "ZPrimeToTT_M1600_W480", "ZPrimeToTT_M1800_W540", "ZPrimeToTT_M2000_W600", "ZPrimeToTT_M2500_W750", "ZPrimeToTT_M3000_W900", "ZPrimeToTT_M3500_W1050", "ZPrimeToTT_M4000_W1200", "ZPrimeToTT_M4500_W1350", "ZPrimeToTT_M5000_W1500", "ZPrimeToTT_M6000_W1800", "ZPrimeToTT_M7000_W2100", "ZPrimeToTT_M8000_W2400", "ZPrimeToTT_M9000_W2700", "ZPrimeToTT_M400_W4", "ZPrimeToTT_M500_W5", "ZPrimeToTT_M600_W6", "ZPrimeToTT_M700_W7", "ZPrimeToTT_M800_W8", "ZPrimeToTT_M900_W9", "ZPrimeToTT_M1000_W10", "ZPrimeToTT_M1200_W12", "ZPrimeToTT_M1400_W14", "ZPrimeToTT_M1600_W16", "ZPrimeToTT_M1800_W18", "ZPrimeToTT_M2000_W20", "ZPrimeToTT_M2500_W25", "ZPrimeToTT_M3000_W30", "ZPrimeToTT_M3500_W35", "ZPrimeToTT_M4000_W40", "ZPrimeToTT_M4500_W45", "ZPrimeToTT_M5000_W50", "ZPrimeToTT_M6000_W60", "ZPrimeToTT_M7000_W70", "ZPrimeToTT_M8000_W80", "ZPrimeToTT_M9000_W90"};

    for(unsigned int i=0; i<names.size(); i++){
      if( ctx.get("dataset_version").find(names.at(i)) != std::string::npos ) sample_name = names.at(i);
    }
    if( (ctx.get("dataset_version").find("TTToHadronic") != std::string::npos) || (ctx.get("dataset_version").find("TTToSemiLeptonic") != std::string::npos) || (ctx.get("dataset_version").find("TTTo2L2Nu") != std::string::npos) ) sample_name = "TTbar";
    if( (ctx.get("dataset_version").find("WW") != std::string::npos) || (ctx.get("dataset_version").find("ZZ") != std::string::npos) || (ctx.get("dataset_version").find("WZ") != std::string::npos) ) sample_name = "Diboson";

    if(isMuon){
      TFile* f_btag2Dsf_muon = new TFile("/nfs/dust/cms/user/jabuschh/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_muon_"+year+".root");
      ratio_hist_muon = (TH2F*)f_btag2Dsf_muon->Get("N_Jets_vs_HT_" + sample_name);
      ratio_hist_muon->SetDirectory(0);
    }
    else if(!isMuon){
      TFile* f_btag2Dsf_ele = new TFile("/nfs/dust/cms/user/jabuschh/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_electron_"+year+".root");
      ratio_hist_ele = (TH2F*)f_btag2Dsf_ele->Get("N_Jets_vs_HT_" + sample_name);
      ratio_hist_ele->SetDirectory(0);
    }
  }
}

/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool ZprimeAnalysisModule_AzCorr_GenStudy::process(uhh2::Event& event){

  if(debug) cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug) cout << " run.event: " << event.run << ". " << event.event << endl;
  // Initialize reco flags with false
  event.set(h_is_zprime_reconstructed_chi2, false);
  event.set(h_is_zprime_reconstructed_correctmatch, false);
  event.set(h_chi2,-100);
  event.set(h_MET,-100);
  event.set(h_Mttbar,-100);
  event.set(h_lep1_pt,-100);
  event.set(h_lep1_eta,-100);
  event.set(h_ak4jet1_pt,-100);
  event.set(h_ak4jet1_eta,-100);
  event.set(h_ak8jet1_pt,-100);
  event.set(h_ak8jet1_eta,-100);
  event.set(h_NPV,-100);
  event.set(h_weight,-100);

  if(debug) cout<<"Initializing Spin Correlation set"<<endl;
  event.set(h_pt_hadTop, -10);      // pt of hadronic top
  // Plotting mass of ttbar system
  event.set(h_ttbar_mass_LabFrame, -10);

  // Plotting longitudinal boost of ttbar system
  event.set(h_ttbar_boost_LabFrame, -10);

  // Phi of lepton from leptonic leg
  event.set(h_phi_lep_LabFrame, -10);    
  event.set(h_phi_lep_CoMFrame, -10);
  event.set(h_phi_lep_helicityFrame, -10);
  event.set(h_phi_lep, -10);
  event.set(h_phi_lep_high, -10);
  event.set(h_phi_lep_low, -10); 
  event.set(h_phi_lep_highMass, -10);
  event.set(h_phi_lep_lowMass, -10);

  // Phi of b-jet from hadronic leg
  // event.set(h_phi_b_LabFrame, -10);
  // event.set(h_phi_b_CoMFrame, -10);
  // event.set(h_phi_b_helicityFrame, -10);
  // event.set(h_phi_b, -10);     
  // event.set(h_phi_b_high, -10);
  // event.set(h_phi_b_low, -10);
  // event.set(h_phi_b_highMass, -10);
  // event.set(h_phi_b_lowMass, -10);

  // Phi of q1 from W-decay of hadronic leg
  event.set(h_phi_hadTop_q1_LabFrame, -10);
  event.set(h_phi_hadTop_q1_CoMFrame, -10);
  event.set(h_phi_hadTop_q1_helicityFrame, -10);

  // Phi of q2 from W-decay of hadronic leg
  event.set(h_phi_hadTop_q2_LabFrame, -10);
  event.set(h_phi_hadTop_q2_CoMFrame, -10);
  event.set(h_phi_hadTop_q2_helicityFrame, -10);

  // Phi of less-energetic (in top's rest-frame) W daughter from hadronic leg
  event.set(h_phi_qlow, -10);     
  event.set(h_phi_qlow_high, -10);
  event.set(h_phi_qlow_low, -10);
  event.set(h_phi_qlow_highMass, -10);
  event.set(h_phi_qlow_lowMass, -10);

  // Sum of phi-coordinates
  event.set(h_sphi, -10);     
  event.set(h_sphi_low, -10); 
  event.set(h_sphi_high, -10);
  event.set(h_sphi_Mass1, -10);
  event.set(h_sphi_Mass2, -10);
  event.set(h_sphi_Mass3, -10);
  event.set(h_sphi_Mass4, -10);
  event.set(h_sphi_boost1, -10);
  event.set(h_sphi_boost2 , -10);
  event.set(h_sphi_boost3 , -10);
  event.set(h_sphi_boost4 , -10);
  
  // Difference of phi-coordinates
  event.set(h_dphi, -10);     
  event.set(h_dphi_low, -10); 
  event.set(h_dphi_high, -10);
  event.set(h_dphi_Mass1, -10);
  event.set(h_dphi_Mass2, -10);
  event.set(h_dphi_Mass3, -10);
  event.set(h_dphi_Mass4, -10);
  event.set(h_dphi_boost1 , -10);
  event.set(h_dphi_boost2 , -10);
  event.set(h_dphi_boost3 , -10);
  event.set(h_dphi_boost4 , -10);

  // Above variables now separated by charge of lepton in system
  event.set(h_phi_lepPlus, -10);
  event.set(h_phi_lepMinus, -10);
  // sphi_plus  is defined as (lepTop + hadTop) for positive leptons
  event.set(h_sphi_plus, -10);         
  event.set(h_sphi_plus_low, -10);     
  event.set(h_sphi_plus_high, -10);
  // dphi_plus is defined as (lepTop - hadTop) for positive leptons
  event.set(h_dphi_plus, -10);         
  event.set(h_dphi_plus_low, -10);     
  event.set(h_dphi_plus_high, -10);
  // sphi_minus is defined as (hadTop + lepTop) for negative leptons
  event.set(h_sphi_minus, -10);         
  event.set(h_sphi_minus_low, -10);     
  event.set(h_sphi_minus_high, -10);
  // dphi_minus is defined as (hadTop - lepTop) for negative leptons
  event.set(h_dphi_minus, -10);         
  event.set(h_dphi_minus_low, -10);     
  event.set(h_dphi_minus_high, -10);


  if(!event.isRealData){
    if(!SignSplit->passes(event)) return false;
  }

  // Run top-tagging
  if(ishotvr){
    TopTaggerHOTVR->process(event);
    hadronic_top->process(event);
  }
  else if(isdeepAK8){
    TopTaggerDeepAK8->process(event);
    hadronic_top->process(event);
  }

  fill_histograms(event, "Weights_Init");
  lumihists_Weights_Init->fill(event);

  if(!HEM_selection->passes(event)){
    if(!isMC) return false;
    else event.weight = event.weight*(1-0.64774715284); // calculated following instructions at https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis
  }
  fill_histograms(event, "Weights_HEM");

  // pileup weight
  PUWeight_module->process(event);
  if(debug) cout << "PUWeight: ok" << endl;
  fill_histograms(event, "Weights_PU");
  lumihists_Weights_PU->fill(event);

  // lumi weight
  LumiWeight_module->process(event);
  if(debug) cout << "LumiWeight: ok" << endl;
  fill_histograms(event, "Weights_Lumi");
  lumihists_Weights_Lumi->fill(event);

  // top pt reweighting
  //TopPtReweight_module->process(event);
  //if(debug) cout << "TopPtReweight: ok" << endl;
  //fill_histograms(event, "Weights_TopPt");
  //lumihists_Weights_TopPt->fill(event);

  // MC scale
  MCScale_module->process(event);
  if(debug) cout << "MCScale: ok" << endl;
  fill_histograms(event, "Weights_MCScale");
  lumihists_Weights_MCScale->fill(event);

  // Prefiring weights
  if(isMC){
    if (Prefiring_direction == "nominal") event.weight *= event.prefiringWeight;
    else if (Prefiring_direction == "up") event.weight *= event.prefiringWeightUp;
    else if (Prefiring_direction == "down") event.weight *= event.prefiringWeightDown;
  }
  fill_histograms(event, "Weights_Prefiring");

  // Write PSWeights from genInfo to own branch in output tree
  ps_weights->process(event);
  if(debug) cout << "Weights_PS: ok" << endl;
  fill_histograms(event, "Weights_PS");
  lumihists_Weights_PS->fill(event);

  // DeepAK8 TopTag SFs
  if(isdeepAK8) sf_toptag->process(event);
  if(debug) cout << "Weights_TopTag_SF: ok" << endl;
  fill_histograms(event, "Weights_TopTag_SF");

  //Clean muon collection with ID based on muon pT
  double muon_pt_high(55.);
  bool muon_is_low = false;
  bool muon_is_high = false;

  if(isMuon){
    vector<Muon>* muons = event.muons;
    for(unsigned int i=0; i<muons->size(); i++){
      if(event.muons->at(i).pt()<=muon_pt_high){
        muon_is_low = true;
      }
      else{
        muon_is_high = true;
      }
    }
  }
  sort_by_pt<Muon>(*event.muons);

  // Select exactly 1 muon and 0 electrons
  if(isMuon){
    if(!isEleTriggerMeasurement){ // For ele trigger SF do not veto additional electrons in muon channel
      if(!EleVeto_selection->passes(event)) return false;
    }
    if(muon_is_low){
      if(!NMuon1_selection->passes(event)) return false;
      muon_cleaner_low->process(event);
      if(!NMuon1_selection->passes(event)) return false;
      fill_histograms(event, "Muon1_LowPt");
    }
    if(muon_is_high){
      if(!NMuon1_selection->passes(event)) return false;
      muon_cleaner_high->process(event);
      if(!NMuon1_selection->passes(event)) return false;
      fill_histograms(event, "Muon1_HighPt");
    }
    if( !(muon_is_high || muon_is_low) ) return false;
    fill_histograms(event, "Muon1_Tot");
  }

  //Clean ele collection with ID based on ele pT
  double electron_pt_high(120.);
  bool ele_is_low = false;
  bool ele_is_high = false;

  if(isElectron || (isMuon && isEleTriggerMeasurement)){
    vector<Electron>* electrons = event.electrons;
    for(unsigned int i=0; i<electrons->size(); i++){
      if(abs(event.electrons->at(i).eta()) > 1.44 && abs(event.electrons->at(i).eta()) < 1.57) return false; // remove gap electrons in transition region between the barrel and endcaps of ECAL
      if(event.electrons->at(i).pt()<=electron_pt_high){
        ele_is_low = true;
      }
      else{
        ele_is_high = true;
      }
    }
  }
  sort_by_pt<Electron>(*event.electrons);

  // For ele trigger SF measurement
  if(isMuon && isEleTriggerMeasurement){
    if(ele_is_low){
      if(!NEle1_selection->passes(event)) return false;
      electron_cleaner_low->process(event);
      if(!NEle1_selection->passes(event)) return false;
      fill_histograms(event, "1Mu1Ele_LowPt");
    }
    if(ele_is_high){
      if(!NEle1_selection->passes(event)) return false;
      electron_cleaner_high->process(event);
      if(!NEle1_selection->passes(event)) return false;
      fill_histograms(event, "1Mu1Ele_HighPt");
    }
    if( !(ele_is_high || ele_is_low) ) return false;
    fill_histograms(event, "1Mu1Ele_Tot");
  }

  // Select exactly 1 electron and 0 muons
  if(isElectron){
    if(!MuonVeto_selection->passes(event)) return false;
    if(ele_is_low){
      if(!NEle1_selection->passes(event)) return false;
      electron_cleaner_low->process(event);
      if(!NEle1_selection->passes(event)) return false;
      fill_histograms(event, "Ele1_LowPt");
    }
    if(ele_is_high){
      if(!NEle1_selection->passes(event)) return false;
      electron_cleaner_high->process(event);
      if(!NEle1_selection->passes(event)) return false;
      fill_histograms(event, "Ele1_HighPt");
    }
    if( !(ele_is_high || ele_is_low)) return false;
    fill_histograms(event, "Ele1_Tot");
  }


  // apply electron id scale factors
  if(isMuon && !isEleTriggerMeasurement){
    sf_ele_id_dummy->process(event);
  }
  if(isMuon && isEleTriggerMeasurement){
    if(ele_is_low){
      sf_ele_id_low->process(event);
    }
    else if(ele_is_high){
      sf_ele_id_high->process(event);
    }
  }
  if(isElectron){
    if(ele_is_low){
      sf_ele_id_low->process(event);
    }
    else if(ele_is_high){
      sf_ele_id_high->process(event);
    }
    fill_histograms(event, "IdEle_SF");
  }


  // apply muon isolation scale factors (low pT only)
  if(isMuon){
    if(muon_is_low){
      sf_muon_iso_low->process(event);
    }
    else if(muon_is_high){
      sf_muon_iso_low_dummy->process(event);
    }
    fill_histograms(event, "IsoMuon_SF");
  }
  if(isElectron){
    sf_muon_iso_low_dummy->process(event);
  }
  // apply muon id scale factors
  if(isMuon){
    if(muon_is_low){
      sf_muon_id_low->process(event);
    }
    else if(muon_is_high){
      sf_muon_id_high->process(event);
    }
    fill_histograms(event, "IdMuon_SF");
  }
  if(isElectron){
    sf_muon_id_dummy->process(event);
  }

  // apply electron reco scale factors
  if(isMuon && !isEleTriggerMeasurement){
    sf_ele_reco_dummy->process(event);
  }
  if(isMuon && isEleTriggerMeasurement){
    sf_ele_reco->process(event);
  }
  if(isElectron){
    sf_ele_reco->process(event);
    fill_histograms(event, "RecoEle_SF");
  }

  // apply muon reco scale factors
  sf_muon_reco->process(event);
  fill_histograms(event, "MuonReco_SF");

  // Trigger MUON channel
  if(isMuon){
    // low pt
    if(muon_is_low){
      if(isUL16preVFP || isUL16postVFP){
        if(!(Trigger_mu_A_selection->passes(event) || Trigger_mu_B_selection->passes(event))) return false;
      }
      else{
        if(!Trigger_mu_A_selection->passes(event)) return false;
      }
    }
    // high pt
    if(muon_is_high){
      if(isUL16preVFP || isUL16postVFP){ // 2016
        if(!isMC){ //2016 DATA RunB
          if( event.run < 274889){
            if(!Trigger_mu_C_selection->passes(event)) return false;
          }else{ // DATA above Run B
            if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_D_selection->passes(event))) return false;
          }
        }else{ // 2016 MC RunB
          float runB_UL16_mu = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
          if( runB_UL16_mu < 0.1429){
            if(!Trigger_mu_C_selection->passes(event)) return false;
          }else{  // 2016 MC above RunB
            if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_D_selection->passes(event))) return false;
          }
        }
      }
      if(isUL17 && !isMC){ //2017 DATA
        if(event.run <= 299329){ //RunB
          if(!Trigger_mu_C_selection->passes(event)) return false;
        }else{
          if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_E_selection->passes(event) || Trigger_mu_F_selection->passes(event))) return false;
        }
      }
      if(isUL17 && isMC){ // 2017 MC
        float runB_mu = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if(runB_mu <= 0.1158){
          if(!Trigger_mu_C_selection->passes(event)) return false;
        }else{
          if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_E_selection->passes(event) || Trigger_mu_F_selection->passes(event))) return false;
        }
      }
      if(isUL18){ //2018
        if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_E_selection->passes(event) || Trigger_mu_F_selection->passes(event))) return false;
      }
    }
    fill_histograms(event, "TriggerMuon");
    lumihists_TriggerMuon->fill(event);
  }

  // Trigger ELECTRON channel
  if(isElectron){
    // low pt
    if(ele_is_low){
      if(isPhoton) return false;
      if(!Trigger_ele_A_selection->passes(event)) return false;
    }
    // high pt
    if(ele_is_high){
      // MC except UL17 that needs special treatment
      if(isMC && (isUL16preVFP || isUL16postVFP || isUL18) ){
        if(!(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event) || Trigger_ph_A_selection->passes(event))) return false;
      }
      if(isMC && isUL17){
        float runB_ele = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if(runB_ele <= 0.1158){ // in RunB (below runnumb 299329) Ele115 does not exist, use Ele35 instead. To apply randomly in MC if random numb < RunB percetage (11.58%, calculated by Christopher Matthies)
          if(!(Trigger_ele_A_selection->passes(event) || Trigger_ph_A_selection->passes(event))) return false;
        }else{
          if(!(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event) || Trigger_ph_A_selection->passes(event))) return false;
        }
      }
      if(!isMC){
        //DATA
        // 2016
        if(isUL16preVFP || isUL16postVFP){
          if(isPhoton){ // photon stream
            if(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event)) return false;
            if(!Trigger_ph_A_selection->passes(event)) return false;
          }else{ // electron stream
            if(!(Trigger_ele_A_selection->passes(event) ||  Trigger_ele_B_selection->passes(event) )) return false;
          }
        }
        // 2017
        if(isUL17){
          // below runnumb trigger Ele115 does not exist
          if(event.run <= 299329){
            if(isPhoton){ // photon stream
              if(Trigger_ele_A_selection->passes(event)) return false;
              if(!Trigger_ph_A_selection->passes(event)) return false;
            }else{ // electron stream
              if(!Trigger_ele_A_selection->passes(event)) return false;
            }
          }else{ // above runnumb with Ele115
            if(isPhoton){ // photon stream
              if(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event)) return false;
              if(!Trigger_ph_A_selection->passes(event)) return false;
            }else{ // electron stream
              if(! (Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event))) return false;
            }
          }
        }
        // 2018 photon & electron together: EGamma DATA
        if(isUL18){
          if(!(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event)|| Trigger_ph_A_selection->passes(event))) return false;
        }

      }
    }
    fill_histograms(event, "TriggerEle");
    lumihists_TriggerEle->fill(event);
  }


  // apply lepton trigger scale factors
  if(isMuon){
    if(muon_is_low){
      sf_muon_trigger_low->process(event);
    }
    if(muon_is_high){
      sf_muon_trigger_high->process(event);
    }
    fill_histograms(event, "TriggerMuon_SF");
  }
  if(isElectron){
    sf_muon_trigger_dummy->process(event);
  }


  if(!isEleTriggerMeasurement){
    if((event.muons->size()+event.electrons->size()) != 1) return false; //veto events without leptons or with too many
  }else{
    if( event.muons->size() != 1 && event.electrons->size() != 1) return false; //for ele trigger SF measurement: 1 ele and 1 mu
  }
  if(debug) cout<<"N leptons ok: Nelectrons="<<event.electrons->size()<<" Nmuons="<<event.muons->size()<<endl;

  if(isMuon && !isEleTriggerMeasurement && muon_is_high){
    if(!TwoDCut_selection->passes(event)) return false;
  }
  fill_histograms(event, "TwoDCut_Muon");
  lumihists_TwoDCut_Muon->fill(event);
  if(isElectron && ele_is_high){
    if(!TwoDCut_selection->passes(event)) return false;
  }
  fill_histograms(event, "TwoDCut_Ele");
  lumihists_TwoDCut_Ele->fill(event);

  if(isMuon && isEleTriggerMeasurement && (muon_is_high || ele_is_high)){
    if(!TwoDCut_selection->passes(event)) return false;
  }

  // match AK4 PUPPI to CHS
  AK4PuppiCHS_matching->process(event);
  h_CHSMatchHists->fill(event);

  if(!Jet1_selection->passes(event)) return false;
  if(debug) cout << "Jet1_selection: ok" << endl;
  fill_histograms(event, "Jet1");
  lumihists_Jet1->fill(event);

  if(!Jet2_selection->passes(event)) return false;
  if(debug) cout << "Jet2_selection: is ok" << endl;
  fill_histograms(event, "Jet2");
  lumihists_Jet2->fill(event);

  // MET selection
  if(!met_sel->passes(event)) return false;
  if(debug) cout << "MET: ok" << endl;
  fill_histograms(event, "MET");
  lumihists_MET->fill(event);
  if(isMuon){
    if(!htlep_sel->passes(event)) return false;
    fill_histograms(event, "HTlep");
    lumihists_HTlep->fill(event);
    if(debug) cout << "HTlep: ok" << endl;
  }

  //Fill histograms before BTagging SF - used to extract Custom BTag SF in (NJets,HT)
  h_CHSMatchHists_beforeBTagSF->fill(event);
  fill_histograms(event, "BeforeBtagSF");

  // btag shape sf (Ak4 chs jets)
  // new: using new modules, with PUPPI-CHS matching
  sf_btagging->process(event);

  h_CHSMatchHists_afterBTagSF->fill(event);
  fill_histograms(event, "AfterBtagSF");

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
  h_CHSMatchHists_after2DBTagSF->fill(event);
  fill_histograms(event, "AfterCustomBtagSF");

  // b-tagging: >= 1 b-tag medium WP (on matched CHS jet)
  if(!AK4PuppiCHS_BTagging->passes(event)) return false;
  fill_histograms(event, "Btags1");
  h_CHSMatchHists_afterBTag->fill(event);

  // Higher order corrections - EWK & QCD NLO
  NLOCorrections_module->process(event);
  fill_histograms(event, "NLOCorrections");

  //// apply ele trigger SF
  if(!isEleTriggerMeasurement) sf_ele_trigger->process(event);
  fill_histograms(event, "TriggerEle_SF");


  // build all possible ttbar candidates
  CandidateBuilder->process(event);
  fill_histograms(event, "TTbarCandidate");
  if(debug) cout << "CandidateBuilder: ok" << endl;

  // matching to gen-level ttbar - to extract chi2 parameters
  CorrectMatchDiscriminatorZprime->process(event);
  fill_histograms(event, "CorrectMatchDiscriminator");
  if(debug) cout << "CorrectMatchDiscriminatorZprime: ok" << endl;
  // if(sample.Contains("_blinded")){
  //   if(!BlindData_selection->passes(event)) return false;
  // }

  // select ttbar candidate with smallest chi2, fill Mtt hists
  Chi2DiscriminatorZprime->process(event);
  fill_histograms(event, "Chi2Discriminator");
  if(debug) cout << "Chi2DiscriminatorZprime: ok" << endl;

  // Variables for NN
  sort_by_pt<Jet>(*event.jets);
  Variables_module->process(event);
  fill_histograms(event, "NNInputsBeforeReweight");
  if(debug) cout << "NNInputsBeforeReweight: ok" << endl;

  // histograms for systematics
  //if(!isEleTriggerMeasurement) SystematicsModule->process(event);

  // filling Analysis Tree histos
  event.set(h_weight, event.weight);
  event.set(h_NPV, event.pvs->size());
  event.set(h_MET, event.met->pt());

  // Everything below this line is for the Spin Correlation studies -----------------------------------------------

  // fill ttbargen information
  ttgenprod->process(event);

  // Variable to access gen-particles
  TTbarGen ttbargen = event.get(h_ttbargen);

  // Make sure ttbar decays semileptonically
  if(ttbargen.IsSemiLeptonicDecay()){
    float pt_hadTop_thresh = 150;         // Define cut-variable as pt of hadTop for low/high regions
    // float ttbar_mass_thresh = 750;         // Define cut-variable as mass of ttbar-system for low/high regions

    //-------------------------------------------------- Start in LAB-frame --------------------------------------------------//

    // Plot pt of hadronic Top 
    LorentzVector Gen_HadTop = ttbargen.TopHad().v4();
    float pt_hadTop = Gen_HadTop.pt();
    if(pt_hadTop != -10) event.set(h_pt_hadTop, pt_hadTop);

    // // Define 4vectors of hadronic b
    // LorentzVector Gen_b = ttbargen.BHad().v4();
    // TLorentzVector hadTop_b;
    // // Set hadronic b 4vector using components
    // hadTop_b.SetPtEtaPhiE(Gen_b.pt(),Gen_b.eta(),Gen_b.phi(),Gen_b.energy());

    // Define 4vector of less energetic (in mother top's rest-frame) W daughter
    TLorentzVector hadTop_qlow;

    // Above vector TO BE SET AFTER TOP REST-FRAME

    // Define 4vectors of W jets
    LorentzVector Gen_q1 = ttbargen.Q1().v4();
    TLorentzVector hadTop_q1;
    LorentzVector Gen_q2 = ttbargen.Q2().v4();
    TLorentzVector hadTop_q2;

    // Set 4vectors of W jets using components
    hadTop_q1.SetPtEtaPhiE(Gen_q1.pt(),Gen_q1.eta(),Gen_q1.phi(),Gen_q1.energy());
    hadTop_q2.SetPtEtaPhiE(Gen_q2.pt(),Gen_q2.eta(),Gen_q2.phi(),Gen_q2.energy());

    // Define 4vectors of lepton
    LorentzVector Gen_Lep = ttbargen.ChargedLepton().v4();
    TLorentzVector lepTop_lep;
    // Set lepton 4vector using components
    lepTop_lep.SetPtEtaPhiE(Gen_Lep.pt(), Gen_Lep.eta(), Gen_Lep.phi(), Gen_Lep.energy());

    // Plot LAB FRAME phi-coordinates
    // float phi_b_LabFrame = hadTop_b.Phi();
    // if(phi_b_LabFrame != -10) event.set(h_phi_b_LabFrame, phi_b_LabFrame);
    // float phi_qlow_LabFrame = hadTop_qlow.Phi();// NOT SET UNTIL TOP REST-FRAME
    // if(phi_qlow_LabFrame != -10) event.set(h_phi_qlow_LabFrame, phi_qlow_LabFrame);// NOT SET UNTIL TOP REST-FRAME
    float phi_hadTop_q1_LabFrame = hadTop_q1.Phi();
    if(phi_hadTop_q1_LabFrame != -10) event.set(h_phi_hadTop_q1_LabFrame, phi_hadTop_q1_LabFrame);
    float phi_hadTop_q2_LabFrame = hadTop_q2.Phi();
    if(phi_hadTop_q2_LabFrame != -10) event.set(h_phi_hadTop_q2_LabFrame, phi_hadTop_q2_LabFrame);
    float phi_lep_LabFrame = lepTop_lep.Phi();
    if(phi_lep_LabFrame != -10) event.set(h_phi_lep_LabFrame, phi_lep_LabFrame);

    // Top vectors
    TLorentzVector PosTop;
    TLorentzVector NegTop;
    // Set top 4vectors using components
    PosTop.SetPtEtaPhiE(ttbargen.Top().v4().pt(), ttbargen.Top().v4().eta(), ttbargen.Top().v4().phi(), ttbargen.Top().v4().energy());
    NegTop.SetPtEtaPhiE(ttbargen.Antitop().v4().pt(), ttbargen.Antitop().v4().eta(), ttbargen.Antitop().v4().phi(), ttbargen.Antitop().v4().energy());

    // 4vector to represent ttbar system
    TLorentzVector ttbar(PosTop + NegTop);

    // Plotting mass of ttbar system
    event.set(h_ttbar_mass_LabFrame, ttbar.M());

    // Constructing longitudinal boost of ttbar system
    float numerator = fabs(PosTop.Pz() + NegTop.Pz());
    float denominator = PosTop.E() + NegTop.E();
    float boost = numerator/denominator;

    // Plotting longitudinal boost of ttbar system
    event.set(h_ttbar_boost_LabFrame, boost);

    //---------------------------------------------------------- Boost into CoM-frame ----------------------------------------------------------//
    if(debug) cout<<" Start First Boost"<<endl;

    // Boost into ttbar Center of Momentum configuration 
    lepTop_lep.Boost(-ttbar.BoostVector());
    // hadTop_b.Boost(-ttbar.BoostVector());
    // hadTop_qlow.Boost(-ttbar.BoostVector());// NOT SET UNTIL TOP REST-FRAME
    hadTop_q1.Boost(-ttbar.BoostVector());
    hadTop_q2.Boost(-ttbar.BoostVector());
    PosTop.Boost(-ttbar.BoostVector());
    NegTop.Boost(-ttbar.BoostVector());

    // Plot phi-coordinates in CoM frame
    // float phi_b_CoMFrame = hadTop_b.Phi();
    // event.set(h_phi_b_CoMFrame, phi_b_CoMFrame);
    // float phi_qlow_CoMFrame = hadTop_qlow.Phi();// NOT SET UNTIL TOP REST-FRAME
    // event.set(h_phi_qlow_CoMFrame, phi_qlow_CoMFrame);// NOT SET UNTIL TOP REST-FRAME
    float phi_hadTop_q1_CoMFrame = hadTop_q1.Phi();
    event.set(h_phi_hadTop_q1_CoMFrame, phi_hadTop_q1_CoMFrame);
    float phi_hadTop_q2_CoMFrame = hadTop_q2.Phi();
    event.set(h_phi_hadTop_q2_CoMFrame, phi_hadTop_q2_CoMFrame);
    float phi_lep_CoMFrame = lepTop_lep.Phi();
    event.set(h_phi_lep_CoMFrame, phi_lep_CoMFrame);

    // Cross-check back-to-back configuration of ttbar system
    if(debug) cout<<" Angle between tops is: "<< PosTop.Angle(NegTop.Vect())<<endl;
    //---------------------------------------------------------- Boost into CoM-frame ----------------------------------------------------------//

    //-------------- Rotate into Helicity Frame --------------//
    if(debug) cout<<"  Start Rotating"<<endl;

    // Rotation angles
    float angle1 = PosTop.Phi();
    float angle2 = PosTop.Theta();

    // Rotate tops and decay products about beam-line
    lepTop_lep.RotateZ(-1.*angle1);
    // hadTop_b.RotateZ(-1.*angle1);
    // hadTop_qlow.RotateZ(-1.*angle1);// NOT SET UNTIL TOP REST-FRAME
    hadTop_q1.RotateZ(-1.*angle1);
    hadTop_q2.RotateZ(-1.*angle1);
    PosTop.RotateZ(-1.*angle1);
    NegTop.RotateZ(-1.*angle1);
    // Rotate tops and decay products about y-zxis
    lepTop_lep.RotateY(-1.*angle2);
    // hadTop_b.RotateY(-1.*angle2);
    // hadTop_qlow.RotateY(-1.*angle2);// NOT SET UNTIL TOP REST-FRAME
    hadTop_q1.RotateY(-1.*angle2);
    hadTop_q2.RotateY(-1.*angle2);
    PosTop.RotateY(-1.*angle2);
    NegTop.RotateY(-1.*angle2);

    // Plot phi-coordinates in helicity-frame
    // float phi_b_helicityFrame = hadTop_b.Phi();
    // event.set(h_phi_b_helicityFrame, phi_b_helicityFrame);
    // float phi_qlow_helicityFrame = hadTop_qlow.Phi(); // NOT SET UNTIL TOP REST-FRAME
    // event.set(h_phi_qlow_helicityFrame, phi_qlow_helicityFrame); // NOT SET UNTIL TOP REST-FRAME
    float phi_hadTop_q1_helicityFrame = hadTop_q1.Phi();
    event.set(h_phi_hadTop_q1_helicityFrame, phi_hadTop_q1_helicityFrame);
    float phi_hadTop_q2_helicityFrame = hadTop_q2.Phi();
    event.set(h_phi_hadTop_q2_helicityFrame, phi_hadTop_q2_helicityFrame);
    float phi_lep_helicityFrame = lepTop_lep.Phi();
    event.set(h_phi_lep_helicityFrame, phi_lep_helicityFrame);

    //-------------- Rotate into Helicity Frame --------------//

    //--------------------------- Boost into ttbar rest-frame ---------------------------//
    if(debug) cout<<"   Start Second Boost"<<endl;
    // // Boost the top's to rest individually, bringing their children with them
    // // Decay products get boosted in opposite directions depending on their mother top

    // POSITIVE LEPTON CONFIGURATION
    if(ttbargen.ChargedLepton().charge() > 0){ // Positively charged lepton
      lepTop_lep.Boost(-PosTop.BoostVector()); // means lepton has Positive Top mother
      // hadTop_b.Boost(-NegTop.BoostVector()); // means b-jet has Negative Top mother
      // hadTop_qlow.Boost(-NegTop.BoostVector()); // means W-daughter has Negative Top mother// NOT SET UNTIL TOP REST-FRAME
      hadTop_q1.Boost(-NegTop.BoostVector());
      hadTop_q2.Boost(-NegTop.BoostVector());
    }

    // NEGATIVE LEPTON CONFIGURATION
    if(ttbargen.ChargedLepton().charge() < 0){ // Negatively charged lepton
      lepTop_lep.Boost(-NegTop.BoostVector()); // means lepton has Negative Top mother
      // hadTop_b.Boost(-PosTop.BoostVector()); // means b-jet has Positive Top mother
      // hadTop_qlow.Boost(-PosTop.BoostVector()); // means W-daughter has Positive Top mother// NOT SET UNTIL TOP REST-FRAME
      hadTop_q1.Boost(-PosTop.BoostVector());
      hadTop_q2.Boost(-PosTop.BoostVector());
    }
    //--------------------------- Boost into ttbar rest-frame ---------------------------//

    // Set 4vector of less-energetic W-daughter // Now we can compare their energies to determine which is less energetic in this rest-frame
    if(hadTop_q1.E() < hadTop_q2.E()){hadTop_qlow = hadTop_q1;}
    else{hadTop_qlow = hadTop_q2;}

    if(debug) cout<<"   Finished results are:"<<endl;
    // Define phi-coordinates from these twice boosted and once rotated 4vectors
    float phi_lep = lepTop_lep.Phi();
    event.set(h_phi_lep, phi_lep);
    if(debug) cout<<"   phi_lep is: "<< phi_lep <<endl;
    // float phi_b = hadTop_b.Phi();
    // event.set(h_phi_b, phi_b);
    // if(debug) cout<<"   phi_b is: "<< phi_b <<endl;
    float phi_qlow = hadTop_qlow.Phi();
    event.set(h_phi_qlow, phi_qlow);
    if(debug) cout<<"   phi_qlow is: "<< phi_qlow <<endl;

    // // Define sphi as sum and dphi as difference of phi's and plot (mixed charges)
    // Also apply mapping to both to keep original domain of [-pi, pi]
    float sphi = phi_lep + phi_qlow;
    // Map back into original domain if necessary
    if(sphi > TMath::Pi()) sphi = sphi - 2*TMath::Pi();
    if(sphi < -TMath::Pi()) sphi = sphi + 2*TMath::Pi();
    event.set(h_sphi, sphi);
    if(debug) cout<<"    sphi is: "<< sphi <<endl;
    float dphi = phi_lep - phi_qlow;
    // Map back into original domain if necessary
    if(dphi > TMath::Pi()) dphi = dphi - 2*TMath::Pi();
    if(dphi < -TMath::Pi()) dphi = dphi + 2*TMath::Pi();
    event.set(h_dphi, dphi);
    if(debug) cout<<"    dphi is: "<< dphi <<endl;

    // Plot dphi and sphi for high-pt ranges
    if(pt_hadTop > pt_hadTop_thresh){
      event.set(h_phi_lep_high, phi_lep);
      event.set(h_phi_qlow_high, phi_qlow);
      event.set(h_sphi_high, sphi);
      event.set(h_dphi_high, dphi);
    }
    // Plot dphi and sphi for low-pt ranges
    if(pt_hadTop < pt_hadTop_thresh){
      event.set(h_phi_lep_low, phi_lep);
      event.set(h_phi_qlow_low, phi_qlow);
      event.set(h_sphi_low, sphi);
      event.set(h_dphi_low, dphi);
    }

    // Plot dphi and sphi for various ranges of ttbar mass
    // 0 < Mass1 < 500
    if(ttbar.M() < 500){
      event.set(h_sphi_Mass1, sphi);
      event.set(h_dphi_Mass1, dphi);
    }
    // 500 < Mass2 < 750
    if(ttbar.M() > 500 &&  ttbar.M() < 750){
      event.set(h_sphi_Mass2, sphi);
      event.set(h_dphi_Mass2, dphi);
    }
    // 750 < Mass3 < 1500
    if(ttbar.M() > 750 &&  ttbar.M() < 1500){
      event.set(h_sphi_Mass3, sphi);
      event.set(h_dphi_Mass3, dphi);
    }
    // 1500 < Mass4
    if(ttbar.M() > 1500){
      event.set(h_sphi_Mass4, sphi);
      event.set(h_dphi_Mass4, dphi);
    }

    // Plot dphi and sphi for various ranges of longitudinal boost of ttbar system
    // 0 < boost1 < 0.3
    if(boost < 0.3){
      event.set(h_sphi_boost1, sphi);
      event.set(h_dphi_boost1, dphi);
    }
    // 0.3 < boost2 < 0.6
    if(boost > 0.3 &&  boost < 0.6){
      event.set(h_sphi_boost2, sphi);
      event.set(h_dphi_boost2, dphi);
    }
    // 0.6 < boost3 < 0.8
    if(boost > 0.6 &&  boost < 0.8){
      event.set(h_sphi_boost3, sphi);
      event.set(h_dphi_boost3, dphi);
    }
    // 0.8 < boost4 < 1.0
    if(boost > 0.8 &&  boost < 1.0){
      event.set(h_sphi_boost4, sphi);
      event.set(h_dphi_boost4, dphi);
    }

    // Positively charged leptons
    // if(ttbargen.ChargedLepton().charge() > 0){
    //   // Plot phi for positive leptons
    //   event.set(h_phi_lepPlus, phi_lep);
    //   // sphi_plus  is defined as (lep + had) for positive leptons
    //   // float sphi_plus = phi_lep + phi_b;
    //   float sphi_plus = phi_lep + phi_qlow;
    //   // Map back into original domain if necessary
    //   if(sphi_plus > TMath::Pi()) sphi_plus = sphi_plus - 2*TMath::Pi();
    //   if(sphi_plus < -TMath::Pi()) sphi_plus = sphi_plus + 2*TMath::Pi();
    //   event.set(h_sphi_plus, sphi_plus);
    //   if(pt_hadTop > pt_hadTop_thresh){event.set(h_sphi_plus_high, sphi_plus);}
    //   if(pt_hadTop < pt_hadTop_thresh){event.set(h_sphi_plus_low, sphi_plus);}
    //   // dphi_plus is defined as (lep - had) for positive leptons
    //   // float dphi_plus = phi_lep - phi_b;
    //   float dphi_plus = phi_lep - phi_qlow;
    //   // Map back into original domain if necessary
    //   if(dphi_plus > TMath::Pi()) dphi_plus = dphi_plus - 2*TMath::Pi();
    //   if(dphi_plus < -TMath::Pi()) dphi_plus = dphi_plus + 2*TMath::Pi();
    //   event.set(h_dphi_plus, dphi_plus);         
    //   if(pt_hadTop > pt_hadTop_thresh){event.set(h_dphi_plus_high, dphi_plus);}     
    //   if(pt_hadTop < pt_hadTop_thresh){event.set(h_dphi_plus_low, dphi_plus);}
    // }

    //Negatively charged leptons
  //   if(ttbargen.ChargedLepton().charge() < 0){
  //     // Plot phi for negative leptons
  //     event.set(h_phi_lepMinus, phi_lep);
  //     // sphi_minus is defined as (had + lep) for negative leptons
  //     // float sphi_minus = phi_b + phi_lep;
  //     float sphi_minus = phi_qlow + phi_lep;
  //     // Map back into original domain if necessary
  //     if(sphi_minus > TMath::Pi()) sphi_minus = sphi_minus - 2*TMath::Pi();
  //     if(sphi_minus < -TMath::Pi()) sphi_minus = sphi_minus + 2*TMath::Pi();
  //     event.set(h_sphi_minus, sphi_minus);         
  //     if(pt_hadTop > pt_hadTop_thresh){event.set(h_sphi_minus_high, sphi_minus);}     
  //     if(pt_hadTop < pt_hadTop_thresh){event.set(h_sphi_minus_low, sphi_minus);}
  //     // dphi_minus is defined as (had - lep) for negative leptons
  //     // float dphi_minus = phi_b - phi_lep;
  //     float dphi_minus = phi_qlow - phi_lep;
  //     // Map back into original domain if necessary
  //     if(dphi_minus > TMath::Pi()) dphi_minus = dphi_minus - 2*TMath::Pi();
  //     if(dphi_minus < -TMath::Pi()) dphi_minus = dphi_minus + 2*TMath::Pi();
  //     event.set(h_dphi_minus, dphi_minus);         
  //     if(pt_hadTop > pt_hadTop_thresh){event.set(h_dphi_minus_high, dphi_minus);}
  //     if(pt_hadTop < pt_hadTop_thresh){event.set(h_dphi_minus_low, dphi_minus);}
  //   }
 }

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule_AzCorr_GenStudy)
