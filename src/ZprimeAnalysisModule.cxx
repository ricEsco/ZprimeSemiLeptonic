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
#include <UHH2/ZprimeSemiLeptonic/include/TopMistagScaleFactor.h>

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

class ZprimeAnalysisModule : public ModuleBASE {

public:

  explicit ZprimeAnalysisModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:

  bool debug; // enable debug flag
  bool isEFT; // running EFTsample flag
  
  
  ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Declare unique_ptr's /////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // Top-tagging Scale Factors
  unique_ptr<AnalysisModule> sf_toptag;
  unique_ptr<AnalysisModule> sf_topmistag;

  // Reweighting modules
  unique_ptr<DeepAK8TopTagger> TopTaggerDeepAK8;
  unique_ptr<AnalysisModule>   PUWeight_module, LumiWeight_module, MCScale_module;
  unique_ptr<TopPtReweighting> TopPtReweight_module;
  // Prefiring weights applied manually
  unique_ptr<PSWeights> ps_weights;
  unique_ptr<AnalysisModule> NLOCorrections_module;

  // Lepton ID Cleaners
  unique_ptr<ElectronCleaner> electron_cleaner_low,electron_cleaner_high;
  unique_ptr<MuonCleaner> muon_cleaner_low, muon_cleaner_high;
  
  // Lepton *_Scale Factors
  unique_ptr<AnalysisModule> sf_ele_id_low, sf_ele_id_high;                       // electron ID
  unique_ptr<AnalysisModule> sf_muon_iso_stat_low;                                // muon ISO stat (low-pT ONLY)
  unique_ptr<AnalysisModule> sf_muon_iso_syst_low;                                // muon ISO syst (low-pT ONLY)
  unique_ptr<AnalysisModule> sf_muon_id_stat_low,  sf_muon_id_stat_high;          // muon ID stat
  unique_ptr<AnalysisModule> sf_muon_id_syst_low, sf_muon_id_syst_high;           // muon ID syst
  unique_ptr<AnalysisModule> sf_ele_reco;                                         // electron RECO
  unique_ptr<MuonRecoSF>     sf_muon_reco;                                        // muon RECO
  unique_ptr<AnalysisModule> sf_muon_trigger_stat_low, sf_muon_trigger_stat_high; // muon Trigger stat
  unique_ptr<AnalysisModule> sf_muon_trigger_syst_low, sf_muon_trigger_syst_high; // muon Trigger syst
  unique_ptr<AnalysisModule> sf_ele_trigger;                                      // electron Trigger

  // b-tagging SFs
  unique_ptr<AnalysisModule> sf_btagging;                                     

  // dummy modules to avoid set SF value errors for events where SFs are not applied (e.g. muon SFs in electron channel and vice versa)
  unique_ptr<AnalysisModule> sf_ele_id_dummy, sf_ele_reco_dummy;
  unique_ptr<AnalysisModule> sf_muon_iso_stat_low_dummy, sf_muon_iso_syst_low_dummy, sf_muon_id_stat_dummy, sf_muon_id_syst_dummy, sf_muon_trigger_stat_dummy, sf_muon_trigger_syst_dummy;
  
  // *_Selections
  unique_ptr<Selection> HEM_selection;
  unique_ptr<Selection> MuonVeto_selection, EleVeto_selection;
  unique_ptr<Selection> NMuon1_selection, NEle1_selection;
  unique_ptr<Selection> Trigger_ele_A_selection, Trigger_ele_B_selection, Trigger_ph_A_selection;
  unique_ptr<Selection> Trigger_mu_A_selection, Trigger_mu_B_selection, Trigger_mu_C_selection, Trigger_mu_D_selection, Trigger_mu_E_selection, Trigger_mu_F_selection;
  unique_ptr<Selection> TwoDCut_selection,TwoDCut_selection_low1, TwoDCut_selection_low2; 
  unique_ptr<Selection> Jet1_selection, Jet2_selection;
  unique_ptr<Selection> met_sel;
  unique_ptr<Selection> htlep_sel;
  unique_ptr<Selection> sel_1btag, sel_2btag;
  unique_ptr<Selection> TopTagVetoSelection;
  unique_ptr<Selection> DeltaEta_selection;

  // ttbar reconstruction
  unique_ptr<ZprimeCandidateBuilder> CandidateBuilder;                         // creates all possible permutations (candidates) of the ttbar system
  unique_ptr<ZprimeChi2Discriminator> Chi2DiscriminatorZprime;                 // extracts chi2(reco, avg) from reconstructed candidates and sets bestCandidate pointer to candidate with lowest chi2
  unique_ptr<ZprimeCorrectMatchDiscriminator> CorrectMatchDiscriminatorZprime; // extracts dr(gen, reco) from reconstructed candidates and sets bestCandidate pointer to candidate with lowest dr
  unique_ptr<Selection> Chi2_selection;                 // selects candidates with chi2 < chi2_max (30.)
  unique_ptr<Selection> Chi2CandidateMatched_selection; // selects candidates with chi2 < chi2_max and CorrectMatch discriminant dr < 10.
  unique_ptr<Selection> TTbarMatchable_selection;       // selects events where ttbar gen particles can be deltaR matched to corresponding reco objects
  unique_ptr<Selection> ZprimeTopTag_selection;         // selects events with AK8 top-tag was found with lepton outside cone
  
  // declare Variables_NN pointer
  unique_ptr<Variables_NN> Variables_module;
  // declare SprinCorrelations pointer
  unique_ptr<SpinCorrelations> SpinCorrelations_module;


  // Do not use the following modules-- might delete these
  unique_ptr<HOTVRTopTagger> TopTaggerHOTVR; // don't use HOTVR-- might delete
  unique_ptr<AnalysisModule> hadronic_top;   // only used for HOTVR-- might delete


  ///////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Declare handles ///////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////

  Event::Handle< std::vector<Jet> > h_CHSjets_matched;  // Collection of CHS matched jets
  Event::Handle< std::vector<TopJet> > h_DeepAK8TopTags;  // Collection of DeepAK8TopTagged jets

  // ttbar reconstruction
  Event::Handle<float> h_chi2;   
  Event::Handle<bool> h_is_zprime_reconstructed_chi2;
  Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;

  Event::Handle<bool> h_is_zprime_reconstructed_correctmatch;
  Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateCorrectMatch;

  // Gen-level variables
  Event::Handle<TTbarGen> h_ttbargen;
  unique_ptr<TTbarGenProducer> ttgenprod;

  // Lumi histograms
  std::unique_ptr<Hists> 
  lumihists_Weights_Init, lumihists_Weights_PU, lumihists_Weights_Lumi, lumihists_Weights_TopPt, lumihists_Weights_MCScale, lumihists_Weights_PS, 
  lumihists_Muon1_LowPt, lumihists_Muon1_HighPt, lumihists_Ele1_LowPt, lumihists_Ele1_HighPt, 
  lumihists_TriggerMuon, lumihists_TriggerEle, 
  lumihists_TwoDCut_Muon, lumihists_TwoDCut_Ele, 
  lumihists_Jet1, lumihists_Jet2, 
  lumihists_MET, 
  lumihists_HTlep,
  lumihists_TwoDCut_Muon_LowPt,
  lumihists_Chi2;

  // PUPPI CHS match module
  unique_ptr<PuppiCHS_matching> AK4PuppiCHS_matching;
  // b-tagging (on CHS matched jets) module
  unique_ptr<Selection> AK4PuppiCHS_BTagging;

  // Histograms with matched CHS jets
  unique_ptr<Hists> h_CHSMatchHists;
  unique_ptr<Hists> h_CHSMatchHists_beforeBTagSF;
  unique_ptr<Hists> h_CHSMatchHists_afterBTagSF;
  unique_ptr<Hists> h_CHSMatchHists_after2DBTagSF;
  unique_ptr<Hists> h_CHSMatchHists_afterBTag;


  /////////////////////////////////////////////////////////////
  /////////////////////// Configuration ///////////////////////
  /////////////////////////////////////////////////////////////
  bool isMC, ishotvr, isdeepAK8;
  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;
  bool isMuon, isElectron;
  bool isPhoton;
  bool isEleTriggerMeasurement;

  int runnr_oldtriggers = 299368;

  string Sys_PU, Prefiring_direction, Sys_TopPt_a, Sys_TopPt_b;
  TString sample;
  TString year, channel;

  TH2F *ratio_hist_muon;
  TH2F *ratio_hist_ele;
};

void ZprimeAnalysisModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag + "_Skimming";
    mytag = tag + "_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicHists(ctx,mytag));
  }
}

void ZprimeAnalysisModule::fill_histograms(uhh2::Event& event, string tag){
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

ZprimeAnalysisModule::ZprimeAnalysisModule(uhh2::Context& ctx){
  debug = false;

  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }

  /////////////////////////////////////////////////
  ///////////////// Configuration /////////////////

  // dataset type
  isMC = (ctx.get("dataset_type") == "MC");
  isPhoton = (ctx.get("dataset_version").find("SinglePhoton") != std::string::npos);

  // Access to gen-level particles
  if(isMC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  // Jet collection
  ishotvr = (ctx.get("is_hotvr") == "true");
  isdeepAK8 = (ctx.get("is_deepAK8") == "true");
  TString mode = "hotvr";
  if(isdeepAK8) mode = "deepAK8";

  // Year
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

  // Running script as electron trigger measurement
  isEleTriggerMeasurement = (ctx.get("is_EleTriggerMeasurement") == "true");

  // Lepton channel
  isMuon = false; isElectron = false;
  channel = ctx.get("channel");
  if(channel == "muon") isMuon = true;
  if(channel == "electron") isElectron = true;

  // Systematics
  Sys_PU = ctx.get("Sys_PU");                     // Pileup systematic
  Prefiring_direction = ctx.get("Sys_prefiring"); // L1 prefiring systematic, mainly for high-energy and forward jets in Forward ECAL
  Sys_TopPt_a = ctx.get("Systematic_TopPt_a");    // a parameter for TopPt reweighting systematic
  Sys_TopPt_b = ctx.get("Systematic_TopPt_b");    // b parameter for TopPt reweighting systematic

  /// Important selection values
  double MET_cut, HT_lep_cut; // based on lepton channel, set below
  double TwoD_dr = 0.4;       // 2D-cut deltaR(lep, jet)
  double TwoD_ptrel = 25.;    // 2D-cut relative pT(lep, jet)
  double jet1_pt(50.);        // pT threshold for leading AK4 jet
  double jet2_pt;             // based on lepton channel, set below
  double chi2_max(30.);       // max chi2 for reconstructed candidates to be accepted


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// LEPTON STUFF ///////////////////////////////////////////

  // set lepton pT thresholds for high/low categories
  double electron_pt_low;
  if(isUL17){electron_pt_low = 38.;} // UL17 ele trigger (HLT_Ele35WPTight _Gsf) threshold is 35 GeV
  else{electron_pt_low = 35.;}
  double electron_pt_high(120.);
  double muon_pt_low(30.);
  double muon_pt_high(55.);

  // set IDs for lepton cleaning
  ElectronId eleID_low  = ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp80);
  const ElectronId electronID_low(AndId<Electron>(PtEtaSCCut(electron_pt_low, 2.5), eleID_low));
  electron_cleaner_low.reset(new ElectronCleaner(electronID_low));
  ElectronId eleID_high = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp80);
  const ElectronId electronID_high(AndId<Electron>(PtEtaSCCut(electron_pt_high, 2.5), eleID_high));
  electron_cleaner_high.reset(new ElectronCleaner(electronID_high));

  MuonId     muID_low   = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  const MuonId muonID_low(AndId<Muon>(PtEtaCut(muon_pt_low, 2.4), muID_low));
  muon_cleaner_low.reset(new MuonCleaner(muonID_low));
  MuonId     muID_high  = MuonID(Muon::CutBasedIdGlobalHighPt);
  const MuonId muonID_high(AndId<Muon>(PtEtaCut(muon_pt_high, 2.4), muID_high));
  muon_cleaner_high.reset(new MuonCleaner(muonID_high));
  
  // define Triggers
  string trigger_mu_A, trigger_mu_B, trigger_mu_C, trigger_mu_D, trigger_mu_E, trigger_mu_F;
  string trigger_ele_A, trigger_ele_B;
  string trigger_ph_A;

  // set Triggers and jet2 pT, MET, and HTlep cuts
  if(isElectron){
    trigger_ele_B = "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*";
    if(isUL16preVFP || isUL16postVFP){trigger_ele_A = "HLT_Ele27_WPTight_Gsf_v*";}
    if(isUL17){trigger_ele_A = "HLT_Ele35_WPTight_Gsf_v*";}
    if(isUL18){trigger_ele_A = "HLT_Ele32_WPTight_Gsf_v*";}
    if(isUL16preVFP || isUL16postVFP){trigger_ph_A = "HLT_Photon175_v*";}
    else{trigger_ph_A = "HLT_Photon200_v*";}
    jet2_pt = 40;
    MET_cut = 60;
    HT_lep_cut = 0;
  }
  if(isMuon){
    if(isUL17){trigger_mu_A = "HLT_IsoMu27_v*";}
    else{trigger_mu_A = "HLT_IsoMu24_v*";}
    trigger_mu_B = "HLT_IsoTkMu24_v*";
    trigger_mu_C = "HLT_Mu50_v*";
    trigger_mu_D = "HLT_TkMu50_v*";
    trigger_mu_E = "HLT_OldMu100_v*";
    trigger_mu_F = "HLT_TkMu100_v*";
    jet2_pt = 50;
    MET_cut = 70;
    HT_lep_cut = 0;
  }
  
  //////////////////////////////////////////// Jets ////////////////////////////////////////////
  // JetId used by PuppiCHS_matching()
  const JetPFID jetID_CHS(JetPFID::WP_TIGHT_CHS);

  // b-tagging parameters
  BTag::algo btag_algo = BTag::DEEPJET;     // algorithm used, e.g. DeepCSV, DeepJet, CSVv2, etc.
  BTag::wp btag_wp = BTag::WP_MEDIUM;       // working point
  JetId id_btag = BTag(btag_algo, btag_wp); // b-tagging JetId(algorithm, working point)

  // Top pT Reweighting parameters
  double a_toppt = 0.0615;  // par a TopPt Reweighting
  double b_toppt = -0.0005; // par b TopPt Reweighting
  // const TopJetId toptagID = AndId<TopJet>(HOTVRTopTag(0.8, 140.0, 220.0, 50.0), Tau32Groomed(0.56));


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////// POINTER REASSIGNING //////////////////////////////////////////////////////////

  /////////////////////////////////////// MC Reweighting modules ///////////////////////////////////////
  LumiWeight_module.reset(new MCLumiWeight(ctx));
  PUWeight_module.reset(new MCPileupReweight(ctx, Sys_PU));
  TopPtReweight_module.reset(new TopPtReweighting(ctx, a_toppt, b_toppt, Sys_TopPt_a, Sys_TopPt_b, ""));
  MCScale_module.reset(new MCScaleVariation(ctx));
  ps_weights.reset(new PSWeights(ctx));
  NLOCorrections_module.reset(new NLOCorrections(ctx));

  ///////////////////////////////// Lepton scale factors: see UHH2/common/include/LeptonScaleFactors.h /////////////////////////////////
  // muon_iso: stat/sys x low pT
  sf_muon_iso_stat_low.reset(new uhh2::MuonIsoScaleFactors_stat(ctx, Muon::Selector::PFIsoTight, Muon::Selector::CutBasedIdTight, true));
  sf_muon_iso_syst_low.reset(new uhh2::MuonIsoScaleFactors_syst(ctx, Muon::Selector::PFIsoTight, Muon::Selector::CutBasedIdTight, true));
  // muon_id: stat/sys x low/high pT
  sf_muon_id_stat_low.reset(new uhh2::MuonIdScaleFactors_stat(ctx, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_stat_high.reset(new uhh2::MuonIdScaleFactors_stat(ctx, Muon::Selector::CutBasedIdGlobalHighPt, true));
  sf_muon_id_syst_low.reset(new uhh2::MuonIdScaleFactors_syst(ctx, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_syst_high.reset(new uhh2::MuonIdScaleFactors_syst(ctx, Muon::Selector::CutBasedIdGlobalHighPt, true));
  // muon_trigger: stat/sys x low/high pT
  sf_muon_trigger_stat_low.reset(new uhh2::MuonTriggerScaleFactors_stat(ctx, false, true));
  sf_muon_trigger_stat_high.reset(new uhh2::MuonTriggerScaleFactors_stat(ctx, true, false));
  sf_muon_trigger_syst_low.reset(new uhh2::MuonTriggerScaleFactors_syst(ctx, false, true));
  sf_muon_trigger_syst_high.reset(new uhh2::MuonTriggerScaleFactors_syst(ctx, true, false));
  // muon reco
  sf_muon_reco.reset(new MuonRecoSF(ctx));

  // electron_id: stat/sys x low/high pT
  sf_ele_id_low.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_iso_V2_wp80, true));
  sf_ele_id_high.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_noIso_V2_wp80, true));
  // electron reco
  sf_ele_reco.reset(new uhh2::ElectronRecoScaleFactors(ctx, false, true));

  // apply electron trigger SFs only for non-trigger-measurement analyses
  if(!isEleTriggerMeasurement) sf_ele_trigger.reset( new uhh2::ElecTriggerSF(ctx, "central", "eta_ptbins", year) );

  // dummies (needed to aviod set value errors)
  sf_muon_iso_stat_low_dummy.reset(new uhh2::MuonIsoScaleFactors_stat(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_id_stat_dummy.reset(new uhh2::MuonIdScaleFactors_stat(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_trigger_stat_dummy.reset(new uhh2::MuonTriggerScaleFactors_stat(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_iso_syst_low_dummy.reset(new uhh2::MuonIsoScaleFactors_syst(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_id_syst_dummy.reset(new uhh2::MuonIdScaleFactors_syst(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_trigger_syst_dummy.reset(new uhh2::MuonTriggerScaleFactors_syst(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_ele_id_dummy.reset(new uhh2::ElectronIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_ele_reco_dummy.reset(new uhh2::ElectronRecoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));

  //////////////////////////////////////////// Jets ////////////////////////////////////////////
    
  // Top Taggers
  TopTaggerHOTVR.reset(new HOTVRTopTagger(ctx));
  TopTaggerDeepAK8.reset(new DeepAK8TopTagger(ctx));

  // Top-tags
  hadronic_top.reset(new HadronicTop(ctx)); // only used for HOTVR-- might delete
  // sf_toptag.reset(new HOTVRScaleFactor(ctx, toptagID, ctx.get("Sys_TopTag", "nominal"), "HadronicTop", "TopTagSF", "HOTVRTopTagSFs")); // don't use HOTVR-- might delete
  // Top-tag SF's
  sf_toptag.reset(new TopTagScaleFactor(ctx));
  sf_topmistag.reset(new TopMistagScaleFactor(ctx));

  // b-tag SF
  sf_btagging.reset(new MCBTagDiscriminantReweighting(ctx, BTag::algo::DEEPJET, "CHS_matched"));

  ////////////////////////////// Selection modules //////////////////////////////

  // NLeptons
  MuonVeto_selection.reset(new NMuonSelection(0, 0));    // ==0 muons
  EleVeto_selection.reset(new NElectronSelection(0, 0)); // ==0 electrons
  NMuon1_selection.reset(new NMuonSelection(1, 1));      // ==1 muon
  NEle1_selection.reset(new NElectronSelection(1, 1));   // ==1 electron

  // Lepton triggers
  Trigger_ele_A_selection.reset(new TriggerSelection(trigger_ele_A));
  Trigger_ele_B_selection.reset(new TriggerSelection(trigger_ele_B));
  Trigger_ph_A_selection.reset(new TriggerSelection(trigger_ph_A));

  Trigger_mu_A_selection.reset(new TriggerSelection(trigger_mu_A));
  Trigger_mu_B_selection.reset(new TriggerSelection(trigger_mu_B));
  Trigger_mu_C_selection.reset(new TriggerSelection(trigger_mu_C));
  Trigger_mu_D_selection.reset(new TriggerSelection(trigger_mu_D));
  Trigger_mu_E_selection.reset(new TriggerSelection(trigger_mu_E));
  Trigger_mu_F_selection.reset(new TriggerSelection(trigger_mu_F));

  // 2D-cuts
  TwoDCut_selection.reset(new TwoDCut(TwoD_dr, TwoD_ptrel));
  TwoDCut_selection_low1.reset(new TwoDCut(0.3, 10.));
  TwoDCut_selection_low2.reset(new TwoDCut(0.3, 15.));

  // Jets
  if(isUL16preVFP || isUL16postVFP){ 
    Jet1_selection.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.4))));
    Jet2_selection.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.4))));
  }
  else{
    Jet1_selection.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.5))));
    Jet2_selection.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.5))));
  }

  // b-tagging modules (require >= 1 or >=2 btagged jets in event)
  sel_1btag.reset(new NJetSelection(1, -1, id_btag));
  sel_2btag.reset(new NJetSelection(2,-1, id_btag));

  // MET and HTlep
  met_sel.reset(new METCut  (MET_cut   , uhh2::infinity));
  htlep_sel.reset(new HTlepCut(HT_lep_cut, uhh2::infinity));

  // top-tag veto (require <=1 top-tagged jet in event)
  TopTagVetoSelection.reset(new TopTag_VetoSelection(ctx, mode));

  // Cut on DeltaEta(j1,j2) < 3 to reduce QCD spikes
  DeltaEta_selection.reset(new DeltaEtaSelection());

 // Chi2 for reconstructed ttbar candidates
  Chi2_selection.reset(new Chi2Cut(ctx, 0., chi2_max));
  TTbarMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
  Chi2CandidateMatched_selection.reset(new Chi2CandidateMatchedSelection(ctx));
  ZprimeTopTag_selection.reset(new ZprimeTopTagSelection(ctx));

  // HEM issue in 2018, veto on leptons and jets
  HEM_selection.reset(new HEMSelection(ctx));



  // Variables for DNN
  Variables_module.reset(new Variables_NN(ctx, mode)); // NEED TO RUN Variables_module->process() else job will crash
  SpinCorrelations_module.reset(new SpinCorrelations(ctx, mode)); // module to calculate spin correlation variables
  // if(!isEleTriggerMeasurement) SystematicsModule.reset(new ZprimeSemiLeptonicSystematicsModule(ctx));

  

  // ttbar reconstruction candidate builder
  CandidateBuilder.reset(new ZprimeCandidateBuilder(ctx, mode));

  // ttbar reconstruction discriminators
  Chi2DiscriminatorZprime.reset(new ZprimeChi2Discriminator(ctx)); // chi2 discriminator
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  CorrectMatchDiscriminatorZprime.reset(new ZprimeCorrectMatchDiscriminator(ctx)); // correctMatch discriminator
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
  h_BestZprimeCandidateCorrectMatch = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestCorrectMatch");


  // PUPPI CHS match modules & hists
  AK4PuppiCHS_matching.reset(new PuppiCHS_matching(ctx)); // match AK4 PUPPI jets to AK4 CHS jets for b-tagging
  AK4PuppiCHS_BTagging.reset(new PuppiCHS_BTagging(ctx)); // b-tagging on matched CHS jets
  h_CHSMatchHists.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch"));
  h_CHSMatchHists_beforeBTagSF.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_beforeBTagSF"));
  h_CHSMatchHists_afterBTagSF.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_afterBTagSF"));
  h_CHSMatchHists_after2DBTagSF.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_after2DBTagSF"));
  h_CHSMatchHists_afterBTag.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_afterBTag"));

  // Lumi hists
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
  lumihists_TwoDCut_Muon_LowPt.reset(new LuminosityHists(ctx, "Lumi_TwoDCut_Muon_LowPt"));
  lumihists_Chi2.reset(new LuminosityHists(ctx, "Lumi_Chi2"));



  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////// Defining event output handles ////////////////////////////////////////

  // Related to Spin Correlation Variables
  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");       // Collection of CHS matched jets
  h_DeepAK8TopTags = ctx.get_handle< std::vector<TopJet>>("DeepAK8TopTags"); // Collection of DeepAK8TopTagged jets

  
  // Book histograms
  vector<string> histogram_tags = {
    "Weights_Init", 
    "Weights_HEM", 
    "Weights_PU", "Weights_Lumi", "Weights_TopPt", "Weights_MCScale", "Weights_Prefiring", "Weights_PS", 
    "Weights_TopTag_SF", "Weights_TopMistag_SF", 

    "Muon1_LowPt", "Muon1_HighPt", "Muon1_Tot", 
    "1Mu1Ele_LowPt", "1Mu1Ele_HighPt", "1Mu1Ele_Tot", // electronTriggerSF measurement in muon channel
    "Ele1_LowPt", "Ele1_HighPt", "Ele1_Tot", 
    
    "IdEle_SF", 
    "IsoMuon_SF", 
    "IdMuon_SF", 
    "RecoEle_SF", 
    "MuonReco_SF", 
    "TriggerMuon", "TriggerEle", 
    "TriggerMuon_SF", 

    "TwoDCut_Muon", "TwoDCut_Ele",
    "CHS_Before", "CHS_After", 
    "Jet1", "Jet2", 
    "MET", "HTlep", 
    "BeforeBtagSF", "AfterBtagSF", "AfterCustomBtagSF", "Btags1", 
    "NLOCorrections", 
    "TriggerEle_SF", 
    "TopTagVeto", 
    "DeltaEtaCut_FAIL", "DeltaEtaCut_PASS",
    "LowPtMuons_Before2Dcut", "LowPtMuons_FAIL2Dcut", "LowPtMuons_PASS2Dcut",
    "EndOfBaselineSelection",

    "MatchableBeforeChi2Cut", "CorrectMatchBeforeChi2Cut",
    "Chi2cut_FAIL", "Chi2cut_PASS", 
    "Matchable", "CorrectMatch",
    "Merged", "Matchable_Merged","CorrectMatch_Merged",
    "Resolved", "Matchable_Resolved","CorrectMatch_Resolved"
  };
  book_histograms(ctx, histogram_tags);
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Handling dataset sample name, includig EFT samples //////////////////////////////////////////////////////////////
  bool isEFT = false;

  if(isMC){
    TString sample_name = "";
    vector<TString> names = {"MC_EFT_Mttbar_0-700_UL17", "MC_EFT_Mttbar_700-900_UL17", "MC_EFT_Mttbar_900-Inf_UL17", "ST", "WJets", "DY", "QCD"};

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

    // Get b-tagging SF histograms (2D histograms of N_jets vs HT) for muon and electron channels
    if(isMuon){
      TFile* f_btag2Dsf = new TFile("/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_muon_"+year+".root");
      if(isEFT){
        // *** CHANGED *** For EFT muon: always use TTbar
        ratio_hist_muon = (TH2F*)f_btag2Dsf->Get("N_Jets_vs_HT_TTbar");
      }
      else {
        ratio_hist_muon = (TH2F*)f_btag2Dsf->Get("N_Jets_vs_HT_" + sample_name);
      }
      ratio_hist_muon->SetDirectory(0);
    }
    else if(!isMuon){
      TFile* f_btag2Dsf = new TFile("/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_electron_"+year+".root");
      if(isEFT){
        ratio_hist_ele = (TH2F*)f_btag2Dsf->Get("N_Jets_vs_HT_TTbar");
      } 
      else{
        ratio_hist_ele = (TH2F*)f_btag2Dsf->Get("N_Jets_vs_HT_" + sample_name);
      }
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

bool ZprimeAnalysisModule::process(uhh2::Event& event){

  if(debug) cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug) cout << " run.event: " << event.run << "." << event.event << endl;

  // Initialize reco flags with false
  event.set(h_is_zprime_reconstructed_chi2, false);
  event.set(h_is_zprime_reconstructed_correctmatch, false);
  
  // Process TTbarGen to set gen-level ttbar system variables
  if(isMC && ttgenprod) ttgenprod->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] TTbarGen processed: ok" << endl;

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
    if(debug) cout << "[ZprimeAnalysisModule] TopTaggerDeepAK8: ok" << endl;
  }

  // Initial weight (gen weight * pileup weight from preselection)
  fill_histograms(event, "Weights_Init");
  lumihists_Weights_Init->fill(event);

  // HEM veto for 2018 data and MC on eta and phi of jets and leptons in specific runs
  if(!HEM_selection->passes(event)){
    if(!isMC) return false;
    else event.weight = event.weight*(1-0.64774715284); // see https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis for calculation
  }
  if(debug) cout << "[ZprimeAnalysisModule] HEM Selection: passed" << endl;
  fill_histograms(event, "Weights_HEM");

  // pileup weight
  PUWeight_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] PUWeight: ok" << endl;
  fill_histograms(event, "Weights_PU");
  lumihists_Weights_PU->fill(event);
  
  // lumi weight
  LumiWeight_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] LumiWeight: ok" << endl;
  fill_histograms(event, "Weights_Lumi");
  lumihists_Weights_Lumi->fill(event);

  // top pt reweighting
  TopPtReweight_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] TopPtReweight: ok" << endl;
  fill_histograms(event, "Weights_TopPt");
  lumihists_Weights_TopPt->fill(event);

  // MC scale
  MCScale_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] MCScale: ok" << endl;
  fill_histograms(event, "Weights_MCScale");
  lumihists_Weights_MCScale->fill(event);

  // Prefiring weights
  if(isMC){
    if (Prefiring_direction == "nominal") event.weight *= event.prefiringWeight;
    else if (Prefiring_direction == "up") event.weight *= event.prefiringWeightUp;
    else if (Prefiring_direction == "down") event.weight *= event.prefiringWeightDown;
  }
  if(debug) cout << "[ZprimeAnalysisModule] Prefiring: ok" << endl;
  fill_histograms(event, "Weights_Prefiring");

  // Write PSWeights from genInfo to own branch in output tree
  ps_weights->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] Weights_PS: ok" << endl;
  fill_histograms(event, "Weights_PS");
  lumihists_Weights_PS->fill(event);

  // DeepAK8 TopTag SFs
  if(isdeepAK8) sf_toptag->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] TopTag_SF: ok" << endl;
  fill_histograms(event, "Weights_TopTag_SF");

  // DeepAK8 TopMistag SFs
  if(isdeepAK8) sf_topmistag->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] TopMistag_SF: ok" << endl;
  fill_histograms(event, "Weights_TopMistag_SF");




  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Start: Cleaning lepton collections with ID based on pT //////////////////////////////

  // muons
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

  // Selection: exactly 1 muon and 0 electrons
  if(isMuon){
    // For electronTrigger SF do not veto additional electrons in muon channel
    if(!isEleTriggerMeasurement){                         // skip for electronTrigger SF measurement
      if(!EleVeto_selection->passes(event)) return false; // require ==0 electrons
    }
    // low-pT muons
    if(muon_is_low){
      if(!NMuon1_selection->passes(event)) return false; // require ==1 muon
        muon_cleaner_low->process(event);                // clean low-pT muons
      if(!NMuon1_selection->passes(event)) return false; // require ==1 muon after cleaning
        fill_histograms(event, "Muon1_LowPt");
        lumihists_Muon1_LowPt->fill(event);
    }
    if(muon_is_high){
      if(!NMuon1_selection->passes(event)) return false; // require ==1 muon
      muon_cleaner_high->process(event);                 // clean high-pT muons
      if(!NMuon1_selection->passes(event)) return false; // require ==1 muon after cleaning
      fill_histograms(event, "Muon1_HighPt");
      lumihists_Muon1_HighPt->fill(event);
    }
    if( !(muon_is_high || muon_is_low) ) return false; // require ==1 muon in event
    if(debug) cout << "[ZprimeAnalysisModule] muonID cleaner: ok" << endl;
    fill_histograms(event, "Muon1_Tot");               // we have ==1 muon and ==0 electrons in event
  }


  // electrons 
  // (includes some isMuon statements for eleTrigger SF measurement in muon channel)
  double electron_pt_high(120.); // electron category pT threshold
  bool ele_is_low = false;       // boolean for low-pT electron category
  bool ele_is_high = false;      // boolean for high-pT electron category

  // Set electron category booleans based on pT
  // but first remove ECAL-gap electrons
  if(isElectron || (isMuon && isEleTriggerMeasurement)){
    vector<Electron>* electrons = event.electrons;
    for(unsigned int i=0; i<electrons->size(); i++){
      // transition region between the barrel and endcaps of ECAL defined as (1.44< eta <1.57)
      if(abs(event.electrons->at(i).eta()) > 1.44 && abs(event.electrons->at(i).eta()) < 1.57) return false;
      if(event.electrons->at(i).pt()<=electron_pt_high){
        ele_is_low = true;}
      else{ele_is_high = true;}
    }
  }
  sort_by_pt<Electron>(*event.electrons);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////// ElectronTrigger SF measurement DONE IN MUON CHANNEL ///////////////////////
  if(isMuon && isEleTriggerMeasurement){
    if(ele_is_low){
      if(!NEle1_selection->passes(event)) return false; // require ==1 electron
      electron_cleaner_low->process(event);             // clean low-pT electrons
      if(!NEle1_selection->passes(event)) return false; // require ==1 electron after cleaning
      fill_histograms(event, "1Mu1Ele_LowPt");
    }
    if(ele_is_high){
      if(!NEle1_selection->passes(event)) return false; // require ==1 electron
      electron_cleaner_high->process(event);            // clean high-pT electrons
      if(!NEle1_selection->passes(event)) return false; // require ==1 electron after cleaning
      fill_histograms(event, "1Mu1Ele_HighPt");
    }
    if( !(ele_is_high || ele_is_low) ) return false;    // require ==1 electron in event
    if(debug) cout << "[ZprimeAnalysisModule] electronID cleaner in EleTriggerMeasurement: ok" << endl;
    fill_histograms(event, "1Mu1Ele_Tot");              // we have ==1 electron and ==1 muon in event
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  // Selection: exactly 1 electron and 0 muons
  if(isElectron){
    if(!MuonVeto_selection->passes(event)) return false; // require ==0 muons
    if(ele_is_low){
      if(!NEle1_selection->passes(event)) return false; // require ==1 electron
      electron_cleaner_low->process(event);             // clean low-pT electrons
      if(!NEle1_selection->passes(event)) return false; // require ==1 electron after cleaning
      fill_histograms(event, "Ele1_LowPt");
      lumihists_Ele1_LowPt->fill(event);
    }
    if(ele_is_high){
      if(!NEle1_selection->passes(event)) return false; // require ==1 electron
      electron_cleaner_high->process(event);            // clean high-pT electrons
      if(!NEle1_selection->passes(event)) return false; // require ==1 electron after cleaning
      fill_histograms(event, "Ele1_HighPt");
      lumihists_Ele1_HighPt->fill(event);
    }
    if( !(ele_is_high || ele_is_low)) return false;     // require ==1 electron in event
    if(debug) cout << "[ZprimeAnalysisModule] electronID cleaner: ok" << endl;
    fill_histograms(event, "Ele1_Tot");                 // we have ==1 electron and ==0 muons in event
  }
  ////////////////////// End: Cleaning lepton collections with ID based on pT///////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////




  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// lepton ID_SF, ISO_SF (low-pT muon only), RECO_SF, Trigger Selection & SF, Nlepton Selection, and 2D Selection //////

  // electronID SF
  if(isMuon && !isEleTriggerMeasurement){ // muon channel withOUT eleTriggerSF measurement
    sf_ele_id_dummy->process(event);      // dummy eleID SF
  }
  if(isMuon && isEleTriggerMeasurement){  // muon channel with eleTriggerSF measurement 
    if(ele_is_low){sf_ele_id_low->process(event);}        // low-pT eleID SF
    else if(ele_is_high){sf_ele_id_high->process(event);} // high-pT eleID SF
    if(debug) cout << "[ZprimeAnalysisModule] electronID SF in EleTriggerMeasurement: ok" << endl;
  }
  if(isElectron){
    if(ele_is_low){sf_ele_id_low->process(event);}        // low-pT eleID SF
    else if(ele_is_high){sf_ele_id_high->process(event);} // high-pT eleID SF
    if(debug) cout << "[ZprimeAnalysisModule] electronID SF: ok" << endl;
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
    if(debug) cout << "[ZprimeAnalysisModule] muonISO SF: ok" << endl;
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
    if(debug) cout << "[ZprimeAnalysisModule] muonID SF: ok" << endl;
    fill_histograms(event, "IdMuon_SF");
  }
  if(isElectron){ // dummy muonID SFs for electron channel
    sf_muon_id_stat_dummy->process(event);
    sf_muon_id_syst_dummy->process(event);
  }

  // electronRECO SF
  if(isMuon && !isEleTriggerMeasurement){ // muon channel withOUT eleTriggerSF measurement
    sf_ele_reco_dummy->process(event);    // dummy electronRECO SF
  }
  if(isMuon && isEleTriggerMeasurement){  // muon channel with eleTriggerSF measurement
    sf_ele_reco->process(event);          // electronRECO SF
    if(debug) cout << "[ZprimeAnalysisModule] electronRECO SF in EleTriggerMeasurement: ok" << endl;
  }
  if(isElectron){                         // electron channel
    sf_ele_reco->process(event);          // electronRECO SF
    if(debug) cout << "[ZprimeAnalysisModule] electronRECO SF: ok" << endl;
    fill_histograms(event, "RecoEle_SF");
  }

  // muonRECO SF
  sf_muon_reco->process(event);
  fill_histograms(event, "MuonReco_SF");

  //////////////////////////////////////////////////////////////////// muonTrigger SELECTION ////////////////////////////////////////////////////////////////////
  if(isMuon){
    // low pt
    if(muon_is_low){ 
      if(isUL16preVFP || isUL16postVFP){ // UL16 pre&post eras
        if(!(Trigger_mu_A_selection->passes(event) || Trigger_mu_B_selection->passes(event))) return false;}
      else{ // UL17 and UL18 eras
        if(!Trigger_mu_A_selection->passes(event)) return false;}}

    // high pt
    if(muon_is_high){
      // UL16 pre&post eras
      if(isUL16preVFP || isUL16postVFP){
        if(!isMC){ // DATA RunB
          if( event.run < 274889){if(!Trigger_mu_C_selection->passes(event)) return false;} 
          else{ if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_D_selection->passes(event))) return false;}
        }
        else{ // MC RunB
          float runB_UL16_mu = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
          if( runB_UL16_mu < 0.1429){if(!Trigger_mu_C_selection->passes(event)) return false;}     
          else{ if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_D_selection->passes(event))) return false;}
        }
      }
      //UL17 DATA RunB
      if(isUL17 && !isMC){ 
        if(event.run <= 299329){ if(!Trigger_mu_C_selection->passes(event)) return false;}
        else{ if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_E_selection->passes(event) || Trigger_mu_F_selection->passes(event))) return false;}
      }
      // UL17 MC RunB
      if(isUL17 && isMC){ 
        float runB_mu = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if(runB_mu <= 0.1158){
          if( !Trigger_mu_C_selection->passes(event)) return false;}
          else{ if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_E_selection->passes(event) || Trigger_mu_F_selection->passes(event))) return false;
        }
      }
      // UL18
      if(isUL18){
        if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_E_selection->passes(event) || Trigger_mu_F_selection->passes(event))) return false;
      }
    }
    if(debug) cout << "[ZprimeAnalysisModule] muonTrigger Selection: passed" << endl;
    fill_histograms(event, "TriggerMuon");
    lumihists_TriggerMuon->fill(event);
  }

  ///////////////////////////////////////////////////////////////////// electronTrigger SELECTION /////////////////////////////////////////////////////////////////////
  if(isElectron){
    // low pt
    if(ele_is_low){
      if(isPhoton) return false;
      if(!Trigger_ele_A_selection->passes(event)) return false;
    }

    // high pt
    if(ele_is_high){
      // MC //
      // UL16 pre&post and UL18 eras
      if(isMC && (isUL16preVFP || isUL16postVFP || isUL18) ){ 
        if(!(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event) || Trigger_ph_A_selection->passes(event))) return false;
      }
      // UL17 era
      if(isMC && isUL17){ 
        float runB_ele = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if(runB_ele <= 0.1158){ // in RunB (below runnumb 299329) Ele115 does not exist, use Ele35 instead. To apply randomly in MC if random numb < RunB percetage (11.58%, calculated by Christopher Matthies)
          if(!(Trigger_ele_A_selection->passes(event) || Trigger_ph_A_selection->passes(event))) return false;}
          else{ if(!(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event) || Trigger_ph_A_selection->passes(event))) return false;}
      }
      // DATA //
      if(!isMC){
        // UL16 pre&post eras
        if(isUL16preVFP || isUL16postVFP){ 
          if(isPhoton){ // photon stream
            if(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event)) return false;
            if(!Trigger_ph_A_selection->passes(event)) return false;}
          else{ // electron stream
            if(!(Trigger_ele_A_selection->passes(event) ||  Trigger_ele_B_selection->passes(event) )) return false;
          }
        }
        // UL17 era
        if(isUL17){ 
          // below runnumb trigger Ele115 does not exist
          if(event.run <= 299329){ 
            if(isPhoton){ // photon stream
              if(Trigger_ele_A_selection->passes(event)) return false;
              if(!Trigger_ph_A_selection->passes(event)) return false;}
            else{ // electron stream
              if(!Trigger_ele_A_selection->passes(event)) return false;}
            }
          // above runnumb with Ele115
          else{ 
            if(isPhoton){ // photon stream
              if(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event)) return false;
              if(!Trigger_ph_A_selection->passes(event)) return false;}
            else{ // electron stream
              if(! (Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event))) return false;}
          }
        }
        // UL18 era
        if(isUL18){ 
          if(!(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event)|| Trigger_ph_A_selection->passes(event))) return false;
        }
      }
    }
    if(debug) cout << "[ZprimeAnalysisModule] electronTrigger Selection: passed" << endl;
    fill_histograms(event, "TriggerEle");
    lumihists_TriggerEle->fill(event);
  }

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
    if(debug) cout << "[ZprimeAnalysisModule] muonTrigger SF: ok" << endl;
    fill_histograms(event, "TriggerMuon_SF");
  }
  if(isElectron){ // dummy muonTrigger SFs for electron channel
    sf_muon_trigger_stat_dummy->process(event);
    sf_muon_trigger_syst_dummy->process(event);
  }


  /////////////////////////////////////////////////////////////// Nlepton selection ///////////////////////////////////////////////////////////////
  // when NOT doing eleTrigger SF measurement
  if(!isEleTriggerMeasurement){ 
    if((event.muons->size() + event.electrons->size()) != 1) return false; // require (==1electron & ==0muon) or (==0electron & ==1muon) in event
  }
  // when DOING eleTrigger SF measurement
  else{ 
    if( event.muons->size() != 1 && event.electrons->size() != 1) return false; // require (==1electron & ==1muon) in event
  }
  if(debug) cout << "[ZprimeAnalysisModule] N leptons ok: Nelectrons= " << event.electrons->size() << ", Nmuons= " << event.muons->size() << endl;




  //////////////////////////////// 2D-cut selection [dr(lep,jet) OR pTrel(lep,jet)] on high-pT leptons ONLY ///////////////////////////////////
  // muon channel (withOUT eleTrigger SF measurement) in high-pT category
  if(isMuon && !isEleTriggerMeasurement && muon_is_high){ 
    if(!TwoDCut_selection->passes(event)) return false; // dr >0.4 OR pTrel >25 GeV
    if(debug) cout << "[ZprimeAnalysisModule] 2D-cut on high-pT muon Selection: passed" << endl;
  }
  fill_histograms(event, "TwoDCut_Muon");
  lumihists_TwoDCut_Muon->fill(event);

  // electron channel in high-pT category
  if(isElectron && ele_is_high){
    if(!TwoDCut_selection->passes(event)) return false;
    if(debug) cout << "[ZprimeAnalysisModule] 2D-cut on high-pT electron Selection: passed" << endl;
  }
  fill_histograms(event, "TwoDCut_Ele");
  lumihists_TwoDCut_Ele->fill(event);

  // muon channel WITH eleTrigger SF measurement in both lepton high-pT categories
  if(isMuon && isEleTriggerMeasurement && (muon_is_high || ele_is_high)){
    if(!TwoDCut_selection->passes(event)) return false;
    if(debug) cout << "[ZprimeAnalysisModule] 2D-cut on high-pT lepton Selection in EleTriggerMeasurement: passed" << endl;
  }
  ///////////////////////////// 2D-cut selection is applied on low-pT muons at end of baseline selection //////////////////////////////////////




  if(debug) cout << "[ZprimeAnalysisModule] before matching" << endl;
  fill_histograms(event, "CHS_Before");

  ///////// Jet Matching: CHS to PUPPI /////////
  AK4PuppiCHS_matching->process(event);
  // if(!AK4PuppiCHS_matching->process(event)) return false; // use if we want to veto events where not matching CHS is found for any PUPPI jet
  if(debug) cout << "[ZprimeAnalysisModule] AK4PuppiCHS_matching: ok" << endl;
  fill_histograms(event, "CHS_After");
  h_CHSMatchHists->fill(event);
  


  /////// Jet Selection: >=1 AK4 PUPPI jets ///////
  // pT >50 GeV
  if(!Jet1_selection->passes(event)) return false;
  if(debug) cout << "[ZprimeAnalysisModule] Jet1 Selection: passed" << endl;
  fill_histograms(event, "Jet1");
  lumihists_Jet1->fill(event); // Fill lumihist after requiring >=1 AK4 PUPPI jet


  /////// Jet Selection: >=2 AK4 PUPPI jets ///////
  // pT >40 GeV for electron channel
  // pT >50 GeV for muon channel
  if(!Jet2_selection->passes(event)) return false;
  if(debug) cout << "[ZprimeAnalysisModule] Jet2 Selection: passed" << endl;
  fill_histograms(event, "Jet2");
  lumihists_Jet2->fill(event); // Fill lumihist after requiring >=2 AK4 PUPPI jets


  ////////////// MET Selection //////////////
  // >60 GeV for electron channel
  // >70 GeV for muon channel
  if(!met_sel->passes(event)) return false;
  if(debug) cout << "[ZprimeAnalysisModule] MET Selection: passed" << endl;
  fill_histograms(event, "MET");
  lumihists_MET->fill(event);  // Fill lumihist after MET cut


  // HT-lep Selection: >0 GeV //
  if(isMuon){
    if(!htlep_sel->passes(event)) return false;
    if(debug) cout << "[ZprimeAnalysisModule] HTlep Selection: passed" << endl;
    fill_histograms(event, "HTlep");
    lumihists_HTlep->fill(event); // Fill lumihist after HTlep cut
  }



  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// Begin b-tagging section ///////////////////////////////

  // Fill histograms before b-tagging SF - used to extract Custom BTag SF in (NJets,HT)
  if(debug) cout << "[ZprimeAnalysisModule] before b-tagging" << endl;
  fill_histograms(event, "BeforeBtagSF");
  h_CHSMatchHists_beforeBTagSF->fill(event);

  ////// Apply b-tag Shape SF //////
  sf_btagging->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] b-tag SF: ok" << endl;
  fill_histograms(event, "AfterBtagSF");
  h_CHSMatchHists_afterBTagSF->fill(event);

  // Apply custom SF to correct for b-tag SF shape effects on NJets/HT (HT := sum of jet pTs)
  if(isMC && isMuon){
    float custom_sf;
    vector<Jet>* jets = event.jets;
    int Njets = jets->size();
    double st_jets = 0.;
    for(const auto & jet : *jets) st_jets += jet.pt();
    custom_sf = ratio_hist_muon->GetBinContent( ratio_hist_muon->GetXaxis()->FindBin(Njets), ratio_hist_muon->GetYaxis()->FindBin(st_jets) );
    event.weight *= custom_sf;
    if(debug) cout << "[ZprimeAnalysisModule] custom SF in muon channel: ok" << endl;
  }
  if(isMC && !isMuon){
    float custom_sf;
    vector<Jet>* jets = event.jets;
    int Njets = jets->size();
    double st_jets = 0.;
    for(const auto & jet : *jets) st_jets += jet.pt();
    custom_sf = ratio_hist_ele->GetBinContent( ratio_hist_ele->GetXaxis()->FindBin(Njets), ratio_hist_ele->GetYaxis()->FindBin(st_jets) );
    event.weight *= custom_sf;
    if(debug) cout << "[ZprimeAnalysisModule] custom SF in electron channel: ok" << endl;
  }
  h_CHSMatchHists_after2DBTagSF->fill(event);
  fill_histograms(event, "AfterCustomBtagSF");

  // b-tag Selection: >= 1 b-tag medium WP 
  if(!AK4PuppiCHS_BTagging->passes(event)) return false;
  if(debug) cout << "[ZprimeAnalysisModule] b-tag1 Selection: passed" << endl;
  fill_histograms(event, "Btags1");
  h_CHSMatchHists_afterBTag->fill(event);
  /////////////////////////////// End b-tagging section ////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////




  // Higher Order Corrections: EWK & QCD NLO //
  NLOCorrections_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] NLO Correction: ok" << endl;
  fill_histograms(event, "NLOCorrections");

  // electronTrigger SF 
  if(!isEleTriggerMeasurement) sf_ele_trigger->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] electronTrigger SF: ok" << endl;
  fill_histograms(event, "TriggerEle_SF");

  // Veto events with >= 2 TopTagged large-R jets
  if(!TopTagVetoSelection->passes(event)) return false;
  if(debug) cout << "[ZprimeAnalysisModule] TopTag Veto Selection: passed" << endl;
  fill_histograms(event, "TopTagVeto");

  // Veto events with DeltaEta(j1, j2) > 3 to suppress multijet background
  if(!DeltaEta_selection->passes(event)){
    if(debug) cout << "[ZprimeAnalysisModule] DeltaEta(j1, j2) Selection: failed" << endl;
    fill_histograms(event, "DeltaEtaCut_FAIL");
    return false;
  }
  if(debug) cout << "[ZprimeAnalysisModule] DeltaEta(j1, j2) Selection: passed" << endl;
  fill_histograms(event, "DeltaEtaCut_PASS");            // all (high&low) muons before 2D cut

  // TwoD for low-pT muons (dr >0.3 OR pTrel >10 GeV)
  if(isMuon && muon_is_low){
    if(debug) cout << "[ZprimeAnalysisModule] 2D-cut on low-pT muons: before"<<endl;
    fill_histograms(event, "LowPtMuons_Before2Dcut"); // low-pT muons before 2D cut

    // Events that FAIL the low-pT muon 2D cut
    if(!TwoDCut_selection_low1->passes(event)){ // QCD enriched region
      if(debug) cout <<"[ZprimeAnalysisModule] 2D-cut on low-pT muons: failed"<<endl;   
      fill_histograms(event, "LowPtMuons_FAIL2Dcut"); // low-pT muons that fail
      return false;
    }
    if(debug) cout << "[ZprimeAnalysisModule] 2D-cut on low-pT muons: after"<<endl;
    fill_histograms(event, "LowPtMuons_PASS2Dcut");  // low-pT muons that pass
  }

  if(debug)cout <<"[ZprimeAnalysisModule] 2D-cut on low-pT muons: passed"<<endl;  
  fill_histograms(event, "EndOfBaselineSelection");     // all (high&low) muons that pass
  lumihists_TwoDCut_Muon_LowPt->fill(event);

  /////////////////////////////////////////////// End Baseline section ///////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// ttbar system reconstruction ///////////////////////////////

  // Build all possible ttbar candidates and extract discriminators
  CandidateBuilder->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] CandidateBuilder: ok" << endl;

  if(isMC) CorrectMatchDiscriminatorZprime->process(event);
  if(isMC && debug) cout << "[ZprimeAnalysisModule] CorrectMatchDiscriminator: ok" << endl;

  Chi2DiscriminatorZprime->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] Chi2Discriminator: ok" << endl;

  // Save Spin Correlation variables to TTree
  SpinCorrelations_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] SpinCorrelations_module: ok" << endl;



  //////////////////////// check matchable and correct match selections BEFORE Chi2 cut //////////////////////////
  if(isMC && TTbarMatchable_selection->passes(event)){                                                          //
    if(debug) cout << "[ZprimeAnalysisModule] TTbarMatchable Selection before chi2: passed" << endl;            //
    fill_histograms(event, "MatchableBeforeChi2Cut");}                                                          //
  if(isMC && CorrectMatchDiscriminatorZprime->process(event)){                                                  //
    if(debug) cout << "[ZprimeAnalysisModule] CorrectMatchDiscriminator selection before chi2: passed" << endl; //
    fill_histograms(event, "CorrectMatchBeforeChi2Cut");}                                                       //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Select events whose chi2 candidates have chi2 < 30
  if(!Chi2_selection->passes(event)){ // <---------------------------------------------------------------- chi2 cut
    if(debug) cout << "[ZprimeAnalysisModule] Chi2 Selection: failed" << endl;
    fill_histograms(event, "Chi2cut_FAIL");
    return false;
  }

  if(debug) cout << "[ZprimeAnalysisModule] Chi2 Selection: passed" << endl;
  fill_histograms(event, "Chi2cut_PASS");
  lumihists_Chi2->fill(event);

  //////////////////////// check matchable and correct match selections AFTER Passing Chi2 cut /////////////////
  if(isMC && TTbarMatchable_selection->passes(event)){                                                        //
    if(debug) cout << "[ZprimeAnalysisModule] TTbarMatchable Selection after chi2 selection: passed" << endl; //
    fill_histograms(event, "Matchable");}                                                                     //
  if(isMC && CorrectMatchDiscriminatorZprime->process(event)){                                                //
    if(debug) cout << "[ZprimeAnalysisModule] CorrectMatchDiscriminator after chi2 selection: ok" << endl;    //
    fill_histograms(event, "CorrectMatch");}                                                                  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// Define Event Topologies //////////////////////////////////

  if(ZprimeTopTag_selection->passes(event)){ // Merged
    if(debug) cout << "[ZprimeAnalysisModule] TopTag Selection: passed" << endl; 
    fill_histograms(event, "Merged");

    ////////////// check matchable and correct match selections AFTER Chi2 cut in Merged ////
    if(isMC && TTbarMatchable_selection->passes(event)){                                   //
      if(debug) cout << "[ZprimeAnalysisModule] TTbarMatchable Selection: passed" << endl; //
      fill_histograms(event, "Matchable_Merged");}                                         //
    if(isMC && CorrectMatchDiscriminatorZprime->process(event)){                           //
      if(debug) cout << "[ZprimeAnalysisModule] CorrectMatchDiscriminator: ok" << endl;    //
      fill_histograms(event, "CorrectMatch_Merged");}                                      //
    /////////////////////////////////////////////////////////////////////////////////////////
  }
  else{ // Resolved
    if(debug) cout << "[ZprimeAnalysisModule] TopTag Selection: failed" << endl; 
    fill_histograms(event, "Resolved");

    ////////////// check matchable and correct match selections AFTER Chi2 cut in Resolved //
    if(isMC && TTbarMatchable_selection->passes(event)){                                   //
      if(debug) cout << "[ZprimeAnalysisModule] TTbarMatchable Selection: passed" << endl; //
      fill_histograms(event, "Matchable_Resolved");}                                       //
    if(isMC && CorrectMatchDiscriminatorZprime->process(event)){                           //
      if(debug) cout << "[ZprimeAnalysisModule] CorrectMatchDiscriminator: ok" << endl;    //
      fill_histograms(event, "CorrectMatch_Resolved");}                                    //
    /////////////////////////////////////////////////////////////////////////////////////////
  }

  // Variables for DNN
  sort_by_pt<Jet>(*event.jets);
  Variables_module->process(event);
  if(debug) cout << "[ZprimeAnalysisModule] Variables_module: ok" << endl;

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule)