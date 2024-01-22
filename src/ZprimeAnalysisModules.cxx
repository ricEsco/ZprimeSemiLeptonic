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
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSystematicsModule.h>
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
  unique_ptr<AnalysisModule> sf_topmistag;
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
  // unique_ptr<ZprimeSemiLeptonicSystematicsModule> SystematicsModule;

  //Handles
  Event::Handle<bool> h_is_zprime_reconstructed_chi2, h_is_zprime_reconstructed_correctmatch;
  Event::Handle<float> h_chi2;   Event::Handle<float> h_weight;
  Event::Handle<float> h_MET;   Event::Handle<int> h_NPV;
  Event::Handle<float> h_lep1_pt; Event::Handle<float> h_lep1_eta;
  Event::Handle<float> h_ak4jet1_pt; Event::Handle<float> h_ak4jet1_eta;
  Event::Handle<float> h_ak8jet1_pt; Event::Handle<float> h_ak8jet1_eta;
  Event::Handle<float> h_Mttbar;

  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateCorrectMatch;

  // Reconstruction Optimization variables
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<TTbarGenProducer> ttgenprod;

  uhh2::Event::Handle< std::vector<Jet> > h_CHSjets_matched;  // Collection of CHS matched jets
  uhh2::Event::Handle< std::vector<TopJet> > h_DeepAK8TopTags;  // Collection of DeepAK8TopTagged jets

  Event::Handle<float> h_res_jet_bscore;
  Event::Handle<float> h_mer_subjet_bscore;
  Event::Handle<float> h_bscore_max;

  Event::Handle<float> h_res_jet_bscore_matchable;
  Event::Handle<float> h_mer_subjet_bscore_matchable;
  Event::Handle<float> h_bscore_max_matchable;

  Event::Handle<float> h_BestChi2;                            // Total Chi2 of "best" ttbar candidate
  Event::Handle<float> h_BestChi2_resolved;
  Event::Handle<float> h_BestChi2_merged;

  Event::Handle<float> h_BestChi2_passedMatching;             // Total Chi2 of "best" ttbar candidate that also passed matching requirements
  Event::Handle<float> h_BestChi2_passedMatching_resolved;
  Event::Handle<float> h_BestChi2_passedMatching_merged;

  Event::Handle<float> h_BestChi2_failedMatching;             // Total Chi2 of "best" ttbar candidate that did NOT matching requirements
  Event::Handle<float> h_BestChi2_failedMatching_resolved;
  Event::Handle<float> h_BestChi2_failedMatching_merged;

  Event::Handle<float> h_BestChi2hadronicleg;                 // Chi2 of hadronic leg of "best" ttbar candidate
  Event::Handle<float> h_BestChi2hadronicleg_resolved;
  Event::Handle<float> h_BestChi2hadronicleg_merged;

  Event::Handle<float> h_BestChi2hadronicleg_passedMatching;  // Chi2 of hadronic leg of "best" ttbar candidate that also passed matching requirements
  Event::Handle<float> h_BestChi2hadronicleg_passedMatching_resolved;
  Event::Handle<float> h_BestChi2hadronicleg_passedMatching_merged;

  Event::Handle<float> h_BestChi2hadronicleg_failedMatching;  // Chi2 of hadronic leg of "best" ttbar candidate that did NOT pass matching requirements
  Event::Handle<float> h_BestChi2hadronicleg_failedMatching_resolved;
  Event::Handle<float> h_BestChi2hadronicleg_failedMatching_merged;

  Event::Handle<float> h_BestChi2leptonicleg;                 // Chi2 of leptonic leg of "best" ttbar candidate
  Event::Handle<float> h_BestChi2leptonicleg_resolved;
  Event::Handle<float> h_BestChi2leptonicleg_merged;

  Event::Handle<float> h_BestChi2leptonicleg_passedMatching;  // Chi2 of leptonic leg of "best" ttbar candidate that also passed matching requirements
  Event::Handle<float> h_BestChi2leptonicleg_passedMatching_resolved;
  Event::Handle<float> h_BestChi2leptonicleg_passedMatching_merged;

  Event::Handle<float> h_BestChi2leptonicleg_failedMatching;  // Chi2 of leptonic leg of "best" ttbar candidate that did NOT pass matching requirements
  Event::Handle<float> h_BestChi2leptonicleg_failedMatching_resolved;
  Event::Handle<float> h_BestChi2leptonicleg_failedMatching_merged;

  Event::Handle<float> h_DrHadB;                           // DeltaR between hadronic b's (gen vs reco)
  Event::Handle<float> h_DrHadB_resolved;
  Event::Handle<float> h_DrHadB_merged;

  Event::Handle<float> h_DrHadB_passedChi2cut;             // DeltaR between hadronic b's (gen vs reco) that ALSO passed matching requirements
  Event::Handle<float> h_DrHadB_passedChi2cut_resolved;
  Event::Handle<float> h_DrHadB_passedChi2cut_merged;

  Event::Handle<float> h_DrHadB_failedChi2cut;             // DeltaR between hadronic b's (gen vs reco) that did NOT pass matching requirements
  Event::Handle<float> h_DrHadB_failedChi2cut_resolved;
  Event::Handle<float> h_DrHadB_failedChi2cut_merged;

  Event::Handle<float> h_DrLep;                          // DeltaR between leptons (gen vs reco)
  Event::Handle<float> h_DrLep_resolved;
  Event::Handle<float> h_DrLep_merged;

  Event::Handle<float> h_DrLep_passedChi2cut;            // DeltaR between leptons (gen vs reco) that ALSO passed matching requirements
  Event::Handle<float> h_DrLep_passedChi2cut_resolved;
  Event::Handle<float> h_DrLep_passedChi2cut_merged;

  Event::Handle<float> h_DrLep_failedChi2cut;            // DeltaR between leptons (gen vs reco) that did NOT pass matching requirements
  Event::Handle<float> h_DrLep_failedChi2cut_resolved;
  Event::Handle<float> h_DrLep_failedChi2cut_merged;

  Event::Handle<float> h_correctDr;                          // correct_dr of matcehd ttbar candidates
  Event::Handle<float> h_correctDr_resolved;
  Event::Handle<float> h_correctDr_merged;

  Event::Handle<float> h_correctDr_passedChi2cut;            // correct_dr of ttbar candidates that ALSO passed matching requirements
  Event::Handle<float> h_correctDr_passedChi2cut_resolved;
  Event::Handle<float> h_correctDr_passedChi2cut_merged;

  Event::Handle<float> h_correctDr_failedChi2cut;            // correct_dr of ttbar candidates that did NOT pass matching requirements
  Event::Handle<float> h_correctDr_failedChi2cut_resolved;
  Event::Handle<float> h_correctDr_failedChi2cut_merged;

  Event::Handle<float> h_DrHadB_matchable;                           // DeltaR between hadronic b's (gen vs reco) that are matchable
  Event::Handle<float> h_DrHadB_matchable_resolved;
  Event::Handle<float> h_DrHadB_matchable_merged;

  Event::Handle<float> h_DrHadB_matchable_passedChi2cut;             // DeltaR between hadronic b's (gen vs reco) that are matchable ALSO passed matching requirements
  Event::Handle<float> h_DrHadB_matchable_passedChi2cut_resolved;
  Event::Handle<float> h_DrHadB_matchable_passedChi2cut_merged;

  Event::Handle<float> h_DrHadB_matchable_failedChi2cut;             // DeltaR between hadronic b's (gen vs reco) that are matchable did NOT pass matching requirements
  Event::Handle<float> h_DrHadB_matchable_failedChi2cut_resolved;
  Event::Handle<float> h_DrHadB_matchable_failedChi2cut_merged;

  Event::Handle<float> h_DrLep_matchable;                          // DeltaR between leptons (gen vs reco) are matchable
  Event::Handle<float> h_DrLep_matchable_resolved;
  Event::Handle<float> h_DrLep_matchable_merged;

  Event::Handle<float> h_DrLep_matchable_passedChi2cut;            // DeltaR between leptons (gen vs reco) that are matchable ALSO passed matching requirements
  Event::Handle<float> h_DrLep_matchable_passedChi2cut_resolved;
  Event::Handle<float> h_DrLep_matchable_passedChi2cut_merged;

  Event::Handle<float> h_DrLep_matchable_failedChi2cut;            // DeltaR between leptons (gen vs reco) that are matchable did NOT pass matching requirements
  Event::Handle<float> h_DrLep_matchable_failedChi2cut_resolved;
  Event::Handle<float> h_DrLep_matchable_failedChi2cut_merged;


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
  sf_topmistag.reset(new TopMistagScaleFactor(ctx));
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

  // if(!isEleTriggerMeasurement) SystematicsModule.reset(new ZprimeSemiLeptonicSystematicsModule(ctx));


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
  h_BestZprimeCandidateCorrectMatch = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestCorrectMatch");
  h_chi2 = ctx.declare_event_output<float> ("rec_chi2");
  h_MET = ctx.declare_event_output<float> ("met_pt");
  h_Mttbar = ctx.declare_event_output<float> ("Mttbar");
  h_lep1_pt = ctx.declare_event_output<float> ("lep1_pt");
  h_lep1_eta = ctx.declare_event_output<float> ("lep1_eta");
  h_ak4jet1_pt = ctx.declare_event_output<float> ("ak4jet1_pt");
  h_ak4jet1_eta = ctx.declare_event_output<float> ("ak4jet1_eta");
  h_ak8jet1_pt = ctx.declare_event_output<float> ("ak8jet1_pt");
  h_ak8jet1_eta = ctx.declare_event_output<float> ("ak8jet1_eta");

  // Reconstruction Optimization variables
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");                         // Access to gen-level particles

  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");       // Collection of CHS matched jets
  h_DeepAK8TopTags = ctx.get_handle< std::vector<TopJet>>("DeepAK8TopTags"); // Collection of DeepAK8TopTagged jets

  h_res_jet_bscore = ctx.declare_event_output<float> ("res_jet_bscore");
  h_mer_subjet_bscore = ctx.declare_event_output<float> ("mer_subjet_bscore");
  h_bscore_max = ctx.declare_event_output<float> ("bscore_max");

  h_res_jet_bscore_matchable = ctx.declare_event_output<float> ("res_jet_bscore_matchable");
  h_mer_subjet_bscore_matchable = ctx.declare_event_output<float> ("mer_subjet_bscore_matchable");
  h_bscore_max_matchable = ctx.declare_event_output<float> ("bscore_max_matchable");

  h_BestChi2 = ctx.declare_event_output<float> ("BestChi2");                               // Total Chi2 of "best" ttbar candidate
  h_BestChi2_resolved = ctx.declare_event_output<float> ("BestChi2_resolved");
  h_BestChi2_merged = ctx.declare_event_output<float> ("BestChi2_merged");

  h_BestChi2_passedMatching = ctx.declare_event_output<float> ("BestChi2_passedMatching"); // Total Chi2 of "best" ttbar candidate that ALSO passed matching requirements
  h_BestChi2_passedMatching_resolved = ctx.declare_event_output<float> ("BestChi2_passedMatching_resolved");
  h_BestChi2_passedMatching_merged = ctx.declare_event_output<float> ("BestChi2_passedMatching_merged");

  h_BestChi2_failedMatching = ctx.declare_event_output<float> ("BestChi2_failedMatching"); // Total Chi2 of "best" ttbar candidate that did NOT pass matching requirements
  h_BestChi2_failedMatching_resolved = ctx.declare_event_output<float> ("BestChi2_failedMatching_resolved");
  h_BestChi2_failedMatching_merged = ctx.declare_event_output<float> ("BestChi2_failedMatching_merged");

  h_BestChi2hadronicleg = ctx.declare_event_output<float> ("BestChi2hadronicleg");                               // Chi2 of hadronic leg of "best" ttbar candidate
  h_BestChi2hadronicleg_resolved = ctx.declare_event_output<float> ("BestChi2hadronicleg_resolved");
  h_BestChi2hadronicleg_merged = ctx.declare_event_output<float> ("BestChi2hadronicleg_merged");

  h_BestChi2hadronicleg_passedMatching = ctx.declare_event_output<float> ("BestChi2hadronicleg_passedMatching"); // Chi2 of hadronic leg of "best" ttbar candidate that ALSO passed matching requirements
  h_BestChi2hadronicleg_passedMatching_resolved = ctx.declare_event_output<float> ("BestChi2hadronicleg_passedMatching_resolved");
  h_BestChi2hadronicleg_passedMatching_merged = ctx.declare_event_output<float> ("BestChi2hadronicleg_passedMatching_merged");

  h_BestChi2hadronicleg_failedMatching = ctx.declare_event_output<float> ("BestChi2hadronicleg_failedMatching"); // Chi2 of hadronic leg of "best" ttbar candidate that did NOT pass matching requirements
  h_BestChi2hadronicleg_failedMatching_resolved = ctx.declare_event_output<float> ("BestChi2hadronicleg_failedMatching_resolved");
  h_BestChi2hadronicleg_failedMatching_merged = ctx.declare_event_output<float> ("BestChi2hadronicleg_failedMatching_merged");

  h_BestChi2leptonicleg = ctx.declare_event_output<float> ("BestChi2leptonicleg");                              // Chi2 of leptonic leg of "best" ttbar candidate
  h_BestChi2leptonicleg_resolved = ctx.declare_event_output<float> ("BestChi2leptonicleg_resolved");
  h_BestChi2leptonicleg_merged = ctx.declare_event_output<float> ("BestChi2leptonicleg_merged");

  h_BestChi2leptonicleg_passedMatching = ctx.declare_event_output<float> ("BestChi2leptonicleg_passedMatching"); // Chi2 of leptonic leg of "best" ttbar candidate that ALSO passed matching requirements
  h_BestChi2leptonicleg_passedMatching_resolved = ctx.declare_event_output<float> ("BestChi2leptonicleg_passedMatching_resolved");
  h_BestChi2leptonicleg_passedMatching_merged = ctx.declare_event_output<float> ("BestChi2leptonicleg_passedMatching_merged");

  h_BestChi2leptonicleg_failedMatching = ctx.declare_event_output<float> ("BestChi2leptonicleg_failedMatching"); // Chi2 of leptonic leg of "best" ttbar candidate that did NOT pass matching requirements
  h_BestChi2leptonicleg_failedMatching_resolved = ctx.declare_event_output<float> ("BestChi2leptonicleg_failedMatching_resolved");
  h_BestChi2leptonicleg_failedMatching_merged = ctx.declare_event_output<float> ("BestChi2leptonicleg_failedMatching_merged");

  h_DrHadB = ctx.declare_event_output<float> ("DrHadB");                              // DeltaR between hadronic b's (gen vs reco)
  h_DrHadB_resolved = ctx.declare_event_output<float> ("DrHadB_resolved");
  h_DrHadB_merged = ctx.declare_event_output<float> ("DrHadB_merged");

  h_DrHadB_passedChi2cut = ctx.declare_event_output<float> ("DrHadB_passedChi2cut");  // DeltaR between hadronic b's (gen vs reco) that ALSO passed the chi2 cut
  h_DrHadB_passedChi2cut_resolved = ctx.declare_event_output<float> ("DrHadB_passedChi2cut_resolved");
  h_DrHadB_passedChi2cut_merged = ctx.declare_event_output<float> ("DrHadB_passedChi2cut_merged");

  h_DrHadB_failedChi2cut = ctx.declare_event_output<float> ("DrHadB_failedChi2cut");  // DeltaR between hadronic b's (gen vs reco) that did NOT pass the chi2 cut
  h_DrHadB_failedChi2cut_resolved = ctx.declare_event_output<float> ("DrHadB_failedChi2cut_resolved");
  h_DrHadB_failedChi2cut_merged = ctx.declare_event_output<float> ("DrHadB_failedChi2cut_merged");

  h_DrLep = ctx.declare_event_output<float> ("DrLep");                             // DeltaR between lepton's (gen vs reco)
  h_DrLep_resolved = ctx.declare_event_output<float> ("DrLep_resolved");
  h_DrLep_merged = ctx.declare_event_output<float> ("DrLep_merged");

  h_DrLep_passedChi2cut = ctx.declare_event_output<float> ("DrLep_passedChi2cut"); // DeltaR between lepton's (gen vs reco) that ALSO passed the chi2 cut
  h_DrLep_passedChi2cut_resolved = ctx.declare_event_output<float> ("DrLep_passedChi2cut_resolved");
  h_DrLep_passedChi2cut_merged = ctx.declare_event_output<float> ("DrLep_passedChi2cut_merged");

  h_DrLep_failedChi2cut = ctx.declare_event_output<float> ("DrLep_failedChi2cut"); // DeltaR between lepton's (gen vs reco) that did NOT pass the chi2 cut
  h_DrLep_failedChi2cut_resolved = ctx.declare_event_output<float> ("DrLep_failedChi2cut_resolved");
  h_DrLep_failedChi2cut_merged = ctx.declare_event_output<float> ("DrLep_failedChi2cut_merged");

  h_correctDr = ctx.declare_event_output<float> ("correctDr");                             // correct_dr of ttbar candidates that pass matching requirements
  h_correctDr_resolved = ctx.declare_event_output<float> ("correctDr_resolved");
  h_correctDr_merged = ctx.declare_event_output<float> ("correctDr_merged");

  h_correctDr_passedChi2cut = ctx.declare_event_output<float> ("correctDr_passedChi2cut"); // correct_dr of ttbar candidates that ALSO passed the chi2 cut
  h_correctDr_passedChi2cut_resolved = ctx.declare_event_output<float> ("correctDr_passedChi2cut_resolved");
  h_correctDr_passedChi2cut_merged = ctx.declare_event_output<float> ("correctDr_passedChi2cut_merged");

  h_correctDr_failedChi2cut = ctx.declare_event_output<float> ("correctDr_failedChi2cut"); // correct_dr of ttbar candidates that did NOT pass the chi2 cut
  h_correctDr_failedChi2cut_resolved = ctx.declare_event_output<float> ("correctDr_failedChi2cut_resolved");
  h_correctDr_failedChi2cut_merged = ctx.declare_event_output<float> ("correctDr_failedChi2cut_merged");

  h_DrHadB_matchable = ctx.declare_event_output<float> ("DrHadB_matchable");                             // DeltaR between hadronic b's (gen vs reco) that are matchable
  h_DrHadB_matchable_resolved = ctx.declare_event_output<float> ("DrHadB_matchable_resolved");
  h_DrHadB_matchable_merged = ctx.declare_event_output<float> ("DrHadB_matchable_merged");

  h_DrHadB_matchable_passedChi2cut = ctx.declare_event_output<float> ("DrHadB_matchable_passedChi2cut"); // DeltaR between hadronic b's (gen vs reco) that are matchable ALSO passed the chi2 cut
  h_DrHadB_matchable_passedChi2cut_resolved = ctx.declare_event_output<float> ("DrHadB_matchable_passedChi2cut_resolved");
  h_DrHadB_matchable_passedChi2cut_merged = ctx.declare_event_output<float> ("DrHadB_matchable_passedChi2cut_merged");

  h_DrHadB_matchable_failedChi2cut = ctx.declare_event_output<float> ("DrHadB_matchable_failedChi2cut"); // DeltaR between hadronic b's (gen vs reco) that are matchable did NOT pass the chi2 cut
  h_DrHadB_matchable_failedChi2cut_resolved = ctx.declare_event_output<float> ("DrHadB_matchable_failedChi2cut_resolved");
  h_DrHadB_matchable_failedChi2cut_merged = ctx.declare_event_output<float> ("DrHadB_matchable_failedChi2cut_merged");

  h_DrLep_matchable = ctx.declare_event_output<float> ("DrLep_matchable");                             // DeltaR between lepton's (gen vs reco) are matchable
  h_DrLep_matchable_resolved = ctx.declare_event_output<float> ("DrLep_matchable_resolved");
  h_DrLep_matchable_merged = ctx.declare_event_output<float> ("DrLep_matchable_merged");

  h_DrLep_matchable_passedChi2cut = ctx.declare_event_output<float> ("DrLep_matchable_passedChi2cut"); // DeltaR between lepton's (gen vs reco) that are matchable ALSO passed the chi2 cut
  h_DrLep_matchable_passedChi2cut_resolved = ctx.declare_event_output<float> ("DrLep_matchable_passedChi2cut_resolved");
  h_DrLep_matchable_passedChi2cut_merged = ctx.declare_event_output<float> ("DrLep_matchable_passedChi2cut_merged");

  h_DrLep_matchable_failedChi2cut = ctx.declare_event_output<float> ("DrLep_matchable_failedChi2cut");  // DeltaR between lepton's (gen vs reco) that are matchable did NOT pass the chi2 cut
  h_DrLep_matchable_failedChi2cut_resolved = ctx.declare_event_output<float> ("DrLep_matchable_failedChi2cut_resolved");
  h_DrLep_matchable_failedChi2cut_merged = ctx.declare_event_output<float> ("DrLep_matchable_failedChi2cut_merged");


  h_NPV = ctx.declare_event_output<int> ("NPV");
  h_weight = ctx.declare_event_output<float> ("weight");

  sel_1btag.reset(new NJetSelection(1, -1, id_btag));
  sel_2btag.reset(new NJetSelection(2,-1, id_btag));

  // PUPPI CHS match modules & hists
  AK4PuppiCHS_matching.reset(new PuppiCHS_matching(ctx)); // match AK4 PUPPI jets to AK4 CHS jets for b-tagging
  AK4PuppiCHS_BTagging.reset(new PuppiCHS_BTagging(ctx)); // b-tagging on matched CHS jets
  h_CHSMatchHists.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch"));
  h_CHSMatchHists_beforeBTagSF.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_beforeBTagSF"));
  h_CHSMatchHists_afterBTagSF.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_afterBTagSF"));
  h_CHSMatchHists_after2DBTagSF.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_after2DBTagSF"));
  h_CHSMatchHists_afterBTag.reset(new ZprimeSemiLeptonicCHSMatchHists(ctx, "CHSMatch_afterBTag"));

  // Book histograms
  vector<string> histogram_tags = {"Weights_Init", "Weights_HEM", "Weights_PU", "Weights_Lumi", "Weights_TopPt", "Weights_MCScale", "Weights_Prefiring", "Weights_PS", "Weights_TopTag_SF", "Weights_TopMistag_SF", "Muon1_LowPt", "Muon1_HighPt", "Muon1_Tot", "Ele1_LowPt", "Ele1_HighPt", "Ele1_Tot", "1Mu1Ele_LowPt", "1Mu1Ele_HighPt", "1Mu1Ele_Tot", "IdMuon_SF", "IdEle_SF", "IsoMuon_SF", "RecoEle_SF", "MuonReco_SF", "TriggerMuon", "TriggerEle", "TriggerMuon_SF", "TwoDCut_Muon", "TwoDCut_Ele", "Jet1", "Jet2", "MET", "HTlep", "BeforeBtagSF", "AfterBtagSF", "AfterCustomBtagSF", "Btags1", "NLOCorrections", "TriggerEle_SF", "TTbarCandidate", "CorrectMatchDiscriminator", "Chi2Discriminator", "NNInputsBeforeReweight"};
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
    vector<TString> names = {"ST", "WJets", "DY", "QCD", "ALP_ttbar_signal", "ALP_ttbar_interference", "HscalarToTTTo1L1Nu2J_m365_w36p5_res", "HscalarToTTTo1L1Nu2J_m400_w40p0_res", "HscalarToTTTo1L1Nu2J_m500_w50p0_res", "HscalarToTTTo1L1Nu2J_m600_w60p0_res", "HscalarToTTTo1L1Nu2J_m800_w80p0_res", "HscalarToTTTo1L1Nu2J_m1000_w100p0_res", "HscalarToTTTo1L1Nu2J_m365_w36p5_int_pos", "HscalarToTTTo1L1Nu2J_m400_w40p0_int_pos", "HscalarToTTTo1L1Nu2J_m500_w50p0_int_pos", "HscalarToTTTo1L1Nu2J_m600_w60p0_int_pos", "HscalarToTTTo1L1Nu2J_m800_w80p0_int_pos", "HscalarToTTTo1L1Nu2J_m1000_w100p0_int_pos", "HscalarToTTTo1L1Nu2J_m365_w36p5_int_neg", "HscalarToTTTo1L1Nu2J_m400_w40p0_int_neg", "HscalarToTTTo1L1Nu2J_m500_w50p0_int_neg", "HscalarToTTTo1L1Nu2J_m600_w60p0_int_neg", "HscalarToTTTo1L1Nu2J_m800_w80p0_int_neg", "HscalarToTTTo1L1Nu2J_m1000_w100p0_int_neg", "HpseudoToTTTo1L1Nu2J_m365_w36p5_res", "HpseudoToTTTo1L1Nu2J_m400_w40p0_res", "HpseudoToTTTo1L1Nu2J_m500_w50p0_res", "HpseudoToTTTo1L1Nu2J_m600_w60p0_res", "HpseudoToTTTo1L1Nu2J_m800_w80p0_res", "HpseudoToTTTo1L1Nu2J_m1000_w100p0_res", "HpseudoToTTTo1L1Nu2J_m365_w36p5_int_pos", "HpseudoToTTTo1L1Nu2J_m400_w40p0_int_pos", "HpseudoToTTTo1L1Nu2J_m500_w50p0_int_pos", "HpseudoToTTTo1L1Nu2J_m600_w60p0_int_pos", "HpseudoToTTTo1L1Nu2J_m800_w80p0_int_pos", "HpseudoToTTTo1L1Nu2J_m1000_w100p0_int_pos", "HpseudoToTTTo1L1Nu2J_m365_w36p5_int_neg", "HpseudoToTTTo1L1Nu2J_m400_w40p0_int_neg", "HpseudoToTTTo1L1Nu2J_m500_w50p0_int_neg", "HpseudoToTTTo1L1Nu2J_m600_w60p0_int_neg", "HpseudoToTTTo1L1Nu2J_m800_w80p0_int_neg", "HpseudoToTTTo1L1Nu2J_m1000_w100p0_int_neg", "HscalarToTTTo1L1Nu2J_m365_w91p25_res", "HscalarToTTTo1L1Nu2J_m400_w100p0_res", "HscalarToTTTo1L1Nu2J_m500_w125p0_res", "HscalarToTTTo1L1Nu2J_m600_w150p0_res", "HscalarToTTTo1L1Nu2J_m800_w200p0_res", "HscalarToTTTo1L1Nu2J_m1000_w250p0_res", "HscalarToTTTo1L1Nu2J_m365_w91p25_int_pos", "HscalarToTTTo1L1Nu2J_m400_w100p0_int_pos", "HscalarToTTTo1L1Nu2J_m500_w125p0_int_pos", "HscalarToTTTo1L1Nu2J_m600_w150p0_int_pos", "HscalarToTTTo1L1Nu2J_m800_w200p0_int_pos", "HscalarToTTTo1L1Nu2J_m1000_w250p0_int_pos", "HscalarToTTTo1L1Nu2J_m365_w91p25_int_neg", "HscalarToTTTo1L1Nu2J_m400_w100p0_int_neg", "HscalarToTTTo1L1Nu2J_m500_w125p0_int_neg", "HscalarToTTTo1L1Nu2J_m600_w150p0_int_neg", "HscalarToTTTo1L1Nu2J_m800_w200p0_int_neg", "HscalarToTTTo1L1Nu2J_m1000_w250p0_int_neg", "HpseudoToTTTo1L1Nu2J_m365_w91p25_res", "HpseudoToTTTo1L1Nu2J_m400_w100p0_res", "HpseudoToTTTo1L1Nu2J_m500_w125p0_res", "HpseudoToTTTo1L1Nu2J_m600_w150p0_res", "HpseudoToTTTo1L1Nu2J_m800_w200p0_res", "HpseudoToTTTo1L1Nu2J_m1000_w250p0_res", "HpseudoToTTTo1L1Nu2J_m365_w91p25_int_pos", "HpseudoToTTTo1L1Nu2J_m400_w100p0_int_pos", "HpseudoToTTTo1L1Nu2J_m500_w125p0_int_pos", "HpseudoToTTTo1L1Nu2J_m600_w150p0_int_pos", "HpseudoToTTTo1L1Nu2J_m800_w200p0_int_pos", "HpseudoToTTTo1L1Nu2J_m1000_w250p0_int_pos", "HpseudoToTTTo1L1Nu2J_m365_w91p25_int_neg", "HpseudoToTTTo1L1Nu2J_m400_w100p0_int_neg", "HpseudoToTTTo1L1Nu2J_m500_w125p0_int_neg", "HpseudoToTTTo1L1Nu2J_m600_w150p0_int_neg", "HpseudoToTTTo1L1Nu2J_m800_w200p0_int_neg", "HpseudoToTTTo1L1Nu2J_m1000_w250p0_int_neg", "HscalarToTTTo1L1Nu2J_m365_w9p125_res", "HscalarToTTTo1L1Nu2J_m400_w10p0_res", "HscalarToTTTo1L1Nu2J_m500_w12p5_res", "HscalarToTTTo1L1Nu2J_m600_w15p0_res", "HscalarToTTTo1L1Nu2J_m800_w20p0_res", "HscalarToTTTo1L1Nu2J_m1000_w25p0_res", "HscalarToTTTo1L1Nu2J_m365_w9p125_int_pos", "HscalarToTTTo1L1Nu2J_m400_w10p0_int_pos", "HscalarToTTTo1L1Nu2J_m500_w12p5_int_pos", "HscalarToTTTo1L1Nu2J_m600_w15p0_int_pos", "HscalarToTTTo1L1Nu2J_m800_w20p0_int_pos", "HscalarToTTTo1L1Nu2J_m1000_w25p0_int_pos", "HscalarToTTTo1L1Nu2J_m365_w9p125_int_neg", "HscalarToTTTo1L1Nu2J_m400_w10p0_int_neg", "HscalarToTTTo1L1Nu2J_m500_w12p5_int_neg", "HscalarToTTTo1L1Nu2J_m600_w15p0_int_neg", "HscalarToTTTo1L1Nu2J_m800_w20p0_int_neg", "HscalarToTTTo1L1Nu2J_m1000_w25p0_int_neg", "HpseudoToTTTo1L1Nu2J_m365_w9p125_res", "HpseudoToTTTo1L1Nu2J_m400_w10p0_res", "HpseudoToTTTo1L1Nu2J_m500_w12p5_res", "HpseudoToTTTo1L1Nu2J_m600_w15p0_res", "HpseudoToTTTo1L1Nu2J_m800_w20p0_res", "HpseudoToTTTo1L1Nu2J_m1000_w25p0_res", "HpseudoToTTTo1L1Nu2J_m365_w9p125_int_pos", "HpseudoToTTTo1L1Nu2J_m400_w10p0_int_pos", "HpseudoToTTTo1L1Nu2J_m500_w12p5_int_pos", "HpseudoToTTTo1L1Nu2J_m600_w15p0_int_pos", "HpseudoToTTTo1L1Nu2J_m800_w20p0_int_pos", "HpseudoToTTTo1L1Nu2J_m1000_w25p0_int_pos", "HpseudoToTTTo1L1Nu2J_m365_w9p125_int_neg", "HpseudoToTTTo1L1Nu2J_m400_w10p0_int_neg", "HpseudoToTTTo1L1Nu2J_m500_w12p5_int_neg", "HpseudoToTTTo1L1Nu2J_m600_w15p0_int_neg", "HpseudoToTTTo1L1Nu2J_m800_w20p0_int_neg", "HpseudoToTTTo1L1Nu2J_m1000_w25p0_int_neg", "RSGluonToTT_M-500", "RSGluonToTT_M-1000", "RSGluonToTT_M-1500", "RSGluonToTT_M-2000", "RSGluonToTT_M-2500", "RSGluonToTT_M-3000", "RSGluonToTT_M-3500", "RSGluonToTT_M-4000", "RSGluonToTT_M-4500", "RSGluonToTT_M-5000", "RSGluonToTT_M-5500", "RSGluonToTT_M-6000", "ZPrimeToTT_M400_W40", "ZPrimeToTT_M500_W50", "ZPrimeToTT_M600_W60", "ZPrimeToTT_M700_W70", "ZPrimeToTT_M800_W80", "ZPrimeToTT_M900_W90", "ZPrimeToTT_M1000_W100", "ZPrimeToTT_M1200_W120", "ZPrimeToTT_M1400_W140", "ZPrimeToTT_M1600_W160", "ZPrimeToTT_M1800_W180", "ZPrimeToTT_M2000_W200", "ZPrimeToTT_M2500_W250", "ZPrimeToTT_M3000_W300", "ZPrimeToTT_M3500_W350", "ZPrimeToTT_M4000_W400", "ZPrimeToTT_M4500_W450", "ZPrimeToTT_M5000_W500", "ZPrimeToTT_M6000_W600",  "ZPrimeToTT_M7000_W700", "ZPrimeToTT_M8000_W800", "ZPrimeToTT_M9000_W900", "ZPrimeToTT_M400_W120", "ZPrimeToTT_M500_W150", "ZPrimeToTT_M600_W180", "ZPrimeToTT_M700_W210", "ZPrimeToTT_M800_W240", "ZPrimeToTT_M900_W270", "ZPrimeToTT_M1000_W300", "ZPrimeToTT_M1200_W360", "ZPrimeToTT_M1400_W420", "ZPrimeToTT_M1600_W480", "ZPrimeToTT_M1800_W540", "ZPrimeToTT_M2000_W600", "ZPrimeToTT_M2500_W750", "ZPrimeToTT_M3000_W900", "ZPrimeToTT_M3500_W1050", "ZPrimeToTT_M4000_W1200", "ZPrimeToTT_M4500_W1350", "ZPrimeToTT_M5000_W1500", "ZPrimeToTT_M6000_W1800", "ZPrimeToTT_M7000_W2100", "ZPrimeToTT_M8000_W2400", "ZPrimeToTT_M9000_W2700", "ZPrimeToTT_M400_W4", "ZPrimeToTT_M500_W5", "ZPrimeToTT_M600_W6", "ZPrimeToTT_M700_W7", "ZPrimeToTT_M800_W8", "ZPrimeToTT_M900_W9", "ZPrimeToTT_M1000_W10", "ZPrimeToTT_M1200_W12", "ZPrimeToTT_M1400_W14", "ZPrimeToTT_M1600_W16", "ZPrimeToTT_M1800_W18", "ZPrimeToTT_M2000_W20", "ZPrimeToTT_M2500_W25", "ZPrimeToTT_M3000_W30", "ZPrimeToTT_M3500_W35", "ZPrimeToTT_M4000_W40", "ZPrimeToTT_M4500_W45", "ZPrimeToTT_M5000_W50", "ZPrimeToTT_M6000_W60", "ZPrimeToTT_M7000_W70", "ZPrimeToTT_M8000_W80", "ZPrimeToTT_M9000_W90", "ZprimeDMToTTbarResoIncl_MZp1000_Mchi10_V1", "ZprimeDMToTTbarResoIncl_MZp1500_Mchi10_V1", "ZprimeDMToTTbarResoIncl_MZp2000_Mchi10_V1", "ZprimeDMToTTbarResoIncl_MZp2500_Mchi10_V1", "ZprimeDMToTTbarResoIncl_MZp3000_Mchi10_V1", "ZprimeDMToTTbarResoIncl_MZp3500_Mchi10_V1", "ZprimeDMToTTbarResoIncl_MZp4000_Mchi10_V1", "ZprimeDMToTTbarResoIncl_MZp4500_Mchi10_V1", "ZprimeDMToTTbarResoIncl_MZp5000_Mchi10_V1", "ZprimeDMToTTbarResoIncl_MZp2500_Mchi1000_A1", "ZprimeDMToTTbarResoIncl_MZp2500_Mchi1000_V1", "ZprimeDMToTTbarResoIncl_MZp2500_Mchi10_A1"};

    for(unsigned int i=0; i<names.size(); i++){
      if( ctx.get("dataset_version").find(names.at(i)) != std::string::npos ) sample_name = names.at(i);
    }
    if( (ctx.get("dataset_version").find("TTToHadronic") != std::string::npos) || (ctx.get("dataset_version").find("TTToSemiLeptonic") != std::string::npos) || (ctx.get("dataset_version").find("TTTo2L2Nu") != std::string::npos) ) sample_name = "TTbar";
    if( (ctx.get("dataset_version").find("WW") != std::string::npos) || (ctx.get("dataset_version").find("ZZ") != std::string::npos) || (ctx.get("dataset_version").find("WZ") != std::string::npos) ) sample_name = "Diboson";

    if(isMuon){
      TFile* f_btag2Dsf_muon = new TFile("/nfs/dust/cms/user/deleokse/RunII_106_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_muon_"+year+".root");
      ratio_hist_muon = (TH2F*)f_btag2Dsf_muon->Get("N_Jets_vs_HT_" + sample_name);
      ratio_hist_muon->SetDirectory(0);
    }
    else if(!isMuon){
      TFile* f_btag2Dsf_ele = new TFile("/nfs/dust/cms/user/deleokse/RunII_106_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/macros/src/files_BTagSF/customBtagSF_electron_"+year+".root");
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

bool ZprimeAnalysisModule::process(uhh2::Event& event){

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

  event.set(h_res_jet_bscore, -10);
  event.set(h_mer_subjet_bscore, -10);
  event.set(h_bscore_max, -10);

  event.set(h_res_jet_bscore_matchable, -10);
  event.set(h_mer_subjet_bscore_matchable, -10);
  event.set(h_bscore_max_matchable, -10);

  // Chi2 values are between 0 and ~10,000
  event.set(h_BestChi2, -10);
  event.set(h_BestChi2_resolved, -10);
  event.set(h_BestChi2_merged, -10);

  event.set(h_BestChi2_passedMatching, -10);
  event.set(h_BestChi2_passedMatching_resolved, -10);
  event.set(h_BestChi2_passedMatching_merged, -10);

  event.set(h_BestChi2_failedMatching, -10);
  event.set(h_BestChi2_failedMatching_resolved, -10);
  event.set(h_BestChi2_failedMatching_merged, -10);

  event.set(h_BestChi2hadronicleg, -10);
  event.set(h_BestChi2hadronicleg_resolved, -10);
  event.set(h_BestChi2hadronicleg_merged, -10);

  event.set(h_BestChi2hadronicleg_passedMatching, -10);
  event.set(h_BestChi2hadronicleg_passedMatching_resolved, -10);
  event.set(h_BestChi2hadronicleg_passedMatching_merged, -10);

  event.set(h_BestChi2hadronicleg_failedMatching, -10);
  event.set(h_BestChi2hadronicleg_failedMatching_resolved, -10);
  event.set(h_BestChi2hadronicleg_failedMatching_merged, -10);

  event.set(h_BestChi2leptonicleg, -10);
  event.set(h_BestChi2leptonicleg_resolved, -10);
  event.set(h_BestChi2leptonicleg_merged, -10);

  event.set(h_BestChi2leptonicleg_passedMatching, -10);
  event.set(h_BestChi2leptonicleg_passedMatching_resolved, -10);
  event.set(h_BestChi2leptonicleg_passedMatching_merged, -10);

  event.set(h_BestChi2leptonicleg_failedMatching, -10);
  event.set(h_BestChi2leptonicleg_failedMatching_resolved, -10);
  event.set(h_BestChi2leptonicleg_failedMatching_merged, -10);

  event.set(h_DrHadB, -10);
  event.set(h_DrHadB_resolved, -10);
  event.set(h_DrHadB_merged, -10);

  event.set(h_DrHadB_passedChi2cut, -10);
  event.set(h_DrHadB_passedChi2cut_resolved, -10);
  event.set(h_DrHadB_passedChi2cut_merged, -10);

  event.set(h_DrHadB_failedChi2cut, -10);
  event.set(h_DrHadB_failedChi2cut_resolved, -10);
  event.set(h_DrHadB_failedChi2cut_merged, -10);

  event.set(h_DrLep, -10);
  event.set(h_DrLep_resolved, -10);
  event.set(h_DrLep_merged, -10);

  event.set(h_DrLep_passedChi2cut, -10);
  event.set(h_DrLep_passedChi2cut_resolved, -10);
  event.set(h_DrLep_passedChi2cut_merged, -10);

  event.set(h_DrLep_failedChi2cut, -10);
  event.set(h_DrLep_failedChi2cut_resolved, -10);
  event.set(h_DrLep_failedChi2cut_merged, -10);

  // Require to pass CorrectMatch selection
  event.set(h_correctDr, -10);
  event.set(h_correctDr_resolved, -10);
  event.set(h_correctDr_merged, -10);

  event.set(h_correctDr_passedChi2cut, -10);
  event.set(h_correctDr_passedChi2cut_resolved, -10);
  event.set(h_correctDr_passedChi2cut_merged, -10);

  event.set(h_correctDr_failedChi2cut, -10);
  event.set(h_correctDr_failedChi2cut_resolved, -10);
  event.set(h_correctDr_failedChi2cut_merged, -10);

  event.set(h_DrHadB_matchable, -10);
  event.set(h_DrHadB_matchable_resolved, -10);
  event.set(h_DrHadB_matchable_merged, -10);

  event.set(h_DrHadB_matchable_passedChi2cut, -10);
  event.set(h_DrHadB_matchable_passedChi2cut_resolved, -10);
  event.set(h_DrHadB_matchable_passedChi2cut_merged, -10);

  event.set(h_DrHadB_matchable_failedChi2cut, -10);
  event.set(h_DrHadB_matchable_failedChi2cut_resolved, -10);
  event.set(h_DrHadB_matchable_failedChi2cut_merged, -10);

  event.set(h_DrLep_matchable, -10);
  event.set(h_DrLep_matchable_resolved, -10);
  event.set(h_DrLep_matchable_merged, -10);

  event.set(h_DrLep_matchable_passedChi2cut, -10);
  event.set(h_DrLep_matchable_passedChi2cut_resolved, -10);
  event.set(h_DrLep_matchable_passedChi2cut_merged, -10);

  event.set(h_DrLep_matchable_failedChi2cut, -10);
  event.set(h_DrLep_matchable_failedChi2cut_resolved, -10);
  event.set(h_DrLep_matchable_failedChi2cut_merged, -10);

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
  if(isdeepAK8) sf_topmistag->process(event);
  if(debug) cout << "Weights_TopMistag_SF: ok" << endl;
  fill_histograms(event, "Weights_TopMistag_SF");

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
  // if(!isEleTriggerMeasurement) SystematicsModule->process(event);


  // Begin Reconstruction Efficiency Study ------------------------------------------------------------------------------------------------------

  // Variable to access gen-ttbar
  TTbarGen ttbargen = event.get(h_ttbargen);

  // Check if ttbar decayed semileptonically
  if(ttbargen.DecayChannel() == TTbarGen::e_muhad || ttbargen.DecayChannel() == TTbarGen::e_ehad){
  
    // Tells me that ttbar was reconstructed using Chi2DiscriminatorZprime
    bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2);

    // Tells me that ttbar was reconstructed using CorrectMatchDiscriminatorZprime
    bool is_zprime_reconstructed_correctmatch = event.get(h_is_zprime_reconstructed_correctmatch);

    // Best Chi2 Candidate variables
    ZprimeCandidate* BestCandidate_chi2 = event.get(h_BestZprimeCandidateChi2);     // best ttbar candidate based on Chi2
    float bestChi2hadronicleg = BestCandidate_chi2->discriminator("chi2_hadronic"); // = pow((mhad - mtophad_) / sigmatophad_,2)
    float bestChi2leptonicleg = BestCandidate_chi2->discriminator("chi2_leptonic"); // = pow((mlep - mtoplep_) / sigmatoplep_,2)
    float bestChi2 =            BestCandidate_chi2->discriminator("chi2_total");    // = chi2_hadronic + chi2_leptonic

    // Variables for identifying btagged hadronic jet
    vector <Jet> AK4CHSjets_matched = event.get(h_CHSjets_matched);                 // AK4Puppijets that have been matched to CHSjets
    vector <float> jets_hadronic_bscores;                                           // bScores-vector for resolved hadronic jets
    float btag_medWP = 0.2783;    // medium WP for UL18
    // float btag_medWP = 0.3040; // medium WP for UL17

    // DeltaR between gen- and reco-level decay products
    float DrHadB, DrLep;


    // Begin Chi2 section --- 
    if(is_zprime_reconstructed_chi2){ // Events whose besCandidates were assigned a Chi2 value by Chi2DiscriminatorZprime

      bool is_toptag_reconstruction = BestCandidate_chi2->is_toptag_reconstruction();  // Reconstruction process id: to distinguish toplogy of ttbar reconstruction (0=resolved or 1=merged)

      // Variables for identifying btagged hadronic jet
      vector<Particle> jets_had_resolved = BestCandidate_chi2->jets_hadronic();        // vector of resolved hadronic jets

      // Define 4vector of hadronic b-jet
      LorentzVector HadTop_b;

      // Extract Hadronic bjet ------------------------------------------------------------------------------
      float bscore_max = -2;

      // Extract highest bscore of Resolved jets
      if(!is_toptag_reconstruction){
        // Loop over resolved hadronic jets to find their bscore via CHS jets
        for(unsigned int i=0; i<jets_had_resolved.size(); i++){
          double deltaR_min = 99;
          // Match resolved hadronic jets to CHS jets (which have bscores)
          for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
            double deltaR_CHS = deltaR(jets_had_resolved.at(i), AK4CHSjets_matched.at(j));
            if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;}
          // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
          for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
            if(deltaR(jets_had_resolved.at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
            jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet());}
        }
        // Loop over bScores-vector to extract highest bscore
        for(unsigned int i=0; i<jets_hadronic_bscores.size(); i++){
          float bscore = jets_hadronic_bscores.at(i);
          if(bscore > bscore_max) bscore_max = bscore;
        }
        // Plot all bscores of resolved top's jets
        for(unsigned int j=0; j<jets_hadronic_bscores.size(); j++){
          float res_jet_bscore = jets_hadronic_bscores.at(j);
          event.set(h_res_jet_bscore, res_jet_bscore);
        }
      }

      // Extract highest bscore of Merged jet
      if(is_toptag_reconstruction){
        // Loop over hadronic top's subjets to extract highest bscore
        for(unsigned int i=0; i < BestCandidate_chi2->tophad_topjet_ptr()->subjets().size(); i++){
          float bscore = BestCandidate_chi2->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet();
          if(bscore > bscore_max) bscore_max = bscore;
        }
        // Plot all bscores of merged top's subjets
        for(unsigned int j=0; j < BestCandidate_chi2->tophad_topjet_ptr()->subjets().size(); j++){
          float mer_subjet_bscore = BestCandidate_chi2->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
          event.set(h_mer_subjet_bscore, mer_subjet_bscore);
        }
      }

      // Cut on bscore_max >= btag_WP
      if(bscore_max >= btag_medWP){
        event.set(h_bscore_max, bscore_max); // Plot max bscores

        // Assign Hadronic bjet in Resolved topology
        if(!is_toptag_reconstruction){
          for(unsigned int i=0; i< jets_had_resolved.size(); i++){
            float bscore = jets_hadronic_bscores.at(i);
            if(bscore == bscore_max) HadTop_b = jets_had_resolved.at(i).v4();
          }
        }
        // Assign Hadronic bjet in Merged topology
        if(is_toptag_reconstruction){
          for(unsigned int j=0; j < BestCandidate_chi2->tophad_topjet_ptr()->subjets().size(); j++){
            float bscore = BestCandidate_chi2->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
            if(bscore == bscore_max) HadTop_b = BestCandidate_chi2->tophad_topjet_ptr()->subjets().at(j).v4();
          }
        }
      }
      // Hadronic bjet Extracted ------------------------------------------------------------------------------
      DrHadB = deltaR(ttbargen.BHad().v4(), HadTop_b);

      // Define 4vector of lepton
      LorentzVector LepTop_lepton = BestCandidate_chi2->lepton().v4();
      DrLep = deltaR(ttbargen.ChargedLepton().v4(), LepTop_lepton);


      if(bestChi2 != -10){
        // Fill Chi2 for
        event.set(h_BestChi2, bestChi2);   // all ttbar
        event.set(h_BestChi2hadronicleg, bestChi2hadronicleg);
        event.set(h_BestChi2leptonicleg, bestChi2leptonicleg);
        if(!is_toptag_reconstruction){   // resolved ttbar
          event.set(h_BestChi2_resolved, bestChi2); 
          event.set(h_BestChi2hadronicleg_resolved, bestChi2hadronicleg);
          event.set(h_BestChi2leptonicleg_resolved, bestChi2leptonicleg);
        }
        if(is_toptag_reconstruction){   // merged ttbar
          event.set(h_BestChi2_merged, bestChi2);   
          event.set(h_BestChi2hadronicleg_merged, bestChi2hadronicleg);
          event.set(h_BestChi2leptonicleg_merged, bestChi2leptonicleg);
        }

        if(is_zprime_reconstructed_correctmatch){ // Events that ALSO passed the matching requirements

          // Fill Chi2 for
          event.set(h_BestChi2_passedMatching, bestChi2);   // all ttbar
          event.set(h_BestChi2hadronicleg_passedMatching, bestChi2hadronicleg);
          event.set(h_BestChi2leptonicleg_passedMatching, bestChi2leptonicleg);
          if(!is_toptag_reconstruction){   // resolved ttbar
            event.set(h_BestChi2_passedMatching_resolved, bestChi2); 
            event.set(h_BestChi2hadronicleg_passedMatching_resolved, bestChi2hadronicleg);
            event.set(h_BestChi2leptonicleg_passedMatching_resolved, bestChi2leptonicleg);
          }
          if(is_toptag_reconstruction){   // merged ttbar
            event.set(h_BestChi2_passedMatching_merged, bestChi2);   
            event.set(h_BestChi2hadronicleg_passedMatching_merged, bestChi2hadronicleg);
            event.set(h_BestChi2leptonicleg_passedMatching_merged, bestChi2leptonicleg);
          }
        }
        else{ // Events that did NOT pass the matching requirements
          event.set(h_BestChi2_failedMatching, bestChi2);   // all ttbar
          event.set(h_BestChi2hadronicleg_failedMatching, bestChi2hadronicleg);
          event.set(h_BestChi2leptonicleg_failedMatching, bestChi2leptonicleg);
          if(!is_toptag_reconstruction){   // resolved ttbar
            event.set(h_BestChi2_failedMatching_resolved, bestChi2); 
            event.set(h_BestChi2hadronicleg_failedMatching_resolved, bestChi2hadronicleg);
            event.set(h_BestChi2leptonicleg_failedMatching_resolved, bestChi2leptonicleg);
          }
          if(is_toptag_reconstruction){   // merged ttbar
            event.set(h_BestChi2_failedMatching_merged, bestChi2);   
            event.set(h_BestChi2hadronicleg_failedMatching_merged, bestChi2hadronicleg);
            event.set(h_BestChi2leptonicleg_failedMatching_merged, bestChi2leptonicleg);
          }
        }
        
      }

      // Fill DeltaR for
      if(DrHadB != -10){
        event.set(h_DrHadB, DrHadB);   // all hadronic b-quarks
        if(!is_toptag_reconstruction) event.set(h_DrHadB_resolved, DrHadB);   // resolved hadronic b-quarks
        if(is_toptag_reconstruction) event.set(h_DrHadB_merged, DrHadB);   // merged hadronic b-quarks

        if(bestChi2 != -10 && bestChi2 < 30){ // Events that ALSO passed the chi2 cut
          event.set(h_DrHadB_passedChi2cut, DrHadB); // all hadronic b-quarks
          if(!is_toptag_reconstruction) event.set(h_DrHadB_passedChi2cut_resolved, DrHadB); // resolved hadronic b-quarks
          if(is_toptag_reconstruction) event.set(h_DrHadB_passedChi2cut_merged, DrHadB); // merged hadronic b-quarks
        }
        if(bestChi2 != -10 && bestChi2 >= 30){ // Events that did NOT pass the chi2 cut
          event.set(h_DrHadB_failedChi2cut, DrHadB); // all hadronic b-quarks
          if(!is_toptag_reconstruction) event.set(h_DrHadB_failedChi2cut_resolved, DrHadB); // resolved hadronic b-quarks
          if(is_toptag_reconstruction) event.set(h_DrHadB_failedChi2cut_merged, DrHadB); // merged hadronic b-quarks
        }
      }

      if(DrLep != -10){ 
        event.set(h_DrLep, DrLep);   // all leptons
        if(!is_toptag_reconstruction) event.set(h_DrLep_resolved, DrLep);   // resolved leptons
        if(is_toptag_reconstruction) event.set(h_DrLep_merged, DrLep);   // merged leptons

        if(is_zprime_reconstructed_chi2 && bestChi2 != -10 && bestChi2 < 30){ // Events that ALSO passed the chi2 cut
          event.set(h_DrLep_passedChi2cut, DrLep); // all leptons
          if(!is_toptag_reconstruction) event.set(h_DrLep_passedChi2cut_resolved, DrLep); // resolved leptons
          if(is_toptag_reconstruction) event.set(h_DrLep_passedChi2cut_merged, DrLep); // merged leptons
        }
        if(is_zprime_reconstructed_chi2 && bestChi2 != -10 && bestChi2 >= 30){ // Events that did NOT pass the chi2 cut
          event.set(h_DrLep_failedChi2cut, DrLep); // all hadronic b-quarks
          if(!is_toptag_reconstruction) event.set(h_DrLep_failedChi2cut_resolved, DrLep); // resolved leptons
          if(is_toptag_reconstruction) event.set(h_DrLep_failedChi2cut_merged, DrLep); // merged leptons
        }
      }

    }
    // End Chi2 section ---



    // Begin CorrectMatch section ---
    if(is_zprime_reconstructed_correctmatch){ // Events whose besCandidates passed matching requirements and were assigned a correct_dr value by CorrectMatchDiscriminatorZprime

      // Best CorrectMatch Candidate variables
      ZprimeCandidate* BestCandidate_correctMatch = event.get(h_BestZprimeCandidateCorrectMatch);   // best ttbar candidate based on Matching parameter dr
      float correct_dr = BestCandidate_correctMatch->discriminator("correct_match");                // = min(sum(dr)) between all pairs of gen-lvl decay products and corresponding reco-lvl objects

      // Variables for identifying btagged hadronic jet
      vector<Particle> jets_had_resolved = BestCandidate_correctMatch->jets_hadronic();             // vector of resolved hadronic jets
      bool is_toptag_reconstruction = BestCandidate_correctMatch->is_toptag_reconstruction();       // Reconstruction process id: to distinguish toplogy of ttbar reconstruction (0=resolved or 1=merged)

      // Define 4vector of hadronic b-jet
      LorentzVector HadTop_b;

      // Extract Hadronic bjet ------------------------------------------------------------------------------
      float bscore_max = -2;

      // Extract highest bscore of Resolved jets
      if(!is_toptag_reconstruction){
        // Loop over resolved hadronic jets to find their bscore via CHS jets
        for(unsigned int i=0; i<jets_had_resolved.size(); i++){
          double deltaR_min = 99;
          // Match resolved hadronic jets to CHS jets (which have bscores)
          for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
            double deltaR_CHS = deltaR(jets_had_resolved.at(i), AK4CHSjets_matched.at(j));
            if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;}
          // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
          for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
            if(deltaR(jets_had_resolved.at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
            jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet());}
        }
        // Loop over bScores-vector to extract highest bscore
        for(unsigned int i=0; i<jets_hadronic_bscores.size(); i++){
          float bscore = jets_hadronic_bscores.at(i);
          if(bscore > bscore_max) bscore_max = bscore;
        }
        // Plot all bscores of resolved top's jets
        for(unsigned int j=0; j<jets_hadronic_bscores.size(); j++){
          float res_jet_bscore = jets_hadronic_bscores.at(j);
          event.set(h_res_jet_bscore_matchable, res_jet_bscore);
        }
      }

      // Extract highest bscore of Merged jet
      if(is_toptag_reconstruction){
        // Loop over hadronic top's subjets to extract highest bscore
        for(unsigned int i=0; i < BestCandidate_correctMatch->tophad_topjet_ptr()->subjets().size(); i++){
          float bscore = BestCandidate_correctMatch->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet();
          if(bscore > bscore_max) bscore_max = bscore;
        }
        // Plot all bscores of merged top's subjets
        for(unsigned int j=0; j < BestCandidate_correctMatch->tophad_topjet_ptr()->subjets().size(); j++){
          float mer_subjet_bscore = BestCandidate_correctMatch->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
          event.set(h_mer_subjet_bscore_matchable, mer_subjet_bscore);
        }
      }

      // Cut on bscore_max >= btag_WP
      if(bscore_max >= btag_medWP){
        event.set(h_bscore_max_matchable, bscore_max); // Plot max bscores

        // Assign Hadronic bjet in Resolved topology
        if(!is_toptag_reconstruction){
          for(unsigned int i=0; i< jets_had_resolved.size(); i++){
            float bscore = jets_hadronic_bscores.at(i);
            if(bscore == bscore_max) HadTop_b = jets_had_resolved.at(i).v4();
          }
        }
        // Assign Hadronic bjet in Merged topology
        if(is_toptag_reconstruction){
          for(unsigned int j=0; j < BestCandidate_correctMatch->tophad_topjet_ptr()->subjets().size(); j++){
            float bscore = BestCandidate_correctMatch->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
            if(bscore == bscore_max) HadTop_b = BestCandidate_correctMatch->tophad_topjet_ptr()->subjets().at(j).v4();
          }
        }
      }
      // Hadronic bjet Extracted ------------------------------------------------------------------------------
      DrHadB = deltaR(ttbargen.BHad().v4(), HadTop_b);

      // Define 4vector of lepton
      LorentzVector LepTop_lepton = BestCandidate_correctMatch->lepton().v4();
      DrLep = deltaR(ttbargen.ChargedLepton().v4(), LepTop_lepton);


      if(correct_dr != -10 && correct_dr != 9999999){

        // Fill correct_dr for
        event.set(h_correctDr, correct_dr); // all ttbar
        if(!is_toptag_reconstruction) event.set(h_correctDr_resolved, correct_dr); // resolved ttbar
        if(is_toptag_reconstruction) event.set(h_correctDr_merged, correct_dr); // merged ttbar

        if(is_zprime_reconstructed_chi2 && bestChi2 != -10 && bestChi2 < 30){ // Events that ALSO passed the chi2 cut
          event.set(h_correctDr_passedChi2cut, correct_dr); // all ttbar
          if(!is_toptag_reconstruction) event.set(h_correctDr_passedChi2cut_resolved, correct_dr); // resolved ttbar
          if(is_toptag_reconstruction) event.set(h_correctDr_passedChi2cut_merged, correct_dr); // merged ttbar
        }
        if(is_zprime_reconstructed_chi2 && bestChi2 != -10 && bestChi2 >= 30){ // Events that failed the chi2 cut
          event.set(h_correctDr_failedChi2cut, correct_dr); // all ttbar
          if(!is_toptag_reconstruction) event.set(h_correctDr_failedChi2cut_resolved, correct_dr); // resolved ttbar
          if(is_toptag_reconstruction) event.set(h_correctDr_failedChi2cut_merged, correct_dr); // merged ttbar
        }

      }

      // Fill DeltaR for
      if(DrHadB != -10){
        event.set(h_DrHadB_matchable, DrHadB);   // all hadronic b-quarks
        if(!is_toptag_reconstruction) event.set(h_DrHadB_matchable_resolved, DrHadB);   // resolved hadronic b-quarks
        if(is_toptag_reconstruction) event.set(h_DrHadB_matchable_merged, DrHadB);   // merged hadronic b-quarks

        if(is_zprime_reconstructed_chi2 && bestChi2 != -10 && bestChi2 < 30){ // Events that ALSO passed the chi2 cut
          event.set(h_DrHadB_matchable_passedChi2cut, DrHadB); // all hadronic b-quarks
          if(!is_toptag_reconstruction) event.set(h_DrHadB_matchable_passedChi2cut_resolved, DrHadB); // resolved hadronic b-quarks
          if(is_toptag_reconstruction) event.set(h_DrHadB_matchable_passedChi2cut_merged, DrHadB); // merged hadronic b-quarks
        }
        if(is_zprime_reconstructed_chi2 && bestChi2 != -10 && bestChi2 >= 30){ // Events that did NOT pass the chi2 cut
          event.set(h_DrHadB_matchable_failedChi2cut, DrHadB); // all hadronic b-quarks
          if(!is_toptag_reconstruction) event.set(h_DrHadB_matchable_failedChi2cut_resolved, DrHadB); // resolved hadronic b-quarks
          if(is_toptag_reconstruction) event.set(h_DrHadB_matchable_failedChi2cut_merged, DrHadB); // merged hadronic b-quarks
        }
      }

      if(DrLep != -10){ 
        event.set(h_DrLep_matchable, DrLep);   // all leptons
        if(!is_toptag_reconstruction) event.set(h_DrLep_matchable_resolved, DrLep);   // resolved leptons
        if(is_toptag_reconstruction) event.set(h_DrLep_matchable_merged, DrLep);   // merged leptons

        if(is_zprime_reconstructed_chi2 && bestChi2 != -10 && bestChi2 < 30){ // Events that ALSO passed the chi2 cut
          event.set(h_DrLep_matchable_passedChi2cut, DrLep); // all leptons
          if(!is_toptag_reconstruction) event.set(h_DrLep_matchable_passedChi2cut_resolved, DrLep); // resolved leptons
          if(is_toptag_reconstruction) event.set(h_DrLep_matchable_passedChi2cut_merged, DrLep); // merged leptons
        }
        if(is_zprime_reconstructed_chi2 && bestChi2 != -10 && bestChi2 >= 30){ // Events that did NOT pass the chi2 cut
          event.set(h_DrLep_matchable_failedChi2cut, DrLep); // all hadronic b-quarks
          if(!is_toptag_reconstruction) event.set(h_DrLep_matchable_failedChi2cut_resolved, DrLep); // resolved leptons
          if(is_toptag_reconstruction) event.set(h_DrLep_matchable_failedChi2cut_merged, DrLep); // merged leptons
        }
      }

    }
    // End CorrectMatch section ---

  }

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule)
