#include <iostream>
#include <memory>

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

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/TTbarLJHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>

//#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

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
  
  // Cleaners (to make sure the pre-selection did everything right)
  std::unique_ptr<MuonCleaner> muon_cleaner;
  std::unique_ptr<ElectronCleaner> electron_cleaner;

  // Scale Factors -- Systematics
  unique_ptr<MCMuonScaleFactor> MuonID_module, MuonTrigger_module;
  unique_ptr<MCElecScaleFactor> EleID_module, EleTrigger_module, EleRec_module;
  unique_ptr<MCToptaggSF> ToptagSF_module;
  unique_ptr<MCNjetsHTScaleFactor> HT_module; 
   
  // AnalysisModules
  unique_ptr<AnalysisModule> LumiWeight_module, PUWeight_module, printer_genparticles, BTagWeight_module, TopPtReweight_module, MCScale_module;

  // Taggers
  unique_ptr<AK8PuppiTopTagger> TopTaggerPuppi;
  unique_ptr<AK8PuppiWTagger> WTaggerPuppi;

  // Mass reconstruction
  unique_ptr<ZprimeCandidateBuilder> CandidateBuilder;

  // Chi2 discriminator
  unique_ptr<ZprimeChi2Discriminator> Chi2DiscriminatorZprime;
  unique_ptr<ZprimeCorrectMatchDiscriminator> CorrectMatchDiscriminatorZprime;

  // Selections
  unique_ptr<Selection> Trigger1_selection, Trigger2_selection, NMuon1_selection, NMuon2_selection, NElectron_selection, TwoDCut_selection, Jet1_selection, Jet2_selection, Met_selection, Chi2_selection, TTbarMatchable_selection, Chi2CandidateMatched_selection, ZprimeTopTag_selection, BlindData_selection;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  std::unique_ptr<Selection> sel_1btag, sel_2btag;
  std::unique_ptr<Selection> TopJetBtagSubjet_selection;
  std::unique_ptr<Selection> HEM_selection;

  // NN variables handles
  unique_ptr<Variables_NN> Variables_module;

  //Handles
  Event::Handle<bool> h_is_zprime_reconstructed_chi2, h_is_zprime_reconstructed_correctmatch;
  Event::Handle<float> h_chi2;   Event::Handle<float> h_weight;
  Event::Handle<float> h_MET;   Event::Handle<int> h_NPV;
  Event::Handle<float> h_lep1_pt; Event::Handle<float> h_lep1_eta; 
  Event::Handle<float> h_ak4jet1_pt; Event::Handle<float> h_ak4jet1_eta; 
  Event::Handle<float> h_ak8jet1_pt; Event::Handle<float> h_ak8jet1_eta; 
  Event::Handle<float> h_Mttbar;
  Event::Handle<float> h_TH_M;
  Event::Handle<float> h_TL_M; 
  Event::Handle<int> h_ttagN;
  Event::Handle<int> h_btagN;
  Event::Handle<int> h_wtagN;
  Event::Handle<float> h_matched_ttbar; 
  Event::Handle<float> h_DeltaY;
  Event::Handle<float> h_ISRup;
  Event::Handle<float> h_ISRdown;
  Event::Handle<float> h_FSRup;
  Event::Handle<float> h_FSRdown;
  Event::Handle<float> h_pTttbar;
  Event::Handle<float> h_Rapidityttbar;
  Event::Handle<float> h_pThad;
  Event::Handle<float> h_pTlep;
  Event::Handle<float> h_Rapidityhad;
  Event::Handle<float> h_Rapiditylep;
  Event::Handle<float> h_DeltaPhi_thad_lepton;
  Event::Handle<float> h_Cos_theta_thad_tlep;
  Event::Handle<float> h_tlead_pT;
  Event::Handle<float> h_tsublead_pT;
  Event::Handle<float> h_muonrecSF_nominal;
  Event::Handle<float> h_muonrecSF_up;
  Event::Handle<float> h_muonrecSF_down;
  Event::Handle<std::vector<float> > h_wgtMC__PDF;

  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
  uhh2::Event::Handle<std::vector<TopJet>> h_AK8TopTags;
  uhh2::Event::Handle<std::vector<TopJet>> h_AK8WTags;
 
  // Configuration
  bool isMC, ispuppi, islooserselection;
  string Sys_MuonID, Sys_MuonTrigger, Sys_PU, Sys_btag, Sys_EleID, Sys_EleTrigger;
  TString sample;
  int runnr_oldtriggers = 299368;

  bool is2016v2, is2016v3, is2017v2, is2018;
  bool isMuon, isElectron;
};

void ZprimeAnalysisModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag + "_Skimming";
    mytag = tag+"_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicHists(ctx,mytag));
    mytag = tag+"_Generator";
    book_HFolder(mytag, new ZprimeSemiLeptonicGeneratorHists(ctx,mytag));
  }
}

void ZprimeAnalysisModule::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag + "_Skimming";
  mytag = tag+"_General";
  HFolder(mytag)->fill(event);
  mytag = tag+"_Generator";
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
  //debug = true;
  debug = false;
  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }
  // Configuration
  isMC = (ctx.get("dataset_type") == "MC");
  ispuppi = (ctx.get("is_puppi") == "true");
  TString mode = "chs";
  if(ispuppi) mode = "puppi";
  string tmp = ctx.get("dataset_version");
  sample = tmp;
  is2016v2 = (ctx.get("dataset_version").find("2016v2") != std::string::npos);
  is2016v3 = (ctx.get("dataset_version").find("2016v3") != std::string::npos);
  is2017v2 = (ctx.get("dataset_version").find("2017v2") != std::string::npos);
  is2018 = (ctx.get("dataset_version").find("2018") != std::string::npos);

  // Important selection values
  islooserselection = (ctx.get("is_looser_selection") == "true");
  double muon_pt(30.);
  double elec_pt(30.);
  double jet1_pt(150.);
  double jet2_pt(50.);
  double chi2_max(300.);
  double mtt_blind(3000.);
  int nmuon_min1, nmuon_max1;
  int nmuon_min2, nmuon_max2;
  int nele_min, nele_max;
  string trigger1,trigger2;
  double MET_cut, HT_lep_cut;
  isMuon = false; isElectron = false;
  if(ctx.get("channel") == "muon") isMuon = true;
  if(ctx.get("channel") == "electron") isElectron = true;

  if(isMuon){//semileptonic muon channel
    trigger1 = "HLT_Mu50_v*";
    if(is2016v2 || is2016v3)
      //trigger2 = "HLT_TkMu50_v*";
      trigger2 = "HLT_Mu50_v*"; //TkMu path does not exist in 2017/2018 and RunB 2016
    else
      trigger2 = "HLT_Mu50_v*"; //TkMu path does not exist in 2017/2018
    nmuon_min1 = 1, nmuon_max1 = -1;
    nmuon_min2 = 1, nmuon_max2 = 1;
    nele_min = 0; nele_max = 0;
    MET_cut = 50;
    jet1_pt = 150.;
    HT_lep_cut = 150;
  }
  if(isElectron){//semileptonic electron channel
    nmuon_min1 = 0; nmuon_max1 = 0;
    nmuon_min2 = 0; nmuon_max2 = 0;
    nele_min = 1; nele_max = 1;
    trigger1 = "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*";
    trigger2 = "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*";
    MET_cut = 120;
    jet1_pt = 185.;
    HT_lep_cut = 0;
  }

 
  double TwoD_dr = 0.4, TwoD_ptrel = 25.;
  if(islooserselection){
    jet1_pt = 100.;
    TwoD_dr = 0.2;
    TwoD_ptrel = 10.;
  }
  const MuonId muonID(PtEtaCut(muon_pt, 2.4));
  const ElectronId electronID(PtEtaSCCut(elec_pt, 2.5));

  

  Sys_MuonID = ctx.get("Sys_MuonID");
  Sys_MuonTrigger = ctx.get("Sys_MuonTrigger");
  Sys_EleID = ctx.get("Sys_EleID");
  Sys_EleTrigger = ctx.get("Sys_EleTrigger");
  Sys_PU = ctx.get("Sys_PU");
  Sys_btag = ctx.get("Sys_BTagSF");

  BTag::algo btag_algo = BTag::DEEPJET;
  BTag::wp btag_wp_tight = BTag::WP_TIGHT;
  JetId id_btag = BTag(btag_algo, btag_wp_tight);

  double a_toppt = 0.0615; // par a TopPt Reweighting
  double b_toppt = -0.0005; // par b TopPt Reweighting 

  // Modules
  printer_genparticles.reset(new GenParticlesPrinter(ctx));
  muon_cleaner.reset(new MuonCleaner(muonID));
  electron_cleaner.reset(new ElectronCleaner(electronID));
  LumiWeight_module.reset(new MCLumiWeight(ctx));
  PUWeight_module.reset(new MCPileupReweight(ctx, Sys_PU));
  BTagWeight_module.reset(new MCBTagDiscriminantReweighting(ctx, btag_algo, "jets", Sys_btag));
  TopPtReweight_module.reset(new TopPtReweight(ctx, a_toppt, b_toppt));
  MCScale_module.reset(new MCScaleVariation(ctx));


  if((is2016v3 || is2016v2) && isMuon){
    MuonID_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root", "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt", 0., "HighPtID", false, Sys_MuonID));
    MuonTrigger_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu50_OR_IsoTkMu50_PtEtaBins", 0.5, "MuonTrigger", false, Sys_MuonTrigger));
    HT_module.reset(new MCNjetsHTScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/HT_SF_muon.root", 1.0, "HT", Sys_EleID));
  }
  if(is2017v2 && isMuon){
    MuonID_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/RunBCDEF_SF_ID.root", "NUM_HighPtID_DEN_genTracks_pair_newTuneP_probe_pt_abseta", 0., "HighPtID", true, Sys_MuonID));
    MuonTrigger_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/MuonTrigger_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root", "Mu50_PtEtaBins/pt_abseta_ratio", 0.5, "Trigger", true, Sys_MuonTrigger));
    HT_module.reset(new MCNjetsHTScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/HT_SF_muon.root", 1.0, "HT", Sys_EleID));
  }
  if(is2018 && isMuon){
    MuonID_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/Muon_ID_SF_RunABCD.root", "NUM_HighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta", 0., "HighPtID", true, Sys_MuonID));
    MuonTrigger_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/Muon_Trigger_Eff_SF_AfterMuonHLTUpdate.root", "Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/pt_abseta_ratio", 0.5, "Trigger", true, Sys_MuonTrigger)); 
    HT_module.reset(new MCNjetsHTScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/HT_SF_muon.root", 1.0, "HT", Sys_EleID));
  }

  if((is2016v3 || is2016v2) && isElectron){
    EleID_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root", 1.0, "TightID", Sys_EleID));
    EleTrigger_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/SF_Ele50_Ele115_2016.root", 0.5, "Trigger", Sys_EleTrigger, "electrons", "abseta_pt_ratio"));
    EleRec_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1.0, "Rec", Sys_EleID));
    HT_module.reset(new MCNjetsHTScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/HT_SF_electron.root", 1.0, "HT", Sys_EleID));
  }

  if(is2017v2 && isElectron){
    EleID_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/2017_ElectronTight.root", 1.0, "TightID", Sys_EleID));
    EleTrigger_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/SF_Ele50_Ele115_2017.root", 0.5, "Trigger", Sys_EleTrigger, "electrons", "abseta_pt_ratio"));
    EleRec_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", 1.0, "Rec", Sys_EleID));
    HT_module.reset(new MCNjetsHTScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/HT_SF_electron.root", 1.0, "HT", Sys_EleID));
  }


  if(is2018 && isElectron){
    cout << "is electron" << endl;
    EleID_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/RunII_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/2018_ElectronTight.root", 1.0, "TightID", Sys_EleID));
    EleTrigger_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/SF_Ele50_Ele115_2018.root", 0.5, "Trigger", Sys_EleTrigger, "electrons", "abseta_pt_ratio"));
    EleRec_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/RunII_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/egammaEffi.txt_EGM2D_updatedAll.root", 1.0, "Rec", Sys_EleID));
    HT_module.reset(new MCNjetsHTScaleFactor(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/HT_SF_electron.root", 1.0, "HT", Sys_EleID));
  }


  if(is2018 || is2017v2 || is2016v3){
      ToptagSF_module.reset(new MCToptaggSF(ctx, "/nfs/dust/cms/user/hugobg/UHH2_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/2018TopTaggingScaleFactors.root", "PUPPI_wp4"));   
  }

  // Selection modules
  Trigger1_selection.reset(new TriggerSelection(trigger1));
  Trigger2_selection.reset(new TriggerSelection(trigger2));
  NMuon1_selection.reset(new NMuonSelection(nmuon_min1, nmuon_max1));
  NMuon2_selection.reset(new NMuonSelection(nmuon_min2, nmuon_max2));
  NElectron_selection.reset(new NElectronSelection(nele_min, nele_max));
  TwoDCut_selection.reset(new TwoDCut(TwoD_dr, TwoD_ptrel));
  Jet1_selection.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.4))));
  Jet2_selection.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.4))));
  met_sel.reset(new METCut  (MET_cut   , uhh2::infinity));
  htlep_sel.reset(new HTlepCut(HT_lep_cut, uhh2::infinity));

  Chi2_selection.reset(new Chi2Cut(ctx, 0., chi2_max));
  TTbarMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
  Chi2CandidateMatched_selection.reset(new Chi2CandidateMatchedSelection(ctx));
  ZprimeTopTag_selection.reset(new ZprimeTopTagSelection(ctx));
  BlindData_selection.reset(new BlindDataSelection(ctx, mtt_blind));

  HEM_selection.reset(new HEMSelection(ctx)); // HEM issue in 2018, veto on leptons and jets

  Variables_module.reset(new Variables_NN(ctx)); // variables for NN

  // Taggers
  TopTaggerPuppi.reset(new AK8PuppiTopTagger(ctx));
  WTaggerPuppi.reset(new AK8PuppiWTagger(ctx));

  // Zprime candidate builder
  CandidateBuilder.reset(new ZprimeCandidateBuilder(ctx, mode));

  // Zprime discriminators
  Chi2DiscriminatorZprime.reset(new ZprimeChi2Discriminator(ctx));
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  CorrectMatchDiscriminatorZprime.reset(new ZprimeCorrectMatchDiscriminator(ctx));
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
  h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("AK8PuppiTopTags");
  h_AK8WTags = ctx.get_handle<std::vector<TopJet>>("AK8PuppiWTags");

  h_chi2 = ctx.declare_event_output<float> ("rec_chi2");
  h_MET = ctx.declare_event_output<float> ("met_pt");
  h_Mttbar = ctx.declare_event_output<float> ("Mttbar");
  h_TH_M = ctx.declare_event_output<float> ("TH_M");
  h_TL_M = ctx.declare_event_output<float> ("TL_M");
  h_ttagN = ctx.declare_event_output<int>("ttagN");
  h_btagN = ctx.declare_event_output<int>("btagN");
  h_wtagN = ctx.declare_event_output<int>("wtagN");
  h_matched_ttbar = ctx.declare_event_output<float>("matched_ttbar");
  h_DeltaY = ctx.declare_event_output<float> ("DeltaY");
  h_ISRup = ctx.declare_event_output<float> ("ISRup");
  h_ISRdown = ctx.declare_event_output<float> ("ISRdown");
  h_FSRup = ctx.declare_event_output<float> ("FSRup");
  h_FSRdown = ctx.declare_event_output<float> ("FSRdown");
  h_pTttbar = ctx.declare_event_output<float> ("pT_ttbar");
  h_Rapidityttbar= ctx.declare_event_output<float> ("Rapidity_ttbar");
  h_pThad= ctx.declare_event_output<float> ("pT_had");
  h_pTlep= ctx.declare_event_output<float> ("pT_lep");
  h_Rapidityhad= ctx.declare_event_output<float> ("Rapidity_had");
  h_Rapiditylep= ctx.declare_event_output<float> ("Rapidity_lep");
  h_DeltaPhi_thad_lepton = ctx.declare_event_output<float> ("DeltaPhi_thad_lepton");
  h_Cos_theta_thad_tlep = ctx.declare_event_output<float> ("Cos_theta_thad_tlep");
  h_tlead_pT = ctx.declare_event_output<float> ("tlead_pT");
  h_tsublead_pT = ctx.declare_event_output<float> ("tsublead_pT");
  h_lep1_pt = ctx.declare_event_output<float> ("lep1_pt");
  h_lep1_eta = ctx.declare_event_output<float> ("lep1_eta");
  h_ak4jet1_pt = ctx.declare_event_output<float> ("ak4jet1_pt");
  h_ak4jet1_eta = ctx.declare_event_output<float> ("ak4jet1_eta");
  h_ak8jet1_pt = ctx.declare_event_output<float> ("ak8jet1_pt");
  h_ak8jet1_eta = ctx.declare_event_output<float> ("ak8jet1_eta");

  h_NPV = ctx.declare_event_output<int> ("NPV");
  h_weight = ctx.declare_event_output<float> ("weight");
  h_muonrecSF_nominal = ctx.declare_event_output<float> ("muonrecSF_nominal");
  h_muonrecSF_up = ctx.declare_event_output<float> ("muonrecSF_up");
  h_muonrecSF_down = ctx.declare_event_output<float> ("muonrecSF_down");
  h_wgtMC__PDF = ctx.declare_event_output<std::vector<float> >("wgtMC__PDF");


  sel_1btag.reset(new NJetSelection(1, 1, id_btag));
  sel_2btag.reset(new NJetSelection(2,-1, id_btag));


  TopJetBtagSubjet_selection.reset(new ZprimeBTagFatSubJetSelection(ctx));

  // Book histograms
  vector<string> histogram_tags = {"Weights", "Weights_MuonID", "Weights_PU", "Weights_Lumi", "Weights_TopPt", "Weights_MCScale", "Muon1", "TriggerMuon", "Muon2", "Electron1", "TriggerEle", "TwoDCut", "Jet1", "Jet2", "MET", "HTlep", "NNInputsBeforeReweight", "MatchableBeforeChi2Cut", "NotMatchableBeforeChi2Cut", "CorrectMatchBeforeChi2Cut", "NotCorrectMatchBeforeChi2Cut", "Chi2", "Matchable", "NotMatchable", "CorrectMatch", "NotCorrectMatch", "TopTagReconstruction", "NotTopTagReconstruction", "Btags2", "Btags1","TopJetBtagSubjet"};
  book_histograms(ctx, histogram_tags);
}

/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool ZprimeAnalysisModule::process(uhh2::Event& event){
  if(debug)   cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug)   cout<<" run.event: "<<event.run<<". "<<event.event<<endl;
  // Initialize reco flags with false
  event.set(h_is_zprime_reconstructed_chi2, false);
  event.set(h_is_zprime_reconstructed_correctmatch, false);
  event.set(h_chi2,-100);
  event.set(h_MET,-100);
  event.set(h_Mttbar,-100);
  event.set(h_TH_M,-100);
  event.set(h_TL_M,-100);
  event.set(h_ttagN,-100);
  event.set(h_wtagN,-100);
  event.set(h_btagN,-100);
  event.set(h_matched_ttbar,-100);
  event.set(h_DeltaY,-100);
  event.set(h_FSRup,-100);
  event.set(h_FSRdown,-100);
  event.set(h_ISRup,-100);
  event.set(h_ISRdown,-100);
  event.set(h_pTttbar,-100);
  event.set(h_Rapidityttbar,-100);
  event.set(h_pThad,-100);
  event.set(h_pTlep,-100);
  event.set(h_Rapidityhad,-100);
  event.set(h_Rapiditylep,-100);
  event.set(h_DeltaPhi_thad_lepton,-100);
  event.set(h_Cos_theta_thad_tlep,-100);
  event.set(h_tlead_pT,-100);
  event.set(h_tsublead_pT,-100);
  event.set(h_lep1_pt,-100);
  event.set(h_lep1_eta,-100);
  event.set(h_ak4jet1_pt,-100);
  event.set(h_ak4jet1_eta,-100);
  event.set(h_ak8jet1_pt,-100);
  event.set(h_ak8jet1_eta,-100);
  event.set(h_NPV,-100);
  event.set(h_weight,-100);
  event.set(h_muonrecSF_nominal,-100);  
  event.set(h_muonrecSF_up,-100);
  event.set(h_muonrecSF_down,-100);

  std::vector<float> w_PDF;
  w_PDF.clear();

  // Printing
  // if(!event.isRealData) printer_genparticles->process(event);

  // TODO Apply things that should've been done in the pre-selection already... Fix pre-selection and then remove these steps


  if(isMuon) muon_cleaner->process(event);
  if(isElectron) electron_cleaner->process(event);


  if(isMuon){
    if(!NMuon1_selection->passes(event)) return false;
    if(!NMuon2_selection->passes(event)) return false;
  }

 
  if(isElectron){
    if(!NElectron_selection->passes(event)) return false;
  }

  if(isMuon){
  if(event.muons->at(0).pt() < 55) return false;
  }

  if(isElectron){
  if(event.electrons->at(0).pt() < 80) return false;
  }



//  if(debug)  cout<<"Muon and Electron cleaner ok"<<endl;


  if(!HEM_selection->passes(event)){
    if(!isMC) return false;
    else event.weight = event.weight*(1-0.64774715284); // calculated following instructions ar https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis
  }

  // Run top-tagging
  TopTaggerPuppi->process(event);
  if(debug) cout<<"Top Tagger ok"<<endl;
  //Run W-tagging
  WTaggerPuppi->process(event); 
  if(debug) cout << "W Tagger ok" << endl;
 
  fill_histograms(event, "Weights");

  // Weight modules
  if(isMuon){
    MuonID_module->process(event);
    if(debug)  cout<<"MuonID ok"<<endl;
  }
  fill_histograms(event, "Weights_MuonID");
  if(isElectron){
    EleID_module->process(event);
    if(debug)  cout<<"EleID ok"<<endl;
  }

  if(isElectron){
    EleRec_module->process(event);
    if(debug)  cout<<"EleID ok"<<endl;
  }

  HT_module->process(event);
  ToptagSF_module->process(event);     

  PUWeight_module->process(event);
  if(debug)  cout<<"PUWeight ok"<<endl;
  fill_histograms(event, "Weights_PU");

  LumiWeight_module->process(event);
  if(debug)  cout<<"LumiWeight ok"<<endl;
  fill_histograms(event, "Weights_Lumi");

//  if(isMC){
  TopPtReweight_module->process(event);
  if(debug)  cout<<"TopPt ok"<<endl;
  fill_histograms(event, "Weights_TopPt");
//  }

  MCScale_module->process(event);
    if(debug)  cout<<"MCsCWeight ok"<<endl;
  fill_histograms(event, "Weights_MCScale");
 
  BTagWeight_module->process(event);

  if(!(Trigger1_selection->passes(event)|| Trigger2_selection->passes(event))) return false;
  if(isMuon){
    if(!NMuon1_selection->passes(event)) return false;
    fill_histograms(event, "Muon1");
    MuonTrigger_module->process_onemuon(event, 0);
    fill_histograms(event, "TriggerMuon");
    if(!NMuon2_selection->passes(event)) return false;
    fill_histograms(event, "Muon2");
  }
  if(isElectron){
    if(!NElectron_selection->passes(event)) return false;
    fill_histograms(event, "Electron1");
    EleTrigger_module->process(event);
    fill_histograms(event, "TriggerEle");
  }

  if(isMuon){
  if(event.muons->at(0).pt() < 55) return false;
  }

  if(isElectron){
  if(event.electrons->at(0).pt() < 80) return false;
  }

  if((event.muons->size()+event.electrons->size()) != 1) return false; //veto events without leptons or with too many 
  if(debug) cout<<"N leptons ok: Nelectrons="<<event.electrons->size()<<" Nmuons="<<event.muons->size()<<endl;


  if(!TwoDCut_selection->passes(event)) return false;
  fill_histograms(event, "TwoDCut");

  // Here, the Zprime must be reconstructed (we ensured to have >= 2 AK4 jets, >= 1 muon)
  // Only consider well-separated AK4 jets
  // //FixME: 
  // for(const TopJet & chsjet : *event.topjets){
  //   int n_separatedjets = 0;
  //   for(unsigned int k = 0; k < event.jets->size(); k++){
  //     if(deltaR(event.jets->at(k), chsjet) > 1.2) n_separatedjets++;
  //   }
  //   //    cout<<"n_separatedjets = "<<n_separatedjets<<endl;
  //   if(n_separatedjets<1) return false;
  // }

  CandidateBuilder->process(event);
  if(debug) cout<<"CandidateBuilder is ok"<<endl;
  Chi2DiscriminatorZprime->process(event);
  if(debug)  cout<<"Chi2DiscriminatorZprime is ok"<<endl;
  CorrectMatchDiscriminatorZprime->process(event);
  if(debug) cout<<"CorrectMatchDiscriminatorZprime is ok"<<endl;
  if(sample.Contains("_blinded")){
    if(!BlindData_selection->passes(event)) return false;
  }
  if(!Jet1_selection->passes(event)) return false;
  if(debug) cout<<"Jet1_selection is ok"<<endl;
  fill_histograms(event, "Jet1");

  if(!Jet2_selection->passes(event)) return false;
  if(debug) cout<<"Jet2_selection is ok"<<endl;
  fill_histograms(event, "Jet2");

  // MET selection
  if(!met_sel->passes(event)) return false;
  if(debug) cout<<"MET is ok"<<endl;
  fill_histograms(event, "MET");
  if(isMuon){
    if(!htlep_sel->passes(event)) return false;
    fill_histograms(event, "HTlep");
    if(debug) cout<<"HTlep is ok"<<endl;
  }

  //  Variables for NN 
  Variables_module->process(event);
  fill_histograms(event, "NNInputsBeforeReweight");

  if(TTbarMatchable_selection->passes(event)) fill_histograms(event, "MatchableBeforeChi2Cut");
  else fill_histograms(event, "NotMatchableBeforeChi2Cut");
  if(debug) cout<<"TTbarMatchable_selection is ok"<<endl;

  if(Chi2CandidateMatched_selection->passes(event)) fill_histograms(event, "CorrectMatchBeforeChi2Cut");
  else fill_histograms(event, "NotCorrectMatchBeforeChi2Cut");
  if(debug) cout<<"Chi2CandidateMatched_selection is ok"<<endl;

  if(!Chi2_selection->passes(event)) return false;

  fill_histograms(event, "Chi2");

  if(debug) cout<<"Chi2_selection is ok"<<endl;

  if(TTbarMatchable_selection->passes(event)) fill_histograms(event, "Matchable");
  else fill_histograms(event, "NotMatchable");
  if(debug) cout<<"TTbarMatchable_selection is ok"<<endl;

  if(Chi2CandidateMatched_selection->passes(event)) fill_histograms(event, "CorrectMatch");
  else fill_histograms(event, "NotCorrectMatch");
  if(debug) cout<<"Chi2CandidateMatched_selection is ok"<<endl;

  if(ZprimeTopTag_selection->passes(event)) fill_histograms(event, "TopTagReconstruction");
  else fill_histograms(event, "NotTopTagReconstruction");
  if(debug) cout<<"ZprimeTopTag_selection is ok"<<endl;

  //Test with b-tagging
  if(sel_2btag->passes(event)) fill_histograms(event, "Btags2");
  if(sel_1btag->passes(event)) fill_histograms(event, "Btags1");
  if(debug) cout<<"Btags1 is ok"<<endl;
  if(TopJetBtagSubjet_selection->passes(event)) fill_histograms(event, "TopJetBtagSubjet");
  if(debug) cout<<"TopJetBtagSubjet_selection is ok"<<endl;

  vector<TopJet> WTags = event.get(h_AK8WTags); 
  event.set(h_wtagN,WTags.size());

  vector<TopJet> TopTags = event.get(h_AK8TopTags);
  event.set(h_ttagN,TopTags.size());

  int btagN = 0;
  int jetbtagN = 0;
  int btag_subjet=0;

  if(TopTags.size() == 0){
     if(sel_2btag->passes(event)) jetbtagN = jetbtagN + 2;
     if(sel_1btag->passes(event)) jetbtagN = jetbtagN + 1;
  }

  if(TopTags.size() == 1 && WTags.size() != 1){
     for(auto & topjet : *event.topjets){
        auto subjets = topjet.subjets();
        for (auto & subjet : subjets) {
          if(sel_1btag->passes(event)) btag_subjet++;
        }


        for(unsigned int j=0; j<event.jets->size(); j++){
           double dR_jet = deltaR(topjet, event.jets->at(j));
           if(dR_jet > 1.2) continue;
           if(sel_1btag->passes(event)) jetbtagN++;
        }
     }
  }

  btagN = jetbtagN + btag_subjet;
  event.set(h_btagN,btagN);


  //Fill some vars for monitoring and comparison
  bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2); 
  if(is_zprime_reconstructed_chi2){
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
    //    float chi2 = BestZprimeCandidate->discriminator("chi2_total");
    event.set(h_chi2,BestZprimeCandidate->discriminator("chi2_total"));
    event.set(h_Mttbar,BestZprimeCandidate->Zprime_v4().M());

    if (isMuon){
        if(event.muons->at(0).charge() == 1){
                event.set(h_DeltaY,TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()));
        }
        if(event.muons->at(0).charge() == -1){
                event.set(h_DeltaY,TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()));
        }

        event.set(h_DeltaPhi_thad_lepton,deltaPhi(event.muons->at(0),BestZprimeCandidate->top_hadronic_v4()));
    }

    if (isElectron){
        if(event.electrons->at(0).charge() == 1){
                event.set(h_DeltaY,TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()));
        }
        if(event.electrons->at(0).charge() == -1){
                event.set(h_DeltaY,TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()));
        }

        event.set(h_DeltaPhi_thad_lepton,deltaPhi(event.electrons->at(0),BestZprimeCandidate->top_hadronic_v4()));
    }


    event.set(h_pTttbar,BestZprimeCandidate->Zprime_v4().Pt());
    event.set(h_Rapidityttbar,BestZprimeCandidate->Zprime_v4().Rapidity());
    event.set(h_pThad,BestZprimeCandidate->top_hadronic_v4().Pt());
    event.set(h_pTlep,BestZprimeCandidate->top_leptonic_v4().Pt());
    event.set(h_Rapidityhad,BestZprimeCandidate->top_hadronic_v4().Rapidity());
    event.set(h_Rapiditylep,BestZprimeCandidate->top_leptonic_v4().Rapidity());

    double mag_lep = TMath::Sqrt(BestZprimeCandidate->top_leptonic_v4().Px()*BestZprimeCandidate->top_leptonic_v4().Px() + BestZprimeCandidate->top_leptonic_v4().Py()*BestZprimeCandidate->top_leptonic_v4().Py() + BestZprimeCandidate->top_leptonic_v4().Pz()*BestZprimeCandidate->top_leptonic_v4().Pz());
    double mag_had = TMath::Sqrt(BestZprimeCandidate->top_hadronic_v4().Px()*BestZprimeCandidate->top_hadronic_v4().Px() + BestZprimeCandidate->top_hadronic_v4().Py()*BestZprimeCandidate->top_hadronic_v4().Py() + BestZprimeCandidate->top_hadronic_v4().Pz()*BestZprimeCandidate->top_hadronic_v4().Pz());
    double dot_prod = BestZprimeCandidate->top_hadronic_v4().Px()*BestZprimeCandidate->top_leptonic_v4().Px() + BestZprimeCandidate->top_hadronic_v4().Py()*BestZprimeCandidate->top_leptonic_v4().Py() + BestZprimeCandidate->top_hadronic_v4().Pz()*BestZprimeCandidate->top_leptonic_v4().Pz();
    event.set(h_Cos_theta_thad_tlep,dot_prod/(mag_had*mag_lep));

    if( BestZprimeCandidate->top_leptonic_v4().Pt() > BestZprimeCandidate->top_hadronic_v4().Pt()){
        event.set(h_tlead_pT,BestZprimeCandidate->top_leptonic_v4().Pt());
        event.set(h_tsublead_pT,BestZprimeCandidate->top_hadronic_v4().Pt());
    }else{
        event.set(h_tlead_pT,BestZprimeCandidate->top_hadronic_v4().Pt());
        event.set(h_tsublead_pT,BestZprimeCandidate->top_leptonic_v4().Pt());
    }

    event.set(h_TH_M,BestZprimeCandidate->top_hadronic_v4().M());
    event.set(h_TL_M,BestZprimeCandidate->top_leptonic_v4().M());
  }
  if(debug) cout<<"Set ttbar reconstruction vars for monitoring"<<endl;

  event.set(h_weight,event.weight);
  event.set(h_MET,event.met->pt());
  if(isMuon){
    event.set(h_lep1_pt,event.muons->at(0).pt());
    event.set(h_lep1_eta,event.muons->at(0).eta());
  }
  if(isElectron){
    event.set(h_lep1_pt,event.electrons->at(0).pt());
    event.set(h_lep1_eta,event.electrons->at(0).eta());
  }
  if(event.jets->size()>0){
    event.set(h_ak4jet1_pt,event.jets->at(0).pt());
    event.set(h_ak4jet1_eta,event.jets->at(0).eta());
  }
  if(event.topjets->size()>0){
    event.set(h_ak8jet1_pt,event.topjets->at(0).pt());
    event.set(h_ak8jet1_eta,event.topjets->at(0).eta());
  }
  event.set(h_NPV,event.pvs->size());
  if(debug) cout<<"Set some vars for monitoring"<<endl;

  const auto & sys_weights = event.genInfo->systweights();
  const auto & weights = event.genInfo->weights();

  if(isMC){
      if(!sys_weights.empty()){
             event.set(h_ISRup, weights[6]/event.genInfo->originalXWGTUP() );
             event.set(h_FSRup, weights[7]/event.genInfo->originalXWGTUP() );
             event.set(h_ISRdown, weights[8]/event.genInfo->originalXWGTUP() );
             event.set(h_FSRdown, weights[9]/event.genInfo->originalXWGTUP() );
      }
  }

  if(isMC){
      if(!sys_weights.empty()){ 
          float orig_weight = event.genInfo->originalXWGTUP();
          int MY_FIRST_INDEX = 9;
          for (unsigned i=0; i < 100; ++i) {
                const float pdf_w(sys_weights[i+MY_FIRST_INDEX]/orig_weight);
                w_PDF.push_back(pdf_w);
          }
       }
       event.set(h_wgtMC__PDF, std::move(w_PDF));
   }
   if(!isMC){
          w_PDF.push_back(1.);
          event.set(h_wgtMC__PDF, std::move(w_PDF));
          event.set(h_ISRup, 1.);
          event.set(h_FSRup, 1.);
          event.set(h_ISRdown, 1.);
          event.set(h_FSRdown, 1.); 
   }
  

  if(isMC){
      if(isMuon){ 
         float Tot_P = event.muons->at(0).pt()*cosh(event.muons->at(0).eta());
         if(is2018){
             if(-1.6 < event.muons->at(0).eta() && event.muons->at(0).eta() < 1.6){ 
                 if( 50 < Tot_P && Tot_P <= 100){ event.set(h_muonrecSF_nominal, 0.9943); event.set(h_muonrecSF_up, 0.995); event.set(h_muonrecSF_down, 0.9936); }
                 if( 100 < Tot_P && Tot_P <= 150){ event.set(h_muonrecSF_nominal, 0.9948); event.set(h_muonrecSF_up, 0.9955); event.set(h_muonrecSF_down, 0.9941); }
                 if( 150 < Tot_P && Tot_P <= 200){ event.set(h_muonrecSF_nominal, 0.9950); event.set(h_muonrecSF_up, 0.9959); event.set(h_muonrecSF_down, 0.9941); }
                 if( 200 < Tot_P && Tot_P <= 300){ event.set(h_muonrecSF_nominal, 0.994); event.set(h_muonrecSF_up, 0.9941); event.set(h_muonrecSF_down, 0.9939); }
                 if( 300 < Tot_P && Tot_P <= 400){ event.set(h_muonrecSF_nominal, 0.9914); event.set(h_muonrecSF_up, 0.9923); event.set(h_muonrecSF_down, 0.9905); }
                 if( 400 < Tot_P && Tot_P <= 600){ event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.9932); event.set(h_muonrecSF_down, 0.9928); }
                 if( 600 < Tot_P && Tot_P <= 1500){ event.set(h_muonrecSF_nominal, 0.9914); event.set(h_muonrecSF_up, 0.9914); event.set(h_muonrecSF_down, 0.9906); }
                 if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 1.0); event.set(h_muonrecSF_up, 1.1); event.set(h_muonrecSF_down, 0.9); }
             }

             if( (-1.6 > event.muons->at(0).eta() && event.muons->at(0).eta() > -2.4) || (1.6 < event.muons->at(0).eta() && event.muons->at(0).eta() < 2.4) ){
                 if( 100 < Tot_P && Tot_P <= 150){ event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.995); event.set(h_muonrecSF_down, 0.994); }
                 if( 150 < Tot_P && Tot_P <= 200){ event.set(h_muonrecSF_nominal, 0.990); event.set(h_muonrecSF_up, 0.991); event.set(h_muonrecSF_down, 0.989); }
                 if( 200 < Tot_P && Tot_P <= 300){ event.set(h_muonrecSF_nominal, 0.988); event.set(h_muonrecSF_up, 0.989); event.set(h_muonrecSF_down, 0.987); }
                 if( 300 < Tot_P && Tot_P <= 400){ event.set(h_muonrecSF_nominal, 0.981); event.set(h_muonrecSF_up, 0.983); event.set(h_muonrecSF_down, 0.979); }
                 if( 400 < Tot_P && Tot_P <= 600){ event.set(h_muonrecSF_nominal, 0.983); event.set(h_muonrecSF_up, 0.986); event.set(h_muonrecSF_down, 0.980); }
                 if( 600 < Tot_P && Tot_P <= 1500){ event.set(h_muonrecSF_nominal, 0.978); event.set(h_muonrecSF_up, 0.984); event.set(h_muonrecSF_down, 0.972);}
                 if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 0.98); event.set(h_muonrecSF_up, 1.01); event.set(h_muonrecSF_down, 0.95);  }
             }
         }
         if(is2017v2){
             if(-1.6 < event.muons->at(0).eta() && event.muons->at(0).eta() < 1.6){
                 if( 50 < Tot_P && Tot_P <= 100){ event.set(h_muonrecSF_nominal, 0.9938); event.set(h_muonrecSF_up, 0.9944); event.set(h_muonrecSF_down, 0.9932); }
                 if( 100 < Tot_P && Tot_P <= 150){ event.set(h_muonrecSF_nominal, 0.995); event.set(h_muonrecSF_up, 0.9957); event.set(h_muonrecSF_down, 0.9943); }
                 if( 150 < Tot_P && Tot_P <= 200){ event.set(h_muonrecSF_nominal, 0.996); event.set(h_muonrecSF_up, 0.997); event.set(h_muonrecSF_down, 0.995); }
                 if( 200 < Tot_P && Tot_P <= 300){ event.set(h_muonrecSF_nominal, 0.996); event.set(h_muonrecSF_up, 0.997); event.set(h_muonrecSF_down, 0.995); }
                 if( 300 < Tot_P && Tot_P <= 400){ event.set(h_muonrecSF_nominal, 0.994); event.set(h_muonrecSF_up, 0.995); event.set(h_muonrecSF_down, 0.993); }
                 if( 400 < Tot_P && Tot_P <= 600){ event.set(h_muonrecSF_nominal, 1.003); event.set(h_muonrecSF_up, 1.009); event.set(h_muonrecSF_down, 0.997); }
                 if( 600 < Tot_P && Tot_P <= 1500){ event.set(h_muonrecSF_nominal, 0.987); event.set(h_muonrecSF_up, 0.984); event.set(h_muonrecSF_down, 0.99); }
                 if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 0.9); event.set(h_muonrecSF_up, 1.); event.set(h_muonrecSF_down, 0.8); }
             }

             if( (-1.6 > event.muons->at(0).eta() && event.muons->at(0).eta() > -2.4) || (1.6 < event.muons->at(0).eta() && event.muons->at(0).eta() < 2.4) ){
                 if( 100 < Tot_P && Tot_P <= 150){ event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.994); event.set(h_muonrecSF_down, 0.992); }
                 if( 150 < Tot_P && Tot_P <= 200){ event.set(h_muonrecSF_nominal, 0.989); event.set(h_muonrecSF_up, 0.990); event.set(h_muonrecSF_down, 0.991); }
                 if( 200 < Tot_P && Tot_P <= 300){ event.set(h_muonrecSF_nominal, 0.986); event.set(h_muonrecSF_up, 0.987); event.set(h_muonrecSF_down, 0.985); }
                 if( 300 < Tot_P && Tot_P <= 400){ event.set(h_muonrecSF_nominal, 0.989); event.set(h_muonrecSF_up, 0.990); event.set(h_muonrecSF_down, 0.988); }
                 if( 400 < Tot_P && Tot_P <= 600){ event.set(h_muonrecSF_nominal, 0.983); event.set(h_muonrecSF_up, 0.986); event.set(h_muonrecSF_down, 0.980); }
                 if( 600 < Tot_P && Tot_P <= 1500){ event.set(h_muonrecSF_nominal, 0.986); event.set(h_muonrecSF_up, 0.992); event.set(h_muonrecSF_down, 0.980);}
                 if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 1.01); event.set(h_muonrecSF_up, 1.02); event.set(h_muonrecSF_down, 1.00);  }
             }
         }
         if(is2016v3 || is2016v2){
             if(-1.6 < event.muons->at(0).eta() && event.muons->at(0).eta() < 1.6){
                 if( 50 < Tot_P && Tot_P <= 100){ event.set(h_muonrecSF_nominal, 0.9914); event.set(h_muonrecSF_up, 0.9922); event.set(h_muonrecSF_down, 0.9906); }
                 if( 100 < Tot_P && Tot_P <= 150){ event.set(h_muonrecSF_nominal, 0.9936); event.set(h_muonrecSF_up, 0.9945); event.set(h_muonrecSF_down, 0.9927); }
                 if( 150 < Tot_P && Tot_P <= 200){ event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.994); event.set(h_muonrecSF_down, 0.992); }
                 if( 200 < Tot_P && Tot_P <= 300){ event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.995); event.set(h_muonrecSF_down, 0.991); }
                 if( 300 < Tot_P && Tot_P <= 400){ event.set(h_muonrecSF_nominal, 0.990); event.set(h_muonrecSF_up, 0.994); event.set(h_muonrecSF_down, 0.986); }
                 if( 400 < Tot_P && Tot_P <= 600){ event.set(h_muonrecSF_nominal, 0.990); event.set(h_muonrecSF_up, 0.993); event.set(h_muonrecSF_down, 0.987); }
                 if( 600 < Tot_P && Tot_P <= 1500){ event.set(h_muonrecSF_nominal, 0.989); event.set(h_muonrecSF_up, 0.985); event.set(h_muonrecSF_down, 0.993); }
                 if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 0.8); event.set(h_muonrecSF_up, 0.83); event.set(h_muonrecSF_down, 0.77); }
             }

             if( (-1.6 > event.muons->at(0).eta() && event.muons->at(0).eta() > -2.4) || (1.6 < event.muons->at(0).eta() && event.muons->at(0).eta() < 2.4) ){
                 if( 100 < Tot_P && Tot_P <= 150){ event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.994); event.set(h_muonrecSF_down, 0.992); }
                 if( 150 < Tot_P && Tot_P <= 200){ event.set(h_muonrecSF_nominal, 0.991); event.set(h_muonrecSF_up, 0.992); event.set(h_muonrecSF_down, 0.990); }
                 if( 200 < Tot_P && Tot_P <= 300){ event.set(h_muonrecSF_nominal, 0.985); event.set(h_muonrecSF_up, 0.986); event.set(h_muonrecSF_down, 0.984); }
                 if( 300 < Tot_P && Tot_P <= 400){ event.set(h_muonrecSF_nominal, 0.981); event.set(h_muonrecSF_up, 0.983); event.set(h_muonrecSF_down, 0.979); }
                 if( 400 < Tot_P && Tot_P <= 600){ event.set(h_muonrecSF_nominal, 0.979); event.set(h_muonrecSF_up, 0.983); event.set(h_muonrecSF_down, 0.975); }
                 if( 600 < Tot_P && Tot_P <= 1500){ event.set(h_muonrecSF_nominal, 0.978); event.set(h_muonrecSF_up, 0.983); event.set(h_muonrecSF_down, 0.973);}
                 if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 0.9); event.set(h_muonrecSF_up, 0.7); event.set(h_muonrecSF_down, 1.1);  }
             }
           

         }
      }
  }

 
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule)
