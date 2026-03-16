#include <iostream>
#include <limits>
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
#include <UHH2/common/include/JetCorrectionSets.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include "UHH2/common/include/LuminosityHists.h"
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/CommonModules.h>
// #include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicPreselectionHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/CHSJetCorrections.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopPuppiJetCorrections.h>
#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"

using namespace std;
using namespace uhh2;

class ZprimePreselectionModule : public ModuleBASE {

public:
  explicit ZprimePreselectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:
  bool debug;
  
  // mttbar mass bin edges
  const std::vector<double> mttbar_bin_edges = {0., 350., 500., 750., 1000., 1500., 20000.};
  // Bins: [0-350), [350-500), [500-750), [750-1000), [1000-1500), [1500-20000)
  
  // Helper function to find mttbar bin
  inline int find_mtt_bin(double mtt) {
    for (size_t i = 0; i+1 < mttbar_bin_edges.size(); ++i) {
      if (mtt >= mttbar_bin_edges[i] && mtt < mttbar_bin_edges[i+1]) return static_cast<int>(i);
    }
    return -1;
  }

  // Corrections
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<AnalysisModule> hotvrjetCorr;
  std::unique_ptr<TopPuppiJetCorrections> toppuppijetCorr;
  std::unique_ptr<CHSJetCorrections> CHSjetCorr;

  // Cleaners
  std::unique_ptr<JetCleaner>      jet_IDcleaner, jet_cleaner1, jet_cleaner2;
  std::unique_ptr<AnalysisModule>  hotvrjet_cleaner;
  std::unique_ptr<TopJetCleaner>   topjet_puppi_IDcleaner, topjet_puppi_cleaner;

  // Selections
  std::unique_ptr<uhh2::Selection> genflavor_sel;
  std::unique_ptr<uhh2::Selection> jet1_sel;
  std::unique_ptr<uhh2::Selection> jet2_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  unique_ptr<Selection> SignSplit;

  bool isMC, isHOTVR;
  string Sys_PU;

  std::unique_ptr<Hists> lumihists;
  TString METcollection;

  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;

  // additional branch with AK4 CHS jets -> for b-tagging
  Event::Handle<vector<Jet>> h_CHSjets;
  
  // TTbarGen handle for mttbar calculation
  Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<TTbarGenProducer> ttgenprod;
  uhh2::Event::Handle<float> h_xi_gen;
  uhh2::Event::Handle<float> h_mtt_gen;
  uhh2::Event::Handle<float> h_DeltaY_gen;

};

void ZprimePreselectionModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag+"_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicPreselectionHists(ctx,mytag));
  }
}

void ZprimePreselectionModule::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag+"_General";
  HFolder(mytag)->fill(event);
}

ZprimePreselectionModule::ZprimePreselectionModule(uhh2::Context& ctx) {

  debug = false; // true/false

  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }

  //// CONFIGURATION
  const TString METcollection = ctx.get("METName");
  isMC    = ctx.get("dataset_type") == "MC";
  isHOTVR = ctx.get("is_HOTVR") == "true";
  Sys_PU  = ctx.get("Sys_PU");

  isUL16preVFP  = (ctx.get("dataset_version").find("UL16preVFP")  != std::string::npos);
  isUL16postVFP = (ctx.get("dataset_version").find("UL16postVFP") != std::string::npos);
  isUL17        = (ctx.get("dataset_version").find("UL17")        != std::string::npos);
  isUL18        = (ctx.get("dataset_version").find("UL18")        != std::string::npos);
  
  // lepton IDs
  ElectronId eleID_veto = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90);
  MuonId     muID_veto  = MuonID(Muon::CutBasedIdTight);

  double electron_pt(25.);
  double muon_pt(25.);
  double jet1_pt(30.);
  double jet2_pt(30.);
  double MET(20.);

  // GEN Flavor selection [W+jets flavor-splitting]
  if(ctx.get("dataset_version").find("WJets") != std::string::npos){
    if     (ctx.get("dataset_version").find("_B") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("b"));
    else if(ctx.get("dataset_version").find("_C") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("c"));
    else if(ctx.get("dataset_version").find("_L") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("l"));
    else genflavor_sel.reset(new uhh2::AndSelection(ctx));
  }
  else genflavor_sel.reset(new uhh2::AndSelection(ctx));

  // Cleaning: Mu, Ele, Jets
  const MuonId muonID_veto(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID_veto));
  const ElectronId electronID_veto(AndId<Electron>(PtEtaSCCut(electron_pt, 2.5), eleID_veto));
  const JetPFID jetID_CHS(JetPFID::WP_TIGHT_CHS);
  const JetPFID jetID_PUPPI(JetPFID::WP_TIGHT_PUPPI);

  jet_IDcleaner.reset(new JetCleaner(ctx, jetID_PUPPI));
  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 20., 2.5));
  hotvrjet_cleaner.reset(new TopJetCleaner(ctx, PtEtaCut(200., 2.5)));
  topjet_puppi_IDcleaner.reset(new TopJetCleaner(ctx, jetID_PUPPI, "toppuppijets"));
  topjet_puppi_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(200., 2.5)), "toppuppijets"));

  // Split interference signal samples by sign
  if(ctx.get("dataset_version").find("_int") != std::string::npos){
    if     (ctx.get("dataset_version").find("_pos") != std::string::npos) SignSplit.reset(new SignSelection("pos"));
    else if(ctx.get("dataset_version").find("_neg") != std::string::npos) SignSplit.reset(new SignSelection("neg"));
    else SignSplit.reset(new uhh2::AndSelection(ctx));
  }
  else SignSplit.reset(new uhh2::AndSelection(ctx));

  // common modules
  common.reset(new CommonModules());
  common->switch_jetlepcleaner(true);
  // common->disable_pvfilter();
  common->disable_jetpfidfilter();
  common->switch_jetPtSorter(true);
  common->switch_metcorrection(true);
  common->set_muon_id(muonID_veto);
  common->set_electron_id(electronID_veto);
  common->init(ctx, Sys_PU);

  hotvrjetCorr.reset(new HOTVRJetCorrectionModule(ctx));

  toppuppijetCorr.reset(new TopPuppiJetCorrections());
  toppuppijetCorr->init(ctx);

  CHSjetCorr.reset(new CHSJetCorrections());
  CHSjetCorr->init(ctx);

  // TTbarGen producer
  if(isMC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", true));

  //// EVENT SELECTION
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.5))));
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.5))));
  met_sel.reset(new METCut(MET, uhh2::infinity));

  // additional branch with Ak4 CHS jets
  h_CHSjets = ctx.get_handle<vector<Jet>>("jetsAk4CHS");
  
  // TTbarGen handle for mttbar calculation
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  // GEN-level outputs (so they exist in the event and output tree)
  if (isMC) {
    h_xi_gen     = ctx.declare_event_output<float>("xi_gen");
    h_mtt_gen    = ctx.declare_event_output<float>("mtt_gen");
    h_DeltaY_gen = ctx.declare_event_output<float>("DeltaY_gen");
  }

  // Book histograms
  vector<string> histogram_tags = {"Input", "mtt_gen_inclusive", "CommonModules", "HOTVRCorrections", "PUPPICorrections", "Lepton1", "JetID", "JetCleaner1", "JetCleaner2", "TopjetCleaner", "Jet1", "Jet2", "MET"};
    
  // Add mttbar bin tags
  for (size_t i = 0; i+1 < mttbar_bin_edges.size(); ++i) {
    const double low = mttbar_bin_edges[i];
    const double high = mttbar_bin_edges[i+1];
    const string bin_tag = "mtt_gen_" + to_string((int)low) + "_" + to_string((int)high);
    histogram_tags.push_back(bin_tag);
  }
  
  book_histograms(ctx, histogram_tags);

  lumihists.reset(new LuminosityHists(ctx, "lumi"));
}

bool ZprimePreselectionModule::process(uhh2::Event& event){
  
  // Process TTbarGen first
  // if (isMC && ttgenprod) {
  //   ttgenprod->process(event);
  // }
  
  if(debug) cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug) cout << " run.event: " << event.run << ". " << event.event << endl;

  if(debug) cout << " event.year: " << event.year << ". " << event.event << endl;
 
  if(!event.isRealData){
    if(!SignSplit->passes(event)) return false;
  }
  if(debug) cout << "beginning: ok" << endl;

  fill_histograms(event, "Input");
  if(debug) cout << "first plots input: ok" << endl;

  // // Calculate mttbar and fill appropriate bin histograms
  // if (isMC && event.is_valid(h_ttbargen)) {
  //   // set defaults first, every event
  //   event.set(h_xi_gen,     std::numeric_limits<float>::quiet_NaN());
  //   event.set(h_mtt_gen,    std::numeric_limits<float>::quiet_NaN());
  //   event.set(h_DeltaY_gen, std::numeric_limits<float>::quiet_NaN());
  //   const auto& ttbargen = event.get(h_ttbargen);
  //   if (ttbargen.IsSemiLeptonicDecay()) {
  //     int lepId = std::abs(ttbargen.ChargedLepton().pdgId());
  //     if (lepId == 11 || lepId == 13) { 
  //       const auto& top  = ttbargen.Top();
  //       const auto& atop = ttbargen.Antitop();
  //       double mtt = (top.v4() + atop.v4()).M();
  //       double dy  = std::abs(top.v4().Rapidity()) - std::abs(atop.v4().Rapidity());
  //       event.set(h_xi_gen,     std::tanh(dy));
  //       event.set(h_mtt_gen,    static_cast<float>(mtt));
  //       event.set(h_DeltaY_gen, static_cast<float>(dy));

  //       // Only fill histograms if e/muon semileptonic
  //       fill_histograms(event, "mtt_gen_inclusive");
  //       const int ibin = find_mtt_bin(mtt);
  //       if (ibin >= 0) {
  //         const string bin_tag = "mtt_gen_" + to_string((int)mttbar_bin_edges[ibin]) + "_" + to_string((int)mttbar_bin_edges[ibin+1]);
  //         fill_histograms(event, bin_tag);
  //       }
  //     }
  //   }
  // }

  bool commonResult = common->process(event);
  if (!commonResult) return false;
  if(debug) cout << "CommonModules: ok" << endl;
  // fill_histograms(event, "CommonModules");

  sort_by_pt<Muon>(*event.muons);
  sort_by_pt<Electron>(*event.electrons);

  // Correct AK4 CHS jets
  CHSjetCorr->process(event);

  if(isHOTVR){
    hotvrjetCorr->process(event);
    // fill_histograms(event, "HOTVRCorrections");
  }

  toppuppijetCorr->process(event);
  if(debug) cout << "TopPuppiJetCorrections: ok" << endl;
  // fill_histograms(event, "PUPPICorrections");

  // GEN ME quark-flavor selection
  if(!event.isRealData){
    if(!genflavor_sel->passes(event)) return false;
  }
  if(debug) cout << "GenFlavorSelection: ok" << endl;

  // cout << "event.muons->size(): " << event.muons->size() << endl;
  // cout << "event.electrons->size(): " << event.electrons->size() << endl;
  const bool pass_lep1 = ((event.muons->size() >= 1) || (event.electrons->size() >= 1));
  // cout << "pass_lep1: " << pass_lep1 << endl;
  if(!pass_lep1) return false;
  if(debug) cout << "≥1 leptons: ok" << endl;
  // fill_histograms(event, "Lepton1");

  jet_IDcleaner->process(event);
  // fill_histograms(event, "JetID");
  if(debug) cout << "JetCleaner ID: ok" << endl;

  jet_cleaner1->process(event);
  sort_by_pt<Jet>(*event.jets);
  // fill_histograms(event, "JetCleaner1");
  if(debug) cout << "JetCleaner1: ok" << endl;

  // Lepton-2Dcut variables
  for(auto& muo : *event.muons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
    muo.set_tag(Muon::twodcut_dRmin, dRmin);
    muo.set_tag(Muon::twodcut_pTrel, pTrel);
  }

  for(auto& ele : *event.electrons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);
    ele.set_tag(Electron::twodcut_dRmin, dRmin);
    ele.set_tag(Electron::twodcut_pTrel, pTrel);
  }

  jet_cleaner2->process(event);
  sort_by_pt<Jet>(*event.jets);
  // fill_histograms(event, "JetCleaner2");
  if(debug) cout << "JetCleaner2: ok" << endl;

  hotvrjet_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);

  topjet_puppi_IDcleaner->process(event);
  topjet_puppi_cleaner->process(event);
  sort_by_pt<TopJet>(*event.toppuppijets);

  // fill_histograms(event, "TopjetCleaner");
  if(debug) cout << "TopJetCleaner: ok" << endl;

  // 1st AK4 jet selection
  const bool pass_jet1 = jet1_sel->passes(event);
  if(!pass_jet1) return false;
  if(debug) cout << "NJetSelection1: ok" << endl;
  // fill_histograms(event, "Jet1");

  // 2nd AK4 jet selection
  const bool pass_jet2 = jet2_sel->passes(event);
  if(!pass_jet2) return false;
  if(debug) cout << "NJetSelection2: ok" << endl;
  // fill_histograms(event, "Jet2");

  // MET selection
  const bool pass_met = met_sel->passes(event);
  if(!pass_met) return false;
  if(debug) cout << "METCut: ok" << endl;
  fill_histograms(event, "MET");


  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimePreselectionModule)