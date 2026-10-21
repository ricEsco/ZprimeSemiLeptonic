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
#include <UHH2/common/include/JetCorrectionSets.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/LuminosityHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/CommonModules.h>

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicPreselectionHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/EntanglementStrengthMapHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/CHSJetCorrections.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopPuppiJetCorrections.h>
#include <UHH2/HOTVR/include/HOTVRJetCorrectionModule.h>

using namespace std;
using namespace uhh2;

/*
Preselection: first step of the analysis, run on the inclusive UHH2 ntuples.
  1) Inclusive gen-level alpha strength maps (nominal TTToSemiLeptonic only), filled BEFORE any selection
  2) Object corrections (CommonModules, AK4 CHS, HOTVR, AK8 PUPPI) and cleaning
  3) Loose event selection: >=1 lepton, >=2 AK4 jets, MET
The output AnalysisTree is the input of ZprimeAnalysisModule.
*/

/*
██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class ZprimePreselectionModule : public ModuleBASE {

public:

  explicit ZprimePreselectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:

  bool debug; // enable debug flag


  ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Declare unique_ptr's /////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // Gen-level ttbar system and inclusive strength maps
  std::unique_ptr<TTbarGenProducer> ttgenprod;   // builds the TTbarGen object ("ttbargen" handle) from GenParticles
  std::unique_ptr<Hists> strength_maps;          // inclusive alpha(m_tt, cos theta*) maps, see EntanglementStrengthMapHists.h

  // Corrections
  std::unique_ptr<CommonModules> common;                    // lumi selection, MET filters, PV, lepton IDs, AK4 JEC/JER, MET correction, PU weights
  std::unique_ptr<AnalysisModule> hotvrjetCorr;             // HOTVR jet corrections
  std::unique_ptr<TopPuppiJetCorrections> toppuppijetCorr;  // AK8 PUPPI jet corrections
  std::unique_ptr<CHSJetCorrections> CHSjetCorr;            // AK4 CHS jet corrections (CHS jets are used for b-tagging)

  // Cleaners
  std::unique_ptr<JetCleaner>      jet_IDcleaner, jet_cleaner1, jet_cleaner2;   // AK4 PUPPI jet ID and pT/eta cleaners
  std::unique_ptr<AnalysisModule>  hotvrjet_cleaner;                            // HOTVR jet pT/eta cleaner
  std::unique_ptr<TopJetCleaner>   topjet_puppi_IDcleaner, topjet_puppi_cleaner; // AK8 PUPPI jet ID and pT/eta cleaners

  // Selections
  std::unique_ptr<uhh2::Selection> genflavor_sel; // W+jets gen-flavor splitting (b, c, light)
  std::unique_ptr<uhh2::Selection> jet1_sel;      // >=1 AK4 jet
  std::unique_ptr<uhh2::Selection> jet2_sel;      // >=2 AK4 jets
  std::unique_ptr<uhh2::Selection> met_sel;       // MET cut
  std::unique_ptr<uhh2::Selection> SignSplit;     // splits interference signal samples by weight sign

  std::unique_ptr<Hists> lumihists;


  ///////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Declare handles ///////////////////////////////////////

  Event::Handle<vector<Jet>> h_CHSjets; // additional branch with AK4 CHS jets (b-tagging)


  /////////////////////////////////////////////////////////////
  /////////////////////// Configuration ///////////////////////
  bool isMC, isHOTVR;
  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;
  bool fill_strength_maps; // true only for the nominal TTToSemiLeptonic sample of each era
  string Sys_PU;
};

void ZprimePreselectionModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag + "_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicPreselectionHists(ctx, mytag));
  }
}

void ZprimePreselectionModule::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag + "_General";
  HFolder(mytag)->fill(event);
}


/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

ZprimePreselectionModule::ZprimePreselectionModule(uhh2::Context& ctx){

  debug = false;

  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }

  /////////////////////////////////////////////////
  ///////////////// Configuration /////////////////

  isMC    = (ctx.get("dataset_type") == "MC");
  isHOTVR = (ctx.get("is_HOTVR") == "true");
  Sys_PU  = ctx.get("Sys_PU");

  const string dataset_version = ctx.get("dataset_version");
  isUL16preVFP  = (dataset_version.find("UL16preVFP")  != std::string::npos);
  isUL16postVFP = (dataset_version.find("UL16postVFP") != std::string::npos);
  isUL17        = (dataset_version.find("UL17")        != std::string::npos);
  isUL18        = (dataset_version.find("UL18")        != std::string::npos);

  string year = "";
  if(isUL16preVFP)  year = "UL16preVFP";
  if(isUL16postVFP) year = "UL16postVFP";
  if(isUL17)        year = "UL17";
  if(isUL18)        year = "UL18";

  // Strength maps are gen-level only: fill them once per era, on the nominal TTToSemiLeptonic sample
  // (skipped in the JEC/JER-varied preselection runs, which would only produce identical copies)
  const bool is_nominal_jets = (ctx.get("jecsmear_direction", "nominal") == "nominal") && (ctx.get("jersmear_direction", "nominal") == "nominal");
  fill_strength_maps = isMC && (dataset_version == "TTToSemiLeptonic_" + year) && is_nominal_jets;

  /// Important selection values
  double electron_pt(25.); // veto-lepton pT thresholds
  double muon_pt(25.);
  double jet1_pt(30.);     // pT threshold for leading AK4 jet
  double jet2_pt(30.);     // pT threshold for subleading AK4 jet
  double MET(20.);         // MET threshold


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////// GEN-LEVEL /////////////////////////////////////////////

  // TTbarGen: does not throw for non-ttbar samples (decay channel is then simply "not found")
  if(isMC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));

  // Inclusive alpha strength maps -> output folder "StrengthMaps"
  if(fill_strength_maps) strength_maps.reset(new EntanglementStrengthMapHists(ctx, "StrengthMaps"));

  // W+jets gen-flavor splitting
  if(dataset_version.find("WJets") != std::string::npos){
    if     (dataset_version.find("_B") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("b"));
    else if(dataset_version.find("_C") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("c"));
    else if(dataset_version.find("_L") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("l"));
    else genflavor_sel.reset(new uhh2::AndSelection(ctx));
  }
  else genflavor_sel.reset(new uhh2::AndSelection(ctx));

  // Split interference signal samples by sign
  if(dataset_version.find("_int") != std::string::npos){
    if     (dataset_version.find("_pos") != std::string::npos) SignSplit.reset(new SignSelection("pos"));
    else if(dataset_version.find("_neg") != std::string::npos) SignSplit.reset(new SignSelection("neg"));
    else SignSplit.reset(new uhh2::AndSelection(ctx));
  }
  else SignSplit.reset(new uhh2::AndSelection(ctx));


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// OBJECT IDs AND CLEANING //////////////////////////////////////

  // veto-lepton IDs
  ElectronId eleID_veto = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90);
  MuonId     muID_veto  = MuonID(Muon::CutBasedIdTight);
  const MuonId muonID_veto(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID_veto));
  const ElectronId electronID_veto(AndId<Electron>(PtEtaSCCut(electron_pt, 2.5), eleID_veto));

  // jet IDs
  const JetPFID jetID_CHS(JetPFID::WP_TIGHT_CHS);
  const JetPFID jetID_PUPPI(JetPFID::WP_TIGHT_PUPPI);

  // jet cleaners
  jet_IDcleaner.reset(new JetCleaner(ctx, jetID_PUPPI));
  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 20., 2.5));
  hotvrjet_cleaner.reset(new TopJetCleaner(ctx, PtEtaCut(200., 2.5)));
  topjet_puppi_IDcleaner.reset(new TopJetCleaner(ctx, jetID_PUPPI, "toppuppijets"));
  topjet_puppi_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(200., 2.5)), "toppuppijets"));


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// CORRECTIONS ////////////////////////////////////////////

  // common modules
  common.reset(new CommonModules());
  common->switch_jetlepcleaner(true);
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


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// EVENT SELECTION //////////////////////////////////////////

  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.5))));
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.5))));
  met_sel.reset(new METCut(MET, uhh2::infinity));

  // additional branch with AK4 CHS jets
  h_CHSjets = ctx.get_handle<vector<Jet>>("jetsAk4CHS");

  // Book histograms
  // dropped "HOTVRCorrections" since we don't use them
  vector<string> histogram_tags = {"Input", "CommonModules", "PUPPICorrections", "Lepton1", "JetID", "JetCleaner1", "JetCleaner2", "TopjetCleaner", "Jet1", "Jet2", "MET"};
  book_histograms(ctx, histogram_tags);

  lumihists.reset(new LuminosityHists(ctx, "lumi"));
}


/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool ZprimePreselectionModule::process(uhh2::Event& event){

  if(debug) cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug) cout << " run.event: " << event.run << ". " << event.event << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Gen-level: TTbarGen and inclusive strength maps /////////////////////////

  // Build the gen-level ttbar system
  if(isMC && ttgenprod) ttgenprod->process(event);
  if(debug) cout << "[ZprimePreselectionModule] TTbarGen: ok" << endl;

  // Inclusive strength maps: filled for every event BEFORE any selection
  if(fill_strength_maps) strength_maps->fill(event);
  if(debug) cout << "[ZprimePreselectionModule] strength maps: ok" << endl;

  // Split interference signal samples by sign
  if(!event.isRealData){
    if(!SignSplit->passes(event)) return false;
  }
  fill_histograms(event, "Input");


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// Corrections /////////////////////////////////////////////

  // CommonModules (can reject events, e.g. lumi selection and MET filters)
  if(!common->process(event)) return false;
  if(debug) cout << "[ZprimePreselectionModule] CommonModules: ok" << endl;
  // fill_histograms(event, "CommonModules");

  sort_by_pt<Muon>(*event.muons);
  sort_by_pt<Electron>(*event.electrons);

  // AK4 CHS jets
  CHSjetCorr->process(event);

  // HOTVR jets
  if(isHOTVR){
    hotvrjetCorr->process(event);
    // fill_histograms(event, "HOTVRCorrections");
  }

  // AK8 PUPPI jets
  toppuppijetCorr->process(event);
  if(debug) cout << "[ZprimePreselectionModule] TopPuppiJetCorrections: ok" << endl;
  // fill_histograms(event, "PUPPICorrections");


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// Selection and cleaning //////////////////////////////////////

  // W+jets gen-flavor selection
  if(!event.isRealData){
    if(!genflavor_sel->passes(event)) return false;
  }
  if(debug) cout << "[ZprimePreselectionModule] GenFlavorSelection: ok" << endl;

  // >=1 lepton (muon or electron)
  const bool pass_lep1 = ((event.muons->size() >= 1) || (event.electrons->size() >= 1));
  if(!pass_lep1) return false;
  if(debug) cout << "[ZprimePreselectionModule] >=1 lepton: ok" << endl;
  // fill_histograms(event, "Lepton1");

  // AK4 PUPPI jet ID
  jet_IDcleaner->process(event);
  if(debug) cout << "[ZprimePreselectionModule] JetCleaner ID: ok" << endl;
  // fill_histograms(event, "JetID");

  // AK4 PUPPI jets: pT > 15 GeV, |eta| < 3.0
  jet_cleaner1->process(event);
  sort_by_pt<Jet>(*event.jets);
  if(debug) cout << "[ZprimePreselectionModule] JetCleaner1: ok" << endl;
  // fill_histograms(event, "JetCleaner1");

  // Lepton 2D-cut variables (dRmin, pTrel to the closest jet), used by the 2D-cut in the Analysis
  for(auto& muo : *event.muons){
    float dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
    muo.set_tag(Muon::twodcut_dRmin, dRmin);
    muo.set_tag(Muon::twodcut_pTrel, pTrel);
  }
  for(auto& ele : *event.electrons){
    float dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);
    ele.set_tag(Electron::twodcut_dRmin, dRmin);
    ele.set_tag(Electron::twodcut_pTrel, pTrel);
  }

  // AK4 PUPPI jets: pT > 20 GeV, |eta| < 2.5
  jet_cleaner2->process(event);
  sort_by_pt<Jet>(*event.jets);
  if(debug) cout << "[ZprimePreselectionModule] JetCleaner2: ok" << endl;
  // fill_histograms(event, "JetCleaner2");

  // HOTVR jets: pT > 200 GeV, |eta| < 2.5
  hotvrjet_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);

  // AK8 PUPPI jets: ID and pT > 200 GeV, |eta| < 2.5
  topjet_puppi_IDcleaner->process(event);
  topjet_puppi_cleaner->process(event);
  sort_by_pt<TopJet>(*event.toppuppijets);
  if(debug) cout << "[ZprimePreselectionModule] TopJetCleaner: ok" << endl;
  // fill_histograms(event, "TopjetCleaner");

  // >=1 AK4 jet with pT > 30 GeV
  if(!jet1_sel->passes(event)) return false;
  if(debug) cout << "[ZprimePreselectionModule] NJetSelection1: ok" << endl;
  // fill_histograms(event, "Jet1");

  // >=2 AK4 jets with pT > 30 GeV
  if(!jet2_sel->passes(event)) return false;
  if(debug) cout << "[ZprimePreselectionModule] NJetSelection2: ok" << endl;
  // fill_histograms(event, "Jet2");

  // MET > 20 GeV
  if(!met_sel->passes(event)) return false;
  if(debug) cout << "[ZprimePreselectionModule] METCut: ok" << endl;
  fill_histograms(event, "MET");

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimePreselectionModule)
