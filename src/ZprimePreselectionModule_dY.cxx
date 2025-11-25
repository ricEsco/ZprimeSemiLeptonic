#include <iostream>
#include <memory>
#include <cmath>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include "UHH2/common/include/PrintingModules.h"
#include <UHH2/core/include/Hists.h>

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/core/include/Utils.h>
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

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicPreselectionHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/CHSJetCorrections.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopPuppiJetCorrections.h>
#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"

#include <TMath.h>

using namespace std;
using namespace uhh2;

//dY_gen LINES
namespace {
  // 6 gen-level mttbar bins + overflow edge
  constexpr double mt_edges[7] = {0.0, 350.0, 500.0, 750.0, 1000.0, 1500.0, 1e9};
  inline int mt_bin(double mtt){
    for(int i=0; i<6; ++i){
      if(mtt >= mt_edges[i] && mtt < mt_edges[i+1]) return i;
    }
    return -1;
  }
  // Returns +1 for e+/μ+ (PDG -11/-13), −1 for e-/μ- (PDG +11/+13)
  inline int lepton_charge_from_pdg(int pdg){
    if(pdg == -11 || pdg == -13) return +1; // positive lepton
    if(pdg ==  11 || pdg ==  13) return -1; // negative lepton
    return 0;
  }
}

// === Two-bin delta|y| histogram for gen-level counting (x-axis is Δ|y| from -2.5..+2.5) ===
class DYGenHists : public uhh2::Hists {
public:
  explicit DYGenHists(uhh2::Context & ctx, const std::string & dirname)
    : Hists(ctx, dirname)
  {
    // Event counter (sum of weights)
    h_events = book<TH1F>("Events", "Events", 1, 0.5, 1.5);

    // Two-bin delta|y| hist with desired axis range
    // Bins: [-2.5, 0) → "neg" ; [0, 2.5] → "pos" (ROOT puts exactly 0 into the second bin)
    h_dy = book<TH1F>("dYgen", "gen #Delta|y|;#Delta|y|;events", 2, -2.5, 2.5);
    h_dy->GetXaxis()->SetBinLabel(1, "neg");
    h_dy->GetXaxis()->SetBinLabel(2, "pos");
  }

  virtual void fill(const uhh2::Event & event) override {
    if(!event.is_valid(h_DeltaY_gen_handle)) return;
    const double w  = event.weight;
    const double dY = event.get(h_DeltaY_gen_handle);

    h_events->Fill(1, w);

    // Clamp to histogram range so no event goes to overflow (we only care about sign)
    double dYc = dY;
    if (dYc < -2.5) dYc = -2.5;
    else if (dYc > 2.5) dYc = 2.5;

    h_dy->Fill(dYc, w);
  }

  // Handle to read delta|y| from the event (set by the module)
  static uhh2::Event::Handle<float> h_DeltaY_gen_handle;

private:
  TH1F *h_events{nullptr};
  TH1F *h_dy{nullptr};
};
uhh2::Event::Handle<float> DYGenHists::h_DeltaY_gen_handle;

class ZprimePreselectionModule_dY : public ModuleBASE {
public:
  explicit ZprimePreselectionModule_dY(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, const vector<string>&);
  void fill_histograms(uhh2::Event&, const string&);

protected:
  bool debug;

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
  std::unique_ptr<uhh2::Selection> SignSplit;

  // Gen outputs
  Event::Handle<float> h_DeltaY_gen;
  Event::Handle<float> h_mttbar_gen;

  bool isMC, isHOTVR;
  string Sys_PU;

  std::unique_ptr<Hists> lumihists;
  TString METcollection;

  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;

  // additional branch with AK4 CHS jets -> for b-tagging
  Event::Handle<vector<Jet>> h_CHSjets;
};

// void ZprimePreselectionModule_dY::book_histograms(uhh2::Context& ctx, const vector<string>& tags){
//   for(const auto & tag : tags){
//     const string mytag = tag + "_General";
//     book_HFolder(mytag, new ZprimeSemiLeptonicPreselectionHists(ctx, mytag));
//   }
// }

void ZprimePreselectionModule_dY::book_histograms(uhh2::Context& ctx, const vector<string>& tags){
  for (const auto & tag : tags){
    const string folder = tag + "_General";
    const bool isDY =
      (tag.rfind("DY_", 0) == 0)    ||
      (tag.rfind("SL_DY_", 0) == 0) ||
      (tag.rfind("DL_DY_", 0) == 0);

    if (isDY) {
      book_HFolder(folder, new DYGenHists(ctx, folder));
    } else {
      // everything else uses your usual preselection hists
      book_HFolder(folder, new ZprimeSemiLeptonicPreselectionHists(ctx, folder));
    }
  }
}


void ZprimePreselectionModule_dY::fill_histograms(uhh2::Event& event, const string& tag){
  const string mytag = tag + "_General";
  HFolder(mytag)->fill(event);
}

ZprimePreselectionModule_dY::ZprimePreselectionModule_dY(uhh2::Context& ctx){

  debug = false;

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
  // ElectronId eleID_veto = ElectronID_Fall17_tight_noIso; - eleID_veto tag changed
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

  //// EVENT SELECTION
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.5))));
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.5))));
  met_sel.reset(new METCut(MET, uhh2::infinity));

  // additional branch with Ak4 CHS jets
  h_CHSjets = ctx.get_handle<vector<Jet>>("jetsAk4CHS");

  //dY_gen LINES
  // DeltaY & Mttbar outputs
  h_DeltaY_gen = ctx.declare_event_output<float>("DeltaY_gen");
  h_mttbar_gen = ctx.declare_event_output<float>("mttbar_gen");
  DYGenHists::h_DeltaY_gen_handle = h_DeltaY_gen; // pass handle to hist class

  // Book usual diagnostic folders
  vector<string> histogram_tags = {"Input",
// semileptonic Δ|y| by mtt bin
    "SL_DY_0_350", "SL_DY_350_500", "SL_DY_500_750",
    "SL_DY_750_1000", "SL_DY_1000_1500", "SL_DY_1500Inf",

    // dileptonic Δ|y| by mtt bin
    "DL_DY_0_350", "DL_DY_350_500", "DL_DY_500_750",
    "DL_DY_750_1000", "DL_DY_1000_1500", "DL_DY_1500Inf",   
    "CommonModules", "HOTVRCorrections", "PUPPICorrections", "Lepton1", "JetID", "JetCleaner1", "JetCleaner2", "TopjetCleaner", "Jet1", "Jet2", "MET"};
  book_histograms(ctx, histogram_tags);

  lumihists.reset(new LuminosityHists(ctx, "lumi"));
}

bool ZprimePreselectionModule_dY::process(uhh2::Event& event){

  if(debug) cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug) cout << " run.event: " << event.run << ". " << event.event << endl;

  //dY_gen LINES
  // Init outputs
  event.set(h_DeltaY_gen, -999.f);
  event.set(h_mttbar_gen, -999.f);
  //dY_gen LINES

  // Fill input folder (for bookkeeping)
  fill_histograms(event, "Input");

  // Interference sign split if applicable
  if(!event.isRealData){
    if(!SignSplit->passes(event)) return false;
  }

  // ===== GENERATOR LEVEL (BEFORE ANY CUTS) =====
  if(!event.isRealData){
    TTbarGen ttbargen(*event.genparticles);

    const auto ch = ttbargen.DecayChannel();

    const bool is_semilep = (ch == TTbarGen::e_ehad) || (ch == TTbarGen::e_muhad);
    const bool is_dilep   = (ch == TTbarGen::e_ee)   || (ch == TTbarGen::e_mumu) || (ch == TTbarGen::e_emu);

    // We fill delta|y| and mttbar for either semileptonic or dileptonic channels.
    if(is_semilep || is_dilep){
        const GenParticle top     = ttbargen.Top();
        const GenParticle antitop = ttbargen.Antitop();
        
        const double mtt = (top.v4() + antitop.v4()).M();
        event.set(h_mttbar_gen, float(mtt));

      double dY = 1e9;

      // semileptonic case
      if (is_semilep) {
        // semilep -sign logic (d|Y| = |y(t)| - |y(tbar)|)
        const bool top_is_lep     = !ttbargen.IsTopHadronicDecay();
        const bool antitop_is_lep = !ttbargen.IsAntiTopHadronicDecay();

        // one leptonic top and one hadronic top
        if ( (top_is_lep && !antitop_is_lep) || (!top_is_lep && antitop_is_lep) ) {

          // truth tops
          const GenParticle top     = ttbargen.Top();
          const GenParticle antitop = ttbargen.Antitop();

          // identify which truth top is the leptonic one
          GenParticle t_lep, t_had;
          if (top_is_lep) {
            t_lep = top;
            t_had = antitop;
          } else {
            t_lep = antitop;
            t_had = top;
          }

          // charged lepton and its charge (+1 for e+/mu+, -1 for e-/mu-)
          const GenParticle lep = ttbargen.ChargedLepton();
          const int pdg     = lep.pdgId();
          const int abs_pdg = std::abs(pdg);

          if (abs_pdg == 11 || abs_pdg == 13) {
            int ql = 0;
            if (pdg == -11 || pdg == -13) ql = +1;
            else if (pdg ==  11 || pdg ==  13) ql = -1;

            const double yL = std::abs(t_lep.v4().Rapidity());
            const double yH = std::abs(t_had.v4().Rapidity());

            if (ql > 0) {
              dY = yL - yH;   // positive lepton: |y_lep| - |y_had|
            } else {
              dY = yH - yL;   // negative lepton: |y_had| - |y_lep|
            }
          } else {
            dY = std::abs(top.v4().Rapidity()) - std::abs(antitop.v4().Rapidity());
          }
        }

        // only set if computed
        if (std::abs(dY) < 1e8) {
          event.set(h_DeltaY_gen, float(dY));
        }
      }

      // dilepton case
      else if(is_dilep){
        dY = std::abs(top.v4().Rapidity()) - std::abs(antitop.v4().Rapidity());
      }

      if (std::abs(dY) < 1e8) {
        event.set(h_DeltaY_gen, float(dY));

        if (is_semilep) {
          if      (mtt >= 0   && mtt < 350)   fill_histograms(event, "SL_DY_0_350");
          else if (mtt >= 350 && mtt < 500)   fill_histograms(event, "SL_DY_350_500");
          else if (mtt >= 500 && mtt < 750)   fill_histograms(event, "SL_DY_500_750");
          else if (mtt >= 750 && mtt < 1000)  fill_histograms(event, "SL_DY_750_1000");
          else if (mtt >= 1000 && mtt < 1500) fill_histograms(event, "SL_DY_1000_1500");
          else if (mtt >= 1500)               fill_histograms(event, "SL_DY_1500Inf");
        }
        else if (is_dilep) {
          if      (mtt >= 0   && mtt < 350)   fill_histograms(event, "DL_DY_0_350");
          else if (mtt >= 350 && mtt < 500)   fill_histograms(event, "DL_DY_350_500");
          else if (mtt >= 500 && mtt < 750)   fill_histograms(event, "DL_DY_500_750");
          else if (mtt >= 750 && mtt < 1000)  fill_histograms(event, "DL_DY_750_1000");
          else if (mtt >= 1000 && mtt < 1500) fill_histograms(event, "DL_DY_1000_1500");
          else if (mtt >= 1500)               fill_histograms(event, "DL_DY_1500Inf");
        }
      }

      }
    }
  //dY_gen LINES

  // ===== STANDARD PRESELECTION CUTS =====
  
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
  // fill_histograms(event, "MET");


  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimePreselectionModule_dY)