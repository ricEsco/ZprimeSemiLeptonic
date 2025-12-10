#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSystematicsHists.h"
#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include "UHH2/common/include/JetIds.h"
#include <math.h>

#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <UHH2/core/include/LorentzVector.h>

#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <glob.h>
#include <cstring>

using namespace std;
using namespace uhh2;

//template method start
// STEP 1: Build the NoAC weights from the generator histogram

// --- NoAC helpers (mirror/S+A/lookup) ---
std::unique_ptr<TH1D> ZprimeSemiLeptonicSystematicsHists::mirror_hist_1d(const TH1D &src){
  auto H = std::unique_ptr<TH1D>(static_cast<TH1D*>(src.Clone("H_mirror_tmp")));
  H->SetDirectory(0);
  H->Reset("ICES");

  // mirror the histogram around 0
  const TAxis *ax = src.GetXaxis();
  const int nb = ax->GetNbins();
  for(int i=1;i<=nb;++i){
    const double xc = ax->GetBinCenter(i);
    int j = ax->FindBin(-xc);
    if(j < 1) j = 1;
    if(j > nb) j = nb;
    H->SetBinContent(i, src.GetBinContent(j));
    H->SetBinError(i, src.GetBinError(j));
  }
  return H;
}

// build the NoAC weights from the generator histogram
// STEP 3: Decompose the generator histogram into S and A components
std::unique_ptr<TH1D> ZprimeSemiLeptonicSystematicsHists::build_noac_weights_from_gen(const TH1D &Hgen_in, double f_noac){
  auto Hmir = mirror_hist_1d(Hgen_in);
  auto S = std::unique_ptr<TH1D>(static_cast<TH1D*>(Hgen_in.Clone("H_S")));
  auto A = std::unique_ptr<TH1D>(static_cast<TH1D*>(Hgen_in.Clone("H_A")));
  S->SetDirectory(0); A->SetDirectory(0);
  S->Reset(); A->Reset();
  S->Add(&Hgen_in, Hmir.get(), 0.5,  0.5); //0.5 is the weight for the original and mirrored histogram
  A->Add(&Hgen_in, Hmir.get(), 0.5, -0.5); //0.5 is the weight for the original and mirrored (-0.5)histogram

  // STEP 4: Build the NoAC weights from the S and A components
  auto W = std::unique_ptr<TH1D>(static_cast<TH1D*>(Hgen_in.Clone("NoAC_W_sys")));
  W->SetDirectory(0); // detach the histogram from the file
  W->Reset(); // reset the histogram
  const int nb = W->GetXaxis()->GetNbins();
  for(int i=1;i<=nb;++i){
    const double s = S->GetBinContent(i);
    const double a = A->GetBinContent(i);
    const double denom = s + a;
    const double numer = s + (1.0 - f_noac) * a; 
    const double w = denom > 0.0 ? (numer / denom) : 1.0; //if the denominator is greater than 0, divide the numerator by the denominator, otherwise set to 1.0
    W->SetBinContent(i, w); //set the content of the bin to the weight
  }
  
  // STEP 5: Normalize the NoAC weights to the sum of the content of the histogram
  // The fit in Combine assumes that the shape variations (NoAC_up and down) have the same total number of events as the nominal shape.
  // normalize <W>_Hgen = 1
  // CHECKED: Before normalization, W(f=+1) + W(f=-1) = 2 in each bin.
  // After normalization, ⟨W⟩ = 1 for each f, which is what Combine wants for pure shape variations.
  double sumW = 0.0, sumH = 0.0;
  for(int i=1;i<=nb;++i){
    const double wi = W->GetBinContent(i);
    const double hi = Hgen_in.GetBinContent(i);
    sumW += wi * hi; //sum of the weights * the content of the histogram. Calulcates the total weighted number of events. 
    sumH += hi; //sum of the content of the histogram
  }
  if(sumW > 0.0 && sumH > 0.0){
    const double norm = sumW / sumH; //normalize the weights to the sum of the content of the histogram
    for(int i=1;i<=nb;++i){
      W->SetBinContent(i, W->GetBinContent(i) / norm);
    }
  }
  return W;
}

// lookup the NoAC weight for a given xi
// Clamps xi to be inside the histogram range to avoid edge effects
// Clamps xi_gen_evt to [-0.999999, 0.999999] to avoid edge effects
double ZprimeSemiLeptonicSystematicsHists::lookup_noac_weight(const TH1 *W, double xi){
  if(!W || !std::isfinite(xi)) return 1.0;
  double xmin = W->GetXaxis()->GetXmin();
  double xmax = W->GetXaxis()->GetXmax();
  // Clamp xi to be inside the histogram range to avoid edge effects
  if(xi <= xmin) xi = std::nextafter(xmin, xmax);
  if(xi >= xmax) xi = std::nextafter(xmax, xmin);
  int bin = W->GetXaxis()->FindFixBin(xi);
  if(bin < 1) bin = 1;
  if(bin > W->GetNbinsX()) bin = W->GetNbinsX();
  const double w = W->GetBinContent(bin);
  if(!std::isfinite(w) || w <= 0.0 || w > 100.0) return 1.0;
  return w;
}

static inline double clamp_xi(double x){
  // keep inside histogram range (-1,1) to avoid edge effects
  if (x <= -1.0) return -0.999999;
  if (x >= +1.0) return +0.999999;
  return x;
}
//template method end

// Static member definitions - shared across all instances
std::map<float, std::unique_ptr<TH1D>> ZprimeSemiLeptonicSystematicsHists::noac_weights_map;
bool ZprimeSemiLeptonicSystematicsHists::noac_weights_initialized = false;

ZprimeSemiLeptonicSystematicsHists::ZprimeSemiLeptonicSystematicsHists(uhh2::Context& ctx, const std::string& dirname):
Hists(ctx, dirname) {

  is_mc = ctx.get("dataset_type") == "MC";
  is_Muon = ctx.get("channel") == "muon";
  std::string dataset_version = ctx.get("dataset_version");
  is_tt = (dataset_version.find("TTTo") == 0) || (dataset_version.find("EFT") != std::string::npos);

  isMuon = false; isElectron = false;
  if(ctx.get("channel") == "muon") isMuon = true;
  if(ctx.get("channel") == "electron") isElectron = true;
  ishotvr = (ctx.get("is_hotvr") == "true");
  isdeepAK8 = (ctx.get("is_deepAK8") == "true");
  if(isdeepAK8){
    h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("DeepAK8TopTags");
  }else if(ishotvr){
    h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("HOTVRTopTags");
  }
  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");

  // electron systematics: reco, id, trigger
  h_ele_reco           = ctx.get_handle<float>("weight_sfelec_reco");
  h_ele_reco_up        = ctx.get_handle<float>("weight_sfelec_reco_up");
  h_ele_reco_down      = ctx.get_handle<float>("weight_sfelec_reco_down");
  h_ele_id             = ctx.get_handle<float>("weight_sfelec_id");
  h_ele_id_up          = ctx.get_handle<float>("weight_sfelec_id_up");
  h_ele_id_down        = ctx.get_handle<float>("weight_sfelec_id_down");
  h_ele_trigger        = ctx.get_handle<float>("weight_sfelec_trigger");
  h_ele_trigger_up     = ctx.get_handle<float>("weight_sfelec_trigger_up");
  h_ele_trigger_down   = ctx.get_handle<float>("weight_sfelec_trigger_down");
  // muon systematics: reco, id_stat, id_syst, trigger_stat, trigger_syst, iso_stat, iso_syst
  h_mu_reco            = ctx.get_handle<float>("weight_sfmu_reco");
  h_mu_reco_up         = ctx.get_handle<float>("weight_sfmu_reco_up");
  h_mu_reco_down       = ctx.get_handle<float>("weight_sfmu_reco_down");
  h_mu_id_stat         = ctx.get_handle<float>("weight_sfmu_id_stat");
  h_mu_id_stat_up      = ctx.get_handle<float>("weight_sfmu_id_stat_up");
  h_mu_id_stat_down    = ctx.get_handle<float>("weight_sfmu_id_stat_down");
  h_mu_id_syst         = ctx.get_handle<float>("weight_sfmu_id_syst");
  h_mu_id_syst_up      = ctx.get_handle<float>("weight_sfmu_id_syst_up");
  h_mu_id_syst_down    = ctx.get_handle<float>("weight_sfmu_id_syst_down");
  h_mu_trigger_stat    = ctx.get_handle<float>("weight_sfmu_trigger_stat");
  h_mu_trigger_stat_up = ctx.get_handle<float>("weight_sfmu_trigger_stat_up");
  h_mu_trigger_stat_down = ctx.get_handle<float>("weight_sfmu_trigger_stat_down");
  h_mu_trigger_syst    = ctx.get_handle<float>("weight_sfmu_trigger_syst");
  h_mu_trigger_syst_up = ctx.get_handle<float>("weight_sfmu_trigger_syst_up");
  h_mu_trigger_syst_down = ctx.get_handle<float>("weight_sfmu_trigger_syst_down");
  h_mu_iso_stat        = ctx.get_handle<float>("weight_sfmu_iso_stat");
  h_mu_iso_stat_up     = ctx.get_handle<float>("weight_sfmu_iso_stat_up");
  h_mu_iso_stat_down   = ctx.get_handle<float>("weight_sfmu_iso_stat_down");
  h_mu_iso_syst        = ctx.get_handle<float>("weight_sfmu_iso_syst");
  h_mu_iso_syst_up     = ctx.get_handle<float>("weight_sfmu_iso_syst_up");
  h_mu_iso_syst_down   = ctx.get_handle<float>("weight_sfmu_iso_syst_down");
  // Pileup reweighting systematics
  h_pu                 = ctx.get_handle<float>("weight_pu");
  h_pu_up              = ctx.get_handle<float>("weight_pu_up");
  h_pu_down            = ctx.get_handle<float>("weight_pu_down");
  // Prefiring systematics
  h_prefiring          = ctx.get_handle<float>("prefiringWeight");
  h_prefiring_up       = ctx.get_handle<float>("prefiringWeightUp");
  h_prefiring_down     = ctx.get_handle<float>("prefiringWeightDown");
  // b-tagging systematics: central, cferr1, cferr2, hf, hfstats1, hfstats2, lf, lfstats1, lfstats2
  h_btag               = ctx.get_handle<float>("weight_btagdisc_central");
  h_btag_cferr1_up     = ctx.get_handle<float>("weight_btagdisc_cferr1_up");
  h_btag_cferr1_down   = ctx.get_handle<float>("weight_btagdisc_cferr1_down");
  h_btag_cferr2_up     = ctx.get_handle<float>("weight_btagdisc_cferr2_up");
  h_btag_cferr2_down   = ctx.get_handle<float>("weight_btagdisc_cferr2_down");
  h_btag_hf_up         = ctx.get_handle<float>("weight_btagdisc_hf_up");
  h_btag_hf_down       = ctx.get_handle<float>("weight_btagdisc_hf_down");
  h_btag_hfstats1_up   = ctx.get_handle<float>("weight_btagdisc_hfstats1_up");
  h_btag_hfstats1_down = ctx.get_handle<float>("weight_btagdisc_hfstats1_down");
  h_btag_hfstats2_up   = ctx.get_handle<float>("weight_btagdisc_hfstats2_up");
  h_btag_hfstats2_down = ctx.get_handle<float>("weight_btagdisc_hfstats2_down");
  h_btag_lf_up         = ctx.get_handle<float>("weight_btagdisc_lf_up");
  h_btag_lf_down       = ctx.get_handle<float>("weight_btagdisc_lf_down");
  h_btag_lfstats1_up   = ctx.get_handle<float>("weight_btagdisc_lfstats1_up");
  h_btag_lfstats1_down = ctx.get_handle<float>("weight_btagdisc_lfstats1_down");
  h_btag_lfstats2_up   = ctx.get_handle<float>("weight_btagdisc_lfstats2_up");
  h_btag_lfstats2_down = ctx.get_handle<float>("weight_btagdisc_lfstats2_down");
  // Top tagging systematics
  h_ttag               = ctx.get_handle<float>("weight_toptagsf");
  h_ttag_corr_up       = ctx.get_handle<float>("weight_toptagsf_corr_up");
  h_ttag_corr_down     = ctx.get_handle<float>("weight_toptagsf_corr_down");
  h_ttag_uncorr_up     = ctx.get_handle<float>("weight_toptagsf_uncorr_up");
  h_ttag_uncorr_down   = ctx.get_handle<float>("weight_toptagsf_uncorr_down");
  // Top mistagging systematics
  h_tmistag            = ctx.get_handle<float>("weight_topmistagsf");
  h_tmistag_up         = ctx.get_handle<float>("weight_topmistagsf_up");
  h_tmistag_down       = ctx.get_handle<float>("weight_topmistagsf_down");
  // Top pT reweighting systematics
  h_toppt_a_up         = ctx.get_handle<float>("weight_toppt_a_up");
  h_toppt_a_down       = ctx.get_handle<float>("weight_toppt_a_down");
  h_toppt_b_up         = ctx.get_handle<float>("weight_toppt_b_up");
  h_toppt_b_down       = ctx.get_handle<float>("weight_toppt_b_down");
  // muR and muF at ME level systematics
  h_murmuf_upup        = ctx.get_handle<float>("weight_murmuf_upup");
  h_murmuf_upnone      = ctx.get_handle<float>("weight_murmuf_upnone");
  h_murmuf_noneup      = ctx.get_handle<float>("weight_murmuf_noneup");
  h_murmuf_nonedown    = ctx.get_handle<float>("weight_murmuf_nonedown");
  h_murmuf_downnone    = ctx.get_handle<float>("weight_murmuf_downnone");
  h_murmuf_downdown    = ctx.get_handle<float>("weight_murmuf_downdown");
  // ISR and FSR systematics
  h_isr_up             = ctx.get_handle<float>("weight_isr_2_up");
  h_isr_down           = ctx.get_handle<float>("weight_isr_2_down");
  h_fsr_up             = ctx.get_handle<float>("weight_fsr_2_up");
  h_fsr_down           = ctx.get_handle<float>("weight_fsr_2_down");
  
  // STEP 6: Build weight maps for the NoAC weights. 
  // Creates a weight histogram for each f value, using the NoAC weights and the sumH histogram.
  // Stores them in the noac_weights_map for event by event.
  //template method start
  if(is_mc && is_tt){
    // --- NoAC setup for xi = tanh(DeltaY) (multi-f) ---
    h_xi_gen = ctx.get_handle<float>("xi_gen");
    use_noac_evtweights_ = (ctx.get("noac_apply_event_weight") == string("true"));
    noac_gen_file_ = ctx.get("noac_gen_file");
    noac_gen_hist_ = ctx.get("noac_gen_hist");

    // f values for the NoAC weights (all f values from kNoACSpecs)
    f_values = {-100.0f, -12.0f, -8.0f, -4.0f, -2.0f, -1.0f, -0.8f, -0.6f, -0.4f, -0.2f, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 2.0f, 4.0f, 8.0f, 12.0f, 100.0f};

    if(use_noac_evtweights_ && !noac_gen_file_.empty()){
      // Only initialize once - weights are shared across all instances
      if(!noac_weights_initialized){
        // reading from TTree only
        TH1D *sumH = nullptr;
        
        glob_t gl; memset(&gl, 0, sizeof(gl));
        int r = glob(noac_gen_file_.c_str(), 0, nullptr, &gl);
        int files_processed = 0;
        int files_with_tree = 0;
        int files_with_branch = 0;
        Long64_t total_events = 0;
        Long64_t total_nan_count = 0;
        Long64_t total_valid_count = 0;

        // STEP 1: Sum all of the DeltaY_xi_gen histograms from ttree which was carried from preselection
        // Reads gen-level tanh(delta|y|) histograms from all TTbar files
        if(r == 0 && gl.gl_pathc > 0){
          // Create 300-bin histogram for TTree reading
          sumH = new TH1D("Hgen_sum_from_tree", "GEN histogram from TTree xi_gen branch", 
                          300, -1.0, 1.0);
          sumH->SetDirectory(0);
          
          if (debug) cout << "INFO: Initializing NoAC weights from " << gl.gl_pathc << " preselection files..." << endl;          
          for(size_t i=0; i<gl.gl_pathc; ++i){
            const char *fp = gl.gl_pathv[i];
            std::unique_ptr<TFile> f(TFile::Open(fp));
            if(!f || f->IsZombie()) continue;
            files_processed++;
            // Try to get the TTree ("AnalysisTree" )
            TTree *tree = dynamic_cast<TTree*>(f->Get("AnalysisTree"));
            
            if(tree){
              files_with_tree++;
              // Check if xi_gen branch exists
              TBranch *br = tree->GetBranch("xi_gen");
              if(br){
                files_with_branch++;
                
                Long64_t nentries = tree->GetEntries();
                total_events += nentries;
                
                // Create a temporary histogram with unique name for this file
                TString tempH_name = TString::Format("tempH_%zu", i);
                TH1D *tempH = new TH1D(tempH_name.Data(), "temp", 300, -1.0, 1.0);
                tempH->SetDirectory(0);
                
                // Check if weight branch exists
                TBranch *wbr = tree->GetBranch("weight");
                bool has_weight = (wbr != nullptr);
                
                // Build selection string: finite xi_gen and in range [-1, 1]
                // (xi_gen == xi_gen) is false for NaN, true for finite numbers
                // (xi_gen*xi_gen < 1e10) filters out Inf values
                TString selection = "(xi_gen == xi_gen) && (xi_gen*xi_gen < 1e10) && (xi_gen >= -1.0) && (xi_gen <= 1.0)";
                TString weight_expr = has_weight ? "weight" : "1.0";
                
                // Project directly into histogram
                Long64_t nselected = tree->Project(tempH_name.Data(), "xi_gen", selection.Data(), weight_expr.Data());
                
                // Add to sum histogram
                sumH->Add(tempH);
                
                total_valid_count += nselected;
                total_nan_count += (nentries - nselected);  // Approximate
                
                delete tempH;
              } else {
                cout << "WARNING: TTree found in " << fp << " but 'xi_gen' branch not found!" << endl;
              }
            } else {
              cout << "WARNING: No TTree found in " << fp << " (tried AnalysisTree, tree, Tree)" << endl;
            }
            f->Close();
          }
          globfree(&gl);
        }  
        if(sumH){
          // Verify we have 300 bins
          if(sumH->GetNbinsX() != 300){
            // cout << "ERROR: Expected 300 bins but got " << sumH->GetNbinsX() << " bins!" << endl;
            delete sumH;
            sumH = nullptr;
          } else {
            // build the NoAC weights from the generator histogram (300 bins from TTree)
            noac_weights_map.clear();
            for(const float fv : f_values){
              try{ 
                noac_weights_map[fv] = build_noac_weights_from_gen(*sumH, fv);
              } catch(...){ 
                // Skip failed weights
                cout << "WARNING: Failed to build NoAC weights for f=" << fv << endl;
              }
            }
            noac_weights_initialized = true;  // Mark as initialized so other instances skip
            delete sumH;
          }
        }
      }  // End of initialization block
    }
  }
  //template method end


  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  init();
}

void ZprimeSemiLeptonicSystematicsHists::init(){

  //------------------------------------------------------------------- Legacy variables -------------------------------------------------------------------//
  DeltaY                    = book<TH1F>("DeltaY", "#DeltaY_{t#bar{t}} ",                                       2, -2.5, 2.5);
  DeltaY_ele_reco_up        = book<TH1F>("DeltaY_ele_reco_up",   "#DeltaY_{t#bar{t}} ele_reco_up",              2, -2.5, 2.5);// electron systematics
  DeltaY_ele_reco_down      = book<TH1F>("DeltaY_ele_reco_down", "#DeltaY_{t#bar{t}} ele_reco_down",            2, -2.5, 2.5);
  DeltaY_ele_id_up          = book<TH1F>("DeltaY_ele_id_up",   "#DeltaY_{t#bar{t}} ele_id_up",                  2, -2.5, 2.5);
  DeltaY_ele_id_down        = book<TH1F>("DeltaY_ele_id_down", "#DeltaY_{t#bar{t}} ele_id_down",                2, -2.5, 2.5);
  DeltaY_ele_trigger_up     = book<TH1F>("DeltaY_ele_trigger_up",   "#DeltaY_{t#bar{t}} ele_trigger_up",        2, -2.5, 2.5);
  DeltaY_ele_trigger_down   = book<TH1F>("DeltaY_ele_trigger_down", "#DeltaY_{t#bar{t}} ele_trigger_down",      2, -2.5, 2.5);
  DeltaY_mu_reco_up         = book<TH1F>("DeltaY_mu_reco_up",   "#DeltaY_{t#bar{t}} mu_reco_up",                2, -2.5, 2.5);// muon systematics
  DeltaY_mu_reco_down       = book<TH1F>("DeltaY_mu_reco_down", "#DeltaY_{t#bar{t}} mu_reco_down",              2, -2.5, 2.5);
  DeltaY_mu_id_stat_up      = book<TH1F>("DeltaY_mu_id_stat_up",   "#DeltaY_{t#bar{t}} mu_id_stat_up",          2, -2.5, 2.5);
  DeltaY_mu_id_stat_down    = book<TH1F>("DeltaY_mu_id_stat_down",   "#DeltaY_{t#bar{t}} mu_id_stat_down",      2, -2.5, 2.5);
  DeltaY_mu_id_syst_up      = book<TH1F>("DeltaY_mu_id_syst_up",   "#DeltaY_{t#bar{t}} mu_id_syst_up",          2, -2.5, 2.5);
  DeltaY_mu_id_syst_down    = book<TH1F>("DeltaY_mu_id_syst_down",   "#DeltaY_{t#bar{t}} mu_id_syst_down",      2, -2.5, 2.5);
  DeltaY_mu_trigger_stat_up   = book<TH1F>("DeltaY_mu_trigger_stat_up",   "#DeltaY_{t#bar{t}} mu_trigger_stat_up",     2, -2.5, 2.5);
  DeltaY_mu_trigger_stat_down = book<TH1F>("DeltaY_mu_trigger_stat_down",   "#DeltaY_{t#bar{t}} mu_trigger_stat_down", 2, -2.5, 2.5);
  DeltaY_mu_trigger_syst_up   = book<TH1F>("DeltaY_mu_trigger_syst_up",   "#DeltaY_{t#bar{t}} mu_trigger_syst_up",     2, -2.5, 2.5);
  DeltaY_mu_trigger_syst_down = book<TH1F>("DeltaY_mu_trigger_syst_down",   "#DeltaY_{t#bar{t}} mu_trigger_syst_down", 2, -2.5, 2.5);
  DeltaY_mu_iso_stat_up     = book<TH1F>("DeltaY_mu_iso_stat_up",   "#DeltaY_{t#bar{t}} mu_iso_stat_up",        2, -2.5, 2.5);
  DeltaY_mu_iso_stat_down   = book<TH1F>("DeltaY_mu_iso_stat_down",   "#DeltaY_{t#bar{t}} mu_iso_stat_down",    2, -2.5, 2.5);
  DeltaY_mu_iso_syst_up     = book<TH1F>("DeltaY_mu_iso_syst_up",   "#DeltaY_{t#bar{t}} mu_iso_syst_up",        2, -2.5, 2.5);
  DeltaY_mu_iso_syst_down   = book<TH1F>("DeltaY_mu_iso_syst_down",   "#DeltaY_{t#bar{t}} mu_iso_syst_down",    2, -2.5, 2.5);
  DeltaY_pu_up              = book<TH1F>("DeltaY_pu_up",   "#DeltaY_{t#bar{t}} pu_up",                          2, -2.5, 2.5);// Pileup reweighting systematics
  DeltaY_pu_down            = book<TH1F>("DeltaY_pu_down", "#DeltaY_{t#bar{t}} pu_down",                        2, -2.5, 2.5);
  DeltaY_prefiring_up       = book<TH1F>("DeltaY_prefiring_up",   "#DeltaY_{t#bar{t}} prefiring_up",            2, -2.5, 2.5);// Prefiring systematics
  DeltaY_prefiring_down     = book<TH1F>("DeltaY_prefiring_down", "#DeltaY_{t#bar{t}} prefiring_down",          2, -2.5, 2.5);
  DeltaY_btag_cferr1_up     = book<TH1F>("DeltaY_btag_cferr1_up", "#DeltaY_{t#bar{t}} btag_cferr1_up",          2, -2.5, 2.5);// b-tagging systematics
  DeltaY_btag_cferr1_down   = book<TH1F>("DeltaY_btag_cferr1_down", "#DeltaY_{t#bar{t}} btag_cferr1_down",      2, -2.5, 2.5);
  DeltaY_btag_cferr2_up     = book<TH1F>("DeltaY_btag_cferr2_up", "#DeltaY_{t#bar{t}} btag_cferr2_up",          2, -2.5, 2.5);
  DeltaY_btag_cferr2_down   = book<TH1F>("DeltaY_btag_cferr2_down", "#DeltaY_{t#bar{t}} btag_cferr2_down",      2, -2.5, 2.5);
  DeltaY_btag_hf_up         = book<TH1F>("DeltaY_btag_hf_up", "#DeltaY_{t#bar{t}} btag_hf_up",                  2, -2.5, 2.5);
  DeltaY_btag_hf_down       = book<TH1F>("DeltaY_btag_hf_down", "#DeltaY_{t#bar{t}} btag_hf_down",              2, -2.5, 2.5);
  DeltaY_btag_hfstats1_up   = book<TH1F>("DeltaY_btag_hfstats1_up", "#DeltaY_{t#bar{t}} btag_hfstats1_up",      2, -2.5, 2.5);
  DeltaY_btag_hfstats1_down = book<TH1F>("DeltaY_btag_hfstats1_down", "#DeltaY_{t#bar{t}} btag_hfstats1_down",  2, -2.5, 2.5);
  DeltaY_btag_hfstats2_up   = book<TH1F>("DeltaY_btag_hfstats2_up", "#DeltaY_{t#bar{t}} btag_hfstats2_up",      2, -2.5, 2.5);
  DeltaY_btag_hfstats2_down = book<TH1F>("DeltaY_btag_hfstats2_down", "#DeltaY_{t#bar{t}} btag_hfstats2_down",  2, -2.5, 2.5);
  DeltaY_btag_lf_up         = book<TH1F>("DeltaY_btag_lf_up", "#DeltaY_{t#bar{t}} btag_lf_up",                  2, -2.5, 2.5);
  DeltaY_btag_lf_down       = book<TH1F>("DeltaY_btag_lf_down", "#DeltaY_{t#bar{t}} btag_lf_down",              2, -2.5, 2.5);
  DeltaY_btag_lfstats1_up   = book<TH1F>("DeltaY_btag_lfstats1_up", "#DeltaY_{t#bar{t}} btag_lfstats1_up",      2, -2.5, 2.5);
  DeltaY_btag_lfstats1_down = book<TH1F>("DeltaY_btag_lfstats1_down", "#DeltaY_{t#bar{t}} btag_lfstats1_down",  2, -2.5, 2.5);
  DeltaY_btag_lfstats2_up   = book<TH1F>("DeltaY_btag_lfstats2_up", "#DeltaY_{t#bar{t}} btag_lfstats2_up",      2, -2.5, 2.5);
  DeltaY_btag_lfstats2_down = book<TH1F>("DeltaY_btag_lfstats2_down", "#DeltaY_{t#bar{t}} btag_lfstats2_down",  2, -2.5, 2.5);
  DeltaY_ttag_corr_up       = book<TH1F>("DeltaY_ttag_corr_up", "#DeltaY_{t#bar{t}} ttag_corr_up",              2, -2.5, 2.5);// Top tagging systematics
  DeltaY_ttag_corr_down     = book<TH1F>("DeltaY_ttag_corr_down", "#DeltaY_{t#bar{t}} ttag_corr_down",          2, -2.5, 2.5);
  DeltaY_ttag_uncorr_up     = book<TH1F>("DeltaY_ttag_uncorr_up", "#DeltaY_{t#bar{t}} ttag_uncorr_up",          2, -2.5, 2.5);
  DeltaY_ttag_uncorr_down   = book<TH1F>("DeltaY_ttag_uncorr_down", "#DeltaY_{t#bar{t}} ttag_counrr_down",      2, -2.5, 2.5);
  DeltaY_tmistag_up         = book<TH1F>("DeltaY_tmistag_up", "#DeltaY_{t#bar{t}} [GeV] tmistag_up",            2, -2.5, 2.5);// Top mistagging systematics
  DeltaY_tmistag_down       = book<TH1F>("DeltaY_tmistag_down", "#DeltaY_{t#bar{t}} [GeV] tmistag_down",        2, -2.5, 2.5);
  DeltaY_toppt_a_up         = book<TH1F>("DeltaY_toppt_a_up", "#DeltaY_{t#bar{t}} [GeV] toppt_a_up",                2, -2.5, 2.5);// Top pT reweighting systematics
  DeltaY_toppt_a_down       = book<TH1F>("DeltaY_toppt_a_down", "#DeltaY_{t#bar{t}} [GeV] toppt_a_down",            2, -2.5, 2.5);
  DeltaY_toppt_b_up         = book<TH1F>("DeltaY_toppt_b_up", "#DeltaY_{t#bar{t}} [GeV] toppt_b_up",                2, -2.5, 2.5);
  DeltaY_toppt_b_down       = book<TH1F>("DeltaY_toppt_b_down", "#DeltaY_{t#bar{t}} [GeV] toppt_b_down",            2, -2.5, 2.5);
  DeltaY_murmuf_upup        = book<TH1F>("DeltaY_murmuf_upup", "#DeltaY_{t#bar{t}} murmuf_upup",                2, -2.5, 2.5);// muR and muF at ME level systematics
  DeltaY_murmuf_upnone      = book<TH1F>("DeltaY_murmuf_upnone", "#DeltaY_{t#bar{t}} murmuf_upnone",            2, -2.5, 2.5);
  DeltaY_murmuf_noneup      = book<TH1F>("DeltaY_murmuf_noneup", "#DeltaY_{t#bar{t}} murmuf_noneup",            2, -2.5, 2.5);
  DeltaY_murmuf_nonedown    = book<TH1F>("DeltaY_murmuf_nonedown", "#DeltaY_{t#bar{t}} murmuf_nonedown",        2, -2.5, 2.5);
  DeltaY_murmuf_downnone    = book<TH1F>("DeltaY_murmuf_downnone", "#DeltaY_{t#bar{t}} murmuf_downnone",        2, -2.5, 2.5);
  DeltaY_murmuf_downdown    = book<TH1F>("DeltaY_murmuf_downdown", "#DeltaY_{t#bar{t}} murmuf_downdown",        2, -2.5, 2.5);
  DeltaY_isr_up             = book<TH1F>("DeltaY_isr_up", "#DeltaY_{t#bar{t}} isr_up",                          2, -2.5, 2.5);// ISR and FSR systematics
  DeltaY_isr_down           = book<TH1F>("DeltaY_isr_down", "#DeltaY_{t#bar{t}} isr_down",                      2, -2.5, 2.5);
  DeltaY_fsr_up             = book<TH1F>("DeltaY_fsr_up", "#DeltaY_{t#bar{t}} fsr_up",                          2, -2.5, 2.5);
  DeltaY_fsr_down           = book<TH1F>("DeltaY_fsr_down", "#DeltaY_{t#bar{t}} fsr_down",                      2, -2.5, 2.5);

  // Response matrixs for unfolding
  DeltaY_tt                    = book<TH2F>("DeltaY_tt", "#DeltaY_{t#bar{t}} ",                                       2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_reco_up_tt         = book<TH2F>("DeltaY_mu_reco_up_tt",   "#DeltaY_{t#bar{t}} mu_reco_up",                2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_reco_down_tt       = book<TH2F>("DeltaY_mu_reco_down_tt", "#DeltaY_{t#bar{t}} mu_reco_down",              2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_pu_up_tt              = book<TH2F>("DeltaY_pu_up_tt",   "#DeltaY_{t#bar{t}} pu_up",                          2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_pu_down_tt            = book<TH2F>("DeltaY_pu_down_tt", "#DeltaY_{t#bar{t}} pu_down",                        2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_prefiring_up_tt       = book<TH2F>("DeltaY_prefiring_up_tt",   "#DeltaY_{t#bar{t}} prefiring_up",            2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_prefiring_down_tt     = book<TH2F>("DeltaY_prefiring_down_tt", "#DeltaY_{t#bar{t}} prefiring_down",          2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_id_stat_up_tt      = book<TH2F>("DeltaY_mu_id_stat_up_tt",   "#DeltaY_{t#bar{t}} mu_id_stat_up",         2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_id_stat_down_tt    = book<TH2F>("DeltaY_mu_id_stat_down_tt",   "#DeltaY_{t#bar{t}} mu_id_stat_down",     2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_id_syst_up_tt      = book<TH2F>("DeltaY_mu_id_syst_up_tt",   "#DeltaY_{t#bar{t}} mu_id_syst_up",         2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_id_syst_down_tt    = book<TH2F>("DeltaY_mu_id_syst_down_tt",   "#DeltaY_{t#bar{t}} mu_id_syst_down",     2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_iso_stat_up_tt     = book<TH2F>("DeltaY_mu_iso_stat_up_tt",   "#DeltaY_{t#bar{t}} mu_iso_stat_up",        2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_iso_stat_down_tt   = book<TH2F>("DeltaY_mu_iso_stat_down_tt",   "#DeltaY_{t#bar{t}} mu_iso_stat_down",    2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_iso_syst_up_tt     = book<TH2F>("DeltaY_mu_iso_syst_up_tt",   "#DeltaY_{t#bar{t}} mu_iso_syst_up",       2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_iso_syst_down_tt   = book<TH2F>("DeltaY_mu_iso_syst_down_tt",   "#DeltaY_{t#bar{t}} mu_iso_syst_down",    2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_trigger_stat_up_tt     = book<TH2F>("DeltaY_mu_trigger_stat_up_tt",   "#DeltaY_{t#bar{t}} mu_trigger_stat_up",        2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_trigger_stat_down_tt   = book<TH2F>("DeltaY_mu_trigger_stat_down_tt",   "#DeltaY_{t#bar{t}} mu_trigger_stat_down",    2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_trigger_syst_up_tt     = book<TH2F>("DeltaY_mu_trigger_syst_up_tt",   "#DeltaY_{t#bar{t}} mu_trigger_syst_up",        2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_mu_trigger_syst_down_tt   = book<TH2F>("DeltaY_mu_trigger_syst_down_tt",   "#DeltaY_{t#bar{t}} mu_trigger_syst_down",    2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ele_id_up_tt          = book<TH2F>("DeltaY_ele_id_up_tt",   "#DeltaY_{t#bar{t}} ele_id_up",                  2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ele_id_down_tt        = book<TH2F>("DeltaY_ele_id_down_tt", "#DeltaY_{t#bar{t}} ele_id_down",                2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ele_trigger_up_tt     = book<TH2F>("DeltaY_ele_trigger_up_tt",   "#DeltaY_{t#bar{t}} ele_trigger_up",        2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ele_trigger_down_tt   = book<TH2F>("DeltaY_ele_trigger_down_tt", "#DeltaY_{t#bar{t}} ele_trigger_down",      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ele_reco_up_tt        = book<TH2F>("DeltaY_ele_reco_up_tt",   "#DeltaY_{t#bar{t}} ele_reco_up",              2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ele_reco_down_tt      = book<TH2F>("DeltaY_ele_reco_down_tt", "#DeltaY_{t#bar{t}} ele_reco_down",            2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_murmuf_upup_tt        = book<TH2F>("DeltaY_murmuf_upup_tt", "#DeltaY_{t#bar{t}} murmuf_upup",                2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_murmuf_upnone_tt      = book<TH2F>("DeltaY_murmuf_upnone_tt", "#DeltaY_{t#bar{t}} murmuf_upnone",            2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_murmuf_noneup_tt      = book<TH2F>("DeltaY_murmuf_noneup_tt", "#DeltaY_{t#bar{t}} murmuf_noneup",            2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_murmuf_nonedown_tt    = book<TH2F>("DeltaY_murmuf_nonedown_tt", "#DeltaY_{t#bar{t}} murmuf_nonedown",        2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_murmuf_downnone_tt    = book<TH2F>("DeltaY_murmuf_downnone_tt", "#DeltaY_{t#bar{t}} murmuf_downnone",        2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_murmuf_downdown_tt    = book<TH2F>("DeltaY_murmuf_downdown_tt", "#DeltaY_{t#bar{t}} murmuf_downdown",        2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_isr_up_tt             = book<TH2F>("DeltaY_isr_up_tt", "#DeltaY_{t#bar{t}} isr_up",                          2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_isr_down_tt           = book<TH2F>("DeltaY_isr_down_tt", "#DeltaY_{t#bar{t}} isr_down",                      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_fsr_up_tt             = book<TH2F>("DeltaY_fsr_up_tt", "#DeltaY_{t#bar{t}} fsr_up",                          2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_fsr_down_tt           = book<TH2F>("DeltaY_fsr_down_tt", "#DeltaY_{t#bar{t}} fsr_down",                      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_cferr1_up_tt     = book<TH2F>("DeltaY_btag_cferr1_up_tt", "#DeltaY_{t#bar{t}} btag_cferr1_up",          2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_cferr1_down_tt   = book<TH2F>("DeltaY_btag_cferr1_down_tt", "#DeltaY_{t#bar{t}} btag_cferr1_down",      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_cferr2_up_tt     = book<TH2F>("DeltaY_btag_cferr2_up_tt", "#DeltaY_{t#bar{t}} btag_cferr2_up",          2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_cferr2_down_tt   = book<TH2F>("DeltaY_btag_cferr2_down_tt", "#DeltaY_{t#bar{t}} btag_cferr2_down",      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_hf_up_tt         = book<TH2F>("DeltaY_btag_hf_up_tt", "#DeltaY_{t#bar{t}} btag_hf_up",                  2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_hf_down_tt       = book<TH2F>("DeltaY_btag_hf_down_tt", "#DeltaY_{t#bar{t}} btag_hf_down",              2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_hfstats1_up_tt   = book<TH2F>("DeltaY_btag_hfstats1_up_tt", "#DeltaY_{t#bar{t}} btag_hfstats1_up",      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_hfstats1_down_tt = book<TH2F>("DeltaY_btag_hfstats1_down_tt", "#DeltaY_{t#bar{t}} btag_hfstats1_down",  2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_hfstats2_up_tt   = book<TH2F>("DeltaY_btag_hfstats2_up_tt", "#DeltaY_{t#bar{t}} btag_hfstats2_up",      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_hfstats2_down_tt = book<TH2F>("DeltaY_btag_hfstats2_down_tt", "#DeltaY_{t#bar{t}} btag_hfstats2_down",  2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_lf_up_tt         = book<TH2F>("DeltaY_btag_lf_up_tt", "#DeltaY_{t#bar{t}} btag_lf_up",                  2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_lf_down_tt       = book<TH2F>("DeltaY_btag_lf_down_tt", "#DeltaY_{t#bar{t}} btag_lf_down",              2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_lfstats1_up_tt   = book<TH2F>("DeltaY_btag_lfstats1_up_tt", "#DeltaY_{t#bar{t}} btag_lfstats1_up",      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_lfstats1_down_tt = book<TH2F>("DeltaY_btag_lfstats1_down_tt", "#DeltaY_{t#bar{t}} btag_lfstats1_down",  2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_lfstats2_up_tt   = book<TH2F>("DeltaY_btag_lfstats2_up_tt", "#DeltaY_{t#bar{t}} btag_lfstats2_up",      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_btag_lfstats2_down_tt = book<TH2F>("DeltaY_btag_lfstats2_down_tt", "#DeltaY_{t#bar{t}} btag_lfstats2_down",  2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ttag_corr_up_tt       = book<TH2F>("DeltaY_ttag_corr_up_tt", "#DeltaY_{t#bar{t}} ttag_corr_up",              2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ttag_corr_down_tt     = book<TH2F>("DeltaY_ttag_corr_down_tt", "#DeltaY_{t#bar{t}} ttag_corr_down",          2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ttag_uncorr_up_tt     = book<TH2F>("DeltaY_ttag_uncorr_up_tt", "#DeltaY_{t#bar{t}} ttag_uncorr_up",          2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_ttag_uncorr_down_tt   = book<TH2F>("DeltaY_ttag_uncorr_down_tt", "#DeltaY_{t#bar{t}} ttag_counrr_down",      2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_tmistag_up_tt         = book<TH2F>("DeltaY_tmistag_up_tt", "#DeltaY_{t#bar{t}} [GeV] tmistag_up",            2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_tmistag_down_tt       = book<TH2F>("DeltaY_tmistag_down_tt", "#DeltaY_{t#bar{t}} [GeV] tmistag_down",        2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_toppt_a_up_tt         = book<TH2F>("DeltaY_toppt_a_up_tt", "#DeltaY_{t#bar{t}} [GeV] toppt_a_up",                2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_toppt_a_down_tt       = book<TH2F>("DeltaY_toppt_a_down_tt", "#DeltaY_{t#bar{t}} [GeV] toppt_a_down",            2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_toppt_b_up_tt         = book<TH2F>("DeltaY_toppt_b_up_tt", "#DeltaY_{t#bar{t}} [GeV] toppt_b_up",                2, -2.5, 2.5, 2, -2.5, 2.5);
  DeltaY_toppt_b_down_tt       = book<TH2F>("DeltaY_toppt_b_down_tt", "#DeltaY_{t#bar{t}} [GeV] toppt_b_down",            2, -2.5, 2.5, 2, -2.5, 2.5);
  //------------------------------------------------------------------- END legacy variables ---------------------------------------------------------------//



  //---------------------------------------------------------------------- EFT variables ----------------------------------------------------------------------//
  // DeltaY_reco_d1                    = book<TH1F>("DeltaY_reco_d1",   "#DeltaY_reco_d1_{t#bar{t}} ",                                     2, -2.5, 2.5);
  DeltaY_reco_d1_mu_reco_up         = book<TH1F>("DeltaY_reco_d1_mu_reco_up",   "#DeltaY_reco_d1_{t#bar{t}} mu_reco_up",                2, -2.5, 2.5);
  DeltaY_reco_d1_mu_reco_down       = book<TH1F>("DeltaY_reco_d1_mu_reco_down", "#DeltaY_reco_d1_{t#bar{t}} mu_reco_down",              2, -2.5, 2.5);
  DeltaY_reco_d1_pu_up              = book<TH1F>("DeltaY_reco_d1_pu_up",   "#DeltaY_reco_d1_{t#bar{t}} pu_up",                          2, -2.5, 2.5);
  DeltaY_reco_d1_pu_down            = book<TH1F>("DeltaY_reco_d1_pu_down", "#DeltaY_reco_d1_{t#bar{t}} pu_down",                        2, -2.5, 2.5);
  DeltaY_reco_d1_prefiring_up       = book<TH1F>("DeltaY_reco_d1_prefiring_up",   "#DeltaY_reco_d1_{t#bar{t}} prefiring_up",            2, -2.5, 2.5);
  DeltaY_reco_d1_prefiring_down     = book<TH1F>("DeltaY_reco_d1_prefiring_down", "#DeltaY_reco_d1_{t#bar{t}} prefiring_down",          2, -2.5, 2.5);
  DeltaY_reco_d1_mu_id_stat_up      = book<TH1F>("DeltaY_reco_d1_mu_id_stat_up",   "#DeltaY_reco_d1_{t#bar{t}} mu_id_stat_up",          2, -2.5, 2.5);
  DeltaY_reco_d1_mu_id_stat_down    = book<TH1F>("DeltaY_reco_d1_mu_id_stat_down",   "#DeltaY_reco_d1_{t#bar{t}} mu_id_stat_down",      2, -2.5, 2.5);
  DeltaY_reco_d1_mu_id_syst_up      = book<TH1F>("DeltaY_reco_d1_mu_id_syst_up",   "#DeltaY_reco_d1_{t#bar{t}} mu_id_syst_up",          2, -2.5, 2.5);
  DeltaY_reco_d1_mu_id_syst_down    = book<TH1F>("DeltaY_reco_d1_mu_id_syst_down",   "#DeltaY_reco_d1_{t#bar{t}} mu_id_syst_down",      2, -2.5, 2.5);
  DeltaY_reco_d1_mu_iso_stat_up     = book<TH1F>("DeltaY_reco_d1_mu_iso_stat_up",   "#DeltaY_reco_d1_{t#bar{t}} mu_iso_stat_up",        2, -2.5, 2.5);
  DeltaY_reco_d1_mu_iso_stat_down   = book<TH1F>("DeltaY_reco_d1_mu_iso_stat_down",   "#DeltaY_reco_d1_{t#bar{t}} mu_iso_stat_down",    2, -2.5, 2.5);
  DeltaY_reco_d1_mu_iso_syst_up     = book<TH1F>("DeltaY_reco_d1_mu_iso_syst_up",   "#DeltaY_reco_d1_{t#bar{t}} mu_iso_syst_up",        2, -2.5, 2.5);
  DeltaY_reco_d1_mu_iso_syst_down   = book<TH1F>("DeltaY_reco_d1_mu_iso_syst_down",   "#DeltaY_reco_d1_{t#bar{t}} mu_iso_syst_down",    2, -2.5, 2.5);
  DeltaY_reco_d1_mu_trigger_stat_up     = book<TH1F>("DeltaY_reco_d1_mu_trigger_stat_up",   "#DeltaY_reco_d1_{t#bar{t}} mu_trigger_stat_up",        2, -2.5, 2.5);
  DeltaY_reco_d1_mu_trigger_stat_down   = book<TH1F>("DeltaY_reco_d1_mu_trigger_stat_down",   "#DeltaY_reco_d1_{t#bar{t}} mu_trigger_stat_down",    2, -2.5, 2.5);
  DeltaY_reco_d1_mu_trigger_syst_up     = book<TH1F>("DeltaY_reco_d1_mu_trigger_syst_up",   "#DeltaY_reco_d1_{t#bar{t}} mu_trigger_syst_up",        2, -2.5, 2.5);
  DeltaY_reco_d1_mu_trigger_syst_down   = book<TH1F>("DeltaY_reco_d1_mu_trigger_syst_down",   "#DeltaY_reco_d1_{t#bar{t}} mu_trigger_syst_down",    2, -2.5, 2.5);
  DeltaY_reco_d1_ele_id_up          = book<TH1F>("DeltaY_reco_d1_ele_id_up",   "#DeltaY_reco_d1_{t#bar{t}} ele_id_up",                  2, -2.5, 2.5);
  DeltaY_reco_d1_ele_id_down        = book<TH1F>("DeltaY_reco_d1_ele_id_down", "#DeltaY_reco_d1_{t#bar{t}} ele_id_down",                2, -2.5, 2.5);
  DeltaY_reco_d1_ele_trigger_up     = book<TH1F>("DeltaY_reco_d1_ele_trigger_up",   "#DeltaY_reco_d1_{t#bar{t}} ele_trigger_up",        2, -2.5, 2.5);
  DeltaY_reco_d1_ele_trigger_down   = book<TH1F>("DeltaY_reco_d1_ele_trigger_down", "#DeltaY_reco_d1_{t#bar{t}} ele_trigger_down",      2, -2.5, 2.5);
  DeltaY_reco_d1_ele_reco_up        = book<TH1F>("DeltaY_reco_d1_ele_reco_up",   "#DeltaY_reco_d1_{t#bar{t}} ele_reco_up",              2, -2.5, 2.5);
  DeltaY_reco_d1_ele_reco_down      = book<TH1F>("DeltaY_reco_d1_ele_reco_down", "#DeltaY_reco_d1_{t#bar{t}} ele_reco_down",            2, -2.5, 2.5);
  DeltaY_reco_d1_murmuf_upup        = book<TH1F>("DeltaY_reco_d1_murmuf_upup", "#DeltaY_reco_d1_{t#bar{t}} murmuf_upup",                2, -2.5, 2.5);
  DeltaY_reco_d1_murmuf_upnone      = book<TH1F>("DeltaY_reco_d1_murmuf_upnone", "#DeltaY_reco_d1_{t#bar{t}} murmuf_upnone",            2, -2.5, 2.5);
  DeltaY_reco_d1_murmuf_noneup      = book<TH1F>("DeltaY_reco_d1_murmuf_noneup", "#DeltaY_reco_d1_{t#bar{t}} murmuf_noneup",            2, -2.5, 2.5);
  DeltaY_reco_d1_murmuf_nonedown    = book<TH1F>("DeltaY_reco_d1_murmuf_nonedown", "#DeltaY_reco_d1_{t#bar{t}} murmuf_nonedown",        2, -2.5, 2.5);
  DeltaY_reco_d1_murmuf_downnone    = book<TH1F>("DeltaY_reco_d1_murmuf_downnone", "#DeltaY_reco_d1_{t#bar{t}} murmuf_downnone",        2, -2.5, 2.5);
  DeltaY_reco_d1_murmuf_downdown    = book<TH1F>("DeltaY_reco_d1_murmuf_downdown", "#DeltaY_reco_d1_{t#bar{t}} murmuf_downdown",        2, -2.5, 2.5);
  DeltaY_reco_d1_isr_up             = book<TH1F>("DeltaY_reco_d1_isr_up", "#DeltaY_reco_d1_{t#bar{t}} isr_up",                          2, -2.5, 2.5);
  DeltaY_reco_d1_isr_down           = book<TH1F>("DeltaY_reco_d1_isr_down", "#DeltaY_reco_d1_{t#bar{t}} isr_down",                      2, -2.5, 2.5);
  DeltaY_reco_d1_fsr_up             = book<TH1F>("DeltaY_reco_d1_fsr_up", "#DeltaY_reco_d1_{t#bar{t}} fsr_up",                          2, -2.5, 2.5);
  DeltaY_reco_d1_fsr_down           = book<TH1F>("DeltaY_reco_d1_fsr_down", "#DeltaY_reco_d1_{t#bar{t}} fsr_down",                      2, -2.5, 2.5);
  DeltaY_reco_d1_btag_cferr1_up     = book<TH1F>("DeltaY_reco_d1_btag_cferr1_up", "#DeltaY_reco_d1_{t#bar{t}} btag_cferr1_up",          2, -2.5, 2.5);
  DeltaY_reco_d1_btag_cferr1_down   = book<TH1F>("DeltaY_reco_d1_btag_cferr1_down", "#DeltaY_reco_d1_{t#bar{t}} btag_cferr1_down",      2, -2.5, 2.5);
  DeltaY_reco_d1_btag_cferr2_up     = book<TH1F>("DeltaY_reco_d1_btag_cferr2_up", "#DeltaY_reco_d1_{t#bar{t}} btag_cferr2_up",          2, -2.5, 2.5);
  DeltaY_reco_d1_btag_cferr2_down   = book<TH1F>("DeltaY_reco_d1_btag_cferr2_down", "#DeltaY_reco_d1_{t#bar{t}} btag_cferr2_down",      2, -2.5, 2.5);
  DeltaY_reco_d1_btag_hf_up         = book<TH1F>("DeltaY_reco_d1_btag_hf_up", "#DeltaY_reco_d1_{t#bar{t}} btag_hf_up",                  2, -2.5, 2.5);
  DeltaY_reco_d1_btag_hf_down       = book<TH1F>("DeltaY_reco_d1_btag_hf_down", "#DeltaY_reco_d1_{t#bar{t}} btag_hf_down",              2, -2.5, 2.5);
  DeltaY_reco_d1_btag_hfstats1_up   = book<TH1F>("DeltaY_reco_d1_btag_hfstats1_up", "#DeltaY_reco_d1_{t#bar{t}} btag_hfstats1_up",      2, -2.5, 2.5);
  DeltaY_reco_d1_btag_hfstats1_down = book<TH1F>("DeltaY_reco_d1_btag_hfstats1_down", "#DeltaY_reco_d1_{t#bar{t}} btag_hfstats1_down",  2, -2.5, 2.5);
  DeltaY_reco_d1_btag_hfstats2_up   = book<TH1F>("DeltaY_reco_d1_btag_hfstats2_up", "#DeltaY_reco_d1_{t#bar{t}} btag_hfstats2_up",      2, -2.5, 2.5);
  DeltaY_reco_d1_btag_hfstats2_down = book<TH1F>("DeltaY_reco_d1_btag_hfstats2_down", "#DeltaY_reco_d1_{t#bar{t}} btag_hfstats2_down",  2, -2.5, 2.5);
  DeltaY_reco_d1_btag_lf_up         = book<TH1F>("DeltaY_reco_d1_btag_lf_up", "#DeltaY_reco_d1_{t#bar{t}} btag_lf_up",                  2, -2.5, 2.5);
  DeltaY_reco_d1_btag_lf_down       = book<TH1F>("DeltaY_reco_d1_btag_lf_down", "#DeltaY_reco_d1_{t#bar{t}} btag_lf_down",              2, -2.5, 2.5);
  DeltaY_reco_d1_btag_lfstats1_up   = book<TH1F>("DeltaY_reco_d1_btag_lfstats1_up", "#DeltaY_reco_d1_{t#bar{t}} btag_lfstats1_up",      2, -2.5, 2.5);
  DeltaY_reco_d1_btag_lfstats1_down = book<TH1F>("DeltaY_reco_d1_btag_lfstats1_down", "#DeltaY_reco_d1_{t#bar{t}} btag_lfstats1_down",  2, -2.5, 2.5);
  DeltaY_reco_d1_btag_lfstats2_up   = book<TH1F>("DeltaY_reco_d1_btag_lfstats2_up", "#DeltaY_reco_d1_{t#bar{t}} btag_lfstats2_up",      2, -2.5, 2.5);
  DeltaY_reco_d1_btag_lfstats2_down = book<TH1F>("DeltaY_reco_d1_btag_lfstats2_down", "#DeltaY_reco_d1_{t#bar{t}} btag_lfstats2_down",  2, -2.5, 2.5);
  DeltaY_reco_d1_ttag_corr_up       = book<TH1F>("DeltaY_reco_d1_ttag_corr_up", "#DeltaY_reco_d1_{t#bar{t}} ttag_corr_up",              2, -2.5, 2.5);
  DeltaY_reco_d1_ttag_corr_down     = book<TH1F>("DeltaY_reco_d1_ttag_corr_down", "#DeltaY_reco_d1_{t#bar{t}} ttag_corr_down",          2, -2.5, 2.5);
  DeltaY_reco_d1_ttag_uncorr_up     = book<TH1F>("DeltaY_reco_d1_ttag_uncorr_up", "#DeltaY_reco_d1_{t#bar{t}} ttag_uncorr_up",          2, -2.5, 2.5);
  DeltaY_reco_d1_ttag_uncorr_down   = book<TH1F>("DeltaY_reco_d1_ttag_uncorr_down", "#DeltaY_reco_d1_{t#bar{t}} ttag_counrr_down",      2, -2.5, 2.5);
  DeltaY_reco_d1_tmistag_up         = book<TH1F>("DeltaY_reco_d1_tmistag_up", "#DeltaY_reco_d1_{t#bar{t}} [GeV] tmistag_up",            2, -2.5, 2.5);
  DeltaY_reco_d1_tmistag_down       = book<TH1F>("DeltaY_reco_d1_tmistag_down", "#DeltaY_reco_d1_{t#bar{t}} [GeV] tmistag_down",        2, -2.5, 2.5);
  DeltaY_reco_d1_toppt_a_up         = book<TH1F>("DeltaY_reco_d1_toppt_a_up", "#DeltaY_reco_d1_{t#bar{t}} [GeV] toppt_a_up",                2, -2.5, 2.5);
  DeltaY_reco_d1_toppt_a_down       = book<TH1F>("DeltaY_reco_d1_toppt_a_down", "#DeltaY_reco_d1_{t#bar{t}} [GeV] toppt_a_down",            2, -2.5, 2.5);
  DeltaY_reco_d1_toppt_b_up         = book<TH1F>("DeltaY_reco_d1_toppt_b_up", "#DeltaY_reco_d1_{t#bar{t}} [GeV] toppt_b_up",                2, -2.5, 2.5);
  DeltaY_reco_d1_toppt_b_down       = book<TH1F>("DeltaY_reco_d1_toppt_b_down", "#DeltaY_reco_d1_{t#bar{t}} [GeV] toppt_b_down",            2, -2.5, 2.5);

  // DeltaY_reco_d2                    = book<TH1F>("DeltaY_reco_d2",   "#DeltaY_reco_d2_{t#bar{t}} ",                                     2, -2.5, 2.5);
  DeltaY_reco_d2_mu_reco_up         = book<TH1F>("DeltaY_reco_d2_mu_reco_up",   "#DeltaY_reco_d2_{t#bar{t}} mu_reco_up",                2, -2.5, 2.5);
  DeltaY_reco_d2_mu_reco_down       = book<TH1F>("DeltaY_reco_d2_mu_reco_down", "#DeltaY_reco_d2_{t#bar{t}} mu_reco_down",              2, -2.5, 2.5);
  DeltaY_reco_d2_pu_up              = book<TH1F>("DeltaY_reco_d2_pu_up",   "#DeltaY_reco_d2_{t#bar{t}} pu_up",                          2, -2.5, 2.5);
  DeltaY_reco_d2_pu_down            = book<TH1F>("DeltaY_reco_d2_pu_down", "#DeltaY_reco_d2_{t#bar{t}} pu_down",                        2, -2.5, 2.5);
  DeltaY_reco_d2_prefiring_up       = book<TH1F>("DeltaY_reco_d2_prefiring_up",   "#DeltaY_reco_d2_{t#bar{t}} prefiring_up",            2, -2.5, 2.5);
  DeltaY_reco_d2_prefiring_down     = book<TH1F>("DeltaY_reco_d2_prefiring_down", "#DeltaY_reco_d2_{t#bar{t}} prefiring_down",          2, -2.5, 2.5);
  DeltaY_reco_d2_mu_id_stat_up      = book<TH1F>("DeltaY_reco_d2_mu_id_stat_up",   "#DeltaY_reco_d2_{t#bar{t}} mu_id_stat_up",          2, -2.5, 2.5);
  DeltaY_reco_d2_mu_id_stat_down    = book<TH1F>("DeltaY_reco_d2_mu_id_stat_down",   "#DeltaY_reco_d2_{t#bar{t}} mu_id_stat_down",      2, -2.5, 2.5);
  DeltaY_reco_d2_mu_id_syst_up      = book<TH1F>("DeltaY_reco_d2_mu_id_syst_up",   "#DeltaY_reco_d2_{t#bar{t}} mu_id_syst_up",          2, -2.5, 2.5);
  DeltaY_reco_d2_mu_id_syst_down    = book<TH1F>("DeltaY_reco_d2_mu_id_syst_down",   "#DeltaY_reco_d2_{t#bar{t}} mu_id_syst_down",      2, -2.5, 2.5);
  DeltaY_reco_d2_mu_iso_stat_up     = book<TH1F>("DeltaY_reco_d2_mu_iso_stat_up",   "#DeltaY_reco_d2_{t#bar{t}} mu_iso_stat_up",        2, -2.5, 2.5);
  DeltaY_reco_d2_mu_iso_stat_down   = book<TH1F>("DeltaY_reco_d2_mu_iso_stat_down",   "#DeltaY_reco_d2_{t#bar{t}} mu_iso_stat_down",    2, -2.5, 2.5);
  DeltaY_reco_d2_mu_iso_syst_up     = book<TH1F>("DeltaY_reco_d2_mu_iso_syst_up",   "#DeltaY_reco_d2_{t#bar{t}} mu_iso_syst_up",        2, -2.5, 2.5);
  DeltaY_reco_d2_mu_iso_syst_down   = book<TH1F>("DeltaY_reco_d2_mu_iso_syst_down",   "#DeltaY_reco_d2_{t#bar{t}} mu_iso_syst_down",    2, -2.5, 2.5);
  DeltaY_reco_d2_mu_trigger_stat_up     = book<TH1F>("DeltaY_reco_d2_mu_trigger_stat_up",   "#DeltaY_reco_d2_{t#bar{t}} mu_trigger_stat_up",        2, -2.5, 2.5);
  DeltaY_reco_d2_mu_trigger_stat_down   = book<TH1F>("DeltaY_reco_d2_mu_trigger_stat_down",   "#DeltaY_reco_d2_{t#bar{t}} mu_trigger_stat_down",    2, -2.5, 2.5);
  DeltaY_reco_d2_mu_trigger_syst_up     = book<TH1F>("DeltaY_reco_d2_mu_trigger_syst_up",   "#DeltaY_reco_d2_{t#bar{t}} mu_trigger_syst_up",        2, -2.5, 2.5);
  DeltaY_reco_d2_mu_trigger_syst_down   = book<TH1F>("DeltaY_reco_d2_mu_trigger_syst_down",   "#DeltaY_reco_d2_{t#bar{t}} mu_trigger_syst_down",    2, -2.5, 2.5);
  DeltaY_reco_d2_ele_id_up          = book<TH1F>("DeltaY_reco_d2_ele_id_up",   "#DeltaY_reco_d2_{t#bar{t}} ele_id_up",                  2, -2.5, 2.5);
  DeltaY_reco_d2_ele_id_down        = book<TH1F>("DeltaY_reco_d2_ele_id_down", "#DeltaY_reco_d2_{t#bar{t}} ele_id_down",                2, -2.5, 2.5);
  DeltaY_reco_d2_ele_trigger_up     = book<TH1F>("DeltaY_reco_d2_ele_trigger_up",   "#DeltaY_reco_d2_{t#bar{t}} ele_trigger_up",        2, -2.5, 2.5);
  DeltaY_reco_d2_ele_trigger_down   = book<TH1F>("DeltaY_reco_d2_ele_trigger_down", "#DeltaY_reco_d2_{t#bar{t}} ele_trigger_down",      2, -2.5, 2.5);
  DeltaY_reco_d2_ele_reco_up        = book<TH1F>("DeltaY_reco_d2_ele_reco_up",   "#DeltaY_reco_d2_{t#bar{t}} ele_reco_up",              2, -2.5, 2.5);
  DeltaY_reco_d2_ele_reco_down      = book<TH1F>("DeltaY_reco_d2_ele_reco_down", "#DeltaY_reco_d2_{t#bar{t}} ele_reco_down",            2, -2.5, 2.5);
  DeltaY_reco_d2_murmuf_upup        = book<TH1F>("DeltaY_reco_d2_murmuf_upup", "#DeltaY_reco_d2_{t#bar{t}} murmuf_upup",                2, -2.5, 2.5);
  DeltaY_reco_d2_murmuf_upnone      = book<TH1F>("DeltaY_reco_d2_murmuf_upnone", "#DeltaY_reco_d2_{t#bar{t}} murmuf_upnone",            2, -2.5, 2.5);
  DeltaY_reco_d2_murmuf_noneup      = book<TH1F>("DeltaY_reco_d2_murmuf_noneup", "#DeltaY_reco_d2_{t#bar{t}} murmuf_noneup",            2, -2.5, 2.5);
  DeltaY_reco_d2_murmuf_nonedown    = book<TH1F>("DeltaY_reco_d2_murmuf_nonedown", "#DeltaY_reco_d2_{t#bar{t}} murmuf_nonedown",        2, -2.5, 2.5);
  DeltaY_reco_d2_murmuf_downnone    = book<TH1F>("DeltaY_reco_d2_murmuf_downnone", "#DeltaY_reco_d2_{t#bar{t}} murmuf_downnone",        2, -2.5, 2.5);
  DeltaY_reco_d2_murmuf_downdown    = book<TH1F>("DeltaY_reco_d2_murmuf_downdown", "#DeltaY_reco_d2_{t#bar{t}} murmuf_downdown",        2, -2.5, 2.5);
  DeltaY_reco_d2_isr_up             = book<TH1F>("DeltaY_reco_d2_isr_up", "#DeltaY_reco_d2_{t#bar{t}} isr_up",                          2, -2.5, 2.5);
  DeltaY_reco_d2_isr_down           = book<TH1F>("DeltaY_reco_d2_isr_down", "#DeltaY_reco_d2_{t#bar{t}} isr_down",                      2, -2.5, 2.5);
  DeltaY_reco_d2_fsr_up             = book<TH1F>("DeltaY_reco_d2_fsr_up", "#DeltaY_reco_d2_{t#bar{t}} fsr_up",                          2, -2.5, 2.5);
  DeltaY_reco_d2_fsr_down           = book<TH1F>("DeltaY_reco_d2_fsr_down", "#DeltaY_reco_d2_{t#bar{t}} fsr_down",                      2, -2.5, 2.5);
  DeltaY_reco_d2_btag_cferr1_up     = book<TH1F>("DeltaY_reco_d2_btag_cferr1_up", "#DeltaY_reco_d2_{t#bar{t}} btag_cferr1_up",          2, -2.5, 2.5);
  DeltaY_reco_d2_btag_cferr1_down   = book<TH1F>("DeltaY_reco_d2_btag_cferr1_down", "#DeltaY_reco_d2_{t#bar{t}} btag_cferr1_down",      2, -2.5, 2.5);
  DeltaY_reco_d2_btag_cferr2_up     = book<TH1F>("DeltaY_reco_d2_btag_cferr2_up", "#DeltaY_reco_d2_{t#bar{t}} btag_cferr2_up",          2, -2.5, 2.5);
  DeltaY_reco_d2_btag_cferr2_down   = book<TH1F>("DeltaY_reco_d2_btag_cferr2_down", "#DeltaY_reco_d2_{t#bar{t}} btag_cferr2_down",      2, -2.5, 2.5);
  DeltaY_reco_d2_btag_hf_up         = book<TH1F>("DeltaY_reco_d2_btag_hf_up", "#DeltaY_reco_d2_{t#bar{t}} btag_hf_up",                  2, -2.5, 2.5);
  DeltaY_reco_d2_btag_hf_down       = book<TH1F>("DeltaY_reco_d2_btag_hf_down", "#DeltaY_reco_d2_{t#bar{t}} btag_hf_down",              2, -2.5, 2.5);
  DeltaY_reco_d2_btag_hfstats1_up   = book<TH1F>("DeltaY_reco_d2_btag_hfstats1_up", "#DeltaY_reco_d2_{t#bar{t}} btag_hfstats1_up",      2, -2.5, 2.5);
  DeltaY_reco_d2_btag_hfstats1_down = book<TH1F>("DeltaY_reco_d2_btag_hfstats1_down", "#DeltaY_reco_d2_{t#bar{t}} btag_hfstats1_down",  2, -2.5, 2.5);
  DeltaY_reco_d2_btag_hfstats2_up   = book<TH1F>("DeltaY_reco_d2_btag_hfstats2_up", "#DeltaY_reco_d2_{t#bar{t}} btag_hfstats2_up",      2, -2.5, 2.5);
  DeltaY_reco_d2_btag_hfstats2_down = book<TH1F>("DeltaY_reco_d2_btag_hfstats2_down", "#DeltaY_reco_d2_{t#bar{t}} btag_hfstats2_down",  2, -2.5, 2.5);
  DeltaY_reco_d2_btag_lf_up         = book<TH1F>("DeltaY_reco_d2_btag_lf_up", "#DeltaY_reco_d2_{t#bar{t}} btag_lf_up",                  2, -2.5, 2.5);
  DeltaY_reco_d2_btag_lf_down       = book<TH1F>("DeltaY_reco_d2_btag_lf_down", "#DeltaY_reco_d2_{t#bar{t}} btag_lf_down",              2, -2.5, 2.5);
  DeltaY_reco_d2_btag_lfstats1_up   = book<TH1F>("DeltaY_reco_d2_btag_lfstats1_up", "#DeltaY_reco_d2_{t#bar{t}} btag_lfstats1_up",      2, -2.5, 2.5);
  DeltaY_reco_d2_btag_lfstats1_down = book<TH1F>("DeltaY_reco_d2_btag_lfstats1_down", "#DeltaY_reco_d2_{t#bar{t}} btag_lfstats1_down",  2, -2.5, 2.5);
  DeltaY_reco_d2_btag_lfstats2_up   = book<TH1F>("DeltaY_reco_d2_btag_lfstats2_up", "#DeltaY_reco_d2_{t#bar{t}} btag_lfstats2_up",      2, -2.5, 2.5);
  DeltaY_reco_d2_btag_lfstats2_down = book<TH1F>("DeltaY_reco_d2_btag_lfstats2_down", "#DeltaY_reco_d2_{t#bar{t}} btag_lfstats2_down",  2, -2.5, 2.5);
  DeltaY_reco_d2_ttag_corr_up       = book<TH1F>("DeltaY_reco_d2_ttag_corr_up", "#DeltaY_reco_d2_{t#bar{t}} ttag_corr_up",              2, -2.5, 2.5);
  DeltaY_reco_d2_ttag_corr_down     = book<TH1F>("DeltaY_reco_d2_ttag_corr_down", "#DeltaY_reco_d2_{t#bar{t}} ttag_corr_down",          2, -2.5, 2.5);
  DeltaY_reco_d2_ttag_uncorr_up     = book<TH1F>("DeltaY_reco_d2_ttag_uncorr_up", "#DeltaY_reco_d2_{t#bar{t}} ttag_uncorr_up",          2, -2.5, 2.5);
  DeltaY_reco_d2_ttag_uncorr_down   = book<TH1F>("DeltaY_reco_d2_ttag_uncorr_down", "#DeltaY_reco_d2_{t#bar{t}} ttag_counrr_down",      2, -2.5, 2.5);
  DeltaY_reco_d2_tmistag_up         = book<TH1F>("DeltaY_reco_d2_tmistag_up", "#DeltaY_reco_d2_{t#bar{t}} [GeV] tmistag_up",            2, -2.5, 2.5);
  DeltaY_reco_d2_tmistag_down       = book<TH1F>("DeltaY_reco_d2_tmistag_down", "#DeltaY_reco_d2DeltaY_{t#bar{t}} [GeV] tmistag_down",        2, -2.5, 2.5);
  DeltaY_reco_d2_toppt_a_up         = book<TH1F>("DeltaY_reco_d2_toppt_a_up", "#DeltaY_reco_d2_{t#bar{t}} [GeV] toppt_a_up",                2, -2.5, 2.5);
  DeltaY_reco_d2_toppt_a_down       = book<TH1F>("DeltaY_reco_d2_toppt_a_down", "#DeltaY_reco_d2_{t#bar{t}} [GeV] toppt_a_down",            2, -2.5, 2.5);
  DeltaY_reco_d2_toppt_b_up         = book<TH1F>("DeltaY_reco_d2_toppt_b_up", "#DeltaY_reco_d2_{t#bar{t}} [GeV] toppt_b_up",                2, -2.5, 2.5);
  DeltaY_reco_d2_toppt_b_down       = book<TH1F>("DeltaY_reco_d2_toppt_b_down", "#DeltaY_reco_d2_{t#bar{t}} [GeV] toppt_b_down",            2, -2.5, 2.5);

  // // cos_theta1k_antiLep = book<TH1F>("cos_theta1k_antiLep", "cos(#theta_{antilep}^{k})",24, -1, 1);----------------------------------------------------------------------//
  // // electron systematics: reco, id, trigger
  // cos_theta1k_antiLep_ele_reco_up          = book<TH1F>("cos_theta1k_antiLep_ele_reco_up",   "cos(#theta_{antilep}^{k}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta1k_antiLep_ele_reco_down        = book<TH1F>("cos_theta1k_antiLep_ele_reco_down", "cos(#theta_{antilep}^{k}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta1k_antiLep_ele_id_up            = book<TH1F>("cos_theta1k_antiLep_ele_id_up",   "cos(#theta_{antilep}^{k}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta1k_antiLep_ele_id_down          = book<TH1F>("cos_theta1k_antiLep_ele_id_down", "cos(#theta_{antilep}^{k}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta1k_antiLep_ele_trigger_up       = book<TH1F>("cos_theta1k_antiLep_ele_trigger_up",   "cos(#theta_{antilep}^{k}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta1k_antiLep_ele_trigger_down     = book<TH1F>("cos_theta1k_antiLep_ele_trigger_down", "cos(#theta_{antilep}^{k}) ele_trigger_down",           24, -1.0, 1.0);
  // // muon systematics: reco, id_stat, id_syst, trigger_stat, trigger_syst, iso_stat, iso_syst
  // cos_theta1k_antiLep_mu_reco_up           = book<TH1F>("cos_theta1k_antiLep_mu_reco_up",   "cos(#theta_{antilep}^{k}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_reco_down         = book<TH1F>("cos_theta1k_antiLep_mu_reco_down", "cos(#theta_{antilep}^{k}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_id_stat_up        = book<TH1F>("cos_theta1k_antiLep_mu_id_stat_up",   "cos(#theta_{antilep}^{k}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_id_stat_down      = book<TH1F>("cos_theta1k_antiLep_mu_id_stat_down",   "cos(#theta_{antilep}^{k}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_id_syst_up        = book<TH1F>("cos_theta1k_antiLep_mu_id_syst_up",   "cos(#theta_{antilep}^{k}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_id_syst_down      = book<TH1F>("cos_theta1k_antiLep_mu_id_syst_down",   "cos(#theta_{antilep}^{k}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_trigger_stat_up   = book<TH1F>("cos_theta1k_antiLep_mu_trigger_stat_up",   "cos(#theta_{antilep}^{k}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_trigger_stat_down = book<TH1F>("cos_theta1k_antiLep_mu_trigger_stat_down",   "cos(#theta_{antilep}^{k}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_trigger_syst_up   = book<TH1F>("cos_theta1k_antiLep_mu_trigger_syst_up",   "cos(#theta_{antilep}^{k}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_trigger_syst_down = book<TH1F>("cos_theta1k_antiLep_mu_trigger_syst_down",   "cos(#theta_{antilep}^{k}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_iso_stat_up       = book<TH1F>("cos_theta1k_antiLep_mu_iso_stat_up",   "cos(#theta_{antilep}^{k}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_iso_stat_down     = book<TH1F>("cos_theta1k_antiLep_mu_iso_stat_down",   "cos(#theta_{antilep}^{k}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_iso_syst_up       = book<TH1F>("cos_theta1k_antiLep_mu_iso_syst_up",   "cos(#theta_{antilep}^{k}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta1k_antiLep_mu_iso_syst_down     = book<TH1F>("cos_theta1k_antiLep_mu_iso_syst_down",   "cos(#theta_{antilep}^{k}) mu_iso_syst_down",         24, -1.0, 1.0);
  // // Pileup reweighting systematics
  // cos_theta1k_antiLep_pu_up                = book<TH1F>("cos_theta1k_antiLep_pu_up",   "cos(#theta_{antilep}^{k}) pu_up",                               24, -1.0, 1.0);
  // cos_theta1k_antiLep_pu_down              = book<TH1F>("cos_theta1k_antiLep_pu_down", "cos(#theta_{antilep}^{k}) pu_down",                             24, -1.0, 1.0);
  // // Prefiring systematics
  // cos_theta1k_antiLep_prefiring_up         = book<TH1F>("cos_theta1k_antiLep_prefiring_up",   "cos(#theta_{antilep}^{k}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta1k_antiLep_prefiring_down       = book<TH1F>("cos_theta1k_antiLep_prefiring_down", "cos(#theta_{antilep}^{k}) prefiring_down",               24, -1.0, 1.0);
  // // b-tagging systematics: cferr1, cferr2, hf, hfstats1, hfstats2, lf, lfstats1, lfstats2
  // cos_theta1k_antiLep_btag_cferr1_up       = book<TH1F>("cos_theta1k_antiLep_btag_cferr1_up", "cos(#theta_{antilep}^{k}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_cferr1_down     = book<TH1F>("cos_theta1k_antiLep_btag_cferr1_down", "cos(#theta_{antilep}^{k}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_cferr2_up       = book<TH1F>("cos_theta1k_antiLep_btag_cferr2_up", "cos(#theta_{antilep}^{k}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_cferr2_down     = book<TH1F>("cos_theta1k_antiLep_btag_cferr2_down", "cos(#theta_{antilep}^{k}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_hf_up           = book<TH1F>("cos_theta1k_antiLep_btag_hf_up", "cos(#theta_{antilep}^{k}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_hf_down         = book<TH1F>("cos_theta1k_antiLep_btag_hf_down", "cos(#theta_{antilep}^{k}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_hfstats1_up     = book<TH1F>("cos_theta1k_antiLep_btag_hfstats1_up", "cos(#theta_{antilep}^{k}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_hfstats1_down   = book<TH1F>("cos_theta1k_antiLep_btag_hfstats1_down", "cos(#theta_{antilep}^{k}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_hfstats2_up     = book<TH1F>("cos_theta1k_antiLep_btag_hfstats2_up", "cos(#theta_{antilep}^{k}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_hfstats2_down   = book<TH1F>("cos_theta1k_antiLep_btag_hfstats2_down", "cos(#theta_{antilep}^{k}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_lf_up           = book<TH1F>("cos_theta1k_antiLep_btag_lf_up", "cos(#theta_{antilep}^{k}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_lf_down         = book<TH1F>("cos_theta1k_antiLep_btag_lf_down", "cos(#theta_{antilep}^{k}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_lfstats1_up     = book<TH1F>("cos_theta1k_antiLep_btag_lfstats1_up", "cos(#theta_{antilep}^{k}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_lfstats1_down   = book<TH1F>("cos_theta1k_antiLep_btag_lfstats1_down", "cos(#theta_{antilep}^{k}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_lfstats2_up     = book<TH1F>("cos_theta1k_antiLep_btag_lfstats2_up", "cos(#theta_{antilep}^{k}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta1k_antiLep_btag_lfstats2_down   = book<TH1F>("cos_theta1k_antiLep_btag_lfstats2_down", "cos(#theta_{antilep}^{k}) btag_lfstats2_down",       24, -1.0, 1.0);
  // // Top tagging systematics
  // cos_theta1k_antiLep_ttag_corr_up         = book<TH1F>("cos_theta1k_antiLep_ttag_corr_up", "cos(#theta_{antilep}^{k}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta1k_antiLep_ttag_corr_down       = book<TH1F>("cos_theta1k_antiLep_ttag_corr_down", "cos(#theta_{antilep}^{k}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta1k_antiLep_ttag_uncorr_up       = book<TH1F>("cos_theta1k_antiLep_ttag_uncorr_up", "cos(#theta_{antilep}^{k}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta1k_antiLep_ttag_uncorr_down     = book<TH1F>("cos_theta1k_antiLep_ttag_uncorr_down", "cos(#theta_{antilep}^{k}) ttag_counrr_down",           24, -1.0, 1.0);
  // // Top mistagging systematics
  // cos_theta1k_antiLep_tmistag_up           = book<TH1F>("cos_theta1k_antiLep_tmistag_up", "cos(#theta_{antilep}^{k}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta1k_antiLep_tmistag_down         = book<TH1F>("cos_theta1k_antiLep_tmistag_down", "cos(#theta_{antilep}^{k}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // // Top pT reweighting systematics
  // cos_theta1k_antiLep_toppt_a_up           = book<TH1F>("cos_theta1k_antiLep_toppt_a_up", "cos(#theta_{antilep}^{k}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta1k_antiLep_toppt_a_down         = book<TH1F>("cos_theta1k_antiLep_toppt_a_down", "cos(#theta_{antilep}^{k}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta1k_antiLep_toppt_b_up           = book<TH1F>("cos_theta1k_antiLep_toppt_b_up", "cos(#theta_{antilep}^{k}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta1k_antiLep_toppt_b_down         = book<TH1F>("cos_theta1k_antiLep_toppt_b_down", "cos(#theta_{antilep}^{k}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // // muR and muF at ME level systematics
  // cos_theta1k_antiLep_murmuf_upup          = book<TH1F>("cos_theta1k_antiLep_murmuf_upup", "cos(#theta_{antilep}^{k}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta1k_antiLep_murmuf_upnone        = book<TH1F>("cos_theta1k_antiLep_murmuf_upnone", "cos(#theta_{antilep}^{k}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta1k_antiLep_murmuf_noneup        = book<TH1F>("cos_theta1k_antiLep_murmuf_noneup", "cos(#theta_{antilep}^{k}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta1k_antiLep_murmuf_nonedown      = book<TH1F>("cos_theta1k_antiLep_murmuf_nonedown", "cos(#theta_{antilep}^{k}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta1k_antiLep_murmuf_downnone      = book<TH1F>("cos_theta1k_antiLep_murmuf_downnone", "cos(#theta_{antilep}^{k}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta1k_antiLep_murmuf_downdown      = book<TH1F>("cos_theta1k_antiLep_murmuf_downdown", "cos(#theta_{antilep}^{k}) murmuf_downdown",             24, -1.0, 1.0);
  // // ISR and FSR systematics
  // cos_theta1k_antiLep_isr_up               = book<TH1F>("cos_theta1k_antiLep_isr_up", "cos(#theta_{antilep}^{k}) isr_up",                               24, -1.0, 1.0);
  // cos_theta1k_antiLep_isr_down             = book<TH1F>("cos_theta1k_antiLep_isr_down", "cos(#theta_{antilep}^{k}) isr_down",                           24, -1.0, 1.0);
  // cos_theta1k_antiLep_fsr_up               = book<TH1F>("cos_theta1k_antiLep_fsr_up", "cos(#theta_{antilep}^{k}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta1k_antiLep_fsr_down             = book<TH1F>("cos_theta1k_antiLep_fsr_down", "cos(#theta_{antilep}^{k}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta1r_antiLep = book<TH1F>("cos_theta1r_antiLep", "cos(#theta_{antilep}^{r})",24, -1, 1);----------------------------------------------------------------------//
  // cos_theta1r_antiLep_ele_reco_up          = book<TH1F>("cos_theta1r_antiLep_ele_reco_up",   "cos(#theta_{antilep}^{r}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta1r_antiLep_ele_reco_down        = book<TH1F>("cos_theta1r_antiLep_ele_reco_down", "cos(#theta_{antilep}^{r}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta1r_antiLep_ele_id_up            = book<TH1F>("cos_theta1r_antiLep_ele_id_up",   "cos(#theta_{antilep}^{r}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta1r_antiLep_ele_id_down          = book<TH1F>("cos_theta1r_antiLep_ele_id_down", "cos(#theta_{antilep}^{r}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta1r_antiLep_ele_trigger_up       = book<TH1F>("cos_theta1r_antiLep_ele_trigger_up",   "cos(#theta_{antilep}^{r}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta1r_antiLep_ele_trigger_down     = book<TH1F>("cos_theta1r_antiLep_ele_trigger_down", "cos(#theta_{antilep}^{r}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_reco_up           = book<TH1F>("cos_theta1r_antiLep_mu_reco_up",   "cos(#theta_{antilep}^{r}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_reco_down         = book<TH1F>("cos_theta1r_antiLep_mu_reco_down", "cos(#theta_{antilep}^{r}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_id_stat_up        = book<TH1F>("cos_theta1r_antiLep_mu_id_stat_up",   "cos(#theta_{antilep}^{r}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_id_stat_down      = book<TH1F>("cos_theta1r_antiLep_mu_id_stat_down",   "cos(#theta_{antilep}^{r}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_id_syst_up        = book<TH1F>("cos_theta1r_antiLep_mu_id_syst_up",   "cos(#theta_{antilep}^{r}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_id_syst_down      = book<TH1F>("cos_theta1r_antiLep_mu_id_syst_down",   "cos(#theta_{antilep}^{r}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_trigger_stat_up   = book<TH1F>("cos_theta1r_antiLep_mu_trigger_stat_up",   "cos(#theta_{antilep}^{r}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_trigger_stat_down = book<TH1F>("cos_theta1r_antiLep_mu_trigger_stat_down",   "cos(#theta_{antilep}^{r}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_trigger_syst_up   = book<TH1F>("cos_theta1r_antiLep_mu_trigger_syst_up",   "cos(#theta_{antilep}^{r}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_trigger_syst_down = book<TH1F>("cos_theta1r_antiLep_mu_trigger_syst_down",   "cos(#theta_{antilep}^{r}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_iso_stat_up       = book<TH1F>("cos_theta1r_antiLep_mu_iso_stat_up",   "cos(#theta_{antilep}^{r}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_iso_stat_down     = book<TH1F>("cos_theta1r_antiLep_mu_iso_stat_down",   "cos(#theta_{antilep}^{r}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_iso_syst_up       = book<TH1F>("cos_theta1r_antiLep_mu_iso_syst_up",   "cos(#theta_{antilep}^{r}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta1r_antiLep_mu_iso_syst_down     = book<TH1F>("cos_theta1r_antiLep_mu_iso_syst_down",   "cos(#theta_{antilep}^{r}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta1r_antiLep_pu_up                = book<TH1F>("cos_theta1r_antiLep_pu_up",   "cos(#theta_{antilep}^{r}) pu_up",                               24, -1.0, 1.0);
  // cos_theta1r_antiLep_pu_down              = book<TH1F>("cos_theta1r_antiLep_pu_down", "cos(#theta_{antilep}^{r}) pu_down",                             24, -1.0, 1.0);
  // cos_theta1r_antiLep_prefiring_up         = book<TH1F>("cos_theta1r_antiLep_prefiring_up",   "cos(#theta_{antilep}^{r}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta1r_antiLep_prefiring_down       = book<TH1F>("cos_theta1r_antiLep_prefiring_down", "cos(#theta_{antilep}^{r}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_cferr1_up       = book<TH1F>("cos_theta1r_antiLep_btag_cferr1_up", "cos(#theta_{antilep}^{r}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_cferr1_down     = book<TH1F>("cos_theta1r_antiLep_btag_cferr1_down", "cos(#theta_{antilep}^{r}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_cferr2_up       = book<TH1F>("cos_theta1r_antiLep_btag_cferr2_up", "cos(#theta_{antilep}^{r}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_cferr2_down     = book<TH1F>("cos_theta1r_antiLep_btag_cferr2_down", "cos(#theta_{antilep}^{r}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_hf_up           = book<TH1F>("cos_theta1r_antiLep_btag_hf_up", "cos(#theta_{antilep}^{r}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_hf_down         = book<TH1F>("cos_theta1r_antiLep_btag_hf_down", "cos(#theta_{antilep}^{r}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_hfstats1_up     = book<TH1F>("cos_theta1r_antiLep_btag_hfstats1_up", "cos(#theta_{antilep}^{r}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_hfstats1_down   = book<TH1F>("cos_theta1r_antiLep_btag_hfstats1_down", "cos(#theta_{antilep}^{r}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_hfstats2_up     = book<TH1F>("cos_theta1r_antiLep_btag_hfstats2_up", "cos(#theta_{antilep}^{r}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_hfstats2_down   = book<TH1F>("cos_theta1r_antiLep_btag_hfstats2_down", "cos(#theta_{antilep}^{r}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_lf_up           = book<TH1F>("cos_theta1r_antiLep_btag_lf_up", "cos(#theta_{antilep}^{r}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_lf_down         = book<TH1F>("cos_theta1r_antiLep_btag_lf_down", "cos(#theta_{antilep}^{r}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_lfstats1_up     = book<TH1F>("cos_theta1r_antiLep_btag_lfstats1_up", "cos(#theta_{antilep}^{r}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_lfstats1_down   = book<TH1F>("cos_theta1r_antiLep_btag_lfstats1_down", "cos(#theta_{antilep}^{r}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_lfstats2_up     = book<TH1F>("cos_theta1r_antiLep_btag_lfstats2_up", "cos(#theta_{antilep}^{r}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_btag_lfstats2_down   = book<TH1F>("cos_theta1r_antiLep_btag_lfstats2_down", "cos(#theta_{antilep}^{r}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta1r_antiLep_ttag_corr_up         = book<TH1F>("cos_theta1r_antiLep_ttag_corr_up", "cos(#theta_{antilep}^{r}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta1r_antiLep_ttag_corr_down       = book<TH1F>("cos_theta1r_antiLep_ttag_corr_down", "cos(#theta_{antilep}^{r}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta1r_antiLep_ttag_uncorr_up       = book<TH1F>("cos_theta1r_antiLep_ttag_uncorr_up", "cos(#theta_{antilep}^{r}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta1r_antiLep_ttag_uncorr_down     = book<TH1F>("cos_theta1r_antiLep_ttag_uncorr_down", "cos(#theta_{antilep}^{r}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta1r_antiLep_tmistag_up           = book<TH1F>("cos_theta1r_antiLep_tmistag_up", "cos(#theta_{antilep}^{r}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta1r_antiLep_tmistag_down         = book<TH1F>("cos_theta1r_antiLep_tmistag_down", "cos(#theta_{antilep}^{r}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta1r_antiLep_toppt_a_up           = book<TH1F>("cos_theta1r_antiLep_toppt_a_up", "cos(#theta_{antilep}^{r}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta1r_antiLep_toppt_a_down         = book<TH1F>("cos_theta1r_antiLep_toppt_a_down", "cos(#theta_{antilep}^{r}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta1r_antiLep_toppt_b_up           = book<TH1F>("cos_theta1r_antiLep_toppt_b_up", "cos(#theta_{antilep}^{r}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta1r_antiLep_toppt_b_down         = book<TH1F>("cos_theta1r_antiLep_toppt_b_down", "cos(#theta_{antilep}^{r}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta1r_antiLep_murmuf_upup          = book<TH1F>("cos_theta1r_antiLep_murmuf_upup", "cos(#theta_{antilep}^{r}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta1r_antiLep_murmuf_upnone        = book<TH1F>("cos_theta1r_antiLep_murmuf_upnone", "cos(#theta_{antilep}^{r}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta1r_antiLep_murmuf_noneup        = book<TH1F>("cos_theta1r_antiLep_murmuf_noneup", "cos(#theta_{antilep}^{r}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta1r_antiLep_murmuf_nonedown      = book<TH1F>("cos_theta1r_antiLep_murmuf_nonedown", "cos(#theta_{antilep}^{r}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta1r_antiLep_murmuf_downnone      = book<TH1F>("cos_theta1r_antiLep_murmuf_downnone", "cos(#theta_{antilep}^{r}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta1r_antiLep_murmuf_downdown      = book<TH1F>("cos_theta1r_antiLep_murmuf_downdown", "cos(#theta_{antilep}^{r}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta1r_antiLep_isr_up               = book<TH1F>("cos_theta1r_antiLep_isr_up", "cos(#theta_{antilep}^{r}) isr_up",                               24, -1.0, 1.0);
  // cos_theta1r_antiLep_isr_down             = book<TH1F>("cos_theta1r_antiLep_isr_down", "cos(#theta_{antilep}^{r}) isr_down",                           24, -1.0, 1.0);
  // cos_theta1r_antiLep_fsr_up               = book<TH1F>("cos_theta1r_antiLep_fsr_up", "cos(#theta_{antilep}^{r}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta1r_antiLep_fsr_down             = book<TH1F>("cos_theta1r_antiLep_fsr_down", "cos(#theta_{antilep}^{r}) fsr_down",                           24, -1.0, 1.0);
  
  // // cos_theta1n_antiLep = book<TH1F>("cos_theta1n_antiLep", "cos(#theta_{antilep}^{n})",24, -1, 1); -----------------------------------------------------------------//
  // cos_theta1n_antiLep_ele_reco_up          = book<TH1F>("cos_theta1n_antiLep_ele_reco_up",   "cos(#theta_{antilep}^{n}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta1n_antiLep_ele_reco_down        = book<TH1F>("cos_theta1n_antiLep_ele_reco_down", "cos(#theta_{antilep}^{n}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta1n_antiLep_ele_id_up            = book<TH1F>("cos_theta1n_antiLep_ele_id_up",   "cos(#theta_{antilep}^{n}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta1n_antiLep_ele_id_down          = book<TH1F>("cos_theta1n_antiLep_ele_id_down", "cos(#theta_{antilep}^{n}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta1n_antiLep_ele_trigger_up       = book<TH1F>("cos_theta1n_antiLep_ele_trigger_up",   "cos(#theta_{antilep}^{n}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta1n_antiLep_ele_trigger_down     = book<TH1F>("cos_theta1n_antiLep_ele_trigger_down", "cos(#theta_{antilep}^{n}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_reco_up           = book<TH1F>("cos_theta1n_antiLep_mu_reco_up",   "cos(#theta_{antilep}^{n}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_reco_down         = book<TH1F>("cos_theta1n_antiLep_mu_reco_down", "cos(#theta_{antilep}^{n}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_id_stat_up        = book<TH1F>("cos_theta1n_antiLep_mu_id_stat_up",   "cos(#theta_{antilep}^{n}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_id_stat_down      = book<TH1F>("cos_theta1n_antiLep_mu_id_stat_down",   "cos(#theta_{antilep}^{n}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_id_syst_up        = book<TH1F>("cos_theta1n_antiLep_mu_id_syst_up",   "cos(#theta_{antilep}^{n}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_id_syst_down      = book<TH1F>("cos_theta1n_antiLep_mu_id_syst_down",   "cos(#theta_{antilep}^{n}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_trigger_stat_up   = book<TH1F>("cos_theta1n_antiLep_mu_trigger_stat_up",   "cos(#theta_{antilep}^{n}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_trigger_stat_down = book<TH1F>("cos_theta1n_antiLep_mu_trigger_stat_down",   "cos(#theta_{antilep}^{n}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_trigger_syst_up   = book<TH1F>("cos_theta1n_antiLep_mu_trigger_syst_up",   "cos(#theta_{antilep}^{n}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_trigger_syst_down = book<TH1F>("cos_theta1n_antiLep_mu_trigger_syst_down",   "cos(#theta_{antilep}^{n}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_iso_stat_up       = book<TH1F>("cos_theta1n_antiLep_mu_iso_stat_up",   "cos(#theta_{antilep}^{n}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_iso_stat_down     = book<TH1F>("cos_theta1n_antiLep_mu_iso_stat_down",   "cos(#theta_{antilep}^{n}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_iso_syst_up       = book<TH1F>("cos_theta1n_antiLep_mu_iso_syst_up",   "cos(#theta_{antilep}^{n}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta1n_antiLep_mu_iso_syst_down     = book<TH1F>("cos_theta1n_antiLep_mu_iso_syst_down",   "cos(#theta_{antilep}^{n}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta1n_antiLep_pu_up                = book<TH1F>("cos_theta1n_antiLep_pu_up",   "cos(#theta_{antilep}^{n}) pu_up",                               24, -1.0, 1.0);
  // cos_theta1n_antiLep_pu_down              = book<TH1F>("cos_theta1n_antiLep_pu_down", "cos(#theta_{antilep}^{n}) pu_down",                             24, -1.0, 1.0);
  // cos_theta1n_antiLep_prefiring_up         = book<TH1F>("cos_theta1n_antiLep_prefiring_up",   "cos(#theta_{antilep}^{n}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta1n_antiLep_prefiring_down       = book<TH1F>("cos_theta1n_antiLep_prefiring_down", "cos(#theta_{antilep}^{n}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_cferr1_up       = book<TH1F>("cos_theta1n_antiLep_btag_cferr1_up", "cos(#theta_{antilep}^{n}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_cferr1_down     = book<TH1F>("cos_theta1n_antiLep_btag_cferr1_down", "cos(#theta_{antilep}^{n}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_cferr2_up       = book<TH1F>("cos_theta1n_antiLep_btag_cferr2_up", "cos(#theta_{antilep}^{n}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_cferr2_down     = book<TH1F>("cos_theta1n_antiLep_btag_cferr2_down", "cos(#theta_{antilep}^{n}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_hf_up           = book<TH1F>("cos_theta1n_antiLep_btag_hf_up", "cos(#theta_{antilep}^{n}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_hf_down         = book<TH1F>("cos_theta1n_antiLep_btag_hf_down", "cos(#theta_{antilep}^{n}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_hfstats1_up     = book<TH1F>("cos_theta1n_antiLep_btag_hfstats1_up", "cos(#theta_{antilep}^{n}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_hfstats1_down   = book<TH1F>("cos_theta1n_antiLep_btag_hfstats1_down", "cos(#theta_{antilep}^{n}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_hfstats2_up     = book<TH1F>("cos_theta1n_antiLep_btag_hfstats2_up", "cos(#theta_{antilep}^{n}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_hfstats2_down   = book<TH1F>("cos_theta1n_antiLep_btag_hfstats2_down", "cos(#theta_{antilep}^{n}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_lf_up           = book<TH1F>("cos_theta1n_antiLep_btag_lf_up", "cos(#theta_{antilep}^{n}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_lf_down         = book<TH1F>("cos_theta1n_antiLep_btag_lf_down", "cos(#theta_{antilep}^{n}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_lfstats1_up     = book<TH1F>("cos_theta1n_antiLep_btag_lfstats1_up", "cos(#theta_{antilep}^{n}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_lfstats1_down   = book<TH1F>("cos_theta1n_antiLep_btag_lfstats1_down", "cos(#theta_{antilep}^{n}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_lfstats2_up     = book<TH1F>("cos_theta1n_antiLep_btag_lfstats2_up", "cos(#theta_{antilep}^{n}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_btag_lfstats2_down   = book<TH1F>("cos_theta1n_antiLep_btag_lfstats2_down", "cos(#theta_{antilep}^{n}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta1n_antiLep_ttag_corr_up         = book<TH1F>("cos_theta1n_antiLep_ttag_corr_up", "cos(#theta_{antilep}^{n}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta1n_antiLep_ttag_corr_down       = book<TH1F>("cos_theta1n_antiLep_ttag_corr_down", "cos(#theta_{antilep}^{n}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta1n_antiLep_ttag_uncorr_up       = book<TH1F>("cos_theta1n_antiLep_ttag_uncorr_up", "cos(#theta_{antilep}^{n}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta1n_antiLep_ttag_uncorr_down     = book<TH1F>("cos_theta1n_antiLep_ttag_uncorr_down", "cos(#theta_{antilep}^{n}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta1n_antiLep_tmistag_up           = book<TH1F>("cos_theta1n_antiLep_tmistag_up", "cos(#theta_{antilep}^{n}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta1n_antiLep_tmistag_down         = book<TH1F>("cos_theta1n_antiLep_tmistag_down", "cos(#theta_{antilep}^{n}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta1n_antiLep_toppt_a_up           = book<TH1F>("cos_theta1n_antiLep_toppt_a_up", "cos(#theta_{antilep}^{n}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta1n_antiLep_toppt_a_down         = book<TH1F>("cos_theta1n_antiLep_toppt_a_down", "cos(#theta_{antilep}^{n}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta1n_antiLep_toppt_b_up           = book<TH1F>("cos_theta1n_antiLep_toppt_b_up", "cos(#theta_{antilep}^{n}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta1n_antiLep_toppt_b_down         = book<TH1F>("cos_theta1n_antiLep_toppt_b_down", "cos(#theta_{antilep}^{n}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta1n_antiLep_murmuf_upup          = book<TH1F>("cos_theta1n_antiLep_murmuf_upup", "cos(#theta_{antilep}^{n}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta1n_antiLep_murmuf_upnone        = book<TH1F>("cos_theta1n_antiLep_murmuf_upnone", "cos(#theta_{antilep}^{n}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta1n_antiLep_murmuf_noneup        = book<TH1F>("cos_theta1n_antiLep_murmuf_noneup", "cos(#theta_{antilep}^{n}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta1n_antiLep_murmuf_nonedown      = book<TH1F>("cos_theta1n_antiLep_murmuf_nonedown", "cos(#theta_{antilep}^{n}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta1n_antiLep_murmuf_downnone      = book<TH1F>("cos_theta1n_antiLep_murmuf_downnone", "cos(#theta_{antilep}^{n}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta1n_antiLep_murmuf_downdown      = book<TH1F>("cos_theta1n_antiLep_murmuf_downdown", "cos(#theta_{antilep}^{n}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta1n_antiLep_isr_up               = book<TH1F>("cos_theta1n_antiLep_isr_up", "cos(#theta_{antilep}^{n}) isr_up",                               24, -1.0, 1.0);
  // cos_theta1n_antiLep_isr_down             = book<TH1F>("cos_theta1n_antiLep_isr_down", "cos(#theta_{antilep}^{n}) isr_down",                           24, -1.0, 1.0);
  // cos_theta1n_antiLep_fsr_up               = book<TH1F>("cos_theta1n_antiLep_fsr_up", "cos(#theta_{antilep}^{n}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta1n_antiLep_fsr_down             = book<TH1F>("cos_theta1n_antiLep_fsr_down", "cos(#theta_{antilep}^{n}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta1kStar_antiLep = book<TH1F>("cos_theta1kStar_antiLep", "cos(#theta_{antilep}^{k*})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta1kStar_antiLep_ele_reco_up          = book<TH1F>("cos_theta1kStar_antiLep_ele_reco_up",   "cos(#theta_{antilep}^{k*}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_ele_reco_down        = book<TH1F>("cos_theta1kStar_antiLep_ele_reco_down", "cos(#theta_{antilep}^{k*}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_ele_id_up            = book<TH1F>("cos_theta1kStar_antiLep_ele_id_up",   "cos(#theta_{antilep}^{k*}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_ele_id_down          = book<TH1F>("cos_theta1kStar_antiLep_ele_id_down", "cos(#theta_{antilep}^{k*}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_ele_trigger_up       = book<TH1F>("cos_theta1kStar_antiLep_ele_trigger_up",   "cos(#theta_{antilep}^{k*}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_ele_trigger_down     = book<TH1F>("cos_theta1kStar_antiLep_ele_trigger_down", "cos(#theta_{antilep}^{k*}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_reco_up           = book<TH1F>("cos_theta1kStar_antiLep_mu_reco_up",   "cos(#theta_{antilep}^{k*}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_reco_down         = book<TH1F>("cos_theta1kStar_antiLep_mu_reco_down", "cos(#theta_{antilep}^{k*}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_id_stat_up        = book<TH1F>("cos_theta1kStar_antiLep_mu_id_stat_up",   "cos(#theta_{antilep}^{k*}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_id_stat_down      = book<TH1F>("cos_theta1kStar_antiLep_mu_id_stat_down",   "cos(#theta_{antilep}^{k*}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_id_syst_up        = book<TH1F>("cos_theta1kStar_antiLep_mu_id_syst_up",   "cos(#theta_{antilep}^{k*}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_id_syst_down      = book<TH1F>("cos_theta1kStar_antiLep_mu_id_syst_down",   "cos(#theta_{antilep}^{k*}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_trigger_stat_up   = book<TH1F>("cos_theta1kStar_antiLep_mu_trigger_stat_up",   "cos(#theta_{antilep}^{k*}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_trigger_stat_down = book<TH1F>("cos_theta1kStar_antiLep_mu_trigger_stat_down",   "cos(#theta_{antilep}^{k*}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_trigger_syst_up   = book<TH1F>("cos_theta1kStar_antiLep_mu_trigger_syst_up",   "cos(#theta_{antilep}^{k*}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_trigger_syst_down = book<TH1F>("cos_theta1kStar_antiLep_mu_trigger_syst_down",   "cos(#theta_{antilep}^{k*}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_iso_stat_up       = book<TH1F>("cos_theta1kStar_antiLep_mu_iso_stat_up",   "cos(#theta_{antilep}^{k*}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_iso_stat_down     = book<TH1F>("cos_theta1kStar_antiLep_mu_iso_stat_down",   "cos(#theta_{antilep}^{k*}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_iso_syst_up       = book<TH1F>("cos_theta1kStar_antiLep_mu_iso_syst_up",   "cos(#theta_{antilep}^{k*}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_mu_iso_syst_down     = book<TH1F>("cos_theta1kStar_antiLep_mu_iso_syst_down",   "cos(#theta_{antilep}^{k*}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_pu_up                = book<TH1F>("cos_theta1kStar_antiLep_pu_up",   "cos(#theta_{antilep}^{k*}) pu_up",                               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_pu_down              = book<TH1F>("cos_theta1kStar_antiLep_pu_down", "cos(#theta_{antilep}^{k*}) pu_down",                             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_prefiring_up         = book<TH1F>("cos_theta1kStar_antiLep_prefiring_up",   "cos(#theta_{antilep}^{k*}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_prefiring_down       = book<TH1F>("cos_theta1kStar_antiLep_prefiring_down", "cos(#theta_{antilep}^{k*}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_cferr1_up       = book<TH1F>("cos_theta1kStar_antiLep_btag_cferr1_up", "cos(#theta_{antilep}^{k*}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_cferr1_down     = book<TH1F>("cos_theta1kStar_antiLep_btag_cferr1_down", "cos(#theta_{antilep}^{k*}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_cferr2_up       = book<TH1F>("cos_theta1kStar_antiLep_btag_cferr2_up", "cos(#theta_{antilep}^{k*}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_cferr2_down     = book<TH1F>("cos_theta1kStar_antiLep_btag_cferr2_down", "cos(#theta_{antilep}^{k*}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_hf_up           = book<TH1F>("cos_theta1kStar_antiLep_btag_hf_up", "cos(#theta_{antilep}^{k*}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_hf_down         = book<TH1F>("cos_theta1kStar_antiLep_btag_hf_down", "cos(#theta_{antilep}^{k*}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_hfstats1_up     = book<TH1F>("cos_theta1kStar_antiLep_btag_hfstats1_up", "cos(#theta_{antilep}^{k*}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_hfstats1_down   = book<TH1F>("cos_theta1kStar_antiLep_btag_hfstats1_down", "cos(#theta_{antilep}^{k*}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_hfstats2_up     = book<TH1F>("cos_theta1kStar_antiLep_btag_hfstats2_up", "cos(#theta_{antilep}^{k*}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_hfstats2_down   = book<TH1F>("cos_theta1kStar_antiLep_btag_hfstats2_down", "cos(#theta_{antilep}^{k*}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_lf_up           = book<TH1F>("cos_theta1kStar_antiLep_btag_lf_up", "cos(#theta_{antilep}^{k*}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_lf_down         = book<TH1F>("cos_theta1kStar_antiLep_btag_lf_down", "cos(#theta_{antilep}^{k*}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_lfstats1_up     = book<TH1F>("cos_theta1kStar_antiLep_btag_lfstats1_up", "cos(#theta_{antilep}^{k*}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_lfstats1_down   = book<TH1F>("cos_theta1kStar_antiLep_btag_lfstats1_down", "cos(#theta_{antilep}^{k*}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_lfstats2_up     = book<TH1F>("cos_theta1kStar_antiLep_btag_lfstats2_up", "cos(#theta_{antilep}^{k*}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_btag_lfstats2_down   = book<TH1F>("cos_theta1kStar_antiLep_btag_lfstats2_down", "cos(#theta_{antilep}^{k*}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_ttag_corr_up         = book<TH1F>("cos_theta1kStar_antiLep_ttag_corr_up", "cos(#theta_{antilep}^{k*}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_ttag_corr_down       = book<TH1F>("cos_theta1kStar_antiLep_ttag_corr_down", "cos(#theta_{antilep}^{k*}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_ttag_uncorr_up       = book<TH1F>("cos_theta1kStar_antiLep_ttag_uncorr_up", "cos(#theta_{antilep}^{k*}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_ttag_uncorr_down     = book<TH1F>("cos_theta1kStar_antiLep_ttag_uncorr_down", "cos(#theta_{antilep}^{k*}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_tmistag_up           = book<TH1F>("cos_theta1kStar_antiLep_tmistag_up", "cos(#theta_{antilep}^{k*}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_tmistag_down         = book<TH1F>("cos_theta1kStar_antiLep_tmistag_down", "cos(#theta_{antilep}^{k*}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_toppt_a_up           = book<TH1F>("cos_theta1kStar_antiLep_toppt_a_up", "cos(#theta_{antilep}^{k*}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_toppt_a_down         = book<TH1F>("cos_theta1kStar_antiLep_toppt_a_down", "cos(#theta_{antilep}^{k*}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_toppt_b_up           = book<TH1F>("cos_theta1kStar_antiLep_toppt_b_up", "cos(#theta_{antilep}^{k*}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_toppt_b_down         = book<TH1F>("cos_theta1kStar_antiLep_toppt_b_down", "cos(#theta_{antilep}^{k*}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_murmuf_upup          = book<TH1F>("cos_theta1kStar_antiLep_murmuf_upup", "cos(#theta_{antilep}^{k*}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_murmuf_upnone        = book<TH1F>("cos_theta1kStar_antiLep_murmuf_upnone", "cos(#theta_{antilep}^{k*}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_murmuf_noneup        = book<TH1F>("cos_theta1kStar_antiLep_murmuf_noneup", "cos(#theta_{antilep}^{k*}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_murmuf_nonedown      = book<TH1F>("cos_theta1kStar_antiLep_murmuf_nonedown", "cos(#theta_{antilep}^{k*}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_murmuf_downnone      = book<TH1F>("cos_theta1kStar_antiLep_murmuf_downnone", "cos(#theta_{antilep}^{k*}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_murmuf_downdown      = book<TH1F>("cos_theta1kStar_antiLep_murmuf_downdown", "cos(#theta_{antilep}^{k*}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_isr_up               = book<TH1F>("cos_theta1kStar_antiLep_isr_up", "cos(#theta_{antilep}^{k*}) isr_up",                               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_isr_down             = book<TH1F>("cos_theta1kStar_antiLep_isr_down", "cos(#theta_{antilep}^{k*}) isr_down",                           24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_fsr_up               = book<TH1F>("cos_theta1kStar_antiLep_fsr_up", "cos(#theta_{antilep}^{k*}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta1kStar_antiLep_fsr_down             = book<TH1F>("cos_theta1kStar_antiLep_fsr_down", "cos(#theta_{antilep}^{k*}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta1rStar_antiLep = book<TH1F>("cos_theta1rStar_antiLep", "cos(#theta_{antilep}^{r*})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta1rStar_antiLep_ele_reco_up          = book<TH1F>("cos_theta1rStar_antiLep_ele_reco_up",   "cos(#theta_{antilep}^{r*}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_ele_reco_down        = book<TH1F>("cos_theta1rStar_antiLep_ele_reco_down", "cos(#theta_{antilep}^{r*}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_ele_id_up            = book<TH1F>("cos_theta1rStar_antiLep_ele_id_up",   "cos(#theta_{antilep}^{r*}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_ele_id_down          = book<TH1F>("cos_theta1rStar_antiLep_ele_id_down", "cos(#theta_{antilep}^{r*}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_ele_trigger_up       = book<TH1F>("cos_theta1rStar_antiLep_ele_trigger_up",   "cos(#theta_{antilep}^{r*}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_ele_trigger_down     = book<TH1F>("cos_theta1rStar_antiLep_ele_trigger_down", "cos(#theta_{antilep}^{r*}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_reco_up           = book<TH1F>("cos_theta1rStar_antiLep_mu_reco_up",   "cos(#theta_{antilep}^{r*}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_reco_down         = book<TH1F>("cos_theta1rStar_antiLep_mu_reco_down", "cos(#theta_{antilep}^{r*}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_id_stat_up        = book<TH1F>("cos_theta1rStar_antiLep_mu_id_stat_up",   "cos(#theta_{antilep}^{r*}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_id_stat_down      = book<TH1F>("cos_theta1rStar_antiLep_mu_id_stat_down",   "cos(#theta_{antilep}^{r*}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_id_syst_up        = book<TH1F>("cos_theta1rStar_antiLep_mu_id_syst_up",   "cos(#theta_{antilep}^{r*}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_id_syst_down      = book<TH1F>("cos_theta1rStar_antiLep_mu_id_syst_down",   "cos(#theta_{antilep}^{r*}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_trigger_stat_up   = book<TH1F>("cos_theta1rStar_antiLep_mu_trigger_stat_up",   "cos(#theta_{antilep}^{r*}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_trigger_stat_down = book<TH1F>("cos_theta1rStar_antiLep_mu_trigger_stat_down",   "cos(#theta_{antilep}^{r*}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_trigger_syst_up   = book<TH1F>("cos_theta1rStar_antiLep_mu_trigger_syst_up",   "cos(#theta_{antilep}^{r*}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_trigger_syst_down = book<TH1F>("cos_theta1rStar_antiLep_mu_trigger_syst_down",   "cos(#theta_{antilep}^{r*}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_iso_stat_up       = book<TH1F>("cos_theta1rStar_antiLep_mu_iso_stat_up",   "cos(#theta_{antilep}^{r*}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_iso_stat_down     = book<TH1F>("cos_theta1rStar_antiLep_mu_iso_stat_down",   "cos(#theta_{antilep}^{r*}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_iso_syst_up       = book<TH1F>("cos_theta1rStar_antiLep_mu_iso_syst_up",   "cos(#theta_{antilep}^{r*}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_mu_iso_syst_down     = book<TH1F>("cos_theta1rStar_antiLep_mu_iso_syst_down",   "cos(#theta_{antilep}^{r*}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_pu_up                = book<TH1F>("cos_theta1rStar_antiLep_pu_up",   "cos(#theta_{antilep}^{r*}) pu_up",                               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_pu_down              = book<TH1F>("cos_theta1rStar_antiLep_pu_down", "cos(#theta_{antilep}^{r*}) pu_down",                             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_prefiring_up         = book<TH1F>("cos_theta1rStar_antiLep_prefiring_up",   "cos(#theta_{antilep}^{r*}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_prefiring_down       = book<TH1F>("cos_theta1rStar_antiLep_prefiring_down", "cos(#theta_{antilep}^{r*}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_cferr1_up       = book<TH1F>("cos_theta1rStar_antiLep_btag_cferr1_up", "cos(#theta_{antilep}^{r*}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_cferr1_down     = book<TH1F>("cos_theta1rStar_antiLep_btag_cferr1_down", "cos(#theta_{antilep}^{r*}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_cferr2_up       = book<TH1F>("cos_theta1rStar_antiLep_btag_cferr2_up", "cos(#theta_{antilep}^{r*}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_cferr2_down     = book<TH1F>("cos_theta1rStar_antiLep_btag_cferr2_down", "cos(#theta_{antilep}^{r*}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_hf_up           = book<TH1F>("cos_theta1rStar_antiLep_btag_hf_up", "cos(#theta_{antilep}^{r*}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_hf_down         = book<TH1F>("cos_theta1rStar_antiLep_btag_hf_down", "cos(#theta_{antilep}^{r*}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_hfstats1_up     = book<TH1F>("cos_theta1rStar_antiLep_btag_hfstats1_up", "cos(#theta_{antilep}^{r*}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_hfstats1_down   = book<TH1F>("cos_theta1rStar_antiLep_btag_hfstats1_down", "cos(#theta_{antilep}^{r*}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_hfstats2_up     = book<TH1F>("cos_theta1rStar_antiLep_btag_hfstats2_up", "cos(#theta_{antilep}^{r*}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_hfstats2_down   = book<TH1F>("cos_theta1rStar_antiLep_btag_hfstats2_down", "cos(#theta_{antilep}^{r*}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_lf_up           = book<TH1F>("cos_theta1rStar_antiLep_btag_lf_up", "cos(#theta_{antilep}^{r*}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_lf_down         = book<TH1F>("cos_theta1rStar_antiLep_btag_lf_down", "cos(#theta_{antilep}^{r*}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_lfstats1_up     = book<TH1F>("cos_theta1rStar_antiLep_btag_lfstats1_up", "cos(#theta_{antilep}^{r*}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_lfstats1_down   = book<TH1F>("cos_theta1rStar_antiLep_btag_lfstats1_down", "cos(#theta_{antilep}^{r*}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_lfstats2_up     = book<TH1F>("cos_theta1rStar_antiLep_btag_lfstats2_up", "cos(#theta_{antilep}^{r*}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_btag_lfstats2_down   = book<TH1F>("cos_theta1rStar_antiLep_btag_lfstats2_down", "cos(#theta_{antilep}^{r*}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_ttag_corr_up         = book<TH1F>("cos_theta1rStar_antiLep_ttag_corr_up", "cos(#theta_{antilep}^{r*}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_ttag_corr_down       = book<TH1F>("cos_theta1rStar_antiLep_ttag_corr_down", "cos(#theta_{antilep}^{r*}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_ttag_uncorr_up       = book<TH1F>("cos_theta1rStar_antiLep_ttag_uncorr_up", "cos(#theta_{antilep}^{r*}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_ttag_uncorr_down     = book<TH1F>("cos_theta1rStar_antiLep_ttag_uncorr_down", "cos(#theta_{antilep}^{r*}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_tmistag_up           = book<TH1F>("cos_theta1rStar_antiLep_tmistag_up", "cos(#theta_{antilep}^{r*}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_tmistag_down         = book<TH1F>("cos_theta1rStar_antiLep_tmistag_down", "cos(#theta_{antilep}^{r*}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_toppt_a_up           = book<TH1F>("cos_theta1rStar_antiLep_toppt_a_up", "cos(#theta_{antilep}^{r*}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_toppt_a_down         = book<TH1F>("cos_theta1rStar_antiLep_toppt_a_down", "cos(#theta_{antilep}^{r*}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_toppt_b_up           = book<TH1F>("cos_theta1rStar_antiLep_toppt_b_up", "cos(#theta_{antilep}^{r*}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_toppt_b_down         = book<TH1F>("cos_theta1rStar_antiLep_toppt_b_down", "cos(#theta_{antilep}^{r*}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_murmuf_upup          = book<TH1F>("cos_theta1rStar_antiLep_murmuf_upup", "cos(#theta_{antilep}^{r*}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_murmuf_upnone        = book<TH1F>("cos_theta1rStar_antiLep_murmuf_upnone", "cos(#theta_{antilep}^{r*}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_murmuf_noneup        = book<TH1F>("cos_theta1rStar_antiLep_murmuf_noneup", "cos(#theta_{antilep}^{r*}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_murmuf_nonedown      = book<TH1F>("cos_theta1rStar_antiLep_murmuf_nonedown", "cos(#theta_{antilep}^{r*}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_murmuf_downnone      = book<TH1F>("cos_theta1rStar_antiLep_murmuf_downnone", "cos(#theta_{antilep}^{r*}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_murmuf_downdown      = book<TH1F>("cos_theta1rStar_antiLep_murmuf_downdown", "cos(#theta_{antilep}^{r*}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_isr_up               = book<TH1F>("cos_theta1rStar_antiLep_isr_up", "cos(#theta_{antilep}^{r*}) isr_up",                               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_isr_down             = book<TH1F>("cos_theta1rStar_antiLep_isr_down", "cos(#theta_{antilep}^{r*}) isr_down",                           24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_fsr_up               = book<TH1F>("cos_theta1rStar_antiLep_fsr_up", "cos(#theta_{antilep}^{r*}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta1rStar_antiLep_fsr_down             = book<TH1F>("cos_theta1rStar_antiLep_fsr_down", "cos(#theta_{antilep}^{r*}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta2k_Lep = book<TH1F>("cos_theta2k_Lep", "cos(#theta_{lep}^{k})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta2k_Lep_ele_reco_up          = book<TH1F>("cos_theta2k_Lep_ele_reco_up",   "cos(#theta_{lep}^{k}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta2k_Lep_ele_reco_down        = book<TH1F>("cos_theta2k_Lep_ele_reco_down", "cos(#theta_{lep}^{k}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta2k_Lep_ele_id_up            = book<TH1F>("cos_theta2k_Lep_ele_id_up",   "cos(#theta_{lep}^{k}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta2k_Lep_ele_id_down          = book<TH1F>("cos_theta2k_Lep_ele_id_down", "cos(#theta_{lep}^{k}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta2k_Lep_ele_trigger_up       = book<TH1F>("cos_theta2k_Lep_ele_trigger_up",   "cos(#theta_{lep}^{k}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta2k_Lep_ele_trigger_down     = book<TH1F>("cos_theta2k_Lep_ele_trigger_down", "cos(#theta_{lep}^{k}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_reco_up           = book<TH1F>("cos_theta2k_Lep_mu_reco_up",   "cos(#theta_{lep}^{k}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_reco_down         = book<TH1F>("cos_theta2k_Lep_mu_reco_down", "cos(#theta_{lep}^{k}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_id_stat_up        = book<TH1F>("cos_theta2k_Lep_mu_id_stat_up",   "cos(#theta_{lep}^{k}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_id_stat_down      = book<TH1F>("cos_theta2k_Lep_mu_id_stat_down",   "cos(#theta_{lep}^{k}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_id_syst_up        = book<TH1F>("cos_theta2k_Lep_mu_id_syst_up",   "cos(#theta_{lep}^{k}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_id_syst_down      = book<TH1F>("cos_theta2k_Lep_mu_id_syst_down",   "cos(#theta_{lep}^{k}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_trigger_stat_up   = book<TH1F>("cos_theta2k_Lep_mu_trigger_stat_up",   "cos(#theta_{lep}^{k}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_trigger_stat_down = book<TH1F>("cos_theta2k_Lep_mu_trigger_stat_down",   "cos(#theta_{lep}^{k}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_trigger_syst_up   = book<TH1F>("cos_theta2k_Lep_mu_trigger_syst_up",   "cos(#theta_{lep}^{k}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_trigger_syst_down = book<TH1F>("cos_theta2k_Lep_mu_trigger_syst_down",   "cos(#theta_{lep}^{k}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_iso_stat_up       = book<TH1F>("cos_theta2k_Lep_mu_iso_stat_up",   "cos(#theta_{lep}^{k}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_iso_stat_down     = book<TH1F>("cos_theta2k_Lep_mu_iso_stat_down",   "cos(#theta_{lep}^{k}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_iso_syst_up       = book<TH1F>("cos_theta2k_Lep_mu_iso_syst_up",   "cos(#theta_{lep}^{k}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta2k_Lep_mu_iso_syst_down     = book<TH1F>("cos_theta2k_Lep_mu_iso_syst_down",   "cos(#theta_{lep}^{k}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta2k_Lep_pu_up                = book<TH1F>("cos_theta2k_Lep_pu_up",   "cos(#theta_{lep}^{k}) pu_up",                               24, -1.0, 1.0);
  // cos_theta2k_Lep_pu_down              = book<TH1F>("cos_theta2k_Lep_pu_down", "cos(#theta_{lep}^{k}) pu_down",                             24, -1.0, 1.0);
  // cos_theta2k_Lep_prefiring_up         = book<TH1F>("cos_theta2k_Lep_prefiring_up",   "cos(#theta_{lep}^{k}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta2k_Lep_prefiring_down       = book<TH1F>("cos_theta2k_Lep_prefiring_down", "cos(#theta_{lep}^{k}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_cferr1_up       = book<TH1F>("cos_theta2k_Lep_btag_cferr1_up", "cos(#theta_{lep}^{k}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_cferr1_down     = book<TH1F>("cos_theta2k_Lep_btag_cferr1_down", "cos(#theta_{lep}^{k}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_cferr2_up       = book<TH1F>("cos_theta2k_Lep_btag_cferr2_up", "cos(#theta_{lep}^{k}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_cferr2_down     = book<TH1F>("cos_theta2k_Lep_btag_cferr2_down", "cos(#theta_{lep}^{k}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_hf_up           = book<TH1F>("cos_theta2k_Lep_btag_hf_up", "cos(#theta_{lep}^{k}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_hf_down         = book<TH1F>("cos_theta2k_Lep_btag_hf_down", "cos(#theta_{lep}^{k}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_hfstats1_up     = book<TH1F>("cos_theta2k_Lep_btag_hfstats1_up", "cos(#theta_{lep}^{k}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_hfstats1_down   = book<TH1F>("cos_theta2k_Lep_btag_hfstats1_down", "cos(#theta_{lep}^{k}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_hfstats2_up     = book<TH1F>("cos_theta2k_Lep_btag_hfstats2_up", "cos(#theta_{lep}^{k}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_hfstats2_down   = book<TH1F>("cos_theta2k_Lep_btag_hfstats2_down", "cos(#theta_{lep}^{k}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_lf_up           = book<TH1F>("cos_theta2k_Lep_btag_lf_up", "cos(#theta_{lep}^{k}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_lf_down         = book<TH1F>("cos_theta2k_Lep_btag_lf_down", "cos(#theta_{lep}^{k}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_lfstats1_up     = book<TH1F>("cos_theta2k_Lep_btag_lfstats1_up", "cos(#theta_{lep}^{k}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_lfstats1_down   = book<TH1F>("cos_theta2k_Lep_btag_lfstats1_down", "cos(#theta_{lep}^{k}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_lfstats2_up     = book<TH1F>("cos_theta2k_Lep_btag_lfstats2_up", "cos(#theta_{lep}^{k}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta2k_Lep_btag_lfstats2_down   = book<TH1F>("cos_theta2k_Lep_btag_lfstats2_down", "cos(#theta_{lep}^{k}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta2k_Lep_ttag_corr_up         = book<TH1F>("cos_theta2k_Lep_ttag_corr_up", "cos(#theta_{lep}^{k}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta2k_Lep_ttag_corr_down       = book<TH1F>("cos_theta2k_Lep_ttag_corr_down", "cos(#theta_{lep}^{k}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta2k_Lep_ttag_uncorr_up       = book<TH1F>("cos_theta2k_Lep_ttag_uncorr_up", "cos(#theta_{lep}^{k}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta2k_Lep_ttag_uncorr_down     = book<TH1F>("cos_theta2k_Lep_ttag_uncorr_down", "cos(#theta_{lep}^{k}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta2k_Lep_tmistag_up           = book<TH1F>("cos_theta2k_Lep_tmistag_up", "cos(#theta_{lep}^{k}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta2k_Lep_tmistag_down         = book<TH1F>("cos_theta2k_Lep_tmistag_down", "cos(#theta_{lep}^{k}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta2k_Lep_toppt_a_up           = book<TH1F>("cos_theta2k_Lep_toppt_a_up", "cos(#theta_{lep}^{k}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta2k_Lep_toppt_a_down         = book<TH1F>("cos_theta2k_Lep_toppt_a_down", "cos(#theta_{lep}^{k}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta2k_Lep_toppt_b_up           = book<TH1F>("cos_theta2k_Lep_toppt_b_up", "cos(#theta_{lep}^{k}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta2k_Lep_toppt_b_down         = book<TH1F>("cos_theta2k_Lep_toppt_b_down", "cos(#theta_{lep}^{k}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta2k_Lep_murmuf_upup          = book<TH1F>("cos_theta2k_Lep_murmuf_upup", "cos(#theta_{lep}^{k}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta2k_Lep_murmuf_upnone        = book<TH1F>("cos_theta2k_Lep_murmuf_upnone", "cos(#theta_{lep}^{k}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta2k_Lep_murmuf_noneup        = book<TH1F>("cos_theta2k_Lep_murmuf_noneup", "cos(#theta_{lep}^{k}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta2k_Lep_murmuf_nonedown      = book<TH1F>("cos_theta2k_Lep_murmuf_nonedown", "cos(#theta_{lep}^{k}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta2k_Lep_murmuf_downnone      = book<TH1F>("cos_theta2k_Lep_murmuf_downnone", "cos(#theta_{lep}^{k}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta2k_Lep_murmuf_downdown      = book<TH1F>("cos_theta2k_Lep_murmuf_downdown", "cos(#theta_{lep}^{k}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta2k_Lep_isr_up               = book<TH1F>("cos_theta2k_Lep_isr_up", "cos(#theta_{lep}^{k}) isr_up",                               24, -1.0, 1.0);
  // cos_theta2k_Lep_isr_down             = book<TH1F>("cos_theta2k_Lep_isr_down", "cos(#theta_{lep}^{k}) isr_down",                           24, -1.0, 1.0);
  // cos_theta2k_Lep_fsr_up               = book<TH1F>("cos_theta2k_Lep_fsr_up", "cos(#theta_{lep}^{k}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta2k_Lep_fsr_down             = book<TH1F>("cos_theta2k_Lep_fsr_down", "cos(#theta_{lep}^{k}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta2r_Lep = book<TH1F>("cos_theta2r_Lep", "cos(#theta_{lep}^{r})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta2r_Lep_ele_reco_up          = book<TH1F>("cos_theta2r_Lep_ele_reco_up",   "cos(#theta_{lep}^{r}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta2r_Lep_ele_reco_down        = book<TH1F>("cos_theta2r_Lep_ele_reco_down", "cos(#theta_{lep}^{r}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta2r_Lep_ele_id_up            = book<TH1F>("cos_theta2r_Lep_ele_id_up",   "cos(#theta_{lep}^{r}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta2r_Lep_ele_id_down          = book<TH1F>("cos_theta2r_Lep_ele_id_down", "cos(#theta_{lep}^{r}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta2r_Lep_ele_trigger_up       = book<TH1F>("cos_theta2r_Lep_ele_trigger_up",   "cos(#theta_{lep}^{r}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta2r_Lep_ele_trigger_down     = book<TH1F>("cos_theta2r_Lep_ele_trigger_down", "cos(#theta_{lep}^{r}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_reco_up           = book<TH1F>("cos_theta2r_Lep_mu_reco_up",   "cos(#theta_{lep}^{r}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_reco_down         = book<TH1F>("cos_theta2r_Lep_mu_reco_down", "cos(#theta_{lep}^{r}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_id_stat_up        = book<TH1F>("cos_theta2r_Lep_mu_id_stat_up",   "cos(#theta_{lep}^{r}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_id_stat_down      = book<TH1F>("cos_theta2r_Lep_mu_id_stat_down",   "cos(#theta_{lep}^{r}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_id_syst_up        = book<TH1F>("cos_theta2r_Lep_mu_id_syst_up",   "cos(#theta_{lep}^{r}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_id_syst_down      = book<TH1F>("cos_theta2r_Lep_mu_id_syst_down",   "cos(#theta_{lep}^{r}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_trigger_stat_up   = book<TH1F>("cos_theta2r_Lep_mu_trigger_stat_up",   "cos(#theta_{lep}^{r}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_trigger_stat_down = book<TH1F>("cos_theta2r_Lep_mu_trigger_stat_down",   "cos(#theta_{lep}^{r}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_trigger_syst_up   = book<TH1F>("cos_theta2r_Lep_mu_trigger_syst_up",   "cos(#theta_{lep}^{r}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_trigger_syst_down = book<TH1F>("cos_theta2r_Lep_mu_trigger_syst_down",   "cos(#theta_{lep}^{r}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_iso_stat_up       = book<TH1F>("cos_theta2r_Lep_mu_iso_stat_up",   "cos(#theta_{lep}^{r}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_iso_stat_down     = book<TH1F>("cos_theta2r_Lep_mu_iso_stat_down",   "cos(#theta_{lep}^{r}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_iso_syst_up       = book<TH1F>("cos_theta2r_Lep_mu_iso_syst_up",   "cos(#theta_{lep}^{r}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta2r_Lep_mu_iso_syst_down     = book<TH1F>("cos_theta2r_Lep_mu_iso_syst_down",   "cos(#theta_{lep}^{r}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta2r_Lep_pu_up                = book<TH1F>("cos_theta2r_Lep_pu_up",   "cos(#theta_{lep}^{r}) pu_up",                               24, -1.0, 1.0);
  // cos_theta2r_Lep_pu_down              = book<TH1F>("cos_theta2r_Lep_pu_down", "cos(#theta_{lep}^{r}) pu_down",                             24, -1.0, 1.0);
  // cos_theta2r_Lep_prefiring_up         = book<TH1F>("cos_theta2r_Lep_prefiring_up",   "cos(#theta_{lep}^{r}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta2r_Lep_prefiring_down       = book<TH1F>("cos_theta2r_Lep_prefiring_down", "cos(#theta_{lep}^{r}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_cferr1_up       = book<TH1F>("cos_theta2r_Lep_btag_cferr1_up", "cos(#theta_{lep}^{r}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_cferr1_down     = book<TH1F>("cos_theta2r_Lep_btag_cferr1_down", "cos(#theta_{lep}^{r}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_cferr2_up       = book<TH1F>("cos_theta2r_Lep_btag_cferr2_up", "cos(#theta_{lep}^{r}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_cferr2_down     = book<TH1F>("cos_theta2r_Lep_btag_cferr2_down", "cos(#theta_{lep}^{r}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_hf_up           = book<TH1F>("cos_theta2r_Lep_btag_hf_up", "cos(#theta_{lep}^{r}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_hf_down         = book<TH1F>("cos_theta2r_Lep_btag_hf_down", "cos(#theta_{lep}^{r}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_hfstats1_up     = book<TH1F>("cos_theta2r_Lep_btag_hfstats1_up", "cos(#theta_{lep}^{r}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_hfstats1_down   = book<TH1F>("cos_theta2r_Lep_btag_hfstats1_down", "cos(#theta_{lep}^{r}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_hfstats2_up     = book<TH1F>("cos_theta2r_Lep_btag_hfstats2_up", "cos(#theta_{lep}^{r}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_hfstats2_down   = book<TH1F>("cos_theta2r_Lep_btag_hfstats2_down", "cos(#theta_{lep}^{r}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_lf_up           = book<TH1F>("cos_theta2r_Lep_btag_lf_up", "cos(#theta_{lep}^{r}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_lf_down         = book<TH1F>("cos_theta2r_Lep_btag_lf_down", "cos(#theta_{lep}^{r}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_lfstats1_up     = book<TH1F>("cos_theta2r_Lep_btag_lfstats1_up", "cos(#theta_{lep}^{r}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_lfstats1_down   = book<TH1F>("cos_theta2r_Lep_btag_lfstats1_down", "cos(#theta_{lep}^{r}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_lfstats2_up     = book<TH1F>("cos_theta2r_Lep_btag_lfstats2_up", "cos(#theta_{lep}^{r}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta2r_Lep_btag_lfstats2_down   = book<TH1F>("cos_theta2r_Lep_btag_lfstats2_down", "cos(#theta_{lep}^{r}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta2r_Lep_ttag_corr_up         = book<TH1F>("cos_theta2r_Lep_ttag_corr_up", "cos(#theta_{lep}^{r}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta2r_Lep_ttag_corr_down       = book<TH1F>("cos_theta2r_Lep_ttag_corr_down", "cos(#theta_{lep}^{r}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta2r_Lep_ttag_uncorr_up       = book<TH1F>("cos_theta2r_Lep_ttag_uncorr_up", "cos(#theta_{lep}^{r}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta2r_Lep_ttag_uncorr_down     = book<TH1F>("cos_theta2r_Lep_ttag_uncorr_down", "cos(#theta_{lep}^{r}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta2r_Lep_tmistag_up           = book<TH1F>("cos_theta2r_Lep_tmistag_up", "cos(#theta_{lep}^{r}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta2r_Lep_tmistag_down         = book<TH1F>("cos_theta2r_Lep_tmistag_down", "cos(#theta_{lep}^{r}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta2r_Lep_toppt_a_up           = book<TH1F>("cos_theta2r_Lep_toppt_a_up", "cos(#theta_{lep}^{r}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta2r_Lep_toppt_a_down         = book<TH1F>("cos_theta2r_Lep_toppt_a_down", "cos(#theta_{lep}^{r}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta2r_Lep_toppt_b_up           = book<TH1F>("cos_theta2r_Lep_toppt_b_up", "cos(#theta_{lep}^{r}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta2r_Lep_toppt_b_down         = book<TH1F>("cos_theta2r_Lep_toppt_b_down", "cos(#theta_{lep}^{r}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta2r_Lep_murmuf_upup          = book<TH1F>("cos_theta2r_Lep_murmuf_upup", "cos(#theta_{lep}^{r}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta2r_Lep_murmuf_upnone        = book<TH1F>("cos_theta2r_Lep_murmuf_upnone", "cos(#theta_{lep}^{r}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta2r_Lep_murmuf_noneup        = book<TH1F>("cos_theta2r_Lep_murmuf_noneup", "cos(#theta_{lep}^{r}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta2r_Lep_murmuf_nonedown      = book<TH1F>("cos_theta2r_Lep_murmuf_nonedown", "cos(#theta_{lep}^{r}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta2r_Lep_murmuf_downnone      = book<TH1F>("cos_theta2r_Lep_murmuf_downnone", "cos(#theta_{lep}^{r}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta2r_Lep_murmuf_downdown      = book<TH1F>("cos_theta2r_Lep_murmuf_downdown", "cos(#theta_{lep}^{r}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta2r_Lep_isr_up               = book<TH1F>("cos_theta2r_Lep_isr_up", "cos(#theta_{lep}^{r}) isr_up",                               24, -1.0, 1.0);
  // cos_theta2r_Lep_isr_down             = book<TH1F>("cos_theta2r_Lep_isr_down", "cos(#theta_{lep}^{r}) isr_down",                           24, -1.0, 1.0);
  // cos_theta2r_Lep_fsr_up               = book<TH1F>("cos_theta2r_Lep_fsr_up", "cos(#theta_{lep}^{r}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta2r_Lep_fsr_down             = book<TH1F>("cos_theta2r_Lep_fsr_down", "cos(#theta_{lep}^{r}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta2n_Lep = book<TH1F>("cos_theta2n_Lep", "cos(#theta_{lep}^{n})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta2n_Lep_ele_reco_up          = book<TH1F>("cos_theta2n_Lep_ele_reco_up",   "cos(#theta_{lep}^{n}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta2n_Lep_ele_reco_down        = book<TH1F>("cos_theta2n_Lep_ele_reco_down", "cos(#theta_{lep}^{n}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta2n_Lep_ele_id_up            = book<TH1F>("cos_theta2n_Lep_ele_id_up",   "cos(#theta_{lep}^{n}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta2n_Lep_ele_id_down          = book<TH1F>("cos_theta2n_Lep_ele_id_down", "cos(#theta_{lep}^{n}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta2n_Lep_ele_trigger_up       = book<TH1F>("cos_theta2n_Lep_ele_trigger_up",   "cos(#theta_{lep}^{n}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta2n_Lep_ele_trigger_down     = book<TH1F>("cos_theta2n_Lep_ele_trigger_down", "cos(#theta_{lep}^{n}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_reco_up           = book<TH1F>("cos_theta2n_Lep_mu_reco_up",   "cos(#theta_{lep}^{n}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_reco_down         = book<TH1F>("cos_theta2n_Lep_mu_reco_down", "cos(#theta_{lep}^{n}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_id_stat_up        = book<TH1F>("cos_theta2n_Lep_mu_id_stat_up",   "cos(#theta_{lep}^{n}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_id_stat_down      = book<TH1F>("cos_theta2n_Lep_mu_id_stat_down",   "cos(#theta_{lep}^{n}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_id_syst_up        = book<TH1F>("cos_theta2n_Lep_mu_id_syst_up",   "cos(#theta_{lep}^{n}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_id_syst_down      = book<TH1F>("cos_theta2n_Lep_mu_id_syst_down",   "cos(#theta_{lep}^{n}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_trigger_stat_up   = book<TH1F>("cos_theta2n_Lep_mu_trigger_stat_up",   "cos(#theta_{lep}^{n}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_trigger_stat_down = book<TH1F>("cos_theta2n_Lep_mu_trigger_stat_down",   "cos(#theta_{lep}^{n}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_trigger_syst_up   = book<TH1F>("cos_theta2n_Lep_mu_trigger_syst_up",   "cos(#theta_{lep}^{n}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_trigger_syst_down = book<TH1F>("cos_theta2n_Lep_mu_trigger_syst_down",   "cos(#theta_{lep}^{n}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_iso_stat_up       = book<TH1F>("cos_theta2n_Lep_mu_iso_stat_up",   "cos(#theta_{lep}^{n}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_iso_stat_down     = book<TH1F>("cos_theta2n_Lep_mu_iso_stat_down",   "cos(#theta_{lep}^{n}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_iso_syst_up       = book<TH1F>("cos_theta2n_Lep_mu_iso_syst_up",   "cos(#theta_{lep}^{n}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta2n_Lep_mu_iso_syst_down     = book<TH1F>("cos_theta2n_Lep_mu_iso_syst_down",   "cos(#theta_{lep}^{n}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta2n_Lep_pu_up                = book<TH1F>("cos_theta2n_Lep_pu_up",   "cos(#theta_{lep}^{n}) pu_up",                               24, -1.0, 1.0);
  // cos_theta2n_Lep_pu_down              = book<TH1F>("cos_theta2n_Lep_pu_down", "cos(#theta_{lep}^{n}) pu_down",                             24, -1.0, 1.0);
  // cos_theta2n_Lep_prefiring_up         = book<TH1F>("cos_theta2n_Lep_prefiring_up",   "cos(#theta_{lep}^{n}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta2n_Lep_prefiring_down       = book<TH1F>("cos_theta2n_Lep_prefiring_down", "cos(#theta_{lep}^{n}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_cferr1_up       = book<TH1F>("cos_theta2n_Lep_btag_cferr1_up", "cos(#theta_{lep}^{n}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_cferr1_down     = book<TH1F>("cos_theta2n_Lep_btag_cferr1_down", "cos(#theta_{lep}^{n}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_cferr2_up       = book<TH1F>("cos_theta2n_Lep_btag_cferr2_up", "cos(#theta_{lep}^{n}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_cferr2_down     = book<TH1F>("cos_theta2n_Lep_btag_cferr2_down", "cos(#theta_{lep}^{n}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_hf_up           = book<TH1F>("cos_theta2n_Lep_btag_hf_up", "cos(#theta_{lep}^{n}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_hf_down         = book<TH1F>("cos_theta2n_Lep_btag_hf_down", "cos(#theta_{lep}^{n}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_hfstats1_up     = book<TH1F>("cos_theta2n_Lep_btag_hfstats1_up", "cos(#theta_{lep}^{n}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_hfstats1_down   = book<TH1F>("cos_theta2n_Lep_btag_hfstats1_down", "cos(#theta_{lep}^{n}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_hfstats2_up     = book<TH1F>("cos_theta2n_Lep_btag_hfstats2_up", "cos(#theta_{lep}^{n}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_hfstats2_down   = book<TH1F>("cos_theta2n_Lep_btag_hfstats2_down", "cos(#theta_{lep}^{n}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_lf_up           = book<TH1F>("cos_theta2n_Lep_btag_lf_up", "cos(#theta_{lep}^{n}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_lf_down         = book<TH1F>("cos_theta2n_Lep_btag_lf_down", "cos(#theta_{lep}^{n}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_lfstats1_up     = book<TH1F>("cos_theta2n_Lep_btag_lfstats1_up", "cos(#theta_{lep}^{n}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_lfstats1_down   = book<TH1F>("cos_theta2n_Lep_btag_lfstats1_down", "cos(#theta_{lep}^{n}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_lfstats2_up     = book<TH1F>("cos_theta2n_Lep_btag_lfstats2_up", "cos(#theta_{lep}^{n}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta2n_Lep_btag_lfstats2_down   = book<TH1F>("cos_theta2n_Lep_btag_lfstats2_down", "cos(#theta_{lep}^{n}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta2n_Lep_ttag_corr_up         = book<TH1F>("cos_theta2n_Lep_ttag_corr_up", "cos(#theta_{lep}^{n}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta2n_Lep_ttag_corr_down       = book<TH1F>("cos_theta2n_Lep_ttag_corr_down", "cos(#theta_{lep}^{n}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta2n_Lep_ttag_uncorr_up       = book<TH1F>("cos_theta2n_Lep_ttag_uncorr_up", "cos(#theta_{lep}^{n}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta2n_Lep_ttag_uncorr_down     = book<TH1F>("cos_theta2n_Lep_ttag_uncorr_down", "cos(#theta_{lep}^{n}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta2n_Lep_tmistag_up           = book<TH1F>("cos_theta2n_Lep_tmistag_up", "cos(#theta_{lep}^{n}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta2n_Lep_tmistag_down         = book<TH1F>("cos_theta2n_Lep_tmistag_down", "cos(#theta_{lep}^{n}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta2n_Lep_toppt_a_up           = book<TH1F>("cos_theta2n_Lep_toppt_a_up", "cos(#theta_{lep}^{n}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta2n_Lep_toppt_a_down         = book<TH1F>("cos_theta2n_Lep_toppt_a_down", "cos(#theta_{lep}^{n}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta2n_Lep_toppt_b_up           = book<TH1F>("cos_theta2n_Lep_toppt_b_up", "cos(#theta_{lep}^{n}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta2n_Lep_toppt_b_down         = book<TH1F>("cos_theta2n_Lep_toppt_b_down", "cos(#theta_{lep}^{n}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta2n_Lep_murmuf_upup          = book<TH1F>("cos_theta2n_Lep_murmuf_upup", "cos(#theta_{lep}^{n}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta2n_Lep_murmuf_upnone        = book<TH1F>("cos_theta2n_Lep_murmuf_upnone", "cos(#theta_{lep}^{n}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta2n_Lep_murmuf_noneup        = book<TH1F>("cos_theta2n_Lep_murmuf_noneup", "cos(#theta_{lep}^{n}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta2n_Lep_murmuf_nonedown      = book<TH1F>("cos_theta2n_Lep_murmuf_nonedown", "cos(#theta_{lep}^{n}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta2n_Lep_murmuf_downnone      = book<TH1F>("cos_theta2n_Lep_murmuf_downnone", "cos(#theta_{lep}^{n}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta2n_Lep_murmuf_downdown      = book<TH1F>("cos_theta2n_Lep_murmuf_downdown", "cos(#theta_{lep}^{n}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta2n_Lep_isr_up               = book<TH1F>("cos_theta2n_Lep_isr_up", "cos(#theta_{lep}^{n}) isr_up",                               24, -1.0, 1.0);
  // cos_theta2n_Lep_isr_down             = book<TH1F>("cos_theta2n_Lep_isr_down", "cos(#theta_{lep}^{n}) isr_down",                           24, -1.0, 1.0);
  // cos_theta2n_Lep_fsr_up               = book<TH1F>("cos_theta2n_Lep_fsr_up", "cos(#theta_{lep}^{n}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta2n_Lep_fsr_down             = book<TH1F>("cos_theta2n_Lep_fsr_down", "cos(#theta_{lep}^{n}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta2kStar_Lep = book<TH1F>("cos_theta2kStar_Lep", "cos(#theta_{lep}^{k*})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta2kStar_Lep_ele_reco_up          = book<TH1F>("cos_theta2kStar_Lep_ele_reco_up",   "cos(#theta_{lep}^{k*}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta2kStar_Lep_ele_reco_down        = book<TH1F>("cos_theta2kStar_Lep_ele_reco_down", "cos(#theta_{lep}^{k*}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta2kStar_Lep_ele_id_up            = book<TH1F>("cos_theta2kStar_Lep_ele_id_up",   "cos(#theta_{lep}^{k*}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta2kStar_Lep_ele_id_down          = book<TH1F>("cos_theta2kStar_Lep_ele_id_down", "cos(#theta_{lep}^{k*}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta2kStar_Lep_ele_trigger_up       = book<TH1F>("cos_theta2kStar_Lep_ele_trigger_up",   "cos(#theta_{lep}^{k*}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_ele_trigger_down     = book<TH1F>("cos_theta2kStar_Lep_ele_trigger_down", "cos(#theta_{lep}^{k*}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_reco_up           = book<TH1F>("cos_theta2kStar_Lep_mu_reco_up",   "cos(#theta_{lep}^{k*}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_reco_down         = book<TH1F>("cos_theta2kStar_Lep_mu_reco_down", "cos(#theta_{lep}^{k*}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_id_stat_up        = book<TH1F>("cos_theta2kStar_Lep_mu_id_stat_up",   "cos(#theta_{lep}^{k*}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_id_stat_down      = book<TH1F>("cos_theta2kStar_Lep_mu_id_stat_down",   "cos(#theta_{lep}^{k*}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_id_syst_up        = book<TH1F>("cos_theta2kStar_Lep_mu_id_syst_up",   "cos(#theta_{lep}^{k*}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_id_syst_down      = book<TH1F>("cos_theta2kStar_Lep_mu_id_syst_down",   "cos(#theta_{lep}^{k*}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_trigger_stat_up   = book<TH1F>("cos_theta2kStar_Lep_mu_trigger_stat_up",   "cos(#theta_{lep}^{k*}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_trigger_stat_down = book<TH1F>("cos_theta2kStar_Lep_mu_trigger_stat_down",   "cos(#theta_{lep}^{k*}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_trigger_syst_up   = book<TH1F>("cos_theta2kStar_Lep_mu_trigger_syst_up",   "cos(#theta_{lep}^{k*}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_trigger_syst_down = book<TH1F>("cos_theta2kStar_Lep_mu_trigger_syst_down",   "cos(#theta_{lep}^{k*}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_iso_stat_up       = book<TH1F>("cos_theta2kStar_Lep_mu_iso_stat_up",   "cos(#theta_{lep}^{k*}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_iso_stat_down     = book<TH1F>("cos_theta2kStar_Lep_mu_iso_stat_down",   "cos(#theta_{lep}^{k*}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_iso_syst_up       = book<TH1F>("cos_theta2kStar_Lep_mu_iso_syst_up",   "cos(#theta_{lep}^{k*}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_mu_iso_syst_down     = book<TH1F>("cos_theta2kStar_Lep_mu_iso_syst_down",   "cos(#theta_{lep}^{k*}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta2kStar_Lep_pu_up                = book<TH1F>("cos_theta2kStar_Lep_pu_up",   "cos(#theta_{lep}^{k*}) pu_up",                               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_pu_down              = book<TH1F>("cos_theta2kStar_Lep_pu_down", "cos(#theta_{lep}^{k*}) pu_down",                             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_prefiring_up         = book<TH1F>("cos_theta2kStar_Lep_prefiring_up",   "cos(#theta_{lep}^{k*}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta2kStar_Lep_prefiring_down       = book<TH1F>("cos_theta2kStar_Lep_prefiring_down", "cos(#theta_{lep}^{k*}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_cferr1_up       = book<TH1F>("cos_theta2kStar_Lep_btag_cferr1_up", "cos(#theta_{lep}^{k*}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_cferr1_down     = book<TH1F>("cos_theta2kStar_Lep_btag_cferr1_down", "cos(#theta_{lep}^{k*}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_cferr2_up       = book<TH1F>("cos_theta2kStar_Lep_btag_cferr2_up", "cos(#theta_{lep}^{k*}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_cferr2_down     = book<TH1F>("cos_theta2kStar_Lep_btag_cferr2_down", "cos(#theta_{lep}^{k*}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_hf_up           = book<TH1F>("cos_theta2kStar_Lep_btag_hf_up", "cos(#theta_{lep}^{k*}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_hf_down         = book<TH1F>("cos_theta2kStar_Lep_btag_hf_down", "cos(#theta_{lep}^{k*}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_hfstats1_up     = book<TH1F>("cos_theta2kStar_Lep_btag_hfstats1_up", "cos(#theta_{lep}^{k*}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_hfstats1_down   = book<TH1F>("cos_theta2kStar_Lep_btag_hfstats1_down", "cos(#theta_{lep}^{k*}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_hfstats2_up     = book<TH1F>("cos_theta2kStar_Lep_btag_hfstats2_up", "cos(#theta_{lep}^{k*}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_hfstats2_down   = book<TH1F>("cos_theta2kStar_Lep_btag_hfstats2_down", "cos(#theta_{lep}^{k*}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_lf_up           = book<TH1F>("cos_theta2kStar_Lep_btag_lf_up", "cos(#theta_{lep}^{k*}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_lf_down         = book<TH1F>("cos_theta2kStar_Lep_btag_lf_down", "cos(#theta_{lep}^{k*}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_lfstats1_up     = book<TH1F>("cos_theta2kStar_Lep_btag_lfstats1_up", "cos(#theta_{lep}^{k*}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_lfstats1_down   = book<TH1F>("cos_theta2kStar_Lep_btag_lfstats1_down", "cos(#theta_{lep}^{k*}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_lfstats2_up     = book<TH1F>("cos_theta2kStar_Lep_btag_lfstats2_up", "cos(#theta_{lep}^{k*}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_btag_lfstats2_down   = book<TH1F>("cos_theta2kStar_Lep_btag_lfstats2_down", "cos(#theta_{lep}^{k*}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta2kStar_Lep_ttag_corr_up         = book<TH1F>("cos_theta2kStar_Lep_ttag_corr_up", "cos(#theta_{lep}^{k*}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta2kStar_Lep_ttag_corr_down       = book<TH1F>("cos_theta2kStar_Lep_ttag_corr_down", "cos(#theta_{lep}^{k*}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_ttag_uncorr_up       = book<TH1F>("cos_theta2kStar_Lep_ttag_uncorr_up", "cos(#theta_{lep}^{k*}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_ttag_uncorr_down     = book<TH1F>("cos_theta2kStar_Lep_ttag_uncorr_down", "cos(#theta_{lep}^{k*}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_tmistag_up           = book<TH1F>("cos_theta2kStar_Lep_tmistag_up", "cos(#theta_{lep}^{k*}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta2kStar_Lep_tmistag_down         = book<TH1F>("cos_theta2kStar_Lep_tmistag_down", "cos(#theta_{lep}^{k*}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_toppt_a_up           = book<TH1F>("cos_theta2kStar_Lep_toppt_a_up", "cos(#theta_{lep}^{k*}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta2kStar_Lep_toppt_a_down         = book<TH1F>("cos_theta2kStar_Lep_toppt_a_down", "cos(#theta_{lep}^{k*}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_toppt_b_up           = book<TH1F>("cos_theta2kStar_Lep_toppt_b_up", "cos(#theta_{lep}^{k*}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta2kStar_Lep_toppt_b_down         = book<TH1F>("cos_theta2kStar_Lep_toppt_b_down", "cos(#theta_{lep}^{k*}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_murmuf_upup          = book<TH1F>("cos_theta2kStar_Lep_murmuf_upup", "cos(#theta_{lep}^{k*}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta2kStar_Lep_murmuf_upnone        = book<TH1F>("cos_theta2kStar_Lep_murmuf_upnone", "cos(#theta_{lep}^{k*}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta2kStar_Lep_murmuf_noneup        = book<TH1F>("cos_theta2kStar_Lep_murmuf_noneup", "cos(#theta_{lep}^{k*}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta2kStar_Lep_murmuf_nonedown      = book<TH1F>("cos_theta2kStar_Lep_murmuf_nonedown", "cos(#theta_{lep}^{k*}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_murmuf_downnone      = book<TH1F>("cos_theta2kStar_Lep_murmuf_downnone", "cos(#theta_{lep}^{k*}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_murmuf_downdown      = book<TH1F>("cos_theta2kStar_Lep_murmuf_downdown", "cos(#theta_{lep}^{k*}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta2kStar_Lep_isr_up               = book<TH1F>("cos_theta2kStar_Lep_isr_up", "cos(#theta_{lep}^{k*}) isr_up",                               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_isr_down             = book<TH1F>("cos_theta2kStar_Lep_isr_down", "cos(#theta_{lep}^{k*}) isr_down",                           24, -1.0, 1.0);
  // cos_theta2kStar_Lep_fsr_up               = book<TH1F>("cos_theta2kStar_Lep_fsr_up", "cos(#theta_{lep}^{k*}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta2kStar_Lep_fsr_down             = book<TH1F>("cos_theta2kStar_Lep_fsr_down", "cos(#theta_{lep}^{k*}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta2rStar_Lep = book<TH1F>("cos_theta2rStar_Lep", "cos(#theta_{lep}^{r*})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta2rStar_Lep_ele_reco_up          = book<TH1F>("cos_theta2rStar_Lep_ele_reco_up",   "cos(#theta_{lep}^{r*}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta2rStar_Lep_ele_reco_down        = book<TH1F>("cos_theta2rStar_Lep_ele_reco_down", "cos(#theta_{lep}^{r*}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta2rStar_Lep_ele_id_up            = book<TH1F>("cos_theta2rStar_Lep_ele_id_up",   "cos(#theta_{lep}^{r*}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta2rStar_Lep_ele_id_down          = book<TH1F>("cos_theta2rStar_Lep_ele_id_down", "cos(#theta_{lep}^{r*}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta2rStar_Lep_ele_trigger_up       = book<TH1F>("cos_theta2rStar_Lep_ele_trigger_up",   "cos(#theta_{lep}^{r*}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_ele_trigger_down     = book<TH1F>("cos_theta2rStar_Lep_ele_trigger_down", "cos(#theta_{lep}^{r*}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_reco_up           = book<TH1F>("cos_theta2rStar_Lep_mu_reco_up",   "cos(#theta_{lep}^{r*}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_reco_down         = book<TH1F>("cos_theta2rStar_Lep_mu_reco_down", "cos(#theta_{lep}^{r*}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_id_stat_up        = book<TH1F>("cos_theta2rStar_Lep_mu_id_stat_up",   "cos(#theta_{lep}^{r*}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_id_stat_down      = book<TH1F>("cos_theta2rStar_Lep_mu_id_stat_down",   "cos(#theta_{lep}^{r*}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_id_syst_up        = book<TH1F>("cos_theta2rStar_Lep_mu_id_syst_up",   "cos(#theta_{lep}^{r*}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_id_syst_down      = book<TH1F>("cos_theta2rStar_Lep_mu_id_syst_down",   "cos(#theta_{lep}^{r*}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_trigger_stat_up   = book<TH1F>("cos_theta2rStar_Lep_mu_trigger_stat_up",   "cos(#theta_{lep}^{r*}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_trigger_stat_down = book<TH1F>("cos_theta2rStar_Lep_mu_trigger_stat_down",   "cos(#theta_{lep}^{r*}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_trigger_syst_up   = book<TH1F>("cos_theta2rStar_Lep_mu_trigger_syst_up",   "cos(#theta_{lep}^{r*}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_trigger_syst_down = book<TH1F>("cos_theta2rStar_Lep_mu_trigger_syst_down",   "cos(#theta_{lep}^{r*}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_iso_stat_up       = book<TH1F>("cos_theta2rStar_Lep_mu_iso_stat_up",   "cos(#theta_{lep}^{r*}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_iso_stat_down     = book<TH1F>("cos_theta2rStar_Lep_mu_iso_stat_down",   "cos(#theta_{lep}^{r*}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_iso_syst_up       = book<TH1F>("cos_theta2rStar_Lep_mu_iso_syst_up",   "cos(#theta_{lep}^{r*}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_mu_iso_syst_down     = book<TH1F>("cos_theta2rStar_Lep_mu_iso_syst_down",   "cos(#theta_{lep}^{r*}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta2rStar_Lep_pu_up                = book<TH1F>("cos_theta2rStar_Lep_pu_up",   "cos(#theta_{lep}^{r*}) pu_up",                               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_pu_down              = book<TH1F>("cos_theta2rStar_Lep_pu_down", "cos(#theta_{lep}^{r*}) pu_down",                             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_prefiring_up         = book<TH1F>("cos_theta2rStar_Lep_prefiring_up",   "cos(#theta_{lep}^{r*}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta2rStar_Lep_prefiring_down       = book<TH1F>("cos_theta2rStar_Lep_prefiring_down", "cos(#theta_{lep}^{r*}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_cferr1_up       = book<TH1F>("cos_theta2rStar_Lep_btag_cferr1_up", "cos(#theta_{lep}^{r*}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_cferr1_down     = book<TH1F>("cos_theta2rStar_Lep_btag_cferr1_down", "cos(#theta_{lep}^{r*}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_cferr2_up       = book<TH1F>("cos_theta2rStar_Lep_btag_cferr2_up", "cos(#theta_{lep}^{r*}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_cferr2_down     = book<TH1F>("cos_theta2rStar_Lep_btag_cferr2_down", "cos(#theta_{lep}^{r*}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_hf_up           = book<TH1F>("cos_theta2rStar_Lep_btag_hf_up", "cos(#theta_{lep}^{r*}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_hf_down         = book<TH1F>("cos_theta2rStar_Lep_btag_hf_down", "cos(#theta_{lep}^{r*}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_hfstats1_up     = book<TH1F>("cos_theta2rStar_Lep_btag_hfstats1_up", "cos(#theta_{lep}^{r*}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_hfstats1_down   = book<TH1F>("cos_theta2rStar_Lep_btag_hfstats1_down", "cos(#theta_{lep}^{r*}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_hfstats2_up     = book<TH1F>("cos_theta2rStar_Lep_btag_hfstats2_up", "cos(#theta_{lep}^{r*}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_hfstats2_down   = book<TH1F>("cos_theta2rStar_Lep_btag_hfstats2_down", "cos(#theta_{lep}^{r*}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_lf_up           = book<TH1F>("cos_theta2rStar_Lep_btag_lf_up", "cos(#theta_{lep}^{r*}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_lf_down         = book<TH1F>("cos_theta2rStar_Lep_btag_lf_down", "cos(#theta_{lep}^{r*}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_lfstats1_up     = book<TH1F>("cos_theta2rStar_Lep_btag_lfstats1_up", "cos(#theta_{lep}^{r*}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_lfstats1_down   = book<TH1F>("cos_theta2rStar_Lep_btag_lfstats1_down", "cos(#theta_{lep}^{r*}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_lfstats2_up     = book<TH1F>("cos_theta2rStar_Lep_btag_lfstats2_up", "cos(#theta_{lep}^{r*}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_btag_lfstats2_down   = book<TH1F>("cos_theta2rStar_Lep_btag_lfstats2_down", "cos(#theta_{lep}^{r*}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta2rStar_Lep_ttag_corr_up         = book<TH1F>("cos_theta2rStar_Lep_ttag_corr_up", "cos(#theta_{lep}^{r*}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta2rStar_Lep_ttag_corr_down       = book<TH1F>("cos_theta2rStar_Lep_ttag_corr_down", "cos(#theta_{lep}^{r*}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_ttag_uncorr_up       = book<TH1F>("cos_theta2rStar_Lep_ttag_uncorr_up", "cos(#theta_{lep}^{r*}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_ttag_uncorr_down     = book<TH1F>("cos_theta2rStar_Lep_ttag_uncorr_down", "cos(#theta_{lep}^{r*}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_tmistag_up           = book<TH1F>("cos_theta2rStar_Lep_tmistag_up", "cos(#theta_{lep}^{r*}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta2rStar_Lep_tmistag_down         = book<TH1F>("cos_theta2rStar_Lep_tmistag_down", "cos(#theta_{lep}^{r*}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_toppt_a_up           = book<TH1F>("cos_theta2rStar_Lep_toppt_a_up", "cos(#theta_{lep}^{r*}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta2rStar_Lep_toppt_a_down         = book<TH1F>("cos_theta2rStar_Lep_toppt_a_down", "cos(#theta_{lep}^{r*}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_toppt_b_up           = book<TH1F>("cos_theta2rStar_Lep_toppt_b_up", "cos(#theta_{lep}^{r*}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta2rStar_Lep_toppt_b_down         = book<TH1F>("cos_theta2rStar_Lep_toppt_b_down", "cos(#theta_{lep}^{r*}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_murmuf_upup          = book<TH1F>("cos_theta2rStar_Lep_murmuf_upup", "cos(#theta_{lep}^{r*}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta2rStar_Lep_murmuf_upnone        = book<TH1F>("cos_theta2rStar_Lep_murmuf_upnone", "cos(#theta_{lep}^{r*}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta2rStar_Lep_murmuf_noneup        = book<TH1F>("cos_theta2rStar_Lep_murmuf_noneup", "cos(#theta_{lep}^{r*}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta2rStar_Lep_murmuf_nonedown      = book<TH1F>("cos_theta2rStar_Lep_murmuf_nonedown", "cos(#theta_{lep}^{r*}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_murmuf_downnone      = book<TH1F>("cos_theta2rStar_Lep_murmuf_downnone", "cos(#theta_{lep}^{r*}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_murmuf_downdown      = book<TH1F>("cos_theta2rStar_Lep_murmuf_downdown", "cos(#theta_{lep}^{r*}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta2rStar_Lep_isr_up               = book<TH1F>("cos_theta2rStar_Lep_isr_up", "cos(#theta_{lep}^{r*}) isr_up",                               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_isr_down             = book<TH1F>("cos_theta2rStar_Lep_isr_down", "cos(#theta_{lep}^{r*}) isr_down",                           24, -1.0, 1.0);
  // cos_theta2rStar_Lep_fsr_up               = book<TH1F>("cos_theta2rStar_Lep_fsr_up", "cos(#theta_{lep}^{r*}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta2rStar_Lep_fsr_down             = book<TH1F>("cos_theta2rStar_Lep_fsr_down", "cos(#theta_{lep}^{r*}) fsr_down",                           24, -1.0, 1.0);

  // cos_theta1k = book<TH1F>("cos_theta1k", "cos(#theta_{1}^{k})",24, -1, 1);-----------------------------------------------------------------//
  cos_theta1k_ele_reco_up          = book<TH1F>("cos_theta1k_ele_reco_up",   "cos(#theta_{1}^{k}) ele_reco_up",                   24, -1.0, 1.0);
  cos_theta1k_ele_reco_down        = book<TH1F>("cos_theta1k_ele_reco_down", "cos(#theta_{1}^{k}) ele_reco_down",                 24, -1.0, 1.0);
  cos_theta1k_ele_id_up            = book<TH1F>("cos_theta1k_ele_id_up",   "cos(#theta_{1}^{k}) ele_id_up",                       24, -1.0, 1.0);
  cos_theta1k_ele_id_down          = book<TH1F>("cos_theta1k_ele_id_down", "cos(#theta_{1}^{k}) ele_id_down",                     24, -1.0, 1.0);
  cos_theta1k_ele_trigger_up       = book<TH1F>("cos_theta1k_ele_trigger_up",   "cos(#theta_{1}^{k}) ele_trigger_up",             24, -1.0, 1.0);
  cos_theta1k_ele_trigger_down     = book<TH1F>("cos_theta1k_ele_trigger_down", "cos(#theta_{1}^{k}) ele_trigger_down",           24, -1.0, 1.0);
  cos_theta1k_mu_reco_up           = book<TH1F>("cos_theta1k_mu_reco_up",   "cos(#theta_{1}^{k}) mu_reco_up",                     24, -1.0, 1.0);
  cos_theta1k_mu_reco_down         = book<TH1F>("cos_theta1k_mu_reco_down", "cos(#theta_{1}^{k}) mu_reco_down",                   24, -1.0, 1.0);
  cos_theta1k_mu_id_stat_up        = book<TH1F>("cos_theta1k_mu_id_stat_up",   "cos(#theta_{1}^{k}) mu_id_stat_up",               24, -1.0, 1.0);
  cos_theta1k_mu_id_stat_down      = book<TH1F>("cos_theta1k_mu_id_stat_down",   "cos(#theta_{1}^{k}) mu_id_stat_down",           24, -1.0, 1.0);
  cos_theta1k_mu_id_syst_up        = book<TH1F>("cos_theta1k_mu_id_syst_up",   "cos(#theta_{1}^{k}) mu_id_syst_up",               24, -1.0, 1.0);
  cos_theta1k_mu_id_syst_down      = book<TH1F>("cos_theta1k_mu_id_syst_down",   "cos(#theta_{1}^{k}) mu_id_syst_down",           24, -1.0, 1.0);
  cos_theta1k_mu_trigger_stat_up   = book<TH1F>("cos_theta1k_mu_trigger_stat_up",   "cos(#theta_{1}^{k}) mu_trigger_stat_up",     24, -1.0, 1.0);
  cos_theta1k_mu_trigger_stat_down = book<TH1F>("cos_theta1k_mu_trigger_stat_down",   "cos(#theta_{1}^{k}) mu_trigger_stat_down", 24, -1.0, 1.0);
  cos_theta1k_mu_trigger_syst_up   = book<TH1F>("cos_theta1k_mu_trigger_syst_up",   "cos(#theta_{1}^{k}) mu_trigger_syst_up",     24, -1.0, 1.0);
  cos_theta1k_mu_trigger_syst_down = book<TH1F>("cos_theta1k_mu_trigger_syst_down",   "cos(#theta_{1}^{k}) mu_trigger_syst_down", 24, -1.0, 1.0);
  cos_theta1k_mu_iso_stat_up       = book<TH1F>("cos_theta1k_mu_iso_stat_up",   "cos(#theta_{1}^{k}) mu_iso_stat_up",             24, -1.0, 1.0);
  cos_theta1k_mu_iso_stat_down     = book<TH1F>("cos_theta1k_mu_iso_stat_down",   "cos(#theta_{1}^{k}) mu_iso_stat_down",         24, -1.0, 1.0);
  cos_theta1k_mu_iso_syst_up       = book<TH1F>("cos_theta1k_mu_iso_syst_up",   "cos(#theta_{1}^{k}) mu_iso_syst_up",             24, -1.0, 1.0);
  cos_theta1k_mu_iso_syst_down     = book<TH1F>("cos_theta1k_mu_iso_syst_down",   "cos(#theta_{1}^{k}) mu_iso_syst_down",         24, -1.0, 1.0);
  cos_theta1k_pu_up                = book<TH1F>("cos_theta1k_pu_up",   "cos(#theta_{1}^{k}) pu_up",                               24, -1.0, 1.0);
  cos_theta1k_pu_down              = book<TH1F>("cos_theta1k_pu_down", "cos(#theta_{1}^{k}) pu_down",                             24, -1.0, 1.0);
  cos_theta1k_prefiring_up         = book<TH1F>("cos_theta1k_prefiring_up",   "cos(#theta_{1}^{k}) prefiring_up",                 24, -1.0, 1.0);
  cos_theta1k_prefiring_down       = book<TH1F>("cos_theta1k_prefiring_down", "cos(#theta_{1}^{k}) prefiring_down",               24, -1.0, 1.0);
  cos_theta1k_btag_cferr1_up       = book<TH1F>("cos_theta1k_btag_cferr1_up", "cos(#theta_{1}^{k}) btag_cferr1_up",               24, -1.0, 1.0);
  cos_theta1k_btag_cferr1_down     = book<TH1F>("cos_theta1k_btag_cferr1_down", "cos(#theta_{1}^{k}) btag_cferr1_down",           24, -1.0, 1.0);
  cos_theta1k_btag_cferr2_up       = book<TH1F>("cos_theta1k_btag_cferr2_up", "cos(#theta_{1}^{k}) btag_cferr2_up",               24, -1.0, 1.0);
  cos_theta1k_btag_cferr2_down     = book<TH1F>("cos_theta1k_btag_cferr2_down", "cos(#theta_{1}^{k}) btag_cferr2_down",           24, -1.0, 1.0);
  cos_theta1k_btag_hf_up           = book<TH1F>("cos_theta1k_btag_hf_up", "cos(#theta_{1}^{k}) btag_hf_up",                       24, -1.0, 1.0);
  cos_theta1k_btag_hf_down         = book<TH1F>("cos_theta1k_btag_hf_down", "cos(#theta_{1}^{k}) btag_hf_down",                   24, -1.0, 1.0);
  cos_theta1k_btag_hfstats1_up     = book<TH1F>("cos_theta1k_btag_hfstats1_up", "cos(#theta_{1}^{k}) btag_hfstats1_up",           24, -1.0, 1.0);
  cos_theta1k_btag_hfstats1_down   = book<TH1F>("cos_theta1k_btag_hfstats1_down", "cos(#theta_{1}^{k}) btag_hfstats1_down",       24, -1.0, 1.0);
  cos_theta1k_btag_hfstats2_up     = book<TH1F>("cos_theta1k_btag_hfstats2_up", "cos(#theta_{1}^{k}) btag_hfstats2_up",           24, -1.0, 1.0);
  cos_theta1k_btag_hfstats2_down   = book<TH1F>("cos_theta1k_btag_hfstats2_down", "cos(#theta_{1}^{k}) btag_hfstats2_down",       24, -1.0, 1.0);
  cos_theta1k_btag_lf_up           = book<TH1F>("cos_theta1k_btag_lf_up", "cos(#theta_{1}^{k}) btag_lf_up",                       24, -1.0, 1.0);
  cos_theta1k_btag_lf_down         = book<TH1F>("cos_theta1k_btag_lf_down", "cos(#theta_{1}^{k}) btag_lf_down",                   24, -1.0, 1.0);
  cos_theta1k_btag_lfstats1_up     = book<TH1F>("cos_theta1k_btag_lfstats1_up", "cos(#theta_{1}^{k}) btag_lfstats1_up",           24, -1.0, 1.0);
  cos_theta1k_btag_lfstats1_down   = book<TH1F>("cos_theta1k_btag_lfstats1_down", "cos(#theta_{1}^{k}) btag_lfstats1_down",       24, -1.0, 1.0);
  cos_theta1k_btag_lfstats2_up     = book<TH1F>("cos_theta1k_btag_lfstats2_up", "cos(#theta_{1}^{k}) btag_lfstats2_up",           24, -1.0, 1.0);
  cos_theta1k_btag_lfstats2_down   = book<TH1F>("cos_theta1k_btag_lfstats2_down", "cos(#theta_{1}^{k}) btag_lfstats2_down",       24, -1.0, 1.0);
  cos_theta1k_ttag_corr_up         = book<TH1F>("cos_theta1k_ttag_corr_up", "cos(#theta_{1}^{k}) ttag_corr_up",                   24, -1.0, 1.0);
  cos_theta1k_ttag_corr_down       = book<TH1F>("cos_theta1k_ttag_corr_down", "cos(#theta_{1}^{k}) ttag_corr_down",               24, -1.0, 1.0);
  cos_theta1k_ttag_uncorr_up       = book<TH1F>("cos_theta1k_ttag_uncorr_up", "cos(#theta_{1}^{k}) ttag_uncorr_up",               24, -1.0, 1.0);
  cos_theta1k_ttag_uncorr_down     = book<TH1F>("cos_theta1k_ttag_uncorr_down", "cos(#theta_{1}^{k}) ttag_counrr_down",           24, -1.0, 1.0);
  cos_theta1k_tmistag_up           = book<TH1F>("cos_theta1k_tmistag_up", "cos(#theta_{1}^{k}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cos_theta1k_tmistag_down         = book<TH1F>("cos_theta1k_tmistag_down", "cos(#theta_{1}^{k}) [GeV] tmistag_down",             24, -1.0, 1.0);
  cos_theta1k_toppt_a_up           = book<TH1F>("cos_theta1k_toppt_a_up", "cos(#theta_{1}^{k}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cos_theta1k_toppt_a_down         = book<TH1F>("cos_theta1k_toppt_a_down", "cos(#theta_{1}^{k}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cos_theta1k_toppt_b_up           = book<TH1F>("cos_theta1k_toppt_b_up", "cos(#theta_{1}^{k}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cos_theta1k_toppt_b_down         = book<TH1F>("cos_theta1k_toppt_b_down", "cos(#theta_{1}^{k}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cos_theta1k_murmuf_upup          = book<TH1F>("cos_theta1k_murmuf_upup", "cos(#theta_{1}^{k}) murmuf_upup",                     24, -1.0, 1.0);
  cos_theta1k_murmuf_upnone        = book<TH1F>("cos_theta1k_murmuf_upnone", "cos(#theta_{1}^{k}) murmuf_upnone",                 24, -1.0, 1.0);
  cos_theta1k_murmuf_noneup        = book<TH1F>("cos_theta1k_murmuf_noneup", "cos(#theta_{1}^{k}) murmuf_noneup",                 24, -1.0, 1.0);
  cos_theta1k_murmuf_nonedown      = book<TH1F>("cos_theta1k_murmuf_nonedown", "cos(#theta_{1}^{k}) murmuf_nonedown",             24, -1.0, 1.0);
  cos_theta1k_murmuf_downnone      = book<TH1F>("cos_theta1k_murmuf_downnone", "cos(#theta_{1}^{k}) murmuf_downnone",             24, -1.0, 1.0);
  cos_theta1k_murmuf_downdown      = book<TH1F>("cos_theta1k_murmuf_downdown", "cos(#theta_{1}^{k}) murmuf_downdown",             24, -1.0, 1.0);
  cos_theta1k_isr_up               = book<TH1F>("cos_theta1k_isr_up", "cos(#theta_{1}^{k}) isr_up",                               24, -1.0, 1.0);
  cos_theta1k_isr_down             = book<TH1F>("cos_theta1k_isr_down", "cos(#theta_{1}^{k}) isr_down",                           24, -1.0, 1.0);
  cos_theta1k_fsr_up               = book<TH1F>("cos_theta1k_fsr_up", "cos(#theta_{1}^{k}) fsr_up",                               24, -1.0, 1.0);
  cos_theta1k_fsr_down             = book<TH1F>("cos_theta1k_fsr_down", "cos(#theta_{1}^{k}) fsr_down",                           24, -1.0, 1.0);

  // cos_theta1r = book<TH1F>("cos_theta1r", "cos(#theta_{1}^{r})",24, -1, 1);-----------------------------------------------------------------//
  cos_theta1r_ele_reco_up          = book<TH1F>("cos_theta1r_ele_reco_up",   "cos(#theta_{1}^{r}) ele_reco_up",                   24, -1.0, 1.0);
  cos_theta1r_ele_reco_down        = book<TH1F>("cos_theta1r_ele_reco_down", "cos(#theta_{1}^{r}) ele_reco_down",                 24, -1.0, 1.0);
  cos_theta1r_ele_id_up            = book<TH1F>("cos_theta1r_ele_id_up",   "cos(#theta_{1}^{r}) ele_id_up",                       24, -1.0, 1.0);
  cos_theta1r_ele_id_down          = book<TH1F>("cos_theta1r_ele_id_down", "cos(#theta_{1}^{r}) ele_id_down",                     24, -1.0, 1.0);
  cos_theta1r_ele_trigger_up       = book<TH1F>("cos_theta1r_ele_trigger_up",   "cos(#theta_{1}^{r}) ele_trigger_up",             24, -1.0, 1.0);
  cos_theta1r_ele_trigger_down     = book<TH1F>("cos_theta1r_ele_trigger_down", "cos(#theta_{1}^{r}) ele_trigger_down",           24, -1.0, 1.0);
  cos_theta1r_mu_reco_up           = book<TH1F>("cos_theta1r_mu_reco_up",   "cos(#theta_{1}^{r}) mu_reco_up",                     24, -1.0, 1.0);
  cos_theta1r_mu_reco_down         = book<TH1F>("cos_theta1r_mu_reco_down", "cos(#theta_{1}^{r}) mu_reco_down",                   24, -1.0, 1.0);
  cos_theta1r_mu_id_stat_up        = book<TH1F>("cos_theta1r_mu_id_stat_up",   "cos(#theta_{1}^{r}) mu_id_stat_up",               24, -1.0, 1.0);
  cos_theta1r_mu_id_stat_down      = book<TH1F>("cos_theta1r_mu_id_stat_down",   "cos(#theta_{1}^{r}) mu_id_stat_down",           24, -1.0, 1.0);
  cos_theta1r_mu_id_syst_up        = book<TH1F>("cos_theta1r_mu_id_syst_up",   "cos(#theta_{1}^{r}) mu_id_syst_up",               24, -1.0, 1.0);
  cos_theta1r_mu_id_syst_down      = book<TH1F>("cos_theta1r_mu_id_syst_down",   "cos(#theta_{1}^{r}) mu_id_syst_down",           24, -1.0, 1.0);
  cos_theta1r_mu_trigger_stat_up   = book<TH1F>("cos_theta1r_mu_trigger_stat_up",   "cos(#theta_{1}^{r}) mu_trigger_stat_up",     24, -1.0, 1.0);
  cos_theta1r_mu_trigger_stat_down = book<TH1F>("cos_theta1r_mu_trigger_stat_down",   "cos(#theta_{1}^{r}) mu_trigger_stat_down", 24, -1.0, 1.0);
  cos_theta1r_mu_trigger_syst_up   = book<TH1F>("cos_theta1r_mu_trigger_syst_up",   "cos(#theta_{1}^{r}) mu_trigger_syst_up",     24, -1.0, 1.0);
  cos_theta1r_mu_trigger_syst_down = book<TH1F>("cos_theta1r_mu_trigger_syst_down",   "cos(#theta_{1}^{r}) mu_trigger_syst_down", 24, -1.0, 1.0);
  cos_theta1r_mu_iso_stat_up       = book<TH1F>("cos_theta1r_mu_iso_stat_up",   "cos(#theta_{1}^{r}) mu_iso_stat_up",             24, -1.0, 1.0);
  cos_theta1r_mu_iso_stat_down     = book<TH1F>("cos_theta1r_mu_iso_stat_down",   "cos(#theta_{1}^{r}) mu_iso_stat_down",         24, -1.0, 1.0);
  cos_theta1r_mu_iso_syst_up       = book<TH1F>("cos_theta1r_mu_iso_syst_up",   "cos(#theta_{1}^{r}) mu_iso_syst_up",             24, -1.0, 1.0);
  cos_theta1r_mu_iso_syst_down     = book<TH1F>("cos_theta1r_mu_iso_syst_down",   "cos(#theta_{1}^{r}) mu_iso_syst_down",         24, -1.0, 1.0);
  cos_theta1r_pu_up                = book<TH1F>("cos_theta1r_pu_up",   "cos(#theta_{1}^{r}) pu_up",                               24, -1.0, 1.0);
  cos_theta1r_pu_down              = book<TH1F>("cos_theta1r_pu_down", "cos(#theta_{1}^{r}) pu_down",                             24, -1.0, 1.0);
  cos_theta1r_prefiring_up         = book<TH1F>("cos_theta1r_prefiring_up",   "cos(#theta_{1}^{r}) prefiring_up",                 24, -1.0, 1.0);
  cos_theta1r_prefiring_down       = book<TH1F>("cos_theta1r_prefiring_down", "cos(#theta_{1}^{r}) prefiring_down",               24, -1.0, 1.0);
  cos_theta1r_btag_cferr1_up       = book<TH1F>("cos_theta1r_btag_cferr1_up", "cos(#theta_{1}^{r}) btag_cferr1_up",               24, -1.0, 1.0);
  cos_theta1r_btag_cferr1_down     = book<TH1F>("cos_theta1r_btag_cferr1_down", "cos(#theta_{1}^{r}) btag_cferr1_down",           24, -1.0, 1.0);
  cos_theta1r_btag_cferr2_up       = book<TH1F>("cos_theta1r_btag_cferr2_up", "cos(#theta_{1}^{r}) btag_cferr2_up",               24, -1.0, 1.0);
  cos_theta1r_btag_cferr2_down     = book<TH1F>("cos_theta1r_btag_cferr2_down", "cos(#theta_{1}^{r}) btag_cferr2_down",           24, -1.0, 1.0);
  cos_theta1r_btag_hf_up           = book<TH1F>("cos_theta1r_btag_hf_up", "cos(#theta_{1}^{r}) btag_hf_up",                       24, -1.0, 1.0);
  cos_theta1r_btag_hf_down         = book<TH1F>("cos_theta1r_btag_hf_down", "cos(#theta_{1}^{r}) btag_hf_down",                   24, -1.0, 1.0);
  cos_theta1r_btag_hfstats1_up     = book<TH1F>("cos_theta1r_btag_hfstats1_up", "cos(#theta_{1}^{r}) btag_hfstats1_up",           24, -1.0, 1.0);
  cos_theta1r_btag_hfstats1_down   = book<TH1F>("cos_theta1r_btag_hfstats1_down", "cos(#theta_{1}^{r}) btag_hfstats1_down",       24, -1.0, 1.0);
  cos_theta1r_btag_hfstats2_up     = book<TH1F>("cos_theta1r_btag_hfstats2_up", "cos(#theta_{1}^{r}) btag_hfstats2_up",           24, -1.0, 1.0);
  cos_theta1r_btag_hfstats2_down   = book<TH1F>("cos_theta1r_btag_hfstats2_down", "cos(#theta_{1}^{r}) btag_hfstats2_down",       24, -1.0, 1.0);
  cos_theta1r_btag_lf_up           = book<TH1F>("cos_theta1r_btag_lf_up", "cos(#theta_{1}^{r}) btag_lf_up",                       24, -1.0, 1.0);
  cos_theta1r_btag_lf_down         = book<TH1F>("cos_theta1r_btag_lf_down", "cos(#theta_{1}^{r}) btag_lf_down",                   24, -1.0, 1.0);
  cos_theta1r_btag_lfstats1_up     = book<TH1F>("cos_theta1r_btag_lfstats1_up", "cos(#theta_{1}^{r}) btag_lfstats1_up",           24, -1.0, 1.0);
  cos_theta1r_btag_lfstats1_down   = book<TH1F>("cos_theta1r_btag_lfstats1_down", "cos(#theta_{1}^{r}) btag_lfstats1_down",       24, -1.0, 1.0);
  cos_theta1r_btag_lfstats2_up     = book<TH1F>("cos_theta1r_btag_lfstats2_up", "cos(#theta_{1}^{r}) btag_lfstats2_up",           24, -1.0, 1.0);
  cos_theta1r_btag_lfstats2_down   = book<TH1F>("cos_theta1r_btag_lfstats2_down", "cos(#theta_{1}^{r}) btag_lfstats2_down",       24, -1.0, 1.0);
  cos_theta1r_ttag_corr_up         = book<TH1F>("cos_theta1r_ttag_corr_up", "cos(#theta_{1}^{r}) ttag_corr_up",                   24, -1.0, 1.0);
  cos_theta1r_ttag_corr_down       = book<TH1F>("cos_theta1r_ttag_corr_down", "cos(#theta_{1}^{r}) ttag_corr_down",               24, -1.0, 1.0);
  cos_theta1r_ttag_uncorr_up       = book<TH1F>("cos_theta1r_ttag_uncorr_up", "cos(#theta_{1}^{r}) ttag_uncorr_up",               24, -1.0, 1.0);
  cos_theta1r_ttag_uncorr_down     = book<TH1F>("cos_theta1r_ttag_uncorr_down", "cos(#theta_{1}^{r}) ttag_counrr_down",           24, -1.0, 1.0);
  cos_theta1r_tmistag_up           = book<TH1F>("cos_theta1r_tmistag_up", "cos(#theta_{1}^{r}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cos_theta1r_tmistag_down         = book<TH1F>("cos_theta1r_tmistag_down", "cos(#theta_{1}^{r}) [GeV] tmistag_down",             24, -1.0, 1.0);
  cos_theta1r_toppt_a_up           = book<TH1F>("cos_theta1r_toppt_a_up", "cos(#theta_{1}^{r}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cos_theta1r_toppt_a_down         = book<TH1F>("cos_theta1r_toppt_a_down", "cos(#theta_{1}^{r}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cos_theta1r_toppt_b_up           = book<TH1F>("cos_theta1r_toppt_b_up", "cos(#theta_{1}^{r}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cos_theta1r_toppt_b_down         = book<TH1F>("cos_theta1r_toppt_b_down", "cos(#theta_{1}^{r}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cos_theta1r_murmuf_upup          = book<TH1F>("cos_theta1r_murmuf_upup", "cos(#theta_{1}^{r}) murmuf_upup",                     24, -1.0, 1.0);
  cos_theta1r_murmuf_upnone        = book<TH1F>("cos_theta1r_murmuf_upnone", "cos(#theta_{1}^{r}) murmuf_upnone",                 24, -1.0, 1.0);
  cos_theta1r_murmuf_noneup        = book<TH1F>("cos_theta1r_murmuf_noneup", "cos(#theta_{1}^{r}) murmuf_noneup",                 24, -1.0, 1.0);
  cos_theta1r_murmuf_nonedown      = book<TH1F>("cos_theta1r_murmuf_nonedown", "cos(#theta_{1}^{r}) murmuf_nonedown",             24, -1.0, 1.0);
  cos_theta1r_murmuf_downnone      = book<TH1F>("cos_theta1r_murmuf_downnone", "cos(#theta_{1}^{r}) murmuf_downnone",             24, -1.0, 1.0);
  cos_theta1r_murmuf_downdown      = book<TH1F>("cos_theta1r_murmuf_downdown", "cos(#theta_{1}^{r}) murmuf_downdown",             24, -1.0, 1.0);
  cos_theta1r_isr_up               = book<TH1F>("cos_theta1r_isr_up", "cos(#theta_{1}^{r}) isr_up",                               24, -1.0, 1.0);
  cos_theta1r_isr_down             = book<TH1F>("cos_theta1r_isr_down", "cos(#theta_{1}^{r}) isr_down",                           24, -1.0, 1.0);
  cos_theta1r_fsr_up               = book<TH1F>("cos_theta1r_fsr_up", "cos(#theta_{1}^{r}) fsr_up",                               24, -1.0, 1.0);
  cos_theta1r_fsr_down             = book<TH1F>("cos_theta1r_fsr_down", "cos(#theta_{1}^{r}) fsr_down",                           24, -1.0, 1.0);

  // cos_theta1n = book<TH1F>("cos_theta1n", "cos(#theta_{1}^{n})",24, -1, 1);-----------------------------------------------------------------//
  cos_theta1n_ele_reco_up          = book<TH1F>("cos_theta1n_ele_reco_up",   "cos(#theta_{1}^{n}) ele_reco_up",                   24, -1.0, 1.0);
  cos_theta1n_ele_reco_down        = book<TH1F>("cos_theta1n_ele_reco_down", "cos(#theta_{1}^{n}) ele_reco_down",                 24, -1.0, 1.0);
  cos_theta1n_ele_id_up            = book<TH1F>("cos_theta1n_ele_id_up",   "cos(#theta_{1}^{n}) ele_id_up",                       24, -1.0, 1.0);
  cos_theta1n_ele_id_down          = book<TH1F>("cos_theta1n_ele_id_down", "cos(#theta_{1}^{n}) ele_id_down",                     24, -1.0, 1.0);
  cos_theta1n_ele_trigger_up       = book<TH1F>("cos_theta1n_ele_trigger_up",   "cos(#theta_{1}^{n}) ele_trigger_up",             24, -1.0, 1.0);
  cos_theta1n_ele_trigger_down     = book<TH1F>("cos_theta1n_ele_trigger_down", "cos(#theta_{1}^{n}) ele_trigger_down",           24, -1.0, 1.0);
  cos_theta1n_mu_reco_up           = book<TH1F>("cos_theta1n_mu_reco_up",   "cos(#theta_{1}^{n}) mu_reco_up",                     24, -1.0, 1.0);
  cos_theta1n_mu_reco_down         = book<TH1F>("cos_theta1n_mu_reco_down", "cos(#theta_{1}^{n}) mu_reco_down",                   24, -1.0, 1.0);
  cos_theta1n_mu_id_stat_up        = book<TH1F>("cos_theta1n_mu_id_stat_up",   "cos(#theta_{1}^{n}) mu_id_stat_up",               24, -1.0, 1.0);
  cos_theta1n_mu_id_stat_down      = book<TH1F>("cos_theta1n_mu_id_stat_down",   "cos(#theta_{1}^{n}) mu_id_stat_down",           24, -1.0, 1.0);
  cos_theta1n_mu_id_syst_up        = book<TH1F>("cos_theta1n_mu_id_syst_up",   "cos(#theta_{1}^{n}) mu_id_syst_up",               24, -1.0, 1.0);
  cos_theta1n_mu_id_syst_down      = book<TH1F>("cos_theta1n_mu_id_syst_down",   "cos(#theta_{1}^{n}) mu_id_syst_down",           24, -1.0, 1.0);
  cos_theta1n_mu_trigger_stat_up   = book<TH1F>("cos_theta1n_mu_trigger_stat_up",   "cos(#theta_{1}^{n}) mu_trigger_stat_up",     24, -1.0, 1.0);
  cos_theta1n_mu_trigger_stat_down = book<TH1F>("cos_theta1n_mu_trigger_stat_down",   "cos(#theta_{1}^{n}) mu_trigger_stat_down", 24, -1.0, 1.0);
  cos_theta1n_mu_trigger_syst_up   = book<TH1F>("cos_theta1n_mu_trigger_syst_up",   "cos(#theta_{1}^{n}) mu_trigger_syst_up",     24, -1.0, 1.0);
  cos_theta1n_mu_trigger_syst_down = book<TH1F>("cos_theta1n_mu_trigger_syst_down",   "cos(#theta_{1}^{n}) mu_trigger_syst_down", 24, -1.0, 1.0);
  cos_theta1n_mu_iso_stat_up       = book<TH1F>("cos_theta1n_mu_iso_stat_up",   "cos(#theta_{1}^{n}) mu_iso_stat_up",             24, -1.0, 1.0);
  cos_theta1n_mu_iso_stat_down     = book<TH1F>("cos_theta1n_mu_iso_stat_down",   "cos(#theta_{1}^{n}) mu_iso_stat_down",         24, -1.0, 1.0);
  cos_theta1n_mu_iso_syst_up       = book<TH1F>("cos_theta1n_mu_iso_syst_up",   "cos(#theta_{1}^{n}) mu_iso_syst_up",             24, -1.0, 1.0);
  cos_theta1n_mu_iso_syst_down     = book<TH1F>("cos_theta1n_mu_iso_syst_down",   "cos(#theta_{1}^{n}) mu_iso_syst_down",         24, -1.0, 1.0);
  cos_theta1n_pu_up                = book<TH1F>("cos_theta1n_pu_up",   "cos(#theta_{1}^{n}) pu_up",                               24, -1.0, 1.0);
  cos_theta1n_pu_down              = book<TH1F>("cos_theta1n_pu_down", "cos(#theta_{1}^{n}) pu_down",                             24, -1.0, 1.0);
  cos_theta1n_prefiring_up         = book<TH1F>("cos_theta1n_prefiring_up",   "cos(#theta_{1}^{n}) prefiring_up",                 24, -1.0, 1.0);
  cos_theta1n_prefiring_down       = book<TH1F>("cos_theta1n_prefiring_down", "cos(#theta_{1}^{n}) prefiring_down",               24, -1.0, 1.0);
  cos_theta1n_btag_cferr1_up       = book<TH1F>("cos_theta1n_btag_cferr1_up", "cos(#theta_{1}^{n}) btag_cferr1_up",               24, -1.0, 1.0);
  cos_theta1n_btag_cferr1_down     = book<TH1F>("cos_theta1n_btag_cferr1_down", "cos(#theta_{1}^{n}) btag_cferr1_down",           24, -1.0, 1.0);
  cos_theta1n_btag_cferr2_up       = book<TH1F>("cos_theta1n_btag_cferr2_up", "cos(#theta_{1}^{n}) btag_cferr2_up",               24, -1.0, 1.0);
  cos_theta1n_btag_cferr2_down     = book<TH1F>("cos_theta1n_btag_cferr2_down", "cos(#theta_{1}^{n}) btag_cferr2_down",           24, -1.0, 1.0);
  cos_theta1n_btag_hf_up           = book<TH1F>("cos_theta1n_btag_hf_up", "cos(#theta_{1}^{n}) btag_hf_up",                       24, -1.0, 1.0);
  cos_theta1n_btag_hf_down         = book<TH1F>("cos_theta1n_btag_hf_down", "cos(#theta_{1}^{n}) btag_hf_down",                   24, -1.0, 1.0);
  cos_theta1n_btag_hfstats1_up     = book<TH1F>("cos_theta1n_btag_hfstats1_up", "cos(#theta_{1}^{n}) btag_hfstats1_up",           24, -1.0, 1.0);
  cos_theta1n_btag_hfstats1_down   = book<TH1F>("cos_theta1n_btag_hfstats1_down", "cos(#theta_{1}^{n}) btag_hfstats1_down",       24, -1.0, 1.0);
  cos_theta1n_btag_hfstats2_up     = book<TH1F>("cos_theta1n_btag_hfstats2_up", "cos(#theta_{1}^{n}) btag_hfstats2_up",           24, -1.0, 1.0);
  cos_theta1n_btag_hfstats2_down   = book<TH1F>("cos_theta1n_btag_hfstats2_down", "cos(#theta_{1}^{n}) btag_hfstats2_down",       24, -1.0, 1.0);
  cos_theta1n_btag_lf_up           = book<TH1F>("cos_theta1n_btag_lf_up", "cos(#theta_{1}^{n}) btag_lf_up",                       24, -1.0, 1.0);
  cos_theta1n_btag_lf_down         = book<TH1F>("cos_theta1n_btag_lf_down", "cos(#theta_{1}^{n}) btag_lf_down",                   24, -1.0, 1.0);
  cos_theta1n_btag_lfstats1_up     = book<TH1F>("cos_theta1n_btag_lfstats1_up", "cos(#theta_{1}^{n}) btag_lfstats1_up",           24, -1.0, 1.0);
  cos_theta1n_btag_lfstats1_down   = book<TH1F>("cos_theta1n_btag_lfstats1_down", "cos(#theta_{1}^{n}) btag_lfstats1_down",       24, -1.0, 1.0);
  cos_theta1n_btag_lfstats2_up     = book<TH1F>("cos_theta1n_btag_lfstats2_up", "cos(#theta_{1}^{n}) btag_lfstats2_up",           24, -1.0, 1.0);
  cos_theta1n_btag_lfstats2_down   = book<TH1F>("cos_theta1n_btag_lfstats2_down", "cos(#theta_{1}^{n}) btag_lfstats2_down",       24, -1.0, 1.0);
  cos_theta1n_ttag_corr_up         = book<TH1F>("cos_theta1n_ttag_corr_up", "cos(#theta_{1}^{n}) ttag_corr_up",                   24, -1.0, 1.0);
  cos_theta1n_ttag_corr_down       = book<TH1F>("cos_theta1n_ttag_corr_down", "cos(#theta_{1}^{n}) ttag_corr_down",               24, -1.0, 1.0);
  cos_theta1n_ttag_uncorr_up       = book<TH1F>("cos_theta1n_ttag_uncorr_up", "cos(#theta_{1}^{n}) ttag_uncorr_up",               24, -1.0, 1.0);
  cos_theta1n_ttag_uncorr_down     = book<TH1F>("cos_theta1n_ttag_uncorr_down", "cos(#theta_{1}^{n}) ttag_counrr_down",           24, -1.0, 1.0);
  cos_theta1n_tmistag_up           = book<TH1F>("cos_theta1n_tmistag_up", "cos(#theta_{1}^{n}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cos_theta1n_tmistag_down         = book<TH1F>("cos_theta1n_tmistag_down", "cos(#theta_{1}^{n}) [GeV] tmistag_down",             24, -1.0, 1.0);
  cos_theta1n_toppt_a_up           = book<TH1F>("cos_theta1n_toppt_a_up", "cos(#theta_{1}^{n}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cos_theta1n_toppt_a_down         = book<TH1F>("cos_theta1n_toppt_a_down", "cos(#theta_{1}^{n}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cos_theta1n_toppt_b_up           = book<TH1F>("cos_theta1n_toppt_b_up", "cos(#theta_{1}^{n}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cos_theta1n_toppt_b_down         = book<TH1F>("cos_theta1n_toppt_b_down", "cos(#theta_{1}^{n}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cos_theta1n_murmuf_upup          = book<TH1F>("cos_theta1n_murmuf_upup", "cos(#theta_{1}^{n}) murmuf_upup",                     24, -1.0, 1.0);
  cos_theta1n_murmuf_upnone        = book<TH1F>("cos_theta1n_murmuf_upnone", "cos(#theta_{1}^{n}) murmuf_upnone",                 24, -1.0, 1.0);
  cos_theta1n_murmuf_noneup        = book<TH1F>("cos_theta1n_murmuf_noneup", "cos(#theta_{1}^{n}) murmuf_noneup",                 24, -1.0, 1.0);
  cos_theta1n_murmuf_nonedown      = book<TH1F>("cos_theta1n_murmuf_nonedown", "cos(#theta_{1}^{n}) murmuf_nonedown",             24, -1.0, 1.0);
  cos_theta1n_murmuf_downnone      = book<TH1F>("cos_theta1n_murmuf_downnone", "cos(#theta_{1}^{n}) murmuf_downnone",             24, -1.0, 1.0);
  cos_theta1n_murmuf_downdown      = book<TH1F>("cos_theta1n_murmuf_downdown", "cos(#theta_{1}^{n}) murmuf_downdown",             24, -1.0, 1.0);
  cos_theta1n_isr_up               = book<TH1F>("cos_theta1n_isr_up", "cos(#theta_{1}^{n}) isr_up",                               24, -1.0, 1.0);
  cos_theta1n_isr_down             = book<TH1F>("cos_theta1n_isr_down", "cos(#theta_{1}^{n}) isr_down",                           24, -1.0, 1.0);
  cos_theta1n_fsr_up               = book<TH1F>("cos_theta1n_fsr_up", "cos(#theta_{1}^{n}) fsr_up",                               24, -1.0, 1.0);
  cos_theta1n_fsr_down             = book<TH1F>("cos_theta1n_fsr_down", "cos(#theta_{1}^{n}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta1kStar = book<TH1F>("cos_theta1kStar", "cos(#theta_{1}^{k*})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta1kStar_ele_reco_up          = book<TH1F>("cos_theta1kStar_ele_reco_up",   "cos(#theta_{1}^{k*}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta1kStar_ele_reco_down        = book<TH1F>("cos_theta1kStar_ele_reco_down", "cos(#theta_{1}^{k*}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta1kStar_ele_id_up            = book<TH1F>("cos_theta1kStar_ele_id_up",   "cos(#theta_{1}^{k*}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta1kStar_ele_id_down          = book<TH1F>("cos_theta1kStar_ele_id_down", "cos(#theta_{1}^{k*}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta1kStar_ele_trigger_up       = book<TH1F>("cos_theta1kStar_ele_trigger_up",   "cos(#theta_{1}^{k*}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta1kStar_ele_trigger_down     = book<TH1F>("cos_theta1kStar_ele_trigger_down", "cos(#theta_{1}^{k*}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta1kStar_mu_reco_up           = book<TH1F>("cos_theta1kStar_mu_reco_up",   "cos(#theta_{1}^{k*}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta1kStar_mu_reco_down         = book<TH1F>("cos_theta1kStar_mu_reco_down", "cos(#theta_{1}^{k*}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta1kStar_mu_id_stat_up        = book<TH1F>("cos_theta1kStar_mu_id_stat_up",   "cos(#theta_{1}^{k*}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta1kStar_mu_id_stat_down      = book<TH1F>("cos_theta1kStar_mu_id_stat_down",   "cos(#theta_{1}^{k*}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta1kStar_mu_id_syst_up        = book<TH1F>("cos_theta1kStar_mu_id_syst_up",   "cos(#theta_{1}^{k*}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta1kStar_mu_id_syst_down      = book<TH1F>("cos_theta1kStar_mu_id_syst_down",   "cos(#theta_{1}^{k*}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta1kStar_mu_trigger_stat_up   = book<TH1F>("cos_theta1kStar_mu_trigger_stat_up",   "cos(#theta_{1}^{k*}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta1kStar_mu_trigger_stat_down = book<TH1F>("cos_theta1kStar_mu_trigger_stat_down",   "cos(#theta_{1}^{k*}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta1kStar_mu_trigger_syst_up   = book<TH1F>("cos_theta1kStar_mu_trigger_syst_up",   "cos(#theta_{1}^{k*}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta1kStar_mu_trigger_syst_down = book<TH1F>("cos_theta1kStar_mu_trigger_syst_down",   "cos(#theta_{1}^{k*}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta1kStar_mu_iso_stat_up       = book<TH1F>("cos_theta1kStar_mu_iso_stat_up",   "cos(#theta_{1}^{k*}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta1kStar_mu_iso_stat_down     = book<TH1F>("cos_theta1kStar_mu_iso_stat_down",   "cos(#theta_{1}^{k*}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta1kStar_mu_iso_syst_up       = book<TH1F>("cos_theta1kStar_mu_iso_syst_up",   "cos(#theta_{1}^{k*}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta1kStar_mu_iso_syst_down     = book<TH1F>("cos_theta1kStar_mu_iso_syst_down",   "cos(#theta_{1}^{k*}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta1kStar_pu_up                = book<TH1F>("cos_theta1kStar_pu_up",   "cos(#theta_{1}^{k*}) pu_up",                               24, -1.0, 1.0);
  // cos_theta1kStar_pu_down              = book<TH1F>("cos_theta1kStar_pu_down", "cos(#theta_{1}^{k*}) pu_down",                             24, -1.0, 1.0);
  // cos_theta1kStar_prefiring_up         = book<TH1F>("cos_theta1kStar_prefiring_up",   "cos(#theta_{1}^{k*}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta1kStar_prefiring_down       = book<TH1F>("cos_theta1kStar_prefiring_down", "cos(#theta_{1}^{k*}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta1kStar_btag_cferr1_up       = book<TH1F>("cos_theta1kStar_btag_cferr1_up", "cos(#theta_{1}^{k*}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta1kStar_btag_cferr1_down     = book<TH1F>("cos_theta1kStar_btag_cferr1_down", "cos(#theta_{1}^{k*}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta1kStar_btag_cferr2_up       = book<TH1F>("cos_theta1kStar_btag_cferr2_up", "cos(#theta_{1}^{k*}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta1kStar_btag_cferr2_down     = book<TH1F>("cos_theta1kStar_btag_cferr2_down", "cos(#theta_{1}^{k*}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta1kStar_btag_hf_up           = book<TH1F>("cos_theta1kStar_btag_hf_up", "cos(#theta_{1}^{k*}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta1kStar_btag_hf_down         = book<TH1F>("cos_theta1kStar_btag_hf_down", "cos(#theta_{1}^{k*}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta1kStar_btag_hfstats1_up     = book<TH1F>("cos_theta1kStar_btag_hfstats1_up", "cos(#theta_{1}^{k*}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta1kStar_btag_hfstats1_down   = book<TH1F>("cos_theta1kStar_btag_hfstats1_down", "cos(#theta_{1}^{k*}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta1kStar_btag_hfstats2_up     = book<TH1F>("cos_theta1kStar_btag_hfstats2_up", "cos(#theta_{1}^{k*}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta1kStar_btag_hfstats2_down   = book<TH1F>("cos_theta1kStar_btag_hfstats2_down", "cos(#theta_{1}^{k*}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta1kStar_btag_lf_up           = book<TH1F>("cos_theta1kStar_btag_lf_up", "cos(#theta_{1}^{k*}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta1kStar_btag_lf_down         = book<TH1F>("cos_theta1kStar_btag_lf_down", "cos(#theta_{1}^{k*}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta1kStar_btag_lfstats1_up     = book<TH1F>("cos_theta1kStar_btag_lfstats1_up", "cos(#theta_{1}^{k*}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta1kStar_btag_lfstats1_down   = book<TH1F>("cos_theta1kStar_btag_lfstats1_down", "cos(#theta_{1}^{k*}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta1kStar_btag_lfstats2_up     = book<TH1F>("cos_theta1kStar_btag_lfstats2_up", "cos(#theta_{1}^{k*}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta1kStar_btag_lfstats2_down   = book<TH1F>("cos_theta1kStar_btag_lfstats2_down", "cos(#theta_{1}^{k*}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta1kStar_ttag_corr_up         = book<TH1F>("cos_theta1kStar_ttag_corr_up", "cos(#theta_{1}^{k*}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta1kStar_ttag_corr_down       = book<TH1F>("cos_theta1kStar_ttag_corr_down", "cos(#theta_{1}^{k*}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta1kStar_ttag_uncorr_up       = book<TH1F>("cos_theta1kStar_ttag_uncorr_up", "cos(#theta_{1}^{k*}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta1kStar_ttag_uncorr_down     = book<TH1F>("cos_theta1kStar_ttag_uncorr_down", "cos(#theta_{1}^{k*}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta1kStar_tmistag_up           = book<TH1F>("cos_theta1kStar_tmistag_up", "cos(#theta_{1}^{k*}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta1kStar_tmistag_down         = book<TH1F>("cos_theta1kStar_tmistag_down", "cos(#theta_{1}^{k*}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta1kStar_toppt_a_up           = book<TH1F>("cos_theta1kStar_toppt_a_up", "cos(#theta_{1}^{k*}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta1kStar_toppt_a_down         = book<TH1F>("cos_theta1kStar_toppt_a_down", "cos(#theta_{1}^{k*}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta1kStar_toppt_b_up           = book<TH1F>("cos_theta1kStar_toppt_b_up", "cos(#theta_{1}^{k*}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta1kStar_toppt_b_down         = book<TH1F>("cos_theta1kStar_toppt_b_down", "cos(#theta_{1}^{k*}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta1kStar_murmuf_upup          = book<TH1F>("cos_theta1kStar_murmuf_upup", "cos(#theta_{1}^{k*}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta1kStar_murmuf_upnone        = book<TH1F>("cos_theta1kStar_murmuf_upnone", "cos(#theta_{1}^{k*}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta1kStar_murmuf_noneup        = book<TH1F>("cos_theta1kStar_murmuf_noneup", "cos(#theta_{1}^{k*}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta1kStar_murmuf_nonedown      = book<TH1F>("cos_theta1kStar_murmuf_nonedown", "cos(#theta_{1}^{k*}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta1kStar_murmuf_downnone      = book<TH1F>("cos_theta1kStar_murmuf_downnone", "cos(#theta_{1}^{k*}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta1kStar_murmuf_downdown      = book<TH1F>("cos_theta1kStar_murmuf_downdown", "cos(#theta_{1}^{k*}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta1kStar_isr_up               = book<TH1F>("cos_theta1kStar_isr_up", "cos(#theta_{1}^{k*}) isr_up",                               24, -1.0, 1.0);
  // cos_theta1kStar_isr_down             = book<TH1F>("cos_theta1kStar_isr_down", "cos(#theta_{1}^{k*}) isr_down",                           24, -1.0, 1.0);
  // cos_theta1kStar_fsr_up               = book<TH1F>("cos_theta1kStar_fsr_up", "cos(#theta_{1}^{k*}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta1kStar_fsr_down             = book<TH1F>("cos_theta1kStar_fsr_down", "cos(#theta_{1}^{k*}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta1rStar = book<TH1F>("cos_theta1rStar", "cos(#theta_{1}^{r*})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta1rStar_ele_reco_up          = book<TH1F>("cos_theta1rStar_ele_reco_up",   "cos(#theta_{1}^{r*}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta1rStar_ele_reco_down        = book<TH1F>("cos_theta1rStar_ele_reco_down", "cos(#theta_{1}^{r*}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta1rStar_ele_id_up            = book<TH1F>("cos_theta1rStar_ele_id_up",   "cos(#theta_{1}^{r*}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta1rStar_ele_id_down          = book<TH1F>("cos_theta1rStar_ele_id_down", "cos(#theta_{1}^{r*}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta1rStar_ele_trigger_up       = book<TH1F>("cos_theta1rStar_ele_trigger_up",   "cos(#theta_{1}^{r*}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta1rStar_ele_trigger_down     = book<TH1F>("cos_theta1rStar_ele_trigger_down", "cos(#theta_{1}^{r*}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta1rStar_mu_reco_up           = book<TH1F>("cos_theta1rStar_mu_reco_up",   "cos(#theta_{1}^{r*}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta1rStar_mu_reco_down         = book<TH1F>("cos_theta1rStar_mu_reco_down", "cos(#theta_{1}^{r*}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta1rStar_mu_id_stat_up        = book<TH1F>("cos_theta1rStar_mu_id_stat_up",   "cos(#theta_{1}^{r*}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta1rStar_mu_id_stat_down      = book<TH1F>("cos_theta1rStar_mu_id_stat_down",   "cos(#theta_{1}^{r*}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta1rStar_mu_id_syst_up        = book<TH1F>("cos_theta1rStar_mu_id_syst_up",   "cos(#theta_{1}^{r*}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta1rStar_mu_id_syst_down      = book<TH1F>("cos_theta1rStar_mu_id_syst_down",   "cos(#theta_{1}^{r*}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta1rStar_mu_trigger_stat_up   = book<TH1F>("cos_theta1rStar_mu_trigger_stat_up",   "cos(#theta_{1}^{r*}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta1rStar_mu_trigger_stat_down = book<TH1F>("cos_theta1rStar_mu_trigger_stat_down",   "cos(#theta_{1}^{r*}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta1rStar_mu_trigger_syst_up   = book<TH1F>("cos_theta1rStar_mu_trigger_syst_up",   "cos(#theta_{1}^{r*}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta1rStar_mu_trigger_syst_down = book<TH1F>("cos_theta1rStar_mu_trigger_syst_down",   "cos(#theta_{1}^{r*}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta1rStar_mu_iso_stat_up       = book<TH1F>("cos_theta1rStar_mu_iso_stat_up",   "cos(#theta_{1}^{r*}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta1rStar_mu_iso_stat_down     = book<TH1F>("cos_theta1rStar_mu_iso_stat_down",   "cos(#theta_{1}^{r*}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta1rStar_mu_iso_syst_up       = book<TH1F>("cos_theta1rStar_mu_iso_syst_up",   "cos(#theta_{1}^{r*}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta1rStar_mu_iso_syst_down     = book<TH1F>("cos_theta1rStar_mu_iso_syst_down",   "cos(#theta_{1}^{r*}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta1rStar_pu_up                = book<TH1F>("cos_theta1rStar_pu_up",   "cos(#theta_{1}^{r*}) pu_up",                               24, -1.0, 1.0);
  // cos_theta1rStar_pu_down              = book<TH1F>("cos_theta1rStar_pu_down", "cos(#theta_{1}^{r*}) pu_down",                             24, -1.0, 1.0);
  // cos_theta1rStar_prefiring_up         = book<TH1F>("cos_theta1rStar_prefiring_up",   "cos(#theta_{1}^{r*}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta1rStar_prefiring_down       = book<TH1F>("cos_theta1rStar_prefiring_down", "cos(#theta_{1}^{r*}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta1rStar_btag_cferr1_up       = book<TH1F>("cos_theta1rStar_btag_cferr1_up", "cos(#theta_{1}^{r*}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta1rStar_btag_cferr1_down     = book<TH1F>("cos_theta1rStar_btag_cferr1_down", "cos(#theta_{1}^{r*}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta1rStar_btag_cferr2_up       = book<TH1F>("cos_theta1rStar_btag_cferr2_up", "cos(#theta_{1}^{r*}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta1rStar_btag_cferr2_down     = book<TH1F>("cos_theta1rStar_btag_cferr2_down", "cos(#theta_{1}^{r*}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta1rStar_btag_hf_up           = book<TH1F>("cos_theta1rStar_btag_hf_up", "cos(#theta_{1}^{r*}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta1rStar_btag_hf_down         = book<TH1F>("cos_theta1rStar_btag_hf_down", "cos(#theta_{1}^{r*}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta1rStar_btag_hfstats1_up     = book<TH1F>("cos_theta1rStar_btag_hfstats1_up", "cos(#theta_{1}^{r*}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta1rStar_btag_hfstats1_down   = book<TH1F>("cos_theta1rStar_btag_hfstats1_down", "cos(#theta_{1}^{r*}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta1rStar_btag_hfstats2_up     = book<TH1F>("cos_theta1rStar_btag_hfstats2_up", "cos(#theta_{1}^{r*}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta1rStar_btag_hfstats2_down   = book<TH1F>("cos_theta1rStar_btag_hfstats2_down", "cos(#theta_{1}^{r*}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta1rStar_btag_lf_up           = book<TH1F>("cos_theta1rStar_btag_lf_up", "cos(#theta_{1}^{r*}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta1rStar_btag_lf_down         = book<TH1F>("cos_theta1rStar_btag_lf_down", "cos(#theta_{1}^{r*}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta1rStar_btag_lfstats1_up     = book<TH1F>("cos_theta1rStar_btag_lfstats1_up", "cos(#theta_{1}^{r*}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta1rStar_btag_lfstats1_down   = book<TH1F>("cos_theta1rStar_btag_lfstats1_down", "cos(#theta_{1}^{r*}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta1rStar_btag_lfstats2_up     = book<TH1F>("cos_theta1rStar_btag_lfstats2_up", "cos(#theta_{1}^{r*}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta1rStar_btag_lfstats2_down   = book<TH1F>("cos_theta1rStar_btag_lfstats2_down", "cos(#theta_{1}^{r*}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta1rStar_ttag_corr_up         = book<TH1F>("cos_theta1rStar_ttag_corr_up", "cos(#theta_{1}^{r*}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta1rStar_ttag_corr_down       = book<TH1F>("cos_theta1rStar_ttag_corr_down", "cos(#theta_{1}^{r*}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta1rStar_ttag_uncorr_up       = book<TH1F>("cos_theta1rStar_ttag_uncorr_up", "cos(#theta_{1}^{r*}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta1rStar_ttag_uncorr_down     = book<TH1F>("cos_theta1rStar_ttag_uncorr_down", "cos(#theta_{1}^{r*}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta1rStar_tmistag_up           = book<TH1F>("cos_theta1rStar_tmistag_up", "cos(#theta_{1}^{r*}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta1rStar_tmistag_down         = book<TH1F>("cos_theta1rStar_tmistag_down", "cos(#theta_{1}^{r*}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta1rStar_toppt_a_up           = book<TH1F>("cos_theta1rStar_toppt_a_up", "cos(#theta_{1}^{r*}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta1rStar_toppt_a_down         = book<TH1F>("cos_theta1rStar_toppt_a_down", "cos(#theta_{1}^{r*}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta1rStar_toppt_b_up           = book<TH1F>("cos_theta1rStar_toppt_b_up", "cos(#theta_{1}^{r*}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta1rStar_toppt_b_down         = book<TH1F>("cos_theta1rStar_toppt_b_down", "cos(#theta_{1}^{r*}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta1rStar_murmuf_upup          = book<TH1F>("cos_theta1rStar_murmuf_upup", "cos(#theta_{1}^{r*}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta1rStar_murmuf_upnone        = book<TH1F>("cos_theta1rStar_murmuf_upnone", "cos(#theta_{1}^{r*}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta1rStar_murmuf_noneup        = book<TH1F>("cos_theta1rStar_murmuf_noneup", "cos(#theta_{1}^{r*}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta1rStar_murmuf_nonedown      = book<TH1F>("cos_theta1rStar_murmuf_nonedown", "cos(#theta_{1}^{r*}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta1rStar_murmuf_downnone      = book<TH1F>("cos_theta1rStar_murmuf_downnone", "cos(#theta_{1}^{r*}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta1rStar_murmuf_downdown      = book<TH1F>("cos_theta1rStar_murmuf_downdown", "cos(#theta_{1}^{r*}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta1rStar_isr_up               = book<TH1F>("cos_theta1rStar_isr_up", "cos(#theta_{1}^{r*}) isr_up",                               24, -1.0, 1.0);
  // cos_theta1rStar_isr_down             = book<TH1F>("cos_theta1rStar_isr_down", "cos(#theta_{1}^{r*}) isr_down",                           24, -1.0, 1.0);
  // cos_theta1rStar_fsr_up               = book<TH1F>("cos_theta1rStar_fsr_up", "cos(#theta_{1}^{r*}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta1rStar_fsr_down             = book<TH1F>("cos_theta1rStar_fsr_down", "cos(#theta_{1}^{r*}) fsr_down",                           24, -1.0, 1.0);

  // cos_theta2k = book<TH1F>("cos_theta2k", "cos(#theta_{2}^{k})",24, -1, 1);-----------------------------------------------------------------//
  cos_theta2k_ele_reco_up          = book<TH1F>("cos_theta2k_ele_reco_up",   "cos(#theta_{2}^{k}) ele_reco_up",                   24, -1.0, 1.0);
  cos_theta2k_ele_reco_down        = book<TH1F>("cos_theta2k_ele_reco_down", "cos(#theta_{2}^{k}) ele_reco_down",                 24, -1.0, 1.0);
  cos_theta2k_ele_id_up            = book<TH1F>("cos_theta2k_ele_id_up",   "cos(#theta_{2}^{k}) ele_id_up",                       24, -1.0, 1.0);
  cos_theta2k_ele_id_down          = book<TH1F>("cos_theta2k_ele_id_down", "cos(#theta_{2}^{k}) ele_id_down",                     24, -1.0, 1.0);
  cos_theta2k_ele_trigger_up       = book<TH1F>("cos_theta2k_ele_trigger_up",   "cos(#theta_{2}^{k}) ele_trigger_up",             24, -1.0, 1.0);
  cos_theta2k_ele_trigger_down     = book<TH1F>("cos_theta2k_ele_trigger_down", "cos(#theta_{2}^{k}) ele_trigger_down",           24, -1.0, 1.0);
  cos_theta2k_mu_reco_up           = book<TH1F>("cos_theta2k_mu_reco_up",   "cos(#theta_{2}^{k}) mu_reco_up",                     24, -1.0, 1.0);
  cos_theta2k_mu_reco_down         = book<TH1F>("cos_theta2k_mu_reco_down", "cos(#theta_{2}^{k}) mu_reco_down",                   24, -1.0, 1.0);
  cos_theta2k_mu_id_stat_up        = book<TH1F>("cos_theta2k_mu_id_stat_up",   "cos(#theta_{2}^{k}) mu_id_stat_up",               24, -1.0, 1.0);
  cos_theta2k_mu_id_stat_down      = book<TH1F>("cos_theta2k_mu_id_stat_down",   "cos(#theta_{2}^{k}) mu_id_stat_down",           24, -1.0, 1.0);
  cos_theta2k_mu_id_syst_up        = book<TH1F>("cos_theta2k_mu_id_syst_up",   "cos(#theta_{2}^{k}) mu_id_syst_up",               24, -1.0, 1.0);
  cos_theta2k_mu_id_syst_down      = book<TH1F>("cos_theta2k_mu_id_syst_down",   "cos(#theta_{2}^{k}) mu_id_syst_down",           24, -1.0, 1.0);
  cos_theta2k_mu_trigger_stat_up   = book<TH1F>("cos_theta2k_mu_trigger_stat_up",   "cos(#theta_{2}^{k}) mu_trigger_stat_up",     24, -1.0, 1.0);
  cos_theta2k_mu_trigger_stat_down = book<TH1F>("cos_theta2k_mu_trigger_stat_down",   "cos(#theta_{2}^{k}) mu_trigger_stat_down", 24, -1.0, 1.0);
  cos_theta2k_mu_trigger_syst_up   = book<TH1F>("cos_theta2k_mu_trigger_syst_up",   "cos(#theta_{2}^{k}) mu_trigger_syst_up",     24, -1.0, 1.0);
  cos_theta2k_mu_trigger_syst_down = book<TH1F>("cos_theta2k_mu_trigger_syst_down",   "cos(#theta_{2}^{k}) mu_trigger_syst_down", 24, -1.0, 1.0);
  cos_theta2k_mu_iso_stat_up       = book<TH1F>("cos_theta2k_mu_iso_stat_up",   "cos(#theta_{2}^{k}) mu_iso_stat_up",             24, -1.0, 1.0);
  cos_theta2k_mu_iso_stat_down     = book<TH1F>("cos_theta2k_mu_iso_stat_down",   "cos(#theta_{2}^{k}) mu_iso_stat_down",         24, -1.0, 1.0);
  cos_theta2k_mu_iso_syst_up       = book<TH1F>("cos_theta2k_mu_iso_syst_up",   "cos(#theta_{2}^{k}) mu_iso_syst_up",             24, -1.0, 1.0);
  cos_theta2k_mu_iso_syst_down     = book<TH1F>("cos_theta2k_mu_iso_syst_down",   "cos(#theta_{2}^{k}) mu_iso_syst_down",         24, -1.0, 1.0);
  cos_theta2k_pu_up                = book<TH1F>("cos_theta2k_pu_up",   "cos(#theta_{2}^{k}) pu_up",                               24, -1.0, 1.0);
  cos_theta2k_pu_down              = book<TH1F>("cos_theta2k_pu_down", "cos(#theta_{2}^{k}) pu_down",                             24, -1.0, 1.0);
  cos_theta2k_prefiring_up         = book<TH1F>("cos_theta2k_prefiring_up",   "cos(#theta_{2}^{k}) prefiring_up",                 24, -1.0, 1.0);
  cos_theta2k_prefiring_down       = book<TH1F>("cos_theta2k_prefiring_down", "cos(#theta_{2}^{k}) prefiring_down",               24, -1.0, 1.0);
  cos_theta2k_btag_cferr1_up       = book<TH1F>("cos_theta2k_btag_cferr1_up", "cos(#theta_{2}^{k}) btag_cferr1_up",               24, -1.0, 1.0);
  cos_theta2k_btag_cferr1_down     = book<TH1F>("cos_theta2k_btag_cferr1_down", "cos(#theta_{2}^{k}) btag_cferr1_down",           24, -1.0, 1.0);
  cos_theta2k_btag_cferr2_up       = book<TH1F>("cos_theta2k_btag_cferr2_up", "cos(#theta_{2}^{k}) btag_cferr2_up",               24, -1.0, 1.0);
  cos_theta2k_btag_cferr2_down     = book<TH1F>("cos_theta2k_btag_cferr2_down", "cos(#theta_{2}^{k}) btag_cferr2_down",           24, -1.0, 1.0);
  cos_theta2k_btag_hf_up           = book<TH1F>("cos_theta2k_btag_hf_up", "cos(#theta_{2}^{k}) btag_hf_up",                       24, -1.0, 1.0);
  cos_theta2k_btag_hf_down         = book<TH1F>("cos_theta2k_btag_hf_down", "cos(#theta_{2}^{k}) btag_hf_down",                   24, -1.0, 1.0);
  cos_theta2k_btag_hfstats1_up     = book<TH1F>("cos_theta2k_btag_hfstats1_up", "cos(#theta_{2}^{k}) btag_hfstats1_up",           24, -1.0, 1.0);
  cos_theta2k_btag_hfstats1_down   = book<TH1F>("cos_theta2k_btag_hfstats1_down", "cos(#theta_{2}^{k}) btag_hfstats1_down",       24, -1.0, 1.0);
  cos_theta2k_btag_hfstats2_up     = book<TH1F>("cos_theta2k_btag_hfstats2_up", "cos(#theta_{2}^{k}) btag_hfstats2_up",           24, -1.0, 1.0);
  cos_theta2k_btag_hfstats2_down   = book<TH1F>("cos_theta2k_btag_hfstats2_down", "cos(#theta_{2}^{k}) btag_hfstats2_down",       24, -1.0, 1.0);
  cos_theta2k_btag_lf_up           = book<TH1F>("cos_theta2k_btag_lf_up", "cos(#theta_{2}^{k}) btag_lf_up",                       24, -1.0, 1.0);
  cos_theta2k_btag_lf_down         = book<TH1F>("cos_theta2k_btag_lf_down", "cos(#theta_{2}^{k}) btag_lf_down",                   24, -1.0, 1.0);
  cos_theta2k_btag_lfstats1_up     = book<TH1F>("cos_theta2k_btag_lfstats1_up", "cos(#theta_{2}^{k}) btag_lfstats1_up",           24, -1.0, 1.0);
  cos_theta2k_btag_lfstats1_down   = book<TH1F>("cos_theta2k_btag_lfstats1_down", "cos(#theta_{2}^{k}) btag_lfstats1_down",       24, -1.0, 1.0);
  cos_theta2k_btag_lfstats2_up     = book<TH1F>("cos_theta2k_btag_lfstats2_up", "cos(#theta_{2}^{k}) btag_lfstats2_up",           24, -1.0, 1.0);
  cos_theta2k_btag_lfstats2_down   = book<TH1F>("cos_theta2k_btag_lfstats2_down", "cos(#theta_{2}^{k}) btag_lfstats2_down",       24, -1.0, 1.0);
  cos_theta2k_ttag_corr_up         = book<TH1F>("cos_theta2k_ttag_corr_up", "cos(#theta_{2}^{k}) ttag_corr_up",                   24, -1.0, 1.0);
  cos_theta2k_ttag_corr_down       = book<TH1F>("cos_theta2k_ttag_corr_down", "cos(#theta_{2}^{k}) ttag_corr_down",               24, -1.0, 1.0);
  cos_theta2k_ttag_uncorr_up       = book<TH1F>("cos_theta2k_ttag_uncorr_up", "cos(#theta_{2}^{k}) ttag_uncorr_up",               24, -1.0, 1.0);
  cos_theta2k_ttag_uncorr_down     = book<TH1F>("cos_theta2k_ttag_uncorr_down", "cos(#theta_{2}^{k}) ttag_counrr_down",           24, -1.0, 1.0);
  cos_theta2k_tmistag_up           = book<TH1F>("cos_theta2k_tmistag_up", "cos(#theta_{2}^{k}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cos_theta2k_tmistag_down         = book<TH1F>("cos_theta2k_tmistag_down", "cos(#theta_{2}^{k}) [GeV] tmistag_down",             24, -1.0, 1.0);
  cos_theta2k_toppt_a_up           = book<TH1F>("cos_theta2k_toppt_a_up", "cos(#theta_{2}^{k}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cos_theta2k_toppt_a_down         = book<TH1F>("cos_theta2k_toppt_a_down", "cos(#theta_{2}^{k}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cos_theta2k_toppt_b_up           = book<TH1F>("cos_theta2k_toppt_b_up", "cos(#theta_{2}^{k}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cos_theta2k_toppt_b_down         = book<TH1F>("cos_theta2k_toppt_b_down", "cos(#theta_{2}^{k}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cos_theta2k_murmuf_upup          = book<TH1F>("cos_theta2k_murmuf_upup", "cos(#theta_{2}^{k}) murmuf_upup",                     24, -1.0, 1.0);
  cos_theta2k_murmuf_upnone        = book<TH1F>("cos_theta2k_murmuf_upnone", "cos(#theta_{2}^{k}) murmuf_upnone",                 24, -1.0, 1.0);
  cos_theta2k_murmuf_noneup        = book<TH1F>("cos_theta2k_murmuf_noneup", "cos(#theta_{2}^{k}) murmuf_noneup",                 24, -1.0, 1.0);
  cos_theta2k_murmuf_nonedown      = book<TH1F>("cos_theta2k_murmuf_nonedown", "cos(#theta_{2}^{k}) murmuf_nonedown",             24, -1.0, 1.0);
  cos_theta2k_murmuf_downnone      = book<TH1F>("cos_theta2k_murmuf_downnone", "cos(#theta_{2}^{k}) murmuf_downnone",             24, -1.0, 1.0);
  cos_theta2k_murmuf_downdown      = book<TH1F>("cos_theta2k_murmuf_downdown", "cos(#theta_{2}^{k}) murmuf_downdown",             24, -1.0, 1.0);
  cos_theta2k_isr_up               = book<TH1F>("cos_theta2k_isr_up", "cos(#theta_{2}^{k}) isr_up",                               24, -1.0, 1.0);
  cos_theta2k_isr_down             = book<TH1F>("cos_theta2k_isr_down", "cos(#theta_{2}^{k}) isr_down",                           24, -1.0, 1.0);
  cos_theta2k_fsr_up               = book<TH1F>("cos_theta2k_fsr_up", "cos(#theta_{2}^{k}) fsr_up",                               24, -1.0, 1.0);
  cos_theta2k_fsr_down             = book<TH1F>("cos_theta2k_fsr_down", "cos(#theta_{2}^{k}) fsr_down",                           24, -1.0, 1.0);

  // cos_theta2r = book<TH1F>("cos_theta2r", "cos(#theta_{2}^{r})",24, -1, 1);-----------------------------------------------------------------//
  cos_theta2r_ele_reco_up          = book<TH1F>("cos_theta2r_ele_reco_up",   "cos(#theta_{2}^{r}) ele_reco_up",                   24, -1.0, 1.0);
  cos_theta2r_ele_reco_down        = book<TH1F>("cos_theta2r_ele_reco_down", "cos(#theta_{2}^{r}) ele_reco_down",                 24, -1.0, 1.0);
  cos_theta2r_ele_id_up            = book<TH1F>("cos_theta2r_ele_id_up",   "cos(#theta_{2}^{r}) ele_id_up",                       24, -1.0, 1.0);
  cos_theta2r_ele_id_down          = book<TH1F>("cos_theta2r_ele_id_down", "cos(#theta_{2}^{r}) ele_id_down",                     24, -1.0, 1.0);
  cos_theta2r_ele_trigger_up       = book<TH1F>("cos_theta2r_ele_trigger_up",   "cos(#theta_{2}^{r}) ele_trigger_up",             24, -1.0, 1.0);
  cos_theta2r_ele_trigger_down     = book<TH1F>("cos_theta2r_ele_trigger_down", "cos(#theta_{2}^{r}) ele_trigger_down",           24, -1.0, 1.0);
  cos_theta2r_mu_reco_up           = book<TH1F>("cos_theta2r_mu_reco_up",   "cos(#theta_{2}^{r}) mu_reco_up",                     24, -1.0, 1.0);
  cos_theta2r_mu_reco_down         = book<TH1F>("cos_theta2r_mu_reco_down", "cos(#theta_{2}^{r}) mu_reco_down",                   24, -1.0, 1.0);
  cos_theta2r_mu_id_stat_up        = book<TH1F>("cos_theta2r_mu_id_stat_up",   "cos(#theta_{2}^{r}) mu_id_stat_up",               24, -1.0, 1.0);
  cos_theta2r_mu_id_stat_down      = book<TH1F>("cos_theta2r_mu_id_stat_down",   "cos(#theta_{2}^{r}) mu_id_stat_down",           24, -1.0, 1.0);
  cos_theta2r_mu_id_syst_up        = book<TH1F>("cos_theta2r_mu_id_syst_up",   "cos(#theta_{2}^{r}) mu_id_syst_up",               24, -1.0, 1.0);
  cos_theta2r_mu_id_syst_down      = book<TH1F>("cos_theta2r_mu_id_syst_down",   "cos(#theta_{2}^{r}) mu_id_syst_down",           24, -1.0, 1.0);
  cos_theta2r_mu_trigger_stat_up   = book<TH1F>("cos_theta2r_mu_trigger_stat_up",   "cos(#theta_{2}^{r}) mu_trigger_stat_up",     24, -1.0, 1.0);
  cos_theta2r_mu_trigger_stat_down = book<TH1F>("cos_theta2r_mu_trigger_stat_down",   "cos(#theta_{2}^{r}) mu_trigger_stat_down", 24, -1.0, 1.0);
  cos_theta2r_mu_trigger_syst_up   = book<TH1F>("cos_theta2r_mu_trigger_syst_up",   "cos(#theta_{2}^{r}) mu_trigger_syst_up",     24, -1.0, 1.0);
  cos_theta2r_mu_trigger_syst_down = book<TH1F>("cos_theta2r_mu_trigger_syst_down",   "cos(#theta_{2}^{r}) mu_trigger_syst_down", 24, -1.0, 1.0);
  cos_theta2r_mu_iso_stat_up       = book<TH1F>("cos_theta2r_mu_iso_stat_up",   "cos(#theta_{2}^{r}) mu_iso_stat_up",             24, -1.0, 1.0);
  cos_theta2r_mu_iso_stat_down     = book<TH1F>("cos_theta2r_mu_iso_stat_down",   "cos(#theta_{2}^{r}) mu_iso_stat_down",         24, -1.0, 1.0);
  cos_theta2r_mu_iso_syst_up       = book<TH1F>("cos_theta2r_mu_iso_syst_up",   "cos(#theta_{2}^{r}) mu_iso_syst_up",             24, -1.0, 1.0);
  cos_theta2r_mu_iso_syst_down     = book<TH1F>("cos_theta2r_mu_iso_syst_down",   "cos(#theta_{2}^{r}) mu_iso_syst_down",         24, -1.0, 1.0);
  cos_theta2r_pu_up                = book<TH1F>("cos_theta2r_pu_up",   "cos(#theta_{2}^{r}) pu_up",                               24, -1.0, 1.0);
  cos_theta2r_pu_down              = book<TH1F>("cos_theta2r_pu_down", "cos(#theta_{2}^{r}) pu_down",                             24, -1.0, 1.0);
  cos_theta2r_prefiring_up         = book<TH1F>("cos_theta2r_prefiring_up",   "cos(#theta_{2}^{r}) prefiring_up",                 24, -1.0, 1.0);
  cos_theta2r_prefiring_down       = book<TH1F>("cos_theta2r_prefiring_down", "cos(#theta_{2}^{r}) prefiring_down",               24, -1.0, 1.0);
  cos_theta2r_btag_cferr1_up       = book<TH1F>("cos_theta2r_btag_cferr1_up", "cos(#theta_{2}^{r}) btag_cferr1_up",               24, -1.0, 1.0);
  cos_theta2r_btag_cferr1_down     = book<TH1F>("cos_theta2r_btag_cferr1_down", "cos(#theta_{2}^{r}) btag_cferr1_down",           24, -1.0, 1.0);
  cos_theta2r_btag_cferr2_up       = book<TH1F>("cos_theta2r_btag_cferr2_up", "cos(#theta_{2}^{r}) btag_cferr2_up",               24, -1.0, 1.0);
  cos_theta2r_btag_cferr2_down     = book<TH1F>("cos_theta2r_btag_cferr2_down", "cos(#theta_{2}^{r}) btag_cferr2_down",           24, -1.0, 1.0);
  cos_theta2r_btag_hf_up           = book<TH1F>("cos_theta2r_btag_hf_up", "cos(#theta_{2}^{r}) btag_hf_up",                       24, -1.0, 1.0);
  cos_theta2r_btag_hf_down         = book<TH1F>("cos_theta2r_btag_hf_down", "cos(#theta_{2}^{r}) btag_hf_down",                   24, -1.0, 1.0);
  cos_theta2r_btag_hfstats1_up     = book<TH1F>("cos_theta2r_btag_hfstats1_up", "cos(#theta_{2}^{r}) btag_hfstats1_up",           24, -1.0, 1.0);
  cos_theta2r_btag_hfstats1_down   = book<TH1F>("cos_theta2r_btag_hfstats1_down", "cos(#theta_{2}^{r}) btag_hfstats1_down",       24, -1.0, 1.0);
  cos_theta2r_btag_hfstats2_up     = book<TH1F>("cos_theta2r_btag_hfstats2_up", "cos(#theta_{2}^{r}) btag_hfstats2_up",           24, -1.0, 1.0);
  cos_theta2r_btag_hfstats2_down   = book<TH1F>("cos_theta2r_btag_hfstats2_down", "cos(#theta_{2}^{r}) btag_hfstats2_down",       24, -1.0, 1.0);
  cos_theta2r_btag_lf_up           = book<TH1F>("cos_theta2r_btag_lf_up", "cos(#theta_{2}^{r}) btag_lf_up",                       24, -1.0, 1.0);
  cos_theta2r_btag_lf_down         = book<TH1F>("cos_theta2r_btag_lf_down", "cos(#theta_{2}^{r}) btag_lf_down",                   24, -1.0, 1.0);
  cos_theta2r_btag_lfstats1_up     = book<TH1F>("cos_theta2r_btag_lfstats1_up", "cos(#theta_{2}^{r}) btag_lfstats1_up",           24, -1.0, 1.0);
  cos_theta2r_btag_lfstats1_down   = book<TH1F>("cos_theta2r_btag_lfstats1_down", "cos(#theta_{2}^{r}) btag_lfstats1_down",       24, -1.0, 1.0);
  cos_theta2r_btag_lfstats2_up     = book<TH1F>("cos_theta2r_btag_lfstats2_up", "cos(#theta_{2}^{r}) btag_lfstats2_up",           24, -1.0, 1.0);
  cos_theta2r_btag_lfstats2_down   = book<TH1F>("cos_theta2r_btag_lfstats2_down", "cos(#theta_{2}^{r}) btag_lfstats2_down",       24, -1.0, 1.0);
  cos_theta2r_ttag_corr_up         = book<TH1F>("cos_theta2r_ttag_corr_up", "cos(#theta_{2}^{r}) ttag_corr_up",                   24, -1.0, 1.0);
  cos_theta2r_ttag_corr_down       = book<TH1F>("cos_theta2r_ttag_corr_down", "cos(#theta_{2}^{r}) ttag_corr_down",               24, -1.0, 1.0);
  cos_theta2r_ttag_uncorr_up       = book<TH1F>("cos_theta2r_ttag_uncorr_up", "cos(#theta_{2}^{r}) ttag_uncorr_up",               24, -1.0, 1.0);
  cos_theta2r_ttag_uncorr_down     = book<TH1F>("cos_theta2r_ttag_uncorr_down", "cos(#theta_{2}^{r}) ttag_counrr_down",           24, -1.0, 1.0);
  cos_theta2r_tmistag_up           = book<TH1F>("cos_theta2r_tmistag_up", "cos(#theta_{2}^{r}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cos_theta2r_tmistag_down         = book<TH1F>("cos_theta2r_tmistag_down", "cos(#theta_{2}^{r}) [GeV] tmistag_down",             24, -1.0, 1.0);
  cos_theta2r_toppt_a_up           = book<TH1F>("cos_theta2r_toppt_a_up", "cos(#theta_{2}^{r}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cos_theta2r_toppt_a_down         = book<TH1F>("cos_theta2r_toppt_a_down", "cos(#theta_{2}^{r}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cos_theta2r_toppt_b_up           = book<TH1F>("cos_theta2r_toppt_b_up", "cos(#theta_{2}^{r}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cos_theta2r_toppt_b_down         = book<TH1F>("cos_theta2r_toppt_b_down", "cos(#theta_{2}^{r}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cos_theta2r_murmuf_upup          = book<TH1F>("cos_theta2r_murmuf_upup", "cos(#theta_{2}^{r}) murmuf_upup",                     24, -1.0, 1.0);
  cos_theta2r_murmuf_upnone        = book<TH1F>("cos_theta2r_murmuf_upnone", "cos(#theta_{2}^{r}) murmuf_upnone",                 24, -1.0, 1.0);
  cos_theta2r_murmuf_noneup        = book<TH1F>("cos_theta2r_murmuf_noneup", "cos(#theta_{2}^{r}) murmuf_noneup",                 24, -1.0, 1.0);
  cos_theta2r_murmuf_nonedown      = book<TH1F>("cos_theta2r_murmuf_nonedown", "cos(#theta_{2}^{r}) murmuf_nonedown",             24, -1.0, 1.0);
  cos_theta2r_murmuf_downnone      = book<TH1F>("cos_theta2r_murmuf_downnone", "cos(#theta_{2}^{r}) murmuf_downnone",             24, -1.0, 1.0);
  cos_theta2r_murmuf_downdown      = book<TH1F>("cos_theta2r_murmuf_downdown", "cos(#theta_{2}^{r}) murmuf_downdown",             24, -1.0, 1.0);
  cos_theta2r_isr_up               = book<TH1F>("cos_theta2r_isr_up", "cos(#theta_{2}^{r}) isr_up",                               24, -1.0, 1.0);
  cos_theta2r_isr_down             = book<TH1F>("cos_theta2r_isr_down", "cos(#theta_{2}^{r}) isr_down",                           24, -1.0, 1.0);
  cos_theta2r_fsr_up               = book<TH1F>("cos_theta2r_fsr_up", "cos(#theta_{2}^{r}) fsr_up",                               24, -1.0, 1.0);
  cos_theta2r_fsr_down             = book<TH1F>("cos_theta2r_fsr_down", "cos(#theta_{2}^{r}) fsr_down",                           24, -1.0, 1.0);

  // cos_theta2n = book<TH1F>("cos_theta2n", "cos(#theta_{2}^{n})",24, -1, 1);-----------------------------------------------------------------//
  cos_theta2n_ele_reco_up          = book<TH1F>("cos_theta2n_ele_reco_up",   "cos(#theta_{2}^{n}) ele_reco_up",                   24, -1.0, 1.0);
  cos_theta2n_ele_reco_down        = book<TH1F>("cos_theta2n_ele_reco_down", "cos(#theta_{2}^{n}) ele_reco_down",                 24, -1.0, 1.0);
  cos_theta2n_ele_id_up            = book<TH1F>("cos_theta2n_ele_id_up",   "cos(#theta_{2}^{n}) ele_id_up",                       24, -1.0, 1.0);
  cos_theta2n_ele_id_down          = book<TH1F>("cos_theta2n_ele_id_down", "cos(#theta_{2}^{n}) ele_id_down",                     24, -1.0, 1.0);
  cos_theta2n_ele_trigger_up       = book<TH1F>("cos_theta2n_ele_trigger_up",   "cos(#theta_{2}^{n}) ele_trigger_up",             24, -1.0, 1.0);
  cos_theta2n_ele_trigger_down     = book<TH1F>("cos_theta2n_ele_trigger_down", "cos(#theta_{2}^{n}) ele_trigger_down",           24, -1.0, 1.0);
  cos_theta2n_mu_reco_up           = book<TH1F>("cos_theta2n_mu_reco_up",   "cos(#theta_{2}^{n}) mu_reco_up",                     24, -1.0, 1.0);
  cos_theta2n_mu_reco_down         = book<TH1F>("cos_theta2n_mu_reco_down", "cos(#theta_{2}^{n}) mu_reco_down",                   24, -1.0, 1.0);
  cos_theta2n_mu_id_stat_up        = book<TH1F>("cos_theta2n_mu_id_stat_up",   "cos(#theta_{2}^{n}) mu_id_stat_up",               24, -1.0, 1.0);
  cos_theta2n_mu_id_stat_down      = book<TH1F>("cos_theta2n_mu_id_stat_down",   "cos(#theta_{2}^{n}) mu_id_stat_down",           24, -1.0, 1.0);
  cos_theta2n_mu_id_syst_up        = book<TH1F>("cos_theta2n_mu_id_syst_up",   "cos(#theta_{2}^{n}) mu_id_syst_up",               24, -1.0, 1.0);
  cos_theta2n_mu_id_syst_down      = book<TH1F>("cos_theta2n_mu_id_syst_down",   "cos(#theta_{2}^{n}) mu_id_syst_down",           24, -1.0, 1.0);
  cos_theta2n_mu_trigger_stat_up   = book<TH1F>("cos_theta2n_mu_trigger_stat_up",   "cos(#theta_{2}^{n}) mu_trigger_stat_up",     24, -1.0, 1.0);
  cos_theta2n_mu_trigger_stat_down = book<TH1F>("cos_theta2n_mu_trigger_stat_down",   "cos(#theta_{2}^{n}) mu_trigger_stat_down", 24, -1.0, 1.0);
  cos_theta2n_mu_trigger_syst_up   = book<TH1F>("cos_theta2n_mu_trigger_syst_up",   "cos(#theta_{2}^{n}) mu_trigger_syst_up",     24, -1.0, 1.0);
  cos_theta2n_mu_trigger_syst_down = book<TH1F>("cos_theta2n_mu_trigger_syst_down",   "cos(#theta_{2}^{n}) mu_trigger_syst_down", 24, -1.0, 1.0);
  cos_theta2n_mu_iso_stat_up       = book<TH1F>("cos_theta2n_mu_iso_stat_up",   "cos(#theta_{2}^{n}) mu_iso_stat_up",             24, -1.0, 1.0);
  cos_theta2n_mu_iso_stat_down     = book<TH1F>("cos_theta2n_mu_iso_stat_down",   "cos(#theta_{2}^{n}) mu_iso_stat_down",         24, -1.0, 1.0);
  cos_theta2n_mu_iso_syst_up       = book<TH1F>("cos_theta2n_mu_iso_syst_up",   "cos(#theta_{2}^{n}) mu_iso_syst_up",             24, -1.0, 1.0);
  cos_theta2n_mu_iso_syst_down     = book<TH1F>("cos_theta2n_mu_iso_syst_down",   "cos(#theta_{2}^{n}) mu_iso_syst_down",         24, -1.0, 1.0);
  cos_theta2n_pu_up                = book<TH1F>("cos_theta2n_pu_up",   "cos(#theta_{2}^{n}) pu_up",                               24, -1.0, 1.0);
  cos_theta2n_pu_down              = book<TH1F>("cos_theta2n_pu_down", "cos(#theta_{2}^{n}) pu_down",                             24, -1.0, 1.0);
  cos_theta2n_prefiring_up         = book<TH1F>("cos_theta2n_prefiring_up",   "cos(#theta_{2}^{n}) prefiring_up",                 24, -1.0, 1.0);
  cos_theta2n_prefiring_down       = book<TH1F>("cos_theta2n_prefiring_down", "cos(#theta_{2}^{n}) prefiring_down",               24, -1.0, 1.0);
  cos_theta2n_btag_cferr1_up       = book<TH1F>("cos_theta2n_btag_cferr1_up", "cos(#theta_{2}^{n}) btag_cferr1_up",               24, -1.0, 1.0);
  cos_theta2n_btag_cferr1_down     = book<TH1F>("cos_theta2n_btag_cferr1_down", "cos(#theta_{2}^{n}) btag_cferr1_down",           24, -1.0, 1.0);
  cos_theta2n_btag_cferr2_up       = book<TH1F>("cos_theta2n_btag_cferr2_up", "cos(#theta_{2}^{n}) btag_cferr2_up",               24, -1.0, 1.0);
  cos_theta2n_btag_cferr2_down     = book<TH1F>("cos_theta2n_btag_cferr2_down", "cos(#theta_{2}^{n}) btag_cferr2_down",           24, -1.0, 1.0);
  cos_theta2n_btag_hf_up           = book<TH1F>("cos_theta2n_btag_hf_up", "cos(#theta_{2}^{n}) btag_hf_up",                       24, -1.0, 1.0);
  cos_theta2n_btag_hf_down         = book<TH1F>("cos_theta2n_btag_hf_down", "cos(#theta_{2}^{n}) btag_hf_down",                   24, -1.0, 1.0);
  cos_theta2n_btag_hfstats1_up     = book<TH1F>("cos_theta2n_btag_hfstats1_up", "cos(#theta_{2}^{n}) btag_hfstats1_up",           24, -1.0, 1.0);
  cos_theta2n_btag_hfstats1_down   = book<TH1F>("cos_theta2n_btag_hfstats1_down", "cos(#theta_{2}^{n}) btag_hfstats1_down",       24, -1.0, 1.0);
  cos_theta2n_btag_hfstats2_up     = book<TH1F>("cos_theta2n_btag_hfstats2_up", "cos(#theta_{2}^{n}) btag_hfstats2_up",           24, -1.0, 1.0);
  cos_theta2n_btag_hfstats2_down   = book<TH1F>("cos_theta2n_btag_hfstats2_down", "cos(#theta_{2}^{n}) btag_hfstats2_down",       24, -1.0, 1.0);
  cos_theta2n_btag_lf_up           = book<TH1F>("cos_theta2n_btag_lf_up", "cos(#theta_{2}^{n}) btag_lf_up",                       24, -1.0, 1.0);
  cos_theta2n_btag_lf_down         = book<TH1F>("cos_theta2n_btag_lf_down", "cos(#theta_{2}^{n}) btag_lf_down",                   24, -1.0, 1.0);
  cos_theta2n_btag_lfstats1_up     = book<TH1F>("cos_theta2n_btag_lfstats1_up", "cos(#theta_{2}^{n}) btag_lfstats1_up",           24, -1.0, 1.0);
  cos_theta2n_btag_lfstats1_down   = book<TH1F>("cos_theta2n_btag_lfstats1_down", "cos(#theta_{2}^{n}) btag_lfstats1_down",       24, -1.0, 1.0);
  cos_theta2n_btag_lfstats2_up     = book<TH1F>("cos_theta2n_btag_lfstats2_up", "cos(#theta_{2}^{n}) btag_lfstats2_up",           24, -1.0, 1.0);
  cos_theta2n_btag_lfstats2_down   = book<TH1F>("cos_theta2n_btag_lfstats2_down", "cos(#theta_{2}^{n}) btag_lfstats2_down",       24, -1.0, 1.0);
  cos_theta2n_ttag_corr_up         = book<TH1F>("cos_theta2n_ttag_corr_up", "cos(#theta_{2}^{n}) ttag_corr_up",                   24, -1.0, 1.0);
  cos_theta2n_ttag_corr_down       = book<TH1F>("cos_theta2n_ttag_corr_down", "cos(#theta_{2}^{n}) ttag_corr_down",               24, -1.0, 1.0);
  cos_theta2n_ttag_uncorr_up       = book<TH1F>("cos_theta2n_ttag_uncorr_up", "cos(#theta_{2}^{n}) ttag_uncorr_up",               24, -1.0, 1.0);
  cos_theta2n_ttag_uncorr_down     = book<TH1F>("cos_theta2n_ttag_uncorr_down", "cos(#theta_{2}^{n}) ttag_counrr_down",           24, -1.0, 1.0);
  cos_theta2n_tmistag_up           = book<TH1F>("cos_theta2n_tmistag_up", "cos(#theta_{2}^{n}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cos_theta2n_tmistag_down         = book<TH1F>("cos_theta2n_tmistag_down", "cos(#theta_{2}^{n}) [GeV] tmistag_down",             24, -1.0, 1.0);
  cos_theta2n_toppt_a_up           = book<TH1F>("cos_theta2n_toppt_a_up", "cos(#theta_{2}^{n}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cos_theta2n_toppt_a_down         = book<TH1F>("cos_theta2n_toppt_a_down", "cos(#theta_{2}^{n}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cos_theta2n_toppt_b_up           = book<TH1F>("cos_theta2n_toppt_b_up", "cos(#theta_{2}^{n}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cos_theta2n_toppt_b_down         = book<TH1F>("cos_theta2n_toppt_b_down", "cos(#theta_{2}^{n}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cos_theta2n_murmuf_upup          = book<TH1F>("cos_theta2n_murmuf_upup", "cos(#theta_{2}^{n}) murmuf_upup",                     24, -1.0, 1.0);
  cos_theta2n_murmuf_upnone        = book<TH1F>("cos_theta2n_murmuf_upnone", "cos(#theta_{2}^{n}) murmuf_upnone",                 24, -1.0, 1.0);
  cos_theta2n_murmuf_noneup        = book<TH1F>("cos_theta2n_murmuf_noneup", "cos(#theta_{2}^{n}) murmuf_noneup",                 24, -1.0, 1.0);
  cos_theta2n_murmuf_nonedown      = book<TH1F>("cos_theta2n_murmuf_nonedown", "cos(#theta_{2}^{n}) murmuf_nonedown",             24, -1.0, 1.0);
  cos_theta2n_murmuf_downnone      = book<TH1F>("cos_theta2n_murmuf_downnone", "cos(#theta_{2}^{n}) murmuf_downnone",             24, -1.0, 1.0);
  cos_theta2n_murmuf_downdown      = book<TH1F>("cos_theta2n_murmuf_downdown", "cos(#theta_{2}^{n}) murmuf_downdown",             24, -1.0, 1.0);
  cos_theta2n_isr_up               = book<TH1F>("cos_theta2n_isr_up", "cos(#theta_{2}^{n}) isr_up",                               24, -1.0, 1.0);
  cos_theta2n_isr_down             = book<TH1F>("cos_theta2n_isr_down", "cos(#theta_{2}^{n}) isr_down",                           24, -1.0, 1.0);
  cos_theta2n_fsr_up               = book<TH1F>("cos_theta2n_fsr_up", "cos(#theta_{2}^{n}) fsr_up",                               24, -1.0, 1.0);
  cos_theta2n_fsr_down             = book<TH1F>("cos_theta2n_fsr_down", "cos(#theta_{2}^{n}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta2kStar = book<TH1F>("cos_theta2kStar", "cos(#theta_{2}^{k*})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta2kStar_ele_reco_up          = book<TH1F>("cos_theta2kStar_ele_reco_up",   "cos(#theta_{2}^{k*}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta2kStar_ele_reco_down        = book<TH1F>("cos_theta2kStar_ele_reco_down", "cos(#theta_{2}^{k*}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta2kStar_ele_id_up            = book<TH1F>("cos_theta2kStar_ele_id_up",   "cos(#theta_{2}^{k*}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta2kStar_ele_id_down          = book<TH1F>("cos_theta2kStar_ele_id_down", "cos(#theta_{2}^{k*}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta2kStar_ele_trigger_up       = book<TH1F>("cos_theta2kStar_ele_trigger_up",   "cos(#theta_{2}^{k*}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta2kStar_ele_trigger_down     = book<TH1F>("cos_theta2kStar_ele_trigger_down", "cos(#theta_{2}^{k*}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta2kStar_mu_reco_up           = book<TH1F>("cos_theta2kStar_mu_reco_up",   "cos(#theta_{2}^{k*}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta2kStar_mu_reco_down         = book<TH1F>("cos_theta2kStar_mu_reco_down", "cos(#theta_{2}^{k*}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta2kStar_mu_id_stat_up        = book<TH1F>("cos_theta2kStar_mu_id_stat_up",   "cos(#theta_{2}^{k*}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta2kStar_mu_id_stat_down      = book<TH1F>("cos_theta2kStar_mu_id_stat_down",   "cos(#theta_{2}^{k*}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta2kStar_mu_id_syst_up        = book<TH1F>("cos_theta2kStar_mu_id_syst_up",   "cos(#theta_{2}^{k*}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta2kStar_mu_id_syst_down      = book<TH1F>("cos_theta2kStar_mu_id_syst_down",   "cos(#theta_{2}^{k*}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta2kStar_mu_trigger_stat_up   = book<TH1F>("cos_theta2kStar_mu_trigger_stat_up",   "cos(#theta_{2}^{k*}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta2kStar_mu_trigger_stat_down = book<TH1F>("cos_theta2kStar_mu_trigger_stat_down",   "cos(#theta_{2}^{k*}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta2kStar_mu_trigger_syst_up   = book<TH1F>("cos_theta2kStar_mu_trigger_syst_up",   "cos(#theta_{2}^{k*}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta2kStar_mu_trigger_syst_down = book<TH1F>("cos_theta2kStar_mu_trigger_syst_down",   "cos(#theta_{2}^{k*}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta2kStar_mu_iso_stat_up       = book<TH1F>("cos_theta2kStar_mu_iso_stat_up",   "cos(#theta_{2}^{k*}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta2kStar_mu_iso_stat_down     = book<TH1F>("cos_theta2kStar_mu_iso_stat_down",   "cos(#theta_{2}^{k*}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta2kStar_mu_iso_syst_up       = book<TH1F>("cos_theta2kStar_mu_iso_syst_up",   "cos(#theta_{2}^{k*}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta2kStar_mu_iso_syst_down     = book<TH1F>("cos_theta2kStar_mu_iso_syst_down",   "cos(#theta_{2}^{k*}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta2kStar_pu_up                = book<TH1F>("cos_theta2kStar_pu_up",   "cos(#theta_{2}^{k*}) pu_up",                               24, -1.0, 1.0);
  // cos_theta2kStar_pu_down              = book<TH1F>("cos_theta2kStar_pu_down", "cos(#theta_{2}^{k*}) pu_down",                             24, -1.0, 1.0);
  // cos_theta2kStar_prefiring_up         = book<TH1F>("cos_theta2kStar_prefiring_up",   "cos(#theta_{2}^{k*}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta2kStar_prefiring_down       = book<TH1F>("cos_theta2kStar_prefiring_down", "cos(#theta_{2}^{k*}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta2kStar_btag_cferr1_up       = book<TH1F>("cos_theta2kStar_btag_cferr1_up", "cos(#theta_{2}^{k*}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta2kStar_btag_cferr1_down     = book<TH1F>("cos_theta2kStar_btag_cferr1_down", "cos(#theta_{2}^{k*}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta2kStar_btag_cferr2_up       = book<TH1F>("cos_theta2kStar_btag_cferr2_up", "cos(#theta_{2}^{k*}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta2kStar_btag_cferr2_down     = book<TH1F>("cos_theta2kStar_btag_cferr2_down", "cos(#theta_{2}^{k*}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta2kStar_btag_hf_up           = book<TH1F>("cos_theta2kStar_btag_hf_up", "cos(#theta_{2}^{k*}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta2kStar_btag_hf_down         = book<TH1F>("cos_theta2kStar_btag_hf_down", "cos(#theta_{2}^{k*}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta2kStar_btag_hfstats1_up     = book<TH1F>("cos_theta2kStar_btag_hfstats1_up", "cos(#theta_{2}^{k*}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta2kStar_btag_hfstats1_down   = book<TH1F>("cos_theta2kStar_btag_hfstats1_down", "cos(#theta_{2}^{k*}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta2kStar_btag_hfstats2_up     = book<TH1F>("cos_theta2kStar_btag_hfstats2_up", "cos(#theta_{2}^{k*}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta2kStar_btag_hfstats2_down   = book<TH1F>("cos_theta2kStar_btag_hfstats2_down", "cos(#theta_{2}^{k*}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta2kStar_btag_lf_up           = book<TH1F>("cos_theta2kStar_btag_lf_up", "cos(#theta_{2}^{k*}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta2kStar_btag_lf_down         = book<TH1F>("cos_theta2kStar_btag_lf_down", "cos(#theta_{2}^{k*}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta2kStar_btag_lfstats1_up     = book<TH1F>("cos_theta2kStar_btag_lfstats1_up", "cos(#theta_{2}^{k*}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta2kStar_btag_lfstats1_down   = book<TH1F>("cos_theta2kStar_btag_lfstats1_down", "cos(#theta_{2}^{k*}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta2kStar_btag_lfstats2_up     = book<TH1F>("cos_theta2kStar_btag_lfstats2_up", "cos(#theta_{2}^{k*}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta2kStar_btag_lfstats2_down   = book<TH1F>("cos_theta2kStar_btag_lfstats2_down", "cos(#theta_{2}^{k*}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta2kStar_ttag_corr_up         = book<TH1F>("cos_theta2kStar_ttag_corr_up", "cos(#theta_{2}^{k*}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta2kStar_ttag_corr_down       = book<TH1F>("cos_theta2kStar_ttag_corr_down", "cos(#theta_{2}^{k*}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta2kStar_ttag_uncorr_up       = book<TH1F>("cos_theta2kStar_ttag_uncorr_up", "cos(#theta_{2}^{k*}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta2kStar_ttag_uncorr_down     = book<TH1F>("cos_theta2kStar_ttag_uncorr_down", "cos(#theta_{2}^{k*}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta2kStar_tmistag_up           = book<TH1F>("cos_theta2kStar_tmistag_up", "cos(#theta_{2}^{k*}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta2kStar_tmistag_down         = book<TH1F>("cos_theta2kStar_tmistag_down", "cos(#theta_{2}^{k*}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta2kStar_toppt_a_up           = book<TH1F>("cos_theta2kStar_toppt_a_up", "cos(#theta_{2}^{k*}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta2kStar_toppt_a_down         = book<TH1F>("cos_theta2kStar_toppt_a_down", "cos(#theta_{2}^{k*}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta2kStar_toppt_b_up           = book<TH1F>("cos_theta2kStar_toppt_b_up", "cos(#theta_{2}^{k*}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta2kStar_toppt_b_down         = book<TH1F>("cos_theta2kStar_toppt_b_down", "cos(#theta_{2}^{k*}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta2kStar_murmuf_upup          = book<TH1F>("cos_theta2kStar_murmuf_upup", "cos(#theta_{2}^{k*}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta2kStar_murmuf_upnone        = book<TH1F>("cos_theta2kStar_murmuf_upnone", "cos(#theta_{2}^{k*}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta2kStar_murmuf_noneup        = book<TH1F>("cos_theta2kStar_murmuf_noneup", "cos(#theta_{2}^{k*}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta2kStar_murmuf_nonedown      = book<TH1F>("cos_theta2kStar_murmuf_nonedown", "cos(#theta_{2}^{k*}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta2kStar_murmuf_downnone      = book<TH1F>("cos_theta2kStar_murmuf_downnone", "cos(#theta_{2}^{k*}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta2kStar_murmuf_downdown      = book<TH1F>("cos_theta2kStar_murmuf_downdown", "cos(#theta_{2}^{k*}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta2kStar_isr_up               = book<TH1F>("cos_theta2kStar_isr_up", "cos(#theta_{2}^{k*}) isr_up",                               24, -1.0, 1.0);
  // cos_theta2kStar_isr_down             = book<TH1F>("cos_theta2kStar_isr_down", "cos(#theta_{2}^{k*}) isr_down",                           24, -1.0, 1.0);
  // cos_theta2kStar_fsr_up               = book<TH1F>("cos_theta2kStar_fsr_up", "cos(#theta_{2}^{k*}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta2kStar_fsr_down             = book<TH1F>("cos_theta2kStar_fsr_down", "cos(#theta_{2}^{k*}) fsr_down",                           24, -1.0, 1.0);

  // // cos_theta2rStar = book<TH1F>("cos_theta2rStar", "cos(#theta_{2}^{r*})",24, -1, 1);-----------------------------------------------------------------//
  // cos_theta2rStar_ele_reco_up          = book<TH1F>("cos_theta2rStar_ele_reco_up",   "cos(#theta_{2}^{r*}) ele_reco_up",                   24, -1.0, 1.0);
  // cos_theta2rStar_ele_reco_down        = book<TH1F>("cos_theta2rStar_ele_reco_down", "cos(#theta_{2}^{r*}) ele_reco_down",                 24, -1.0, 1.0);
  // cos_theta2rStar_ele_id_up            = book<TH1F>("cos_theta2rStar_ele_id_up",   "cos(#theta_{2}^{r*}) ele_id_up",                       24, -1.0, 1.0);
  // cos_theta2rStar_ele_id_down          = book<TH1F>("cos_theta2rStar_ele_id_down", "cos(#theta_{2}^{r*}) ele_id_down",                     24, -1.0, 1.0);
  // cos_theta2rStar_ele_trigger_up       = book<TH1F>("cos_theta2rStar_ele_trigger_up",   "cos(#theta_{2}^{r*}) ele_trigger_up",             24, -1.0, 1.0);
  // cos_theta2rStar_ele_trigger_down     = book<TH1F>("cos_theta2rStar_ele_trigger_down", "cos(#theta_{2}^{r*}) ele_trigger_down",           24, -1.0, 1.0);
  // cos_theta2rStar_mu_reco_up           = book<TH1F>("cos_theta2rStar_mu_reco_up",   "cos(#theta_{2}^{r*}) mu_reco_up",                     24, -1.0, 1.0);
  // cos_theta2rStar_mu_reco_down         = book<TH1F>("cos_theta2rStar_mu_reco_down", "cos(#theta_{2}^{r*}) mu_reco_down",                   24, -1.0, 1.0);
  // cos_theta2rStar_mu_id_stat_up        = book<TH1F>("cos_theta2rStar_mu_id_stat_up",   "cos(#theta_{2}^{r*}) mu_id_stat_up",               24, -1.0, 1.0);
  // cos_theta2rStar_mu_id_stat_down      = book<TH1F>("cos_theta2rStar_mu_id_stat_down",   "cos(#theta_{2}^{r*}) mu_id_stat_down",           24, -1.0, 1.0);
  // cos_theta2rStar_mu_id_syst_up        = book<TH1F>("cos_theta2rStar_mu_id_syst_up",   "cos(#theta_{2}^{r*}) mu_id_syst_up",               24, -1.0, 1.0);
  // cos_theta2rStar_mu_id_syst_down      = book<TH1F>("cos_theta2rStar_mu_id_syst_down",   "cos(#theta_{2}^{r*}) mu_id_syst_down",           24, -1.0, 1.0);
  // cos_theta2rStar_mu_trigger_stat_up   = book<TH1F>("cos_theta2rStar_mu_trigger_stat_up",   "cos(#theta_{2}^{r*}) mu_trigger_stat_up",     24, -1.0, 1.0);
  // cos_theta2rStar_mu_trigger_stat_down = book<TH1F>("cos_theta2rStar_mu_trigger_stat_down",   "cos(#theta_{2}^{r*}) mu_trigger_stat_down", 24, -1.0, 1.0);
  // cos_theta2rStar_mu_trigger_syst_up   = book<TH1F>("cos_theta2rStar_mu_trigger_syst_up",   "cos(#theta_{2}^{r*}) mu_trigger_syst_up",     24, -1.0, 1.0);
  // cos_theta2rStar_mu_trigger_syst_down = book<TH1F>("cos_theta2rStar_mu_trigger_syst_down",   "cos(#theta_{2}^{r*}) mu_trigger_syst_down", 24, -1.0, 1.0);
  // cos_theta2rStar_mu_iso_stat_up       = book<TH1F>("cos_theta2rStar_mu_iso_stat_up",   "cos(#theta_{2}^{r*}) mu_iso_stat_up",             24, -1.0, 1.0);
  // cos_theta2rStar_mu_iso_stat_down     = book<TH1F>("cos_theta2rStar_mu_iso_stat_down",   "cos(#theta_{2}^{r*}) mu_iso_stat_down",         24, -1.0, 1.0);
  // cos_theta2rStar_mu_iso_syst_up       = book<TH1F>("cos_theta2rStar_mu_iso_syst_up",   "cos(#theta_{2}^{r*}) mu_iso_syst_up",             24, -1.0, 1.0);
  // cos_theta2rStar_mu_iso_syst_down     = book<TH1F>("cos_theta2rStar_mu_iso_syst_down",   "cos(#theta_{2}^{r*}) mu_iso_syst_down",         24, -1.0, 1.0);
  // cos_theta2rStar_pu_up                = book<TH1F>("cos_theta2rStar_pu_up",   "cos(#theta_{2}^{r*}) pu_up",                               24, -1.0, 1.0);
  // cos_theta2rStar_pu_down              = book<TH1F>("cos_theta2rStar_pu_down", "cos(#theta_{2}^{r*}) pu_down",                             24, -1.0, 1.0);
  // cos_theta2rStar_prefiring_up         = book<TH1F>("cos_theta2rStar_prefiring_up",   "cos(#theta_{2}^{r*}) prefiring_up",                 24, -1.0, 1.0);
  // cos_theta2rStar_prefiring_down       = book<TH1F>("cos_theta2rStar_prefiring_down", "cos(#theta_{2}^{r*}) prefiring_down",               24, -1.0, 1.0);
  // cos_theta2rStar_btag_cferr1_up       = book<TH1F>("cos_theta2rStar_btag_cferr1_up", "cos(#theta_{2}^{r*}) btag_cferr1_up",               24, -1.0, 1.0);
  // cos_theta2rStar_btag_cferr1_down     = book<TH1F>("cos_theta2rStar_btag_cferr1_down", "cos(#theta_{2}^{r*}) btag_cferr1_down",           24, -1.0, 1.0);
  // cos_theta2rStar_btag_cferr2_up       = book<TH1F>("cos_theta2rStar_btag_cferr2_up", "cos(#theta_{2}^{r*}) btag_cferr2_up",               24, -1.0, 1.0);
  // cos_theta2rStar_btag_cferr2_down     = book<TH1F>("cos_theta2rStar_btag_cferr2_down", "cos(#theta_{2}^{r*}) btag_cferr2_down",           24, -1.0, 1.0);
  // cos_theta2rStar_btag_hf_up           = book<TH1F>("cos_theta2rStar_btag_hf_up", "cos(#theta_{2}^{r*}) btag_hf_up",                       24, -1.0, 1.0);
  // cos_theta2rStar_btag_hf_down         = book<TH1F>("cos_theta2rStar_btag_hf_down", "cos(#theta_{2}^{r*}) btag_hf_down",                   24, -1.0, 1.0);
  // cos_theta2rStar_btag_hfstats1_up     = book<TH1F>("cos_theta2rStar_btag_hfstats1_up", "cos(#theta_{2}^{r*}) btag_hfstats1_up",           24, -1.0, 1.0);
  // cos_theta2rStar_btag_hfstats1_down   = book<TH1F>("cos_theta2rStar_btag_hfstats1_down", "cos(#theta_{2}^{r*}) btag_hfstats1_down",       24, -1.0, 1.0);
  // cos_theta2rStar_btag_hfstats2_up     = book<TH1F>("cos_theta2rStar_btag_hfstats2_up", "cos(#theta_{2}^{r*}) btag_hfstats2_up",           24, -1.0, 1.0);
  // cos_theta2rStar_btag_hfstats2_down   = book<TH1F>("cos_theta2rStar_btag_hfstats2_down", "cos(#theta_{2}^{r*}) btag_hfstats2_down",       24, -1.0, 1.0);
  // cos_theta2rStar_btag_lf_up           = book<TH1F>("cos_theta2rStar_btag_lf_up", "cos(#theta_{2}^{r*}) btag_lf_up",                       24, -1.0, 1.0);
  // cos_theta2rStar_btag_lf_down         = book<TH1F>("cos_theta2rStar_btag_lf_down", "cos(#theta_{2}^{r*}) btag_lf_down",                   24, -1.0, 1.0);
  // cos_theta2rStar_btag_lfstats1_up     = book<TH1F>("cos_theta2rStar_btag_lfstats1_up", "cos(#theta_{2}^{r*}) btag_lfstats1_up",           24, -1.0, 1.0);
  // cos_theta2rStar_btag_lfstats1_down   = book<TH1F>("cos_theta2rStar_btag_lfstats1_down", "cos(#theta_{2}^{r*}) btag_lfstats1_down",       24, -1.0, 1.0);
  // cos_theta2rStar_btag_lfstats2_up     = book<TH1F>("cos_theta2rStar_btag_lfstats2_up", "cos(#theta_{2}^{r*}) btag_lfstats2_up",           24, -1.0, 1.0);
  // cos_theta2rStar_btag_lfstats2_down   = book<TH1F>("cos_theta2rStar_btag_lfstats2_down", "cos(#theta_{2}^{r*}) btag_lfstats2_down",       24, -1.0, 1.0);
  // cos_theta2rStar_ttag_corr_up         = book<TH1F>("cos_theta2rStar_ttag_corr_up", "cos(#theta_{2}^{r*}) ttag_corr_up",                   24, -1.0, 1.0);
  // cos_theta2rStar_ttag_corr_down       = book<TH1F>("cos_theta2rStar_ttag_corr_down", "cos(#theta_{2}^{r*}) ttag_corr_down",               24, -1.0, 1.0);
  // cos_theta2rStar_ttag_uncorr_up       = book<TH1F>("cos_theta2rStar_ttag_uncorr_up", "cos(#theta_{2}^{r*}) ttag_uncorr_up",               24, -1.0, 1.0);
  // cos_theta2rStar_ttag_uncorr_down     = book<TH1F>("cos_theta2rStar_ttag_uncorr_down", "cos(#theta_{2}^{r*}) ttag_counrr_down",           24, -1.0, 1.0);
  // cos_theta2rStar_tmistag_up           = book<TH1F>("cos_theta2rStar_tmistag_up", "cos(#theta_{2}^{r*}) [GeV] tmistag_up",                 24, -1.0, 1.0);
  // cos_theta2rStar_tmistag_down         = book<TH1F>("cos_theta2rStar_tmistag_down", "cos(#theta_{2}^{r*}) [GeV] tmistag_down",             24, -1.0, 1.0);
  // cos_theta2rStar_toppt_a_up           = book<TH1F>("cos_theta2rStar_toppt_a_up", "cos(#theta_{2}^{r*}) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  // cos_theta2rStar_toppt_a_down         = book<TH1F>("cos_theta2rStar_toppt_a_down", "cos(#theta_{2}^{r*}) [GeV] toppt_a_down",             24, -1.0, 1.0);
  // cos_theta2rStar_toppt_b_up           = book<TH1F>("cos_theta2rStar_toppt_b_up", "cos(#theta_{2}^{r*}) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  // cos_theta2rStar_toppt_b_down         = book<TH1F>("cos_theta2rStar_toppt_b_down", "cos(#theta_{2}^{r*}) [GeV] toppt_b_down",             24, -1.0, 1.0);
  // cos_theta2rStar_murmuf_upup          = book<TH1F>("cos_theta2rStar_murmuf_upup", "cos(#theta_{2}^{r*}) murmuf_upup",                     24, -1.0, 1.0);
  // cos_theta2rStar_murmuf_upnone        = book<TH1F>("cos_theta2rStar_murmuf_upnone", "cos(#theta_{2}^{r*}) murmuf_upnone",                 24, -1.0, 1.0);
  // cos_theta2rStar_murmuf_noneup        = book<TH1F>("cos_theta2rStar_murmuf_noneup", "cos(#theta_{2}^{r*}) murmuf_noneup",                 24, -1.0, 1.0);
  // cos_theta2rStar_murmuf_nonedown      = book<TH1F>("cos_theta2rStar_murmuf_nonedown", "cos(#theta_{2}^{r*}) murmuf_nonedown",             24, -1.0, 1.0);
  // cos_theta2rStar_murmuf_downnone      = book<TH1F>("cos_theta2rStar_murmuf_downnone", "cos(#theta_{2}^{r*}) murmuf_downnone",             24, -1.0, 1.0);
  // cos_theta2rStar_murmuf_downdown      = book<TH1F>("cos_theta2rStar_murmuf_downdown", "cos(#theta_{2}^{r*}) murmuf_downdown",             24, -1.0, 1.0);
  // cos_theta2rStar_isr_up               = book<TH1F>("cos_theta2rStar_isr_up", "cos(#theta_{2}^{r*}) isr_up",                               24, -1.0, 1.0);
  // cos_theta2rStar_isr_down             = book<TH1F>("cos_theta2rStar_isr_down", "cos(#theta_{2}^{r*}) isr_down",                           24, -1.0, 1.0);
  // cos_theta2rStar_fsr_up               = book<TH1F>("cos_theta2rStar_fsr_up", "cos(#theta_{2}^{r*}) fsr_up",                               24, -1.0, 1.0);
  // cos_theta2rStar_fsr_down             = book<TH1F>("cos_theta2rStar_fsr_down", "cos(#theta_{2}^{r*}) fsr_down",                           24, -1.0, 1.0);

  // Ckk = book<TH1F>("Ckk", "C_{kk}",24, -1, 1);-----------------------------------------------------------------//


  // Crr = book<TH1F>("Crr", "C_{rr}",24, -1, 1);-----------------------------------------------------------------//


  // Cnn = book<TH1F>("Cnn", "C_{nn}",24, -1, 1);-----------------------------------------------------------------//


  // Crk_plus = book<TH1F>("Crk_plus",  "C_{rk} + C_{kr}",24, -1, 1);------------------------------------------------------------------//


  // Crk_minus = book<TH1F>("Crk_minus", "C_{rk} - C_{kr}",24, -1, 1);-----------------------------------------------------------------//


  // Cnr_plus = book<TH1F>("Cnr_plus",  "C_{nr} + C_{rn}",24, -1, 1);------------------------------------------------------------------//


  // Cnr_minus = book<TH1F>("Cnr_minus", "C_{nr} - C_{rn}",24, -1, 1);-----------------------------------------------------------------//


  // Cnk_plus = book<TH1F>("Cnk_plus",  "C_{nk} + C_{kn}",24, -1, 1);------------------------------------------------------------------//


  // Cnk_minus = book<TH1F>("Cnk_minus", "C_{nk} - C_{kn}",24, -1, 1);-----------------------------------------------------------------//


  // cHel_Mtt300_400 = book<TH1F>("cHel_Mtt300_400",           "cos(#phi_{lb}) (M_{tt} [300,400] GeV)",24, -1, 1);-----------------------------------------------------------------//
  cHel_Mtt300_400_ele_reco_up          = book<TH1F>("cHel_Mtt300_400_ele_reco_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ele_reco_up",                   24, -1.0, 1.0);
  cHel_Mtt300_400_ele_reco_down        = book<TH1F>("cHel_Mtt300_400_ele_reco_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ele_reco_down",                 24, -1.0, 1.0);
  cHel_Mtt300_400_ele_id_up            = book<TH1F>("cHel_Mtt300_400_ele_id_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ele_id_up",                       24, -1.0, 1.0);
  cHel_Mtt300_400_ele_id_down          = book<TH1F>("cHel_Mtt300_400_ele_id_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ele_id_down",                     24, -1.0, 1.0);
  cHel_Mtt300_400_ele_trigger_up       = book<TH1F>("cHel_Mtt300_400_ele_trigger_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ele_trigger_up",             24, -1.0, 1.0);
  cHel_Mtt300_400_ele_trigger_down     = book<TH1F>("cHel_Mtt300_400_ele_trigger_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ele_trigger_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_mu_reco_up           = book<TH1F>("cHel_Mtt300_400_mu_reco_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_reco_up",                     24, -1.0, 1.0);
  cHel_Mtt300_400_mu_reco_down         = book<TH1F>("cHel_Mtt300_400_mu_reco_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_reco_down",                   24, -1.0, 1.0);
  cHel_Mtt300_400_mu_id_stat_up        = book<TH1F>("cHel_Mtt300_400_mu_id_stat_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_id_stat_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_mu_id_stat_down      = book<TH1F>("cHel_Mtt300_400_mu_id_stat_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_id_stat_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_mu_id_syst_up        = book<TH1F>("cHel_Mtt300_400_mu_id_syst_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_id_syst_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_mu_id_syst_down      = book<TH1F>("cHel_Mtt300_400_mu_id_syst_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_id_syst_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_mu_trigger_stat_up   = book<TH1F>("cHel_Mtt300_400_mu_trigger_stat_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_trigger_stat_up",     24, -1.0, 1.0);
  cHel_Mtt300_400_mu_trigger_stat_down = book<TH1F>("cHel_Mtt300_400_mu_trigger_stat_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_trigger_stat_down", 24, -1.0, 1.0);
  cHel_Mtt300_400_mu_trigger_syst_up   = book<TH1F>("cHel_Mtt300_400_mu_trigger_syst_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_trigger_syst_up",     24, -1.0, 1.0);
  cHel_Mtt300_400_mu_trigger_syst_down = book<TH1F>("cHel_Mtt300_400_mu_trigger_syst_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_trigger_syst_down", 24, -1.0, 1.0);
  cHel_Mtt300_400_mu_iso_stat_up       = book<TH1F>("cHel_Mtt300_400_mu_iso_stat_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_iso_stat_up",             24, -1.0, 1.0);
  cHel_Mtt300_400_mu_iso_stat_down     = book<TH1F>("cHel_Mtt300_400_mu_iso_stat_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_iso_stat_down",         24, -1.0, 1.0);
  cHel_Mtt300_400_mu_iso_syst_up       = book<TH1F>("cHel_Mtt300_400_mu_iso_syst_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_iso_syst_up",             24, -1.0, 1.0);
  cHel_Mtt300_400_mu_iso_syst_down     = book<TH1F>("cHel_Mtt300_400_mu_iso_syst_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) mu_iso_syst_down",         24, -1.0, 1.0);
  cHel_Mtt300_400_pu_up                = book<TH1F>("cHel_Mtt300_400_pu_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) pu_up",                               24, -1.0, 1.0);
  cHel_Mtt300_400_pu_down              = book<TH1F>("cHel_Mtt300_400_pu_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) pu_down",                             24, -1.0, 1.0);
  cHel_Mtt300_400_prefiring_up         = book<TH1F>("cHel_Mtt300_400_prefiring_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV) prefiring_up",                 24, -1.0, 1.0);
  cHel_Mtt300_400_prefiring_down       = book<TH1F>("cHel_Mtt300_400_prefiring_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) prefiring_down",               24, -1.0, 1.0);
  cHel_Mtt300_400_btag_cferr1_up       = book<TH1F>("cHel_Mtt300_400_btag_cferr1_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_cferr1_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_btag_cferr1_down     = book<TH1F>("cHel_Mtt300_400_btag_cferr1_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_cferr1_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_btag_cferr2_up       = book<TH1F>("cHel_Mtt300_400_btag_cferr2_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_cferr2_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_btag_cferr2_down     = book<TH1F>("cHel_Mtt300_400_btag_cferr2_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_cferr2_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_btag_hf_up           = book<TH1F>("cHel_Mtt300_400_btag_hf_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_hf_up",                       24, -1.0, 1.0);
  cHel_Mtt300_400_btag_hf_down         = book<TH1F>("cHel_Mtt300_400_btag_hf_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_hf_down",                   24, -1.0, 1.0);
  cHel_Mtt300_400_btag_hfstats1_up     = book<TH1F>("cHel_Mtt300_400_btag_hfstats1_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_hfstats1_up",           24, -1.0, 1.0);
  cHel_Mtt300_400_btag_hfstats1_down   = book<TH1F>("cHel_Mtt300_400_btag_hfstats1_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_hfstats1_down",       24, -1.0, 1.0);
  cHel_Mtt300_400_btag_hfstats2_up     = book<TH1F>("cHel_Mtt300_400_btag_hfstats2_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_hfstats2_up",           24, -1.0, 1.0);
  cHel_Mtt300_400_btag_hfstats2_down   = book<TH1F>("cHel_Mtt300_400_btag_hfstats2_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_hfstats2_down",       24, -1.0, 1.0);
  cHel_Mtt300_400_btag_lf_up           = book<TH1F>("cHel_Mtt300_400_btag_lf_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_lf_up",                       24, -1.0, 1.0);
  cHel_Mtt300_400_btag_lf_down         = book<TH1F>("cHel_Mtt300_400_btag_lf_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_lf_down",                   24, -1.0, 1.0);
  cHel_Mtt300_400_btag_lfstats1_up     = book<TH1F>("cHel_Mtt300_400_btag_lfstats1_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_lfstats1_up",           24, -1.0, 1.0);
  cHel_Mtt300_400_btag_lfstats1_down   = book<TH1F>("cHel_Mtt300_400_btag_lfstats1_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_lfstats1_down",       24, -1.0, 1.0);
  cHel_Mtt300_400_btag_lfstats2_up     = book<TH1F>("cHel_Mtt300_400_btag_lfstats2_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_lfstats2_up",           24, -1.0, 1.0);
  cHel_Mtt300_400_btag_lfstats2_down   = book<TH1F>("cHel_Mtt300_400_btag_lfstats2_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) btag_lfstats2_down",       24, -1.0, 1.0);
  cHel_Mtt300_400_ttag_corr_up         = book<TH1F>("cHel_Mtt300_400_ttag_corr_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ttag_corr_up",                   24, -1.0, 1.0);
  cHel_Mtt300_400_ttag_corr_down       = book<TH1F>("cHel_Mtt300_400_ttag_corr_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ttag_corr_down",               24, -1.0, 1.0);
  cHel_Mtt300_400_ttag_uncorr_up       = book<TH1F>("cHel_Mtt300_400_ttag_uncorr_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ttag_uncorr_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_ttag_uncorr_down     = book<TH1F>("cHel_Mtt300_400_ttag_uncorr_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) ttag_counrr_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_tmistag_up           = book<TH1F>("cHel_Mtt300_400_tmistag_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cHel_Mtt300_400_tmistag_down         = book<TH1F>("cHel_Mtt300_400_tmistag_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) [GeV] tmistag_down",             24, -1.0, 1.0);
  cHel_Mtt300_400_toppt_a_up           = book<TH1F>("cHel_Mtt300_400_toppt_a_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cHel_Mtt300_400_toppt_a_down         = book<TH1F>("cHel_Mtt300_400_toppt_a_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cHel_Mtt300_400_toppt_b_up           = book<TH1F>("cHel_Mtt300_400_toppt_b_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cHel_Mtt300_400_toppt_b_down         = book<TH1F>("cHel_Mtt300_400_toppt_b_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cHel_Mtt300_400_murmuf_upup          = book<TH1F>("cHel_Mtt300_400_murmuf_upup", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) murmuf_upup",                     24, -1.0, 1.0);
  cHel_Mtt300_400_murmuf_upnone        = book<TH1F>("cHel_Mtt300_400_murmuf_upnone", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) murmuf_upnone",                 24, -1.0, 1.0);
  cHel_Mtt300_400_murmuf_noneup        = book<TH1F>("cHel_Mtt300_400_murmuf_noneup", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) murmuf_noneup",                 24, -1.0, 1.0);
  cHel_Mtt300_400_murmuf_nonedown      = book<TH1F>("cHel_Mtt300_400_murmuf_nonedown", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) murmuf_nonedown",             24, -1.0, 1.0);
  cHel_Mtt300_400_murmuf_downnone      = book<TH1F>("cHel_Mtt300_400_murmuf_downnone", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) murmuf_downnone",             24, -1.0, 1.0);
  cHel_Mtt300_400_murmuf_downdown      = book<TH1F>("cHel_Mtt300_400_murmuf_downdown", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) murmuf_downdown",             24, -1.0, 1.0);
  cHel_Mtt300_400_isr_up               = book<TH1F>("cHel_Mtt300_400_isr_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) isr_up",                               24, -1.0, 1.0);
  cHel_Mtt300_400_isr_down             = book<TH1F>("cHel_Mtt300_400_isr_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) isr_down",                           24, -1.0, 1.0);
  cHel_Mtt300_400_fsr_up               = book<TH1F>("cHel_Mtt300_400_fsr_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) fsr_up",                               24, -1.0, 1.0);
  cHel_Mtt300_400_fsr_down             = book<TH1F>("cHel_Mtt300_400_fsr_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV) fsr_down",                           24, -1.0, 1.0);

  // cHel_Mtt300_400_betaLT0p9 = book<TH1F>("cHel_Mtt300_400_betaLT0p9", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9)",24, -1, 1);-------------------------------------//
  cHel_Mtt300_400_betaLT0p9_ele_reco_up          = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ele_reco_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ele_reco_up",                   24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_ele_reco_down        = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ele_reco_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ele_reco_down",                 24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_ele_id_up            = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ele_id_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ele_id_up",                       24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_ele_id_down          = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ele_id_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ele_id_down",                     24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_ele_trigger_up       = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ele_trigger_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ele_trigger_up",             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_ele_trigger_down     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ele_trigger_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ele_trigger_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_reco_up           = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_reco_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_reco_up",                     24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_reco_down         = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_reco_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_reco_down",                   24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_id_stat_up        = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_id_stat_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_id_stat_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_id_stat_down      = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_id_stat_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_id_stat_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_id_syst_up        = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_id_syst_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_id_syst_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_id_syst_down      = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_id_syst_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_id_syst_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_trigger_stat_up   = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_trigger_stat_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_trigger_stat_up",     24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_trigger_stat_down = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_trigger_stat_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_trigger_stat_down", 24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_trigger_syst_up   = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_trigger_syst_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_trigger_syst_up",     24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_trigger_syst_down = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_trigger_syst_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_trigger_syst_down", 24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_iso_stat_up       = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_iso_stat_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_iso_stat_up",             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_iso_stat_down     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_iso_stat_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_iso_stat_down",         24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_iso_syst_up       = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_iso_syst_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_iso_syst_up",             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_mu_iso_syst_down     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_mu_iso_syst_down",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) mu_iso_syst_down",         24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_pu_up                = book<TH1F>("cHel_Mtt300_400_betaLT0p9_pu_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) pu_up",                               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_pu_down              = book<TH1F>("cHel_Mtt300_400_betaLT0p9_pu_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) pu_down",                             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_prefiring_up         = book<TH1F>("cHel_Mtt300_400_betaLT0p9_prefiring_up",   "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) prefiring_up",                 24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_prefiring_down       = book<TH1F>("cHel_Mtt300_400_betaLT0p9_prefiring_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) prefiring_down",               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_cferr1_up       = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_cferr1_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_cferr1_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_cferr1_down     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_cferr1_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_cferr1_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_cferr2_up       = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_cferr2_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_cferr2_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_cferr2_down     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_cferr2_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_cferr2_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_hf_up           = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_hf_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_hf_up",                       24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_hf_down         = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_hf_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_hf_down",                   24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_hfstats1_up     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_hfstats1_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_hfstats1_up",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_hfstats1_down   = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_hfstats1_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_hfstats1_down",       24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_hfstats2_up     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_hfstats2_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_hfstats2_up",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_hfstats2_down   = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_hfstats2_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_hfstats2_down",       24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_lf_up           = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_lf_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_lf_up",                       24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_lf_down         = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_lf_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_lf_down",                   24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_lfstats1_up     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_lfstats1_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_lfstats1_up",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_lfstats1_down   = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_lfstats1_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_lfstats1_down",       24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_lfstats2_up     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_lfstats2_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_lfstats2_up",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_btag_lfstats2_down   = book<TH1F>("cHel_Mtt300_400_betaLT0p9_btag_lfstats2_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) btag_lfstats2_down",       24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_ttag_corr_up         = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ttag_corr_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ttag_corr_up",                   24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_ttag_corr_down       = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ttag_corr_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ttag_corr_down",               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_ttag_uncorr_up       = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ttag_uncorr_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ttag_uncorr_up",               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_ttag_uncorr_down     = book<TH1F>("cHel_Mtt300_400_betaLT0p9_ttag_uncorr_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) ttag_counrr_down",           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_tmistag_up           = book<TH1F>("cHel_Mtt300_400_betaLT0p9_tmistag_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_tmistag_down         = book<TH1F>("cHel_Mtt300_400_betaLT0p9_tmistag_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) [GeV] tmistag_down",             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_toppt_a_up           = book<TH1F>("cHel_Mtt300_400_betaLT0p9_toppt_a_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_toppt_a_down         = book<TH1F>("cHel_Mtt300_400_betaLT0p9_toppt_a_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_toppt_b_up           = book<TH1F>("cHel_Mtt300_400_betaLT0p9_toppt_b_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_toppt_b_down         = book<TH1F>("cHel_Mtt300_400_betaLT0p9_toppt_b_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_murmuf_upup          = book<TH1F>("cHel_Mtt300_400_betaLT0p9_murmuf_upup", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) murmuf_upup",                     24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_murmuf_upnone        = book<TH1F>("cHel_Mtt300_400_betaLT0p9_murmuf_upnone", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) murmuf_upnone",                 24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_murmuf_noneup        = book<TH1F>("cHel_Mtt300_400_betaLT0p9_murmuf_noneup", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) murmuf_noneup",                 24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_murmuf_nonedown      = book<TH1F>("cHel_Mtt300_400_betaLT0p9_murmuf_nonedown", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) murmuf_nonedown",             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_murmuf_downnone      = book<TH1F>("cHel_Mtt300_400_betaLT0p9_murmuf_downnone", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) murmuf_downnone",             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_murmuf_downdown      = book<TH1F>("cHel_Mtt300_400_betaLT0p9_murmuf_downdown", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) murmuf_downdown",             24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_isr_up               = book<TH1F>("cHel_Mtt300_400_betaLT0p9_isr_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) isr_up",                               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_isr_down             = book<TH1F>("cHel_Mtt300_400_betaLT0p9_isr_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) isr_down",                           24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_fsr_up               = book<TH1F>("cHel_Mtt300_400_betaLT0p9_fsr_up", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) fsr_up",                               24, -1.0, 1.0);
  cHel_Mtt300_400_betaLT0p9_fsr_down             = book<TH1F>("cHel_Mtt300_400_betaLT0p9_fsr_down", "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) fsr_down",                           24, -1.0, 1.0);

  // cHel_P3n_Mtt800_Inf = book<TH1F>("cHel_P3n_Mtt800_Inf",               "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV)",24, -1, 1);-----------------------------------------------------------------//
  cHel_P3n_Mtt800_Inf_ele_reco_up          = book<TH1F>("cHel_P3n_Mtt800_Inf_ele_reco_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ele_reco_up",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_ele_reco_down        = book<TH1F>("cHel_P3n_Mtt800_Inf_ele_reco_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ele_reco_down",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_ele_id_up            = book<TH1F>("cHel_P3n_Mtt800_Inf_ele_id_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ele_id_up",                       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_ele_id_down          = book<TH1F>("cHel_P3n_Mtt800_Inf_ele_id_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ele_id_down",                     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_ele_trigger_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_ele_trigger_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ele_trigger_up",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_ele_trigger_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_ele_trigger_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ele_trigger_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_reco_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_reco_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_reco_up",                     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_reco_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_reco_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_reco_down",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_id_stat_up        = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_id_stat_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_id_stat_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_id_stat_down      = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_id_stat_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_id_stat_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_id_syst_up        = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_id_syst_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_id_syst_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_id_syst_down      = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_id_syst_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_id_syst_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_trigger_stat_up   = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_trigger_stat_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_trigger_stat_up",     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_trigger_stat_down = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_trigger_stat_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_trigger_stat_down", 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_trigger_syst_up   = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_trigger_syst_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_trigger_syst_up",     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_trigger_syst_down = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_trigger_syst_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_trigger_syst_down", 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_iso_stat_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_iso_stat_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_iso_stat_up",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_iso_stat_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_iso_stat_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_iso_stat_down",         24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_iso_syst_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_iso_syst_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_iso_syst_up",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_mu_iso_syst_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_mu_iso_syst_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) mu_iso_syst_down",         24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_pu_up                = book<TH1F>("cHel_P3n_Mtt800_Inf_pu_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) pu_up",                               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_pu_down              = book<TH1F>("cHel_P3n_Mtt800_Inf_pu_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) pu_down",                             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_prefiring_up         = book<TH1F>("cHel_P3n_Mtt800_Inf_prefiring_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) prefiring_up",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_prefiring_down       = book<TH1F>("cHel_P3n_Mtt800_Inf_prefiring_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) prefiring_down",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_cferr1_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_cferr1_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_cferr1_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_cferr1_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_cferr1_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_cferr1_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_cferr2_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_cferr2_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_cferr2_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_cferr2_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_cferr2_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_cferr2_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_hf_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_hf_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_hf_up",                       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_hf_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_hf_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_hf_down",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_hfstats1_up     = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_hfstats1_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_hfstats1_up",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_hfstats1_down   = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_hfstats1_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_hfstats1_down",       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_hfstats2_up     = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_hfstats2_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_hfstats2_up",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_hfstats2_down   = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_hfstats2_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_hfstats2_down",       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_lf_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_lf_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_lf_up",                       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_lf_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_lf_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_lf_down",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_lfstats1_up     = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_lfstats1_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_lfstats1_up",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_lfstats1_down   = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_lfstats1_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_lfstats1_down",       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_lfstats2_up     = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_lfstats2_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_lfstats2_up",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_btag_lfstats2_down   = book<TH1F>("cHel_P3n_Mtt800_Inf_btag_lfstats2_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) btag_lfstats2_down",       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_ttag_corr_up         = book<TH1F>("cHel_P3n_Mtt800_Inf_ttag_corr_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ttag_corr_up",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_ttag_corr_down       = book<TH1F>("cHel_P3n_Mtt800_Inf_ttag_corr_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ttag_corr_down",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_ttag_uncorr_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_ttag_uncorr_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ttag_uncorr_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_ttag_uncorr_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_ttag_uncorr_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) ttag_counrr_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_tmistag_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_tmistag_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_tmistag_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_tmistag_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) [GeV] tmistag_down",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_toppt_a_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_toppt_a_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_toppt_a_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_toppt_a_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_toppt_b_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_toppt_b_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_toppt_b_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_toppt_b_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_murmuf_upup          = book<TH1F>("cHel_P3n_Mtt800_Inf_murmuf_upup", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) murmuf_upup",                     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_murmuf_upnone        = book<TH1F>("cHel_P3n_Mtt800_Inf_murmuf_upnone", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) murmuf_upnone",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_murmuf_noneup        = book<TH1F>("cHel_P3n_Mtt800_Inf_murmuf_noneup", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) murmuf_noneup",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_murmuf_nonedown      = book<TH1F>("cHel_P3n_Mtt800_Inf_murmuf_nonedown", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) murmuf_nonedown",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_murmuf_downnone      = book<TH1F>("cHel_P3n_Mtt800_Inf_murmuf_downnone", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) murmuf_downnone",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_murmuf_downdown      = book<TH1F>("cHel_P3n_Mtt800_Inf_murmuf_downdown", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) murmuf_downdown",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_isr_up               = book<TH1F>("cHel_P3n_Mtt800_Inf_isr_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) isr_up",                               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_isr_down             = book<TH1F>("cHel_P3n_Mtt800_Inf_isr_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) isr_down",                           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_fsr_up               = book<TH1F>("cHel_P3n_Mtt800_Inf_fsr_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) fsr_up",                               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_fsr_down             = book<TH1F>("cHel_P3n_Mtt800_Inf_fsr_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) fsr_down",                           24, -1.0, 1.0);

  // cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4)",24, -1, 1);------------------------------//
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_reco_up          = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_reco_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ele_reco_up",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_reco_down        = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_reco_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ele_reco_down",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_id_up            = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_id_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ele_id_up",                       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_id_down          = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_id_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ele_id_down",                     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_trigger_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_trigger_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ele_trigger_up",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_trigger_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_trigger_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ele_trigger_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_reco_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_reco_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_reco_up",                     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_reco_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_reco_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_reco_down",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_stat_up        = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_stat_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_id_stat_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_stat_down      = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_stat_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_id_stat_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_syst_up        = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_syst_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_id_syst_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_syst_down      = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_syst_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_id_syst_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_stat_up   = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_stat_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_trigger_stat_up",     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_stat_down = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_stat_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_trigger_stat_down", 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_syst_up   = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_syst_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_trigger_syst_up",     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_syst_down = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_syst_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_trigger_syst_down", 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_stat_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_stat_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_iso_stat_up",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_stat_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_stat_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_iso_stat_down",         24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_syst_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_syst_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_iso_syst_up",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_syst_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_syst_down",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) mu_iso_syst_down",         24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_pu_up                = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_pu_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) pu_up",                               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_pu_down              = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_pu_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) pu_down",                             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_prefiring_up         = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_prefiring_up",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) prefiring_up",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_prefiring_down       = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_prefiring_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) prefiring_down",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr1_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr1_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_cferr1_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr1_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr1_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_cferr1_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr2_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr2_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_cferr2_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr2_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr2_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_cferr2_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hf_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hf_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_hf_up",                       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hf_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hf_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_hf_down",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats1_up     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats1_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_hfstats1_up",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats1_down   = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats1_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_hfstats1_down",       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats2_up     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats2_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_hfstats2_up",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats2_down   = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats2_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_hfstats2_down",       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lf_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lf_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_lf_up",                       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lf_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lf_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_lf_down",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats1_up     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats1_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_lfstats1_up",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats1_down   = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats1_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_lfstats1_down",       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats2_up     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats2_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_lfstats2_up",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats2_down   = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats2_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) btag_lfstats2_down",       24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_corr_up         = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_corr_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ttag_corr_up",                   24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_corr_down       = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_corr_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ttag_corr_down",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_uncorr_up       = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_uncorr_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ttag_uncorr_up",               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_uncorr_down     = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_uncorr_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) ttag_counrr_down",           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_tmistag_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_tmistag_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) [GeV] tmistag_up",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_tmistag_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_tmistag_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) [GeV] tmistag_down",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_a_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_a_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) [GeV] toppt_a_up",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_a_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_a_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) [GeV] toppt_a_down",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_b_up           = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_b_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) [GeV] toppt_b_up",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_b_down         = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_b_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) [GeV] toppt_b_down",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_upup          = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_upup", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) murmuf_upup",                     24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_upnone        = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_upnone", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) murmuf_upnone",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_noneup        = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_noneup", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) murmuf_noneup",                 24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_nonedown      = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_nonedown", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) murmuf_nonedown",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_downnone      = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_downnone", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) murmuf_downnone",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_downdown      = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_downdown", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) murmuf_downdown",             24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_isr_up               = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_isr_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) isr_up",                               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_isr_down             = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_isr_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) isr_down",                           24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_fsr_up               = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_fsr_up", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) fsr_up",                               24, -1.0, 1.0);
  cHel_P3n_Mtt800_Inf_cosThetaLT0p4_fsr_down             = book<TH1F>("cHel_P3n_Mtt800_Inf_cosThetaLT0p4_fsr_down", "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) fsr_down",                           24, -1.0, 1.0);


  // Sigma_phi_1                    = book<TH1F>("Sigma_phi_1",   "#Sigma_phi_1_{t#bar{t}} ",                                     16, -3.2, 3.2);
  Sigma_phi_1_mu_reco_up         = book<TH1F>("Sigma_phi_1_mu_reco_up",   "#Sigma_phi_1_{t#bar{t}} mu_reco_up",                16, -3.2, 3.2);
  Sigma_phi_1_mu_reco_down       = book<TH1F>("Sigma_phi_1_mu_reco_down", "#Sigma_phi_1_{t#bar{t}} mu_reco_down",              16,-3.2,3.2);
  Sigma_phi_1_pu_up              = book<TH1F>("Sigma_phi_1_pu_up",   "#Sigma_phi_1_{t#bar{t}} pu_up",                          16,-3.2,3.2);
  Sigma_phi_1_pu_down            = book<TH1F>("Sigma_phi_1_pu_down", "#Sigma_phi_1_{t#bar{t}} pu_down",                        16,-3.2,3.2);
  Sigma_phi_1_prefiring_up       = book<TH1F>("Sigma_phi_1_prefiring_up",   "#Sigma_phi_1_{t#bar{t}} prefiring_up",            16,-3.2,3.2);
  Sigma_phi_1_prefiring_down     = book<TH1F>("Sigma_phi_1_prefiring_down", "#Sigma_phi_1_{t#bar{t}} prefiring_down",          16,-3.2,3.2);
  Sigma_phi_1_mu_id_stat_up      = book<TH1F>("Sigma_phi_1_mu_id_stat_up",   "#Sigma_phi_1_{t#bar{t}} mu_id_stat_up",          16,-3.2,3.2);
  Sigma_phi_1_mu_id_stat_down    = book<TH1F>("Sigma_phi_1_mu_id_stat_down",   "#Sigma_phi_1_{t#bar{t}} mu_id_stat_down",      16,-3.2,3.2);
  Sigma_phi_1_mu_id_syst_up      = book<TH1F>("Sigma_phi_1_mu_id_syst_up",   "#Sigma_phi_1_{t#bar{t}} mu_id_syst_up",          16,-3.2,3.2);
  Sigma_phi_1_mu_id_syst_down    = book<TH1F>("Sigma_phi_1_mu_id_syst_down",   "#Sigma_phi_1_{t#bar{t}} mu_id_syst_down",      16,-3.2,3.2);
  Sigma_phi_1_mu_iso_stat_up     = book<TH1F>("Sigma_phi_1_mu_iso_stat_up",   "#Sigma_phi_1_{t#bar{t}} mu_iso_stat_up",        16,-3.2,3.2);
  Sigma_phi_1_mu_iso_stat_down   = book<TH1F>("Sigma_phi_1_mu_iso_stat_down",   "#Sigma_phi_1_{t#bar{t}} mu_iso_stat_down",    16,-3.2,3.2);
  Sigma_phi_1_mu_iso_syst_up     = book<TH1F>("Sigma_phi_1_mu_iso_syst_up",   "#Sigma_phi_1_{t#bar{t}} mu_iso_syst_up",        16,-3.2,3.2);
  Sigma_phi_1_mu_iso_syst_down   = book<TH1F>("Sigma_phi_1_mu_iso_syst_down",   "#Sigma_phi_1_{t#bar{t}} mu_iso_syst_down",    16,-3.2,3.2);
  Sigma_phi_1_mu_trigger_stat_up     = book<TH1F>("Sigma_phi_1_mu_trigger_stat_up",   "#Sigma_phi_1_{t#bar{t}} mu_trigger_stat_up",        16,-3.2,3.2);
  Sigma_phi_1_mu_trigger_stat_down   = book<TH1F>("Sigma_phi_1_mu_trigger_stat_down",   "#Sigma_phi_1_{t#bar{t}} mu_trigger_stat_down",    16,-3.2,3.2);
  Sigma_phi_1_mu_trigger_syst_up     = book<TH1F>("Sigma_phi_1_mu_trigger_syst_up",   "#Sigma_phi_1_{t#bar{t}} mu_trigger_syst_up",        16,-3.2,3.2);
  Sigma_phi_1_mu_trigger_syst_down   = book<TH1F>("Sigma_phi_1_mu_trigger_syst_down",   "#Sigma_phi_1_{t#bar{t}} mu_trigger_syst_down",    16,-3.2,3.2);
  Sigma_phi_1_ele_id_up          = book<TH1F>("Sigma_phi_1_ele_id_up",   "#Sigma_phi_1_{t#bar{t}} ele_id_up",                  16,-3.2,3.2);
  Sigma_phi_1_ele_id_down        = book<TH1F>("Sigma_phi_1_ele_id_down", "#Sigma_phi_1_{t#bar{t}} ele_id_down",                16,-3.2,3.2);
  Sigma_phi_1_ele_trigger_up     = book<TH1F>("Sigma_phi_1_ele_trigger_up",   "#Sigma_phi_1_{t#bar{t}} ele_trigger_up",        16,-3.2,3.2);
  Sigma_phi_1_ele_trigger_down   = book<TH1F>("Sigma_phi_1_ele_trigger_down", "#Sigma_phi_1_{t#bar{t}} ele_trigger_down",      16,-3.2,3.2);
  Sigma_phi_1_ele_reco_up        = book<TH1F>("Sigma_phi_1_ele_reco_up",   "#Sigma_phi_1_{t#bar{t}} ele_reco_up",              16,-3.2,3.2);
  Sigma_phi_1_ele_reco_down      = book<TH1F>("Sigma_phi_1_ele_reco_down", "#Sigma_phi_1_{t#bar{t}} ele_reco_down",            16,-3.2,3.2);
  Sigma_phi_1_murmuf_upup        = book<TH1F>("Sigma_phi_1_murmuf_upup", "#Sigma_phi_1_{t#bar{t}} murmuf_upup",                16,-3.2,3.2);
  Sigma_phi_1_murmuf_upnone      = book<TH1F>("Sigma_phi_1_murmuf_upnone", "#Sigma_phi_1_{t#bar{t}} murmuf_upnone",            16,-3.2,3.2);
  Sigma_phi_1_murmuf_noneup      = book<TH1F>("Sigma_phi_1_murmuf_noneup", "#Sigma_phi_1_{t#bar{t}} murmuf_noneup",            16,-3.2,3.2);
  Sigma_phi_1_murmuf_nonedown    = book<TH1F>("Sigma_phi_1_murmuf_nonedown", "#Sigma_phi_1_{t#bar{t}} murmuf_nonedown",        16,-3.2,3.2);
  Sigma_phi_1_murmuf_downnone    = book<TH1F>("Sigma_phi_1_murmuf_downnone", "#Sigma_phi_1_{t#bar{t}} murmuf_downnone",        16,-3.2,3.2);
  Sigma_phi_1_murmuf_downdown    = book<TH1F>("Sigma_phi_1_murmuf_downdown", "#Sigma_phi_1_{t#bar{t}} murmuf_downdown",        16,-3.2,3.2);
  Sigma_phi_1_isr_up             = book<TH1F>("Sigma_phi_1_isr_up", "#Sigma_phi_1_{t#bar{t}} isr_up",                          16,-3.2,3.2);
  Sigma_phi_1_isr_down           = book<TH1F>("Sigma_phi_1_isr_down", "#Sigma_phi_1_{t#bar{t}} isr_down",                      16,-3.2,3.2);
  Sigma_phi_1_fsr_up             = book<TH1F>("Sigma_phi_1_fsr_up", "#Sigma_phi_1_{t#bar{t}} fsr_up",                          16,-3.2,3.2);
  Sigma_phi_1_fsr_down           = book<TH1F>("Sigma_phi_1_fsr_down", "#Sigma_phi_1_{t#bar{t}} fsr_down",                      16,-3.2,3.2);
  Sigma_phi_1_btag_cferr1_up     = book<TH1F>("Sigma_phi_1_btag_cferr1_up", "#Sigma_phi_1_{t#bar{t}} btag_cferr1_up",          16,-3.2,3.2);
  Sigma_phi_1_btag_cferr1_down   = book<TH1F>("Sigma_phi_1_btag_cferr1_down", "#Sigma_phi_1_{t#bar{t}} btag_cferr1_down",      16,-3.2,3.2);
  Sigma_phi_1_btag_cferr2_up     = book<TH1F>("Sigma_phi_1_btag_cferr2_up", "#Sigma_phi_1_{t#bar{t}} btag_cferr2_up",          16,-3.2,3.2);
  Sigma_phi_1_btag_cferr2_down   = book<TH1F>("Sigma_phi_1_btag_cferr2_down", "#Sigma_phi_1_{t#bar{t}} btag_cferr2_down",      16,-3.2,3.2);
  Sigma_phi_1_btag_hf_up         = book<TH1F>("Sigma_phi_1_btag_hf_up", "#Sigma_phi_1_{t#bar{t}} btag_hf_up",                  16,-3.2,3.2);
  Sigma_phi_1_btag_hf_down       = book<TH1F>("Sigma_phi_1_btag_hf_down", "#Sigma_phi_1_{t#bar{t}} btag_hf_down",              16,-3.2,3.2);
  Sigma_phi_1_btag_hfstats1_up   = book<TH1F>("Sigma_phi_1_btag_hfstats1_up", "#Sigma_phi_1_{t#bar{t}} btag_hfstats1_up",      16,-3.2,3.2);
  Sigma_phi_1_btag_hfstats1_down = book<TH1F>("Sigma_phi_1_btag_hfstats1_down", "#Sigma_phi_1_{t#bar{t}} btag_hfstats1_down",  16,-3.2,3.2);
  Sigma_phi_1_btag_hfstats2_up   = book<TH1F>("Sigma_phi_1_btag_hfstats2_up", "#Sigma_phi_1_{t#bar{t}} btag_hfstats2_up",      16,-3.2,3.2);
  Sigma_phi_1_btag_hfstats2_down = book<TH1F>("Sigma_phi_1_btag_hfstats2_down", "#Sigma_phi_1_{t#bar{t}} btag_hfstats2_down",  16,-3.2,3.2);
  Sigma_phi_1_btag_lf_up         = book<TH1F>("Sigma_phi_1_btag_lf_up", "#Sigma_phi_1_{t#bar{t}} btag_lf_up",                  16,-3.2,3.2);
  Sigma_phi_1_btag_lf_down       = book<TH1F>("Sigma_phi_1_btag_lf_down", "#Sigma_phi_1_{t#bar{t}} btag_lf_down",              16,-3.2,3.2);
  Sigma_phi_1_btag_lfstats1_up   = book<TH1F>("Sigma_phi_1_btag_lfstats1_up", "#Sigma_phi_1_{t#bar{t}} btag_lfstats1_up",      16,-3.2,3.2);
  Sigma_phi_1_btag_lfstats1_down = book<TH1F>("Sigma_phi_1_btag_lfstats1_down", "#Sigma_phi_1_{t#bar{t}} btag_lfstats1_down",  16,-3.2,3.2);
  Sigma_phi_1_btag_lfstats2_up   = book<TH1F>("Sigma_phi_1_btag_lfstats2_up", "#Sigma_phi_1_{t#bar{t}} btag_lfstats2_up",      16,-3.2,3.2);
  Sigma_phi_1_btag_lfstats2_down = book<TH1F>("Sigma_phi_1_btag_lfstats2_down", "#Sigma_phi_1_{t#bar{t}} btag_lfstats2_down",  16,-3.2,3.2);
  Sigma_phi_1_ttag_corr_up       = book<TH1F>("Sigma_phi_1_ttag_corr_up", "#Sigma_phi_1_{t#bar{t}} ttag_corr_up",              16,-3.2,3.2);
  Sigma_phi_1_ttag_corr_down     = book<TH1F>("Sigma_phi_1_ttag_corr_down", "#Sigma_phi_1_{t#bar{t}} ttag_corr_down",          16,-3.2,3.2);
  Sigma_phi_1_ttag_uncorr_up     = book<TH1F>("Sigma_phi_1_ttag_uncorr_up", "#Sigma_phi_1_{t#bar{t}} ttag_uncorr_up",          16,-3.2,3.2);
  Sigma_phi_1_ttag_uncorr_down   = book<TH1F>("Sigma_phi_1_ttag_uncorr_down", "#Sigma_phi_1_{t#bar{t}} ttag_counrr_down",      16,-3.2,3.2);
  Sigma_phi_1_tmistag_up         = book<TH1F>("Sigma_phi_1_tmistag_up", "#Sigma_phi_1_{t#bar{t}} [GeV] tmistag_up",            16,-3.2,3.2);
  Sigma_phi_1_tmistag_down       = book<TH1F>("Sigma_phi_1_tmistag_down", "#Sigma_phi_1_{t#bar{t}} [GeV] tmistag_down",        16,-3.2,3.2);
  Sigma_phi_1_toppt_a_up         = book<TH1F>("Sigma_phi_1_toppt_a_up", "#Sigma_phi_1_{t#bar{t}} [GeV] toppt_a_up",                16,-3.2,3.2);
  Sigma_phi_1_toppt_a_down       = book<TH1F>("Sigma_phi_1_toppt_a_down", "#Sigma_phi_1_{t#bar{t}} [GeV] toppt_a_down",            16,-3.2,3.2);
  Sigma_phi_1_toppt_b_up         = book<TH1F>("Sigma_phi_1_toppt_b_up", "#Sigma_phi_1_{t#bar{t}} [GeV] toppt_b_up",                16,-3.2,3.2);
  Sigma_phi_1_toppt_b_down       = book<TH1F>("Sigma_phi_1_toppt_b_down", "#Sigma_phi_1_{t#bar{t}} [GeV] toppt_b_down",            16,-3.2,3.2);
  // Sigma_phi_2                    = book<TH1F>("Sigma_phi_2",   "#Sigma_phi_2_{t#bar{t}} ",                                      16, -3.2, 3.2);
  Sigma_phi_2_mu_reco_up         = book<TH1F>("Sigma_phi_2_mu_reco_up",   "# Sigma_phi_2_{t#bar{t}} mu_reco_up",                16, -3.2, 3.2);
  Sigma_phi_2_mu_reco_down       = book<TH1F>("Sigma_phi_2_mu_reco_down", "# Sigma_phi_2_{t#bar{t}} mu_reco_down",              16,-3.2,3.2);
  Sigma_phi_2_pu_up              = book<TH1F>("Sigma_phi_2_pu_up",   "# Sigma_phi_2_{t#bar{t}} pu_up",                          16,-3.2,3.2);
  Sigma_phi_2_pu_down            = book<TH1F>("Sigma_phi_2_pu_down", "# Sigma_phi_2_{t#bar{t}} pu_down",                        16,-3.2,3.2);
  Sigma_phi_2_prefiring_up       = book<TH1F>("Sigma_phi_2_prefiring_up",   "# Sigma_phi_2_{t#bar{t}} prefiring_up",            16,-3.2,3.2);
  Sigma_phi_2_prefiring_down     = book<TH1F>("Sigma_phi_2_prefiring_down", "#Sigma_phi_2_{t#bar{t}} prefiring_down",          16,-3.2,3.2);
  Sigma_phi_2_mu_id_stat_up      = book<TH1F>("Sigma_phi_2_mu_id_stat_up",   "#Sigma_phi_2_{t#bar{t}} mu_id_stat_up",          16,-3.2,3.2);
  Sigma_phi_2_mu_id_stat_down    = book<TH1F>("Sigma_phi_2_mu_id_stat_down",   "#Sigma_phi_2_{t#bar{t}} mu_id_stat_down",      16,-3.2,3.2);
  Sigma_phi_2_mu_id_syst_up      = book<TH1F>("Sigma_phi_2_mu_id_syst_up",   "#Sigma_phi_2_{t#bar{t}} mu_id_syst_up",          16,-3.2,3.2);
  Sigma_phi_2_mu_id_syst_down    = book<TH1F>("Sigma_phi_2_mu_id_syst_down",   "#Sigma_phi_2_{t#bar{t}} mu_id_syst_down",      16,-3.2,3.2);
  Sigma_phi_2_mu_iso_stat_up     = book<TH1F>("Sigma_phi_2_mu_iso_stat_up",   "#Sigma_phi_2_{t#bar{t}} mu_iso_stat_up",        16,-3.2,3.2);
  Sigma_phi_2_mu_iso_stat_down   = book<TH1F>("Sigma_phi_2_mu_iso_stat_down",   "#Sigma_phi_2_{t#bar{t}} mu_iso_stat_down",    16,-3.2,3.2);
  Sigma_phi_2_mu_iso_syst_up     = book<TH1F>("Sigma_phi_2_mu_iso_syst_up",   "#Sigma_phi_2_{t#bar{t}} mu_iso_syst_up",        16,-3.2,3.2);
  Sigma_phi_2_mu_iso_syst_down   = book<TH1F>("Sigma_phi_2_mu_iso_syst_down",   "#Sigma_phi_2_{t#bar{t}} mu_iso_syst_down",    16,-3.2,3.2);
  Sigma_phi_2_mu_trigger_stat_up     = book<TH1F>("Sigma_phi_2_mu_trigger_stat_up",   "#Sigma_phi_2_{t#bar{t}} mu_trigger_stat_up",        16,-3.2,3.2);
  Sigma_phi_2_mu_trigger_stat_down   = book<TH1F>("Sigma_phi_2_mu_trigger_stat_down",   "#Sigma_phi_2_{t#bar{t}} mu_trigger_stat_down",    16,-3.2,3.2);
  Sigma_phi_2_mu_trigger_syst_up     = book<TH1F>("Sigma_phi_2_mu_trigger_syst_up",   "#Sigma_phi_2_{t#bar{t}} mu_trigger_syst_up",        16,-3.2,3.2);
  Sigma_phi_2_mu_trigger_syst_down   = book<TH1F>("Sigma_phi_2_mu_trigger_syst_down",   "#Sigma_phi_2_{t#bar{t}} mu_trigger_syst_down",    16,-3.2,3.2);
  Sigma_phi_2_ele_id_up          = book<TH1F>("Sigma_phi_2_ele_id_up",   "#Sigma_phi_2_{t#bar{t}} ele_id_up",                  16,-3.2,3.2);
  Sigma_phi_2_ele_id_down        = book<TH1F>("Sigma_phi_2_ele_id_down", "#Sigma_phi_2_{t#bar{t}} ele_id_down",                16,-3.2,3.2);
  Sigma_phi_2_ele_trigger_up     = book<TH1F>("Sigma_phi_2_ele_trigger_up",   "#Sigma_phi_2_{t#bar{t}} ele_trigger_up",        16,-3.2,3.2);
  Sigma_phi_2_ele_trigger_down   = book<TH1F>("Sigma_phi_2_ele_trigger_down", "#Sigma_phi_2_{t#bar{t}} ele_trigger_down",      16,-3.2,3.2);
  Sigma_phi_2_ele_reco_up        = book<TH1F>("Sigma_phi_2_ele_reco_up",   "#Sigma_phi_2_{t#bar{t}} ele_reco_up",              16,-3.2,3.2);
  Sigma_phi_2_ele_reco_down      = book<TH1F>("Sigma_phi_2_ele_reco_down", "#Sigma_phi_2_{t#bar{t}} ele_reco_down",            16,-3.2,3.2);
  Sigma_phi_2_murmuf_upup        = book<TH1F>("Sigma_phi_2_murmuf_upup", "#Sigma_phi_2_{t#bar{t}} murmuf_upup",                16,-3.2,3.2);
  Sigma_phi_2_murmuf_upnone      = book<TH1F>("Sigma_phi_2_murmuf_upnone", "#Sigma_phi_2_{t#bar{t}} murmuf_upnone",            16,-3.2,3.2);
  Sigma_phi_2_murmuf_noneup      = book<TH1F>("Sigma_phi_2_murmuf_noneup", "#Sigma_phi_2_{t#bar{t}} murmuf_noneup",            16,-3.2,3.2);
  Sigma_phi_2_murmuf_nonedown    = book<TH1F>("Sigma_phi_2_murmuf_nonedown", "#Sigma_phi_2_{t#bar{t}} murmuf_nonedown",        16,-3.2,3.2);
  Sigma_phi_2_murmuf_downnone    = book<TH1F>("Sigma_phi_2_murmuf_downnone", "#Sigma_phi_2_{t#bar{t}} murmuf_downnone",        16,-3.2,3.2);
  Sigma_phi_2_murmuf_downdown    = book<TH1F>("Sigma_phi_2_murmuf_downdown", "#Sigma_phi_2_{t#bar{t}} murmuf_downdown",        16,-3.2,3.2);
  Sigma_phi_2_isr_up             = book<TH1F>("Sigma_phi_2_isr_up", "#Sigma_phi_2_{t#bar{t}} isr_up",                          16,-3.2,3.2);
  Sigma_phi_2_isr_down           = book<TH1F>("Sigma_phi_2_isr_down", "#Sigma_phi_2_{t#bar{t}} isr_down",                      16,-3.2,3.2);
  Sigma_phi_2_fsr_up             = book<TH1F>("Sigma_phi_2_fsr_up", "#Sigma_phi_2_{t#bar{t}} fsr_up",                          16,-3.2,3.2);
  Sigma_phi_2_fsr_down           = book<TH1F>("Sigma_phi_2_fsr_down", "#Sigma_phi_2_{t#bar{t}} fsr_down",                      16,-3.2,3.2);
  Sigma_phi_2_btag_cferr1_up     = book<TH1F>("Sigma_phi_2_btag_cferr1_up", "#Sigma_phi_2_{t#bar{t}} btag_cferr1_up",          16,-3.2,3.2);
  Sigma_phi_2_btag_cferr1_down   = book<TH1F>("Sigma_phi_2_btag_cferr1_down", "#Sigma_phi_2_{t#bar{t}} btag_cferr1_down",      16,-3.2,3.2);
  Sigma_phi_2_btag_cferr2_up     = book<TH1F>("Sigma_phi_2_btag_cferr2_up", "#Sigma_phi_2_{t#bar{t}} btag_cferr2_up",          16,-3.2,3.2);
  Sigma_phi_2_btag_cferr2_down   = book<TH1F>("Sigma_phi_2_btag_cferr2_down", "#Sigma_phi_2_{t#bar{t}} btag_cferr2_down",      16,-3.2,3.2);
  Sigma_phi_2_btag_hf_up         = book<TH1F>("Sigma_phi_2_btag_hf_up", "#Sigma_phi_2_{t#bar{t}} btag_hf_up",                  16,-3.2,3.2);
  Sigma_phi_2_btag_hf_down       = book<TH1F>("Sigma_phi_2_btag_hf_down", "#Sigma_phi_2_{t#bar{t}} btag_hf_down",              16,-3.2,3.2);
  Sigma_phi_2_btag_hfstats1_up   = book<TH1F>("Sigma_phi_2_btag_hfstats1_up", "#Sigma_phi_2_{t#bar{t}} btag_hfstats1_up",      16,-3.2,3.2);
  Sigma_phi_2_btag_hfstats1_down = book<TH1F>("Sigma_phi_2_btag_hfstats1_down", "#Sigma_phi_2_{t#bar{t}} btag_hfstats1_down",  16,-3.2,3.2);
  Sigma_phi_2_btag_hfstats2_up   = book<TH1F>("Sigma_phi_2_btag_hfstats2_up", "#Sigma_phi_2_{t#bar{t}} btag_hfstats2_up",      16,-3.2,3.2);
  Sigma_phi_2_btag_hfstats2_down = book<TH1F>("Sigma_phi_2_btag_hfstats2_down", "#Sigma_phi_2_{t#bar{t}} btag_hfstats2_down",  16,-3.2,3.2);
  Sigma_phi_2_btag_lf_up         = book<TH1F>("Sigma_phi_2_btag_lf_up", "#Sigma_phi_2_{t#bar{t}} btag_lf_up",                  16,-3.2,3.2);
  Sigma_phi_2_btag_lf_down       = book<TH1F>("Sigma_phi_2_btag_lf_down", "#Sigma_phi_2_{t#bar{t}} btag_lf_down",              16,-3.2,3.2);
  Sigma_phi_2_btag_lfstats1_up   = book<TH1F>("Sigma_phi_2_btag_lfstats1_up", "#Sigma_phi_2_{t#bar{t}} btag_lfstats1_up",      16,-3.2,3.2);
  Sigma_phi_2_btag_lfstats1_down = book<TH1F>("Sigma_phi_2_btag_lfstats1_down", "#Sigma_phi_2_{t#bar{t}} btag_lfstats1_down",  16,-3.2,3.2);
  Sigma_phi_2_btag_lfstats2_up   = book<TH1F>("Sigma_phi_2_btag_lfstats2_up", "#Sigma_phi_2_{t#bar{t}} btag_lfstats2_up",      16,-3.2,3.2);
  Sigma_phi_2_btag_lfstats2_down = book<TH1F>("Sigma_phi_2_btag_lfstats2_down", "#Sigma_phi_2_{t#bar{t}} btag_lfstats2_down",  16,-3.2,3.2);
  Sigma_phi_2_ttag_corr_up       = book<TH1F>("Sigma_phi_2_ttag_corr_up", "#Sigma_phi_2_{t#bar{t}} ttag_corr_up",              16,-3.2,3.2);
  Sigma_phi_2_ttag_corr_down     = book<TH1F>("Sigma_phi_2_ttag_corr_down", "#Sigma_phi_2_{t#bar{t}} ttag_corr_down",          16,-3.2,3.2);
  Sigma_phi_2_ttag_uncorr_up     = book<TH1F>("Sigma_phi_2_ttag_uncorr_up", "#Sigma_phi_2_{t#bar{t}} ttag_uncorr_up",          16,-3.2,3.2);
  Sigma_phi_2_ttag_uncorr_down   = book<TH1F>("Sigma_phi_2_ttag_uncorr_down", "#Sigma_phi_2_{t#bar{t}} ttag_counrr_down",      16,-3.2,3.2);
  Sigma_phi_2_tmistag_up         = book<TH1F>("Sigma_phi_2_tmistag_up", "#Sigma_phi_2_{t#bar{t}} [GeV] tmistag_up",            16,-3.2,3.2);
  Sigma_phi_2_tmistag_down       = book<TH1F>("Sigma_phi_2_tmistag_down", "#Sigma_phi_2_{t#bar{t}} [GeV] tmistag_down",        16,-3.2,3.2);
  Sigma_phi_2_toppt_a_up         = book<TH1F>("Sigma_phi_2_toppt_a_up", "#Sigma_phi_2_{t#bar{t}} [GeV] toppt_a_up",                16,-3.2,3.2);
  Sigma_phi_2_toppt_a_down       = book<TH1F>("Sigma_phi_2_toppt_a_down", "#Sigma_phi_2_{t#bar{t}} [GeV] toppt_a_down",            16,-3.2,3.2);
  Sigma_phi_2_toppt_b_up         = book<TH1F>("Sigma_phi_2_toppt_b_up", "#Sigma_phi_2_{t#bar{t}} [GeV] toppt_b_up",                16,-3.2,3.2);
  Sigma_phi_2_toppt_b_down       = book<TH1F>("Sigma_phi_2_toppt_b_down", "#Sigma_phi_2_{t#bar{t}} [GeV] toppt_b_down",            16,-3.2,3.2);


  //template method start
  // template method variable: xi = tanh(DeltaY), 6 bins in [-1,1]
  // Nominals
  DeltaY_xi_reco_6          = book<TH1F>("DeltaY_xi_reco_6", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_12         = book<TH1F>("DeltaY_xi_reco_12", ";tanh(#Delta y)_{reco};Events / bin", 12, -1.0, 1.0);
  DeltaY_xi_reco_18         = book<TH1F>("DeltaY_xi_reco_18", ";tanh(#Delta y)_{reco};Events / bin", 18, -1.0, 1.0);
  DeltaY_xi_reco_24         = book<TH1F>("DeltaY_xi_reco_24", ";tanh(#Delta y)_{reco};Events / bin", 24, -1.0, 1.0);
  DeltaY_xi_reco_30         = book<TH1F>("DeltaY_xi_reco_30", ";tanh(#Delta y)_{reco};Events / bin", 30, -1.0, 1.0);
  DeltaY_xi_reco_36         = book<TH1F>("DeltaY_xi_reco_36", ";tanh(#Delta y)_{reco};Events / bin", 36, -1.0, 1.0);
  DeltaY_xi_reco_50         = book<TH1F>("DeltaY_xi_reco_50", ";tanh(#Delta y)_{reco};Events / bin", 50, -1.0, 1.0);
  // NoAC systematics (f=±1)
  DeltaY_xi_reco_6_NoAC_up   = book<TH1F>("DeltaY_xi_reco_6_NoAC_up", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_NoAC_down = book<TH1F>("DeltaY_xi_reco_6_NoAC_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_12_NoAC_up   = book<TH1F>("DeltaY_xi_reco_12_NoAC_up", ";tanh(#Delta y)_{reco};Events / bin", 12, -1.0, 1.0);
  DeltaY_xi_reco_12_NoAC_down = book<TH1F>("DeltaY_xi_reco_12_NoAC_down", ";tanh(#Delta y)_{reco};Events / bin", 12, -1.0, 1.0);
  DeltaY_xi_reco_18_NoAC_up   = book<TH1F>("DeltaY_xi_reco_18_NoAC_up", ";tanh(#Delta y)_{reco};Events / bin", 18, -1.0, 1.0);
  DeltaY_xi_reco_18_NoAC_down = book<TH1F>("DeltaY_xi_reco_18_NoAC_down", ";tanh(#Delta y)_{reco};Events / bin", 18, -1.0, 1.0);
  DeltaY_xi_reco_24_NoAC_up   = book<TH1F>("DeltaY_xi_reco_24_NoAC_up", ";tanh(#Delta y)_{reco};Events / bin", 24, -1.0, 1.0);
  DeltaY_xi_reco_24_NoAC_down = book<TH1F>("DeltaY_xi_reco_24_NoAC_down", ";tanh(#Delta y)_{reco};Events / bin", 24, -1.0, 1.0);
  DeltaY_xi_reco_30_NoAC_up   = book<TH1F>("DeltaY_xi_reco_30_NoAC_up", ";tanh(#Delta y)_{reco};Events / bin", 30, -1.0, 1.0);
  DeltaY_xi_reco_30_NoAC_down = book<TH1F>("DeltaY_xi_reco_30_NoAC_down", ";tanh(#Delta y)_{reco};Events / bin", 30, -1.0, 1.0);
  DeltaY_xi_reco_36_NoAC_up   = book<TH1F>("DeltaY_xi_reco_36_NoAC_up", ";tanh(#Delta y)_{reco};Events / bin", 36, -1.0, 1.0);
  DeltaY_xi_reco_36_NoAC_down = book<TH1F>("DeltaY_xi_reco_36_NoAC_down", ";tanh(#Delta y)_{reco};Events / bin", 36, -1.0, 1.0);
  DeltaY_xi_reco_50_NoAC_up   = book<TH1F>("DeltaY_xi_reco_50_NoAC_up", ";tanh(#Delta y)_{reco};Events / bin", 50, -1.0, 1.0);
  DeltaY_xi_reco_50_NoAC_down = book<TH1F>("DeltaY_xi_reco_50_NoAC_down", ";tanh(#Delta y)_{reco};Events / bin", 50, -1.0, 1.0);
  // NoAC systematics for f=±2, ±8, ±12
  DeltaY_xi_reco_6_NoAC_up_f2   = book<TH1F>("DeltaY_xi_reco_6_NoAC_up_f2", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_NoAC_down_f2 = book<TH1F>("DeltaY_xi_reco_6_NoAC_down_f2", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_NoAC_up_f8   = book<TH1F>("DeltaY_xi_reco_6_NoAC_up_f8", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_NoAC_down_f8 = book<TH1F>("DeltaY_xi_reco_6_NoAC_down_f8", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_NoAC_up_f12   = book<TH1F>("DeltaY_xi_reco_6_NoAC_up_f12", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_NoAC_down_f12 = book<TH1F>("DeltaY_xi_reco_6_NoAC_down_f12", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_12_NoAC_up_f2   = book<TH1F>("DeltaY_xi_reco_12_NoAC_up_f2", ";tanh(#Delta y)_{reco};Events / bin", 12, -1.0, 1.0);
  DeltaY_xi_reco_12_NoAC_down_f2 = book<TH1F>("DeltaY_xi_reco_12_NoAC_down_f2", ";tanh(#Delta y)_{reco};Events / bin", 12, -1.0, 1.0);
  DeltaY_xi_reco_12_NoAC_up_f8   = book<TH1F>("DeltaY_xi_reco_12_NoAC_up_f8", ";tanh(#Delta y)_{reco};Events / bin", 12, -1.0, 1.0);
  DeltaY_xi_reco_12_NoAC_down_f8 = book<TH1F>("DeltaY_xi_reco_12_NoAC_down_f8", ";tanh(#Delta y)_{reco};Events / bin", 12, -1.0, 1.0);
  DeltaY_xi_reco_12_NoAC_up_f12   = book<TH1F>("DeltaY_xi_reco_12_NoAC_up_f12", ";tanh(#Delta y)_{reco};Events / bin", 12, -1.0, 1.0);
  DeltaY_xi_reco_12_NoAC_down_f12 = book<TH1F>("DeltaY_xi_reco_12_NoAC_down_f12", ";tanh(#Delta y)_{reco};Events / bin", 12, -1.0, 1.0);
  DeltaY_xi_reco_18_NoAC_up_f2   = book<TH1F>("DeltaY_xi_reco_18_NoAC_up_f2", ";tanh(#Delta y)_{reco};Events / bin", 18, -1.0, 1.0);
  DeltaY_xi_reco_18_NoAC_down_f2 = book<TH1F>("DeltaY_xi_reco_18_NoAC_down_f2", ";tanh(#Delta y)_{reco};Events / bin", 18, -1.0, 1.0);
  DeltaY_xi_reco_18_NoAC_up_f8   = book<TH1F>("DeltaY_xi_reco_18_NoAC_up_f8", ";tanh(#Delta y)_{reco};Events / bin", 18, -1.0, 1.0);
  DeltaY_xi_reco_18_NoAC_down_f8 = book<TH1F>("DeltaY_xi_reco_18_NoAC_down_f8", ";tanh(#Delta y)_{reco};Events / bin", 18, -1.0, 1.0);
  DeltaY_xi_reco_18_NoAC_up_f12   = book<TH1F>("DeltaY_xi_reco_18_NoAC_up_f12", ";tanh(#Delta y)_{reco};Events / bin", 18, -1.0, 1.0);
  DeltaY_xi_reco_18_NoAC_down_f12 = book<TH1F>("DeltaY_xi_reco_18_NoAC_down_f12", ";tanh(#Delta y)_{reco};Events / bin", 18, -1.0, 1.0);
  DeltaY_xi_reco_24_NoAC_up_f2   = book<TH1F>("DeltaY_xi_reco_24_NoAC_up_f2", ";tanh(#Delta y)_{reco};Events / bin", 24, -1.0, 1.0);
  DeltaY_xi_reco_24_NoAC_down_f2 = book<TH1F>("DeltaY_xi_reco_24_NoAC_down_f2", ";tanh(#Delta y)_{reco};Events / bin", 24, -1.0, 1.0);
  DeltaY_xi_reco_24_NoAC_up_f8   = book<TH1F>("DeltaY_xi_reco_24_NoAC_up_f8", ";tanh(#Delta y)_{reco};Events / bin", 24, -1.0, 1.0);
  DeltaY_xi_reco_24_NoAC_down_f8 = book<TH1F>("DeltaY_xi_reco_24_NoAC_down_f8", ";tanh(#Delta y)_{reco};Events / bin", 24, -1.0, 1.0);
  DeltaY_xi_reco_24_NoAC_up_f12   = book<TH1F>("DeltaY_xi_reco_24_NoAC_up_f12", ";tanh(#Delta y)_{reco};Events / bin", 24, -1.0, 1.0);
  DeltaY_xi_reco_24_NoAC_down_f12 = book<TH1F>("DeltaY_xi_reco_24_NoAC_down_f12", ";tanh(#Delta y)_{reco};Events / bin", 24, -1.0, 1.0);
  DeltaY_xi_reco_30_NoAC_up_f2   = book<TH1F>("DeltaY_xi_reco_30_NoAC_up_f2", ";tanh(#Delta y)_{reco};Events / bin", 30, -1.0, 1.0);
  DeltaY_xi_reco_30_NoAC_down_f2 = book<TH1F>("DeltaY_xi_reco_30_NoAC_down_f2", ";tanh(#Delta y)_{reco};Events / bin", 30, -1.0, 1.0);
  DeltaY_xi_reco_30_NoAC_up_f8   = book<TH1F>("DeltaY_xi_reco_30_NoAC_up_f8", ";tanh(#Delta y)_{reco};Events / bin", 30, -1.0, 1.0);
  DeltaY_xi_reco_30_NoAC_down_f8 = book<TH1F>("DeltaY_xi_reco_30_NoAC_down_f8", ";tanh(#Delta y)_{reco};Events / bin", 30, -1.0, 1.0);
  DeltaY_xi_reco_30_NoAC_up_f12   = book<TH1F>("DeltaY_xi_reco_30_NoAC_up_f12", ";tanh(#Delta y)_{reco};Events / bin", 30, -1.0, 1.0);
  DeltaY_xi_reco_30_NoAC_down_f12 = book<TH1F>("DeltaY_xi_reco_30_NoAC_down_f12", ";tanh(#Delta y)_{reco};Events / bin", 30, -1.0, 1.0);
  DeltaY_xi_reco_36_NoAC_up_f2   = book<TH1F>("DeltaY_xi_reco_36_NoAC_up_f2", ";tanh(#Delta y)_{reco};Events / bin", 36, -1.0, 1.0);
  DeltaY_xi_reco_36_NoAC_down_f2 = book<TH1F>("DeltaY_xi_reco_36_NoAC_down_f2", ";tanh(#Delta y)_{reco};Events / bin", 36, -1.0, 1.0);
  DeltaY_xi_reco_36_NoAC_up_f8   = book<TH1F>("DeltaY_xi_reco_36_NoAC_up_f8", ";tanh(#Delta y)_{reco};Events / bin", 36, -1.0, 1.0);
  DeltaY_xi_reco_36_NoAC_down_f8 = book<TH1F>("DeltaY_xi_reco_36_NoAC_down_f8", ";tanh(#Delta y)_{reco};Events / bin", 36, -1.0, 1.0);
  DeltaY_xi_reco_36_NoAC_up_f12   = book<TH1F>("DeltaY_xi_reco_36_NoAC_up_f12", ";tanh(#Delta y)_{reco};Events / bin", 36, -1.0, 1.0);
  DeltaY_xi_reco_36_NoAC_down_f12 = book<TH1F>("DeltaY_xi_reco_36_NoAC_down_f12", ";tanh(#Delta y)_{reco};Events / bin", 36, -1.0, 1.0);
  DeltaY_xi_reco_50_NoAC_up_f2   = book<TH1F>("DeltaY_xi_reco_50_NoAC_up_f2", ";tanh(#Delta y)_{reco};Events / bin", 50, -1.0, 1.0);
  DeltaY_xi_reco_50_NoAC_down_f2 = book<TH1F>("DeltaY_xi_reco_50_NoAC_down_f2", ";tanh(#Delta y)_{reco};Events / bin", 50, -1.0, 1.0);
  DeltaY_xi_reco_50_NoAC_up_f8   = book<TH1F>("DeltaY_xi_reco_50_NoAC_up_f8", ";tanh(#Delta y)_{reco};Events / bin", 50, -1.0, 1.0);
  DeltaY_xi_reco_50_NoAC_down_f8 = book<TH1F>("DeltaY_xi_reco_50_NoAC_down_f8", ";tanh(#Delta y)_{reco};Events / bin", 50, -1.0, 1.0);
  DeltaY_xi_reco_50_NoAC_up_f12   = book<TH1F>("DeltaY_xi_reco_50_NoAC_up_f12", ";tanh(#Delta y)_{reco};Events / bin", 50, -1.0, 1.0);
  DeltaY_xi_reco_50_NoAC_down_f12 = book<TH1F>("DeltaY_xi_reco_50_NoAC_down_f12", ";tanh(#Delta y)_{reco};Events / bin", 50, -1.0, 1.0);
  // lepton/trigger/pileup/prefiring
  DeltaY_xi_reco_6_ele_reco_up      = book<TH1F>("DeltaY_xi_reco_6_ele_reco_up",      ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_ele_reco_down    = book<TH1F>("DeltaY_xi_reco_6_ele_reco_down",    ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_ele_id_up        = book<TH1F>("DeltaY_xi_reco_6_ele_id_up",        ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_ele_id_down      = book<TH1F>("DeltaY_xi_reco_6_ele_id_down",      ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_ele_trigger_up   = book<TH1F>("DeltaY_xi_reco_6_ele_trigger_up",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_ele_trigger_down = book<TH1F>("DeltaY_xi_reco_6_ele_trigger_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_reco_up       = book<TH1F>("DeltaY_xi_reco_6_mu_reco_up",       ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_reco_down     = book<TH1F>("DeltaY_xi_reco_6_mu_reco_down",     ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_iso_stat_up   = book<TH1F>("DeltaY_xi_reco_6_mu_iso_stat_up",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_iso_stat_down = book<TH1F>("DeltaY_xi_reco_6_mu_iso_stat_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_iso_syst_up   = book<TH1F>("DeltaY_xi_reco_6_mu_iso_syst_up",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_iso_syst_down = book<TH1F>("DeltaY_xi_reco_6_mu_iso_syst_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_id_stat_up    = book<TH1F>("DeltaY_xi_reco_6_mu_id_stat_up",    ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_id_stat_down  = book<TH1F>("DeltaY_xi_reco_6_mu_id_stat_down",  ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_id_syst_up    = book<TH1F>("DeltaY_xi_reco_6_mu_id_syst_up",    ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_id_syst_down  = book<TH1F>("DeltaY_xi_reco_6_mu_id_syst_down",  ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_trigger_stat_up   = book<TH1F>("DeltaY_xi_reco_6_mu_trigger_stat_up",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_trigger_stat_down = book<TH1F>("DeltaY_xi_reco_6_mu_trigger_stat_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_trigger_syst_up   = book<TH1F>("DeltaY_xi_reco_6_mu_trigger_syst_up",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_mu_trigger_syst_down = book<TH1F>("DeltaY_xi_reco_6_mu_trigger_syst_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_pu_up              = book<TH1F>("DeltaY_xi_reco_6_pu_up",              ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_pu_down            = book<TH1F>("DeltaY_xi_reco_6_pu_down",            ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_prefiring_up       = book<TH1F>("DeltaY_xi_reco_6_prefiring_up",       ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_prefiring_down     = book<TH1F>("DeltaY_xi_reco_6_prefiring_down",     ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  // scales, isr fsr
  DeltaY_xi_reco_6_murmuf_upup        = book<TH1F>("DeltaY_xi_reco_6_murmuf_upup",        ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_murmuf_upnone      = book<TH1F>("DeltaY_xi_reco_6_murmuf_upnone",      ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_murmuf_noneup      = book<TH1F>("DeltaY_xi_reco_6_murmuf_noneup",      ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_murmuf_nonedown    = book<TH1F>("DeltaY_xi_reco_6_murmuf_nonedown",    ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_murmuf_downnone    = book<TH1F>("DeltaY_xi_reco_6_murmuf_downnone",    ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_murmuf_downdown    = book<TH1F>("DeltaY_xi_reco_6_murmuf_downdown",    ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_isr_up             = book<TH1F>("DeltaY_xi_reco_6_isr_up",             ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_isr_down           = book<TH1F>("DeltaY_xi_reco_6_isr_down",           ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_fsr_up             = book<TH1F>("DeltaY_xi_reco_6_fsr_up",             ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_fsr_down           = book<TH1F>("DeltaY_xi_reco_6_fsr_down",           ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  // btag
  DeltaY_xi_reco_6_btag_cferr1_up     = book<TH1F>("DeltaY_xi_reco_6_btag_cferr1_up",     ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_cferr1_down   = book<TH1F>("DeltaY_xi_reco_6_btag_cferr1_down",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_cferr2_up     = book<TH1F>("DeltaY_xi_reco_6_btag_cferr2_up",     ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_cferr2_down   = book<TH1F>("DeltaY_xi_reco_6_btag_cferr2_down",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_hf_up         = book<TH1F>("DeltaY_xi_reco_6_btag_hf_up",         ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_hf_down       = book<TH1F>("DeltaY_xi_reco_6_btag_hf_down",       ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_hfstats1_up   = book<TH1F>("DeltaY_xi_reco_6_btag_hfstats1_up",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_hfstats1_down = book<TH1F>("DeltaY_xi_reco_6_btag_hfstats1_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_hfstats2_up   = book<TH1F>("DeltaY_xi_reco_6_btag_hfstats2_up",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_hfstats2_down = book<TH1F>("DeltaY_xi_reco_6_btag_hfstats2_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_lf_up         = book<TH1F>("DeltaY_xi_reco_6_btag_lf_up",         ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_lf_down       = book<TH1F>("DeltaY_xi_reco_6_btag_lf_down",       ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_lfstats1_up   = book<TH1F>("DeltaY_xi_reco_6_btag_lfstats1_up",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_lfstats1_down = book<TH1F>("DeltaY_xi_reco_6_btag_lfstats1_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_lfstats2_up   = book<TH1F>("DeltaY_xi_reco_6_btag_lfstats2_up",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_btag_lfstats2_down = book<TH1F>("DeltaY_xi_reco_6_btag_lfstats2_down", ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  // ttag, mistag, toppt
  DeltaY_xi_reco_6_ttag_corr_up       = book<TH1F>("DeltaY_xi_reco_6_ttag_corr_up",       ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_ttag_corr_down     = book<TH1F>("DeltaY_xi_reco_6_ttag_corr_down",     ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_ttag_uncorr_up     = book<TH1F>("DeltaY_xi_reco_6_ttag_uncorr_up",     ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_ttag_uncorr_down   = book<TH1F>("DeltaY_xi_reco_6_ttag_uncorr_down",   ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_tmistag_up         = book<TH1F>("DeltaY_xi_reco_6_tmistag_up",         ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_tmistag_down       = book<TH1F>("DeltaY_xi_reco_6_tmistag_down",       ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_toppt_a_up         = book<TH1F>("DeltaY_xi_reco_6_toppt_a_up",         ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_toppt_a_down       = book<TH1F>("DeltaY_xi_reco_6_toppt_a_down",       ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_toppt_b_up         = book<TH1F>("DeltaY_xi_reco_6_toppt_b_up",         ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);
  DeltaY_xi_reco_6_toppt_b_down       = book<TH1F>("DeltaY_xi_reco_6_toppt_b_down",       ";tanh(#Delta y)_{reco};Events / bin", 6, -1.0, 1.0);

  // --- NoAC templates removed (now using dedicated NoAC_up/down histograms) ---
  // --- BOOK xi systematics for 12 and 50 bins into name->hist map ---
  {
    auto book_by_name = [&](const std::string &hname, int nb){
      if(h_deltaY_xi_reco_map.find(hname) == h_deltaY_xi_reco_map.end()){
        h_deltaY_xi_reco_map[hname] = book<TH1F>(hname.c_str(), ";tanh(#Delta y)_{reco};Events / bin", nb, -1.0, 1.0);
      }
    };
    // Helper function to get exact suffix matching kNoACSpecs
    auto get_noac_suffix = [](float fv) -> std::string {
      static const std::map<float, std::string> f_to_suffix = {
        {-100.0f, "noacm100"}, {-12.0f, "noacm12"}, {-8.0f, "noacm8"}, {-4.0f, "noacm4"}, {-2.0f, "noacm2"},
        {-1.0f, "noacm1"}, {-0.8f, "noacm08"}, {-0.6f, "noacm06"}, {-0.4f, "noacm04"}, {-0.2f, "noacm02"},
        {0.0f, "noac0"},
        {0.2f, "noac02"}, {0.4f, "noac04"}, {0.6f, "noac06"}, {0.8f, "noac08"},
        {1.0f, "noac1"}, {2.0f, "noac2"}, {4.0f, "noac4"}, {8.0f, "noac8"}, {12.0f, "noac12"}, {100.0f, "noac100"}
      };
      auto it = f_to_suffix.find(fv);
      if(it != f_to_suffix.end()) return it->second;
      // Fallback (shouldn't happen)
      if(std::abs(fv) < 0.01f) return "noac0";
      else if(fv < 0) return "noacm" + std::to_string(static_cast<int>(std::abs(fv)));
      else return "noac" + std::to_string(static_cast<int>(fv));
    };
    // simple up/down tags
    const std::vector<std::pair<std::string,int>> xi_binnings = {
      {"DeltaY_xi_reco_12_", 12},
      {"DeltaY_xi_reco_18_", 18},
      {"DeltaY_xi_reco_24_", 24},
      {"DeltaY_xi_reco_30_", 30},
      {"DeltaY_xi_reco_36_", 36},
      {"DeltaY_xi_reco_50_", 50}
    };
    // Book all NoAC f-value histograms (matching Hists naming: DeltaY_xi_reco_N_noacX)
    const int all_bins[] = {6, 12, 18, 24, 30, 36, 50};
    const std::vector<float> all_f_values = {-100.0f, -12.0f, -8.0f, -4.0f, -2.0f, -1.0f, -0.8f, -0.6f, -0.4f, -0.2f, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 2.0f, 4.0f, 8.0f, 12.0f, 100.0f};
    for(int nb : all_bins){
      for(float fv : all_f_values){
        std::string suffix = get_noac_suffix(fv);
        std::string hname = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + suffix;
        book_by_name(hname, nb);
      }
    }
    std::vector<std::string> simple_tags = {"ele_reco","ele_id","ele_trigger","mu_reco","mu_iso_stat","mu_iso_syst","mu_id_stat","mu_id_syst","mu_trigger_stat","mu_trigger_syst","pu","prefiring"};
    for(const auto &tg : simple_tags){
      for(const auto &cfg : xi_binnings){
        book_by_name(cfg.first + tg + "_up", cfg.second);
        book_by_name(cfg.first + tg + "_down", cfg.second);
      }
    }
    // murmuf 6-point
    std::vector<std::string> mm = {"upup","upnone","noneup","nonedown","downnone","downdown"};
    for(const auto &v : mm){
      for(const auto &cfg : xi_binnings){
        book_by_name(cfg.first + std::string("murmuf_") + v, cfg.second);
      }
    }
    // ps (isr/fsr)
    std::vector<std::string> ps = {"isr_up","isr_down","fsr_up","fsr_down"};
    for(const auto &v : ps){
      for(const auto &cfg : xi_binnings){
        book_by_name(cfg.first + v, cfg.second);
      }
    }
    // btag set
    std::vector<std::string> btag_tags = {
      "btag_cferr1_up","btag_cferr1_down","btag_cferr2_up","btag_cferr2_down",
      "btag_hf_up","btag_hf_down","btag_hfstats1_up","btag_hfstats1_down",
      "btag_hfstats2_up","btag_hfstats2_down","btag_lf_up","btag_lf_down",
      "btag_lfstats1_up","btag_lfstats1_down","btag_lfstats2_up","btag_lfstats2_down"
    };
    for(const auto &v : btag_tags){
      for(const auto &cfg : xi_binnings){
        book_by_name(cfg.first + v, cfg.second);
      }
    }
    // ttag
    std::vector<std::string> ttag_tags = {"ttag_corr_up","ttag_corr_down","ttag_uncorr_up","ttag_uncorr_down"};
    for(const auto &v : ttag_tags){
      for(const auto &cfg : xi_binnings){
        book_by_name(cfg.first + v, cfg.second);
      }
    }
    // tmistag
    for(const auto &v : {std::string("tmistag_up"), std::string("tmistag_down")}){
      for(const auto &cfg : xi_binnings){
        book_by_name(cfg.first + v, cfg.second);
      }
    }
    // toppt
    std::vector<std::string> toppt_tags = {"toppt_a_up","toppt_a_down","toppt_b_up","toppt_b_down"};
    for(const auto &v : toppt_tags){
      for(const auto &cfg : xi_binnings){
        book_by_name(cfg.first + v, cfg.second);
      }
    }
  }
  //template method end
}



void ZprimeSemiLeptonicSystematicsHists::fill(const Event & event){

  double weight = event.weight;
  // electron systematics: reco, id, trigger
  float ele_reco_nominal   = event.get(h_ele_reco);
  float ele_reco_up        = event.get(h_ele_reco_up);
  float ele_reco_down      = event.get(h_ele_reco_down);
  float ele_id_nominal     = event.get(h_ele_id);
  float ele_id_up          = event.get(h_ele_id_up);
  float ele_id_down        = event.get(h_ele_id_down);
  float ele_trigger_nominal= event.get(h_ele_trigger);
  float ele_trigger_up     = event.get(h_ele_trigger_up);
  float ele_trigger_down   = event.get(h_ele_trigger_down);
  // muon systematics: reco, id_stat, id_syst, trigger_stat, trigger_syst, iso_stat, iso_syst
  float mu_reco_nominal    = event.get(h_mu_reco);
  float mu_reco_up         = event.get(h_mu_reco_up);
  float mu_reco_down       = event.get(h_mu_reco_down);
  float mu_id_stat_nominal     = event.get(h_mu_id_stat);
  float mu_id_stat_up     = event.get(h_mu_id_stat_up);
  float mu_id_stat_down   = event.get(h_mu_id_stat_down);
  float mu_id_syst_nominal     = event.get(h_mu_id_syst);
  float mu_id_syst_up     = event.get(h_mu_id_syst_up);
  float mu_id_syst_down   = event.get(h_mu_id_syst_down);
  float mu_trigger_stat_nominal     = event.get(h_mu_trigger_stat);
  float mu_trigger_stat_up     = event.get(h_mu_trigger_stat_up);
  float mu_trigger_stat_down   = event.get(h_mu_trigger_stat_down);
  float mu_trigger_syst_nominal     = event.get(h_mu_trigger_syst);
  float mu_trigger_syst_up     = event.get(h_mu_trigger_syst_up);
  float mu_trigger_syst_down   = event.get(h_mu_trigger_syst_down);
  float mu_iso_stat_nominal     = event.get(h_mu_iso_stat);
  float mu_iso_stat_up     = event.get(h_mu_iso_stat_up);
  float mu_iso_stat_down   = event.get(h_mu_iso_stat_down);
  float mu_iso_syst_nominal     = event.get(h_mu_iso_syst);
  float mu_iso_syst_up     = event.get(h_mu_iso_syst_up);
  float mu_iso_syst_down   = event.get(h_mu_iso_syst_down);
  // Pileup reweighting systematics
  float pu_nominal         = event.get(h_pu);
  float pu_up              = event.get(h_pu_up);
  float pu_down            = event.get(h_pu_down);
  // Prefiring systematics
  float prefiring_nominal  = event.get(h_prefiring);
  float prefiring_up       = event.get(h_prefiring_up);
  float prefiring_down     = event.get(h_prefiring_down);
  // b-tagging systematics: central, cferr1, cferr2, hf, hfstats1, hfstats2, lf, lfstats1, lfstats2
  float btag_nominal       = event.get(h_btag);
  float btag_cferr1_up     = event.get(h_btag_cferr1_up);
  float btag_cferr1_down   = event.get(h_btag_cferr1_down);
  float btag_cferr2_up     = event.get(h_btag_cferr2_up);
  float btag_cferr2_down   = event.get(h_btag_cferr2_down);
  float btag_hf_up         = event.get(h_btag_hf_up);
  float btag_hf_down       = event.get(h_btag_hf_down);
  float btag_hfstats1_up   = event.get(h_btag_hfstats1_up);
  float btag_hfstats1_down = event.get(h_btag_hfstats1_down);
  float btag_hfstats2_up   = event.get(h_btag_hfstats2_up);
  float btag_hfstats2_down = event.get(h_btag_hfstats2_down);
  float btag_lf_up         = event.get(h_btag_lf_up);
  float btag_lf_down       = event.get(h_btag_lf_down);
  float btag_lfstats1_up   = event.get(h_btag_lfstats1_up);
  float btag_lfstats1_down = event.get(h_btag_lfstats1_down);
  float btag_lfstats2_up   = event.get(h_btag_lfstats2_up);
  float btag_lfstats2_down = event.get(h_btag_lfstats2_down);
  // Top tagging systematics
  float ttag_nominal       = event.get(h_ttag);
  float ttag_corr_up       = event.get(h_ttag_corr_up);
  float ttag_corr_down     = event.get(h_ttag_corr_down);
  float ttag_uncorr_up     = event.get(h_ttag_uncorr_up);
  float ttag_uncorr_down   = event.get(h_ttag_uncorr_down);
  // Top mistagging systematics
  float tmistag_nominal    = event.get(h_tmistag);
  float tmistag_up         = event.get(h_tmistag_up);
  float tmistag_down       = event.get(h_tmistag_down);
  // Top pT reweighting systematics
  float toppt_a_up         = event.get(h_toppt_a_up);
  float toppt_a_down       = event.get(h_toppt_a_down);
  float toppt_b_up         = event.get(h_toppt_b_up);
  float toppt_b_down       = event.get(h_toppt_b_down);
  // muR and muF at ME level systematics
  float murmuf_upup        = event.get(h_murmuf_upup);
  float murmuf_upnone      = event.get(h_murmuf_upnone);
  float murmuf_noneup      = event.get(h_murmuf_noneup);
  float murmuf_nonedown    = event.get(h_murmuf_nonedown);
  float murmuf_downnone    = event.get(h_murmuf_downnone);
  float murmuf_downdown    = event.get(h_murmuf_downdown);
  // ISR and FSR systematics
  float isr_up             = event.get(h_isr_up);
  float isr_down           = event.get(h_isr_down);
  float fsr_up             = event.get(h_fsr_up);
  float fsr_down           = event.get(h_fsr_down);

  // up/down variations: electron, muon, pu, prefiring
  vector<string> names       = {"ele_reco", "ele_id", "ele_trigger", "mu_reco", "mu_iso_stat","mu_iso_syst", "mu_id_stat","mu_id_syst","mu_trigger_stat","mu_trigger_syst", "pu", "prefiring"};
  vector<float> syst_nominal = {ele_reco_nominal, ele_id_nominal, ele_trigger_nominal, mu_reco_nominal, mu_iso_stat_nominal, mu_iso_syst_nominal, mu_id_stat_nominal, mu_id_syst_nominal, mu_trigger_stat_nominal, mu_trigger_syst_nominal, pu_nominal, prefiring_nominal};
  vector<float> syst_up      = {ele_reco_up, ele_id_up, ele_trigger_up, mu_reco_up, mu_iso_stat_up, mu_iso_syst_up, mu_id_stat_up, mu_id_syst_up, mu_trigger_stat_up, mu_trigger_syst_up, pu_up, prefiring_up};
  vector<float> syst_down    = {ele_reco_down, ele_id_down, ele_trigger_down, mu_reco_down, mu_iso_stat_down, mu_iso_syst_down, mu_id_stat_down, mu_id_syst_down, mu_trigger_stat_down, mu_trigger_syst_down, pu_down, prefiring_down};
  // up variation histograms
  vector<TH1F*> hists_up     = {DeltaY_ele_reco_up, DeltaY_ele_id_up, DeltaY_ele_trigger_up, DeltaY_mu_reco_up, DeltaY_mu_iso_stat_up,DeltaY_mu_iso_syst_up, DeltaY_mu_id_stat_up, DeltaY_mu_id_syst_up, DeltaY_mu_trigger_stat_up, DeltaY_mu_trigger_syst_up, DeltaY_pu_up, DeltaY_prefiring_up};
  vector<TH1F*> hists_dy_d1_up = {DeltaY_reco_d1_ele_reco_up, DeltaY_reco_d1_ele_id_up, DeltaY_reco_d1_ele_trigger_up, DeltaY_reco_d1_mu_reco_up, DeltaY_reco_d1_mu_iso_stat_up,DeltaY_reco_d1_mu_iso_syst_up, DeltaY_reco_d1_mu_id_stat_up, DeltaY_reco_d1_mu_id_syst_up, DeltaY_reco_d1_mu_trigger_stat_up, DeltaY_reco_d1_mu_trigger_syst_up, DeltaY_reco_d1_pu_up, DeltaY_reco_d1_prefiring_up};
  vector<TH1F*> hists_dy_d2_up = {DeltaY_reco_d2_ele_reco_up, DeltaY_reco_d2_ele_id_up, DeltaY_reco_d2_ele_trigger_up, DeltaY_reco_d2_mu_reco_up, DeltaY_reco_d2_mu_iso_stat_up,DeltaY_reco_d2_mu_iso_syst_up, DeltaY_reco_d2_mu_id_stat_up, DeltaY_reco_d2_mu_id_syst_up, DeltaY_reco_d2_mu_trigger_stat_up, DeltaY_reco_d2_mu_trigger_syst_up, DeltaY_reco_d2_pu_up, DeltaY_reco_d2_prefiring_up};
  vector<TH1F*> hists_sigma_1_up = {Sigma_phi_1_ele_reco_up, Sigma_phi_1_ele_id_up, Sigma_phi_1_ele_trigger_up, Sigma_phi_1_mu_reco_up, Sigma_phi_1_mu_iso_stat_up,Sigma_phi_1_mu_iso_syst_up, Sigma_phi_1_mu_id_stat_up, Sigma_phi_1_mu_id_syst_up, Sigma_phi_1_mu_trigger_stat_up, Sigma_phi_1_mu_trigger_syst_up, Sigma_phi_1_pu_up, Sigma_phi_1_prefiring_up};
  vector<TH1F*> hists_sigma_2_up = {Sigma_phi_2_ele_reco_up, Sigma_phi_2_ele_id_up, Sigma_phi_2_ele_trigger_up, Sigma_phi_2_mu_reco_up, Sigma_phi_2_mu_iso_stat_up,Sigma_phi_2_mu_iso_syst_up, Sigma_phi_2_mu_id_stat_up, Sigma_phi_2_mu_id_syst_up, Sigma_phi_2_mu_trigger_stat_up, Sigma_phi_2_mu_trigger_syst_up, Sigma_phi_2_pu_up, Sigma_phi_2_prefiring_up};
  vector<TH1F*> hists_deltaY_xi_reco_6_up = {DeltaY_xi_reco_6_ele_reco_up, DeltaY_xi_reco_6_ele_id_up, DeltaY_xi_reco_6_ele_trigger_up, DeltaY_xi_reco_6_mu_reco_up, DeltaY_xi_reco_6_mu_iso_stat_up,DeltaY_xi_reco_6_mu_iso_syst_up, DeltaY_xi_reco_6_mu_id_stat_up, DeltaY_xi_reco_6_mu_id_syst_up, DeltaY_xi_reco_6_mu_trigger_stat_up, DeltaY_xi_reco_6_mu_trigger_syst_up, DeltaY_xi_reco_6_pu_up, DeltaY_xi_reco_6_prefiring_up};
  // cos_theta1k_antiLep
  // cos_theta1r_antiLep
  // cos_theta1n_antiLep
  // cos_theta1kStar_antiLep
  // cos_theta1rStar_antiLep
  // cos_theta2k_Lep
  // cos_theta2r_Lep
  // cos_theta2n_Lep
  // cos_theta2kStar_Lep
  // cos_theta2rStar_Lep
  // cos_theta1k
  // cos_theta1r
  // cos_theta1n
  // cos_theta1kStar
  // cos_theta1rStar
  // cos_theta2k
  // cos_theta2r
  // cos_theta2n
  // cos_theta2kStar
  // cos_theta2rStar
  // Ckk
  // Crr
  // Cnn
  // Crk_plus
  // Crk_minus
  // Cnr_plus
  // Cnr_minus
  // Cnk_plus
  // Cnk_minus
  // cHel_Mtt300_400
  vector<TH1F*> hists_cHel_Mtt300_400_up = {cHel_Mtt300_400_ele_reco_up, cHel_Mtt300_400_ele_id_up, cHel_Mtt300_400_ele_trigger_up, cHel_Mtt300_400_mu_reco_up, cHel_Mtt300_400_mu_iso_stat_up,cHel_Mtt300_400_mu_iso_syst_up, cHel_Mtt300_400_mu_id_stat_up, cHel_Mtt300_400_mu_id_syst_up, cHel_Mtt300_400_mu_trigger_stat_up, cHel_Mtt300_400_mu_trigger_syst_up, cHel_Mtt300_400_pu_up, cHel_Mtt300_400_prefiring_up};
  // cHel_Mtt300_400_betaLT0p9
  vector<TH1F*> hists_cHel_Mtt300_400_betaLT0p9_up = {cHel_Mtt300_400_betaLT0p9_ele_reco_up, cHel_Mtt300_400_betaLT0p9_ele_id_up, cHel_Mtt300_400_betaLT0p9_ele_trigger_up, cHel_Mtt300_400_betaLT0p9_mu_reco_up, cHel_Mtt300_400_betaLT0p9_mu_iso_stat_up,cHel_Mtt300_400_betaLT0p9_mu_iso_syst_up, cHel_Mtt300_400_betaLT0p9_mu_id_stat_up, cHel_Mtt300_400_betaLT0p9_mu_id_syst_up, cHel_Mtt300_400_betaLT0p9_mu_trigger_stat_up, cHel_Mtt300_400_betaLT0p9_mu_trigger_syst_up, cHel_Mtt300_400_betaLT0p9_pu_up, cHel_Mtt300_400_betaLT0p9_prefiring_up};
  // cHel_P3n_Mtt800_Inf
  vector<TH1F*> hists_cHel_P3n_Mtt800_Inf_up = {cHel_P3n_Mtt800_Inf_ele_reco_up, cHel_P3n_Mtt800_Inf_ele_id_up, cHel_P3n_Mtt800_Inf_ele_trigger_up, cHel_P3n_Mtt800_Inf_mu_reco_up, cHel_P3n_Mtt800_Inf_mu_iso_stat_up,cHel_P3n_Mtt800_Inf_mu_iso_syst_up, cHel_P3n_Mtt800_Inf_mu_id_stat_up, cHel_P3n_Mtt800_Inf_mu_id_syst_up, cHel_P3n_Mtt800_Inf_mu_trigger_stat_up, cHel_P3n_Mtt800_Inf_mu_trigger_syst_up, cHel_P3n_Mtt800_Inf_pu_up, cHel_P3n_Mtt800_Inf_prefiring_up};
  // cHel_P3n_Mtt800_Inf_cosThetaLT0p4
  vector<TH1F*> hists_cHel_P3n_Mtt800_Inf_cosThetaLT0p4_up = {cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_reco_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_id_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_trigger_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_reco_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_stat_up,cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_syst_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_stat_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_syst_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_stat_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_syst_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_pu_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_prefiring_up};
  // down variation histograms
  vector<TH1F*> hists_down   = {DeltaY_ele_reco_down, DeltaY_ele_id_down, DeltaY_ele_trigger_down, DeltaY_mu_reco_down, DeltaY_mu_iso_stat_down, DeltaY_mu_iso_syst_down, DeltaY_mu_id_stat_down, DeltaY_mu_id_syst_down, DeltaY_mu_trigger_stat_down, DeltaY_mu_trigger_syst_down, DeltaY_pu_down, DeltaY_prefiring_down};
  vector<TH1F*> hists_dy_d1_down = {DeltaY_reco_d1_ele_reco_down, DeltaY_reco_d1_ele_id_down, DeltaY_reco_d1_ele_trigger_down, DeltaY_reco_d1_mu_reco_down, DeltaY_reco_d1_mu_iso_stat_down,DeltaY_reco_d1_mu_iso_syst_down, DeltaY_reco_d1_mu_id_stat_down, DeltaY_reco_d1_mu_id_syst_down, DeltaY_reco_d1_mu_trigger_stat_down, DeltaY_reco_d1_mu_trigger_syst_down, DeltaY_reco_d1_pu_down, DeltaY_reco_d1_prefiring_down};
  vector<TH1F*> hists_dy_d2_down = {DeltaY_reco_d2_ele_reco_down, DeltaY_reco_d2_ele_id_down, DeltaY_reco_d2_ele_trigger_down, DeltaY_reco_d2_mu_reco_down, DeltaY_reco_d2_mu_iso_stat_down,DeltaY_reco_d2_mu_iso_syst_down, DeltaY_reco_d2_mu_id_stat_down, DeltaY_reco_d2_mu_id_syst_down, DeltaY_reco_d2_mu_trigger_stat_down, DeltaY_reco_d2_mu_trigger_syst_down, DeltaY_reco_d2_pu_down, DeltaY_reco_d2_prefiring_down};
  vector<TH1F*> hists_sigma_1_down = {Sigma_phi_1_ele_reco_down, Sigma_phi_1_ele_id_down, Sigma_phi_1_ele_trigger_down, Sigma_phi_1_mu_reco_down, Sigma_phi_1_mu_iso_stat_down, Sigma_phi_1_mu_iso_syst_down, Sigma_phi_1_mu_id_stat_down, Sigma_phi_1_mu_id_syst_down, Sigma_phi_1_mu_trigger_stat_down, Sigma_phi_1_mu_trigger_syst_down, Sigma_phi_1_pu_down, Sigma_phi_1_prefiring_down};
  vector<TH1F*> hists_sigma_2_down = {Sigma_phi_2_ele_reco_down, Sigma_phi_2_ele_id_down, Sigma_phi_2_ele_trigger_down, Sigma_phi_2_mu_reco_down, Sigma_phi_2_mu_iso_stat_down, Sigma_phi_2_mu_iso_syst_down, Sigma_phi_2_mu_id_stat_down, Sigma_phi_2_mu_id_syst_down, Sigma_phi_2_mu_trigger_stat_down, Sigma_phi_2_mu_trigger_syst_down, Sigma_phi_2_pu_down, Sigma_phi_2_prefiring_down};
  vector<TH1F*> hists_deltaY_xi_reco_6_down = {DeltaY_xi_reco_6_ele_reco_down, DeltaY_xi_reco_6_ele_id_down, DeltaY_xi_reco_6_ele_trigger_down, DeltaY_xi_reco_6_mu_reco_down, DeltaY_xi_reco_6_mu_iso_stat_down,DeltaY_xi_reco_6_mu_iso_syst_down, DeltaY_xi_reco_6_mu_id_stat_down, DeltaY_xi_reco_6_mu_id_syst_down, DeltaY_xi_reco_6_mu_trigger_stat_down, DeltaY_xi_reco_6_mu_trigger_syst_down, DeltaY_xi_reco_6_pu_down, DeltaY_xi_reco_6_prefiring_down};
  // cos_theta1k_antiLep
  // cos_theta1r_antiLep
  // cos_theta1n_antiLep
  // cos_theta1kStar_antiLep
  // cos_theta1rStar_antiLep
  // cos_theta2k_Lep
  // cos_theta2r_Lep
  // cos_theta2n_Lep
  // cos_theta2kStar_Lep
  // cos_theta2rStar_Lep
  // cos_theta1k
  // cos_theta1r
  // cos_theta1n
  // cos_theta1kStar
  // cos_theta1rStar
  // cos_theta2k
  // cos_theta2r
  // cos_theta2n
  // cos_theta2kStar
  // cos_theta2rStar
  // Ckk
  // Crr
  // Cnn
  // Crk_plus
  // Crk_minus
  // Cnr_plus
  // Cnr_minus
  // Cnk_plus
  // Cnk_minus
  // cHel_Mtt300_400
  vector<TH1F*> hists_cHel_Mtt300_400_down = {cHel_Mtt300_400_ele_reco_down, cHel_Mtt300_400_ele_id_down, cHel_Mtt300_400_ele_trigger_down, cHel_Mtt300_400_mu_reco_down, cHel_Mtt300_400_mu_iso_stat_down,cHel_Mtt300_400_mu_iso_syst_down, cHel_Mtt300_400_mu_id_stat_down, cHel_Mtt300_400_mu_id_syst_down, cHel_Mtt300_400_mu_trigger_stat_down, cHel_Mtt300_400_mu_trigger_syst_down, cHel_Mtt300_400_pu_down, cHel_Mtt300_400_prefiring_down};
  // cHel_Mtt300_400_betaLT0p9
  vector<TH1F*> hists_cHel_Mtt300_400_betaLT0p9_down = {cHel_Mtt300_400_betaLT0p9_ele_reco_down, cHel_Mtt300_400_betaLT0p9_ele_id_down, cHel_Mtt300_400_betaLT0p9_ele_trigger_down, cHel_Mtt300_400_betaLT0p9_mu_reco_down, cHel_Mtt300_400_betaLT0p9_mu_iso_stat_down,cHel_Mtt300_400_betaLT0p9_mu_iso_syst_down, cHel_Mtt300_400_betaLT0p9_mu_id_stat_down, cHel_Mtt300_400_betaLT0p9_mu_id_syst_down, cHel_Mtt300_400_betaLT0p9_mu_trigger_stat_down, cHel_Mtt300_400_betaLT0p9_mu_trigger_syst_down, cHel_Mtt300_400_betaLT0p9_pu_down, cHel_Mtt300_400_betaLT0p9_prefiring_down};
  // cHel_P3n_Mtt800_Inf
  vector<TH1F*> hists_cHel_P3n_Mtt800_Inf_down = {cHel_P3n_Mtt800_Inf_ele_reco_down, cHel_P3n_Mtt800_Inf_ele_id_down, cHel_P3n_Mtt800_Inf_ele_trigger_down, cHel_P3n_Mtt800_Inf_mu_reco_down, cHel_P3n_Mtt800_Inf_mu_iso_stat_down,cHel_P3n_Mtt800_Inf_mu_iso_syst_down, cHel_P3n_Mtt800_Inf_mu_id_stat_down, cHel_P3n_Mtt800_Inf_mu_id_syst_down, cHel_P3n_Mtt800_Inf_mu_trigger_stat_down, cHel_P3n_Mtt800_Inf_mu_trigger_syst_down, cHel_P3n_Mtt800_Inf_pu_down, cHel_P3n_Mtt800_Inf_prefiring_down};
  // cHel_P3n_Mtt800_Inf_cosThetaLT0p4
  vector<TH1F*> hists_cHel_P3n_Mtt800_Inf_cosThetaLT0p4_down = {cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_reco_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_id_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ele_trigger_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_reco_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_stat_down,cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_iso_syst_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_stat_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_id_syst_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_stat_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_mu_trigger_syst_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_pu_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_prefiring_down};


  // b-tagging variations
  vector<float> syst_btag   = {btag_cferr1_up, btag_cferr1_down, btag_cferr2_up, btag_cferr2_down, btag_hf_up, btag_hf_down, btag_hfstats1_up, btag_hfstats1_down, btag_hfstats2_up, btag_hfstats2_down, btag_lf_up, btag_lf_down, btag_lfstats1_up, btag_lfstats1_down, btag_lfstats2_up, btag_lfstats2_down};
  vector<TH1F*> hists_btag  = {DeltaY_btag_cferr1_up, DeltaY_btag_cferr1_down, DeltaY_btag_cferr2_up, DeltaY_btag_cferr2_down, DeltaY_btag_hf_up, DeltaY_btag_hf_down, DeltaY_btag_hfstats1_up, DeltaY_btag_hfstats1_down, DeltaY_btag_hfstats2_up, DeltaY_btag_hfstats2_down, DeltaY_btag_lf_up, DeltaY_btag_lf_down, DeltaY_btag_lfstats1_up, DeltaY_btag_lfstats1_down, DeltaY_btag_lfstats2_up, DeltaY_btag_lfstats2_down};
  vector<TH1F*> hists_btag_dy_d1 = {DeltaY_reco_d1_btag_cferr1_up, DeltaY_reco_d1_btag_cferr1_down, DeltaY_reco_d1_btag_cferr2_up, DeltaY_reco_d1_btag_cferr2_down, DeltaY_reco_d1_btag_hf_up, DeltaY_reco_d1_btag_hf_down, DeltaY_reco_d1_btag_hfstats1_up, DeltaY_reco_d1_btag_hfstats1_down, DeltaY_reco_d1_btag_hfstats2_up, DeltaY_reco_d1_btag_hfstats2_down, DeltaY_reco_d1_btag_lf_up, DeltaY_reco_d1_btag_lf_down, DeltaY_reco_d1_btag_lfstats1_up, DeltaY_reco_d1_btag_lfstats1_down, DeltaY_reco_d1_btag_lfstats2_up, DeltaY_reco_d1_btag_lfstats2_down};
  vector<TH1F*> hists_btag_dy_d2 = {DeltaY_reco_d2_btag_cferr1_up, DeltaY_reco_d2_btag_cferr1_down, DeltaY_reco_d2_btag_cferr2_up, DeltaY_reco_d2_btag_cferr2_down, DeltaY_reco_d2_btag_hf_up, DeltaY_reco_d2_btag_hf_down, DeltaY_reco_d2_btag_hfstats1_up, DeltaY_reco_d2_btag_hfstats1_down, DeltaY_reco_d2_btag_hfstats2_up, DeltaY_reco_d2_btag_hfstats2_down, DeltaY_reco_d2_btag_lf_up, DeltaY_reco_d2_btag_lf_down, DeltaY_reco_d2_btag_lfstats1_up, DeltaY_reco_d2_btag_lfstats1_down, DeltaY_reco_d2_btag_lfstats2_up, DeltaY_reco_d2_btag_lfstats2_down};
  vector<TH1F*> hists_btag_sigma_1 = {Sigma_phi_1_btag_cferr1_up, Sigma_phi_1_btag_cferr1_down, Sigma_phi_1_btag_cferr2_up, Sigma_phi_1_btag_cferr2_down, Sigma_phi_1_btag_hf_up, Sigma_phi_1_btag_hf_down, Sigma_phi_1_btag_hfstats1_up, Sigma_phi_1_btag_hfstats1_down, Sigma_phi_1_btag_hfstats2_up, Sigma_phi_1_btag_hfstats2_down, Sigma_phi_1_btag_lf_up, Sigma_phi_1_btag_lf_down, Sigma_phi_1_btag_lfstats1_up, Sigma_phi_1_btag_lfstats1_down, Sigma_phi_1_btag_lfstats2_up, Sigma_phi_1_btag_lfstats2_down};
  vector<TH1F*> hists_btag_sigma_2 = {Sigma_phi_2_btag_cferr1_up, Sigma_phi_2_btag_cferr1_down, Sigma_phi_2_btag_cferr2_up, Sigma_phi_2_btag_cferr2_down, Sigma_phi_2_btag_hf_up, Sigma_phi_2_btag_hf_down, Sigma_phi_2_btag_hfstats1_up, Sigma_phi_2_btag_hfstats1_down, Sigma_phi_2_btag_hfstats2_up, Sigma_phi_2_btag_hfstats2_down, Sigma_phi_2_btag_lf_up, Sigma_phi_2_btag_lf_down, Sigma_phi_2_btag_lfstats1_up, Sigma_phi_2_btag_lfstats1_down, Sigma_phi_2_btag_lfstats2_up, Sigma_phi_2_btag_lfstats2_down};
  vector<TH1F*> hists_btag_deltaY_xi_reco_6 = {DeltaY_xi_reco_6_btag_cferr1_up, DeltaY_xi_reco_6_btag_cferr1_down, DeltaY_xi_reco_6_btag_cferr2_up, DeltaY_xi_reco_6_btag_cferr2_down, DeltaY_xi_reco_6_btag_hf_up, DeltaY_xi_reco_6_btag_hf_down, DeltaY_xi_reco_6_btag_hfstats1_up, DeltaY_xi_reco_6_btag_hfstats1_down, DeltaY_xi_reco_6_btag_hfstats2_up, DeltaY_xi_reco_6_btag_hfstats2_down, DeltaY_xi_reco_6_btag_lf_up, DeltaY_xi_reco_6_btag_lf_down, DeltaY_xi_reco_6_btag_lfstats1_up, DeltaY_xi_reco_6_btag_lfstats1_down, DeltaY_xi_reco_6_btag_lfstats2_up, DeltaY_xi_reco_6_btag_lfstats2_down};
  // cHel_Mtt300_400
  vector<TH1F*> hists_btag_cHel_Mtt300_400 = {cHel_Mtt300_400_btag_cferr1_up, cHel_Mtt300_400_btag_cferr1_down, cHel_Mtt300_400_btag_cferr2_up, cHel_Mtt300_400_btag_cferr2_down, cHel_Mtt300_400_btag_hf_up, cHel_Mtt300_400_btag_hf_down, cHel_Mtt300_400_btag_hfstats1_up, cHel_Mtt300_400_btag_hfstats1_down, cHel_Mtt300_400_btag_hfstats2_up, cHel_Mtt300_400_btag_hfstats2_down, cHel_Mtt300_400_btag_lf_up, cHel_Mtt300_400_btag_lf_down, cHel_Mtt300_400_btag_lfstats1_up, cHel_Mtt300_400_btag_lfstats1_down, cHel_Mtt300_400_btag_lfstats2_up, cHel_Mtt300_400_btag_lfstats2_down};
  // cHel_Mtt300_400_betaLT0p9
  vector<TH1F*> hists_btag_cHel_Mtt300_400_betaLT0p9 = {cHel_Mtt300_400_betaLT0p9_btag_cferr1_up, cHel_Mtt300_400_betaLT0p9_btag_cferr1_down, cHel_Mtt300_400_betaLT0p9_btag_cferr2_up, cHel_Mtt300_400_betaLT0p9_btag_cferr2_down, cHel_Mtt300_400_betaLT0p9_btag_hf_up, cHel_Mtt300_400_betaLT0p9_btag_hf_down, cHel_Mtt300_400_betaLT0p9_btag_hfstats1_up, cHel_Mtt300_400_betaLT0p9_btag_hfstats1_down, cHel_Mtt300_400_betaLT0p9_btag_hfstats2_up, cHel_Mtt300_400_betaLT0p9_btag_hfstats2_down, cHel_Mtt300_400_betaLT0p9_btag_lf_up, cHel_Mtt300_400_betaLT0p9_btag_lf_down, cHel_Mtt300_400_betaLT0p9_btag_lfstats1_up, cHel_Mtt300_400_betaLT0p9_btag_lfstats1_down, cHel_Mtt300_400_betaLT0p9_btag_lfstats2_up, cHel_Mtt300_400_betaLT0p9_btag_lfstats2_down};
  // cHel_P3n_Mtt800_Inf
  vector<TH1F*> hists_btag_cHel_P3n_Mtt800_Inf = {cHel_P3n_Mtt800_Inf_btag_cferr1_up, cHel_P3n_Mtt800_Inf_btag_cferr1_down, cHel_P3n_Mtt800_Inf_btag_cferr2_up, cHel_P3n_Mtt800_Inf_btag_cferr2_down, cHel_P3n_Mtt800_Inf_btag_hf_up, cHel_P3n_Mtt800_Inf_btag_hf_down, cHel_P3n_Mtt800_Inf_btag_hfstats1_up, cHel_P3n_Mtt800_Inf_btag_hfstats1_down, cHel_P3n_Mtt800_Inf_btag_hfstats2_up, cHel_P3n_Mtt800_Inf_btag_hfstats2_down, cHel_P3n_Mtt800_Inf_btag_lf_up, cHel_P3n_Mtt800_Inf_btag_lf_down, cHel_P3n_Mtt800_Inf_btag_lfstats1_up, cHel_P3n_Mtt800_Inf_btag_lfstats1_down, cHel_P3n_Mtt800_Inf_btag_lfstats2_up, cHel_P3n_Mtt800_Inf_btag_lfstats2_down};
  // cHel_P3n_Mtt800_Inf_cosThetaLT0p4
  vector<TH1F*> hists_btag_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = {cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr1_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr1_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr2_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_cferr2_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hf_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hf_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats1_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats1_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats2_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_hfstats2_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lf_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lf_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats1_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats1_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats2_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_btag_lfstats2_down};

  
  // Top tagging variations
  vector<float> syst_ttag = {ttag_corr_up, ttag_corr_down, ttag_uncorr_up, ttag_uncorr_down};
  vector<TH1F*> hists_ttag = {DeltaY_ttag_corr_up, DeltaY_ttag_corr_down, DeltaY_ttag_uncorr_up, DeltaY_ttag_uncorr_down};
  vector<TH1F*> hists_ttag_dy_d1 = {DeltaY_reco_d1_ttag_corr_up, DeltaY_reco_d1_ttag_corr_down, DeltaY_reco_d1_ttag_uncorr_up, DeltaY_reco_d1_ttag_uncorr_down};
  vector<TH1F*> hists_ttag_dy_d2 = {DeltaY_reco_d2_ttag_corr_up, DeltaY_reco_d2_ttag_corr_down, DeltaY_reco_d2_ttag_uncorr_up, DeltaY_reco_d2_ttag_uncorr_down};
  vector<TH1F*> hists_ttag_sigma_1 = {Sigma_phi_1_ttag_corr_up, Sigma_phi_1_ttag_corr_down, Sigma_phi_1_ttag_uncorr_up, Sigma_phi_1_ttag_uncorr_down};
  vector<TH1F*> hists_ttag_sigma_2 = {Sigma_phi_2_ttag_corr_up, Sigma_phi_2_ttag_corr_down, Sigma_phi_2_ttag_uncorr_up, Sigma_phi_2_ttag_uncorr_down};
  vector<TH1F*> hists_ttag_deltaY_xi_reco_6 = {DeltaY_xi_reco_6_ttag_corr_up, DeltaY_xi_reco_6_ttag_corr_down, DeltaY_xi_reco_6_ttag_uncorr_up, DeltaY_xi_reco_6_ttag_uncorr_down};
  // cHel_Mtt300_400
  vector<TH1F*> hists_ttag_cHel_Mtt300_400 = {cHel_Mtt300_400_ttag_corr_up, cHel_Mtt300_400_ttag_corr_down, cHel_Mtt300_400_ttag_uncorr_up, cHel_Mtt300_400_ttag_uncorr_down};
  // cHel_Mtt300_400_betaLT0p9
  vector<TH1F*> hists_ttag_cHel_Mtt300_400_betaLT0p9 = {cHel_Mtt300_400_betaLT0p9_ttag_corr_up, cHel_Mtt300_400_betaLT0p9_ttag_corr_down, cHel_Mtt300_400_betaLT0p9_ttag_uncorr_up, cHel_Mtt300_400_betaLT0p9_ttag_uncorr_down};
  // cHel_P3n_Mtt800_Inf
  vector<TH1F*> hists_ttag_cHel_P3n_Mtt800_Inf = {cHel_P3n_Mtt800_Inf_ttag_corr_up, cHel_P3n_Mtt800_Inf_ttag_corr_down, cHel_P3n_Mtt800_Inf_ttag_uncorr_up, cHel_P3n_Mtt800_Inf_ttag_uncorr_down};
  // cHel_P3n_Mtt800_Inf_cosThetaLT0p4
  vector<TH1F*> hists_ttag_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = {cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_corr_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_corr_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_uncorr_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_ttag_uncorr_down};

  
  // tmistag variations
  vector<float> syst_tmistag = {tmistag_up, tmistag_down};
  vector<TH1F*> hists_tmistag = {DeltaY_tmistag_up, DeltaY_tmistag_down};
  vector<TH1F*> hists_tmistag_dy_d1 = {DeltaY_reco_d1_tmistag_up, DeltaY_reco_d1_tmistag_down};
  vector<TH1F*> hists_tmistag_dy_d2 = {DeltaY_reco_d2_tmistag_up, DeltaY_reco_d2_tmistag_down}; 
  vector<TH1F*> hists_tmistag_sigma_1 = {Sigma_phi_1_tmistag_up, Sigma_phi_1_tmistag_down};
  vector<TH1F*> hists_tmistag_sigma_2 = {Sigma_phi_2_tmistag_up, Sigma_phi_2_tmistag_down};
  vector<TH1F*> hists_tmistag_deltaY_xi_reco_6 = {DeltaY_xi_reco_6_tmistag_up, DeltaY_xi_reco_6_tmistag_down};
  // cHel_Mtt300_400
  vector<TH1F*> hists_tmistag_cHel_Mtt300_400 = {cHel_Mtt300_400_tmistag_up, cHel_Mtt300_400_tmistag_down};
  // cHel_Mtt300_400_betaLT0p9
  vector<TH1F*> hists_tmistag_cHel_Mtt300_400_betaLT0p9 = {cHel_Mtt300_400_betaLT0p9_tmistag_up, cHel_Mtt300_400_betaLT0p9_tmistag_down};
  // cHel_P3n_Mtt800_Inf
  vector<TH1F*> hists_tmistag_cHel_P3n_Mtt800_Inf = {cHel_P3n_Mtt800_Inf_tmistag_up, cHel_P3n_Mtt800_Inf_tmistag_down};
  // cHel_P3n_Mtt800_Inf_cosThetaLT0p4
  vector<TH1F*> hists_tmistag_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = {cHel_P3n_Mtt800_Inf_cosThetaLT0p4_tmistag_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_tmistag_down};
  
  
  //Top pt reweighting
  vector<float> syst_toppt = {toppt_a_up, toppt_a_down, toppt_b_up, toppt_b_down};
  vector<TH1F*> hists_toppt = {DeltaY_toppt_a_up, DeltaY_toppt_a_down, DeltaY_toppt_b_up, DeltaY_toppt_b_down};
  vector<TH1F*> hists_toppt_dy_d1 = {DeltaY_reco_d1_toppt_a_up, DeltaY_reco_d1_toppt_a_down, DeltaY_reco_d1_toppt_b_up, DeltaY_reco_d1_toppt_b_down};
  vector<TH1F*> hists_toppt_dy_d2 = {DeltaY_reco_d2_toppt_a_up, DeltaY_reco_d2_toppt_a_down, DeltaY_reco_d2_toppt_b_up, DeltaY_reco_d2_toppt_b_down};
  vector<TH1F*> hists_toppt_sigma_1 = {Sigma_phi_1_toppt_a_up, Sigma_phi_1_toppt_a_down, Sigma_phi_1_toppt_b_up, Sigma_phi_1_toppt_b_down};
  vector<TH1F*> hists_toppt_sigma_2 = {Sigma_phi_2_toppt_a_up, Sigma_phi_2_toppt_a_down, Sigma_phi_2_toppt_b_up, Sigma_phi_2_toppt_b_down};
  vector<TH1F*> hists_toppt_deltaY_xi_reco_6 = {DeltaY_xi_reco_6_toppt_a_up, DeltaY_xi_reco_6_toppt_a_down, DeltaY_xi_reco_6_toppt_b_up, DeltaY_xi_reco_6_toppt_b_down};
  // cHel_Mtt300_400
  vector<TH1F*> hists_toppt_cHel_Mtt300_400 = {cHel_Mtt300_400_toppt_a_up, cHel_Mtt300_400_toppt_a_down, cHel_Mtt300_400_toppt_b_up, cHel_Mtt300_400_toppt_b_down};
  // cHel_Mtt300_400_betaLT0p9
  vector<TH1F*> hists_toppt_cHel_Mtt300_400_betaLT0p9 = {cHel_Mtt300_400_betaLT0p9_toppt_a_up, cHel_Mtt300_400_betaLT0p9_toppt_a_down, cHel_Mtt300_400_betaLT0p9_toppt_b_up, cHel_Mtt300_400_betaLT0p9_toppt_b_down};
  // cHel_P3n_Mtt800_Inf
  vector<TH1F*> hists_toppt_cHel_P3n_Mtt800_Inf = {cHel_P3n_Mtt800_Inf_toppt_a_up, cHel_P3n_Mtt800_Inf_toppt_a_down, cHel_P3n_Mtt800_Inf_toppt_b_up, cHel_P3n_Mtt800_Inf_toppt_b_down};
  // cHel_P3n_Mtt800_Inf_cosThetaLT0p4
  vector<TH1F*> hists_toppt_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = {cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_a_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_a_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_b_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_toppt_b_down};
  
  
  // muR and muF variations
  vector<float> syst_scale  = {murmuf_upup, murmuf_upnone, murmuf_noneup, murmuf_nonedown, murmuf_downnone, murmuf_downdown};
  vector<TH1F*> hists_scale = {DeltaY_murmuf_upup, DeltaY_murmuf_upnone, DeltaY_murmuf_noneup, DeltaY_murmuf_nonedown, DeltaY_murmuf_downnone, DeltaY_murmuf_downdown};
  vector<TH1F*> hists_scale_dy_d1 = {DeltaY_reco_d1_murmuf_upup, DeltaY_reco_d1_murmuf_upnone, DeltaY_reco_d1_murmuf_noneup, DeltaY_reco_d1_murmuf_nonedown, DeltaY_reco_d1_murmuf_downnone, DeltaY_reco_d1_murmuf_downdown};
  vector<TH1F*> hists_scale_dy_d2 = {DeltaY_reco_d2_murmuf_upup, DeltaY_reco_d2_murmuf_upnone, DeltaY_reco_d2_murmuf_noneup, DeltaY_reco_d2_murmuf_nonedown, DeltaY_reco_d2_murmuf_downnone, DeltaY_reco_d2_murmuf_downdown};
  vector<TH1F*> hists_scale_sigma_1 = {Sigma_phi_1_murmuf_upup, Sigma_phi_1_murmuf_upnone, Sigma_phi_1_murmuf_noneup, Sigma_phi_1_murmuf_nonedown, Sigma_phi_1_murmuf_downnone, Sigma_phi_1_murmuf_downdown};
  vector<TH1F*> hists_scale_sigma_2 = {Sigma_phi_2_murmuf_upup, Sigma_phi_2_murmuf_upnone, Sigma_phi_2_murmuf_noneup, Sigma_phi_2_murmuf_nonedown, Sigma_phi_2_murmuf_downnone, Sigma_phi_2_murmuf_downdown};
  vector<TH1F*> hists_scale_deltaY_xi_reco_6 = {DeltaY_xi_reco_6_murmuf_upup, DeltaY_xi_reco_6_murmuf_upnone, DeltaY_xi_reco_6_murmuf_noneup, DeltaY_xi_reco_6_murmuf_nonedown, DeltaY_xi_reco_6_murmuf_downnone, DeltaY_xi_reco_6_murmuf_downdown};
  // cHel_Mtt300_400
  vector<TH1F*> hists_scale_cHel_Mtt300_400 = {cHel_Mtt300_400_murmuf_upup, cHel_Mtt300_400_murmuf_upnone, cHel_Mtt300_400_murmuf_noneup, cHel_Mtt300_400_murmuf_nonedown, cHel_Mtt300_400_murmuf_downnone, cHel_Mtt300_400_murmuf_downdown};
  // cHel_Mtt300_400_betaLT0p9
  vector<TH1F*> hists_scale_cHel_Mtt300_400_betaLT0p9 = {cHel_Mtt300_400_betaLT0p9_murmuf_upup, cHel_Mtt300_400_betaLT0p9_murmuf_upnone, cHel_Mtt300_400_betaLT0p9_murmuf_noneup, cHel_Mtt300_400_betaLT0p9_murmuf_nonedown, cHel_Mtt300_400_betaLT0p9_murmuf_downnone, cHel_Mtt300_400_betaLT0p9_murmuf_downdown};
  // cHel_P3n_Mtt800_Inf
  vector<TH1F*> hists_scale_cHel_P3n_Mtt800_Inf = {cHel_P3n_Mtt800_Inf_murmuf_upup, cHel_P3n_Mtt800_Inf_murmuf_upnone, cHel_P3n_Mtt800_Inf_murmuf_noneup, cHel_P3n_Mtt800_Inf_murmuf_nonedown, cHel_P3n_Mtt800_Inf_murmuf_downnone, cHel_P3n_Mtt800_Inf_murmuf_downdown};
  // cHel_P3n_Mtt800_Inf_cosThetaLT0p4
  vector<TH1F*> hists_scale_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = {cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_upup, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_upnone, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_noneup, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_nonedown, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_downnone, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_murmuf_downdown};
  
  
  // parton shower variations (ISR, FSR)
  vector<float> syst_ps = {isr_up, isr_down, fsr_up, fsr_down};
  vector<TH1F*> hists_ps = {DeltaY_isr_up, DeltaY_isr_down, DeltaY_fsr_up, DeltaY_fsr_down}; 
  vector<TH1F*> hists_ps_dy_d1 = {DeltaY_reco_d1_isr_up, DeltaY_reco_d1_isr_down, DeltaY_reco_d1_fsr_up, DeltaY_reco_d1_fsr_down}; 
  vector<TH1F*> hists_ps_dy_d2 = {DeltaY_reco_d2_isr_up, DeltaY_reco_d2_isr_down, DeltaY_reco_d2_fsr_up, DeltaY_reco_d2_fsr_down};
  vector<TH1F*> hists_ps_sigma_1 = {Sigma_phi_1_isr_up, Sigma_phi_1_isr_down, Sigma_phi_1_fsr_up, Sigma_phi_1_fsr_down};
  vector<TH1F*> hists_ps_sigma_2 = {Sigma_phi_2_isr_up, Sigma_phi_2_isr_down, Sigma_phi_2_fsr_up, Sigma_phi_2_fsr_down};
  vector<TH1F*> hists_ps_deltaY_xi_reco_6 = {DeltaY_xi_reco_6_isr_up, DeltaY_xi_reco_6_isr_down, DeltaY_xi_reco_6_fsr_up, DeltaY_xi_reco_6_fsr_down};
  // cHel_Mtt300_400
  vector<TH1F*> hists_ps_cHel_Mtt300_400 = {cHel_Mtt300_400_isr_up, cHel_Mtt300_400_isr_down, cHel_Mtt300_400_fsr_up, cHel_Mtt300_400_fsr_down}; 
  // cHel_Mtt300_400_betaLT0p9
  vector<TH1F*> hists_ps_cHel_Mtt300_400_betaLT0p9 = {cHel_Mtt300_400_betaLT0p9_isr_up, cHel_Mtt300_400_betaLT0p9_isr_down, cHel_Mtt300_400_betaLT0p9_fsr_up, cHel_Mtt300_400_betaLT0p9_fsr_down}; 
  // cHel_P3n_Mtt800_Inf
  vector<TH1F*> hists_ps_cHel_P3n_Mtt800_Inf = {cHel_P3n_Mtt800_Inf_isr_up, cHel_P3n_Mtt800_Inf_isr_down, cHel_P3n_Mtt800_Inf_fsr_up, cHel_P3n_Mtt800_Inf_fsr_down}; 
  // cHel_P3n_Mtt800_Inf_cosThetaLT0p4
  vector<TH1F*> hists_ps_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = {cHel_P3n_Mtt800_Inf_cosThetaLT0p4_isr_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_isr_down, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_fsr_up, cHel_P3n_Mtt800_Inf_cosThetaLT0p4_fsr_down}; 
  
  
  //-------------------------------------------------------------------------- ttbar 2D --------------------------------------------------------------------------//
  // up/down variations
  vector<TH2F*> hists_up_tt  = {DeltaY_ele_reco_up_tt, DeltaY_ele_id_up_tt, DeltaY_ele_trigger_up_tt, DeltaY_mu_reco_up_tt, DeltaY_mu_iso_stat_up_tt, DeltaY_mu_id_stat_up_tt, DeltaY_mu_trigger_stat_up_tt, DeltaY_mu_iso_syst_up_tt, DeltaY_mu_id_syst_up_tt, DeltaY_mu_trigger_syst_up_tt,  DeltaY_pu_up_tt, DeltaY_prefiring_up_tt};
  vector<TH2F*> hists_down_tt= {DeltaY_ele_reco_down_tt, DeltaY_ele_id_down_tt, DeltaY_ele_trigger_down_tt, DeltaY_mu_reco_down_tt, DeltaY_mu_iso_stat_down_tt, DeltaY_mu_id_stat_down_tt, DeltaY_mu_trigger_stat_down_tt, DeltaY_mu_iso_syst_down_tt, DeltaY_mu_id_syst_down_tt, DeltaY_mu_trigger_syst_down_tt, DeltaY_pu_down_tt, DeltaY_prefiring_down_tt};
  // ttbar 2D btag variations
  vector<TH2F*> hists_btag_tt  = {DeltaY_btag_cferr1_up_tt, DeltaY_btag_cferr1_down_tt, DeltaY_btag_cferr2_up_tt, DeltaY_btag_cferr2_down_tt, DeltaY_btag_hf_up_tt, DeltaY_btag_hf_down_tt, DeltaY_btag_hfstats1_up_tt, DeltaY_btag_hfstats1_down_tt, DeltaY_btag_hfstats2_up_tt, DeltaY_btag_hfstats2_down_tt, DeltaY_btag_lf_up_tt, DeltaY_btag_lf_down_tt, DeltaY_btag_lfstats1_up_tt, DeltaY_btag_lfstats1_down_tt, DeltaY_btag_lfstats2_up_tt, DeltaY_btag_lfstats2_down_tt};
  // ttbar 2D ttag variations
  vector<TH2F*> hists_ttag_tt = {DeltaY_ttag_corr_up_tt, DeltaY_ttag_corr_down_tt, DeltaY_ttag_uncorr_up_tt, DeltaY_ttag_uncorr_down_tt};
  // ttbar 2D tmistag variations
  vector<TH2F*> hists_tmistag_tt = {DeltaY_tmistag_up_tt, DeltaY_tmistag_down_tt};
  // ttbar 2D toppt variations
  vector<TH2F*> hists_toppt_tt = {DeltaY_toppt_a_up_tt, DeltaY_toppt_a_down_tt, DeltaY_toppt_b_up_tt, DeltaY_toppt_b_down_tt};
  // ttbar 2D scale variations
  vector<TH2F*> hists_scale_tt = {DeltaY_murmuf_upup_tt, DeltaY_murmuf_upnone_tt, DeltaY_murmuf_noneup_tt, DeltaY_murmuf_nonedown_tt, DeltaY_murmuf_downnone_tt, DeltaY_murmuf_downdown_tt};
  // ttbar 2D parton shower variations
  vector<TH2F*> hists_ps_tt = {DeltaY_isr_up_tt, DeltaY_isr_down_tt, DeltaY_fsr_up_tt, DeltaY_fsr_down_tt};


  bool debug=false;
  // Zprime reco
  bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2);
  ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
  if(is_zprime_reconstructed_chi2 && is_mc){
    if(is_tt){
      if (debug)cout << "check ttbar all sys for deltay RM" <<endl;
      const auto& genparticles = event.genparticles;
      // ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);

      GenParticle top, antitop;
      for(const GenParticle & gp : *genparticles){
        if(gp.pdgId() == 6){
            top = gp;
        }
        else if(gp.pdgId() == -6){
            antitop = gp;
        }
      }
       if (debug)cout << "ttbar loop over gen" <<endl;
      // The Lorentz vectors represent the 4-momenta (energy, and three spatial momentum components) for the leptonic and hadronic tops from the "BestZprimeCandidate" object
      LorentzVector lep_top = BestZprimeCandidate->top_leptonic_v4();
      LorentzVector had_top = BestZprimeCandidate->top_hadronic_v4();
      if (debug)cout << "ttbar all sys define best cand leg" <<endl;
      // vectors to store the deltaR values for the leptonic and hadronic tops with each gen particle
      // this part initializes vectors to store deltaR values with a default of 99.0 and fills in the actual deltaR values by looping over the gen particles (top)
      std::vector<std::pair<double, int>> deltaR_leptonic_values; // ((dR, index), (dR, index), ...)
      std::vector<std::pair<double, int>> deltaR_hadronic_values;

      double deltaR_min_leptonic = 99.0;
      double deltaR_sec_min_leptonic = 99.0;
      int best_gen_for_leptop = -1;
      int sec_best_gen_for_leptop = -1;
      // bool is_leptop_matched = false;

      double deltaR_min_hadronic = 99.0;
      double deltaR_sec_min_hadronic = 99.0;
      int best_gen_for_hadtop = -1;
      int sec_best_gen_for_hadtop = -1;
      // bool is_hadtop_matched = false;

      for(unsigned int j=0; j<genparticles->size(); ++j) {
        if(abs(genparticles->at(j).pdgId()) == 6 ){
          if (genparticles->at(j).index() == 2 || genparticles->at(j).index() == 3){
            LorentzVector genparticle_p4(genparticles->at(j).pt(), genparticles->at(j).eta(), genparticles->at(j).phi(), genparticles->at(j).energy());
            deltaR_leptonic_values.push_back(std::make_pair(deltaR(lep_top, genparticle_p4),genparticles->at(j).index() ));
            deltaR_hadronic_values.push_back(std::make_pair(deltaR(had_top, genparticle_p4), genparticles->at(j).index()));
            // if (debug)cout << "deltaR: " << deltaR(lep_top, genparticle_p4) << j << endl;
          }
        }
      }  
      if (debug)cout << "ttbar all sys assign gen particles" <<endl;
      for (const auto& pair_lep : deltaR_leptonic_values) {
        if (pair_lep.first > 0 && pair_lep.first < deltaR_min_leptonic) {
          // deltaR_min_leptonic = pair_lep.first;
          deltaR_sec_min_leptonic = deltaR_min_leptonic;
          deltaR_min_leptonic = pair_lep.first;
          sec_best_gen_for_leptop = best_gen_for_leptop;
          best_gen_for_leptop = pair_lep.second;
          // is_leptop_matched = true;
        }
        else if (pair_lep.first > 0 && pair_lep.first < deltaR_sec_min_leptonic && pair_lep.first != deltaR_min_leptonic && pair_lep.second != best_gen_for_leptop) {
          deltaR_sec_min_leptonic = pair_lep.first;
          sec_best_gen_for_leptop = pair_lep.second;
        }
      }
       if (debug)cout << "ttbar all sys assign gen particles matching" <<endl;
      // deltaY values calculation starts:

      // matched gen particles
      for (const auto& pair_had : deltaR_hadronic_values) {
        if (pair_had.first > 0 && pair_had.first < deltaR_min_hadronic) {
          deltaR_sec_min_hadronic = deltaR_min_hadronic;
          deltaR_min_hadronic = pair_had.first;
          sec_best_gen_for_hadtop = best_gen_for_hadtop;
          best_gen_for_hadtop = pair_had.second;
          // is_hadtop_matched = true;
        }
        else if (pair_had.first > 0 && pair_had.first < deltaR_sec_min_hadronic && pair_had.first != deltaR_min_hadronic && pair_had.second != best_gen_for_hadtop) {
          deltaR_sec_min_hadronic = pair_had.first;
          sec_best_gen_for_hadtop = pair_had.second;
        }
      }

      if (best_gen_for_hadtop == best_gen_for_leptop){
        if (deltaR_min_leptonic <= deltaR_min_hadronic){
          best_gen_for_hadtop = sec_best_gen_for_hadtop;
          deltaR_min_hadronic = deltaR_sec_min_hadronic;
        } else {
          best_gen_for_leptop = sec_best_gen_for_leptop;
          deltaR_min_leptonic = deltaR_sec_min_leptonic;
        }
      }

      GenParticle best_matched_gen_leptop;
      GenParticle best_matched_gen_hadtop;

      float_t DeltaY_gen_best = 99.0;
      float_t DeltaY_reco_best = 99.0;
      bool isLeptonPositive = false;
       if (debug)cout << "ttbar all sys calc deltay " <<endl;
      if(isMuon){
        if (event.muons->at(0).charge() == 1){
          isLeptonPositive = true;
        } else {
          isLeptonPositive = false;
        }
      }

      if(isElectron){
        if (event.electrons->at(0).charge() == 1){
          isLeptonPositive = true;
        } else {
          isLeptonPositive = false;
        }
      }

      if (deltaR_min_leptonic < 0.4 && deltaR_min_hadronic < 0.4 && best_gen_for_leptop >= 0 && best_gen_for_hadtop >= 0) {
     
        if(static_cast<std::size_t>(best_gen_for_leptop) < genparticles->size()) {
          best_matched_gen_leptop = genparticles->at(best_gen_for_leptop);
        }
        
        if(static_cast<std::size_t>(best_gen_for_hadtop) < genparticles->size()) {
            best_matched_gen_hadtop = genparticles->at(best_gen_for_hadtop);
            // if( ) if (debug)cout << "hadtop match done" << endl;
        }



        // Calculates the delta y (with reco particles) values for the leptonic and hadronic tops depending on the charge of the lepton
        if (isLeptonPositive) {
          //if (debug)cout<<"lepton is positive"<<endl;
          DeltaY_reco_best = TMath::Abs(0.5*TMath::Log((lep_top.energy() + lep_top.pt()*TMath::SinH(lep_top.eta()))/(lep_top.energy() - lep_top.pt()*TMath::SinH(lep_top.eta())))) - TMath::Abs(0.5*TMath::Log((had_top.energy() + had_top.pt()*TMath::SinH(had_top.eta()))/(had_top.energy() - had_top.pt()*TMath::SinH(had_top.eta()))));
          DeltaY_gen_best = TMath::Abs(0.5*TMath::Log((best_matched_gen_leptop.energy() + best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta()))/(best_matched_gen_leptop.energy() - best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta())))) - TMath::Abs(0.5*TMath::Log((best_matched_gen_hadtop.energy() + best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta()))/(best_matched_gen_hadtop.energy() - best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta()))));

        } else {
          //if (debug)cout<<"lepton is negative"<<endl;
          DeltaY_reco_best = TMath::Abs(0.5*TMath::Log((had_top.energy() + had_top.pt()*TMath::SinH(had_top.eta()))/(had_top.energy() - had_top.pt()*TMath::SinH(had_top.eta())))) - TMath::Abs(0.5*TMath::Log((lep_top.energy() + lep_top.pt()*TMath::SinH(lep_top.eta()))/(lep_top.energy() - lep_top.pt()*TMath::SinH(lep_top.eta()))));
          DeltaY_gen_best = TMath::Abs(0.5*TMath::Log((best_matched_gen_hadtop.energy() + best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta()))/(best_matched_gen_hadtop.energy() - best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta())))) - TMath::Abs(0.5*TMath::Log((best_matched_gen_leptop.energy() + best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta()))/(best_matched_gen_leptop.energy() - best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta()))));
        }
      
      }

      if (debug)cout << "ttbar all sys fill deltay " <<endl;
      DeltaY_tt->Fill(DeltaY_reco_best, DeltaY_gen_best, weight);
      if (debug)cout << "ttbar all sys fill deltay-nominal " <<endl;
      
      // STEP 7: Calculate the reco and gen xi values. Event by event.
      //template method start
      // --- Template method: xi = tanh(DeltaY_reco) built from reco tops (not best/gen-matched)
      double deltay_reco = 0.0;
      if (isLeptonPositive) {
        deltay_reco = TMath::Abs(lep_top.Rapidity()) - TMath::Abs(had_top.Rapidity());
      } else {
        deltay_reco = TMath::Abs(had_top.Rapidity()) - TMath::Abs(lep_top.Rapidity());
      }

      const double deltay_xi = TMath::TanH(deltay_reco);
      double xi_gen_evt = static_cast<double>(event.get(h_xi_gen));
      // Clamp xi_gen_evt to valid range to avoid edge effects in weight lookup
      if(xi_gen_evt <= -1.0) xi_gen_evt = -0.999999;
      if(xi_gen_evt >= +1.0) xi_gen_evt = +0.999999;

      // Fill base xi nominals (matching Hists behavior: with NoAC weights if enabled)
      if(use_noac_evtweights_ && !noac_weights_map.empty() && noac_weights_map.count(0.0f)){
        // Use f=0 weight for nominal (or could use a default configured f value)
        const double w_nom = lookup_noac_weight(noac_weights_map[0.0f].get(), xi_gen_evt);
        DeltaY_xi_reco_6->Fill(deltay_xi, weight * w_nom);
        DeltaY_xi_reco_12->Fill(deltay_xi, weight * w_nom);
        DeltaY_xi_reco_18->Fill(deltay_xi, weight * w_nom);
        DeltaY_xi_reco_24->Fill(deltay_xi, weight * w_nom);
        DeltaY_xi_reco_30->Fill(deltay_xi, weight * w_nom);
        DeltaY_xi_reco_36->Fill(deltay_xi, weight * w_nom);
        DeltaY_xi_reco_50->Fill(deltay_xi, weight * w_nom);
      } else {
        // Fill without NoAC weights
        DeltaY_xi_reco_6->Fill(deltay_xi, weight);
        DeltaY_xi_reco_12->Fill(deltay_xi, weight);
        DeltaY_xi_reco_18->Fill(deltay_xi, weight);
        DeltaY_xi_reco_24->Fill(deltay_xi, weight);
        DeltaY_xi_reco_30->Fill(deltay_xi, weight);
        DeltaY_xi_reco_36->Fill(deltay_xi, weight);
        DeltaY_xi_reco_50->Fill(deltay_xi, weight);
      }

      // STEP 9: Fill NoAC systematics for all f values (matching Hists naming convention exactly)

      // RECO deltay_xi determines which bin is filled
      // GEN xi_gen_evt determines which weight is applied from the weight map
      // Weight w_f depends on the f value

      // Fill NoAC systematics for all f values (matching Hists naming convention exactly)
      if(use_noac_evtweights_ && !noac_weights_map.empty()){
        const int bin_schemes[] = {6, 12, 18, 24, 30, 36, 50};
        // Map f values to exact suffixes from kNoACSpecs (same as in init)
        static const std::map<float, std::string> f_to_suffix = {
          {-100.0f, "noacm100"}, {-12.0f, "noacm12"}, {-8.0f, "noacm8"}, {-4.0f, "noacm4"}, {-2.0f, "noacm2"},
          {-1.0f, "noacm1"}, {-0.8f, "noacm08"}, {-0.6f, "noacm06"}, {-0.4f, "noacm04"}, {-0.2f, "noacm02"},
          {0.0f, "noac0"},
          {0.2f, "noac02"}, {0.4f, "noac04"}, {0.6f, "noac06"}, {0.8f, "noac08"},
          {1.0f, "noac1"}, {2.0f, "noac2"}, {4.0f, "noac4"}, {8.0f, "noac8"}, {12.0f, "noac12"}, {100.0f, "noac100"}
        };
        for(const float fv : f_values){
          if(!noac_weights_map.count(fv)) continue;
          // STEP 8: Lookup the NoAC weight for the current f value
          const double w_f = lookup_noac_weight(noac_weights_map[fv].get(), xi_gen_evt);
          // Get exact suffix matching kNoACSpecs
          std::string suffix;
          auto it_suffix = f_to_suffix.find(fv);
          if(it_suffix != f_to_suffix.end()){
            suffix = it_suffix->second;
          } else {
            // Fallback (shouldn't happen if f_values matches kNoACSpecs)
            if(std::abs(fv) < 0.01f) suffix = "noac0";
            else if(fv < 0) suffix = "noacm" + std::to_string(static_cast<int>(std::abs(fv)));
            else suffix = "noac" + std::to_string(static_cast<int>(fv));
          }
          // Fill all binning schemes for this f value
          for(int nb : bin_schemes){
            std::string hname = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + suffix;
            auto it = h_deltaY_xi_reco_map.find(hname);
            if(it != h_deltaY_xi_reco_map.end()) it->second->Fill(deltay_xi, weight * w_f);
          }
        }
      }
      //template method end
      // up/down variations (multi-f map)
      // Get NoAC weight for f=0 (nominal) if available
      double w_noac_nominal = 1.0;
      if(use_noac_evtweights_ && !noac_weights_map.empty() && noac_weights_map.count(0.0f)){
        w_noac_nominal = lookup_noac_weight(noac_weights_map[0.0f].get(), xi_gen_evt);
      }
      for(unsigned int i=0; i<names.size(); i++){
          const double w_up = weight * syst_up.at(i)/syst_nominal.at(i);
          const double w_dn = weight * syst_down.at(i)/syst_nominal.at(i);
  
        const double nom = syst_nominal.at(i);
        const double up = syst_up.at(i);
        const double dn = syst_down.at(i);
        
        if(!std::isfinite(nom) || nom == 0.0 || !std::isfinite(up) || up == 0.0 || !std::isfinite(dn) || dn == 0.0){
          cout << "DEBUG SYSTEMATIC [" << names.at(i) << "]: nom=" << nom << ", up=" << up << ", dn=" << dn << endl;
        }

        hists_up_tt.at(i)->Fill(DeltaY_reco_best, DeltaY_gen_best, w_up);
        hists_down_tt.at(i)->Fill(DeltaY_reco_best, DeltaY_gen_best, w_dn);
        // xi systematics (with NoAC f=0 weight applied to match nominal)
        hists_deltaY_xi_reco_6_up.at(i)->Fill(deltay_xi, w_up * w_noac_nominal);
        hists_deltaY_xi_reco_6_down.at(i)->Fill(deltay_xi, w_dn * w_noac_nominal);
        // also fill 12/50 bin booked maps
        {
          std::string tag = names.at(i);
          auto fill_by_name = [&](const std::string &nm, double w){
            auto it = h_deltaY_xi_reco_map.find(nm);
            if(it != h_deltaY_xi_reco_map.end()) it->second->Fill(deltay_xi, w);
          };
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string prefix = "DeltaY_xi_reco_" + std::to_string(nb) + "_";
            fill_by_name(prefix + tag + "_up", w_up * w_noac_nominal);
            fill_by_name(prefix + tag + "_down", w_dn * w_noac_nominal);
          }
        }
      }
      if (debug)cout << "ttbar all sys fill deltay-systematics " <<endl;
      // scale variations
      for(unsigned int i=0; i<hists_scale.size(); i++){
        const double w_sc = weight * syst_scale.at(i);
        hists_scale_tt.at(i)->Fill(DeltaY_reco_best, DeltaY_gen_best, w_sc);
        hists_scale_deltaY_xi_reco_6.at(i)->Fill(deltay_xi, w_sc * w_noac_nominal);
        // name-based 12/50
        static const char* mmv[6] = {"upup","upnone","noneup","nonedown","downnone","downdown"};
        if(i<6){
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_murmuf_" + mmv[i];
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(deltay_xi, w_sc * w_noac_nominal);
          }
        }
      }
      if (debug)cout << "ttbar all sys fill deltay-scale variations " <<endl;
      // btag variations
      for(unsigned int i=0; i<hists_btag.size(); i++){
        const double w_bt = weight * (syst_btag.at(i)/btag_nominal);
        hists_btag_tt.at(i)->Fill(DeltaY_reco_best, DeltaY_gen_best, w_bt);
        hists_btag_deltaY_xi_reco_6.at(i)->Fill(deltay_xi, w_bt * w_noac_nominal);
        // derive tag name from order
        static const char* btags[] = {"btag_cferr1_up","btag_cferr1_down","btag_cferr2_up","btag_cferr2_down","btag_hf_up","btag_hf_down","btag_hfstats1_up","btag_hfstats1_down","btag_hfstats2_up","btag_hfstats2_down","btag_lf_up","btag_lf_down","btag_lfstats1_up","btag_lfstats1_down","btag_lfstats2_up","btag_lfstats2_down"};
        if(i<16){
          std::string tag = btags[i];
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + tag;
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(deltay_xi, w_bt * w_noac_nominal);
          }
        }
      }
      if (debug)cout << "ttbar all sys fill deltay - btag" <<endl;
      // ttag variations!
      for(unsigned int i=0; i<hists_ttag.size(); i++){
        const double w_tt = weight * (syst_ttag.at(i)/ttag_nominal);
        hists_ttag_tt.at(i)->Fill(DeltaY_reco_best, DeltaY_gen_best, w_tt);
        hists_ttag_deltaY_xi_reco_6.at(i)->Fill(deltay_xi, w_tt * w_noac_nominal);
        static const char* tags[] = {"ttag_corr_up","ttag_corr_down","ttag_uncorr_up","ttag_uncorr_down"};
        if(i<4){
          std::string t = tags[i];
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + t;
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(deltay_xi, w_tt * w_noac_nominal);
          }
        }
      }
      if (debug)cout << "ttbar all sys fill deltay- ttag " <<endl;
      // tmistag variations
      if (debug)cout <<"size for mistag: "<<hists_tmistag.size()<<endl;
      for(unsigned int i=0; i<hists_tmistag.size(); i++){
        const double w_tm = weight * (syst_tmistag.at(i)/tmistag_nominal);
        hists_tmistag_tt.at(i)->Fill(DeltaY_reco_best, DeltaY_gen_best, w_tm);
        hists_tmistag_deltaY_xi_reco_6.at(i)->Fill(deltay_xi, w_tm * w_noac_nominal);
        const char* t = (i==0? "tmistag_up" : "tmistag_down");
        const int map_bins[] = {12,18,24,30,36,50};
        for(int nb : map_bins){
          std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + t;
          auto it = h_deltaY_xi_reco_map.find(name);
          if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(deltay_xi, w_tm * w_noac_nominal);
        }
      }
      for(unsigned int i=0; i<hists_toppt.size(); i++){
        const double w_tp = weight * syst_toppt.at(i);
        hists_toppt_tt.at(i)->Fill(DeltaY_reco_best, DeltaY_gen_best, w_tp);
        hists_toppt_deltaY_xi_reco_6.at(i)->Fill(deltay_xi, w_tp * w_noac_nominal);
        static const char* tags[] = {"toppt_a_up","toppt_a_down","toppt_b_up","toppt_b_down"};
        if(i<4){
          std::string t = tags[i];
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + t;
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(deltay_xi, w_tp * w_noac_nominal);
          }
        }
      }
      if (debug)cout << "ttbar all sys fill deltay -mistag" <<endl;
      // ps variations
      for(unsigned int i=0; i<hists_ps.size(); i++){
        const double w_ps = weight * syst_ps.at(i);
        hists_ps_tt.at(i)->Fill(DeltaY_reco_best, DeltaY_gen_best, w_ps);
        hists_ps_deltaY_xi_reco_6.at(i)->Fill(deltay_xi, w_ps * w_noac_nominal);
        static const char* tags[] = {"isr_up","isr_down","fsr_up","fsr_down"};
        if(i<4){
          std::string t = tags[i];
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + t;
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(deltay_xi, w_ps * w_noac_nominal);
          }
        }
      }
      if (debug)cout << "done with ttbar RM deltay" <<endl;
    }
    else{
      // ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
      // float Mreco = BestZprimeCandidate->Zprime_v4().M();
      if (debug)cout << "in other MC cat"<<endl;
      if (debug)cout <<"one event loop"<<endl;
      float deltay=99.;
      bool isLeptonPositive = false;
      if (debug)cout << "ttbar all sys calc deltay " <<endl;
      if(isMuon){
        if (event.muons->at(0).charge() == 1){
          isLeptonPositive = true;
        } else {
          isLeptonPositive = false;
        }
      }
     
      if(isElectron){
        if (event.electrons->at(0).charge() == 1){
          isLeptonPositive = true;
        } else {
          isLeptonPositive = false;
        }
      }
      if (isLeptonPositive) {
        deltay=(TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()));
      }
      else {
        deltay=(TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity())); 
      }
      DeltaY->Fill(deltay, weight);
      // Template method xi for non-TT MC (no NoAC). Fill xi nominal only
      const double xi_reco = TMath::TanH(deltay);
      DeltaY_xi_reco_6->Fill(xi_reco, weight);
      DeltaY_xi_reco_12->Fill(xi_reco, weight);
      DeltaY_xi_reco_18->Fill(xi_reco, weight);
      DeltaY_xi_reco_24->Fill(xi_reco, weight);
      DeltaY_xi_reco_30->Fill(xi_reco, weight);
      DeltaY_xi_reco_36->Fill(xi_reco, weight);
      DeltaY_xi_reco_50->Fill(xi_reco, weight);
      if (debug)cout <<"fill nominal "<<endl;
      // Fill explicit xi systematics (no f dependence)
      for(unsigned int i=0; i<names.size(); i++){
        if (debug)cout <<"filling : "<<names[i]<<endl;
        hists_up.at(i)->Fill(deltay, weight * syst_up.at(i)/syst_nominal.at(i));
        hists_down.at(i)->Fill(deltay, weight * syst_down.at(i)/syst_nominal.at(i));
        // only 6-bin per-syst are booked
        hists_deltaY_xi_reco_6_up.at(i)->Fill(xi_reco, weight * (syst_up.at(i)/syst_nominal.at(i)));
        hists_deltaY_xi_reco_6_down.at(i)->Fill(xi_reco, weight * (syst_down.at(i)/syst_nominal.at(i)));
        // also fill 12/50 bin maps
        {
          std::string tag = names.at(i);
          auto fill_by_name = [&](const std::string &nm, double w){
            auto it = h_deltaY_xi_reco_map.find(nm);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(xi_reco, w);
          };
          const double w_up_val = weight * syst_up.at(i)/syst_nominal.at(i);
          const double w_down_val = weight * syst_down.at(i)/syst_nominal.at(i);
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string prefix = "DeltaY_xi_reco_" + std::to_string(nb) + "_";
            fill_by_name(prefix + tag + "_up", w_up_val);
            fill_by_name(prefix + tag + "_down", w_down_val);
          }
        }
        // cout << "filling xi systematics for " << names[i] <<" xi_reco: "<<xi_reco<< " weight: "<<weight * (syst_up.at(i)/syst_nominal.at(i))  << endl; 
      }
      // scale variations
      for(unsigned int i=0; i<hists_scale.size(); i++){
        hists_scale.at(i)->Fill(deltay, weight * syst_scale.at(i));
        hists_scale_deltaY_xi_reco_6.at(i)->Fill(xi_reco, weight * syst_scale.at(i));
        // name-based 12/50
        static const char* mmv[6] = {"upup","upnone","noneup","nonedown","downnone","downdown"};
        if(i<6){
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_murmuf_" + mmv[i];
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(xi_reco, weight * syst_scale.at(i));
          }
        }
      }
      // btag variations
      for(unsigned int i=0; i<hists_btag.size(); i++){
        hists_btag.at(i)->Fill(deltay, weight * syst_btag.at(i)/btag_nominal);
        hists_btag_deltaY_xi_reco_6.at(i)->Fill(xi_reco, weight * (syst_btag.at(i)/btag_nominal));
        // derive tag name from order
        static const char* btags[] = {"btag_cferr1_up","btag_cferr1_down","btag_cferr2_up","btag_cferr2_down","btag_hf_up","btag_hf_down","btag_hfstats1_up","btag_hfstats1_down","btag_hfstats2_up","btag_hfstats2_down","btag_lf_up","btag_lf_down","btag_lfstats1_up","btag_lfstats1_down","btag_lfstats2_up","btag_lfstats2_down"};
        if(i<16){
          std::string tag = btags[i];
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + tag;
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(xi_reco, weight * syst_btag.at(i)/btag_nominal);
          }
        }
      }
      // ttag variations
      for(unsigned int i=0; i<hists_ttag.size(); i++){
        hists_ttag.at(i)->Fill(deltay, weight * syst_ttag.at(i)/ttag_nominal);
        hists_ttag_deltaY_xi_reco_6.at(i)->Fill(xi_reco, weight * (syst_ttag.at(i)/ttag_nominal));
        static const char* tags[] = {"ttag_corr_up","ttag_corr_down","ttag_uncorr_up","ttag_uncorr_down"};
        if(i<4){
          std::string t = tags[i];
          const int map_bins[] = {12,18,24,30,36,50};
          for(int nb : map_bins){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + t;
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(xi_reco, weight * syst_ttag.at(i)/ttag_nominal);
          }
        }
      }
      // tmistag variations
      for(unsigned int i=0; i<hists_tmistag.size(); i++){
        hists_tmistag.at(i)->Fill(deltay, weight * syst_tmistag.at(i)/tmistag_nominal);
        hists_tmistag_deltaY_xi_reco_6.at(i)->Fill(xi_reco, weight * (syst_tmistag.at(i)/tmistag_nominal));
        const char* t = (i==0? "tmistag_up" : "tmistag_down");
        const int map_bins[] = {12,18,24,30,36,50};
        for(int nb : map_bins){
          std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + t;
          auto it = h_deltaY_xi_reco_map.find(name);
          if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(xi_reco, weight * syst_tmistag.at(i)/tmistag_nominal);
        }
      }
      // Top pt xi
      for(unsigned int i=0; i<hists_toppt.size(); i++){
        hists_toppt.at(i)->Fill(deltay, weight * syst_toppt.at(i));
        hists_toppt_deltaY_xi_reco_6.at(i)->Fill(xi_reco, weight * syst_toppt.at(i));
        static const char* tags[] = {"toppt_a_up","toppt_a_down","toppt_b_up","toppt_b_down"};
        if(i<4){
          std::string t = tags[i];
          const int map_bins_mm[] = {12,18,24,30,36,50};
          for(int nb : map_bins_mm){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + t;
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(xi_reco, weight * syst_toppt.at(i));
          }
        }
      }
      // isr fsr xi
      for(unsigned int i=0; i<hists_ps.size(); i++){
        hists_ps.at(i)->Fill(deltay, weight * syst_ps.at(i));
        hists_ps_deltaY_xi_reco_6.at(i)->Fill(xi_reco, weight * syst_ps.at(i));
        static const char* tags[] = {"isr_up","isr_down","fsr_up","fsr_down"};
        if(i<4){
          std::string t = tags[i];
          const int map_bins_mm[] = {12,18,24,30,36,50};
          for(int nb : map_bins_mm){
            std::string name = "DeltaY_xi_reco_" + std::to_string(nb) + "_" + t;
            auto it = h_deltaY_xi_reco_map.find(name);
            if(it!=h_deltaY_xi_reco_map.end()) it->second->Fill(xi_reco, weight * syst_ps.at(i));
          }
        }
      }
     
    }//end loop for all MC but ttbar for Deltay
    
    
  }//end loop for Deltay
  //all MC systematics for EFT
  if(is_zprime_reconstructed_chi2 && is_mc){
    if (debug)cout << "check ttbar & all MC sys for EFT reco variables" <<endl;
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);       // Best Z' candidate from chi2 reconstruction
    bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction();   // Reconstruction process id
    vector <Jet> AK4CHSjets_matched = event.get(h_CHSjets_matched);                    // AK4Puppijets that have been matched to CHSjets
    if (debug)cout << "define Jet for all sys for EFT reco variables" <<endl;          
    vector <TopJet> TopTaggedJets = event.get(h_AK8TopTags);                           // AK8Puppi jets TopTagged by DeepAK8TopTagger
    if (debug)cout << "define Top Jet for all sys for EFT reco variables" <<endl;                  
    vector <float> jets_hadronic_bscores;                                              // bScores vector for resolved hadronic jets
    float pt_hadTop_thresh = 150;                                                      // Define cut-variable as pt of hadTop for low/high regions
    float pt_hadTop = BestZprimeCandidate->top_hadronic_v4().pt();                     // pT of hadronic-top jet


    //-------------- Extracting highest b-tag score in Resolved topology, i.e. no top-tagged jet in event --------------//
    float bscore_max = -2;
    if(!is_toptag_reconstruction){
        // Loop over resolved hadronic jets to find their bscore via CHS jets
      for(unsigned int i=0; i<BestZprimeCandidate->jets_hadronic().size(); i++){
        double deltaR_min = 99;
        // Match resolved hadronic jets to CHS jets (which have bscores)
        for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
          double deltaR_CHS = deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(j));
          if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;}
        // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
        for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
          if(deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
          jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet());} // Using DeepJet btag score
      }
      // Loop over bScores-vector to extract highest bscore
      for(unsigned int i=0; i<jets_hadronic_bscores.size(); i++){
        float bscore = jets_hadronic_bscores.at(i);
        if(bscore > bscore_max) bscore_max = bscore;
      }
    }

    //--------------- Extracting highest b-tag score in Merged topology, i.e. with top-tagged jet in event ---------------//
    if(is_toptag_reconstruction){
        // Loop over hadronic top's subjets to extract highest bscore
      for(unsigned int i=0; i < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); i++){
        float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet(); // Using DeepJet btag score
        if(bscore > bscore_max) bscore_max = bscore;
      }
    }

    //------------------------------------Define 4vectors of top quarks and spin analyzers------------------------------------//
    TLorentzVector had_top_b(0, 0, 0, 0); // b-jet 4-vector
    // Resolved topology
    if(!is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK4-jet with highest bscore
      for(unsigned int i=0; i< BestZprimeCandidate->jets_hadronic().size(); i++){
        float bscore = jets_hadronic_bscores.at(i);
        if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->jets_hadronic().at(i).pt(), 
                                                        BestZprimeCandidate->jets_hadronic().at(i).eta(), 
                                                        BestZprimeCandidate->jets_hadronic().at(i).phi(), 
                                                        BestZprimeCandidate->jets_hadronic().at(i).energy());
      }
    }
    // Merged topology
    if(is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK8-subjet with highest bscore
      for(unsigned int j=0; j < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); j++){
        float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
        if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).pt(), 
                                                        BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).eta(), 
                                                        BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).phi(), 
                                                        BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).energy());
      }
    }

    TLorentzVector lep_top_lep(0, 0, 0, 0);  // Lepton 4-vector
    LorentzVector lep = BestZprimeCandidate->lepton().v4();
    lep_top_lep.SetPtEtaPhiE(lep.pt(), lep.eta(), lep.phi(), lep.E());

    TLorentzVector PosTop(0, 0, 0, 0); // top quark 4vector
    TLorentzVector NegTop(0, 0, 0, 0); // top antiquark 4vector
    if(BestZprimeCandidate->lepton().charge() > 0){ // Positively charged Lepton => Positively charged top QUARK mother
      PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                          BestZprimeCandidate->top_leptonic_v4().eta(), 
                          BestZprimeCandidate->top_leptonic_v4().phi(), 
                          BestZprimeCandidate->top_leptonic_v4().energy());
      NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                          BestZprimeCandidate->top_hadronic_v4().eta(), 
                          BestZprimeCandidate->top_hadronic_v4().phi(), 
                          BestZprimeCandidate->top_hadronic_v4().energy());
    }
    else if (BestZprimeCandidate->lepton().charge() < 0){ // Negatively charged Lepton => Positively charged top ANTIQUARK mother
    PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                          BestZprimeCandidate->top_hadronic_v4().eta(), 
                          BestZprimeCandidate->top_hadronic_v4().phi(), 
                          BestZprimeCandidate->top_hadronic_v4().energy());
    NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                          BestZprimeCandidate->top_leptonic_v4().eta(), 
                          BestZprimeCandidate->top_leptonic_v4().phi(), 
                          BestZprimeCandidate->top_leptonic_v4().energy());
    }
    TLorentzVector ttbar = PosTop + NegTop; // ttbar 4vector
    float beta = TMath::Abs(PosTop.Pz() + NegTop.Pz()) / (PosTop.E() + NegTop.E()); // relativistic beta of ttbar system in lab-frame
    float deltay = 0.0;
    if (BestZprimeCandidate->lepton().charge()>0) 
      {deltay = TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity());} 
    else 
      {deltay = TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity());}

    //--------- Boost into ttbar CoM-Frame ---------//
    // Center of Mass frame copies of 4vectors
    TLorentzVector lep_top_lep_CoM = lep_top_lep;
    TLorentzVector had_top_b_CoM = had_top_b;
    TLorentzVector PosTop_CoM = PosTop;
    TLorentzVector NegTop_CoM = NegTop;
    // Boost with negative of ttbar boost vector
    lep_top_lep_CoM.Boost(-1.*ttbar.BoostVector());
    had_top_b_CoM.Boost(-1.*ttbar.BoostVector());
    PosTop_CoM.Boost(-1.*ttbar.BoostVector());
    NegTop_CoM.Boost(-1.*ttbar.BoostVector());


    //-------------------------------------------------------------- Build Bernreuther basis --------------------------------------------------------------//
    // Required axes from CoM frame
    TVector3 beam_axis(0,0,1);                                                                      // Beam unit vector
    TVector3 k_axis = PosTop_CoM.Vect().Unit();                                                     // direction of top quark momentum in ttbar CoM frame
    double cos_PosTop_beam = PosTop_CoM.Vect().Unit().Dot(beam_axis);                               // Cosine of scattering angle, "y" in Bernreuther et al.
    double abs_sin_PosTop_beam = sqrt(1 - cos_PosTop_beam*cos_PosTop_beam);                         // Sine of scattering angle,   "r" in Bernreuther et al.
    TVector3 r_axis = ( (1./abs_sin_PosTop_beam) * (beam_axis - cos_PosTop_beam * k_axis) ).Unit(); // orthogonal to k_axis and lies in production plane
    TVector3 n_axis = ( (1./abs_sin_PosTop_beam) * beam_axis.Cross(k_axis) ).Unit();                // orthogonal to production plane
    double sign_cos_PosTop_beam = (cos_PosTop_beam > 0.) ? 1. : -1.;                                // Bose symmetry factor
    double sign_rapidity = (deltay > 0.) ? 1. : -1.;                                                // Charge asymmetry (in lab-frame) factor
  
    // Basis vectors
    TVector3 kbase = k_axis;
    TVector3 rbase = sign_cos_PosTop_beam * r_axis;
    TVector3 nbase = sign_cos_PosTop_beam * n_axis;
    // CA-corrected basis vectors
    TVector3 kStar = sign_rapidity * k_axis;
    TVector3 rStar = sign_rapidity * sign_cos_PosTop_beam * r_axis;
  
    
    //------------------- Boost spin-analyzers into top quark Rest-Frames -------------------//
    TLorentzVector lep_top_lep_Rest = lep_top_lep_CoM;
    TLorentzVector had_top_b_Rest = had_top_b_CoM;
    
    if(BestZprimeCandidate->lepton().charge() > 0){
      lep_top_lep_Rest.Boost(-1.*PosTop_CoM.BoostVector()); // lepton has Positive Top mother
      had_top_b_Rest.Boost(-1.*NegTop_CoM.BoostVector());   // b-jet has Negative Top mother
    }
    else if (BestZprimeCandidate->lepton().charge() < 0){
      lep_top_lep_Rest.Boost(-1.*NegTop_CoM.BoostVector()); // lepton has Negative Top mother
      had_top_b_Rest.Boost(-1.*PosTop_CoM.BoostVector());   // b-jet has Positive Top mother
    }


    //------------------------------------------- Spin Correlation variables -------------------------------------------//
    // float cosTheta1k_antiLep = 99.;
    // float cosTheta1r_antiLep = 99.;
    // float cosTheta1n_antiLep = 99.;
    // float cosTheta1kStar_antiLep = 99.;
    // float cosTheta1rStar_antiLep = 99.;
    float cosTheta1k = 99.;
    float cosTheta1r = 99.;
    float cosTheta1n = 99.;
    // float cosTheta1kStar = 99.;
    // float cosTheta1rStar = 99.;
    
    // float cosTheta2k_Lep = 99.;
    // float cosTheta2r_Lep = 99.;
    // float cosTheta2n_Lep = 99.;
    // float cosTheta2kStar_Lep = 99.;
    // float cosTheta2rStar_Lep = 99.;
    float cosTheta2k = 99.;
    float cosTheta2r = 99.;
    float cosTheta2n = 99.;
    // float cosTheta2kStar = 99.;
    // float cosTheta2rStar = 99.;

    // float C_kk = 99.;
    // float C_rr = 99.;
    // float C_nn = 99.;
    
    // float C_rk_plus = 99.;
    // float C_rk_minus = 99.;
    // float C_nr_plus = 99.;
    // float C_nr_minus = 99.;
    // float C_nk_plus = 99.;
    // float C_nk_minus = 99.;

    float CHel = 99.;
    float CHel_Mtt300_400 = 99.;
    float CHel_Mtt300_400_betaLT0p9 = 99.;

    float CHel_P3n = 99.;
    float CHel_P3n_Mtt800_Inf = 99.;
    float CHel_P3n_Mtt800_Inf_cosThetaLT0p4 = 99.;

    // // Use only leptons as spin-analyzers
    // if(BestZprimeCandidate->lepton().charge() > 0){
    //   // anti-lepton is spin-analyzer fo top quark
    //   cosTheta1k_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(kbase);
    //   cosTheta1r_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(rbase);
    //   cosTheta1n_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(nbase);
    //   cosTheta1kStar_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(kStar);
    //   cosTheta1rStar_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(rStar);
    // }
    // else if (BestZprimeCandidate->lepton().charge() < 0){
    //   // lepton is spin-analyzer for antitop quark
    //   cosTheta2k_Lep = lep_top_lep_Rest.Vect().Unit().Dot(kbase);
    //   cosTheta2r_Lep = lep_top_lep_Rest.Vect().Unit().Dot(rbase);
    //   cosTheta2n_Lep = lep_top_lep_Rest.Vect().Unit().Dot(nbase);
    //   cosTheta2kStar_Lep = lep_top_lep_Rest.Vect().Unit().Dot(kStar);
    //   cosTheta2rStar_Lep = lep_top_lep_Rest.Vect().Unit().Dot(rStar);
    // }

    // // Assign spin-analyzers depending on lepton charge
    // if(BestZprimeCandidate->lepton().charge() > 0){
    //   // top quark spin-analyzer is lepton
    //   cosTheta1k = lep_top_lep_Rest.Vect().Unit().Dot(kbase);
    //   cosTheta1r = lep_top_lep_Rest.Vect().Unit().Dot(rbase);
    //   cosTheta1n = lep_top_lep_Rest.Vect().Unit().Dot(nbase);
    //   cosTheta1kStar = lep_top_lep_Rest.Vect().Unit().Dot(kStar);
    //   cosTheta1rStar = lep_top_lep_Rest.Vect().Unit().Dot(rStar);
    //   // antitop spin-analyzer is b-jet
    //   cosTheta2k = had_top_b_Rest.Vect().Unit().Dot(kbase);
    //   cosTheta2r = had_top_b_Rest.Vect().Unit().Dot(rbase);
    //   cosTheta2n = had_top_b_Rest.Vect().Unit().Dot(nbase);
    //   cosTheta2kStar = had_top_b_Rest.Vect().Unit().Dot(kStar);
    //   cosTheta2rStar = had_top_b_Rest.Vect().Unit().Dot(rStar);
    // }
    // else if (BestZprimeCandidate->lepton().charge() < 0){
    //   // top quark spin-analyzer is b-jet
    //   cosTheta1k = had_top_b_Rest.Vect().Unit().Dot(kbase);
    //   cosTheta1r = had_top_b_Rest.Vect().Unit().Dot(rbase);
    //   cosTheta1n = had_top_b_Rest.Vect().Unit().Dot(nbase);
    //   cosTheta1kStar = had_top_b_Rest.Vect().Unit().Dot(kStar);
    //   cosTheta1rStar = had_top_b_Rest.Vect().Unit().Dot(rStar);
    //   // antitop spin-analyzer is lepton
    //   cosTheta2k = lep_top_lep_Rest.Vect().Unit().Dot(kbase);
    //   cosTheta2r = lep_top_lep_Rest.Vect().Unit().Dot(rbase);
    //   cosTheta2n = lep_top_lep_Rest.Vect().Unit().Dot(nbase);
    //   cosTheta2kStar = lep_top_lep_Rest.Vect().Unit().Dot(kStar);
    //   cosTheta2rStar = lep_top_lep_Rest.Vect().Unit().Dot(rStar);
    // }

    // // correlation matrix elements
    // C_nn = cosTheta1n * cosTheta2n;
    // C_rr = cosTheta1r * cosTheta2r;
    // C_kk = cosTheta1k * cosTheta2k;
    // // sum and differences of cross correlations
    // C_rk_plus = cosTheta1r * cosTheta2k + cosTheta1k * cosTheta2r;
    // C_rk_minus = cosTheta1r * cosTheta2k - cosTheta1k * cosTheta2r;
    // C_nr_plus = cosTheta1n * cosTheta2r + cosTheta1r * cosTheta2n;
    // C_nr_minus = cosTheta1n * cosTheta2r - cosTheta1r * cosTheta2n;
    // C_nk_plus = cosTheta1n * cosTheta2k + cosTheta1k * cosTheta2n;
    // C_nk_minus = cosTheta1n * cosTheta2k - cosTheta1k * cosTheta2n;

    // entanglement variables
    CHel = lep_top_lep_Rest.Vect().Unit().Dot(had_top_b_Rest.Vect().Unit());
    CHel_P3n = cosTheta1k * cosTheta2k + cosTheta1r * cosTheta2r - cosTheta1n * cosTheta2n;

    // near-threshold
    if(ttbar.M() < 400.){CHel_Mtt300_400 = CHel;
      // near-threshold and "slow"
      if(beta < 0.9){CHel_Mtt300_400_betaLT0p9 = CHel;}
    }

    // boosted
    if(ttbar.M() > 800.){CHel_P3n_Mtt800_Inf = CHel_P3n;
      // boosted and central
      if(TMath::Abs(cos_PosTop_beam) < 0.4){CHel_P3n_Mtt800_Inf_cosThetaLT0p4 = CHel_P3n;}
    }

    // Baumgart et al. angles depend on phi wrt Bernreuther nbase
    float lep_top_lep_phi = atan2(lep_top_lep_Rest.Vect().Unit().Dot(rbase), lep_top_lep_Rest.Vect().Unit().Dot(nbase));
    float had_top_b_phi   = atan2(had_top_b_Rest.Vect().Unit().Dot(rbase),   had_top_b_Rest.Vect().Unit().Dot(nbase));

    // sphi and dphi = PosTopDecayProd_phi +- NegTopDecayProd_phi
    float sphi = lep_top_lep_phi + had_top_b_phi;   // sum is independent of order
    float dphi = -99.;                              // initialize with dummy value
    if(BestZprimeCandidate->lepton().charge() > 0){ // lepton is Positive Top's Decay Product
      dphi = lep_top_lep_phi - had_top_b_phi;
    }
    if(BestZprimeCandidate->lepton().charge() < 0){ // b-jet is Positive Top's Decay Product
      dphi = had_top_b_phi - lep_top_lep_phi;
    }
    
    // Map back into original domain if necessary
    if(sphi > TMath::Pi()) sphi = sphi - 2*TMath::Pi();
    if(sphi < -TMath::Pi()) sphi = sphi + 2*TMath::Pi();
    if(dphi > TMath::Pi()) dphi = dphi - 2*TMath::Pi();
    if(dphi < -TMath::Pi()) dphi = dphi + 2*TMath::Pi();
    if (debug)cout << "spin corr vars defined"<<endl;


    //------------------------------------- Now fill histograms for systematics variations -------------------------------------//
    if (debug)cout <<" about to fill syst histograms" <<endl;

    if (debug)cout <<" how many sys: "<< names.size() <<endl;
    if (debug)cout <<" nominal up size: "<< hists_up.size() <<endl;
    if (debug)cout <<" sigma 1 up size: "<< hists_sigma_1_up.size() <<endl;
    if (debug)cout <<" sigma 2 up size: "<< hists_sigma_2_up.size() <<endl;
    if (debug)cout <<" dy 1 up size: "<< hists_dy_d1_up.size() <<endl;
    if (debug)cout <<" dy 2 up size: "<< hists_dy_d2_up.size() <<endl;

    if (debug)cout <<" nominal down size: "<< hists_down.size() <<endl;
    if (debug)cout <<" sigma 1 down size: "<< hists_sigma_1_down.size() <<endl;
    if (debug)cout <<" sigma 2 down size: "<< hists_sigma_2_down.size() <<endl;
    if (debug)cout <<" dy 1 down size: "<< hists_dy_d1_down.size() <<endl;
    if (debug)cout <<" dy 2 down size: "<< hists_dy_d2_down.size() <<endl;

    // Fill up/down variations ---------------------------------------------------------//
    for(unsigned int i=0; i<names.size(); i++){
      if (pt_hadTop > pt_hadTop_thresh && deltay>0){
        hists_sigma_1_up.at(i)->Fill(sphi, weight * syst_up.at(i)/syst_nominal.at(i));
        hists_sigma_1_down.at(i)->Fill(sphi, weight * syst_down.at(i)/syst_nominal.at(i));
      }
      if (debug)cout <<" done with sigma1" <<endl;
      if (pt_hadTop > pt_hadTop_thresh && deltay<0){
        hists_sigma_2_up.at(i)->Fill(sphi, weight * syst_up.at(i)/syst_nominal.at(i));
        hists_sigma_2_down.at(i)->Fill(sphi, weight * syst_down.at(i)/syst_nominal.at(i));
      }
      if (debug)cout <<" done with sigma2" <<endl;  
      if(pt_hadTop < pt_hadTop_thresh && dphi>0){
        hists_dy_d1_up.at(i)->Fill(deltay,weight * syst_up.at(i)/syst_nominal.at(i));
        hists_dy_d1_down.at(i)->Fill(deltay,weight * syst_down.at(i)/syst_nominal.at(i));
      }
      if (debug)cout <<" done with d1" <<endl;
      if(pt_hadTop < pt_hadTop_thresh && dphi<0){
        hists_dy_d2_up.at(i)->Fill(deltay,weight * syst_up.at(i)/syst_nominal.at(i));
        hists_dy_d2_down.at(i)->Fill(deltay,weight * syst_down.at(i)/syst_nominal.at(i));
      }
      if (debug)cout <<" done with d2" <<endl;
      // Entanglement variables
      // UP variations
      hists_cHel_Mtt300_400_up.at(i)->Fill(CHel_Mtt300_400, weight * syst_up.at(i)/syst_nominal.at(i));
      hists_cHel_Mtt300_400_betaLT0p9_up.at(i)->Fill(CHel_Mtt300_400_betaLT0p9, weight * syst_up.at(i)/syst_nominal.at(i));
      hists_cHel_P3n_Mtt800_Inf_up.at(i)->Fill(CHel_P3n_Mtt800_Inf, weight * syst_up.at(i)/syst_nominal.at(i));
      hists_cHel_P3n_Mtt800_Inf_cosThetaLT0p4_up.at(i)->Fill(CHel_P3n_Mtt800_Inf_cosThetaLT0p4, weight * syst_up.at(i)/syst_nominal.at(i));
      // DOWN variations
      hists_cHel_Mtt300_400_down.at(i)->Fill(CHel_Mtt300_400, weight * syst_down.at(i)/syst_nominal.at(i));
      hists_cHel_Mtt300_400_betaLT0p9_down.at(i)->Fill(CHel_Mtt300_400_betaLT0p9, weight * syst_down.at(i)/syst_nominal.at(i));
      hists_cHel_P3n_Mtt800_Inf_down.at(i)->Fill(CHel_P3n_Mtt800_Inf, weight * syst_down.at(i)/syst_nominal.at(i));
      hists_cHel_P3n_Mtt800_Inf_cosThetaLT0p4_down.at(i)->Fill(CHel_P3n_Mtt800_Inf_cosThetaLT0p4, weight * syst_down.at(i)/syst_nominal.at(i));
    }

    // btag variations -----------------------------------------------------------//
    if (debug)cout <<"how many in btag? "<<hists_btag.size() <<endl;
    if (debug)cout <<"sigma btag 1 "<<hists_btag_sigma_1.size() <<endl;
    if (debug)cout <<"sigma btag 2 "<<hists_btag_sigma_2.size() <<endl;
    if (debug)cout <<"dy btag 1 "<<hists_btag_dy_d1.size() <<endl;
    if (debug)cout <<"dy btag 2 "<<hists_btag_dy_d2.size() <<endl;
    for(unsigned int i=0; i<hists_btag.size(); i++){
      if (pt_hadTop > pt_hadTop_thresh && deltay>0){
        hists_btag_sigma_1.at(i)->Fill(sphi, weight * syst_btag.at(i)/btag_nominal);
      }
      if (debug)cout <<" done with s1 btag " << i << endl; 
      if (pt_hadTop > pt_hadTop_thresh && deltay<0){
        hists_btag_sigma_2.at(i)->Fill(sphi, weight * syst_btag.at(i)/btag_nominal);
      }
      if (debug)cout <<" done with s2 btag " << i << endl; 
      if(pt_hadTop < pt_hadTop_thresh && dphi>0){
        hists_btag_dy_d1.at(i)->Fill(deltay, weight * syst_btag.at(i)/btag_nominal);
      }
      if (debug)cout <<" done with d1 btag " <<i <<endl; 
      if(pt_hadTop < pt_hadTop_thresh && dphi<0){
        hists_btag_dy_d2.at(i)->Fill(deltay, weight * syst_btag.at(i)/btag_nominal);
      }
      if (debug)cout <<" done with d2 btag " <<i << endl;
      // Entanglement variables
      hists_btag_cHel_Mtt300_400.at(i)->Fill(CHel_Mtt300_400, weight * syst_btag.at(i)/btag_nominal);
      hists_btag_cHel_Mtt300_400_betaLT0p9.at(i)->Fill(CHel_Mtt300_400_betaLT0p9, weight * syst_btag.at(i)/btag_nominal);
      hists_btag_cHel_P3n_Mtt800_Inf.at(i)->Fill(CHel_P3n_Mtt800_Inf, weight * syst_btag.at(i)/btag_nominal);
      hists_btag_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.at(i)->Fill(CHel_P3n_Mtt800_Inf_cosThetaLT0p4, weight * syst_btag.at(i)/btag_nominal);
    }

    // ttag variations -------------------------------------------------------------//
    for(unsigned int i=0; i<hists_ttag.size(); i++){
      if (pt_hadTop > pt_hadTop_thresh && deltay>0){
        hists_ttag_sigma_1.at(i)->Fill(sphi, weight * syst_ttag.at(i)/ttag_nominal);
        }
      if (debug)cout <<" done with s1 ttag" <<endl; 
      if (pt_hadTop > pt_hadTop_thresh && deltay<0){
        hists_ttag_sigma_2.at(i)->Fill(sphi, weight * syst_ttag.at(i)/ttag_nominal);
        }
      if (debug)cout <<" done with s2 ttag" <<endl; 
      if(pt_hadTop < pt_hadTop_thresh && dphi>0){
        hists_ttag_dy_d1.at(i)->Fill(deltay, weight * syst_ttag.at(i)/ttag_nominal);
      }
      if (debug)cout <<" done with d1 ttag" <<endl; 
      if(pt_hadTop < pt_hadTop_thresh && dphi<0){
        hists_ttag_dy_d2.at(i)->Fill(deltay, weight * syst_ttag.at(i)/ttag_nominal);
      }
      if (debug)cout <<" done with d2 ttag" <<endl;
      // Entanglement variables
      hists_ttag_cHel_Mtt300_400.at(i)->Fill(CHel_Mtt300_400, weight * syst_ttag.at(i)/ttag_nominal);
      hists_ttag_cHel_Mtt300_400_betaLT0p9.at(i)->Fill(CHel_Mtt300_400_betaLT0p9, weight * syst_ttag.at(i)/ttag_nominal);
      hists_ttag_cHel_P3n_Mtt800_Inf.at(i)->Fill(CHel_P3n_Mtt800_Inf, weight * syst_ttag.at(i)/ttag_nominal);
      hists_ttag_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.at(i)->Fill(CHel_P3n_Mtt800_Inf_cosThetaLT0p4, weight * syst_ttag.at(i)/ttag_nominal);
    }

    // tmistag variations -------------------------------------------------------------------//
    for(unsigned int i=0; i<hists_tmistag.size(); i++){
      if (pt_hadTop > pt_hadTop_thresh && deltay>0){
        hists_tmistag_sigma_1.at(i)->Fill(sphi, weight * syst_tmistag.at(i)/tmistag_nominal);
        }
      if (debug)cout <<" done with s1 mistag" <<endl; 
      if (pt_hadTop > pt_hadTop_thresh && deltay<0){
        hists_tmistag_sigma_2.at(i)->Fill(sphi, weight * syst_tmistag.at(i)/tmistag_nominal);
        }
      if (debug)cout <<" done with s2 mistag" <<endl; 
      if(pt_hadTop < pt_hadTop_thresh && dphi>0){
        hists_tmistag_dy_d1.at(i)->Fill(deltay, weight * syst_tmistag.at(i)/tmistag_nominal);
      }
      if (debug)cout <<" done with d1 mistag" <<endl; 
      if(pt_hadTop < pt_hadTop_thresh && dphi<0){
        hists_tmistag_dy_d2.at(i)->Fill(deltay, weight * syst_tmistag.at(i)/tmistag_nominal);
      }
      if (debug)cout <<" done with d2 mistag" <<endl;
      // Entanglement variables
      hists_tmistag_cHel_Mtt300_400.at(i)->Fill(CHel_Mtt300_400, weight * syst_tmistag.at(i)/tmistag_nominal);
      hists_tmistag_cHel_Mtt300_400_betaLT0p9.at(i)->Fill(CHel_Mtt300_400_betaLT0p9, weight * syst_tmistag.at(i)/tmistag_nominal);
      hists_tmistag_cHel_P3n_Mtt800_Inf.at(i)->Fill(CHel_P3n_Mtt800_Inf, weight * syst_tmistag.at(i)/tmistag_nominal);
      hists_tmistag_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.at(i)->Fill(CHel_P3n_Mtt800_Inf_cosThetaLT0p4, weight * syst_tmistag.at(i)/tmistag_nominal);
    }

    // Top pt reweighting --------------------------------------------------//
    for(unsigned int i=0; i<hists_toppt.size(); i++){
        if (pt_hadTop > pt_hadTop_thresh && deltay>0){
          hists_toppt_sigma_1.at(i)->Fill(sphi, weight * syst_toppt.at(i));
        }
        if (debug)cout <<" done with s1 toppt" <<endl; 
        if (pt_hadTop > pt_hadTop_thresh && deltay<0){
          hists_toppt_sigma_2.at(i)->Fill(sphi, weight * syst_toppt.at(i));
        }
        if (debug)cout <<" done with s2 toppt" <<endl; 
        if(pt_hadTop < pt_hadTop_thresh && dphi>0){
          hists_toppt_dy_d1.at(i)->Fill(deltay, weight * syst_toppt.at(i));
        }
        if (debug)cout <<" done with d1 toppt" <<endl; 
        if(pt_hadTop < pt_hadTop_thresh && dphi<0){
          hists_toppt_dy_d2.at(i)->Fill(deltay, weight * syst_toppt.at(i));
        }
        if (debug)cout <<" done with d2 toppt" <<endl;
        // Entanglement variables
        hists_toppt_cHel_Mtt300_400.at(i)->Fill(CHel_Mtt300_400, weight * syst_toppt.at(i));
        hists_toppt_cHel_Mtt300_400_betaLT0p9.at(i)->Fill(CHel_Mtt300_400_betaLT0p9, weight * syst_toppt.at(i));
        hists_toppt_cHel_P3n_Mtt800_Inf.at(i)->Fill(CHel_P3n_Mtt800_Inf, weight * syst_toppt.at(i));
        hists_toppt_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.at(i)->Fill(CHel_P3n_Mtt800_Inf_cosThetaLT0p4, weight * syst_toppt.at(i));
      }

    // "scale", i.e. muR and muF, variations ------------------------------------------------//
    if (debug)cout <<"how many in scale? "<<hists_scale.size() <<endl;
    if (debug)cout <<"sigma 1 "<<hists_scale_sigma_1.size() <<endl;
    if (debug)cout <<"sigma 2 "<<hists_scale_sigma_2.size() <<endl;
    if (debug)cout <<"dy 1 "<<hists_scale_dy_d1.size() <<endl;
    if (debug)cout <<"dy 2 "<<hists_scale_dy_d2.size() <<endl;
    for(unsigned int i=0; i<hists_scale.size(); i++){
      if (pt_hadTop > pt_hadTop_thresh && deltay>0){
        hists_scale_sigma_1.at(i)->Fill(sphi, weight * syst_scale.at(i));
        }
      if (debug)cout <<" done with s1" <<endl; 
      if (pt_hadTop > pt_hadTop_thresh && deltay<0){
        hists_scale_sigma_2.at(i)->Fill(sphi, weight * syst_scale.at(i));
        }
      if (debug)cout <<" done with s2" <<endl; 
      if(pt_hadTop < pt_hadTop_thresh && dphi>0){
        hists_scale_dy_d1.at(i)->Fill(deltay, weight * syst_scale.at(i));
      }
      if (debug)cout <<" done with d1" <<endl; 
      if(pt_hadTop < pt_hadTop_thresh && dphi<0){
        hists_scale_dy_d2.at(i)->Fill(deltay, weight * syst_scale.at(i));
      }
      if (debug)cout <<" done with d2" <<endl;
      // Entanglement variables
      hists_scale_cHel_Mtt300_400.at(i)->Fill(CHel_Mtt300_400, weight * syst_scale.at(i));
      hists_scale_cHel_Mtt300_400_betaLT0p9.at(i)->Fill(CHel_Mtt300_400_betaLT0p9, weight * syst_scale.at(i));
      hists_scale_cHel_P3n_Mtt800_Inf.at(i)->Fill(CHel_P3n_Mtt800_Inf, weight * syst_scale.at(i));
      hists_scale_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.at(i)->Fill(CHel_P3n_Mtt800_Inf_cosThetaLT0p4, weight * syst_scale.at(i));
    }
    
    // ps variations ------------------------------------------------//
    for(unsigned int i=0; i<hists_ps.size(); i++){
      if (pt_hadTop > pt_hadTop_thresh && deltay>0){
          hists_ps_sigma_1.at(i)->Fill(sphi, weight * syst_ps.at(i));
      }
      if (debug)cout <<" done with s1 ps" <<endl;   
      if (pt_hadTop > pt_hadTop_thresh && deltay<0){
        hists_ps_sigma_2.at(i)->Fill(sphi, weight * syst_ps.at(i));
      }
      if (debug)cout <<" done with s2 ps" <<endl;  
      if(pt_hadTop < pt_hadTop_thresh && dphi>0){
        hists_ps_dy_d1.at(i)->Fill(deltay, weight * syst_ps.at(i));
      }
      if (debug)cout <<" done with d1 ps" <<endl;  
      if(pt_hadTop < pt_hadTop_thresh && dphi<0){
        hists_ps_dy_d2.at(i)->Fill(deltay, weight * syst_ps.at(i));
      }
      if (debug)cout <<" done with d2 ps" <<endl;
      // Entanglement variables
      hists_ps_cHel_Mtt300_400.at(i)->Fill(CHel_Mtt300_400, weight * syst_ps.at(i));
      hists_ps_cHel_Mtt300_400_betaLT0p9.at(i)->Fill(CHel_Mtt300_400_betaLT0p9, weight * syst_ps.at(i));
      hists_ps_cHel_P3n_Mtt800_Inf.at(i)->Fill(CHel_P3n_Mtt800_Inf, weight * syst_ps.at(i));
      hists_ps_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.at(i)->Fill(CHel_P3n_Mtt800_Inf_cosThetaLT0p4, weight * syst_ps.at(i));
    }

  }//end of all MC systematics for EFT  

}

ZprimeSemiLeptonicSystematicsHists::~ZprimeSemiLeptonicSystematicsHists(){}