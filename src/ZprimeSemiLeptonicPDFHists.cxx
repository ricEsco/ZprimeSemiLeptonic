#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicPDFHists.h"
#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>
#include <sstream>
#include <iomanip>

#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <UHH2/core/include/LorentzVector.h>
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <glob.h>
#include <cstring>

using namespace std;
using namespace uhh2;

// NoAC weight helper functions (static methods)
std::unique_ptr<TH1D> ZprimeSemiLeptonicPDFHists::mirror_hist_1d(const TH1* H) {
  auto H1 = dynamic_cast<const TH1D*>(H);
  if(!H1) throw std::runtime_error("[NoAC] mirror_hist_1d expects TH1D");
  auto M = std::unique_ptr<TH1D>(static_cast<TH1D*>(H1->Clone((std::string(H1->GetName())+"_mir").c_str())));
  M->SetDirectory(nullptr);
  M->Reset("ICES");
  const TAxis* xax = H1->GetXaxis();
  for(int i=1;i<=H1->GetNbinsX();++i){
    const double xc = xax->GetBinCenter(i);
    int j = xax->FindBin(-xc);
    if(j < 1) j = 1;
    if(j > H1->GetNbinsX()) j = H1->GetNbinsX();
    M->SetBinContent(i, H1->GetBinContent(j));
    M->SetBinError  (i, H1->GetBinError(j));
  }
  return M;
}

std::unique_ptr<TH1D> ZprimeSemiLeptonicPDFHists::build_noac_weights_from_gen(const TH1D& Hgen_in, float f_noac) {
  auto Hgen = &Hgen_in;
  if(!Hgen) throw std::runtime_error("[NoAC] Hgen is null");

  auto Hmir = mirror_hist_1d(Hgen);
  auto S = std::unique_ptr<TH1D>(static_cast<TH1D*>(Hgen->Clone("NoAC_S"))); S->SetDirectory(nullptr); S->Reset("ICES");
  auto A = std::unique_ptr<TH1D>(static_cast<TH1D*>(Hgen->Clone("NoAC_A"))); A->SetDirectory(nullptr); A->Reset("ICES");

  S->Add(Hgen, Hmir.get(), 0.5,  0.5);
  A->Add(Hgen, Hmir.get(), 0.5, -0.5);

  auto W = std::unique_ptr<TH1D>(static_cast<TH1D*>(Hgen->Clone("NoAC_W")));
  W->SetDirectory(nullptr); W->Reset("ICES");

  const int nb = Hgen->GetNbinsX();
  int n_bad = 0;
  for(int i=1;i<=nb;++i){
    const double s = S->GetBinContent(i);
    const double a = A->GetBinContent(i);
    const double denom = s + a;
    const double numer = s + (1.0 - static_cast<double>(f_noac)) * a;
    double w = (denom>0.0 ? numer/denom : 1.0);
    if(!(denom>0.0)) ++n_bad;
    W->SetBinContent(i, w);
  }
  if(n_bad){
    std::cout << "[NoAC] warning: " << n_bad << " GEN bins had (S+A)<=0; set w=1 there" << std::endl;
  }

  // renormalize <W> wrt Hgen content to keep total yield
  double sumW=0.0, sum1=0.0;
  for(int i=1;i<=nb;++i){
    const double wi = W->GetBinContent(i);
    const double hi = Hgen->GetBinContent(i);
    sumW += wi*hi;
    sum1 += hi;
  }
  if(sumW>0.0 && sum1>0.0){
    const double norm = sumW/sum1;
    for(int i=1;i<=nb;++i){
      W->SetBinContent(i, W->GetBinContent(i)/norm);
    }
  }
  return W;
}

double ZprimeSemiLeptonicPDFHists::lookup_noac_weight(double xi, const TH1* W){
  if(!W || !std::isfinite(xi)) return 1.0;
  double xmin = W->GetXaxis()->GetXmin();
  double xmax = W->GetXaxis()->GetXmax();
  if(xi <= xmin) xi = std::nextafter(xmin, xmax);
  if(xi >= xmax) xi = std::nextafter(xmax, xmin);
  int bin = W->GetXaxis()->FindFixBin(xi);
  if(bin < 1) bin = 1;
  if(bin > W->GetNbinsX()) bin = W->GetNbinsX();
  const double w = W->GetBinContent(bin);
  if(!std::isfinite(w) || w <= 0.0 || w > 100.0) return 1.0;
  return w;
}

namespace {
  // Helper: clone arbitrary TH1 into a detached TH1D
  static std::unique_ptr<TH1D> clone_as_TH1D(const TH1* src, const std::string &out_name){
    if(!src) return nullptr;
    const TAxis* ax = src->GetXaxis();
    const int nb = ax->GetNbins();
    std::unique_ptr<TH1D> dst;
    const TArrayD* xbins = ax->GetXbins();
    if(xbins && xbins->GetSize() > 0){
      dst.reset(new TH1D(out_name.c_str(), src->GetTitle(), nb, xbins->GetArray()));
    } else {
      dst.reset(new TH1D(out_name.c_str(), src->GetTitle(), nb, ax->GetXmin(), ax->GetXmax()));
    }
    dst->SetDirectory(nullptr);
    for(int i=1;i<=nb;++i){
      dst->SetBinContent(i, src->GetBinContent(i));
      dst->SetBinError  (i, src->GetBinError(i));
    }
    return dst;
  }
  
  // Map f values to exact suffixes from kNoACSpecs
  std::string get_noac_suffix(float fv) {
    static const std::map<float, std::string> f_to_suffix = {
      {-100.0f, "noacm100"}, {-12.0f, "noacm12"}, {-8.0f, "noacm8"}, {-4.0f, "noacm4"}, {-2.0f, "noacm2"},
      {-1.0f, "noacm1"}, {-0.8f, "noacm08"}, {-0.6f, "noacm06"}, {-0.4f, "noacm04"}, {-0.2f, "noacm02"},
      {0.0f, "noac0"},
      {0.2f, "noac02"}, {0.4f, "noac04"}, {0.6f, "noac06"}, {0.8f, "noac08"},
      {1.0f, "noac1"}, {2.0f, "noac2"}, {4.0f, "noac4"}, {8.0f, "noac8"}, {12.0f, "noac12"}, {100.0f, "noac100"}
    };
    auto it = f_to_suffix.find(fv);
    if(it != f_to_suffix.end()) return it->second;
    if(std::abs(fv) < 0.01f) return "noac0";
    else if(fv < 0) return "noacm" + std::to_string(static_cast<int>(std::abs(fv)));
    else return "noac" + std::to_string(static_cast<int>(fv));
  }
}

ZprimeSemiLeptonicPDFHists::ZprimeSemiLeptonicPDFHists(uhh2::Context & ctx, const std::string& dirname): 
Hists(ctx, dirname){
  // dataset type flags
  is_mc = ctx.get("dataset_type") == "MC";
  is_dy = ctx.get("dataset_version").find("DYJets") == 0;
  std::string dataset_version = ctx.get("dataset_version");
  is_tt = (dataset_version.find("TTTo") == 0) || (dataset_version.find("EFT") != std::string::npos);
  is_wjets = ctx.get("dataset_version").find("WJets") == 0;
  is_qcd_HTbinned = ctx.get("dataset_version").find("QCD_HT") == 0;
  is_alps = ctx.get("dataset_version").find("ALP") == 0;
  is_azh = ctx.get("dataset_version").find("AZH") == 0;
  is_htott_scalar = ctx.get("dataset_version").find("HscalarToTTTo") == 0;
  is_htott_pseudo = ctx.get("dataset_version").find("HpseudoToTTTo") == 0;
  is_zprimetott = ctx.get("dataset_version").find("ZPrimeToTT_") == 0;

  // handle for best reconstructed candidate and reconstruction flag
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");

  // top tagging algorithm flag
  ishotvr = (ctx.get("is_hotvr") == "true");
  isdeepAK8 = (ctx.get("is_deepAK8") == "true");
 
  // lepton channel flag
  isMuon = false; isElectron = false;
  if(ctx.get("channel") == "muon") isMuon = true;
  if(ctx.get("channel") == "electron") isElectron = true;
  if(isdeepAK8){
    h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("DeepAK8TopTags");
  }else if(ishotvr){
    h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("HOTVRTopTags");
  }
  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");
  
  //template method - NoAC setup for TTbar
  use_noac_evtweights_ = false;
  noac_gen_file_ = "";
  noac_gen_hist_ = "";
  if(is_mc && is_tt){
    h_xi_gen = ctx.get_handle<float>("xi_gen");
    use_noac_evtweights_ = (ctx.get("noac_apply_event_weight") == string("true"));
    noac_gen_file_ = ctx.get("noac_gen_file");
    noac_gen_hist_ = ctx.get("noac_gen_hist");
    
    // f values for the NoAC weights (all f values from kNoACSpecs)
    f_values = {-100.0f, -12.0f, -8.0f, -4.0f, -2.0f, -1.0f, -0.8f, -0.6f, -0.4f, -0.2f, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 2.0f, 4.0f, 8.0f, 12.0f, 100.0f};
    
    if(use_noac_evtweights_ && !noac_gen_file_.empty() && !noac_gen_hist_.empty()){
      // Glob inputs to get the generator histogram files
      glob_t gl; memset(&gl, 0, sizeof(gl));
      int r = glob(noac_gen_file_.c_str(), 0, nullptr, &gl);
      // sum the generator histograms
      TH1D *sumH = nullptr;
      if(r == 0){
        // loop over the generator histogram files
        for(size_t i=0;i<gl.gl_pathc;++i){
          const char *fp = gl.gl_pathv[i];
          std::unique_ptr<TFile> f(TFile::Open(fp));
          if(!f || f->IsZombie()) continue;
          TH1 *h = dynamic_cast<TH1*>(f->Get(noac_gen_hist_.c_str()));
          if(!h) continue;
          std::unique_ptr<TH1D> hD = clone_as_TH1D(h, "_tmpH");
          if(!sumH){ sumH = static_cast<TH1D*>(hD->Clone("Hgen_sum_pdf")); sumH->SetDirectory(0); }
          else { sumH->Add(hD.get()); }
        }
        globfree(&gl);
      }
      if(sumH){
        // build the NoAC weights from the generator histogram
        noac_weights_map.clear();
        for(const float fv : f_values){
          try{ 
            noac_weights_map[fv] = build_noac_weights_from_gen(*sumH, fv); 
          } catch(...){ }
        }
        delete sumH;
      }
    }
  }
  //template method end

  for(int i=0; i<100; i++){
    std::stringstream ss_name;
    std::stringstream ss_name_tt;
    std::stringstream ss_name_dy_d1;
    std::stringstream ss_name_dy_d2;
    std::stringstream ss_name_sigma_1;
    std::stringstream ss_name_sigma_2;
    std::stringstream ss_name_delta_1;
    std::stringstream ss_name_delta_2;
    std::stringstream ss_name_xi;


    ss_name         << "DeltaY_PDF_" << i+1;
    ss_name_tt      << "DeltaY_PDF_RM_" << i+1;
    ss_name_dy_d1   << "DeltaY_reco_d1_PDF_" << i+1;
    ss_name_dy_d2   << "DeltaY_reco_d2_PDF_" << i+1;
    ss_name_sigma_1 << "Sigma_phi_1_PDF_" << i+1;
    ss_name_sigma_2 << "Sigma_phi_2_PDF_" << i+1;
    ss_name_delta_1 << "Delta_phi_1_PDF_" << i+1;
    ss_name_delta_2 << "Delta_phi_2_PDF_" << i+1;
    ss_name_xi      << "DeltaY_xi_reco_6_PDF_" << i+1;


    stringstream ss_title;
    stringstream ss_title_tt;
    stringstream ss_title_dy_d1;
    stringstream ss_title_dy_d2;
    stringstream ss_title_sigma_1;
    stringstream ss_title_sigma_2;
    stringstream ss_title_xi;



    ss_title << "#DeltaY_{t#bar{t}} for PDF No. "  << i+1 << " out of 100" ;
    ss_title_tt << "#DeltaY_{t#bar{t}} RM for PDF No. "  << i+1 << " out of 100" ;

    ss_title_dy_d1 <<"#DeltaY_{t#bar{t}} for #Delta #phi >0 for PDF No. "<< i+1 << " out of 100" ;
    ss_title_dy_d2 <<"#DeltaY_{t#bar{t}} for #Delta #phi <0 for PDF No. "<< i+1 << " out of 100" ;

    ss_title_sigma_1 <<"#Sigma #phi for #DeltaY >0 for PDF No. "<< i+1 << " out of 100" ;
    ss_title_sigma_2 <<"#Sigma #phi for #DeltaY <0 for PDF No. "<< i+1 << " out of 100" ;

    ss_title_xi    << "tanh(#Delta y)_{reco} for PDF No. " << i+1 << " out of 100";

    std::string s_name = ss_name.str();
    std::string s_name_dy_d1 = ss_name_dy_d1.str();
    std::string s_title_dy_d1 = ss_title_dy_d1.str();
    std::string s_name_dy_d2 = ss_name_dy_d2.str();
    std::string s_title_dy_d2 = ss_title_dy_d2.str();
    std::string s_name_sigma_1 = ss_name_sigma_1.str();
    std::string s_title_sigma_1 = ss_title_sigma_1.str();
    std::string s_name_sigma_2 = ss_name_sigma_2.str();
    std::string s_title_sigma_2 = ss_title_sigma_2.str();

    std::string s_name_xi    = ss_name_xi.str();
    std::string s_title_xi    = ss_title_xi.str();


    std::string s_name_tt = ss_name_tt.str();
    std::string s_title = ss_title.str();
    std::string s_title_tt = ss_title_tt.str();
    const char* char_name = s_name.c_str();
    const char* char_name_dy_d1 = s_name_dy_d1.c_str();
    const char* char_name_dy_d2 = s_name_dy_d2.c_str();
    const char* char_name_sigma_1 = s_name_sigma_1.c_str();
    const char* char_name_sigma_2 = s_name_sigma_2.c_str();
    const char* char_name_tt = s_name_tt.c_str();
    const char* char_title = s_title.c_str();
    const char* char_title_tt = s_title_tt.c_str();
    const char* char_title_dy_d1 = s_title_dy_d1.c_str();
    const char* char_title_dy_d2 = s_title_dy_d2.c_str();
    const char* char_title_sigma_1 = s_title_sigma_1.c_str();
    const char* char_title_sigma_2 = s_title_sigma_2.c_str();


 
    hist_names[i] = s_name;
    hist_names_tt[i] = s_name_tt;
    hist_names_dy_d1[i] = s_name_dy_d1;
    hist_names_dy_d2[i] = s_name_dy_d2;
    hist_names_sigma_1[i] = s_name_sigma_1;
    hist_names_sigma_2[i] = s_name_sigma_2;
    hist_names_xi[i]    = s_name_xi;


    book<TH1F>(char_name_dy_d1, char_title_dy_d1,  2, -2.5, 2.5);
    book<TH1F>(char_name_dy_d2, char_title_dy_d2,  2, -2.5, 2.5);
    book<TH1F>(char_name_sigma_1, char_title_sigma_1,  16, -3.2, 3.2);
    book<TH1F>(char_name_sigma_2, char_title_sigma_2,  16, -3.2, 3.2);
    book<TH1F>(char_name, char_title,  2, -2.5, 2.5);
    book<TH2F>(char_name_tt, char_title_tt,  2, -2.5, 2.5, 2, -2.5, 2.5);

    book<TH1F>(s_name_xi.c_str(),    s_title_xi.c_str(),    /*nbins*/ 6,  -1.0,  1.0);
  }
  
  // Book PDF histograms for all binning schemes (6, 12, 18, 24, 30, 36, 50) - nominal only
  // Note: 6-bin histograms are already booked above, so add them to the map
  for(int pdf_idx = 0; pdf_idx < 100; ++pdf_idx){
    std::stringstream ss6;
    ss6 << "DeltaY_xi_reco_6_PDF_" << (pdf_idx + 1);
    std::string hname6 = ss6.str();
    if(h_pdf_xi_reco_map.find(hname6) == h_pdf_xi_reco_map.end()){
      // The 6-bin histograms are already booked as hist_names_xi, so get them from there
      if(auto* hxi = dynamic_cast<TH1F*>(hist(hist_names_xi[pdf_idx].c_str()))){
        h_pdf_xi_reco_map[hname6] = hxi;
      }
    }
  }
  
  const int pdf_bin_schemes[] = {12, 18, 24, 30, 36, 50};
  for(int pdf_idx = 0; pdf_idx < 100; ++pdf_idx){
    for(int nb : pdf_bin_schemes){
      std::stringstream ss;
      ss << "DeltaY_xi_reco_" << nb << "_PDF_" << (pdf_idx + 1);
      std::string hname = ss.str();
      if(h_pdf_xi_reco_map.find(hname) == h_pdf_xi_reco_map.end()){
        std::stringstream ss_title;
        ss_title << "tanh(#Delta y)_{reco} " << nb << " bins for PDF No. " << (pdf_idx + 1) << " out of 100";
        h_pdf_xi_reco_map[hname] = book<TH1F>(hname.c_str(), ss_title.str().c_str(), nb, -1.0, 1.0);
      }
    }
  }
}

void ZprimeSemiLeptonicPDFHists::fill(const Event & event){

  double weight = event.weight;
  bool debug=false;
  bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2);
  if(is_zprime_reconstructed_chi2 && is_mc){
    if(is_tt){
      if (debug)cout<<" doing ttbar sample" <<endl;

      ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
      // float Mreco = BestZprimeCandidate->Zprime_v4().M();
      const auto& genparticles = event.genparticles;
      
      GenParticle top, antitop;
      for(const GenParticle & gp : *genparticles){
          if(gp.pdgId() == 6){
              top = gp;
          }
          else if(gp.pdgId() == -6){
              antitop = gp;
          }
      }
      // The Lorentz vectors represent the 4-momenta (energy, and three spatial momentum components) for the leptonic and hadronic tops from the "BestZprimeCandidate" object
      LorentzVector lep_top = BestZprimeCandidate->top_leptonic_v4();
      LorentzVector had_top = BestZprimeCandidate->top_hadronic_v4();

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

      // deltaY values calculation starts:

      // matched gen particles
      for (const auto& pair_had : deltaR_hadronic_values) {
        if (pair_had.first > 0 && pair_had.first < deltaR_min_hadronic) {
          // deltaR_min_hadronic = pair_had.first;
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
        }

        if (isLeptonPositive) {
        //if (debug)cout <<"Lepton is positive for PDF"<<endl;
        DeltaY_reco_best = TMath::Abs(0.5*TMath::Log((lep_top.energy() + lep_top.pt()*TMath::SinH(lep_top.eta()))/(lep_top.energy() - lep_top.pt()*TMath::SinH(lep_top.eta())))) - TMath::Abs(0.5*TMath::Log((had_top.energy() + had_top.pt()*TMath::SinH(had_top.eta()))/(had_top.energy() - had_top.pt()*TMath::SinH(had_top.eta()))));
        DeltaY_gen_best = TMath::Abs(0.5*TMath::Log((best_matched_gen_leptop.energy() + best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta()))/(best_matched_gen_leptop.energy() - best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta())))) - TMath::Abs(0.5*TMath::Log((best_matched_gen_hadtop.energy() + best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta()))/(best_matched_gen_hadtop.energy() - best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta()))));
        } else {
          //if (debug)cout <<"Lepton is negative for PDF"<<endl;
          DeltaY_reco_best = TMath::Abs(0.5*TMath::Log((had_top.energy() + had_top.pt()*TMath::SinH(had_top.eta()))/(had_top.energy() - had_top.pt()*TMath::SinH(had_top.eta())))) - TMath::Abs(0.5*TMath::Log((lep_top.energy() + lep_top.pt()*TMath::SinH(lep_top.eta()))/(lep_top.energy() - lep_top.pt()*TMath::SinH(lep_top.eta()))));
          DeltaY_gen_best = TMath::Abs(0.5*TMath::Log((best_matched_gen_hadtop.energy() + best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta()))/(best_matched_gen_hadtop.energy() - best_matched_gen_hadtop.pt()*TMath::SinH(best_matched_gen_hadtop.eta())))) - TMath::Abs(0.5*TMath::Log((best_matched_gen_leptop.energy() + best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta()))/(best_matched_gen_leptop.energy() - best_matched_gen_leptop.pt()*TMath::SinH(best_matched_gen_leptop.eta()))));

        }
      }
      
      
      

      // Fill the histogram

      int MY_FIRST_INDEX = 9;
      if(event.genInfo->systweights().size() > (unsigned int) 100 + MY_FIRST_INDEX){
        float orig_weight = event.genInfo->originalXWGTUP();
        for(int i=0; i<100; i++){
          double pdf_weight = event.genInfo->systweights().at(i + MY_FIRST_INDEX);
          const char* name_tt = hist_names_tt[i].c_str();
          TH1* base_hist = hist(name_tt); // Retrieve histogram using base class method
          TH2F* hist2d = dynamic_cast<TH2F*>(base_hist); // Attempt to cast to TH2F

          if(hist2d) {
              hist2d->Fill(DeltaY_reco_best, DeltaY_gen_best, weight * pdf_weight / orig_weight);
          } else {
              // Handle the error if the histogram is not of type TH2F
              std::cerr << "Histogram casting error for: " << name_tt << std::endl;
          }
        }
      }
      
      // template-method variable
      // xi = tanh(dy_reco)
      double deltay_reco = 0.0;
      if (isLeptonPositive) {
        deltay_reco = std::abs(lep_top.Rapidity()) - std::abs(had_top.Rapidity());
      } else {
        deltay_reco = std::abs(had_top.Rapidity()) - std::abs(lep_top.Rapidity());
      }
      const double deltay_xi = TMath::TanH(deltay_reco);

      // ---- guards (early return) ----
      const auto &systw = event.genInfo->systweights();
      const size_t needed = static_cast<size_t>(MY_FIRST_INDEX) + 100u;

      const float orig_w = event.genInfo->originalXWGTUP();
      if (systw.size() > needed && orig_w != 0.f) {
        // Get NoAC weight for f=0 (nominal) if available
        double w_noac_nominal = 1.0;
        if(use_noac_evtweights_ && !noac_weights_map.empty() && noac_weights_map.count(0.0f) && event.is_valid(h_xi_gen)){
          const double xi_gen_evt = static_cast<double>(event.get(h_xi_gen));
          w_noac_nominal = lookup_noac_weight(xi_gen_evt, noac_weights_map[0.0f].get());
        }

        // ---- fill 100 PDF replica histos for 6-bin (existing) ----
        for (int i = 0; i < 100; ++i) {
          const double pdfw = systw.at(MY_FIRST_INDEX + i);
          if (auto* hxi = dynamic_cast<TH1F*>(hist(hist_names_xi[i].c_str())))
            hxi->Fill(deltay_xi, weight * pdfw / orig_w * w_noac_nominal);
        }
        
        // ---- fill 100 PDF replica histos for all other binning schemes (12, 18, 24, 30, 36, 50) ----
        const int pdf_bin_schemes[] = {12, 18, 24, 30, 36, 50};
        for (int i = 0; i < 100; ++i) {
          const double pdfw = systw.at(MY_FIRST_INDEX + i);
          for(int nb : pdf_bin_schemes){
            std::stringstream ss;
            ss << "DeltaY_xi_reco_" << nb << "_PDF_" << (i + 1);
            std::string hname = ss.str();
            auto it = h_pdf_xi_reco_map.find(hname);
            if(it != h_pdf_xi_reco_map.end()){
              it->second->Fill(deltay_xi, weight * pdfw / orig_w * w_noac_nominal);
            }
          }
        }
        //template method end
      }
    }// loop ending for ttbar only
    //else{
    if (debug)cout << "checking for all MC" << endl;
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
    bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction();
    // h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("DeepAK8TopTags");
    // h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");
    if (debug)cout << "about to define vectors" << endl;  
    vector <TopJet> TopTaggedJets = event.get(h_AK8TopTags); 
    if (debug)cout << "define AK8 CHSJets" << endl;   
    vector <Jet> AK4CHSjets_matched = event.get(h_CHSjets_matched);   
    if (debug)cout << "define AK4 CHSJets" << endl;               // AK4Puppijets that have been matched to CHSjets
    vector <float> jets_hadronic_bscores; 
    if (debug)cout << "define float for bscore" << endl;
    // float Mreco = BestZprimeCandidate->Zprime_v4().M();
    
    // bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction();
    // LorentzVector toplep = BestZprimeCandidate->top_leptonic_v4();
    LorentzVector tophad = BestZprimeCandidate->top_hadronic_v4();
    float pt_hadTop_thresh = 150;   
    float pt_hadTop = tophad.pt();
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
    if(is_toptag_reconstruction){
        // Loop over hadronic top's subjets to extract highest bscore
      for(unsigned int i=0; i < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); i++){
        float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet(); // Using DeepJet btag score
        if(bscore > bscore_max) bscore_max = bscore;
      }
    
    }
    if (debug)cout << "about to define boost"<<endl;
    TLorentzVector had_top_b(0, 0, 0, 0);

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

      // Lepton 4-vector
        TLorentzVector lep_top_lep(0, 0, 0, 0);
        LorentzVector lep = BestZprimeCandidate->lepton().v4();
        lep_top_lep.SetPtEtaPhiE(lep.pt(), lep.eta(), lep.phi(), lep.E());
        //------------------------------------Define 4vectors of hadronic b-jet and lepton------------------------------------//


        //-------------------------------- Begin boosting top quarks and their decay products --------------------------------//
        // Define 4vectors of top quarks
        TLorentzVector PosTop(0, 0, 0, 0);
        TLorentzVector NegTop(0, 0, 0, 0);

        // POSITIVE LEPTON CONFIGURATION => Positive charged lepton has Positive Top mother
        if(BestZprimeCandidate->lepton().charge() > 0){
          // Define ttbar system
          PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                              BestZprimeCandidate->top_leptonic_v4().eta(), 
                              BestZprimeCandidate->top_leptonic_v4().phi(), 
                              BestZprimeCandidate->top_leptonic_v4().energy());
          NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                              BestZprimeCandidate->top_hadronic_v4().eta(), 
                              BestZprimeCandidate->top_hadronic_v4().phi(), 
                              BestZprimeCandidate->top_hadronic_v4().energy());
        }
        else if (BestZprimeCandidate->lepton().charge() < 0)
        {
        PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                              BestZprimeCandidate->top_hadronic_v4().eta(), 
                              BestZprimeCandidate->top_hadronic_v4().phi(), 
                              BestZprimeCandidate->top_hadronic_v4().energy());
        NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                              BestZprimeCandidate->top_leptonic_v4().eta(), 
                              BestZprimeCandidate->top_leptonic_v4().phi(), 
                              BestZprimeCandidate->top_leptonic_v4().energy());
        }
        
          TLorentzVector ttbar = PosTop + NegTop;
          // Boost into ttbar CoM-Frame <<<-------//
          lep_top_lep.Boost(-ttbar.BoostVector());
          had_top_b.Boost(-ttbar.BoostVector());
          PosTop.Boost(-ttbar.BoostVector());
          NegTop.Boost(-ttbar.BoostVector());

          // Rotate vectors into Helicity Frame <<<------//
          // Rotate about beamline
          lep_top_lep.RotateZ(-1.*PosTop.Phi());
          had_top_b.RotateZ(-1.*PosTop.Phi());
          PosTop.RotateZ(-1.*PosTop.Phi());
          NegTop.RotateZ(-1.*PosTop.Phi());
          // Rotate about y-axis
          lep_top_lep.RotateY(-1.*PosTop.Theta());
          had_top_b.RotateY(-1.*PosTop.Theta());
          PosTop.RotateY(-1.*PosTop.Theta());
          NegTop.RotateY(-1.*PosTop.Theta());

          // Boost into ttbar Rest-Frame <<<--------//
          lep_top_lep.Boost(-PosTop.BoostVector()); // Positive charged lepton has Positive Top mother
          had_top_b.Boost(-NegTop.BoostVector());   // Positive charged lepton means b-jet has Negative Top mother

        //-------------------------------- End boosting top quarks and their decay products --------------------------------//
    
        // Define angular variables as sum and difference of decay products' phi-coordinates
        // sphi and dphi = PosTopDecayProd_phi +- NegTopDecayProd_phi
        if (debug)cout << "done with boost"<<endl;
        float dphi=0.;
        float sphi = lep_top_lep.Phi() + had_top_b.Phi();
        if(BestZprimeCandidate->lepton().charge() > 0){ // lepton is Positive Top's Decay Product
          dphi = lep_top_lep.Phi() - had_top_b.Phi();
        }
        if(BestZprimeCandidate->lepton().charge() < 0)
        {
          dphi = had_top_b.Phi() - lep_top_lep.Phi();
        }
        if (debug)cout << "spin corr vars defined"<<endl;
          // Map back into original domain if necessary
        if(sphi > TMath::Pi()) sphi = sphi - 2*TMath::Pi();
        if(sphi < -TMath::Pi()) sphi = sphi + 2*TMath::Pi();
        if(dphi > TMath::Pi()) dphi = dphi - 2*TMath::Pi();
        if(dphi < -TMath::Pi()) dphi = dphi + 2*TMath::Pi();
        
    
      float deltay=(TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()));
      int MY_FIRST_INDEX = 9;
      if ( is_dy || is_wjets || is_qcd_HTbinned || is_alps || is_azh || is_htott_scalar || is_htott_pseudo || is_zprimetott ) MY_FIRST_INDEX = 47;
      if(event.genInfo->systweights().size() > (unsigned int) 100 + MY_FIRST_INDEX){
        float orig_weight = event.genInfo->originalXWGTUP();
        for(int i=0; i<100; i++){
          double pdf_weight = event.genInfo->systweights().at(i+MY_FIRST_INDEX);
          const char* name = hist_names[i].c_str();
          const char* name_dy_d1 = hist_names_dy_d1[i].c_str();
          const char* name_dy_d2 = hist_names_dy_d2[i].c_str();
          const char* name_sigma_1 = hist_names_sigma_1[i].c_str();
          const char* name_sigma_2 = hist_names_sigma_2[i].c_str();
          if (debug)cout <<" about to fill histos" <<endl;
          if (pt_hadTop > pt_hadTop_thresh && deltay >0){
            hist(name_sigma_1)->Fill(sphi,weight * pdf_weight / orig_weight);
            if (debug)cout <<" done with sigma 1" <<endl;
          }
          if (pt_hadTop > pt_hadTop_thresh && deltay <0){
            hist(name_sigma_2)->Fill(sphi,weight * pdf_weight / orig_weight);
            if (debug)cout <<" done with sigma 2" <<endl;
          }
          if(pt_hadTop < pt_hadTop_thresh && dphi >0){
            hist(name_dy_d1)->Fill(deltay,weight * pdf_weight / orig_weight);
            if (debug)cout <<" done with dy 1" <<endl;
          }
          if(pt_hadTop < pt_hadTop_thresh && dphi < 0){
            hist(name_dy_d2)->Fill(deltay,weight * pdf_weight / orig_weight);
            if (debug)cout <<" done with dy 2" <<endl;
          }
          hist(name)->Fill(deltay,weight * pdf_weight / orig_weight);
          if (debug)cout <<" done with dy" <<endl;
        }
      }
    //} get rid of else for ttbar 
  }

  

}

ZprimeSemiLeptonicPDFHists::~ZprimeSemiLeptonicPDFHists(){}