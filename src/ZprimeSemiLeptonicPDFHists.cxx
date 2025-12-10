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
    // Create unique histogram names for each PDF set
    std::stringstream ss_name;
    std::stringstream ss_name_tt;
    std::stringstream ss_name_dy_d1;
    std::stringstream ss_name_dy_d2;
    std::stringstream ss_name_sigma_1;
    std::stringstream ss_name_sigma_2;
    std::stringstream ss_name_xi;
    std::stringstream ss_name_cos_theta1k_antiLep;
    std::stringstream ss_name_cos_theta1r_antiLep;
    std::stringstream ss_name_cos_theta1n_antiLep;
    std::stringstream ss_name_cos_theta1kStar_antiLep;
    std::stringstream ss_name_cos_theta1rStar_antiLep;
    std::stringstream ss_name_cos_theta2k_Lep;
    std::stringstream ss_name_cos_theta2r_Lep;
    std::stringstream ss_name_cos_theta2n_Lep;
    std::stringstream ss_name_cos_theta2kStar_Lep;
    std::stringstream ss_name_cos_theta2rStar_Lep;
    std::stringstream ss_name_cos_theta1k;
    std::stringstream ss_name_cos_theta1r;
    std::stringstream ss_name_cos_theta1n;
    std::stringstream ss_name_cos_theta1kStar;
    std::stringstream ss_name_cos_theta1rStar;
    std::stringstream ss_name_cos_theta2k;
    std::stringstream ss_name_cos_theta2r;
    std::stringstream ss_name_cos_theta2n;
    std::stringstream ss_name_cos_theta2kStar;
    std::stringstream ss_name_cos_theta2rStar;
    std::stringstream ss_name_Ckk;
    std::stringstream ss_name_Crr;
    std::stringstream ss_name_Cnn;
    std::stringstream ss_name_Crk_plus;
    std::stringstream ss_name_Crk_minus;
    std::stringstream ss_name_Cnr_plus;
    std::stringstream ss_name_Cnr_minus;
    std::stringstream ss_name_Cnk_plus;
    std::stringstream ss_name_Cnk_minus;
    std::stringstream ss_name_cHel_Mtt300_400;
    std::stringstream ss_name_cHel_Mtt300_400_betaLT0p9;
    std::stringstream ss_name_cHel_P3n_Mtt800_Inf;
    std::stringstream ss_name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4;

    ss_name         << "DeltaY_PDF_" << i+1;
    ss_name_tt      << "DeltaY_PDF_RM_" << i+1;
    ss_name_dy_d1   << "DeltaY_reco_d1_PDF_" << i+1;
    ss_name_dy_d2   << "DeltaY_reco_d2_PDF_" << i+1;
    ss_name_sigma_1 << "Sigma_phi_1_PDF_" << i+1;
    ss_name_sigma_2 << "Sigma_phi_2_PDF_" << i+1;
    ss_name_xi      << "DeltaY_xi_reco_6_PDF_" << i+1;
    ss_name_cos_theta1k_antiLep << "cos_theta1k_antiLep_PDF_" << i+1;
    ss_name_cos_theta1r_antiLep << "cos_theta1r_antiLep_PDF_" << i+1;
    ss_name_cos_theta1n_antiLep << "cos_theta1n_antiLep_PDF_" << i+1;
    ss_name_cos_theta1kStar_antiLep << "cos_theta1kStar_antiLep_PDF_" << i+1;
    ss_name_cos_theta1rStar_antiLep << "cos_theta1rStar_antiLep_PDF_" << i+1;
    ss_name_cos_theta2k_Lep << "cos_theta2k_Lep_PDF_" << i+1;
    ss_name_cos_theta2r_Lep << "cos_theta2r_Lep_PDF_" << i+1;
    ss_name_cos_theta2n_Lep << "cos_theta2n_Lep_PDF_" << i+1;
    ss_name_cos_theta2kStar_Lep << "cos_theta2kStar_Lep_PDF_" << i+1;
    ss_name_cos_theta2rStar_Lep << "cos_theta2rStar_Lep_PDF_" << i+1;
    ss_name_cos_theta1k << "cos_theta1k_PDF_" << i+1;
    ss_name_cos_theta1r << "cos_theta1r_PDF_" << i+1;
    ss_name_cos_theta1n << "cos_theta1n_PDF_" << i+1;
    ss_name_cos_theta1kStar << "cos_theta1kStar_PDF_" << i+1;
    ss_name_cos_theta1rStar << "cos_theta1rStar_PDF_" << i+1;
    ss_name_cos_theta2k << "cos_theta2k_PDF_" << i+1;
    ss_name_cos_theta2r << "cos_theta2r_PDF_" << i+1;
    ss_name_cos_theta2n << "cos_theta2n_PDF_" << i+1;
    ss_name_cos_theta2kStar << "cos_theta2kStar_PDF_" << i+1;
    ss_name_cos_theta2rStar << "cos_theta2rStar_PDF_" << i+1;
    ss_name_Ckk << "Ckk_PDF_" << i+1;
    ss_name_Crr << "Crr _PDF_" << i+1;
    ss_name_Cnn << "Cnn_PDF_" << i+1;
    ss_name_Crk_plus << "Crk_plus_PDF_" << i+1;
    ss_name_Crk_minus << "Crk_minus_PDF_" << i+1;
    ss_name_Cnr_plus << "Cnr_plus_PDF_" << i+1;
    ss_name_Cnr_minus << "Cnr_minus_PDF_" << i+1;
    ss_name_Cnk_plus << "Cnk_plus_PDF_" << i+1;
    ss_name_Cnk_minus << "Cnk_minus_PDF_" << i+1;
    ss_name_cHel_Mtt300_400 << "cHel_Mtt300_400_PDF_" << i+1;
    ss_name_cHel_Mtt300_400_betaLT0p9 << "cHel_Mtt300_400_betaLT0p9_PDF_" << i+1;
    ss_name_cHel_P3n_Mtt800_Inf << "cHel_P3n_Mtt800_Inf_PDF_" << i+1;
    ss_name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 << "cHel_P3n_Mtt800_Inf_cosThetaLT0p4_PDF_" << i+1;


    // Create histogram titles
    stringstream ss_title;
    stringstream ss_title_tt;
    stringstream ss_title_dy_d1;
    stringstream ss_title_dy_d2;
    stringstream ss_title_sigma_1;
    stringstream ss_title_sigma_2;
    stringstream ss_title_delta_1;
    stringstream ss_title_delta_2;
    stringstream ss_title_xi;
    stringstream ss_title_cos_theta1k_antiLep;
    stringstream ss_title_cos_theta1r_antiLep;
    stringstream ss_title_cos_theta1n_antiLep;
    stringstream ss_title_cos_theta1kStar_antiLep;
    stringstream ss_title_cos_theta1rStar_antiLep;
    stringstream ss_title_cos_theta2k_Lep;
    stringstream ss_title_cos_theta2r_Lep;
    stringstream ss_title_cos_theta2n_Lep;
    stringstream ss_title_cos_theta2kStar_Lep;
    stringstream ss_title_cos_theta2rStar_Lep;
    stringstream ss_title_cos_theta1k;
    stringstream ss_title_cos_theta1r;
    stringstream ss_title_cos_theta1n;
    stringstream ss_title_cos_theta1kStar;
    stringstream ss_title_cos_theta1rStar;
    stringstream ss_title_cos_theta2k;
    stringstream ss_title_cos_theta2r;
    stringstream ss_title_cos_theta2n;
    stringstream ss_title_cos_theta2kStar;
    stringstream ss_title_cos_theta2rStar;
    stringstream ss_title_Ckk;
    stringstream ss_title_Crr;
    stringstream ss_title_Cnn;
    stringstream ss_title_Crk_plus;
    stringstream ss_title_Crk_minus;
    stringstream ss_title_Cnr_plus;
    stringstream ss_title_Cnr_minus;
    stringstream ss_title_Cnk_plus;
    stringstream ss_title_Cnk_minus;
    stringstream ss_title_cHel_Mtt300_400;
    stringstream ss_title_cHel_Mtt300_400_betaLT0p9;
    stringstream ss_title_cHel_P3n_Mtt800_Inf;
    stringstream ss_title_cHel_P3n_Mtt800_Inf_cosThetaLT0p4;

    ss_title         << "#DeltaY_{t#bar{t}} for PDF No. "  << i+1 << " out of 100" ;
    ss_title_tt      << "#DeltaY_{t#bar{t}} RM for PDF No. "  << i+1 << " out of 100" ;
    ss_title_dy_d1   << "#DeltaY_{t#bar{t}} for #Delta #phi >0 for PDF No. "<< i+1 << " out of 100" ;
    ss_title_dy_d2   << "#DeltaY_{t#bar{t}} for #Delta #phi <0 for PDF No. "<< i+1 << " out of 100" ;
    ss_title_sigma_1 << "#Sigma #phi for #DeltaY >0 for PDF No. "<< i+1 << " out of 100" ;
    ss_title_sigma_2 << "#Sigma #phi for #DeltaY <0 for PDF No. "<< i+1 << " out of 100" ;
    ss_title_delta_1 << "#Delta #phi for #DeltaY >0 for PDF No. "<< i+1 << " out of 100" ;
    ss_title_delta_2 << "#Delta #phi for #DeltaY <0 for PDF No. "<< i+1 << " out of 100" ;
    ss_title_xi      << "tanh(#Delta y)_{reco} for PDF No. " << i+1 << " out of 100";
    ss_title_cos_theta1k_antiLep << "cos(#theta_{antilep}^{k}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta1r_antiLep << "cos(#theta_{antilep}^{r}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta1n_antiLep << "cos(#theta_{antilep}^{n}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta1kStar_antiLep << "cos(#theta_{antilep}^{k*}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta1rStar_antiLep << "cos(#theta_{antilep}^{r*}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2k_Lep << "cos(#theta_{lep}^{k}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2r_Lep << "cos(#theta_{lep}^{r}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2n_Lep << "cos(#theta_{lep}^{n}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2kStar_Lep << "cos(#theta_{lep}^{k*}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2rStar_Lep << "cos(#theta_{lep}^{r*}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta1k << "cos(#theta_{1}^{k}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta1r << "cos(#theta_{1}^{r}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta1n << "cos(#theta_{1}^{n}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta1kStar << "cos(#theta_{1}^{k*}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta1rStar << "cos(#theta_{1}^{r*}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2k << "cos(#theta_{2}^{k}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2r << "cos(#theta_{2}^{r}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2n << "cos(#theta_{2}^{n}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2kStar << "cos(#theta_{2}^{k*}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cos_theta2rStar << "cos(#theta_{2}^{r*}) for PDF No. " <<  i+1 << " out of 100";
    ss_title_Ckk << "C_{kk} for PDF No. " <<  i+1 << " out of 100";
    ss_title_Crr << "C_{rr} for PDF No. " <<  i+1 << " out of 100";
    ss_title_Cnn << "C_{nn} for PDF No. " <<  i+1 << " out of 100";
    ss_title_Crk_plus << "C_{rk} + C_{kr} for PDF No. " <<  i+1 << " out of 100";
    ss_title_Crk_minus << "C_{rk} - C_{kr} for PDF No. " <<  i+1 << " out of 100";
    ss_title_Cnr_plus << "C_{nr} + C_{rn} for PDF No. " <<  i+1 << " out of 100";
    ss_title_Cnr_minus << "C_{nr} - C_{rn} for PDF No. " <<  i+1 << " out of 100";
    ss_title_Cnk_plus << "C_{nk} + C_{kn} for PDF No. " <<  i+1 << " out of 100";
    ss_title_Cnk_minus << "C_{nk} - C_{kn} for PDF No. " <<  i+1 << " out of 100";
    ss_title_cHel_Mtt300_400 << "cos(#phi_{lb}) (M_{tt} [300,400] GeV) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cHel_Mtt300_400_betaLT0p9 << "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cHel_P3n_Mtt800_Inf << "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV) for PDF No. " <<  i+1 << " out of 100";
    ss_title_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 << "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4) for PDF No. " <<  i+1 << " out of 100";


    // Convert names and titles to std::string
    std::string s_name = ss_name.str();
    std::string s_name_tt = ss_name_tt.str();
    std::string s_name_dy_d1 = ss_name_dy_d1.str();
    std::string s_name_dy_d2 = ss_name_dy_d2.str();
    std::string s_name_sigma_1 = ss_name_sigma_1.str();
    std::string s_name_sigma_2 = ss_name_sigma_2.str();
    std::string s_name_xi    = ss_name_xi.str();
    std::string s_name_cos_theta1k_antiLep = ss_name_cos_theta1k_antiLep.str();
    std::string s_name_cos_theta1r_antiLep = ss_name_cos_theta1r_antiLep.str();
    std::string s_name_cos_theta1n_antiLep = ss_name_cos_theta1n_antiLep.str();
    std::string s_name_cos_theta1kStar_antiLep = ss_name_cos_theta1kStar_antiLep.str();
    std::string s_name_cos_theta1rStar_antiLep = ss_name_cos_theta1rStar_antiLep.str();
    std::string s_name_cos_theta2k_Lep = ss_name_cos_theta2k_Lep.str();
    std::string s_name_cos_theta2r_Lep = ss_name_cos_theta2r_Lep.str();
    std::string s_name_cos_theta2n_Lep = ss_name_cos_theta2n_Lep.str();
    std::string s_name_cos_theta2kStar_Lep = ss_name_cos_theta2kStar_Lep.str();
    std::string s_name_cos_theta2rStar_Lep = ss_name_cos_theta2rStar_Lep.str();
    std::string s_name_cos_theta1k = ss_name_cos_theta1k.str();
    std::string s_name_cos_theta1r = ss_name_cos_theta1r.str();
    std::string s_name_cos_theta1n = ss_name_cos_theta1n.str();
    std::string s_name_cos_theta1kStar = ss_name_cos_theta1kStar.str();
    std::string s_name_cos_theta1rStar = ss_name_cos_theta1rStar.str();
    std::string s_name_cos_theta2k = ss_name_cos_theta2k.str();
    std::string s_name_cos_theta2r = ss_name_cos_theta2r.str();
    std::string s_name_cos_theta2n = ss_name_cos_theta2n.str();
    std::string s_name_cos_theta2kStar = ss_name_cos_theta2kStar.str();
    std::string s_name_cos_theta2rStar = ss_name_cos_theta2rStar.str();
    std::string s_name_Ckk = ss_name_Ckk.str();
    std::string s_name_Crr = ss_name_Crr.str();
    std::string s_name_Cnn = ss_name_Cnn.str();
    std::string s_name_Crk_plus = ss_name_Crk_plus.str();
    std::string s_name_Crk_minus = ss_name_Crk_minus.str();
    std::string s_name_Cnr_plus = ss_name_Cnr_plus.str();
    std::string s_name_Cnr_minus = ss_name_Cnr_minus.str();
    std::string s_name_Cnk_plus = ss_name_Cnk_plus.str();
    std::string s_name_Cnk_minus = ss_name_Cnk_minus.str();
    std::string s_name_cHel_Mtt300_400 = ss_name_cHel_Mtt300_400.str();
    std::string s_name_cHel_Mtt300_400_betaLT0p9 = ss_name_cHel_Mtt300_400_betaLT0p9.str();
    std::string s_name_cHel_P3n_Mtt800_Inf = ss_name_cHel_P3n_Mtt800_Inf.str();
    std::string s_name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = ss_name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.str();

    std::string s_title = ss_title.str();
    std::string s_title_tt = ss_title_tt.str();
    std::string s_title_dy_d1 = ss_title_dy_d1.str();
    std::string s_title_dy_d2 = ss_title_dy_d2.str();
    std::string s_title_sigma_1 = ss_title_sigma_1.str();
    std::string s_title_sigma_2 = ss_title_sigma_2.str();
    std::string s_title_delta_1 = ss_title_delta_1.str();
    std::string s_title_delta_2 = ss_title_delta_2.str();
    std::string s_title_xi    = ss_title_xi.str();
    std::string s_title_cos_theta1k_antiLep = ss_title_cos_theta1k_antiLep.str();
    std::string s_title_cos_theta1r_antiLep = ss_title_cos_theta1r_antiLep.str();
    std::string s_title_cos_theta1n_antiLep = ss_title_cos_theta1n_antiLep.str();
    std::string s_title_cos_theta1kStar_antiLep = ss_title_cos_theta1kStar_antiLep.str();
    std::string s_title_cos_theta1rStar_antiLep = ss_title_cos_theta1rStar_antiLep.str();
    std::string s_title_cos_theta2k_Lep = ss_title_cos_theta2k_Lep.str();
    std::string s_title_cos_theta2r_Lep = ss_title_cos_theta2r_Lep.str();
    std::string s_title_cos_theta2n_Lep = ss_title_cos_theta2n_Lep.str();
    std::string s_title_cos_theta2kStar_Lep = ss_title_cos_theta2kStar_Lep.str();
    std::string s_title_cos_theta2rStar_Lep = ss_title_cos_theta2rStar_Lep.str();
    std::string s_title_cos_theta1k = ss_title_cos_theta1k.str();
    std::string s_title_cos_theta1r = ss_title_cos_theta1r.str();
    std::string s_title_cos_theta1n = ss_title_cos_theta1n.str();
    std::string s_title_cos_theta1kStar = ss_title_cos_theta1kStar.str();
    std::string s_title_cos_theta1rStar = ss_title_cos_theta1rStar.str();
    std::string s_title_cos_theta2k = ss_title_cos_theta2k.str();
    std::string s_title_cos_theta2r = ss_title_cos_theta2r.str();
    std::string s_title_cos_theta2n = ss_title_cos_theta2n.str();
    std::string s_title_cos_theta2kStar = ss_title_cos_theta2kStar.str();
    std::string s_title_cos_theta2rStar = ss_title_cos_theta2rStar.str();
    std::string s_title_Ckk = ss_title_Ckk.str();
    std::string s_title_Crr = ss_title_Crr.str();
    std::string s_title_Cnn = ss_title_Cnn.str();
    std::string s_title_Crk_plus = ss_title_Crk_plus.str();
    std::string s_title_Crk_minus = ss_title_Crk_minus.str();
    std::string s_title_Cnr_plus = ss_title_Cnr_plus.str();
    std::string s_title_Cnr_minus = ss_title_Cnr_minus.str();
    std::string s_title_Cnk_plus = ss_title_Cnk_plus.str();
    std::string s_title_Cnk_minus = ss_title_Cnk_minus.str();
    std::string s_title_cHel_Mtt300_400 = ss_title_cHel_Mtt300_400.str();
    std::string s_title_cHel_Mtt300_400_betaLT0p9 = ss_title_cHel_Mtt300_400_betaLT0p9.str();
    std::string s_title_cHel_P3n_Mtt800_Inf = ss_title_cHel_P3n_Mtt800_Inf.str();
    std::string s_title_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = ss_title_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.str();


    // Convert name and title strings to const char*
    const char* char_name = s_name.c_str();
    const char* char_name_tt = s_name_tt.c_str();
    const char* char_name_dy_d1 = s_name_dy_d1.c_str();
    const char* char_name_dy_d2 = s_name_dy_d2.c_str();
    const char* char_name_sigma_1 = s_name_sigma_1.c_str();
    const char* char_name_sigma_2 = s_name_sigma_2.c_str();
    const char* char_name_cos_theta1k_antiLep = s_name_cos_theta1k_antiLep.c_str();
    const char* char_name_cos_theta1r_antiLep = s_name_cos_theta1r_antiLep.c_str();
    const char* char_name_cos_theta1n_antiLep = s_name_cos_theta1n_antiLep.c_str();
    const char* char_name_cos_theta1kStar_antiLep = s_name_cos_theta1kStar_antiLep.c_str();
    const char* char_name_cos_theta1rStar_antiLep = s_name_cos_theta1rStar_antiLep.c_str();
    const char* char_name_cos_theta2k_Lep = s_name_cos_theta2k_Lep.c_str();
    const char* char_name_cos_theta2r_Lep = s_name_cos_theta2r_Lep.c_str();
    const char* char_name_cos_theta2n_Lep = s_name_cos_theta2n_Lep.c_str();
    const char* char_name_cos_theta2kStar_Lep = s_name_cos_theta2kStar_Lep.c_str();
    const char* char_name_cos_theta2rStar_Lep = s_name_cos_theta2rStar_Lep.c_str();
    const char* char_name_cos_theta1k = s_name_cos_theta1k.c_str();
    const char* char_name_cos_theta1r = s_name_cos_theta1r.c_str();
    const char* char_name_cos_theta1n = s_name_cos_theta1n.c_str();
    const char* char_name_cos_theta1kStar = s_name_cos_theta1kStar.c_str();
    const char* char_name_cos_theta1rStar = s_name_cos_theta1rStar.c_str();
    const char* char_name_cos_theta2k = s_name_cos_theta2k.c_str();
    const char* char_name_cos_theta2r = s_name_cos_theta2r.c_str();
    const char* char_name_cos_theta2n = s_name_cos_theta2n.c_str();
    const char* char_name_cos_theta2kStar = s_name_cos_theta2kStar.c_str();
    const char* char_name_cos_theta2rStar = s_name_cos_theta2rStar.c_str();
    const char* char_name_Ckk = s_name_Ckk.c_str();
    const char* char_name_Crr = s_name_Crr.c_str();
    const char* char_name_Cnn = s_name_Cnn.c_str();
    const char* char_name_Crk_plus = s_name_Crk_plus.c_str();
    const char* char_name_Crk_minus = s_name_Crk_minus.c_str();
    const char* char_name_Cnr_plus = s_name_Cnr_plus.c_str();
    const char* char_name_Cnr_minus = s_name_Cnr_minus.c_str();
    const char* char_name_Cnk_plus = s_name_Cnk_plus.c_str();
    const char* char_name_Cnk_minus = s_name_Cnk_minus.c_str();
    const char* char_name_cHel_Mtt300_400 = s_name_cHel_Mtt300_400.c_str();
    const char* char_name_cHel_Mtt300_400_betaLT0p9 = s_name_cHel_Mtt300_400_betaLT0p9.c_str();
    const char* char_name_cHel_P3n_Mtt800_Inf = s_name_cHel_P3n_Mtt800_Inf.c_str();
    const char* char_name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = s_name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.c_str();
    
    const char* char_title = s_title.c_str();
    const char* char_title_tt = s_title_tt.c_str();
    const char* char_title_dy_d1 = s_title_dy_d1.c_str();
    const char* char_title_dy_d2 = s_title_dy_d2.c_str();
    const char* char_title_sigma_1 = s_title_sigma_1.c_str();
    const char* char_title_sigma_2 = s_title_sigma_2.c_str();
    const char* char_title_delta_1 = s_title_delta_1.c_str();
    const char* char_title_delta_2 = s_title_delta_2.c_str();
    const char* char_title_cos_theta1k_antiLep = s_title_cos_theta1k_antiLep.c_str();
    const char* char_title_cos_theta1r_antiLep = s_title_cos_theta1r_antiLep.c_str();
    const char* char_title_cos_theta1n_antiLep = s_title_cos_theta1n_antiLep.c_str();
    const char* char_title_cos_theta1kStar_antiLep = s_title_cos_theta1kStar_antiLep.c_str();
    const char* char_title_cos_theta1rStar_antiLep = s_title_cos_theta1rStar_antiLep.c_str();
    const char* char_title_cos_theta2k_Lep = s_title_cos_theta2k_Lep.c_str();
    const char* char_title_cos_theta2r_Lep = s_title_cos_theta2r_Lep.c_str();
    const char* char_title_cos_theta2n_Lep = s_title_cos_theta2n_Lep.c_str();
    const char* char_title_cos_theta2kStar_Lep = s_title_cos_theta2kStar_Lep.c_str();
    const char* char_title_cos_theta2rStar_Lep = s_title_cos_theta2rStar_Lep.c_str();
    const char* char_title_cos_theta1k = s_title_cos_theta1k.c_str();
    const char* char_title_cos_theta1r = s_title_cos_theta1r.c_str();
    const char* char_title_cos_theta1n = s_title_cos_theta1n.c_str();
    const char* char_title_cos_theta1kStar = s_title_cos_theta1kStar.c_str();
    const char* char_title_cos_theta1rStar = s_title_cos_theta1rStar.c_str();
    const char* char_title_cos_theta2k = s_title_cos_theta2k.c_str();
    const char* char_title_cos_theta2r = s_title_cos_theta2r.c_str();
    const char* char_title_cos_theta2n = s_title_cos_theta2n.c_str();
    const char* char_title_cos_theta2kStar = s_title_cos_theta2kStar.c_str();
    const char* char_title_cos_theta2rStar = s_title_cos_theta2rStar.c_str();
    const char* char_title_Ckk = s_title_Ckk.c_str();
    const char* char_title_Crr = s_title_Crr.c_str();
    const char* char_title_Cnn = s_title_Cnn.c_str();
    const char* char_title_Crk_plus = s_title_Crk_plus.c_str();
    const char* char_title_Crk_minus = s_title_Crk_minus.c_str();
    const char* char_title_Cnr_plus = s_title_Cnr_plus.c_str();
    const char* char_title_Cnr_minus = s_title_Cnr_minus.c_str();
    const char* char_title_Cnk_plus = s_title_Cnk_plus.c_str();
    const char* char_title_Cnk_minus = s_title_Cnk_minus.c_str();
    const char* char_title_cHel_Mtt300_400 = s_title_cHel_Mtt300_400.c_str();
    const char* char_title_cHel_Mtt300_400_betaLT0p9 = s_title_cHel_Mtt300_400_betaLT0p9.c_str();
    const char* char_title_cHel_P3n_Mtt800_Inf = s_title_cHel_P3n_Mtt800_Inf.c_str();
    const char* char_title_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = s_title_cHel_P3n_Mtt800_Inf_cosThetaLT0p4.c_str();


    // Store histogram names for later use
    hist_names[i] = s_name;
    hist_names_tt[i] = s_name_tt;
    hist_names_dy_d1[i] = s_name_dy_d1;
    hist_names_dy_d2[i] = s_name_dy_d2;
    hist_names_sigma_1[i] = s_name_sigma_1;
    hist_names_sigma_2[i] = s_name_sigma_2;
    hist_names_xi[i]    = s_name_xi;
    hist_names_cos_theta1k_antiLep[i] = s_name_cos_theta1k_antiLep;
    hist_names_cos_theta1r_antiLep[i] = s_name_cos_theta1r_antiLep;
    hist_names_cos_theta1n_antiLep[i] = s_name_cos_theta1n_antiLep;
    hist_names_cos_theta1kStar_antiLep[i] = s_name_cos_theta1kStar_antiLep;
    hist_names_cos_theta1rStar_antiLep[i] = s_name_cos_theta1rStar_antiLep;
    hist_names_cos_theta2k_Lep[i] = s_name_cos_theta2k_Lep;
    hist_names_cos_theta2r_Lep[i] = s_name_cos_theta2r_Lep;
    hist_names_cos_theta2n_Lep[i] = s_name_cos_theta2n_Lep;
    hist_names_cos_theta2kStar_Lep[i] = s_name_cos_theta2kStar_Lep;
    hist_names_cos_theta2rStar_Lep[i] = s_name_cos_theta2rStar_Lep;
    hist_names_cos_theta1k[i] = s_name_cos_theta1k;
    hist_names_cos_theta1r[i] = s_name_cos_theta1r;
    hist_names_cos_theta1n[i] = s_name_cos_theta1n;
    hist_names_cos_theta1kStar[i] = s_name_cos_theta1kStar;
    hist_names_cos_theta1rStar[i] = s_name_cos_theta1rStar;
    hist_names_cos_theta2k[i] = s_name_cos_theta2k;
    hist_names_cos_theta2r[i] = s_name_cos_theta2r;
    hist_names_cos_theta2n[i] = s_name_cos_theta2n;
    hist_names_cos_theta2kStar[i] = s_name_cos_theta2kStar;
    hist_names_cos_theta2rStar[i] = s_name_cos_theta2rStar;
    hist_names_Ckk[i] = s_name_Ckk;
    hist_names_Crr[i] = s_name_Crr;
    hist_names_Cnn[i] = s_name_Cnn;
    hist_names_Crk_plus[i]  = s_name_Crk_plus;
    hist_names_Crk_minus[i] = s_name_Crk_minus;
    hist_names_Cnr_plus[i]  = s_name_Cnr_plus;
    hist_names_Cnr_minus[i] = s_name_Cnr_minus;
    hist_names_Cnk_plus[i]  = s_name_Cnk_plus;
    hist_names_Cnk_minus[i] = s_name_Cnk_minus;
    hist_names_cHel_Mtt300_400[i]           = s_name_cHel_Mtt300_400;
    hist_names_cHel_Mtt300_400_betaLT0p9[i] = s_name_cHel_Mtt300_400_betaLT0p9;
    hist_names_cHel_P3n_Mtt800_Inf[i]               = s_name_cHel_P3n_Mtt800_Inf;
    hist_names_cHel_P3n_Mtt800_Inf_cosThetaLT0p4[i] = s_name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4;


    // Book histograms
    book<TH1F>(char_name, char_title,  2, -2.5, 2.5);
    book<TH2F>(char_name_tt, char_title_tt,  2, -2.5, 2.5, 2, -2.5, 2.5);
    book<TH1F>(char_name_dy_d1, char_title_dy_d1,  2, -2.5, 2.5);
    book<TH1F>(char_name_dy_d2, char_title_dy_d2,  2, -2.5, 2.5);
    book<TH1F>(char_name_sigma_1, char_title_sigma_1,  16, -3.2, 3.2);
    book<TH1F>(char_name_sigma_2, char_title_sigma_2,  16, -3.2, 3.2);
    book<TH1F>(s_name_xi.c_str(),    s_title_xi.c_str(),    /*nbins*/ 6,  -1.0,  1.0);
    book<TH1F>(char_name_cos_theta1k_antiLep, char_title_cos_theta1k_antiLep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta1r_antiLep, char_title_cos_theta1r_antiLep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta1n_antiLep, char_title_cos_theta1n_antiLep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta1kStar_antiLep, char_title_cos_theta1kStar_antiLep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta1rStar_antiLep, char_title_cos_theta1rStar_antiLep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2k_Lep, char_title_cos_theta2k_Lep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2r_Lep, char_title_cos_theta2r_Lep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2n_Lep, char_title_cos_theta2n_Lep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2kStar_Lep, char_title_cos_theta2kStar_Lep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2rStar_Lep, char_title_cos_theta2rStar_Lep, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta1k, char_title_cos_theta1k, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta1r, char_title_cos_theta1r, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta1n, char_title_cos_theta1n, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta1kStar, char_title_cos_theta1kStar,  24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta1rStar, char_title_cos_theta1rStar,  24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2k, char_title_cos_theta2k, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2r, char_title_cos_theta2r, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2n, char_title_cos_theta2n, 24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2kStar, char_title_cos_theta2kStar,  24, -1.0, 1.0);
    book<TH1F>(char_name_cos_theta2rStar, char_title_cos_theta2rStar,  24, -1.0, 1.0);
    book<TH1F>(char_name_Ckk, char_title_Ckk,  24, -1.0, 1.0);
    book<TH1F>(char_name_Crr, char_title_Crr,  24, -1.0, 1.0);
    book<TH1F>(char_name_Cnn, char_title_Cnn,  24, -1.0, 1.0);
    book<TH1F>(char_name_Crk_plus, char_title_Crk_plus,  24, -1.0, 1.0);
    book<TH1F>(char_name_Crk_minus, char_title_Crk_minus,  24, -1.0, 1.0);
    book<TH1F>(char_name_Cnr_plus, char_title_Cnr_plus,  24, -1.0, 1.0);
    book<TH1F>(char_name_Cnr_minus, char_title_Cnr_minus,  24, -1.0, 1.0);
    book<TH1F>(char_name_Cnk_plus, char_title_Cnk_plus,  24, -1.0, 1.0);
    book<TH1F>(char_name_Cnk_minus, char_title_Cnk_minus,  24, -1.0, 1.0);
    book<TH1F>(char_name_cHel_Mtt300_400, char_title_cHel_Mtt300_400,  24, -1.0, 1.0);
    book<TH1F>(char_name_cHel_Mtt300_400_betaLT0p9, char_title_cHel_Mtt300_400_betaLT0p9,  24, -1.0, 1.0);
    book<TH1F>(char_name_cHel_P3n_Mtt800_Inf, char_title_cHel_P3n_Mtt800_Inf,  24, -1.0, 1.0);
    book<TH1F>(char_name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4, char_title_cHel_P3n_Mtt800_Inf_cosThetaLT0p4,  24, -1.0, 1.0);
    
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

    if (debug)cout << "checking for all MC" << endl;

    // Start angular variable definitions -------------------------------------------------------------------------------------------------------//
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);     // Best Z' candidate from chi2 reconstruction
    bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction(); // Reconstruction process id
    vector <Jet> AK4CHSjets_matched = event.get(h_CHSjets_matched);                  // AK4Puppijets that have been matched to CHSjets
    vector <TopJet> TopTaggedJets = event.get(h_AK8TopTags);                         // AK8Puppi jets TopTagged by DeepAK8TopTagger
    vector <float> jets_hadronic_bscores;                                            // bScores vector for resolved hadronic jets
    float pt_hadTop_thresh = 150;                                                    // Define cut-variable as pt of hadTop for low/high regions
    float pt_hadTop = BestZprimeCandidate->top_hadronic_v4().pt();                   // pT of hadronic-top jet

    
    //-------------- Extracting highest b-tag score in Resolved topology, i.e. no top-tagged jet in event --------------//
    float bscore_max = -2;
    if(!is_toptag_reconstruction){
        // Loop over resolved hadronic jets to find their bscore via CHS jets
      for(unsigned int i=0; i<BestZprimeCandidate->jets_hadronic().size(); i++){
        double deltaR_min = 99;
        // Match resolved hadronic jets to CHS jets (which have bscores)
        for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
          double deltaR_CHS = deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(j));
          if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;
          }
        // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
        for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
          if(deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
          jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet()); // Using DeepJet btag score
          } 
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

    TLorentzVector lep_top_lep(0, 0, 0, 0); // Lepton 4-vector
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
    TLorentzVector PosTop_CoM = PosTop;
    TLorentzVector NegTop_CoM = NegTop;
    TLorentzVector lep_top_lep_CoM = lep_top_lep;
    TLorentzVector had_top_b_CoM = had_top_b;
    // Boost with negative of ttbar boost vector
    PosTop_CoM.Boost(-1*ttbar.BoostVector());
    NegTop_CoM.Boost(-1*ttbar.BoostVector());
    lep_top_lep_CoM.Boost(-1*ttbar.BoostVector());
    had_top_b_CoM.Boost(-1*ttbar.BoostVector());


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
    TLorentzVector had_top_b_Rest   = had_top_b_CoM;
    
    // Mother depends on lepton charge
    if(BestZprimeCandidate->lepton().charge() > 0){
      lep_top_lep_Rest.Boost(-1.*PosTop_CoM.BoostVector()); // lepton has Positive Top mother
      had_top_b_Rest.Boost(-1.*NegTop_CoM.BoostVector());   // b-jet has Negative Top mother
    }
    else if (BestZprimeCandidate->lepton().charge() < 0){
      lep_top_lep_Rest.Boost(-1.*NegTop_CoM.BoostVector()); // lepton has Negative Top mother
      had_top_b_Rest.Boost(-1.*PosTop_CoM.BoostVector());   // b-jet has Positive Top mother
    }

  
    //------------------------------------------- Spin Correlation variables -------------------------------------------//
    float cosTheta1k_antiLep = 99.;
    float cosTheta1r_antiLep = 99.;
    float cosTheta1n_antiLep = 99.;
    float cosTheta1kStar_antiLep = 99.;
    float cosTheta1rStar_antiLep = 99.;
    float cosTheta1k = 99.;
    float cosTheta1r = 99.;
    float cosTheta1n = 99.;
    float cosTheta1kStar = 99.;
    float cosTheta1rStar = 99.;
    
    float cosTheta2k_Lep = 99.;
    float cosTheta2r_Lep = 99.;
    float cosTheta2n_Lep = 99.;
    float cosTheta2kStar_Lep = 99.;
    float cosTheta2rStar_Lep = 99.;
    float cosTheta2k = 99.;
    float cosTheta2r = 99.;
    float cosTheta2n = 99.;
    float cosTheta2kStar = 99.;
    float cosTheta2rStar = 99.;

    float C_kk = 99.;
    float C_rr = 99.;
    float C_nn = 99.;
    
    float C_rk_plus = 99.;
    float C_rk_minus = 99.;
    float C_nr_plus = 99.;
    float C_nr_minus = 99.;
    float C_nk_plus = 99.;
    float C_nk_minus = 99.;

    float CHel = 99.;
    float CHel_Mtt300_400 = 99.;
    float CHel_Mtt300_400_betaLT0p9 = 99.;

    float CHel_P3n = 99.;
    float CHel_P3n_Mtt800_Inf = 99.;
    float CHel_P3n_Mtt800_Inf_cosThetaLT0p4 = 99.;

    // Use only leptons as spin-analyzers
    if(BestZprimeCandidate->lepton().charge() > 0){
      // anti-lepton is spin-analyzer fo top quark
      cosTheta1k_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(kbase);
      cosTheta1r_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(rbase);
      cosTheta1n_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(nbase);
      cosTheta1kStar_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(kStar);
      cosTheta1rStar_antiLep = lep_top_lep_Rest.Vect().Unit().Dot(rStar);
    }
    else if (BestZprimeCandidate->lepton().charge() < 0){
      // lepton is spin-analyzer for antitop quark
      cosTheta2k_Lep = lep_top_lep_Rest.Vect().Unit().Dot(kbase);
      cosTheta2r_Lep = lep_top_lep_Rest.Vect().Unit().Dot(rbase);
      cosTheta2n_Lep = lep_top_lep_Rest.Vect().Unit().Dot(nbase);
      cosTheta2kStar_Lep = lep_top_lep_Rest.Vect().Unit().Dot(kStar);
      cosTheta2rStar_Lep = lep_top_lep_Rest.Vect().Unit().Dot(rStar);
    }

    // Assign spin-analyzers depending on lepton charge
    if(BestZprimeCandidate->lepton().charge() > 0){
      // top quark spin-analyzer is lepton
      cosTheta1k = lep_top_lep_Rest.Vect().Unit().Dot(kbase);
      cosTheta1r = lep_top_lep_Rest.Vect().Unit().Dot(rbase);
      cosTheta1n = lep_top_lep_Rest.Vect().Unit().Dot(nbase);
      cosTheta1kStar = lep_top_lep_Rest.Vect().Unit().Dot(kStar);
      cosTheta1rStar = lep_top_lep_Rest.Vect().Unit().Dot(rStar);
      // antitop spin-analyzer is b-jet
      cosTheta2k = had_top_b_Rest.Vect().Unit().Dot(kbase);
      cosTheta2r = had_top_b_Rest.Vect().Unit().Dot(rbase);
      cosTheta2n = had_top_b_Rest.Vect().Unit().Dot(nbase);
      cosTheta2kStar = had_top_b_Rest.Vect().Unit().Dot(kStar);
      cosTheta2rStar = had_top_b_Rest.Vect().Unit().Dot(rStar);
    }
    else if (BestZprimeCandidate->lepton().charge() < 0){
      // top quark spin-analyzer is b-jet
      cosTheta1k = had_top_b_Rest.Vect().Unit().Dot(kbase);
      cosTheta1r = had_top_b_Rest.Vect().Unit().Dot(rbase);
      cosTheta1n = had_top_b_Rest.Vect().Unit().Dot(nbase);
      cosTheta1kStar = had_top_b_Rest.Vect().Unit().Dot(kStar);
      cosTheta1rStar = had_top_b_Rest.Vect().Unit().Dot(rStar);
      // antitop spin-analyzer is lepton
      cosTheta2k = lep_top_lep_Rest.Vect().Unit().Dot(kbase);
      cosTheta2r = lep_top_lep_Rest.Vect().Unit().Dot(rbase);
      cosTheta2n = lep_top_lep_Rest.Vect().Unit().Dot(nbase);
      cosTheta2kStar = lep_top_lep_Rest.Vect().Unit().Dot(kStar);
      cosTheta2rStar = lep_top_lep_Rest.Vect().Unit().Dot(rStar);
    }

    // correlation matrix elements
    C_nn = cosTheta1n * cosTheta2n;
    C_rr = cosTheta1r * cosTheta2r;
    C_kk = cosTheta1k * cosTheta2k;
    // sum and differences of cross correlations
    C_rk_plus = cosTheta1r * cosTheta2k + cosTheta1k * cosTheta2r;
    C_rk_minus = cosTheta1r * cosTheta2k - cosTheta1k * cosTheta2r;
    C_nr_plus = cosTheta1n * cosTheta2r + cosTheta1r * cosTheta2n;
    C_nr_minus = cosTheta1n * cosTheta2r - cosTheta1r * cosTheta2n;
    C_nk_plus = cosTheta1n * cosTheta2k + cosTheta1k * cosTheta2n;
    C_nk_minus = cosTheta1n * cosTheta2k - cosTheta1k * cosTheta2n;

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

  
    
    // Fill PDF systematic histograms for all MC samples ----------------------------------------------------------------------------------------//
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
        const char* name_cos_theta1k_antiLep = hist_names_cos_theta1k_antiLep[i].c_str();
        const char* name_cos_theta1r_antiLep = hist_names_cos_theta1r_antiLep[i].c_str();
        const char* name_cos_theta1n_antiLep = hist_names_cos_theta1n_antiLep[i].c_str();
        const char* name_cos_theta1kStar_antiLep = hist_names_cos_theta1kStar_antiLep[i].c_str();
        const char* name_cos_theta1rStar_antiLep = hist_names_cos_theta1rStar_antiLep[i].c_str();
        const char* name_cos_theta2k_Lep = hist_names_cos_theta2k_Lep[i].c_str();
        const char* name_cos_theta2r_Lep = hist_names_cos_theta2r_Lep[i].c_str();
        const char* name_cos_theta2n_Lep = hist_names_cos_theta2n_Lep[i].c_str();
        const char* name_cos_theta2kStar_Lep = hist_names_cos_theta2kStar_Lep[i].c_str();
        const char* name_cos_theta2rStar_Lep = hist_names_cos_theta2rStar_Lep[i].c_str();
        const char* name_cos_theta1k = hist_names_cos_theta1k[i].c_str();
        const char* name_cos_theta1r = hist_names_cos_theta1r[i].c_str();
        const char* name_cos_theta1n = hist_names_cos_theta1n[i].c_str();
        const char* name_cos_theta1kStar = hist_names_cos_theta1kStar[i].c_str();
        const char* name_cos_theta1rStar = hist_names_cos_theta1rStar[i].c_str();
        const char* name_cos_theta2k = hist_names_cos_theta2k[i].c_str();
        const char* name_cos_theta2r = hist_names_cos_theta2r[i].c_str();
        const char* name_cos_theta2n = hist_names_cos_theta2n[i].c_str();
        const char* name_cos_theta2kStar = hist_names_cos_theta2kStar[i].c_str();
        const char* name_cos_theta2rStar = hist_names_cos_theta2rStar[i].c_str();
        const char* name_Ckk = hist_names_Ckk[i].c_str();
        const char* name_Crr = hist_names_Crr[i].c_str();
        const char* name_Cnn = hist_names_Cnn[i].c_str();
        const char* name_Crk_plus = hist_names_Crk_plus[i].c_str();
        const char* name_Crk_minus = hist_names_Crk_minus[i].c_str();
        const char* name_Cnr_plus = hist_names_Cnr_plus[i].c_str();
        const char* name_Cnr_minus = hist_names_Cnr_minus[i].c_str();
        const char* name_Cnk_plus = hist_names_Cnk_plus[i].c_str();
        const char* name_Cnk_minus = hist_names_Cnk_minus[i].c_str();
        const char* name_cHel_Mtt300_400 = hist_names_cHel_Mtt300_400[i].c_str();
        const char* name_cHel_Mtt300_400_betaLT0p9 = hist_names_cHel_Mtt300_400_betaLT0p9[i].c_str();
        const char* name_cHel_P3n_Mtt800_Inf = hist_names_cHel_P3n_Mtt800_Inf[i].c_str();
        const char* name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 = hist_names_cHel_P3n_Mtt800_Inf_cosThetaLT0p4[i].c_str();

        if (debug)cout <<" about to fill histos" <<endl;

        hist(name)->Fill(deltay,weight * pdf_weight / orig_weight);
        if (debug)cout <<" done with dy" <<endl;

        if(pt_hadTop < pt_hadTop_thresh && dphi < 0){
          hist(name_dy_d2)->Fill(deltay,weight * pdf_weight / orig_weight);
          if (debug)cout <<" done with dy 2" <<endl;
        }
        if(pt_hadTop < pt_hadTop_thresh && dphi >0){
          hist(name_dy_d1)->Fill(deltay,weight * pdf_weight / orig_weight);
          if (debug)cout <<" done with dy 1" <<endl;
        }
        if (pt_hadTop > pt_hadTop_thresh && deltay <0){
          hist(name_sigma_2)->Fill(sphi,weight * pdf_weight / orig_weight);
          if (debug)cout <<" done with sigma 2" <<endl;
        }
        if (pt_hadTop > pt_hadTop_thresh && deltay >0){
          hist(name_sigma_1)->Fill(sphi,weight * pdf_weight / orig_weight);
          if (debug)cout <<" done with sigma 1" <<endl;
        }
        
        hist(name_cos_theta1k_antiLep)->Fill(cosTheta1k_antiLep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta1r_antiLep)->Fill(cosTheta1r_antiLep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta1n_antiLep)->Fill(cosTheta1n_antiLep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta1kStar_antiLep)->Fill(cosTheta1kStar_antiLep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta1rStar_antiLep)->Fill(cosTheta1rStar_antiLep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2k_Lep)->Fill(cosTheta2k_Lep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2r_Lep)->Fill(cosTheta2r_Lep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2n_Lep)->Fill(cosTheta2n_Lep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2kStar_Lep)->Fill(cosTheta2kStar_Lep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2rStar_Lep)->Fill(cosTheta2rStar_Lep, weight * pdf_weight / orig_weight);
        hist(name_cos_theta1k)->Fill(cosTheta1k, weight * pdf_weight / orig_weight);
        hist(name_cos_theta1r)->Fill(cosTheta1r, weight * pdf_weight / orig_weight);
        hist(name_cos_theta1n)->Fill(cosTheta1n, weight * pdf_weight / orig_weight);
        hist(name_cos_theta1kStar)->Fill(cosTheta1kStar, weight * pdf_weight / orig_weight);
        hist(name_cos_theta1rStar)->Fill(cosTheta1rStar, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2k)->Fill(cosTheta2k, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2r)->Fill(cosTheta2r, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2n)->Fill(cosTheta2n, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2kStar)->Fill(cosTheta2kStar, weight * pdf_weight / orig_weight);
        hist(name_cos_theta2rStar)->Fill(cosTheta2rStar, weight * pdf_weight / orig_weight);
        hist(name_Ckk)->Fill(C_kk, weight * pdf_weight / orig_weight);
        hist(name_Crr)->Fill(C_rr, weight * pdf_weight / orig_weight);
        hist(name_Cnn)->Fill(C_nn, weight * pdf_weight / orig_weight);
        hist(name_Crk_plus)->Fill(C_rk_plus, weight * pdf_weight / orig_weight);
        hist(name_Crk_minus)->Fill(C_rk_minus, weight * pdf_weight / orig_weight);
        hist(name_Cnr_plus)->Fill(C_nr_plus, weight * pdf_weight / orig_weight);
        hist(name_Cnr_minus)->Fill(C_nr_minus, weight * pdf_weight / orig_weight);
        hist(name_Cnk_plus)->Fill(C_nk_plus, weight * pdf_weight / orig_weight);
        hist(name_Cnk_minus)->Fill(C_nk_minus, weight * pdf_weight / orig_weight);
        if(CHel_Mtt300_400 != 99.){
          hist(name_cHel_Mtt300_400)->Fill(CHel_Mtt300_400,weight * pdf_weight / orig_weight);
        }
        if(CHel_Mtt300_400_betaLT0p9 != 99.){
          hist(name_cHel_Mtt300_400_betaLT0p9)->Fill(CHel_Mtt300_400_betaLT0p9,weight * pdf_weight / orig_weight);
        }
        if(CHel_P3n_Mtt800_Inf != 99.){
          hist(name_cHel_P3n_Mtt800_Inf)->Fill(CHel_P3n_Mtt800_Inf,weight * pdf_weight / orig_weight);
        }
        if(CHel_P3n_Mtt800_Inf_cosThetaLT0p4 != 99.){
          hist(name_cHel_P3n_Mtt800_Inf_cosThetaLT0p4)->Fill(CHel_P3n_Mtt800_Inf_cosThetaLT0p4,weight * pdf_weight / orig_weight);
        }

      }
    }// Fill PDF systematic histograms

  }// is_zprime_reconstructed_chi2 && is_mc
}

ZprimeSemiLeptonicPDFHists::~ZprimeSemiLeptonicPDFHists(){}