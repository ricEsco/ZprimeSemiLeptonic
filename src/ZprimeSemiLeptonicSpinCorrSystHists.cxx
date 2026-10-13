#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSpinCorrSystHists.h"

using namespace std;
using namespace uhh2;

namespace {
  struct VarSpec {
    const char* handle_name; // branch name declared by SpinCorrelations (ctx.declare_event_output<float>(...))
    const char* hist_name;   // histogram name (matches ZprimeSemiLeptonicHists naming convention)
    const char* title;
    int nbins;
    double xmin, xmax;
  };

  const VarSpec kVars[] = {
    // spin-analyzer (leptons-only) projections
    {"cosTheta1k_antiLep",      "cos_theta1k_antiLep",      "cos(#theta_{antilep}^{k})",   24, -1, 1},
    {"cosTheta1r_antiLep",      "cos_theta1r_antiLep",      "cos(#theta_{antilep}^{r})",   24, -1, 1},
    {"cosTheta1n_antiLep",      "cos_theta1n_antiLep",      "cos(#theta_{antilep}^{n})",   24, -1, 1},
    {"cosTheta1kStar_antiLep",  "cos_theta1kStar_antiLep",  "cos(#theta_{antilep}^{k*})",  24, -1, 1},
    {"cosTheta1rStar_antiLep",  "cos_theta1rStar_antiLep",  "cos(#theta_{antilep}^{r*})",  24, -1, 1},
    {"cosTheta2k_Lep",          "cos_theta2k_Lep",          "cos(#theta_{lep}^{k})",       24, -1, 1},
    {"cosTheta2r_Lep",          "cos_theta2r_Lep",          "cos(#theta_{lep}^{r})",       24, -1, 1},
    {"cosTheta2n_Lep",          "cos_theta2n_Lep",          "cos(#theta_{lep}^{n})",       24, -1, 1},
    {"cosTheta2kStar_Lep",      "cos_theta2kStar_Lep",      "cos(#theta_{lep}^{k*})",      24, -1, 1},
    {"cosTheta2rStar_Lep",      "cos_theta2rStar_Lep",      "cos(#theta_{lep}^{r*})",      24, -1, 1},
    // spin-analyzer projections
    {"cosTheta1k",              "cos_theta1k",              "cos(#theta_{1}^{k})",         24, -1, 1},
    {"cosTheta1r",              "cos_theta1r",              "cos(#theta_{1}^{r})",         24, -1, 1},
    {"cosTheta1n",              "cos_theta1n",              "cos(#theta_{1}^{n})",         24, -1, 1},
    {"cosTheta1kStar",          "cos_theta1kStar",          "cos(#theta_{1}^{k*})",        24, -1, 1},
    {"cosTheta1rStar",          "cos_theta1rStar",          "cos(#theta_{1}^{r*})",        24, -1, 1},
    {"cosTheta2k",              "cos_theta2k",              "cos(#theta_{2}^{k})",         24, -1, 1},
    {"cosTheta2r",              "cos_theta2r",              "cos(#theta_{2}^{r})",         24, -1, 1},
    {"cosTheta2n",              "cos_theta2n",              "cos(#theta_{2}^{n})",         24, -1, 1},
    {"cosTheta2kStar",          "cos_theta2kStar",          "cos(#theta_{2}^{k*})",        24, -1, 1},
    {"cosTheta2rStar",          "cos_theta2rStar",          "cos(#theta_{2}^{r*})",        24, -1, 1},
    // Correlation elements
    {"Cnn", "Cnn", "C_{nn}", 24, -1, 1},
    {"Cnr", "Cnr", "C_{nr}", 24, -1, 1},
    {"Cnk", "Cnk", "C_{nk}", 24, -1, 1},
    {"Crn", "Crn", "C_{rn}", 24, -1, 1},
    {"Crr", "Crr", "C_{rr}", 24, -1, 1},
    {"Crk", "Crk", "C_{rk}", 24, -1, 1},
    {"Ckn", "Ckn", "C_{kn}", 24, -1, 1},
    {"Ckr", "Ckr", "C_{kr}", 24, -1, 1},
    {"Ckk", "Ckk", "C_{kk}", 24, -1, 1},
    // linear combinations
    {"Crk_plus",  "Crk_plus",  "C_{rk} + C_{kr}", 24, -1, 1},
    {"Crk_minus", "Crk_minus", "C_{rk} - C_{kr}", 24, -1, 1},
    {"Cnr_plus",  "Cnr_plus",  "C_{nr} + C_{rn}", 24, -1, 1},
    {"Cnr_minus", "Cnr_minus", "C_{nr} - C_{rn}", 24, -1, 1},
    {"Cnk_plus",  "Cnk_plus",  "C_{nk} + C_{kn}", 24, -1, 1},
    {"Cnk_minus", "Cnk_minus", "C_{nk} - C_{kn}", 24, -1, 1},
    // Entanglement witnesses
    {"cHel",                                "cHel",                                "cos(#phi_{lb})", 24, -1, 1},
    {"cHel_Mtt300_400",                     "cHel_Mtt300_400",                     "cos(#phi_{lb}) (M_{tt} [300,400] GeV)", 24, -1, 1},
    {"cHel_Mtt300_400_betaLT0p9",           "cHel_Mtt300_400_betaLT0p9",           "cos(#phi_{lb}) (M_{tt} [300,400] GeV, #beta_{tt} < 0.9)", 24, -1, 1},
    {"cHel_P3n",                            "cHel_P3n",                            "cos(#phi_{(P3n)lb})", 24, -1, 1},
    {"cHel_P3n_Mtt800_Inf",                 "cHel_P3n_Mtt800_Inf",                 "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV)", 24, -1, 1},
    {"cHel_P3n_Mtt800_Inf_cosThetaLT0p4",   "cHel_P3n_Mtt800_Inf_cosThetaLT0p4",   "cos(#phi_{(P3n)lb}) (M_{tt} > 800 GeV, |cos(#theta)| < 0.4)", 24, -1, 1},
    // Baumgart et al. variables
    {"Sigma_phi", "Sigma_phi", "#Sigma #phi ", 16, -3.2, 3.2},
    {"Delta_phi", "Delta_phi", "#Delta #phi ", 16, -3.2, 3.2},
  };
  const int kNVars = sizeof(kVars) / sizeof(kVars[0]);
}

ZprimeSemiLeptonicSpinCorrSystHists::ZprimeSemiLeptonicSpinCorrSystHists(uhh2::Context& ctx, const std::string& dirname):
Hists(ctx, dirname) {
  for(int i = 0; i < kNVars; ++i){
    var_names_.push_back(kVars[i].handle_name);
    var_handles_.push_back(ctx.get_handle<float>(kVars[i].handle_name));
  }
  init();
}

void ZprimeSemiLeptonicSpinCorrSystHists::init(){
  for(int i = 0; i < kNVars; ++i){
    var_hists_.push_back(book<TH1F>(kVars[i].hist_name, kVars[i].title, kVars[i].nbins, kVars[i].xmin, kVars[i].xmax));
  }
}

void ZprimeSemiLeptonicSpinCorrSystHists::fill(const uhh2::Event& event){
  double weight = event.weight;
  for(size_t i = 0; i < var_handles_.size(); ++i){
    var_hists_[i]->Fill(event.get(var_handles_[i]), weight);
  }
}
