#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSystematicsModule.h>

#include <stdexcept>

using namespace std;
using namespace uhh2;

ZprimeSemiLeptonicSystematicsModule::Folder
ZprimeSemiLeptonicSystematicsModule::book(Context& ctx, const string& variation, const vector<string>& observables){
  return Folder(new ZprimeSemiLeptonicSpinCorrSystHists(ctx, region_ + "_SC_" + variation, observables));
}

ZprimeSemiLeptonicSystematicsModule::ZprimeSemiLeptonicSystematicsModule(Context& ctx, const string& region,
                                                                         const vector<string>& observables,
                                                                         const map<string, string>& nosc_maps,
                                                                         bool fill_nosc):
  region_(region), fill_nosc_(fill_nosc){

  is_mc = (ctx.get("dataset_type") == "MC");
  const string channel = ctx.get("channel");
  const string version = ctx.get("dataset_version");

  // nominal: MC and data
  nominal = book(ctx, "nominal", observables);
  if(!is_mc) return; // data: no variations

  /////////////////////////////////////////////////////////////////////////////////////////
  //////////////////// variations given as varied / nominal weight ////////////////////////
  // {name, nominal branch, up branch, down branch}
  vector<vector<string>> ratio_table = {
    {"pu",            "weight_pu",                "weight_pu_up",                    "weight_pu_down"},
    {"btag_cferr1",   "weight_btagdisc_central",  "weight_btagdisc_cferr1_up",       "weight_btagdisc_cferr1_down"},
    {"btag_cferr2",   "weight_btagdisc_central",  "weight_btagdisc_cferr2_up",       "weight_btagdisc_cferr2_down"},
    {"btag_hf",       "weight_btagdisc_central",  "weight_btagdisc_hf_up",           "weight_btagdisc_hf_down"},
    {"btag_hfstats1", "weight_btagdisc_central",  "weight_btagdisc_hfstats1_up",     "weight_btagdisc_hfstats1_down"},
    {"btag_hfstats2", "weight_btagdisc_central",  "weight_btagdisc_hfstats2_up",     "weight_btagdisc_hfstats2_down"},
    {"btag_lf",       "weight_btagdisc_central",  "weight_btagdisc_lf_up",           "weight_btagdisc_lf_down"},
    {"btag_lfstats1", "weight_btagdisc_central",  "weight_btagdisc_lfstats1_up",     "weight_btagdisc_lfstats1_down"},
    {"btag_lfstats2", "weight_btagdisc_central",  "weight_btagdisc_lfstats2_up",     "weight_btagdisc_lfstats2_down"},
    {"ttag_corr",     "weight_toptagsf",          "weight_toptagsf_corr_up",         "weight_toptagsf_corr_down"},
    {"ttag_uncorr",   "weight_toptagsf",          "weight_toptagsf_uncorr_up",       "weight_toptagsf_uncorr_down"},
    {"tmistag",       "weight_topmistagsf",       "weight_topmistagsf_up",           "weight_topmistagsf_down"},
    {"nnloQCD",       "weight_ttbar_nnlo_qcd",    "weight_ttbar_nnlo_qcd_up",        "weight_ttbar_nnlo_qcd_down"},
    {"ewk",           "weight_ttbar_ewk_nominal", "weight_ttbar_ewk_up",             "weight_ttbar_ewk_down"},
  };
  if(channel == "muon"){
    ratio_table.insert(ratio_table.end(), {
      {"mu_id_stat",      "weight_sfmu_id_stat",      "weight_sfmu_id_stat_up",      "weight_sfmu_id_stat_down"},
      {"mu_id_syst",      "weight_sfmu_id_syst",      "weight_sfmu_id_syst_up",      "weight_sfmu_id_syst_down"},
      {"mu_iso_stat",     "weight_sfmu_iso_stat",     "weight_sfmu_iso_stat_up",     "weight_sfmu_iso_stat_down"},
      {"mu_iso_syst",     "weight_sfmu_iso_syst",     "weight_sfmu_iso_syst_up",     "weight_sfmu_iso_syst_down"},
      {"mu_reco",         "weight_sfmu_reco",         "weight_sfmu_reco_up",         "weight_sfmu_reco_down"},
      {"mu_trigger_stat", "weight_sfmu_trigger_stat", "weight_sfmu_trigger_stat_up", "weight_sfmu_trigger_stat_down"},
      {"mu_trigger_syst", "weight_sfmu_trigger_syst", "weight_sfmu_trigger_syst_up", "weight_sfmu_trigger_syst_down"},
    });
  }
  if(channel == "electron"){
    ratio_table.insert(ratio_table.end(), {
      {"ele_id",      "weight_sfelec_id",      "weight_sfelec_id_up",      "weight_sfelec_id_down"},
      {"ele_reco",    "weight_sfelec_reco",    "weight_sfelec_reco_up",    "weight_sfelec_reco_down"},
      {"ele_trigger", "weight_sfelec_trigger", "weight_sfelec_trigger_up", "weight_sfelec_trigger_down"},
    });
  }
  for(const auto& row : ratio_table){
    RatioSyst syst;
    syst.name      = row[0];
    syst.h_nominal = ctx.get_handle<float>(row[1]);
    syst.h_up      = ctx.get_handle<float>(row[2]);
    syst.h_down    = ctx.get_handle<float>(row[3]);
    syst.up        = book(ctx, row[0] + "_up",   observables);
    syst.down      = book(ctx, row[0] + "_down", observables);
    ratio_systs.push_back(std::move(syst));
  }

  // L1 prefiring: weights are members of the Event (prefiringWeight, prefiringWeightUp/Down)
  prefiring_up   = book(ctx, "prefiring_up",   observables);
  prefiring_down = book(ctx, "prefiring_down", observables);

  /////////////////////////////////////////////////////////////////////////////////////////
  ////////////////// variations stored relative to the nominal weight /////////////////////
  // {folder name, branch}
  const vector<pair<string, string>> relative_table = {
    {"isr_up",          "weight_isr_2_up"},          // PS weights: ISR/FSR scale x2 / x0.5
    {"isr_down",        "weight_isr_2_down"},
    {"fsr_up",          "weight_fsr_2_up"},
    {"fsr_down",        "weight_fsr_2_down"},
    {"murmuf_upup",     "weight_murmuf_upup"},       // ME scales: (muR, muF) x2 ("up") / x0.5 ("down") / x1 ("none")
    {"murmuf_upnone",   "weight_murmuf_upnone"},
    {"murmuf_noneup",   "weight_murmuf_noneup"},
    {"murmuf_nonedown", "weight_murmuf_nonedown"},
    {"murmuf_downnone", "weight_murmuf_downnone"},
    {"murmuf_downdown", "weight_murmuf_downdown"},
  };
  for(const auto& row : relative_table){
    RelativeSyst syst;
    syst.name     = row.first;
    syst.h_weight = ctx.get_handle<float>(row.second);
    syst.folder   = book(ctx, row.first, observables);
    relative_systs.push_back(std::move(syst));
  }

  // PDF: 100 NNPDF3.1 replicas; samples generated with another PDF set store them at a different index
  const bool other_pdf_set = (version.find("WJetsToLNu_HT-") != string::npos) || (version.find("DYJetsToLL_M-50_HT-") != string::npos) || (version.find("QCD_HT") != string::npos);
  pdf_start_index = other_pdf_set ? 47 : 9;
  for(int i = 0; i < 100; ++i) pdf.push_back(book(ctx, "pdf" + to_string(i), observables));

  /////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////// NoSC signal-shape variation ///////////////////////////////////
  if(fill_nosc_){
    vector<string> nosc_observables;
    for(const auto& obs : observables){
      const auto it = nosc_maps.find(obs);
      if(it == nosc_maps.end()) continue;
      nosc_observables.push_back(obs);
      h_nosc_up.push_back(ctx.get_handle<float>("weight_NoSC_" + it->second + "_up"));
      h_nosc_down.push_back(ctx.get_handle<float>("weight_NoSC_" + it->second + "_down"));
    }
    if(nosc_observables.empty()) fill_nosc_ = false;
    else{
      nosc_up   = book(ctx, "NoSC_up",   nosc_observables);
      nosc_down = book(ctx, "NoSC_down", nosc_observables);
    }
  }
}


void ZprimeSemiLeptonicSystematicsModule::fill(Event& event){

  const double weight_nominal = event.weight;
  nominal->fill(event);
  if(!is_mc) return;

  // varied / nominal weight (a nominal weight of 0 means the event already has weight 0)
  for(auto& syst : ratio_systs){
    const double w_nom = event.get(syst.h_nominal);
    event.weight = (w_nom != 0.) ? weight_nominal * event.get(syst.h_up)   / w_nom : 0.;
    syst.up->fill(event);
    event.weight = (w_nom != 0.) ? weight_nominal * event.get(syst.h_down) / w_nom : 0.;
    syst.down->fill(event);
  }

  // L1 prefiring
  const double w_pref = event.prefiringWeight;
  event.weight = (w_pref != 0.) ? weight_nominal * event.prefiringWeightUp   / w_pref : 0.;
  prefiring_up->fill(event);
  event.weight = (w_pref != 0.) ? weight_nominal * event.prefiringWeightDown / w_pref : 0.;
  prefiring_down->fill(event);

  // weights stored relative to nominal
  for(auto& syst : relative_systs){
    event.weight = weight_nominal * event.get(syst.h_weight);
    syst.folder->fill(event);
  }

  // PDF replicas; events without the replica weights keep the nominal weight (no PDF effect)
  const bool has_pdf = event.genInfo && event.genInfo->systweights().size() > (size_t)(pdf_start_index + 100)
                       && event.genInfo->originalXWGTUP() != 0.;
  for(int i = 0; i < 100; ++i){
    event.weight = has_pdf ? weight_nominal * event.genInfo->systweights().at(pdf_start_index + i) / event.genInfo->originalXWGTUP()
                           : weight_nominal;
    pdf[i]->fill(event);
  }

  // NoSC: each observable with its own per-event weight
  event.weight = weight_nominal;
  if(fill_nosc_){
    vector<double> factors_up, factors_down;
    for(size_t i = 0; i < h_nosc_up.size(); ++i){
      factors_up.push_back(event.get(h_nosc_up[i]));
      factors_down.push_back(event.get(h_nosc_down[i]));
    }
    nosc_up->fill_with_factors(event, factors_up);
    nosc_down->fill_with_factors(event, factors_down);
  }

  event.weight = weight_nominal;
}
