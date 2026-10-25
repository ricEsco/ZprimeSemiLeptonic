#include <UHH2/ZprimeSemiLeptonic/include/EntanglementModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/EntanglementGenObservables.h>

#include <TFile.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace uhh2;


/*
███████ ███████ ███    ██      ██████  ██    ██ ████████ ██████  ██    ██ ████████
██      ██      ████   ██     ██    ██ ██    ██    ██    ██   ██ ██    ██    ██
██  ███ █████   ██ ██  ██     ██    ██ ██    ██    ██    ██████  ██    ██    ██
██   ██ ██      ██  ██ ██     ██    ██ ██    ██    ██    ██      ██    ██    ██
 ██████ ███████ ██   ████      ██████   ██████     ██    ██       ██████     ██
*/

namespace {
  const vector<string> kGenParticles = {"Top", "Antitop", "bHad", "Lepton"};
  const vector<string> kP4Components = {"pt", "eta", "phi", "E"};
}

EntanglementGenOutput::EntanglementGenOutput(Context& ctx, const string& ttgen_name){
  is_mc = (ctx.get("dataset_type") == "MC");
  h_ttbargen = ctx.get_handle<TTbarGen>(ttgen_name);

  // one branch per gen observable, same names as the reco branches with a gen_ prefix
  for(const auto& name : entanglement::gen_observable_names()){
    h_obs.emplace_back(name, ctx.declare_event_output<float>("gen_" + name));
  }
  // lab-frame 4-vectors (last-copy tops)
  for(const auto& p : kGenParticles){
    for(const auto& c : kP4Components) h_p4.push_back(ctx.declare_event_output<float>("gen_" + p + "_" + c));
  }
  h_semilep   = ctx.declare_event_output<int>("gen_semilep");
  h_channel   = ctx.declare_event_output<int>("gen_channel");
  h_genweight = ctx.declare_event_output<float>("genweight");
}

bool EntanglementGenOutput::process(Event& event){

  // generator weight and decay channel (data: neutral values)
  event.set(h_genweight, (is_mc && event.genInfo) ? event.genInfo->originalXWGTUP() : 1.f);
  const bool has_ttbargen = is_mc && event.is_valid(h_ttbargen);
  event.set(h_channel, has_ttbargen ? static_cast<int>(event.get(h_ttbargen).DecayChannel()) : -1);

  // gen observables: SENTINEL (-10) unless the event is a semileptonic ttbar decay
  entanglement::GenSpinCorrelations gen;
  if(has_ttbargen) gen = entanglement::compute_gen_spin_correlations(event.get(h_ttbargen));

  event.set(h_semilep, gen.obs_valid ? 1 : 0);
  for(const auto& entry : h_obs){
    event.set(entry.second, gen.obs_valid ? gen.obs.at(entry.first) : entanglement::SENTINEL);
  }

  // 4-vectors: filled whenever the decay is semileptonic
  const vector<const TLorentzVector*> p4 = {&gen.Top, &gen.Antitop, &gen.bHad, &gen.Lepton};
  for(size_t i = 0; i < p4.size(); ++i){
    const bool ok = gen.is_semileptonic;
    event.set(h_p4[4*i + 0], ok ? p4[i]->Pt()  : entanglement::SENTINEL);
    event.set(h_p4[4*i + 1], ok ? p4[i]->Eta() : entanglement::SENTINEL);
    event.set(h_p4[4*i + 2], ok ? p4[i]->Phi() : entanglement::SENTINEL);
    event.set(h_p4[4*i + 3], ok ? p4[i]->E()   : entanglement::SENTINEL);
  }
  return true;
}


/*
███    ██  ██████  ███████  ██████     ██     ██ ███████ ██  ██████  ██   ██ ████████ ███████
████   ██ ██    ██ ██      ██          ██     ██ ██      ██ ██       ██   ██    ██    ██
██ ██  ██ ██    ██ ███████ ██          ██  █  ██ █████   ██ ██   ███ ███████    ██    ███████
██  ██ ██ ██    ██      ██ ██          ██ ███ ██ ██      ██ ██    ██ ██   ██    ██         ██
██   ████  ██████  ███████  ██████      ███ ███  ███████ ██  ██████  ██   ██    ██    ███████
*/

namespace {
  const double kF_up   = +1.;   // NoSC up:   spin correlation removed
  const double kF_down = -1.;   // NoSC down: spin correlation enhanced
  const double kDensityFloor = 1e-6; // guard on |1 + alpha*O|
}

const vector<string>& EntanglementNoSCReweighting::mapped_observables(){
  // same list as EntanglementStrengthMapHists (Preselection)
  static const vector<string> names = {
    "cosTheta1k", "cosTheta1r", "cosTheta1n", "cosTheta2k", "cosTheta2r", "cosTheta2n", // B1_i, B2_i
    "Cnn", "Cnr", "Cnk", "Crn", "Crr", "Crk", "Ckn", "Ckr", "Ckk",                      // C_ij
    "cHel", "cHel_P3n",                                                                 // D, D-tilde
  };
  return names;
}

EntanglementNoSCReweighting::EntanglementNoSCReweighting(Context& ctx){

  // Active only for the nominal semileptonic ttbar sample of this era
  const bool is_mc = (ctx.get("dataset_type") == "MC");
  const string version = ctx.get("dataset_version");
  string year = "";
  for(const string era : {"UL16preVFP", "UL16postVFP", "UL17", "UL18"}){
    if(version.find(era) != string::npos){ year = era; break; }
  }
  active = is_mc && (version == "TTToSemiLeptonic_" + year);

  h_gen_mtt     = ctx.get_handle<float>("gen_M_tt");
  h_gen_cosTS   = ctx.get_handle<float>("gen_cosThetaStar");
  h_gen_semilep = ctx.get_handle<int>("gen_semilep");

  // Load the maps (only needed for the signal sample)
  unique_ptr<TFile> map_file;
  const string folder = ctx.get("StrengthMapFolder", "StrengthMaps");
  if(active){
    const string path = ctx.get("StrengthMapFile", "");
    if(path.empty()) throw runtime_error("EntanglementNoSCReweighting: config item 'StrengthMapFile' is required for " + version);
    map_file.reset(TFile::Open(path.c_str()));
    if(!map_file || map_file->IsZombie()) throw runtime_error("EntanglementNoSCReweighting: cannot open StrengthMapFile " + path);
  }

  for(const auto& obs : mapped_observables()){
    MapEntry entry;
    entry.obs       = obs;
    entry.h_gen_obs = ctx.get_handle<float>("gen_" + obs);
    entry.h_alpha   = ctx.declare_event_output<float>("alpha_" + obs);
    entry.h_up      = ctx.declare_event_output<float>("weight_NoSC_" + obs + "_up");
    entry.h_down    = ctx.declare_event_output<float>("weight_NoSC_" + obs + "_down");
    if(active){
      TProfile2D* map = dynamic_cast<TProfile2D*>(map_file->Get((folder + "/alpha_" + obs).c_str()));
      if(!map) throw runtime_error("EntanglementNoSCReweighting: map " + folder + "/alpha_" + obs + " not found in StrengthMapFile");
      entry.map.reset(static_cast<TProfile2D*>(map->Clone(("NoSC_alpha_" + obs).c_str())));
      entry.map->SetDirectory(nullptr);
    }
    entries.push_back(std::move(entry));
  }
  if(map_file) map_file->Close();
}

// alpha at (m_tt, cos theta*), clamped into the map range so very boosted / edge events use the nearest valid cell
double EntanglementNoSCReweighting::lookup_alpha(TProfile2D& map, double mtt, double cosTS) const {
  const TAxis* ax = map.GetXaxis();
  const TAxis* ay = map.GetYaxis();
  const double eps_x = 1e-3 * (ax->GetXmax() - ax->GetXmin()) / ax->GetNbins();
  const double eps_y = 1e-3 * (ay->GetXmax() - ay->GetXmin()) / ay->GetNbins();
  const double x = std::min(std::max(mtt,   ax->GetXmin() + eps_x), ax->GetXmax() - eps_x);
  const double y = std::min(std::max(cosTS, ay->GetXmin() + eps_y), ay->GetXmax() - eps_y);
  return map.GetBinContent(map.FindBin(x, y));
}

bool EntanglementNoSCReweighting::process(Event& event){

  const bool reweight = active && event.get(h_gen_semilep) == 1;
  const double mtt   = reweight ? event.get(h_gen_mtt)   : 0.;
  const double cosTS = reweight ? event.get(h_gen_cosTS) : 0.;

  for(auto& entry : entries){
    double alpha = 0., w_up = 1., w_down = 1.;
    if(reweight){
      const double obs = event.get(entry.h_gen_obs);
      alpha = lookup_alpha(*entry.map, mtt, cosTS);
      const double density = 1. + alpha * obs;
      if(std::fabs(density) >= kDensityFloor){
        w_up   = (1. + alpha * obs * (1. - kF_up))   / density;
        w_down = (1. + alpha * obs * (1. - kF_down)) / density;
      }
    }
    event.set(entry.h_alpha, static_cast<float>(alpha));
    event.set(entry.h_up,    static_cast<float>(w_up));
    event.set(entry.h_down,  static_cast<float>(w_down));
  }
  return true;
}
