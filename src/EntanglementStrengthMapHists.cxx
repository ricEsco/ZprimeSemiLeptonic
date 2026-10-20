#include <UHH2/ZprimeSemiLeptonic/include/EntanglementStrengthMapHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/EntanglementGenObservables.h>

#include <cmath>

using namespace std;
using namespace uhh2;

EntanglementStrengthMapHists::EntanglementStrengthMapHists(Context& ctx, const string& dirname): Hists(ctx, dirname){

  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen"); // filled by TTbarGenProducer in the calling module

  // Observables with a strength map and their prefactors
  const vector<pair<string, double>> observables = {
    {"cosTheta1k", 2.}, {"cosTheta1r", 2.}, {"cosTheta1n", 2.},    // B1_i (top polarization)
    {"cosTheta2k", 2.}, {"cosTheta2r", 2.}, {"cosTheta2n", 2.},    // B2_i (antitop polarization)
    {"Cnn", 4.}, {"Cnr", 4.}, {"Cnk", 4.},                         // C_ij (spin correlation matrix)
    {"Crn", 4.}, {"Crr", 4.}, {"Crk", 4.},
    {"Ckn", 4.}, {"Ckr", 4.}, {"Ckk", 4.},
    {"cHel", 2.},                                                  // D       (threshold entanglement witness)
    {"cHel_P3n", 2.},                                              // D-tilde (boosted entanglement witness)
  };

  // Binning: 52 bins of ~32.7 GeV in m_tt [300, 2000] x 48 bins in cos theta* [-1, 1]
  const int    n_mtt = 52;   const double mtt_lo = 300.,  mtt_hi = 2000.;
  const int    n_cts = 48;   const double cts_lo = -1.,   cts_hi = 1.;

  for(const auto& obs : observables){
    const string name  = "alpha_" + obs.first;
    const string title = "#alpha " + obs.first + ";m_{t#bar{t}}^{gen} [GeV];cos#theta*^{gen};#alpha";
    TProfile2D* map = book<TProfile2D>(name, title, n_mtt, mtt_lo, mtt_hi, n_cts, cts_lo, cts_hi);
    maps[obs.first] = make_pair(map, obs.second);
  }
}


void EntanglementStrengthMapHists::fill(const Event& event){

  if(event.isRealData || !event.genInfo || !event.is_valid(h_ttbargen)) return;

  // Gen-level observables (last-copy tops); skip non-semileptonic events and top parallel to the beam
  const auto gen = entanglement::compute_gen_spin_correlations(event.get(h_ttbargen));
  if(!gen.obs_valid) return;

  const double mtt     = gen.obs.at("M_tt");
  const double cosTS   = gen.obs.at("cosThetaStar");
  const double genweight = event.genInfo->originalXWGTUP(); // nominal generator weight only

  // Fill prefactor * sign(Obs); the per-cell gen-weighted mean is alpha
  for(const auto& entry : maps){
    const double value = gen.obs.at(entry.first);
    if(!std::isfinite(value)) continue;
    const double sign = (value >= 0.) ? 1. : -1.;
    entry.second.first->Fill(mtt, cosTS, entry.second.second * sign, genweight);
  }
}
