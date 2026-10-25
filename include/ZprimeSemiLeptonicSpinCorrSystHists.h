#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include <string>
#include <vector>

#include <TH1F.h>

// Slim Hists class for the systematic-variation folders written by ZprimeSemiLeptonicSystematicsModule
// (e.g. "Merged_SC_pdf3", "Chi2cut_PASS_SC_mu_id_stat_up").
// It books entanglement/spin-correlation observables computed by SpinCorrelations (ZprimeSemiLeptonicModules.cxx).
// Values are read back via the event handles SpinCorrelations declared with ctx.declare_event_output<float>(...),
// so this class only re-fills the already-computed value into a differently-weighted histogram.
// Valid only for weight-only systematics (the observable itself does not change).
//
// Unavailable values keep their SpinCorrelations sentinels and land outside [-1, 1]:
//   -10 (not reconstructed / outside a windowed slice) -> underflow, 99 (absent (anti)lepton slot) -> overflow

class ZprimeSemiLeptonicSpinCorrSystHists: public uhh2::Hists {
public:
  // all observables of the internal table
  explicit ZprimeSemiLeptonicSpinCorrSystHists(uhh2::Context&, const std::string&);
  // only the listed observables (SpinCorrelations branch names, e.g. "cHel_P3n_Mtt800_Inf_cosThetaLT0p4")
  ZprimeSemiLeptonicSpinCorrSystHists(uhh2::Context&, const std::string&, const std::vector<std::string>& observables);

  // fill every histogram with event.weight
  virtual void fill(const uhh2::Event&) override;
  // fill histogram i with event.weight * weight_factors[i] (per-observable weights, e.g. NoSC variations)
  void fill_with_factors(const uhh2::Event&, const std::vector<double>& weight_factors);

  const std::vector<std::string>& observables() const { return var_names_; }

protected:
  void add(uhh2::Context&, const std::string& observable);

  std::vector<std::string> var_names_;
  std::vector<uhh2::Event::Handle<float>> var_handles_;
  std::vector<TH1F*> var_hists_;
};
