#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include <string>
#include <vector>

#include <TH1F.h>

// Slim Hists class used ONLY for booking/filling systematic-variation folders (e.g. "Merged_pdf3", "Chi2cut_PASS_mu_id_stat_up").
// It books entanglement/spin-correlation observables computed by SpinCorrelations (ZprimeSemiLeptonicModules.cxx),  -- 43 variables total. 
// Values are read back via the event handles SpinCorrelations already declared with ctx.declare_event_output<float>(...)
// this class only re-fills the already-computed value into a differently-weighted histogram. 

class ZprimeSemiLeptonicSpinCorrSystHists: public uhh2::Hists {
public:
  explicit ZprimeSemiLeptonicSpinCorrSystHists(uhh2::Context&, const std::string&);
  virtual void fill(const uhh2::Event&) override;

protected:
  void init();

  std::vector<std::string> var_names_;
  std::vector<uhh2::Event::Handle<float>> var_handles_;
  std::vector<TH1F*> var_hists_;
};
