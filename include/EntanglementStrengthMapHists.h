#pragma once

#include <UHH2/core/include/Hists.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/common/include/TTbarGen.h>

#include <TProfile2D.h>

#include <map>
#include <string>

/*
Inclusive gen-level alpha strength maps for the template reweighting method (docs/reweighting_method.md).

One TProfile2D per spin-correlation observable, binned in gen-level (m_tt, cos theta*):
    alpha(m_tt, cos theta*) = prefactor * < sign(Obs) >_genweight
  - each event fills prefactor * sign(Obs) with weight = generator weight (originalXWGTUP)
  - the TProfile2D bin content is then the gen-weighted mean = prefactor * asymmetry = alpha
  - prefactor converts the asymmetry into the coefficient of the observable's density:
      2 for polarizations B_i and the witnesses D, D-tilde;  4 for the correlation elements C_ij

Filled in ZprimePreselectionModule on the INCLUSIVE ntuples (before any selection), so the maps describe the full
gen-level phase space. Observables come from entanglement::compute_gen_spin_correlations (last-copy tops), the same
function that writes the gen_<Obs> branches in the AnalysisTree, so map and reweighting observables are identical.

Output: <dirname>/alpha_<Obs> in the Preselection output file (hadd-safe: TProfile2D merges correctly).
*/
class EntanglementStrengthMapHists: public uhh2::Hists {
public:
  explicit EntanglementStrengthMapHists(uhh2::Context&, const std::string& dirname);
  virtual void fill(const uhh2::Event&) override;

private:
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  std::map<std::string, std::pair<TProfile2D*, double>> maps; // observable name -> (map, prefactor)
};
