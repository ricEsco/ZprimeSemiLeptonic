#pragma once

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/common/include/TTbarGen.h>

#include <TProfile2D.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

/*
AnalysisModules that write the gen-level entanglement inputs to the AnalysisTree.
Together with the reco branches of SpinCorrelations they replace the skim files of
src/entanglementScripts/makeAnalysisSkims.py and the per-event part of make_templates_cHel_P3n_LepChannel.py.

  EntanglementGenOutput       : gen_<Obs> branches (last-copy tops, see EntanglementGenObservables.h), gen 4-vectors,
                                gen_semilep, gen_channel, genweight
  EntanglementNoSCReweighting : per-event "NoSC" signal-shape weights from the inclusive alpha strength maps
                                built in ZprimePreselectionModule (StrengthMaps/alpha_<Obs>)

Run EntanglementGenOutput first: EntanglementNoSCReweighting reads the gen_<Obs> branches it sets.
*/


// Gen-level observables and 4-vectors of the semileptonic ttbar system
class EntanglementGenOutput : public uhh2::AnalysisModule {
public:
  explicit EntanglementGenOutput(uhh2::Context& ctx, const std::string& ttgen_name = "ttbargen");
  virtual bool process(uhh2::Event& event) override;

private:
  bool is_mc;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  std::vector<std::pair<std::string, uhh2::Event::Handle<float>>> h_obs; // gen_<Obs>
  std::vector<uhh2::Event::Handle<float>> h_p4;                          // gen_<Top|Antitop|bHad|Lepton>_<pt|eta|phi|E>
  uhh2::Event::Handle<int>   h_semilep;   // 1 if the gen observables are defined (semileptonic decay), else 0
  uhh2::Event::Handle<int>   h_channel;   // TTbarGen decay channel: 0 had, 1-3 semilep (e, mu, tau), 4-9 dilep, 10 not found; -1 data
  uhh2::Event::Handle<float> h_genweight; // nominal generator weight (originalXWGTUP); 1 for data
};


/*
NoSC ("no spin correlation") signal-shape variation of the template reweighting method (docs/reweighting_method.md)
For each mapped observable O, with alpha = alpha_O(gen m_tt, gen cos theta*) looked up in the inclusive strength map:
    density  d   = 1 + alpha * O_gen
    varied   d_f = 1 + alpha * (1 - f) * O_gen
    w(f)         = d_f / d          f = +1 -> NoSC up (spin correlation removed), f = -1 -> NoSC down (enhanced)
Active only for the nominal TTToSemiLeptonic sample; every other sample gets weight 1 (and alpha 0).
Branches per observable: alpha_<Obs>, weight_NoSC_<Obs>_up, weight_NoSC_<Obs>_down

Config items: StrengthMapFile   = hadded Preselection output of TTToSemiLeptonic for this era (contains StrengthMaps/)
              StrengthMapFolder = folder holding the maps (default "StrengthMaps")
*/
class EntanglementNoSCReweighting : public uhh2::AnalysisModule {
public:
  explicit EntanglementNoSCReweighting(uhh2::Context& ctx);
  virtual bool process(uhh2::Event& event) override;

  // observables with a strength map, i.e. with a NoSC variation
  static const std::vector<std::string>& mapped_observables();
  bool is_active() const { return active; }

private:
  struct MapEntry {
    std::string obs;
    std::unique_ptr<TProfile2D> map;
    uhh2::Event::Handle<float> h_gen_obs, h_alpha, h_up, h_down;
  };

  double lookup_alpha(TProfile2D& map, double mtt, double cosTS) const; // non-const map: ROOT FindBin is non-const

  bool active;
  std::vector<MapEntry> entries;
  uhh2::Event::Handle<float> h_gen_mtt, h_gen_cosTS;
  uhh2::Event::Handle<int>   h_gen_semilep;
};
