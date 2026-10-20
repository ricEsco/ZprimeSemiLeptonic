#pragma once

#include <UHH2/common/include/TTbarGen.h>

#include <TLorentzVector.h>

#include <map>
#include <string>
#include <vector>

/*
Generator-level spin-correlation observables of the semileptonic ttbar system.

This is the single C++ definition of the gen-level observables used by the analysis:
  - ZprimePreselectionModule fills the inclusive alpha strength maps with them (EntanglementStrengthMapHists)
  - ZprimeAnalysisModule writes them to the AnalysisTree as gen_<name> branches
Using one function for both guarantees that the map observable and the per-event observable used in the
reweighting are identical by construction.

Recipe (same as SpinCorrelations in ZprimeSemiLeptonicModules.cxx, but with gen-particles):
  1) top/antitop = LAST-COPY tops, built as W + b from TTbarGen (the direct parents of the decay products,
     so the spin analyzers are boosted into the rest frame of their actual parent top)
  2) boost everything into the ttbar CM frame and build the Bose-symmetrized Bernreuther (k, r, n) basis
  3) boost the spin analyzers (charged lepton, hadronic b-quark) into their parent top rest frame
  4) project the spin-analyzer directions onto the basis -> polarizations, C_ij, entanglement witnesses, ...

Sentinel conventions (mirror SpinCorrelations exactly):
  -10 : observable not available this event (not semileptonic, or outside a windowed slice)
   99 : (anti)lepton-only projection when that (anti)lepton is not present in the event
*/

namespace entanglement {

  constexpr float SENTINEL    = -10.f;   // observable not available this event
  constexpr float NOT_PRESENT =  99.f;   // (anti)lepton-only projection where that (anti)lepton is absent
  constexpr float PT_HADTOP_THRESHOLD = 150.f; // [GeV] low/high pT(hadronic top) split of the Baumgart variables

  // Result of the gen-level computation for one event
  struct GenSpinCorrelations {
    bool is_semileptonic = false;   // TTbarGen decay is e/mu/tau + jets
    bool obs_valid       = false;   // observables computed (false if not semileptonic or top parallel to beam)

    // lab-frame 4-vectors (filled whenever is_semileptonic)
    TLorentzVector Top, Antitop;    // last-copy top/antitop (W + b)
    TLorentzVector bHad, Lepton;    // hadronic-side b-quark and charged lepton from the W decay

    // observable name -> value, keys are gen_observable_names(); all SENTINEL unless obs_valid
    std::map<std::string, float> obs;
  };

  // Ordered list of observable names (the gen_ branch names without prefix)
  const std::vector<std::string>& gen_observable_names();

  // Compute all gen-level observables from the TTbarGen object
  GenSpinCorrelations compute_gen_spin_correlations(const TTbarGen& ttbargen);

}
