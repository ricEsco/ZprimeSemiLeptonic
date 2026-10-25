#pragma once

#include <UHH2/core/include/Event.h>
#include <UHH2/ZprimeSemiLeptonic/include/ElecTriggerSF.h>
#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSpinCorrSystHists.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

/*
Weight-only systematic variations of the spin-correlation observables in one analysis region.
Owned by ZprimeAnalysisModule (one instance per region) and called where the region is defined.

For each variation the observables are re-filled into their own folder "<region>_SC_<variation>" with
    event.weight_variation = event.weight * (varied weight / nominal weight)
i.e. the nominal correction already applied to event.weight is swapped for its varied value. Variations:
  nominal                                   MC and data
  pu, prefiring, btag_<8 sources>,
  ttag_corr, ttag_uncorr, tmistag,
  nnloQCD, ewk                              ratio to the applied nominal weight        (_up / _down)
  mu_<id/iso/trigger>_<stat|syst>, mu_reco  muon channel only                          (_up / _down)
  ele_<id/reco/trigger>                     electron channel only                      (_up / _down)
  isr, fsr                                  PS weights (stored relative to nominal)    (_up / _down)
  murmuf_<6 points>                         ME scale weights (stored relative to nominal)
  pdf0 ... pdf99                            NNPDF3.1 replica / nominal LHE weight
  NoSC_up, NoSC_down                        signal-shape variation, TTToSemiLeptonic only, mapped observables only
The input weights are the branches written by the scale-factor and weight modules (all stored in the AnalysisTree).
event.weight is unchanged on return.
*/
class ZprimeSemiLeptonicSystematicsModule {
public:
  // region      : folder prefix, e.g. "Merged" -> "Merged_SC_nominal", "Merged_SC_pu_up", ...
  // observables : SpinCorrelations branch names to histogram
  // nosc_maps   : observable -> strength-map observable used for its NoSC variation
  //               (e.g. "cHel_P3n_Mtt800_Inf_cosThetaLT0p4" -> "cHel_P3n"); observables without an entry get no NoSC histogram
  // fill_nosc   : true only for the reweighted signal sample (EntanglementNoSCReweighting::is_active())
  ZprimeSemiLeptonicSystematicsModule(uhh2::Context& ctx, const std::string& region,
                                      const std::vector<std::string>& observables,
                                      const std::map<std::string, std::string>& nosc_maps = {},
                                      bool fill_nosc = false);

  void fill(uhh2::Event& event);

private:
  using Folder = std::unique_ptr<ZprimeSemiLeptonicSpinCorrSystHists>;
  Folder book(uhh2::Context& ctx, const std::string& variation, const std::vector<std::string>& observables);

  // variation given as varied / nominal weight
  struct RatioSyst {
    std::string name;
    uhh2::Event::Handle<float> h_nominal, h_up, h_down;
    Folder up, down;
  };
  // variation weight stored relative to nominal (PS and ME-scale weights)
  struct RelativeSyst {
    std::string name;
    uhh2::Event::Handle<float> h_weight;
    Folder folder;
  };

  bool is_mc;
  std::string region_;
  int pdf_start_index; // index of the first NNPDF3.1 replica in genInfo->systweights()

  Folder nominal;
  std::vector<RatioSyst> ratio_systs;
  Folder prefiring_up, prefiring_down;
  std::vector<RelativeSyst> relative_systs;
  std::vector<Folder> pdf;

  bool fill_nosc_;
  Folder nosc_up, nosc_down;
  std::vector<uhh2::Event::Handle<float>> h_nosc_up, h_nosc_down; // one per NoSC observable, same order as the folder
};
