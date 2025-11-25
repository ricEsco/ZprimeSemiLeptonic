#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h"

#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <TLorentzVector.h>
#include <string>
#include <map>
#include <memory>
#include <vector>

class ZprimeSemiLeptonicSystematicsHists: public uhh2::Hists {
public:
  explicit ZprimeSemiLeptonicSystematicsHists(uhh2::Context&, const std::string&);
  virtual void fill(const uhh2::Event&) override;

protected:
  void init();
  bool is_mc, is_Muon;
  bool is_tt, ishotvr, isdeepAK8;
  bool isMuon;
  bool isElectron;
  bool debug;

  uhh2::Event::Handle<float> h_ele_reco;
  uhh2::Event::Handle<float> h_ele_reco_up;
  uhh2::Event::Handle<float> h_ele_reco_down;
  uhh2::Event::Handle<float> h_ele_id;
  uhh2::Event::Handle<float> h_ele_id_up;
  uhh2::Event::Handle<float> h_ele_id_down;
  uhh2::Event::Handle<float> h_ele_trigger;
  uhh2::Event::Handle<float> h_ele_trigger_up;
  uhh2::Event::Handle<float> h_ele_trigger_down;
  uhh2::Event::Handle<float> h_mu_reco;
  uhh2::Event::Handle<float> h_mu_reco_up;
  uhh2::Event::Handle<float> h_mu_reco_down;
  // uhh2::Event::Handle<float> h_mu_iso;
  // uhh2::Event::Handle<float> h_mu_iso_up;
  // uhh2::Event::Handle<float> h_mu_iso_down;
  uhh2::Event::Handle<float> h_mu_iso_stat;
  uhh2::Event::Handle<float> h_mu_iso_stat_up;
  uhh2::Event::Handle<float> h_mu_iso_stat_down;
  uhh2::Event::Handle<float> h_mu_iso_syst;
  uhh2::Event::Handle<float> h_mu_iso_syst_up;
  uhh2::Event::Handle<float> h_mu_iso_syst_down;

  // uhh2::Event::Handle<float> h_mu_id;
  // uhh2::Event::Handle<float> h_mu_id_up;
  // uhh2::Event::Handle<float> h_mu_id_down;
  uhh2::Event::Handle<float> h_mu_id_stat;
  uhh2::Event::Handle<float> h_mu_id_stat_up;
  uhh2::Event::Handle<float> h_mu_id_stat_down;
  uhh2::Event::Handle<float> h_mu_id_syst;
  uhh2::Event::Handle<float> h_mu_id_syst_up;
  uhh2::Event::Handle<float> h_mu_id_syst_down;

  // uhh2::Event::Handle<float> h_mu_trigger;
  // uhh2::Event::Handle<float> h_mu_trigger_up;
  // uhh2::Event::Handle<float> h_mu_trigger_down;
  uhh2::Event::Handle<float> h_mu_trigger_stat;
  uhh2::Event::Handle<float> h_mu_trigger_stat_up;
  uhh2::Event::Handle<float> h_mu_trigger_stat_down;
  uhh2::Event::Handle<float> h_mu_trigger_syst;
  uhh2::Event::Handle<float> h_mu_trigger_syst_up;
  uhh2::Event::Handle<float> h_mu_trigger_syst_down;

  uhh2::Event::Handle<float> h_pu;
  uhh2::Event::Handle<float> h_pu_up;
  uhh2::Event::Handle<float> h_pu_down;
  uhh2::Event::Handle<float> h_prefiring;
  uhh2::Event::Handle<float> h_prefiring_up;
  uhh2::Event::Handle<float> h_prefiring_down;
  uhh2::Event::Handle<float> h_murmuf_upup;
  uhh2::Event::Handle<float> h_murmuf_upnone;
  uhh2::Event::Handle<float> h_murmuf_noneup;
  uhh2::Event::Handle<float> h_murmuf_nonedown;
  uhh2::Event::Handle<float> h_murmuf_downnone;
  uhh2::Event::Handle<float> h_murmuf_downdown;
  uhh2::Event::Handle<float> h_isr_up;
  uhh2::Event::Handle<float> h_isr_down;
  uhh2::Event::Handle<float> h_fsr_up;
  uhh2::Event::Handle<float> h_fsr_down;
  uhh2::Event::Handle<float> h_btag;
  uhh2::Event::Handle<float> h_btag_cferr1_up;
  uhh2::Event::Handle<float> h_btag_cferr1_down;
  uhh2::Event::Handle<float> h_btag_cferr2_up;
  uhh2::Event::Handle<float> h_btag_cferr2_down;
  uhh2::Event::Handle<float> h_btag_hf_up;
  uhh2::Event::Handle<float> h_btag_hf_down;
  uhh2::Event::Handle<float> h_btag_hfstats1_up;
  uhh2::Event::Handle<float> h_btag_hfstats1_down;
  uhh2::Event::Handle<float> h_btag_hfstats2_up;
  uhh2::Event::Handle<float> h_btag_hfstats2_down;
  uhh2::Event::Handle<float> h_btag_lf_up;
  uhh2::Event::Handle<float> h_btag_lf_down;
  uhh2::Event::Handle<float> h_btag_lfstats1_up;
  uhh2::Event::Handle<float> h_btag_lfstats1_down;
  uhh2::Event::Handle<float> h_btag_lfstats2_up;
  uhh2::Event::Handle<float> h_btag_lfstats2_down;
  uhh2::Event::Handle<float> h_ttag;
  uhh2::Event::Handle<float> h_ttag_corr_up;
  uhh2::Event::Handle<float> h_ttag_corr_down;
  uhh2::Event::Handle<float> h_ttag_uncorr_up;
  uhh2::Event::Handle<float> h_ttag_uncorr_down;
  uhh2::Event::Handle<float> h_tmistag;
  uhh2::Event::Handle<float> h_tmistag_up;
  uhh2::Event::Handle<float> h_tmistag_down;
  uhh2::Event::Handle<float> h_toppt_a_up;
  uhh2::Event::Handle<float> h_toppt_a_down;
  uhh2::Event::Handle<float> h_toppt_b_up;
  uhh2::Event::Handle<float> h_toppt_b_down;

  

  TH1F *DeltaY;
  TH1F *DeltaY_mu_reco_up;
  TH1F *DeltaY_mu_reco_down;
  TH1F *DeltaY_pu_up;
  TH1F *DeltaY_pu_down;
  TH1F *DeltaY_prefiring_up;
  TH1F *DeltaY_prefiring_down;
  TH1F *DeltaY_mu_iso_stat_up;
  TH1F *DeltaY_mu_iso_stat_down;
  TH1F *DeltaY_mu_iso_syst_up;
  TH1F *DeltaY_mu_iso_syst_down;
  TH1F *DeltaY_mu_id_stat_up;
  TH1F *DeltaY_mu_id_stat_down;
  TH1F *DeltaY_mu_id_syst_up;
  TH1F *DeltaY_mu_id_syst_down;
  TH1F *DeltaY_mu_trigger_stat_up;
  TH1F *DeltaY_mu_trigger_stat_down;
  TH1F *DeltaY_mu_trigger_syst_up;
  TH1F *DeltaY_mu_trigger_syst_down;
  TH1F *DeltaY_ele_id_up;
  TH1F *DeltaY_ele_id_down;
  TH1F *DeltaY_ele_trigger_up;
  TH1F *DeltaY_ele_trigger_down;
  TH1F *DeltaY_ele_reco_up;
  TH1F *DeltaY_ele_reco_down;
  TH1F *DeltaY_murmuf_upup;
  TH1F *DeltaY_murmuf_upnone;
  TH1F *DeltaY_murmuf_noneup;
  TH1F *DeltaY_murmuf_nonedown;
  TH1F *DeltaY_murmuf_downnone;
  TH1F *DeltaY_murmuf_downdown;
  TH1F *DeltaY_isr_up;
  TH1F *DeltaY_isr_down;
  TH1F *DeltaY_fsr_up;
  TH1F *DeltaY_fsr_down;
  TH1F *DeltaY_btag_cferr1_up;
  TH1F *DeltaY_btag_cferr1_down;
  TH1F *DeltaY_btag_cferr2_up;
  TH1F *DeltaY_btag_cferr2_down;
  TH1F *DeltaY_btag_hf_up;
  TH1F *DeltaY_btag_hf_down;
  TH1F *DeltaY_btag_hfstats1_up;
  TH1F *DeltaY_btag_hfstats1_down;
  TH1F *DeltaY_btag_hfstats2_up;
  TH1F *DeltaY_btag_hfstats2_down;
  TH1F *DeltaY_btag_lf_up;
  TH1F *DeltaY_btag_lf_down;
  TH1F *DeltaY_btag_lfstats1_up;
  TH1F *DeltaY_btag_lfstats1_down;
  TH1F *DeltaY_btag_lfstats2_up;
  TH1F *DeltaY_btag_lfstats2_down;
  TH1F *DeltaY_ttag_corr_up;
  TH1F *DeltaY_ttag_corr_down;
  TH1F *DeltaY_ttag_uncorr_up;
  TH1F *DeltaY_ttag_uncorr_down;
  TH1F *DeltaY_tmistag_up;
  TH1F *DeltaY_tmistag_down;
  TH1F *DeltaY_toppt_a_up;
  TH1F *DeltaY_toppt_a_down;
  TH1F *DeltaY_toppt_b_up;
  TH1F *DeltaY_toppt_b_down;
  // --- Template method: xi = tanh(DeltaY), multiple binning schemes ---
  TH1F *DeltaY_xi_reco_6;
  TH1F *DeltaY_xi_reco_12;
  TH1F *DeltaY_xi_reco_18;
  TH1F *DeltaY_xi_reco_24;
  TH1F *DeltaY_xi_reco_30;
  TH1F *DeltaY_xi_reco_36;
  TH1F *DeltaY_xi_reco_50;
  // NoAC systematics (f=+1 is Down, f=-1 is Up) - all binning schemes
  TH1F *DeltaY_xi_reco_6_NoAC_up;
  TH1F *DeltaY_xi_reco_6_NoAC_down;
  TH1F *DeltaY_xi_reco_12_NoAC_up;
  TH1F *DeltaY_xi_reco_12_NoAC_down;
  TH1F *DeltaY_xi_reco_18_NoAC_up;
  TH1F *DeltaY_xi_reco_18_NoAC_down;
  TH1F *DeltaY_xi_reco_24_NoAC_up;
  TH1F *DeltaY_xi_reco_24_NoAC_down;
  TH1F *DeltaY_xi_reco_30_NoAC_up;
  TH1F *DeltaY_xi_reco_30_NoAC_down;
  TH1F *DeltaY_xi_reco_36_NoAC_up;
  TH1F *DeltaY_xi_reco_36_NoAC_down;
  TH1F *DeltaY_xi_reco_50_NoAC_up;
  TH1F *DeltaY_xi_reco_50_NoAC_down;
  // NoAC systematics for f=±2, ±8, ±12 - all binning schemes
  TH1F *DeltaY_xi_reco_6_NoAC_up_f2;   TH1F *DeltaY_xi_reco_6_NoAC_down_f2;
  TH1F *DeltaY_xi_reco_6_NoAC_up_f8;   TH1F *DeltaY_xi_reco_6_NoAC_down_f8;
  TH1F *DeltaY_xi_reco_6_NoAC_up_f12;  TH1F *DeltaY_xi_reco_6_NoAC_down_f12;
  TH1F *DeltaY_xi_reco_12_NoAC_up_f2;  TH1F *DeltaY_xi_reco_12_NoAC_down_f2;
  TH1F *DeltaY_xi_reco_12_NoAC_up_f8;  TH1F *DeltaY_xi_reco_12_NoAC_down_f8;
  TH1F *DeltaY_xi_reco_12_NoAC_up_f12; TH1F *DeltaY_xi_reco_12_NoAC_down_f12;
  TH1F *DeltaY_xi_reco_18_NoAC_up_f2;  TH1F *DeltaY_xi_reco_18_NoAC_down_f2;
  TH1F *DeltaY_xi_reco_18_NoAC_up_f8;  TH1F *DeltaY_xi_reco_18_NoAC_down_f8;
  TH1F *DeltaY_xi_reco_18_NoAC_up_f12; TH1F *DeltaY_xi_reco_18_NoAC_down_f12;
  TH1F *DeltaY_xi_reco_24_NoAC_up_f2;  TH1F *DeltaY_xi_reco_24_NoAC_down_f2;
  TH1F *DeltaY_xi_reco_24_NoAC_up_f8;  TH1F *DeltaY_xi_reco_24_NoAC_down_f8;
  TH1F *DeltaY_xi_reco_24_NoAC_up_f12; TH1F *DeltaY_xi_reco_24_NoAC_down_f12;
  TH1F *DeltaY_xi_reco_30_NoAC_up_f2;  TH1F *DeltaY_xi_reco_30_NoAC_down_f2;
  TH1F *DeltaY_xi_reco_30_NoAC_up_f8;  TH1F *DeltaY_xi_reco_30_NoAC_down_f8;
  TH1F *DeltaY_xi_reco_30_NoAC_up_f12; TH1F *DeltaY_xi_reco_30_NoAC_down_f12;
  TH1F *DeltaY_xi_reco_36_NoAC_up_f2;  TH1F *DeltaY_xi_reco_36_NoAC_down_f2;
  TH1F *DeltaY_xi_reco_36_NoAC_up_f8;  TH1F *DeltaY_xi_reco_36_NoAC_down_f8;
  TH1F *DeltaY_xi_reco_36_NoAC_up_f12; TH1F *DeltaY_xi_reco_36_NoAC_down_f12;
  TH1F *DeltaY_xi_reco_50_NoAC_up_f2;  TH1F *DeltaY_xi_reco_50_NoAC_down_f2;
  TH1F *DeltaY_xi_reco_50_NoAC_up_f8;  TH1F *DeltaY_xi_reco_50_NoAC_down_f8;
  TH1F *DeltaY_xi_reco_50_NoAC_up_f12; TH1F *DeltaY_xi_reco_50_NoAC_down_f12;
  // lepton/trigger/pileup/prefiring - 6 bins
  TH1F *DeltaY_xi_reco_6_ele_reco_up;      TH1F *DeltaY_xi_reco_6_ele_reco_down;
  TH1F *DeltaY_xi_reco_6_ele_id_up;        TH1F *DeltaY_xi_reco_6_ele_id_down;
  TH1F *DeltaY_xi_reco_6_ele_trigger_up;   TH1F *DeltaY_xi_reco_6_ele_trigger_down;
  TH1F *DeltaY_xi_reco_6_mu_reco_up;       TH1F *DeltaY_xi_reco_6_mu_reco_down;
  TH1F *DeltaY_xi_reco_6_mu_iso_stat_up;   TH1F *DeltaY_xi_reco_6_mu_iso_stat_down;
  TH1F *DeltaY_xi_reco_6_mu_iso_syst_up;   TH1F *DeltaY_xi_reco_6_mu_iso_syst_down;
  TH1F *DeltaY_xi_reco_6_mu_id_stat_up;    TH1F *DeltaY_xi_reco_6_mu_id_stat_down;
  TH1F *DeltaY_xi_reco_6_mu_id_syst_up;    TH1F *DeltaY_xi_reco_6_mu_id_syst_down;
  TH1F *DeltaY_xi_reco_6_mu_trigger_stat_up; TH1F *DeltaY_xi_reco_6_mu_trigger_stat_down;
  TH1F *DeltaY_xi_reco_6_mu_trigger_syst_up; TH1F *DeltaY_xi_reco_6_mu_trigger_syst_down;
  TH1F *DeltaY_xi_reco_6_pu_up;            TH1F *DeltaY_xi_reco_6_pu_down;
  TH1F *DeltaY_xi_reco_6_prefiring_up;     TH1F *DeltaY_xi_reco_6_prefiring_down;
  // scales, isr fsr - 6 bins
  TH1F *DeltaY_xi_reco_6_murmuf_upup;      TH1F *DeltaY_xi_reco_6_murmuf_upnone;
  TH1F *DeltaY_xi_reco_6_murmuf_noneup;    TH1F *DeltaY_xi_reco_6_murmuf_nonedown;
  TH1F *DeltaY_xi_reco_6_murmuf_downnone;  TH1F *DeltaY_xi_reco_6_murmuf_downdown;
  TH1F *DeltaY_xi_reco_6_isr_up;           TH1F *DeltaY_xi_reco_6_isr_down;
  TH1F *DeltaY_xi_reco_6_fsr_up;           TH1F *DeltaY_xi_reco_6_fsr_down;
  // btag - 6 bins
  TH1F *DeltaY_xi_reco_6_btag_cferr1_up;   TH1F *DeltaY_xi_reco_6_btag_cferr1_down;
  TH1F *DeltaY_xi_reco_6_btag_cferr2_up;   TH1F *DeltaY_xi_reco_6_btag_cferr2_down;
  TH1F *DeltaY_xi_reco_6_btag_hf_up;       TH1F *DeltaY_xi_reco_6_btag_hf_down;
  TH1F *DeltaY_xi_reco_6_btag_hfstats1_up; TH1F *DeltaY_xi_reco_6_btag_hfstats1_down;
  TH1F *DeltaY_xi_reco_6_btag_hfstats2_up; TH1F *DeltaY_xi_reco_6_btag_hfstats2_down;
  TH1F *DeltaY_xi_reco_6_btag_lf_up;       TH1F *DeltaY_xi_reco_6_btag_lf_down;
  TH1F *DeltaY_xi_reco_6_btag_lfstats1_up; TH1F *DeltaY_xi_reco_6_btag_lfstats1_down;
  TH1F *DeltaY_xi_reco_6_btag_lfstats2_up; TH1F *DeltaY_xi_reco_6_btag_lfstats2_down;
  // ttag, mistag, toppt - 6 bins
  TH1F *DeltaY_xi_reco_6_ttag_corr_up;     TH1F *DeltaY_xi_reco_6_ttag_corr_down;
  TH1F *DeltaY_xi_reco_6_ttag_uncorr_up;   TH1F *DeltaY_xi_reco_6_ttag_uncorr_down;
  TH1F *DeltaY_xi_reco_6_tmistag_up;       TH1F *DeltaY_xi_reco_6_tmistag_down;
  TH1F *DeltaY_xi_reco_6_toppt_a_up;       TH1F *DeltaY_xi_reco_6_toppt_a_down;
  TH1F *DeltaY_xi_reco_6_toppt_b_up;       TH1F *DeltaY_xi_reco_6_toppt_b_down;

  // lepton/trigger/pileup/prefiring - 12 bins
  TH1F *DeltaY_xi_reco_12_ele_reco_up;      TH1F *DeltaY_xi_reco_12_ele_reco_down;
  TH1F *DeltaY_xi_reco_12_ele_id_up;        TH1F *DeltaY_xi_reco_12_ele_id_down;
  TH1F *DeltaY_xi_reco_12_ele_trigger_up;   TH1F *DeltaY_xi_reco_12_ele_trigger_down;
  TH1F *DeltaY_xi_reco_12_mu_reco_up;       TH1F *DeltaY_xi_reco_12_mu_reco_down;
  TH1F *DeltaY_xi_reco_12_mu_iso_stat_up;   TH1F *DeltaY_xi_reco_12_mu_iso_stat_down;
  TH1F *DeltaY_xi_reco_12_mu_iso_syst_up;   TH1F *DeltaY_xi_reco_12_mu_iso_syst_down;
  TH1F *DeltaY_xi_reco_12_mu_id_stat_up;    TH1F *DeltaY_xi_reco_12_mu_id_stat_down;
  TH1F *DeltaY_xi_reco_12_mu_id_syst_up;    TH1F *DeltaY_xi_reco_12_mu_id_syst_down;
  TH1F *DeltaY_xi_reco_12_mu_trigger_stat_up; TH1F *DeltaY_xi_reco_12_mu_trigger_stat_down;
  TH1F *DeltaY_xi_reco_12_mu_trigger_syst_up; TH1F *DeltaY_xi_reco_12_mu_trigger_syst_down;
  TH1F *DeltaY_xi_reco_12_pu_up;            TH1F *DeltaY_xi_reco_12_pu_down;
  TH1F *DeltaY_xi_reco_12_prefiring_up;     TH1F *DeltaY_xi_reco_12_prefiring_down;
  // scales, isr fsr - 12 bins
  TH1F *DeltaY_xi_reco_12_murmuf_upup;      TH1F *DeltaY_xi_reco_12_murmuf_upnone;
  TH1F *DeltaY_xi_reco_12_murmuf_noneup;    TH1F *DeltaY_xi_reco_12_murmuf_nonedown;
  TH1F *DeltaY_xi_reco_12_murmuf_downnone;  TH1F *DeltaY_xi_reco_12_murmuf_downdown;
  TH1F *DeltaY_xi_reco_12_isr_up;           TH1F *DeltaY_xi_reco_12_isr_down;
  TH1F *DeltaY_xi_reco_12_fsr_up;           TH1F *DeltaY_xi_reco_12_fsr_down;
  // btag - 12 bins
  TH1F *DeltaY_xi_reco_12_btag_cferr1_up;   TH1F *DeltaY_xi_reco_12_btag_cferr1_down;
  TH1F *DeltaY_xi_reco_12_btag_cferr2_up;   TH1F *DeltaY_xi_reco_12_btag_cferr2_down;
  TH1F *DeltaY_xi_reco_12_btag_hf_up;       TH1F *DeltaY_xi_reco_12_btag_hf_down;
  TH1F *DeltaY_xi_reco_12_btag_hfstats1_up; TH1F *DeltaY_xi_reco_12_btag_hfstats1_down;
  TH1F *DeltaY_xi_reco_12_btag_hfstats2_up; TH1F *DeltaY_xi_reco_12_btag_hfstats2_down;
  TH1F *DeltaY_xi_reco_12_btag_lf_up;       TH1F *DeltaY_xi_reco_12_btag_lf_down;
  TH1F *DeltaY_xi_reco_12_btag_lfstats1_up; TH1F *DeltaY_xi_reco_12_btag_lfstats1_down;
  TH1F *DeltaY_xi_reco_12_btag_lfstats2_up; TH1F *DeltaY_xi_reco_12_btag_lfstats2_down;
  // ttag, mistag, toppt - 12 bins
  TH1F *DeltaY_xi_reco_12_ttag_corr_up;     TH1F *DeltaY_xi_reco_12_ttag_corr_down;
  TH1F *DeltaY_xi_reco_12_ttag_uncorr_up;   TH1F *DeltaY_xi_reco_12_ttag_uncorr_down;
  TH1F *DeltaY_xi_reco_12_tmistag_up;       TH1F *DeltaY_xi_reco_12_tmistag_down;
  TH1F *DeltaY_xi_reco_12_toppt_a_up;       TH1F *DeltaY_xi_reco_12_toppt_a_down;
  TH1F *DeltaY_xi_reco_12_toppt_b_up;       TH1F *DeltaY_xi_reco_12_toppt_b_down;

  // lepton/trigger/pileup/prefiring - 50 bins
  TH1F *DeltaY_xi_reco_50_ele_reco_up;      TH1F *DeltaY_xi_reco_50_ele_reco_down;
  TH1F *DeltaY_xi_reco_50_ele_id_up;        TH1F *DeltaY_xi_reco_50_ele_id_down;
  TH1F *DeltaY_xi_reco_50_ele_trigger_up;   TH1F *DeltaY_xi_reco_50_ele_trigger_down;
  TH1F *DeltaY_xi_reco_50_mu_reco_up;       TH1F *DeltaY_xi_reco_50_mu_reco_down;
  TH1F *DeltaY_xi_reco_50_mu_iso_stat_up;   TH1F *DeltaY_xi_reco_50_mu_iso_stat_down;
  TH1F *DeltaY_xi_reco_50_mu_iso_syst_up;   TH1F *DeltaY_xi_reco_50_mu_iso_syst_down;
  TH1F *DeltaY_xi_reco_50_mu_id_stat_up;    TH1F *DeltaY_xi_reco_50_mu_id_stat_down;
  TH1F *DeltaY_xi_reco_50_mu_id_syst_up;    TH1F *DeltaY_xi_reco_50_mu_id_syst_down;
  TH1F *DeltaY_xi_reco_50_mu_trigger_stat_up; TH1F *DeltaY_xi_reco_50_mu_trigger_stat_down;
  TH1F *DeltaY_xi_reco_50_mu_trigger_syst_up; TH1F *DeltaY_xi_reco_50_mu_trigger_syst_down;
  TH1F *DeltaY_xi_reco_50_pu_up;            TH1F *DeltaY_xi_reco_50_pu_down;
  TH1F *DeltaY_xi_reco_50_prefiring_up;     TH1F *DeltaY_xi_reco_50_prefiring_down;
  // scales, isr fsr - 50 bins
  TH1F *DeltaY_xi_reco_50_murmuf_upup;      TH1F *DeltaY_xi_reco_50_murmuf_upnone;
  TH1F *DeltaY_xi_reco_50_murmuf_noneup;    TH1F *DeltaY_xi_reco_50_murmuf_nonedown;
  TH1F *DeltaY_xi_reco_50_murmuf_downnone;  TH1F *DeltaY_xi_reco_50_murmuf_downdown;
  TH1F *DeltaY_xi_reco_50_isr_up;           TH1F *DeltaY_xi_reco_50_isr_down;
  TH1F *DeltaY_xi_reco_50_fsr_up;           TH1F *DeltaY_xi_reco_50_fsr_down;
  // btag - 50 bins
  TH1F *DeltaY_xi_reco_50_btag_cferr1_up;   TH1F *DeltaY_xi_reco_50_btag_cferr1_down;
  TH1F *DeltaY_xi_reco_50_btag_cferr2_up;   TH1F *DeltaY_xi_reco_50_btag_cferr2_down;
  TH1F *DeltaY_xi_reco_50_btag_hf_up;       TH1F *DeltaY_xi_reco_50_btag_hf_down;
  TH1F *DeltaY_xi_reco_50_btag_hfstats1_up; TH1F *DeltaY_xi_reco_50_btag_hfstats1_down;
  TH1F *DeltaY_xi_reco_50_btag_hfstats2_up; TH1F *DeltaY_xi_reco_50_btag_hfstats2_down;
  TH1F *DeltaY_xi_reco_50_btag_lf_up;       TH1F *DeltaY_xi_reco_50_btag_lf_down;
  TH1F *DeltaY_xi_reco_50_btag_lfstats1_up; TH1F *DeltaY_xi_reco_50_btag_lfstats1_down;
  TH1F *DeltaY_xi_reco_50_btag_lfstats2_up; TH1F *DeltaY_xi_reco_50_btag_lfstats2_down;
  // ttag, mistag, toppt - 50 bins
  TH1F *DeltaY_xi_reco_50_ttag_corr_up;     TH1F *DeltaY_xi_reco_50_ttag_corr_down;
  TH1F *DeltaY_xi_reco_50_ttag_uncorr_up;   TH1F *DeltaY_xi_reco_50_ttag_uncorr_down;
  TH1F *DeltaY_xi_reco_50_tmistag_up;       TH1F *DeltaY_xi_reco_50_tmistag_down;
  TH1F *DeltaY_xi_reco_50_toppt_a_up;       TH1F *DeltaY_xi_reco_50_toppt_a_down;
  TH1F *DeltaY_xi_reco_50_toppt_b_up;       TH1F *DeltaY_xi_reco_50_toppt_b_down;

  // lepton/trigger/pileup/prefiring - 24 bins
  TH1F *DeltaY_xi_reco_24_ele_reco_up;      TH1F *DeltaY_xi_reco_24_ele_reco_down;
  TH1F *DeltaY_xi_reco_24_ele_id_up;        TH1F *DeltaY_xi_reco_24_ele_id_down;
  TH1F *DeltaY_xi_reco_24_ele_trigger_up;   TH1F *DeltaY_xi_reco_24_ele_trigger_down;
  TH1F *DeltaY_xi_reco_24_mu_reco_up;       TH1F *DeltaY_xi_reco_24_mu_reco_down;
  TH1F *DeltaY_xi_reco_24_mu_iso_stat_up;   TH1F *DeltaY_xi_reco_24_mu_iso_stat_down;
  TH1F *DeltaY_xi_reco_24_mu_iso_syst_up;   TH1F *DeltaY_xi_reco_24_mu_iso_syst_down;
  TH1F *DeltaY_xi_reco_24_mu_id_stat_up;    TH1F *DeltaY_xi_reco_24_mu_id_stat_down;
  TH1F *DeltaY_xi_reco_24_mu_id_syst_up;    TH1F *DeltaY_xi_reco_24_mu_id_syst_down;
  TH1F *DeltaY_xi_reco_24_mu_trigger_stat_up; TH1F *DeltaY_xi_reco_24_mu_trigger_stat_down;
  TH1F *DeltaY_xi_reco_24_mu_trigger_syst_up; TH1F *DeltaY_xi_reco_24_mu_trigger_syst_down;
  TH1F *DeltaY_xi_reco_24_pu_up;            TH1F *DeltaY_xi_reco_24_pu_down;
  TH1F *DeltaY_xi_reco_24_prefiring_up;     TH1F *DeltaY_xi_reco_24_prefiring_down;
  // scales, isr fsr - 24 bins
  TH1F *DeltaY_xi_reco_24_murmuf_upup;      TH1F *DeltaY_xi_reco_24_murmuf_upnone;
  TH1F *DeltaY_xi_reco_24_murmuf_noneup;    TH1F *DeltaY_xi_reco_24_murmuf_nonedown;
  TH1F *DeltaY_xi_reco_24_murmuf_downnone;  TH1F *DeltaY_xi_reco_24_murmuf_downdown;
  TH1F *DeltaY_xi_reco_24_isr_up;           TH1F *DeltaY_xi_reco_24_isr_down;
  TH1F *DeltaY_xi_reco_24_fsr_up;           TH1F *DeltaY_xi_reco_24_fsr_down;
  // btag - 24 bins
  TH1F *DeltaY_xi_reco_24_btag_cferr1_up;   TH1F *DeltaY_xi_reco_24_btag_cferr1_down;
  TH1F *DeltaY_xi_reco_24_btag_cferr2_up;   TH1F *DeltaY_xi_reco_24_btag_cferr2_down;
  TH1F *DeltaY_xi_reco_24_btag_hf_up;       TH1F *DeltaY_xi_reco_24_btag_hf_down;
  TH1F *DeltaY_xi_reco_24_btag_hfstats1_up; TH1F *DeltaY_xi_reco_24_btag_hfstats1_down;
  TH1F *DeltaY_xi_reco_24_btag_hfstats2_up; TH1F *DeltaY_xi_reco_24_btag_hfstats2_down;
  TH1F *DeltaY_xi_reco_24_btag_lf_up;       TH1F *DeltaY_xi_reco_24_btag_lf_down;
  TH1F *DeltaY_xi_reco_24_btag_lfstats1_up; TH1F *DeltaY_xi_reco_24_btag_lfstats1_down;
  TH1F *DeltaY_xi_reco_24_btag_lfstats2_up; TH1F *DeltaY_xi_reco_24_btag_lfstats2_down;
  // ttag, mistag, toppt - 24 bins
  TH1F *DeltaY_xi_reco_24_ttag_corr_up;     TH1F *DeltaY_xi_reco_24_ttag_corr_down;
  TH1F *DeltaY_xi_reco_24_ttag_uncorr_up;   TH1F *DeltaY_xi_reco_24_ttag_uncorr_down;
  TH1F *DeltaY_xi_reco_24_tmistag_up;       TH1F *DeltaY_xi_reco_24_tmistag_down;
  TH1F *DeltaY_xi_reco_24_toppt_a_up;       TH1F *DeltaY_xi_reco_24_toppt_a_down;
  TH1F *DeltaY_xi_reco_24_toppt_b_up;       TH1F *DeltaY_xi_reco_24_toppt_b_down;

  // lepton/trigger/pileup/prefiring - 30 bins
  TH1F *DeltaY_xi_reco_30_ele_reco_up;      TH1F *DeltaY_xi_reco_30_ele_reco_down;
  TH1F *DeltaY_xi_reco_30_ele_id_up;        TH1F *DeltaY_xi_reco_30_ele_id_down;
  TH1F *DeltaY_xi_reco_30_ele_trigger_up;   TH1F *DeltaY_xi_reco_30_ele_trigger_down;
  TH1F *DeltaY_xi_reco_30_mu_reco_up;       TH1F *DeltaY_xi_reco_30_mu_reco_down;
  TH1F *DeltaY_xi_reco_30_mu_iso_stat_up;   TH1F *DeltaY_xi_reco_30_mu_iso_stat_down;
  TH1F *DeltaY_xi_reco_30_mu_iso_syst_up;   TH1F *DeltaY_xi_reco_30_mu_iso_syst_down;
  TH1F *DeltaY_xi_reco_30_mu_id_stat_up;    TH1F *DeltaY_xi_reco_30_mu_id_stat_down;
  TH1F *DeltaY_xi_reco_30_mu_id_syst_up;    TH1F *DeltaY_xi_reco_30_mu_id_syst_down;
  TH1F *DeltaY_xi_reco_30_mu_trigger_stat_up; TH1F *DeltaY_xi_reco_30_mu_trigger_stat_down;
  TH1F *DeltaY_xi_reco_30_mu_trigger_syst_up; TH1F *DeltaY_xi_reco_30_mu_trigger_syst_down;
  TH1F *DeltaY_xi_reco_30_pu_up;            TH1F *DeltaY_xi_reco_30_pu_down;
  TH1F *DeltaY_xi_reco_30_prefiring_up;     TH1F *DeltaY_xi_reco_30_prefiring_down;
  // scales, isr fsr - 30 bins
  TH1F *DeltaY_xi_reco_30_murmuf_upup;      TH1F *DeltaY_xi_reco_30_murmuf_upnone;
  TH1F *DeltaY_xi_reco_30_murmuf_noneup;    TH1F *DeltaY_xi_reco_30_murmuf_nonedown;
  TH1F *DeltaY_xi_reco_30_murmuf_downnone;  TH1F *DeltaY_xi_reco_30_murmuf_downdown;
  TH1F *DeltaY_xi_reco_30_isr_up;           TH1F *DeltaY_xi_reco_30_isr_down;
  TH1F *DeltaY_xi_reco_30_fsr_up;           TH1F *DeltaY_xi_reco_30_fsr_down;
  // btag - 30 bins
  TH1F *DeltaY_xi_reco_30_btag_cferr1_up;   TH1F *DeltaY_xi_reco_30_btag_cferr1_down;
  TH1F *DeltaY_xi_reco_30_btag_cferr2_up;   TH1F *DeltaY_xi_reco_30_btag_cferr2_down;
  TH1F *DeltaY_xi_reco_30_btag_hf_up;       TH1F *DeltaY_xi_reco_30_btag_hf_down;
  TH1F *DeltaY_xi_reco_30_btag_hfstats1_up; TH1F *DeltaY_xi_reco_30_btag_hfstats1_down;
  TH1F *DeltaY_xi_reco_30_btag_hfstats2_up; TH1F *DeltaY_xi_reco_30_btag_hfstats2_down;
  TH1F *DeltaY_xi_reco_30_btag_lf_up;       TH1F *DeltaY_xi_reco_30_btag_lf_down;
  TH1F *DeltaY_xi_reco_30_btag_lfstats1_up; TH1F *DeltaY_xi_reco_30_btag_lfstats1_down;
  TH1F *DeltaY_xi_reco_30_btag_lfstats2_up; TH1F *DeltaY_xi_reco_30_btag_lfstats2_down;
  // ttag, mistag, toppt - 30 bins
  TH1F *DeltaY_xi_reco_30_ttag_corr_up;     TH1F *DeltaY_xi_reco_30_ttag_corr_down;
  TH1F *DeltaY_xi_reco_30_ttag_uncorr_up;   TH1F *DeltaY_xi_reco_30_ttag_uncorr_down;
  TH1F *DeltaY_xi_reco_30_tmistag_up;       TH1F *DeltaY_xi_reco_30_tmistag_down;
  TH1F *DeltaY_xi_reco_30_toppt_a_up;       TH1F *DeltaY_xi_reco_30_toppt_a_down;
  TH1F *DeltaY_xi_reco_30_toppt_b_up;       TH1F *DeltaY_xi_reco_30_toppt_b_down;


  TH1F *DeltaY_reco_d1_mu_reco_up;
  TH1F *DeltaY_reco_d1_mu_reco_down;
  TH1F *DeltaY_reco_d1_pu_up;
  TH1F *DeltaY_reco_d1_pu_down;
  TH1F *DeltaY_reco_d1_prefiring_up;
  TH1F *DeltaY_reco_d1_prefiring_down;
  TH1F *DeltaY_reco_d1_mu_iso_stat_up;
  TH1F *DeltaY_reco_d1_mu_iso_stat_down;
  TH1F *DeltaY_reco_d1_mu_iso_syst_up;
  TH1F *DeltaY_reco_d1_mu_iso_syst_down;
  TH1F *DeltaY_reco_d1_mu_id_stat_up;
  TH1F *DeltaY_reco_d1_mu_id_stat_down;
  TH1F *DeltaY_reco_d1_mu_id_syst_up;
  TH1F *DeltaY_reco_d1_mu_id_syst_down;
  TH1F *DeltaY_reco_d1_mu_trigger_stat_up;
  TH1F *DeltaY_reco_d1_mu_trigger_stat_down;
  TH1F *DeltaY_reco_d1_mu_trigger_syst_up;
  TH1F *DeltaY_reco_d1_mu_trigger_syst_down;
  TH1F *DeltaY_reco_d1_ele_id_up;
  TH1F *DeltaY_reco_d1_ele_id_down;
  TH1F *DeltaY_reco_d1_ele_trigger_up;
  TH1F *DeltaY_reco_d1_ele_trigger_down;
  TH1F *DeltaY_reco_d1_ele_reco_up;
  TH1F *DeltaY_reco_d1_ele_reco_down;
  TH1F *DeltaY_reco_d1_murmuf_upup;
  TH1F *DeltaY_reco_d1_murmuf_upnone;
  TH1F *DeltaY_reco_d1_murmuf_noneup;
  TH1F *DeltaY_reco_d1_murmuf_nonedown;
  TH1F *DeltaY_reco_d1_murmuf_downnone;
  TH1F *DeltaY_reco_d1_murmuf_downdown;
  TH1F *DeltaY_reco_d1_isr_up;
  TH1F *DeltaY_reco_d1_isr_down;
  TH1F *DeltaY_reco_d1_fsr_up;
  TH1F *DeltaY_reco_d1_fsr_down;
  TH1F *DeltaY_reco_d1_btag_cferr1_up;
  TH1F *DeltaY_reco_d1_btag_cferr1_down;
  TH1F *DeltaY_reco_d1_btag_cferr2_up;
  TH1F *DeltaY_reco_d1_btag_cferr2_down;
  TH1F *DeltaY_reco_d1_btag_hf_up;
  TH1F *DeltaY_reco_d1_btag_hf_down;
  TH1F *DeltaY_reco_d1_btag_hfstats1_up;
  TH1F *DeltaY_reco_d1_btag_hfstats1_down;
  TH1F *DeltaY_reco_d1_btag_hfstats2_up;
  TH1F *DeltaY_reco_d1_btag_hfstats2_down;
  TH1F *DeltaY_reco_d1_btag_lf_up;
  TH1F *DeltaY_reco_d1_btag_lf_down;
  TH1F *DeltaY_reco_d1_btag_lfstats1_up;
  TH1F *DeltaY_reco_d1_btag_lfstats1_down;
  TH1F *DeltaY_reco_d1_btag_lfstats2_up;
  TH1F *DeltaY_reco_d1_btag_lfstats2_down;
  TH1F *DeltaY_reco_d1_ttag_corr_up;
  TH1F *DeltaY_reco_d1_ttag_corr_down;
  TH1F *DeltaY_reco_d1_ttag_uncorr_up;
  TH1F *DeltaY_reco_d1_ttag_uncorr_down;
  TH1F *DeltaY_reco_d1_tmistag_up;
  TH1F *DeltaY_reco_d1_tmistag_down;
  TH1F *DeltaY_reco_d1_toppt_a_up;
  TH1F *DeltaY_reco_d1_toppt_a_down;
  TH1F *DeltaY_reco_d1_toppt_b_up;
  TH1F *DeltaY_reco_d1_toppt_b_down;


  TH1F *DeltaY_reco_d2_mu_reco_up;
  TH1F *DeltaY_reco_d2_mu_reco_down;
  TH1F *DeltaY_reco_d2_pu_up;
  TH1F *DeltaY_reco_d2_pu_down;
  TH1F *DeltaY_reco_d2_prefiring_up;
  TH1F *DeltaY_reco_d2_prefiring_down;
  TH1F *DeltaY_reco_d2_mu_iso_stat_up;
  TH1F *DeltaY_reco_d2_mu_iso_stat_down;
  TH1F *DeltaY_reco_d2_mu_iso_syst_up;
  TH1F *DeltaY_reco_d2_mu_iso_syst_down;
  TH1F *DeltaY_reco_d2_mu_id_stat_up;
  TH1F *DeltaY_reco_d2_mu_id_stat_down;
  TH1F *DeltaY_reco_d2_mu_id_syst_up;
  TH1F *DeltaY_reco_d2_mu_id_syst_down;
  TH1F *DeltaY_reco_d2_mu_trigger_stat_up;
  TH1F *DeltaY_reco_d2_mu_trigger_stat_down;
  TH1F *DeltaY_reco_d2_mu_trigger_syst_up;
  TH1F *DeltaY_reco_d2_mu_trigger_syst_down;
  TH1F *DeltaY_reco_d2_ele_id_up;
  TH1F *DeltaY_reco_d2_ele_id_down;
  TH1F *DeltaY_reco_d2_ele_trigger_up;
  TH1F *DeltaY_reco_d2_ele_trigger_down;
  TH1F *DeltaY_reco_d2_ele_reco_up;
  TH1F *DeltaY_reco_d2_ele_reco_down;
  TH1F *DeltaY_reco_d2_murmuf_upup;
  TH1F *DeltaY_reco_d2_murmuf_upnone;
  TH1F *DeltaY_reco_d2_murmuf_noneup;
  TH1F *DeltaY_reco_d2_murmuf_nonedown;
  TH1F *DeltaY_reco_d2_murmuf_downnone;
  TH1F *DeltaY_reco_d2_murmuf_downdown;
  TH1F *DeltaY_reco_d2_isr_up;
  TH1F *DeltaY_reco_d2_isr_down;
  TH1F *DeltaY_reco_d2_fsr_up;
  TH1F *DeltaY_reco_d2_fsr_down;
  TH1F *DeltaY_reco_d2_btag_cferr1_up;
  TH1F *DeltaY_reco_d2_btag_cferr1_down;
  TH1F *DeltaY_reco_d2_btag_cferr2_up;
  TH1F *DeltaY_reco_d2_btag_cferr2_down;
  TH1F *DeltaY_reco_d2_btag_hf_up;
  TH1F *DeltaY_reco_d2_btag_hf_down;
  TH1F *DeltaY_reco_d2_btag_hfstats1_up;
  TH1F *DeltaY_reco_d2_btag_hfstats1_down;
  TH1F *DeltaY_reco_d2_btag_hfstats2_up;
  TH1F *DeltaY_reco_d2_btag_hfstats2_down;
  TH1F *DeltaY_reco_d2_btag_lf_up;
  TH1F *DeltaY_reco_d2_btag_lf_down;
  TH1F *DeltaY_reco_d2_btag_lfstats1_up;
  TH1F *DeltaY_reco_d2_btag_lfstats1_down;
  TH1F *DeltaY_reco_d2_btag_lfstats2_up;
  TH1F *DeltaY_reco_d2_btag_lfstats2_down;
  TH1F *DeltaY_reco_d2_ttag_corr_up;
  TH1F *DeltaY_reco_d2_ttag_corr_down;
  TH1F *DeltaY_reco_d2_ttag_uncorr_up;
  TH1F *DeltaY_reco_d2_ttag_uncorr_down;
  TH1F *DeltaY_reco_d2_tmistag_up;
  TH1F *DeltaY_reco_d2_tmistag_down;
  TH1F *DeltaY_reco_d2_toppt_a_up;
  TH1F *DeltaY_reco_d2_toppt_a_down;
  TH1F *DeltaY_reco_d2_toppt_b_up;
  TH1F *DeltaY_reco_d2_toppt_b_down;


  TH1F *Sigma_phi_1_mu_reco_up;
  TH1F *Sigma_phi_1_mu_reco_down;
  TH1F *Sigma_phi_1_pu_up;
  TH1F *Sigma_phi_1_pu_down;
  TH1F *Sigma_phi_1_prefiring_up;
  TH1F *Sigma_phi_1_prefiring_down;
  TH1F *Sigma_phi_1_mu_iso_stat_up;
  TH1F *Sigma_phi_1_mu_iso_stat_down;
  TH1F *Sigma_phi_1_mu_iso_syst_up;
  TH1F *Sigma_phi_1_mu_iso_syst_down;
  TH1F *Sigma_phi_1_mu_id_stat_up;
  TH1F *Sigma_phi_1_mu_id_stat_down;
  TH1F *Sigma_phi_1_mu_id_syst_up;
  TH1F *Sigma_phi_1_mu_id_syst_down;
  TH1F *Sigma_phi_1_mu_trigger_stat_up;
  TH1F *Sigma_phi_1_mu_trigger_stat_down;
  TH1F *Sigma_phi_1_mu_trigger_syst_up;
  TH1F *Sigma_phi_1_mu_trigger_syst_down;
  TH1F *Sigma_phi_1_ele_id_up;
  TH1F *Sigma_phi_1_ele_id_down;
  TH1F *Sigma_phi_1_ele_trigger_up;
  TH1F *Sigma_phi_1_ele_trigger_down;
  TH1F *Sigma_phi_1_ele_reco_up;
  TH1F *Sigma_phi_1_ele_reco_down;
  TH1F *Sigma_phi_1_murmuf_upup;
  TH1F *Sigma_phi_1_murmuf_upnone;
  TH1F *Sigma_phi_1_murmuf_noneup;
  TH1F *Sigma_phi_1_murmuf_nonedown;
  TH1F *Sigma_phi_1_murmuf_downnone;
  TH1F *Sigma_phi_1_murmuf_downdown;
  TH1F *Sigma_phi_1_isr_up;
  TH1F *Sigma_phi_1_isr_down;
  TH1F *Sigma_phi_1_fsr_up;
  TH1F *Sigma_phi_1_fsr_down;
  TH1F *Sigma_phi_1_btag_cferr1_up;
  TH1F *Sigma_phi_1_btag_cferr1_down;
  TH1F *Sigma_phi_1_btag_cferr2_up;
  TH1F *Sigma_phi_1_btag_cferr2_down;
  TH1F *Sigma_phi_1_btag_hf_up;
  TH1F *Sigma_phi_1_btag_hf_down;
  TH1F *Sigma_phi_1_btag_hfstats1_up;
  TH1F *Sigma_phi_1_btag_hfstats1_down;
  TH1F *Sigma_phi_1_btag_hfstats2_up;
  TH1F *Sigma_phi_1_btag_hfstats2_down;
  TH1F *Sigma_phi_1_btag_lf_up;
  TH1F *Sigma_phi_1_btag_lf_down;
  TH1F *Sigma_phi_1_btag_lfstats1_up;
  TH1F *Sigma_phi_1_btag_lfstats1_down;
  TH1F *Sigma_phi_1_btag_lfstats2_up;
  TH1F *Sigma_phi_1_btag_lfstats2_down;
  TH1F *Sigma_phi_1_ttag_corr_up;
  TH1F *Sigma_phi_1_ttag_corr_down;
  TH1F *Sigma_phi_1_ttag_uncorr_up;
  TH1F *Sigma_phi_1_ttag_uncorr_down;
  TH1F *Sigma_phi_1_tmistag_up;
  TH1F *Sigma_phi_1_tmistag_down;
  TH1F *Sigma_phi_1_toppt_a_up;
  TH1F *Sigma_phi_1_toppt_a_down;
  TH1F *Sigma_phi_1_toppt_b_up;
  TH1F *Sigma_phi_1_toppt_b_down;


  TH1F *Sigma_phi_2_mu_reco_up;
  TH1F *Sigma_phi_2_mu_reco_down;
  TH1F *Sigma_phi_2_pu_up;
  TH1F *Sigma_phi_2_pu_down;
  TH1F *Sigma_phi_2_prefiring_up;
  TH1F *Sigma_phi_2_prefiring_down;
  TH1F *Sigma_phi_2_mu_iso_stat_up;
  TH1F *Sigma_phi_2_mu_iso_stat_down;
  TH1F *Sigma_phi_2_mu_iso_syst_up;
  TH1F *Sigma_phi_2_mu_iso_syst_down;
  TH1F *Sigma_phi_2_mu_id_stat_up;
  TH1F *Sigma_phi_2_mu_id_stat_down;
  TH1F *Sigma_phi_2_mu_id_syst_up;
  TH1F *Sigma_phi_2_mu_id_syst_down;
  TH1F *Sigma_phi_2_mu_trigger_stat_up;
  TH1F *Sigma_phi_2_mu_trigger_stat_down;
  TH1F *Sigma_phi_2_mu_trigger_syst_up;
  TH1F *Sigma_phi_2_mu_trigger_syst_down;
  TH1F *Sigma_phi_2_ele_id_up;
  TH1F *Sigma_phi_2_ele_id_down;
  TH1F *Sigma_phi_2_ele_trigger_up;
  TH1F *Sigma_phi_2_ele_trigger_down;
  TH1F *Sigma_phi_2_ele_reco_up;
  TH1F *Sigma_phi_2_ele_reco_down;
  TH1F *Sigma_phi_2_murmuf_upup;
  TH1F *Sigma_phi_2_murmuf_upnone;
  TH1F *Sigma_phi_2_murmuf_noneup;
  TH1F *Sigma_phi_2_murmuf_nonedown;
  TH1F *Sigma_phi_2_murmuf_downnone;
  TH1F *Sigma_phi_2_murmuf_downdown;
  TH1F *Sigma_phi_2_isr_up;
  TH1F *Sigma_phi_2_isr_down;
  TH1F *Sigma_phi_2_fsr_up;
  TH1F *Sigma_phi_2_fsr_down;
  TH1F *Sigma_phi_2_btag_cferr1_up;
  TH1F *Sigma_phi_2_btag_cferr1_down;
  TH1F *Sigma_phi_2_btag_cferr2_up;
  TH1F *Sigma_phi_2_btag_cferr2_down;
  TH1F *Sigma_phi_2_btag_hf_up;
  TH1F *Sigma_phi_2_btag_hf_down;
  TH1F *Sigma_phi_2_btag_hfstats1_up;
  TH1F *Sigma_phi_2_btag_hfstats1_down;
  TH1F *Sigma_phi_2_btag_hfstats2_up;
  TH1F *Sigma_phi_2_btag_hfstats2_down;
  TH1F *Sigma_phi_2_btag_lf_up;
  TH1F *Sigma_phi_2_btag_lf_down;
  TH1F *Sigma_phi_2_btag_lfstats1_up;
  TH1F *Sigma_phi_2_btag_lfstats1_down;
  TH1F *Sigma_phi_2_btag_lfstats2_up;
  TH1F *Sigma_phi_2_btag_lfstats2_down;
  TH1F *Sigma_phi_2_ttag_corr_up;
  TH1F *Sigma_phi_2_ttag_corr_down;
  TH1F *Sigma_phi_2_ttag_uncorr_up;
  TH1F *Sigma_phi_2_ttag_uncorr_down;
  TH1F *Sigma_phi_2_tmistag_up;
  TH1F *Sigma_phi_2_tmistag_down;
  TH1F *Sigma_phi_2_toppt_a_up;
  TH1F *Sigma_phi_2_toppt_a_down;
  TH1F *Sigma_phi_2_toppt_b_up;
  TH1F *Sigma_phi_2_toppt_b_down;


  TH2F *DeltaY_tt;
  TH2F *DeltaY_mu_reco_up_tt;
  TH2F *DeltaY_mu_reco_down_tt;
  TH2F *DeltaY_pu_up_tt;
  TH2F *DeltaY_pu_down_tt;
  TH2F *DeltaY_prefiring_up_tt;
  TH2F *DeltaY_prefiring_down_tt;
  TH2F *DeltaY_mu_id_stat_up_tt;
  TH2F *DeltaY_mu_id_stat_down_tt;
  TH2F *DeltaY_mu_id_syst_up_tt;
  TH2F *DeltaY_mu_id_syst_down_tt;
  TH2F *DeltaY_mu_iso_stat_up_tt;
  TH2F *DeltaY_mu_iso_stat_down_tt;
  TH2F *DeltaY_mu_iso_syst_up_tt;
  TH2F *DeltaY_mu_iso_syst_down_tt;
  TH2F *DeltaY_mu_trigger_stat_up_tt;
  TH2F *DeltaY_mu_trigger_stat_down_tt;
  TH2F *DeltaY_mu_trigger_syst_up_tt;
  TH2F *DeltaY_mu_trigger_syst_down_tt;
  TH2F *DeltaY_ele_id_up_tt;
  TH2F *DeltaY_ele_id_down_tt;
  TH2F *DeltaY_ele_trigger_up_tt;
  TH2F *DeltaY_ele_trigger_down_tt;
  TH2F *DeltaY_ele_reco_up_tt;
  TH2F *DeltaY_ele_reco_down_tt;
  TH2F *DeltaY_murmuf_upup_tt;
  TH2F *DeltaY_murmuf_upnone_tt;
  TH2F *DeltaY_murmuf_noneup_tt;
  TH2F *DeltaY_murmuf_nonedown_tt;
  TH2F *DeltaY_murmuf_downnone_tt;
  TH2F *DeltaY_murmuf_downdown_tt;
  TH2F *DeltaY_isr_up_tt;
  TH2F *DeltaY_isr_down_tt;
  TH2F *DeltaY_fsr_up_tt;
  TH2F *DeltaY_fsr_down_tt;
  TH2F *DeltaY_btag_cferr1_up_tt;
  TH2F *DeltaY_btag_cferr1_down_tt;
  TH2F *DeltaY_btag_cferr2_up_tt;
  TH2F *DeltaY_btag_cferr2_down_tt;
  TH2F *DeltaY_btag_hf_up_tt;
  TH2F *DeltaY_btag_hf_down_tt;
  TH2F *DeltaY_btag_hfstats1_up_tt;
  TH2F *DeltaY_btag_hfstats1_down_tt;
  TH2F *DeltaY_btag_hfstats2_up_tt;
  TH2F *DeltaY_btag_hfstats2_down_tt;
  TH2F *DeltaY_btag_lf_up_tt;
  TH2F *DeltaY_btag_lf_down_tt;
  TH2F *DeltaY_btag_lfstats1_up_tt;
  TH2F *DeltaY_btag_lfstats1_down_tt;
  TH2F *DeltaY_btag_lfstats2_up_tt;
  TH2F *DeltaY_btag_lfstats2_down_tt;
  TH2F *DeltaY_ttag_corr_up_tt;
  TH2F *DeltaY_ttag_corr_down_tt;
  TH2F *DeltaY_ttag_uncorr_up_tt;
  TH2F *DeltaY_ttag_uncorr_down_tt;
  TH2F *DeltaY_tmistag_up_tt;
  TH2F *DeltaY_tmistag_down_tt;
  TH2F *DeltaY_toppt_a_up_tt;
  TH2F *DeltaY_toppt_a_down_tt;
  TH2F *DeltaY_toppt_b_up_tt;
  TH2F *DeltaY_toppt_b_down_tt;


  uhh2::Event::Handle< std::vector<TopJet> > h_AK8TopTags;
  uhh2::Event::Handle< std::vector<Jet> > h_CHSjets_matched;
  uhh2::Event::Handle<bool> h_is_zprime_reconstructed_chi2;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
  uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;

  // --- NoAC support for tanh systematics (TT-only) ---
  bool use_noac_evtweights_ = false;
  std::string noac_gen_file_;
  std::string noac_gen_hist_;
  double noac_fraction_ = 0.0;
  uhh2::Event::Handle<float> h_xi_gen;
  std::unique_ptr<TH1D> noac_weights_;

  static std::unique_ptr<TH1D> mirror_hist_1d(const TH1D &src);
  static std::unique_ptr<TH1D> build_noac_weights_from_gen(const TH1D &Hgen, double f_noac);
  static double lookup_noac_weight(const TH1 *W, double xi);
  
  // Template-method (xi) multi-f infrastructure
  std::vector<float> f_values;
  static std::map<float, std::unique_ptr<TH1D>> noac_weights_map;  // Static: shared across all instances
  static bool noac_weights_initialized;  // Flag to track if weights have been initialized
  std::map<std::string, TH1F*> h_deltaY_xi_reco_map;
  virtual ~ZprimeSemiLeptonicSystematicsHists();
};