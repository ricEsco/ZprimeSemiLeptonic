#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <map>
#include <cmath>
#include <algorithm>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include "UHH2/common/include/PrintingModules.h"

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/JetCorrectionSets.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/core/include/GenParticle.h>
#include <set>
#include <UHH2/common/include/AdditionalSelections.h>
#include "UHH2/common/include/LuminosityHists.h"
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/MCWeight.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <stdexcept>
// #include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicPreselectionHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/CHSJetCorrections.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopPuppiJetCorrections.h>
#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"

using namespace std;
using namespace uhh2;

namespace {
  constexpr float kBFragNominalInputRb = 0.855f;

  // True if tau has e or mu in its decay chain (W->tau->e/mu). False if tau->hadronic or unknown.
  bool tau_decays_to_emumu(const std::vector<GenParticle>* gps, int tau_index, std::set<int>& visited) {
    if (!gps || gps->empty() || tau_index < 0) return false;
    if (visited.count(tau_index)) return false;
    visited.insert(tau_index);
    const unsigned short uinv = (unsigned short)(-1);
    for (const auto& gp : *gps) {
      const bool from_m1 = (gp.mother1() != uinv && (int)gp.mother1() == tau_index);
      const bool from_m2 = (gp.mother2() != uinv && (int)gp.mother2() == tau_index);
      if (!from_m1 && !from_m2) continue;
      const int id = std::abs(gp.pdgId());
      if (id == 11 || id == 13) return true;
      if (id == 15 && tau_decays_to_emumu(gps, gp.index(), visited)) return true;
    }
    return false;
  }

  // b-fragmentation : check if the particle is a b-hadron, b-quark, or b-hadron-quark.

  bool is_b_hadron_pdgid(int pdgid) {
    // This removes the sign, so both B hadrons and anti-B hadrons are accepted.
    const int id = std::abs(pdgid);
    // This removes simple particles like quarks and leptons.
    if (id <= 100) return false;
    // Then it checks whether the PDG ID contains a 5 in the hadron quark-content digits.
    return ((id / 100) % 10 == 5) || ((id / 1000) % 10 == 5) || ((id / 10000) % 10 == 5);
  }

  const GenParticle* find_genparticle_by_index(const std::vector<GenParticle>* gps, int index) {
    if (!gps || index < 0) return nullptr;
    for (const auto& gp : *gps) {
      if ((int)gp.index() == index) return &gp;
    }
    return nullptr;
  }

  bool has_ancestor_index(const std::vector<GenParticle>* gps, const GenParticle& gp, int ancestor_index, std::set<int>& visited) {
    if (!gps || ancestor_index < 0) return false;
    if (visited.count(gp.index())) return false;
    visited.insert(gp.index());

    const unsigned short uinv = (unsigned short)(-1);
    for (const auto mother_index : {gp.mother1(), gp.mother2()}) {
      if (mother_index == uinv) continue;
      if ((int)mother_index == ancestor_index) return true;
      const GenParticle* mother = find_genparticle_by_index(gps, mother_index);
      if (mother && has_ancestor_index(gps, *mother, ancestor_index, visited)) return true;
    }
    return false;
  }

  bool has_b_hadron_daughter(const std::vector<GenParticle>* gps, const GenParticle& gp) {
    if (!gps) return false;
    const unsigned short uinv = (unsigned short)(-1);
    for (const auto& candidate : *gps) {
      if (!is_b_hadron_pdgid(candidate.pdgId())) continue;
      const bool from_m1 = (candidate.mother1() != uinv && (int)candidate.mother1() == (int)gp.index());
      const bool from_m2 = (candidate.mother2() != uinv && (int)candidate.mother2() == (int)gp.index());
      if (from_m1 || from_m2) return true;
    }
    return false;
  }

  std::vector<const GenParticle*> find_last_b_hadrons_from_b(const std::vector<GenParticle>* gps, const GenParticle& b_quark) {
    std::vector<const GenParticle*> out;
    if (!gps || b_quark.index() < 0) return out;

    for (const auto& gp : *gps) {
      if (!is_b_hadron_pdgid(gp.pdgId())) continue;
      std::set<int> visited;
      if (!has_ancestor_index(gps, gp, b_quark.index(), visited)) continue;
      if (has_b_hadron_daughter(gps, gp)) continue; // MiniAOD fallback: last-copy B hadron
      out.push_back(&gp);
    }
    return out;
  }

  double lorentz_dot(const LorentzVector& a, const LorentzVector& b) {
    return a.E() * b.E() - a.Px() * b.Px() - a.Py() * b.Py() - a.Pz() * b.Pz();
  }

  // b-fragmentation : compute the xb value for the b-hadron: more detail: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MLReweighting?extralog=-%20caching%20topic#Bfragm
  float compute_bfrag_xb(const GenParticle& b_hadron, const GenParticle& top, const GenParticle& w_boson) {
    const double mt = top.v4().M();
    const double mw = w_boson.v4().M();
    if (!std::isfinite(mt) || mt <= 0.0 || !std::isfinite(mw)) return std::numeric_limits<float>::quiet_NaN();

    const double mt2 = mt * mt;
    const double w = (mw * mw) / mt2;
    const double denom = 1.0 - w;
    if (!std::isfinite(denom) || denom <= 0.0) return std::numeric_limits<float>::quiet_NaN();

    const double xE = 2.0 * lorentz_dot(b_hadron.v4(), top.v4()) / mt2;
    const double xb = xE / denom;
    if (!std::isfinite(xb)) return std::numeric_limits<float>::quiet_NaN();
    return static_cast<float>(xb);
  }

  struct BFragInputs {
    float xb_top = std::numeric_limits<float>::quiet_NaN();
    float xb_antitop = std::numeric_limits<float>::quiet_NaN();
    int valid = 0;
  };

  // b-fragmentation
  BFragInputs compute_bfrag_inputs(const std::vector<GenParticle>* gps, const TTbarGen& ttbargen) {
    BFragInputs out;
    if (!gps || gps->empty() || ttbargen.DecayChannel() == TTbarGen::e_notfound) return out;

    // ttbargen.bTop() bquark from t->bW+ 
    const auto b_hadrons_top = find_last_b_hadrons_from_b(gps, ttbargen.bTop());
    const auto b_hadrons_antitop = find_last_b_hadrons_from_b(gps, ttbargen.bAntitop());

    // check if there is only one b-hadron from the top and one b-hadron from the antitop
    if (b_hadrons_top.size() != 1 || b_hadrons_antitop.size() != 1) return out;
    if (b_hadrons_top.front()->index() == b_hadrons_antitop.front()->index()) return out;

    out.xb_top = compute_bfrag_xb(*b_hadrons_top.front(), ttbargen.Top(), ttbargen.WTop());
    out.xb_antitop = compute_bfrag_xb(*b_hadrons_antitop.front(), ttbargen.Antitop(), ttbargen.WAntitop());

    const bool xb_ok = std::isfinite(out.xb_top) && std::isfinite(out.xb_antitop)
                    && out.xb_top >= 0.f && out.xb_antitop >= 0.f
                    && out.xb_top <= 1.2f && out.xb_antitop <= 1.2f;
    out.valid = xb_ok ? 1 : 0;
    return out;
  }
  // EW uncertainty costheta info at gen level
  TLorentzVector to_tlorentzvector(const LorentzVector& v4) {
    TLorentzVector out;
    out.SetPxPyPzE(v4.Px(), v4.Py(), v4.Pz(), v4.E());
    return out;
  }

  float compute_top_costheta_ttbarframe(const TTbarGen& ttbargen) {
    if (ttbargen.DecayChannel() == TTbarGen::e_notfound) return std::numeric_limits<float>::quiet_NaN();

    const TLorentzVector top = to_tlorentzvector(ttbargen.Top().v4());
    const TLorentzVector antitop = to_tlorentzvector(ttbargen.Antitop().v4());
    const TLorentzVector ttbar = top + antitop;
    if (!std::isfinite(ttbar.E()) || ttbar.E() <= 0.0) return std::numeric_limits<float>::quiet_NaN();

    TLorentzVector top_ttbar_frame = top;
    top_ttbar_frame.Boost(-ttbar.BoostVector());
    if (!std::isfinite(top_ttbar_frame.P()) || top_ttbar_frame.P() <= 0.0) return std::numeric_limits<float>::quiet_NaN();

    double costheta = top_ttbar_frame.Vect().Unit().Dot(TVector3(0.0, 0.0, 1.0));
    if (!std::isfinite(costheta)) return std::numeric_limits<float>::quiet_NaN();
    costheta = std::max(-1.0, std::min(1.0, costheta));
    return static_cast<float>(costheta);
  }

  // NNLO QCD SF_3(pT), same parameters as TTbarNNLOQCDReweighting (ZprimeSemiLeptonicModules.cxx).
  // Used only for the no-cap cross-check (FillNNLONoCapGenHists); the capped value built from it
  // is compared on every event with the weight stored by the module.
  double nnlo_qcd_sf3(double pt) {
    return 0.17478303752138802 * std::exp(-0.0045665365705007399 * pt) + 0.89410719420857065;
  }

  // Local copy of the TTbarEWKCorrection lookup (ZprimeSemiLeptonicModules.cxx: clamp_to_axis,
  // evaluate_histo, evaluate_ratio and the cos(theta*) definition of process()), so that the EWK
  // weight can be evaluated on caller-supplied tops (last-copy cross-check). Checked on every event
  // against the module's stored weight on the status-22 tops.
  double ewk_clamp_to_axis(const TAxis* axis, double value) {
    if (!axis || !std::isfinite(value)) return std::numeric_limits<double>::quiet_NaN();
    const double xmin = axis->GetXmin();
    const double xmax = axis->GetXmax();
    const double eps = std::max(1e-9, 1e-9 * std::abs(xmax - xmin));
    if (value <= xmin) return xmin + eps;
    if (value >= xmax) return xmax - eps;
    return value;
  }
  double ewk_evaluate_histo(const TH2D* hist, double x, double y) {
    if (!hist) return 0.0;
    const double x_eval = ewk_clamp_to_axis(hist->GetXaxis(), x);
    const double y_eval = ewk_clamp_to_axis(hist->GetYaxis(), y);
    if (!std::isfinite(x_eval) || !std::isfinite(y_eval)) return 0.0;
    return hist->GetBinContent(hist->GetXaxis()->FindFixBin(x_eval), hist->GetYaxis()->FindFixBin(y_eval));
  }
  // EWYt_10 / EWno_0 at (m_ttbar, cos theta*) of the given tops; 1.0 wherever the module stores 1.0
  double ewk_weight_for(const TH2D* lo_qcd, const TH2D* nlo_ewk, const LorentzVector& top_v4, const LorentzVector& antitop_v4) {
    const double mttbar = (top_v4 + antitop_v4).M();
    const TLorentzVector top = to_tlorentzvector(top_v4);
    const TLorentzVector antitop = to_tlorentzvector(antitop_v4);
    const TLorentzVector ttbar = top + antitop;
    if (!(std::isfinite(mttbar) && mttbar > 0.0 && ttbar.E() > 0.0)) return 1.0;
    TLorentzVector top_ttbar_frame = top;
    top_ttbar_frame.Boost(-ttbar.BoostVector());
    if (!(top_ttbar_frame.P() > 0.0)) return 1.0;
    double costheta = top_ttbar_frame.Vect().Unit().Dot(TVector3(0.0, 0.0, 1.0));
    costheta = std::max(-1.0, std::min(1.0, costheta));
    const double lo = ewk_evaluate_histo(lo_qcd, mttbar, costheta);
    const double nlo = ewk_evaluate_histo(nlo_ewk, mttbar, costheta);
    if (lo == 0.0 || !std::isfinite(lo) || !std::isfinite(nlo)) return 1.0;
    const double ratio = nlo / lo;
    return (std::isfinite(ratio) && ratio > 0.0) ? ratio : 1.0;
  }

  // b-fragmentation : check if the particle is a b-hadron, b-quark, or b-hadron-quark.
  // True if any tau in the ttbar decay goes tau->hadronic. Exclude such events.
  bool has_any_hadronic_tau_decay(const std::vector<GenParticle>* gps, const TTbarGen& ttbargen) {
    if (!gps || gps->empty()) return true;  // conservative: exclude when uncertain
    const auto dc = ttbargen.DecayChannel();
    if (dc != TTbarGen::e_tauhad && dc != TTbarGen::e_tautau && dc != TTbarGen::e_etau && dc != TTbarGen::e_mutau)
      return false;
    std::set<int> visited;
    if (dc == TTbarGen::e_tauhad) {
      const auto& tau = ttbargen.ChargedLepton();
      if (!tau_decays_to_emumu(gps, tau.index(), visited)) return true;
      return false;
    }
    if (dc == TTbarGen::e_tautau) {
      for (const auto& wd : {ttbargen.Wdecay1(), ttbargen.Wdecay2(), ttbargen.WMinusdecay1(), ttbargen.WMinusdecay2()}) {
        if (std::abs(wd.pdgId()) == 15 && !tau_decays_to_emumu(gps, wd.index(), visited)) return true;
      }
      return false;
    }
    if (dc == TTbarGen::e_etau || dc == TTbarGen::e_mutau) {
      for (const auto& wd : {ttbargen.Wdecay1(), ttbargen.Wdecay2(), ttbargen.WMinusdecay1(), ttbargen.WMinusdecay2()}) {
        if (std::abs(wd.pdgId()) == 15) {
          if (!tau_decays_to_emumu(gps, wd.index(), visited)) return true;
          return false;
        }
      }
    }
    return false;
  }
}

class ZprimePreselectionModule : public ModuleBASE {

public:
  explicit ZprimePreselectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:
  bool debug;

  // mttbar mass bin lower edges (must match noac_mtt_edges in SystematicsHists for NoAC gen histograms)
  const std::vector<double> mttbar_bin_edges = {0., 500., 750., 1000., 1500.};
  const std::vector<std::string> mttbar_bin_tags = {"0_500", "500_750", "750_1000", "1000_1500", "1500Inf"};
  // Fine bins: [0-500), [500-750), [750-1000), [1000-1500), [1500, Inf)
  // Additionally book/fill mtt_gen_0_750 = merged [0, 750) (overlaps first two fine bins).

  // Helper function to find mttbar bin
  inline int find_mtt_bin(double mtt) {
    for (size_t i = 0; i+1 < mttbar_bin_edges.size(); ++i) {
      if (mtt >= mttbar_bin_edges[i] && mtt < mttbar_bin_edges[i+1]) return static_cast<int>(i);
    }
    if (mtt >= mttbar_bin_edges.back()) return static_cast<int>(mttbar_bin_edges.size() - 1);
    return -1;
  }

  // Standalone lumi weight (applied between raw GEN hists and lumi-scaled GEN hists)
  std::unique_ptr<MCLumiWeight> lumi_weight;

  // Optional NNLO QCD + EWK correction of the GEN xi histograms used to build NoAC weights.
  // Gated by XML flag ApplyNNLOEWKToGenHists (default false): when off, nothing below is
  // constructed, no extra output branches are declared, and the normal preselection ->
  // analysis -> DNN chain is completely unchanged. When on, a parallel "_lumiscaled_nnloewk"
  // folder set is filled with genWeight * lumi * EWK * NNLO-QCD, intended only for producing
  // the NoAC gen-template file (the actual per-event correction stays at the applyNN stage).
  bool apply_nnloewk_to_gen_hists = false;
  std::unique_ptr<AnalysisModule> ttbar_ewk_gen;
  std::unique_ptr<AnalysisModule> ttbar_nnloqcd_gen;
  uhh2::Event::Handle<float> h_w_ewk_nom;
  uhh2::Event::Handle<float> h_w_nnlo_nom;

  // Optional (XML FillNNLONoCapGenHists=true, needs ApplyNNLOEWKToGenHists=true and SF_3): parallel
  // "_lumiscaled_nnloewk_nocap" GEN folders (and gen_lastcopy/*_nnloewk_nocap) filled with
  // genWeight * lumi * EWK * NNLO-QCD WITHOUT the pT > 500 GeV freeze, to compare with the capped
  // _lumiscaled_nnloewk folders. The capped SF_3 recomputed here must equal the module's stored
  // weight; if not, the event keeps the stored weight and a warning is printed.
  bool fill_nnlo_nocap_gen_hists = false;
  long n_nnlo_nocap_mismatch = 0;

  // Corrections
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<AnalysisModule> hotvrjetCorr;
  std::unique_ptr<TopPuppiJetCorrections> toppuppijetCorr;
  std::unique_ptr<CHSJetCorrections> CHSjetCorr;

  // Cleaners
  std::unique_ptr<JetCleaner>      jet_IDcleaner, jet_cleaner1, jet_cleaner2;
  std::unique_ptr<AnalysisModule>  hotvrjet_cleaner;
  std::unique_ptr<TopJetCleaner>   topjet_puppi_IDcleaner, topjet_puppi_cleaner;

  // Selections
  std::unique_ptr<uhh2::Selection> genflavor_sel;
  std::unique_ptr<uhh2::Selection> jet1_sel;
  std::unique_ptr<uhh2::Selection> jet2_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  unique_ptr<Selection> SignSplit;

  bool isMC, isHOTVR;
  string Sys_PU;

  std::unique_ptr<Hists> lumihists;
  TString METcollection;

  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;

  // additional branch with AK4 CHS jets -> for b-tagging
  Event::Handle<vector<Jet>> h_CHSjets;

  // TTbarGen handle for mttbar calculation
  Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<TTbarGenProducer> ttgenprod;
  uhh2::Event::Handle<float> h_xi_gen;
  uhh2::Event::Handle<float> h_mtt_gen;
  uhh2::Event::Handle<float> h_costheta_gen;
  uhh2::Event::Handle<float> h_DeltaY_gen;
  // b-fragmentation 
  uhh2::Event::Handle<float> h_bfrag_xb_top;
  uhh2::Event::Handle<float> h_bfrag_xb_antitop;
  uhh2::Event::Handle<float> h_bfrag_rb;
  uhh2::Event::Handle<int> h_bfrag_valid;
  uhh2::Event::Handle<float> h_weight_bfrag_nom;
  uhh2::Event::Handle<float> h_weight_bfrag_up;

  // Cutflow histograms: survive hadd across Condor jobs
  // h_cf_raw:         bin content = N_raw (unweighted event count)
  // h_cf_genweight:   bin content = Sum(genWeight) before lumi scaling
  // h_cf_lumiscaled:  bin content = Sum(genWeight * lumiScale) after lumi scaling
  // Bin labels: 0_AllEvents, 1_AfterLumiWeight, 2_GEN_semilep, 3_AfterCommonModules,
  //             4_Lepton1, 5_Jet1, 6_Jet2, 7_MET
  TH1D *h_cf_raw, *h_cf_genweight, *h_cf_lumiscaled;
  static const int NCF = 8; // number of cutflow stages
  double lumi_factor = 1.0; // target_lumi / dataset_lumi, cached for cutflow

  // POWHEG EWK lookup ingredients (sigma_POWHEG vs m_ttbar, costheta*), filled with
  // raw genWeight before any selection and for ALL ttbar decays (no tau veto).
  // Enabled by setting PowhegEWKLookupBinningFile; binning cloned from EWno_0.
  // Harvest offline: sigma_POWHEG_bin = powheg_yield * xsec / powheg_sumw.
  TH2D* h_powheg_lookup_yield = nullptr;
  TH1D* h_powheg_lookup_sumw = nullptr;

  // Top-definition cross-check (XML flag FillLastCopyGenHists, default false).
  // Same GEN events and weight (genWeight * lumi) as the _lumiscaled folders, with xi and
  // m_tt computed twice: from the status-22 (hard-process) tops, and from the last top copy
  // (after ISR recoil / FSR, before decay; what TopDecaySubset kStable uses). The last copy
  // is pruned from UHH2 ntuples, but the stored W (status 22) and b (status 23) daughters of
  // the top come from its decay, so t_last = W + b.
  // Output directory gen_lastcopy/:
  //   mtt_22, mtt_last                     gen m_tt for each definition
  //   xi_22_<bin>, xi_last_<bin>           xi per gen-m_tt bin of the SAME definition
  //                                        (<bin> = inclusive, 0_500, ..., 1500Inf)
  //   migration_22_vs_last                 x = 2*mttbin_22 + (xi_22 > 0),
  //                                        y = 2*mttbin_last + (xi_last > 0); Sumw2 per cell
  //                                        gives the paired (same-event) error on A_last - A_22
  // Every histogram exists twice: without suffix (genWeight * lumi) and, when
  // ApplyNNLOEWKToGenHists is also on, with suffix _nnloewk (genWeight * lumi * EWK * NNLO-QCD,
  // the same per-event correction for both definitions; it is evaluated on the status-22 tops).
  struct LastCopyGenHists {
    TH1D* mtt_22 = nullptr;
    TH1D* mtt_last = nullptr;
    std::vector<TH1D*> xi_22;       // [0] = inclusive, [1 + i] = gen m_tt bin i
    std::vector<TH1D*> xi_last;
    TH2D* migration = nullptr;        // filled with the last-copy weight
    TH2D* migration_cross = nullptr;  // only *_ownkin sets: filled with w_22 * w_last
  };
  // Optional (XML FillLastCopyOwnKinGenHists=true; needs FillLastCopyGenHists and
  // ApplyNNLOEWKToGenHists): gen_lastcopy/*_nnloewk_ownkin (and *_nnloewk_nocap_ownkin when
  // FillNNLONoCapGenHists is on). The status-22 quantities keep the module weights; the last-copy
  // quantities get EWK and NNLO-QCD evaluated on the last-copy (W+b) tops, i.e. the dilepton
  // TopAnalysis recipe. The migration matrix is filled with the last-copy weight and
  // migration_*_cross with w_22 * w_last (paired error in compare_lastcopy_genac.py).
  bool fill_lastcopy_ownkin = false;
  std::unique_ptr<TH2D> h_ewk_lo_qcd;
  std::unique_ptr<TH2D> h_ewk_nlo_ewk;
  long n_ownkin_ewk_mismatch = 0;
  // DEBUG (XML DebugTauGen=true, default false): print the stored gen record of the first
  // DebugTauGenMaxEvents tau-channel events and running counts, to check whether tau decay
  // products are stored and whether leptonic-tau events pass fill_gen.
  bool debug_tau_gen = false;
  int  debug_tau_gen_max = 10;
  long n_dbg_tau_events = 0, n_dbg_tau_with_child = 0, n_dbg_tau_kept = 0, n_dbg_lep_printed = 0;

  bool fill_lastcopy_gen_hists = false;
  LastCopyGenHists lc_hists;           // genWeight * lumi
  LastCopyGenHists lc_hists_nnloewk;   // genWeight * lumi * EWK * NNLO-QCD
  LastCopyGenHists lc_hists_nnloewk_nocap;  // same, NNLO-QCD without the pT > 500 GeV freeze
  LastCopyGenHists lc_hists_ownkin;         // corrections on each definition's own tops (capped NNLO)
  LastCopyGenHists lc_hists_ownkin_nocap;   // same, uncapped NNLO

};

void ZprimePreselectionModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag+"_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicPreselectionHists(ctx,mytag));
  }
}

void ZprimePreselectionModule::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag+"_General";
  HFolder(mytag)->fill(event);
}

ZprimePreselectionModule::ZprimePreselectionModule(uhh2::Context& ctx) {

  debug = false; // true/false

  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }

  //// CONFIGURATION
  const TString METcollection = ctx.get("METName");
  isMC    = ctx.get("dataset_type") == "MC";
  isHOTVR = ctx.get("is_HOTVR") == "true";
  Sys_PU  = ctx.get("Sys_PU");

  isUL16preVFP  = (ctx.get("dataset_version").find("UL16preVFP")  != std::string::npos);
  isUL16postVFP = (ctx.get("dataset_version").find("UL16postVFP") != std::string::npos);
  isUL17        = (ctx.get("dataset_version").find("UL17")        != std::string::npos);
  isUL18        = (ctx.get("dataset_version").find("UL18")        != std::string::npos);

  // lepton IDs
  ElectronId eleID_veto = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90);
  MuonId     muID_veto  = MuonID(Muon::CutBasedIdTight);

  double electron_pt(25.);
  double muon_pt(25.);
  double jet1_pt(30.);
  double jet2_pt(30.);
  double MET(20.);

  // GEN Flavor selection [W+jets flavor-splitting]
  if(ctx.get("dataset_version").find("WJets") != std::string::npos){
    if     (ctx.get("dataset_version").find("_B") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("b"));
    else if(ctx.get("dataset_version").find("_C") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("c"));
    else if(ctx.get("dataset_version").find("_L") != std::string::npos) genflavor_sel.reset(new GenFlavorSelection("l"));
    else genflavor_sel.reset(new uhh2::AndSelection(ctx));
  }
  else genflavor_sel.reset(new uhh2::AndSelection(ctx));

  // Cleaning: Mu, Ele, Jets
  const MuonId muonID_veto(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID_veto));
  const ElectronId electronID_veto(AndId<Electron>(PtEtaSCCut(electron_pt, 2.5), eleID_veto));
  const JetPFID jetID_CHS(JetPFID::WP_TIGHT_CHS);
  const JetPFID jetID_PUPPI(JetPFID::WP_TIGHT_PUPPI);

  jet_IDcleaner.reset(new JetCleaner(ctx, jetID_PUPPI));
  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 20., 2.5));
  hotvrjet_cleaner.reset(new TopJetCleaner(ctx, PtEtaCut(200., 2.5)));
  topjet_puppi_IDcleaner.reset(new TopJetCleaner(ctx, jetID_PUPPI, "toppuppijets"));
  topjet_puppi_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(200., 2.5)), "toppuppijets"));

  // Split interference signal samples by sign
  if(ctx.get("dataset_version").find("_int") != std::string::npos){
    if     (ctx.get("dataset_version").find("_pos") != std::string::npos) SignSplit.reset(new SignSelection("pos"));
    else if(ctx.get("dataset_version").find("_neg") != std::string::npos) SignSplit.reset(new SignSelection("neg"));
    else SignSplit.reset(new uhh2::AndSelection(ctx));
  }
  else SignSplit.reset(new uhh2::AndSelection(ctx));

  // Standalone lumi weight: applied AFTER raw GEN hists, BEFORE lumi-scaled GEN hists
  if(isMC) lumi_weight.reset(new MCLumiWeight(ctx));

  // Optional: NNLO QCD + EWK correction applied (after lumi) to a parallel GEN folder set,
  // used only to build NoAC gen templates. Off by default so normal jobs are untouched.
  // Reuses the "ttbargen" handle produced above by TTbarGenProducer. Both modules apply to
  // event.weight only for ttbar datasets (gated internally on dataset_version).
  apply_nnloewk_to_gen_hists = isMC && (ctx.get("ApplyNNLOEWKToGenHists", "false") == "true");
  if(apply_nnloewk_to_gen_hists){
    ttbar_ewk_gen.reset(new TTbarEWKCorrection(ctx, "ttbargen"));
    ttbar_nnloqcd_gen.reset(new TTbarNNLOQCDReweighting(ctx, "ttbargen"));
    // Read back the stored nominal weights (set by the modules above) to apply only to the
    // GEN _nnloewk folders, without touching the global event.weight.
    h_w_ewk_nom  = ctx.get_handle<float>("weight_ttbar_ewk_nominal");
    h_w_nnlo_nom = ctx.get_handle<float>("weight_ttbar_nnlo_qcd");
  }
  fill_nnlo_nocap_gen_hists = isMC && (ctx.get("FillNNLONoCapGenHists", "false") == "true");
  if (fill_nnlo_nocap_gen_hists && !apply_nnloewk_to_gen_hists) {
    cout << "ZprimePreselectionModule WARNING: FillNNLONoCapGenHists=true needs ApplyNNLOEWKToGenHists=true; "
         << "no _nnloewk_nocap histograms will be filled." << endl;
    fill_nnlo_nocap_gen_hists = false;
  }
  if (fill_nnlo_nocap_gen_hists && ctx.get("TTbarNNLOQCDScaleFactor", "SF_3") != "SF_3") {
    cout << "ZprimePreselectionModule WARNING: FillNNLONoCapGenHists implements SF_3 only, but TTbarNNLOQCDScaleFactor="
         << ctx.get("TTbarNNLOQCDScaleFactor") << "; no _nnloewk_nocap histograms will be filled." << endl;
    fill_nnlo_nocap_gen_hists = false;
  }
  fill_lastcopy_ownkin = isMC && (ctx.get("FillLastCopyOwnKinGenHists", "false") == "true");
  if (fill_lastcopy_ownkin && !(apply_nnloewk_to_gen_hists && ctx.get("FillLastCopyGenHists", "false") == "true"
                                && ctx.get("TTbarNNLOQCDScaleFactor", "SF_3") == "SF_3")) {
    cout << "ZprimePreselectionModule WARNING: FillLastCopyOwnKinGenHists=true needs FillLastCopyGenHists=true, "
         << "ApplyNNLOEWKToGenHists=true and SF_3; no *_ownkin histograms will be filled." << endl;
    fill_lastcopy_ownkin = false;
  }
  if (fill_lastcopy_ownkin) {
    // same table and histograms as TTbarEWKCorrection
    const std::string ewk_file_path = ctx.get(
      "TTbarEWKCorrectionFile",
      "/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/top_pt_reweighting/HATHORMGPY_costheta.root");
    TFile ewk_file(locate_file(ewk_file_path).c_str(), "READ");
    if (ewk_file.IsZombie()) {
      throw std::runtime_error("ZprimePreselectionModule: failed to open TTbarEWKCorrectionFile '" + ewk_file_path + "'");
    }
    TH2D* lo_qcd = dynamic_cast<TH2D*>(ewk_file.Get("EWno_0"));
    TH2D* nlo_ewk = dynamic_cast<TH2D*>(ewk_file.Get("EWYt_10"));
    if (!lo_qcd || !nlo_ewk) {
      throw std::runtime_error("ZprimePreselectionModule: no EWno_0/EWYt_10 in '" + ewk_file_path + "'");
    }
    h_ewk_lo_qcd.reset(dynamic_cast<TH2D*>(lo_qcd->Clone("presel_ewk_lo_qcd")));
    h_ewk_nlo_ewk.reset(dynamic_cast<TH2D*>(nlo_ewk->Clone("presel_ewk_nlo_ewk")));
    h_ewk_lo_qcd->SetDirectory(0);
    h_ewk_nlo_ewk->SetDirectory(0);
    ewk_file.Close();
  }

  // Cache lumi factor for cutflow
  if(isMC){
    double dataset_lumi = std::abs(string2double(ctx.get("dataset_lumi")));
    double target_lumi  = string2double(ctx.get("target_lumi"));
    lumi_factor = target_lumi / dataset_lumi;
    cout << "=== LUMI SCALING INFO ===" << endl;
    cout << "  target_lumi  = " << target_lumi  << " pb^-1" << endl;
    cout << "  dataset_lumi = " << dataset_lumi  << " pb^-1" << endl;
    cout << "  lumi_factor  = " << lumi_factor << endl;
    cout << "=========================" << endl;
  }

  // common modules (lumi weight DISABLED - we apply it standalone between the two GEN hist sets)
  common.reset(new CommonModules());
  common->disable_mclumiweight();   // <-- avoid double lumi scaling
  common->switch_jetlepcleaner(true);
  // common->disable_pvfilter();
  common->disable_jetpfidfilter();
  common->switch_jetPtSorter(true);
  common->switch_metcorrection(true);
  common->set_muon_id(muonID_veto);
  common->set_electron_id(electronID_veto);
  common->init(ctx, Sys_PU);

  hotvrjetCorr.reset(new HOTVRJetCorrectionModule(ctx));

  toppuppijetCorr.reset(new TopPuppiJetCorrections());
  toppuppijetCorr->init(ctx);

  CHSjetCorr.reset(new CHSJetCorrections());
  CHSjetCorr->init(ctx);

  // TTbarGen producer. Never throw: this runs on every MC sample, and anything
  // without exactly one top and one antitop (WJets, ST, DY, QCD, toponium EtaT)
  // must yield DecayChannel == e_notfound instead of killing the job.
  if(isMC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));

  //// EVENT SELECTION
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.5))));
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.5))));
  met_sel.reset(new METCut(MET, uhh2::infinity));

  // additional branch with Ak4 CHS jets
  h_CHSjets = ctx.get_handle<vector<Jet>>("jetsAk4CHS");

  // TTbarGen handle for mttbar calculation
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  // GEN-level outputs (so they exist in the event and output tree)
  if (isMC) {
    h_xi_gen     = ctx.declare_event_output<float>("xi_gen");
    h_mtt_gen    = ctx.declare_event_output<float>("mtt_gen");
    h_costheta_gen = ctx.declare_event_output<float>("costheta_gen");
    h_DeltaY_gen = ctx.declare_event_output<float>("DeltaY_gen");
    // b-fragmentation 
    h_bfrag_xb_top      = ctx.declare_event_output<float>("bfrag_xb_top");
    h_bfrag_xb_antitop  = ctx.declare_event_output<float>("bfrag_xb_antitop");
    h_bfrag_rb          = ctx.declare_event_output<float>("bfrag_rb");
    h_bfrag_valid       = ctx.declare_event_output<int>("bfrag_valid");
    h_weight_bfrag_nom  = ctx.declare_event_output<float>("weight_bfrag_nom");
    h_weight_bfrag_up   = ctx.declare_event_output<float>("weight_bfrag_up");
  }

  // Book cutflow histograms (survive hadd across Condor jobs)
  {
    const vector<string> cf_labels = {
      "0_AllEvents", "1_GEN_raw", "2_AfterLumiWeight", "3_GEN_lumiscaled",
      "4_AfterCommonModules", "5_Lepton1", "6_Jet1", "7_MET"
    };
    auto book_cf = [&](const char* name, const char* title) -> TH1D* {
      TH1D* h = new TH1D(name, title, NCF, 0, NCF);
      h->Sumw2();
      for(int i = 0; i < NCF; ++i) h->GetXaxis()->SetBinLabel(i+1, cf_labels[i].c_str());
      ctx.put("cutflow", h);  // attach to output file under cutflow/ directory
      return h;
    };
    h_cf_raw        = book_cf("cf_raw",        "Cutflow: unweighted event count");
    h_cf_genweight  = book_cf("cf_genweight",  "Cutflow: sum of genWeight (before lumi)");
    h_cf_lumiscaled = book_cf("cf_lumiscaled", "Cutflow: sum of lumi-scaled weight");
  }

  // Book histograms
  // --- Raw GEN folders (genWeight only, no lumi) ---
  vector<string> histogram_tags = {"Input", "mtt_gen_inclusive"};
    // "CommonModules", "HOTVRCorrections", "PUPPICorrections", "Lepton1", "JetID", "JetCleaner1", "JetCleaner2", "TopjetCleaner", "Jet1", "Jet2", "MET"};
    
  // Add mttbar bin tags
  for (const auto & tag : mttbar_bin_tags) {
    histogram_tags.push_back("mtt_gen_" + tag);
  }
  histogram_tags.push_back("mtt_gen_0_750");

  // --- Lumi-scaled GEN folders (genWeight * lumi, no PU) ---
  histogram_tags.push_back("mtt_gen_inclusive_lumiscaled");
  for (const auto & tag : mttbar_bin_tags) {
    histogram_tags.push_back("mtt_gen_" + tag + "_lumiscaled");
  }
  histogram_tags.push_back("mtt_gen_0_750_lumiscaled");
  // --- NNLO QCD + EWK corrected GEN folders (genWeight * lumi * EWK * NNLOQCD) ---
  // Parallel to the _lumiscaled set above; only booked when ApplyNNLOEWKToGenHists is on.
  if(apply_nnloewk_to_gen_hists){
    histogram_tags.push_back("mtt_gen_inclusive_lumiscaled_nnloewk");
    for (const auto & tag : mttbar_bin_tags) {
      histogram_tags.push_back("mtt_gen_" + tag + "_lumiscaled_nnloewk");
    }
    histogram_tags.push_back("mtt_gen_0_750_lumiscaled_nnloewk");
  }
  // --- Same, with the NNLO-QCD SF not frozen above pT = 500 GeV (FillNNLONoCapGenHists) ---
  if(fill_nnlo_nocap_gen_hists){
    histogram_tags.push_back("mtt_gen_inclusive_lumiscaled_nnloewk_nocap");
    for (const auto & tag : mttbar_bin_tags) {
      histogram_tags.push_back("mtt_gen_" + tag + "_lumiscaled_nnloewk_nocap");
    }
    histogram_tags.push_back("mtt_gen_0_750_lumiscaled_nnloewk_nocap");
  }

  // --- Reco-level folders ---
  histogram_tags.push_back("MET");

  book_histograms(ctx, histogram_tags);

  // POWHEG EWK lookup: book yield (HATHOR binning) + sum-of-genWeight histograms.
  // Off unless PowhegEWKLookupBinningFile is set in the XML, so normal runs are unchanged.
  const std::string powheg_lookup_file = ctx.get("PowhegEWKLookupBinningFile", "");
  if (isMC && !powheg_lookup_file.empty()) {
    const std::string powheg_lookup_hist = ctx.get("PowhegEWKLookupBinningHist", "EWno_0");
    TFile lookup_file(powheg_lookup_file.c_str(), "READ");
    if (lookup_file.IsZombie()) {
      throw std::runtime_error("ZprimePreselectionModule: failed to open PowhegEWKLookupBinningFile '" + powheg_lookup_file + "'");
    }
    TH2D* lookup_src = dynamic_cast<TH2D*>(lookup_file.Get(powheg_lookup_hist.c_str()));
    if (!lookup_src) {
      throw std::runtime_error("ZprimePreselectionModule: no TH2 '" + powheg_lookup_hist + "' in '" + powheg_lookup_file + "'");
    }
    h_powheg_lookup_yield = dynamic_cast<TH2D*>(lookup_src->Clone("powheg_yield"));
    h_powheg_lookup_yield->SetDirectory(0);
    h_powheg_lookup_yield->Reset("ICES");
    if (h_powheg_lookup_yield->GetSumw2N() == 0) h_powheg_lookup_yield->Sumw2();
    h_powheg_lookup_yield->SetTitle("POWHEG yield, raw genWeight;M_{t#bar{t}} [GeV];cos#theta*");
    h_powheg_lookup_sumw = new TH1D("powheg_sumw", "sum of genWeight, all events", 1, 0., 1.);
    h_powheg_lookup_sumw->Sumw2();
    ctx.put("powheg_ewk_lookup", h_powheg_lookup_yield);
    ctx.put("powheg_ewk_lookup", h_powheg_lookup_sumw);
    lookup_file.Close();
  }

  debug_tau_gen     = isMC && (ctx.get("DebugTauGen", "false") == "true");
  debug_tau_gen_max = std::stoi(ctx.get("DebugTauGenMaxEvents", "10"));

  // Top-definition cross-check: status-22 vs last-copy tops (see member comment).
  // Off unless FillLastCopyGenHists=true in the XML, so normal runs are unchanged.
  fill_lastcopy_gen_hists = isMC && (ctx.get("FillLastCopyGenHists", "false") == "true");
  if (fill_lastcopy_gen_hists) {
    auto book_lc = [&](TH1* h) {
      h->Sumw2();
      ctx.put("gen_lastcopy", h);  // attach to output file under gen_lastcopy/ directory
    };
    // 300 xi bins: rebinnable to 100 (NoAC templates), 12 and 6 (comparisons)
    std::vector<std::string> xi_tags = {"inclusive"};
    xi_tags.insert(xi_tags.end(), mttbar_bin_tags.begin(), mttbar_bin_tags.end());
    const int ncell = 2 * static_cast<int>(mttbar_bin_tags.size());

    auto book_set = [&](const std::string& sfx, const std::string& wlabel, bool with_cross) {
      LastCopyGenHists s;
      s.mtt_22   = new TH1D(("mtt_22" + sfx).c_str(),   ("GEN m_{t#bar{t}}, status-22 tops, " + wlabel + ";m_{t#bar{t}} [GeV]").c_str(), 500, 0., 5000.);
      s.mtt_last = new TH1D(("mtt_last" + sfx).c_str(), ("GEN m_{t#bar{t}}, last-copy tops (W+b), " + wlabel + ";m_{t#bar{t}} [GeV]").c_str(), 500, 0., 5000.);
      book_lc(s.mtt_22);
      book_lc(s.mtt_last);
      for (const auto& tag : xi_tags) {
        TH1D* h22   = new TH1D(("xi_22_" + tag + sfx).c_str(),   ("#xi GEN, status-22 tops, m_{t#bar{t}} " + tag + ", " + wlabel).c_str(),         300, -1., 1.);
        TH1D* hlast = new TH1D(("xi_last_" + tag + sfx).c_str(), ("#xi GEN, last-copy tops (W+b), m_{t#bar{t}} " + tag + ", " + wlabel).c_str(), 300, -1., 1.);
        book_lc(h22);
        book_lc(hlast);
        s.xi_22.push_back(h22);
        s.xi_last.push_back(hlast);
      }
      s.migration = new TH2D(("migration_22_vs_last" + sfx).c_str(),
                             ("GEN (m_{t#bar{t}} bin, sign #xi): status-22 (x) vs last copy (y), " + wlabel +
                              ";2*bin_{22}+(#xi_{22}>0);2*bin_{last}+(#xi_{last}>0)").c_str(),
                             ncell, 0., ncell, ncell, 0., ncell);
      book_lc(s.migration);
      if (with_cross) {
        s.migration_cross = new TH2D(("migration_22_vs_last" + sfx + "_cross").c_str(),
                                     ("sum(w_{22} w_{last}) per cell, " + wlabel +
                                      ";2*bin_{22}+(#xi_{22}>0);2*bin_{last}+(#xi_{last}>0)").c_str(),
                                     ncell, 0., ncell, ncell, 0., ncell);
        book_lc(s.migration_cross);
      }
      return s;
    };
    lc_hists = book_set("", "genWeight*lumi", false);
    if (apply_nnloewk_to_gen_hists) {
      lc_hists_nnloewk = book_set("_nnloewk", "genWeight*lumi*EWK*NNLO", false);
      if (fill_nnlo_nocap_gen_hists) {
        lc_hists_nnloewk_nocap = book_set("_nnloewk_nocap", "genWeight*lumi*EWK*NNLO (no pT cap)", false);
      }
      if (fill_lastcopy_ownkin) {
        lc_hists_ownkin = book_set("_nnloewk_ownkin", "genWeight*lumi*EWK*NNLO, last copy: corrections on W+b tops", true);
        if (fill_nnlo_nocap_gen_hists) {
          lc_hists_ownkin_nocap = book_set("_nnloewk_nocap_ownkin", "genWeight*lumi*EWK*NNLO (no pT cap), last copy: corrections on W+b tops", true);
        }
      }
    }
    else {
      cout << "ZprimePreselectionModule WARNING: FillLastCopyGenHists=true but ApplyNNLOEWKToGenHists is not true: "
           << "only the UNCORRECTED (genWeight*lumi) gen_lastcopy histograms will be filled." << endl;
    }
  }

  lumihists.reset(new LuminosityHists(ctx, "lumi"));
}

bool ZprimePreselectionModule::process(uhh2::Event& event){

  // Process TTbarGen first
  if (isMC && ttgenprod) {
    ttgenprod->process(event);
  }

  // b-fragmentation
  if (isMC) {
    BFragInputs bfrag_inputs;
    if (event.is_valid(h_ttbargen)) {
      bfrag_inputs = compute_bfrag_inputs(event.genparticles, event.get(h_ttbargen));
    }

    event.set(h_bfrag_xb_top, bfrag_inputs.xb_top);
    event.set(h_bfrag_xb_antitop, bfrag_inputs.xb_antitop);
    event.set(h_bfrag_rb, kBFragNominalInputRb);
    event.set(h_bfrag_valid, bfrag_inputs.valid);

    // Real b-fragmentation weights require the trained ONNX models. Keep the
    // branches neutral until model inference or official plugin output is wired.
    event.set(h_weight_bfrag_nom, 1.0f);
    event.set(h_weight_bfrag_up, 1.0f);

    // NNLO QCD + EWK gen-template correction. These modules declare output branches that SFrame
    // requires to be set on EVERY written event, so process() must run for all MC events here
    // (the modules set their weights to 1.0 for non-ttbar). event.weight is restored immediately;
    // the correction itself is applied only to the dedicated _nnloewk GEN folders below, by reading
    // the stored nominal weights from the handles.
    if (apply_nnloewk_to_gen_hists) {
      const double w_keep = event.weight;
      ttbar_ewk_gen->process(event);
      ttbar_nnloqcd_gen->process(event);
      event.weight = w_keep;
    }
  }

  if(debug) cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug) cout << " run.event: " << event.run << ". " << event.event << endl;
  if(debug) cout << " event.year: " << event.year << ". " << event.event << endl;

  if(!event.isRealData){
    if(!SignSplit->passes(event)) return false;
  }
  if(debug) cout << "beginning: ok" << endl;

  // ==========================================
  // At this point: event.weight = genInfo->weights()[0] (set by framework)
  // NO lumi scaling yet, NO PU, NO reco weights
  // ==========================================

  // --- CUTFLOW: 0_AllEvents ---
  h_cf_raw->Fill(0.5, 1.);
  h_cf_genweight->Fill(0.5, event.weight);
  h_cf_lumiscaled->Fill(0.5, event.weight * lumi_factor);

  // POWHEG EWK lookup: every MC event counts toward the cross-section denominator,
  // independent of decay channel or any later selection (event.weight is still the
  // raw genWeight here).
  if (h_powheg_lookup_sumw && !event.isRealData) {
    h_powheg_lookup_sumw->Fill(0.5, event.weight);
  }

  fill_histograms(event, "Input");
  if(debug) cout << "first plots input: ok" << endl;

  // ==========================================
  // FIRST: fill raw GEN mttbar folders (genWeight only — no lumi)
  // ==========================================
  if (isMC) {
    // Defaults are always set for MC, even if h_ttbargen is not valid
    event.set(h_xi_gen,       std::numeric_limits<float>::quiet_NaN());
    event.set(h_mtt_gen,      std::numeric_limits<float>::quiet_NaN());
    event.set(h_costheta_gen, std::numeric_limits<float>::quiet_NaN());
    event.set(h_DeltaY_gen,   std::numeric_limits<float>::quiet_NaN());

    bool gen_filled = false;
    if (event.is_valid(h_ttbargen)) {
      const auto& ttbargen = event.get(h_ttbargen);
      const auto dc = ttbargen.DecayChannel();

      // POWHEG EWK lookup: fill before any selection, for ALL decay channels
      // (deliberately no hadronic-tau veto — sigma_POWHEG is decay-inclusive).
      // Out-of-range values are clamped into the first/last bin, mirroring how
      // TTbarEWKCorrection::evaluate_histo clamps at evaluation time.
      if (h_powheg_lookup_yield && dc != TTbarGen::e_notfound) {
        const double mtt_lookup = (ttbargen.Top().v4() + ttbargen.Antitop().v4()).M();
        const double ct_lookup = compute_top_costheta_ttbarframe(ttbargen);
        if (std::isfinite(mtt_lookup) && std::isfinite(ct_lookup)) {
          const TAxis* ax = h_powheg_lookup_yield->GetXaxis();
          const TAxis* ay = h_powheg_lookup_yield->GetYaxis();
          h_powheg_lookup_yield->Fill(
            std::min(std::max(mtt_lookup, ax->GetBinCenter(1)), ax->GetBinCenter(ax->GetNbins())),
            std::min(std::max(ct_lookup, ay->GetBinCenter(1)), ay->GetBinCenter(ay->GetNbins())),
            event.weight);
        }
      }
      const bool fill_gen = (dc != TTbarGen::e_notfound && !has_any_hadronic_tau_decay(event.genparticles, ttbargen));

      // DEBUG: tau treatment at gen level (XML DebugTauGen=true)
      if (debug_tau_gen && event.genparticles) {
        static const char* dc_names[] = {"e_had", "e_ehad", "e_muhad", "e_tauhad", "e_ee", "e_mumu",
                                         "e_tautau", "e_emu", "e_etau", "e_mutau", "e_notfound"};
        const auto& gps = *event.genparticles;
        const bool tau_channel = (dc == TTbarGen::e_tauhad || dc == TTbarGen::e_tautau ||
                                  dc == TTbarGen::e_etau   || dc == TTbarGen::e_mutau);
        const bool lep_channel = (dc == TTbarGen::e_ehad || dc == TTbarGen::e_muhad);

        // does any stored particle have a tau as mother? (exactly what tau_decays_to_emumu scans)
        auto is_tau_child = [&](const GenParticle& gp) {
          for (const auto& tau : gps) {
            if (std::abs(tau.pdgId()) != 15) continue;
            if (gp.mother1() == tau.index() || gp.mother2() == tau.index()) return true;
          }
          return false;
        };
        bool any_tau_child = false;
        for (const auto& gp : gps) any_tau_child |= is_tau_child(gp);

        if (tau_channel) {
          ++n_dbg_tau_events;
          if (any_tau_child) ++n_dbg_tau_with_child;
          if (fill_gen) ++n_dbg_tau_kept;
        }
        const bool print_event = (tau_channel && n_dbg_tau_events <= debug_tau_gen_max) ||
                                 (lep_channel && n_dbg_lep_printed < 2);   // two e/mu+jets events for comparison
        if (lep_channel && print_event) ++n_dbg_lep_printed;

        if (print_event) {
          const auto old_flags = cout.flags();
          const auto old_prec  = cout.precision();
          cout << "\n[DebugTauGen] run " << event.run << " event " << event.event
               << "  decay channel = " << dc_names[dc] << " (" << dc << ")"
               << "  fill_gen = " << fill_gen << (fill_gen ? "  (KEPT)" : "  (EXCLUDED)")
               << "  any particle with a tau mother: " << (any_tau_child ? "YES" : "NO") << endl;
          cout << "  (mother/daughter index 65535 or -1 = none)" << endl;
          cout << setw(6) << "idx" << setw(8) << "pdgId" << setw(8) << "status"
               << setw(9) << "mother1" << setw(9) << "mother2" << setw(11) << "daughter1" << setw(11) << "daughter2"
               << setw(10) << "pt" << setw(8) << "eta" << endl;
          for (const auto& gp : gps) {
            cout << setw(6) << gp.index() << setw(8) << gp.pdgId() << setw(8) << gp.status()
                 << setw(9) << gp.mother1() << setw(9) << gp.mother2() << setw(11) << gp.daughter1() << setw(11) << gp.daughter2()
                 << setw(10) << fixed << setprecision(1) << gp.pt();
            if (gp.pt() > 1e-3) cout << setw(8) << setprecision(2) << gp.eta();
            else                cout << setw(8) << "-";     // incoming partons: pt = 0, eta undefined
            if (std::abs(gp.pdgId()) == 15) cout << "   <-- tau";
            if (is_tau_child(gp))           cout << "   <-- CHILD OF A TAU";
            cout << endl;
          }
          cout.flags(old_flags);
          cout.precision(old_prec);
        }
        if (tau_channel && (print_event || n_dbg_tau_events % 500 == 0)) {
          cout << "[DebugTauGen] tau-channel events so far: " << n_dbg_tau_events
               << " | with any stored tau daughter: " << n_dbg_tau_with_child
               << " | passing fill_gen (kept): " << n_dbg_tau_kept << endl;
        }
      }

      if (fill_gen) {
        const auto& top  = ttbargen.Top();
        const auto& atop = ttbargen.Antitop();
        double mtt = (top.v4() + atop.v4()).M();
        double dy  = std::abs(top.v4().Rapidity()) - std::abs(atop.v4().Rapidity());
        event.set(h_xi_gen,     std::tanh(dy));
        event.set(h_mtt_gen,    static_cast<float>(mtt));
        event.set(h_costheta_gen, compute_top_costheta_ttbarframe(ttbargen));
        event.set(h_DeltaY_gen, static_cast<float>(dy));

        // --- CUTFLOW: 1_GEN_raw ---
        h_cf_raw->Fill(1.5, 1.);
        h_cf_genweight->Fill(1.5, event.weight);

        // Fill RAW (no lumi) GEN folders
        fill_histograms(event, "mtt_gen_inclusive");
        const int ibin = find_mtt_bin(mtt);
        if (ibin >= 0) {
          fill_histograms(event, "mtt_gen_" + mttbar_bin_tags.at(ibin));
        }
        // Merged [0, 750) GeV: same events as fine bins 0-500 and 500-750
        if (mtt >= mttbar_bin_edges[0] && mtt < mttbar_bin_edges[2]) {
          fill_histograms(event, "mtt_gen_0_750");
        }

        // ==========================================
        // NOW: apply lumi scaling (same events as raw GEN fills above)
        // ==========================================
        lumi_weight->process(event);
        // event.weight = genWeight * (target_lumi / dataset_lumi)

        // Uncapped NNLO-QCD weight for the optional _nnloewk_nocap histograms (status-22 tops).
        // nnlo_sf3_matches_module: the capped SF_3 recomputed here equals the module's stored weight,
        // i.e. the module applied SF_3 to this event (false for non-ttbar datasets, where it stores 1.0).
        double w_nnlo_nocap = 1.0;
        bool nnlo_sf3_matches_module = false;
        if (fill_nnlo_nocap_gen_hists || fill_lastcopy_ownkin) {
          const double pt_t  = top.v4().Pt();
          const double pt_tb = atop.v4().Pt();
          const double w_nnlo_stored = event.get(h_w_nnlo_nom);
          const double w_nnlo_capped = std::sqrt(nnlo_qcd_sf3(std::min(pt_t, 500.)) * nnlo_qcd_sf3(std::min(pt_tb, 500.)));
          nnlo_sf3_matches_module = (std::abs(w_nnlo_capped - w_nnlo_stored) <= 1e-6 * std::abs(w_nnlo_capped));
          if (nnlo_sf3_matches_module) {
            w_nnlo_nocap = std::sqrt(nnlo_qcd_sf3(pt_t) * nnlo_qcd_sf3(pt_tb));
          }
          else {
            w_nnlo_nocap = w_nnlo_stored;  // module did not apply the capped SF_3 here: keep its weight
            if (++n_nnlo_nocap_mismatch <= 5) {
              cerr << "ZprimePreselectionModule WARNING: capped SF_3 " << w_nnlo_capped << " != stored NNLO weight "
                   << w_nnlo_stored << " (event " << event.event << "); using the stored weight for _nnloewk_nocap / _ownkin." << endl;
            }
          }
        }

        // --- CUTFLOW: 2_AfterLumiWeight, 3_GEN_lumiscaled ---
        h_cf_raw->Fill(2.5, 1.);
        h_cf_lumiscaled->Fill(2.5, event.weight);
        h_cf_raw->Fill(3.5, 1.);
        h_cf_lumiscaled->Fill(3.5, event.weight);

        // Fill LUMI-SCALED GEN folders (same events, now with lumi weight)
        fill_histograms(event, "mtt_gen_inclusive_lumiscaled");
        if (ibin >= 0) {
          fill_histograms(event, "mtt_gen_" + mttbar_bin_tags.at(ibin) + "_lumiscaled");
        }
        if (mtt >= mttbar_bin_edges[0] && mtt < mttbar_bin_edges[2]) {
          fill_histograms(event, "mtt_gen_0_750_lumiscaled");
        }

        // OPTIONAL: top-definition cross-check, status-22 vs last-copy (W+b) tops.
        // Same events as the _lumiscaled folders above; event.weight = genWeight * lumi here.
        if (fill_lastcopy_gen_hists) {
          const auto t_last  = ttbargen.WTop().v4()     + ttbargen.bTop().v4();
          const auto tb_last = ttbargen.WAntitop().v4() + ttbargen.bAntitop().v4();
          const double mtt_last  = (t_last + tb_last).M();
          const double xi_last   = std::tanh(std::abs(t_last.Rapidity()) - std::abs(tb_last.Rapidity()));
          const double xi_22     = std::tanh(dy);
          const int    ibin_last = find_mtt_bin(mtt_last);

          // w22: weight for the status-22 quantities, wlast: weight for the last-copy quantities
          auto fill_set = [&](LastCopyGenHists& s, double w22, double wlast) {
            s.mtt_22->Fill(mtt, w22);
            s.mtt_last->Fill(mtt_last, wlast);
            s.xi_22.at(0)->Fill(xi_22, w22);
            s.xi_last.at(0)->Fill(xi_last, wlast);
            if (ibin >= 0)      s.xi_22.at(ibin + 1)->Fill(xi_22, w22);
            if (ibin_last >= 0) s.xi_last.at(ibin_last + 1)->Fill(xi_last, wlast);
            if (ibin >= 0 && ibin_last >= 0) {
              const double cell_22   = 2 * ibin      + (xi_22   > 0. ? 1 : 0) + 0.5;
              const double cell_last = 2 * ibin_last + (xi_last > 0. ? 1 : 0) + 0.5;
              s.migration->Fill(cell_22, cell_last, wlast);
              if (s.migration_cross) s.migration_cross->Fill(cell_22, cell_last, w22 * wlast);
            }
          };
          fill_set(lc_hists, event.weight, event.weight);
          if (apply_nnloewk_to_gen_hists) {
            // same per-event correction as the _lumiscaled_nnloewk folders (status-22 tops)
            const double w_ewk_22  = event.get(h_w_ewk_nom);
            const double w_nnlo_22 = event.get(h_w_nnlo_nom);
            const double w22_cap   = event.weight * w_ewk_22 * w_nnlo_22;
            const double w22_nocap = event.weight * w_ewk_22 * w_nnlo_nocap;
            fill_set(lc_hists_nnloewk, w22_cap, w22_cap);
            if (fill_nnlo_nocap_gen_hists) {
              fill_set(lc_hists_nnloewk_nocap, w22_nocap, w22_nocap);
            }
            if (fill_lastcopy_ownkin) {
              // EWK on the last-copy tops, after checking that the local lookup reproduces the module
              // weight on the status-22 tops; otherwise keep the module weight (non-ttbar datasets).
              double w_ewk_last = ewk_weight_for(h_ewk_lo_qcd.get(), h_ewk_nlo_ewk.get(), t_last, tb_last);
              const double w_ewk_22_local = ewk_weight_for(h_ewk_lo_qcd.get(), h_ewk_nlo_ewk.get(), top.v4(), atop.v4());
              if (std::abs(w_ewk_22_local - w_ewk_22) > 1e-6 * std::abs(w_ewk_22)) {
                w_ewk_last = w_ewk_22;
                if (++n_ownkin_ewk_mismatch <= 5) {
                  cerr << "ZprimePreselectionModule WARNING: local EWK lookup " << w_ewk_22_local << " != stored EWK weight "
                       << w_ewk_22 << " (event " << event.event << "); using the stored weight for *_ownkin." << endl;
                }
              }
              // NNLO-QCD on the last-copy tops (capped / uncapped); module weight if SF_3 did not apply
              const double pt_tl  = t_last.Pt();
              const double pt_tbl = tb_last.Pt();
              const double w_nnlo_last_cap   = nnlo_sf3_matches_module
                ? std::sqrt(nnlo_qcd_sf3(std::min(pt_tl, 500.)) * nnlo_qcd_sf3(std::min(pt_tbl, 500.))) : w_nnlo_22;
              const double w_nnlo_last_nocap = nnlo_sf3_matches_module
                ? std::sqrt(nnlo_qcd_sf3(pt_tl) * nnlo_qcd_sf3(pt_tbl)) : w_nnlo_22;
              fill_set(lc_hists_ownkin, w22_cap, event.weight * w_ewk_last * w_nnlo_last_cap);
              if (fill_nnlo_nocap_gen_hists) {
                fill_set(lc_hists_ownkin_nocap, w22_nocap, event.weight * w_ewk_last * w_nnlo_last_nocap);
              }
            }
          }
        }

        // ==========================================
        // OPTIONAL: NNLO QCD + EWK corrected GEN folders (applied AFTER lumi scaling).
        // Only used to build NoAC gen templates; event.weight is restored afterwards so the
        // downstream reco-level processing (CommonModules onward) sees the lumi-only weight.
        // ==========================================
        if (apply_nnloewk_to_gen_hists) {
          const double w_lumi = event.weight;          // genWeight * lumi
          const float w_ewk  = event.get(h_w_ewk_nom);   // EWK ratio (1.0 if not ttbar)
          const float w_nnlo = event.get(h_w_nnlo_nom);  // sqrt(SF(pt_top)*SF(pt_atop))
          event.weight = w_lumi * w_ewk * w_nnlo;
          fill_histograms(event, "mtt_gen_inclusive_lumiscaled_nnloewk");
          if (ibin >= 0) {
            fill_histograms(event, "mtt_gen_" + mttbar_bin_tags.at(ibin) + "_lumiscaled_nnloewk");
          }
          if (mtt >= mttbar_bin_edges[0] && mtt < mttbar_bin_edges[2]) {
            fill_histograms(event, "mtt_gen_0_750_lumiscaled_nnloewk");
          }
          if (fill_nnlo_nocap_gen_hists) {
            event.weight = w_lumi * w_ewk * w_nnlo_nocap;
            fill_histograms(event, "mtt_gen_inclusive_lumiscaled_nnloewk_nocap");
            if (ibin >= 0) {
              fill_histograms(event, "mtt_gen_" + mttbar_bin_tags.at(ibin) + "_lumiscaled_nnloewk_nocap");
            }
            if (mtt >= mttbar_bin_edges[0] && mtt < mttbar_bin_edges[2]) {
              fill_histograms(event, "mtt_gen_0_750_lumiscaled_nnloewk_nocap");
            }
          }
          event.weight = w_lumi;                        // restore lumi-only weight
        }
        gen_filled = true;
      }
    }

    if (!gen_filled) {
      // Non-semileptonic, hadronic-tau, or no ttbargen: still apply lumi weight for reco-level
      lumi_weight->process(event);
    }
  }

  // ==========================================
  // CommonModules: PU reweight + JEC/JER + cleaners
  // (lumi already applied above, disabled in CommonModules)
  // ==========================================
  bool commonResult = common->process(event);
  if (!commonResult) return false;

  // --- CUTFLOW: 4_AfterCommonModules (lumi + PU) ---
  h_cf_raw->Fill(4.5, 1.);
  h_cf_lumiscaled->Fill(4.5, event.weight);

  if(debug) cout << "CommonModules: ok" << endl;
  // fill_histograms(event, "CommonModules");

  sort_by_pt<Muon>(*event.muons);
  sort_by_pt<Electron>(*event.electrons);

  // Correct AK4 CHS jets
  CHSjetCorr->process(event);

  if(isHOTVR){
    hotvrjetCorr->process(event);
    // fill_histograms(event, "HOTVRCorrections");
  }

  toppuppijetCorr->process(event);
  if(debug) cout << "TopPuppiJetCorrections: ok" << endl;
  // fill_histograms(event, "PUPPICorrections");

  // GEN ME quark-flavor selection
  if(!event.isRealData){
    if(!genflavor_sel->passes(event)) return false;
  }
  if(debug) cout << "GenFlavorSelection: ok" << endl;

  const bool pass_lep1 = ((event.muons->size() >= 1) || (event.electrons->size() >= 1));
  if(!pass_lep1) return false;
  h_cf_raw->Fill(5.5, 1.); h_cf_lumiscaled->Fill(5.5, event.weight);
  if(debug) cout << "≥1 leptons: ok" << endl;
  // fill_histograms(event, "Lepton1");

  jet_IDcleaner->process(event);
  // fill_histograms(event, "JetID");
  if(debug) cout << "JetCleaner ID: ok" << endl;

  jet_cleaner1->process(event);
  sort_by_pt<Jet>(*event.jets);
  // fill_histograms(event, "JetCleaner1");
  if(debug) cout << "JetCleaner1: ok" << endl;

  // Lepton-2Dcut variables
  for(auto& muo : *event.muons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
    muo.set_tag(Muon::twodcut_dRmin, dRmin);
    muo.set_tag(Muon::twodcut_pTrel, pTrel);
  }

  for(auto& ele : *event.electrons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);
    ele.set_tag(Electron::twodcut_dRmin, dRmin);
    ele.set_tag(Electron::twodcut_pTrel, pTrel);
  }

  jet_cleaner2->process(event);
  sort_by_pt<Jet>(*event.jets);
  // fill_histograms(event, "JetCleaner2");
  if(debug) cout << "JetCleaner2: ok" << endl;

  hotvrjet_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);

  topjet_puppi_IDcleaner->process(event);
  topjet_puppi_cleaner->process(event);
  sort_by_pt<TopJet>(*event.toppuppijets);

  // fill_histograms(event, "TopjetCleaner");
  if(debug) cout << "TopJetCleaner: ok" << endl;

  // 1st AK4 jet selection
  const bool pass_jet1 = jet1_sel->passes(event);
  if(!pass_jet1) return false;
  h_cf_raw->Fill(6.5, 1.); h_cf_lumiscaled->Fill(6.5, event.weight);

  if(debug) cout << "NJetSelection1: ok" << endl;
  // fill_histograms(event, "Jet1");

  // 2nd AK4 jet selection
  const bool pass_jet2 = jet2_sel->passes(event);
  if(!pass_jet2) return false;
  if(debug) cout << "NJetSelection2: ok" << endl;
  // fill_histograms(event, "Jet2");

  // MET selection
  const bool pass_met = met_sel->passes(event);
  if(!pass_met) return false;
  h_cf_raw->Fill(7.5, 1.); h_cf_lumiscaled->Fill(7.5, event.weight);
  if(debug) cout << "METCut: ok" << endl;
  fill_histograms(event, "MET");

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimePreselectionModule)