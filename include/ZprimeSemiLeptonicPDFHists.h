#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h"

#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>



#include <TLorentzVector.h>
#include "TH1F.h"
#include "TH1D.h"
#include <string>
#include <math.h>
#include <sstream>
#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include <glob.h>
#include <cstring>

namespace uhh2 {

  class ZprimeSemiLeptonicPDFHists: public uhh2::Hists {
  public:
    explicit ZprimeSemiLeptonicPDFHists(uhh2::Context&, const std::string&);
    virtual void fill(const uhh2::Event&) override;
    std::string hist_names[100];
    std::string hist_names_tt[100];
    std::string hist_names_dy_d1[100];
    std::string hist_names_dy_d2[100];
    std::string hist_names_sigma_1[100];
    std::string hist_names_sigma_2[100];
    
    std::string hist_names_xi[100];
    std::string hist_names_xi_RM[100];
    
    // Map to store PDF histograms for all binning schemes and f-values
    std::map<std::string, TH1F*> h_pdf_xi_reco_map;

  protected:
    bool is_mc;
    bool is_dy;
    bool is_tt;
    bool is_wjets;
    bool is_qcd_HTbinned;
    bool is_alps;
    bool is_azh;
    bool is_htott_scalar;
    bool is_htott_pseudo;
    bool is_zprimetott;
    bool isMuon;
    bool isElectron;
    bool debug;
    bool ishotvr, isdeepAK8;
  
    uhh2::Event::Handle< std::vector<TopJet> > h_AK8TopTags;
    uhh2::Event::Handle< std::vector<Jet> > h_CHSjets_matched;
    uhh2::Event::Handle<bool> h_is_zprime_reconstructed_chi2;
    uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
    uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;
     
    // NoAC weight support for TTbar
     uhh2::Event::Handle<float> h_xi_gen;
     bool use_noac_evtweights_;
     std::string noac_gen_file_;
     std::string noac_gen_hist_;
     std::map<float, std::unique_ptr<TH1D>> noac_weights_map;
     std::vector<float> f_values;
     
     // Helper functions for NoAC weights
     static std::unique_ptr<TH1D> mirror_hist_1d(const TH1* H);
     static std::unique_ptr<TH1D> build_noac_weights_from_gen(const TH1D& Hgen_in, float f_noac);
     static double lookup_noac_weight(double xi, const TH1* W);
     
    virtual ~ZprimeSemiLeptonicPDFHists();
  };
}