#pragma once

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/NtupleObjects.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>
#include <UHH2/core/include/LorentzVector.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/NSelections.h>
#include "UHH2/common/include/BTagCalibrationStandalone.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/Utils.h"

#include <UHH2/HOTVR/include/HOTVRIds.h>

#include <UHH2/ZprimeSemiLeptonic/include/constants.hpp>

#include <Eigen/Dense>
#include "TH1.h"
#include <TLorentzVector.h>

float inv_mass(const LorentzVector&);

class ZprimeCandidateBuilder : uhh2::AnalysisModule{
public:
  explicit ZprimeCandidateBuilder(uhh2::Context&, TString mode, float minDR = 1.2);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates_;
  uhh2::Event::Handle< std::vector<TopJet> > h_AK8TopTags;
  uhh2::Event::Handle< std::vector<const TopJet*> > h_AK8TopTagsPtr;
  float minDR_;
  TString mode_;
};

class ZprimeChi2Discriminator : uhh2::AnalysisModule{
public:
  explicit ZprimeChi2Discriminator(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates_;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestCandidate_;
  uhh2::Event::Handle<bool> h_is_zprime_reconstructed_;
  float mtoplep_, mtoplep_ttag_;
  float sigmatoplep_, sigmatoplep_ttag_;
  float mtophad_, mtophad_ttag_;
  float sigmatophad_, sigmatophad_ttag_;
};

class ZprimeCorrectMatchDiscriminator : uhh2::AnalysisModule{
public:
  explicit ZprimeCorrectMatchDiscriminator(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates_;
  uhh2::Event::Handle<TTbarGen> h_ttbargen_;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestCandidate_;
  uhh2::Event::Handle<bool> h_is_zprime_reconstructed_;
  bool is_mc;
  std::unique_ptr<TTbarGenProducer> ttgenprod;

};

class AK8PuppiTopTagger : public uhh2::AnalysisModule {
public:
  explicit AK8PuppiTopTagger(uhh2::Context&, int min_num_daughters = 2, float max_dR = 1., float min_mass = 105., float max_mass = 210., float max_tau32 = 0.65);
  virtual bool process(uhh2::Event&) override;

private:
  int min_num_daughters_ = 2;
  float max_dR_;
  float min_mass_;
  float max_mass_;
  float max_tau32_;
  uhh2::Event::Handle< std::vector<TopJet> > h_AK8PuppiTopTags_;
  uhh2::Event::Handle< std::vector<const TopJet*> > h_AK8PuppiTopTagsPtr_;
};

class HOTVRTopTagger : public uhh2::AnalysisModule {
public:
  explicit HOTVRTopTagger(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle< std::vector<TopJet> > h_HOTVRTopTags_;
  uhh2::Event::Handle< std::vector<const TopJet*> > h_HOTVRTopTagsPtr_;
  TopJetId toptag_id = AndId<TopJet>(HOTVRTopTag(0.8, 140.0, 220.0, 50.0), Tau32Groomed(0.56));
};

class DeepAK8TopTagger : public uhh2::AnalysisModule {
public:
  explicit DeepAK8TopTagger(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  Year year;
  uhh2::Event::Handle< std::vector<TopJet> > h_DeepAK8TopTags_;
  uhh2::Event::Handle< std::vector<const TopJet*> > h_DeepAK8TopTagsPtr_;
};

class JetLeptonDeltaRCleaner : public uhh2::AnalysisModule {
public:
  explicit JetLeptonDeltaRCleaner(float mindr=0.8): minDR_(mindr) {}
  virtual bool process(uhh2::Event&) override;

private:
  float minDR_;
};

class TopJetLeptonDeltaRCleaner : public uhh2::AnalysisModule {
public:
  explicit TopJetLeptonDeltaRCleaner(float mindr=0.8): minDR_(mindr) {}
  virtual bool process(uhh2::Event&) override;

private:
  float minDR_;
};

const Particle* leading_lepton(const uhh2::Event&);

float STlep(const uhh2::Event&);
float Muon_pfMINIIso    (const Muon&    , const uhh2::Event&, const std::string&);
float Electron_pfMINIIso(const Electron&, const uhh2::Event&, const std::string&);

bool trigger_bit(const uhh2::Event&, const std::string&);

class GENWToLNuFinder : public uhh2::AnalysisModule {
public:
  explicit GENWToLNuFinder(uhh2::Context&, const std::string&);
  virtual ~GENWToLNuFinder() {}
  virtual bool process(uhh2::Event&) override;

protected:
  uhh2::Event::Handle<GenParticle> h_genWln_W_;
  uhh2::Event::Handle<GenParticle> h_genWln_l_;
  uhh2::Event::Handle<GenParticle> h_genWln_n_;
};

class MEPartonFinder : public uhh2::AnalysisModule {
public:
  explicit MEPartonFinder(uhh2::Context&, const std::string&);
  virtual ~MEPartonFinder() {}
  virtual bool process(uhh2::Event&) override;

protected:
  uhh2::Event::Handle<std::vector<GenParticle> > h_meps_;
};

/////////////////////////////////////////////////////////////////////////
////////////////////////////     NN vars     ////////////////////////////
//////////////////////////// i.e. before DNN ////////////////////////////
/////////////////////////////////////////////////////////////////////////
class Variables_NN : uhh2::AnalysisModule{
public:
  explicit Variables_NN(uhh2::Context&, TString mode);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle<std::vector<Jet> > h_CHSjets_matched;
  uhh2::Event::Handle<float> h_Mu_pt, h_Mu_eta, h_Mu_phi, h_Mu_E;
  uhh2::Event::Handle<float> h_Ele_pt, h_Ele_eta, h_Ele_phi, h_Ele_E;
  uhh2::Event::Handle<float> h_MET_pt, h_MET_phi;
  uhh2::Event::Handle<float> h_N_HOTVR;
  uhh2::Event::Handle<float> h_HOTVR_j1_pt, h_HOTVR_j1_eta, h_HOTVR_j1_phi, h_HOTVR_j1_E, h_HOTVR_j1_mSD, h_HOTVR_j1_tau21, h_HOTVR_j1_tau32;
  uhh2::Event::Handle<float> h_HOTVR_j2_pt, h_HOTVR_j2_eta, h_HOTVR_j2_phi, h_HOTVR_j2_E, h_HOTVR_j2_mSD, h_HOTVR_j2_tau21, h_HOTVR_j2_tau32;
  uhh2::Event::Handle<float> h_HOTVR_j3_pt, h_HOTVR_j3_eta, h_HOTVR_j3_phi, h_HOTVR_j3_E, h_HOTVR_j3_mSD, h_HOTVR_j3_tau21, h_HOTVR_j3_tau32;
  uhh2::Event::Handle<float> h_N_Ak8, h_N_Ak4;
  uhh2::Event::Handle<float> h_Ak8_j1_pt, h_Ak8_j1_eta, h_Ak8_j1_phi, h_Ak8_j1_E, h_Ak8_j1_mSD, h_Ak8_j1_tau21, h_Ak8_j1_tau32, h_Ak8_j1_deepak8tscore;
  uhh2::Event::Handle<float> h_Ak8_j2_pt, h_Ak8_j2_eta, h_Ak8_j2_phi, h_Ak8_j2_E, h_Ak8_j2_mSD, h_Ak8_j2_tau21, h_Ak8_j2_tau32, h_Ak8_j2_deepak8tscore;
  uhh2::Event::Handle<float> h_Ak8_j3_pt, h_Ak8_j3_eta, h_Ak8_j3_phi, h_Ak8_j3_E, h_Ak8_j3_mSD, h_Ak8_j3_tau21, h_Ak8_j3_tau32, h_Ak8_j3_deepak8tscore;
  uhh2::Event::Handle<float> h_Ak4_j1_pt, h_Ak4_j1_eta, h_Ak4_j1_phi, h_Ak4_j1_E, h_Ak4_j1_m, h_Ak4_j1_deepjetbscore;
  uhh2::Event::Handle<float> h_Ak4_j2_pt, h_Ak4_j2_eta, h_Ak4_j2_phi, h_Ak4_j2_E, h_Ak4_j2_m, h_Ak4_j2_deepjetbscore;
  uhh2::Event::Handle<float> h_Ak4_j3_pt, h_Ak4_j3_eta, h_Ak4_j3_phi, h_Ak4_j3_E, h_Ak4_j3_m, h_Ak4_j3_deepjetbscore;
  uhh2::Event::Handle<float> h_Ak4_j4_pt, h_Ak4_j4_eta, h_Ak4_j4_phi, h_Ak4_j4_E, h_Ak4_j4_m, h_Ak4_j4_deepjetbscore;
  uhh2::Event::Handle<float> h_Ak4_j5_pt, h_Ak4_j5_eta, h_Ak4_j5_phi, h_Ak4_j5_E, h_Ak4_j5_m, h_Ak4_j5_deepjetbscore;
  uhh2::Event::Handle<float> h_Ak4_j6_pt, h_Ak4_j6_eta, h_Ak4_j6_phi, h_Ak4_j6_E, h_Ak4_j6_m, h_Ak4_j6_deepjetbscore;
  uhh2::Event::Handle<float> h_uniform_random;
  TString mode_;
};



/////////////////////////////////////////////////////////////////////////
//////////////////     Spin Correlation variabless     //////////////////
/////////////////////////////////////////////////////////////////////////
class SpinCorrelations : uhh2::AnalysisModule{
public:
  explicit SpinCorrelations(uhh2::Context&, TString mode);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle<bool> h_is_zprime_reconstructed_chi2;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
  uhh2::Event::Handle<std::vector<Jet> > h_CHSjets_matched;
  uhh2::Event::Handle<float> h_eventweight;
  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;
  uhh2::Event::Handle<float> h_chi2, h_M_tt, h_beta, h_dyreco;
  uhh2::Event::Handle<float> h_cosTheta1k_antiLep, h_cosTheta1r_antiLep, h_cosTheta1n_antiLep, h_cosTheta1kStar_antiLep, h_cosTheta1rStar_antiLep;
  uhh2::Event::Handle<float> h_cosTheta2k_Lep,     h_cosTheta2r_Lep,     h_cosTheta2n_Lep,     h_cosTheta2kStar_Lep,     h_cosTheta2rStar_Lep;
  uhh2::Event::Handle<float> h_cosTheta1k, h_cosTheta1r, h_cosTheta1n, h_cosTheta1kStar, h_cosTheta1rStar;
  uhh2::Event::Handle<float> h_cosTheta2k, h_cosTheta2r, h_cosTheta2n, h_cosTheta2kStar, h_cosTheta2rStar;
  uhh2::Event::Handle<float> h_Cnn, h_Cnr, h_Cnk, h_Crn, h_Crr, h_Crk, h_Ckn, h_Ckr, h_Ckk;
  uhh2::Event::Handle<float> h_Crk_plus, h_Crk_minus, h_Cnr_plus, h_Cnr_minus, h_Cnk_plus, h_Cnk_minus;
  uhh2::Event::Handle<float> h_cHel,     h_cHel_Mtt300_400,     h_cHel_Mtt300_400_betaLT0p9;
  uhh2::Event::Handle<float> h_cHel_P3n, h_cHel_P3n_Mtt800_Inf, h_cHel_P3n_Mtt800_Inf_cosThetaLT0p4;
  uhh2::Event::Handle<float> h_Sigma_phi, h_Delta_phi;
  uhh2::Event::Handle<float> h_Sigma_phi_1, h_Sigma_phi_2, h_Delta_phi_1, h_Delta_phi_2;
  uhh2::Event::Handle<float> h_dyreco_s1, h_dyreco_s2, h_dyreco_d1, h_dyreco_d2;
  uhh2::Event::Handle<float> h_Sigma_phi_high, h_Sigma_phi_low, h_Delta_phi_high, h_Delta_phi_low, h_dyreco_high, h_dyreco_low;
  TString mode_;
};


/////////////////////////////////////////////////////////////////////////
////////////////////////////     EFT SR     /////////////////////////////
/////////////////////////////////////////////////////////////////////////
class Variables_EFT_SR : uhh2::AnalysisModule{
  public:
    explicit Variables_EFT_SR(uhh2::Context&, TString mode);
    virtual bool process(uhh2::Event&) override;
  
  private:
  uhh2::Event::Handle<bool> h_is_zprime_reconstructed_chi2;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
  uhh2::Event::Handle<std::vector<Jet> > h_CHSjets_matched;
  uhh2::Event::Handle<float> h_eventweight_SR;
  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;
  uhh2::Event::Handle<float> h_chi2_SR, h_M_tt_SR, h_beta_SR, h_dyreco_SR;
  uhh2::Event::Handle<float> h_cosTheta1k_antiLep_SR, h_cosTheta1r_antiLep_SR, h_cosTheta1n_antiLep_SR, h_cosTheta1kStar_antiLep_SR, h_cosTheta1rStar_antiLep_SR;
  uhh2::Event::Handle<float> h_cosTheta2k_Lep_SR, h_cosTheta2r_Lep_SR, h_cosTheta2n_Lep_SR, h_cosTheta2kStar_Lep_SR, h_cosTheta2rStar_Lep_SR;
  uhh2::Event::Handle<float> h_cosTheta1k_SR, h_cosTheta1r_SR, h_cosTheta1n_SR, h_cosTheta1kStar_SR, h_cosTheta1rStar_SR;
  uhh2::Event::Handle<float> h_cosTheta2k_SR, h_cosTheta2r_SR, h_cosTheta2n_SR, h_cosTheta2kStar_SR, h_cosTheta2rStar_SR;
  uhh2::Event::Handle<float> h_Cnn_SR, h_Cnr_SR, h_Cnk_SR, h_Crn_SR, h_Crr_SR, h_Crk_SR, h_Ckn_SR, h_Ckr_SR, h_Ckk_SR;
  uhh2::Event::Handle<float> h_Crk_plus_SR, h_Crk_minus_SR, h_Cnr_plus_SR, h_Cnr_minus_SR, h_Cnk_plus_SR, h_Cnk_minus_SR;
  uhh2::Event::Handle<float> h_cHel_SR, h_cHel_Mtt300_400_SR, h_cHel_Mtt300_400_betaLT0p9_SR;
  uhh2::Event::Handle<float> h_cHel_P3n_SR, h_cHel_P3n_Mtt800_Inf_SR, h_cHel_P3n_Mtt800_Inf_cosThetaLT0p4_SR;
  uhh2::Event::Handle<float> h_Sigma_phi_SR, h_Delta_phi_SR;
  uhh2::Event::Handle<float> h_Sigma_phi_1_SR, h_Sigma_phi_2_SR, h_Delta_phi_1_SR, h_Delta_phi_2_SR;
  uhh2::Event::Handle<float> h_dyreco_s1_SR, h_dyreco_s2_SR, h_dyreco_d1_SR, h_dyreco_d2_SR;
  uhh2::Event::Handle<float> h_Sigma_phi_high_SR, h_Sigma_phi_low_SR, h_Delta_phi_high_SR, h_Delta_phi_low_SR, h_dyreco_high_SR, h_dyreco_low_SR;
  TString mode_;
};


//////////////////////////////////////////////////////////////////////////
////////////////////////////     EFT CR1     /////////////////////////////
//////////////////////////////////////////////////////////////////////////
class Variables_EFT_CR1 : uhh2::AnalysisModule{
  public:
    explicit Variables_EFT_CR1(uhh2::Context&, TString mode);
    virtual bool process(uhh2::Event&) override;
  
  private:
  uhh2::Event::Handle<bool> h_is_zprime_reconstructed_chi2;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
  uhh2::Event::Handle<std::vector<Jet> > h_CHSjets_matched;
  uhh2::Event::Handle<float> h_eventweight_CR1;
  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;
  uhh2::Event::Handle<float> h_chi2_CR1, h_M_tt_CR1, h_beta_CR1, h_dyreco_CR1;
  uhh2::Event::Handle<float> h_cosTheta1k_antiLep_CR1, h_cosTheta1r_antiLep_CR1, h_cosTheta1n_antiLep_CR1, h_cosTheta1kStar_antiLep_CR1, h_cosTheta1rStar_antiLep_CR1;
  uhh2::Event::Handle<float> h_cosTheta2k_Lep_CR1, h_cosTheta2r_Lep_CR1, h_cosTheta2n_Lep_CR1, h_cosTheta2kStar_Lep_CR1, h_cosTheta2rStar_Lep_CR1;
  uhh2::Event::Handle<float> h_cosTheta1k_CR1, h_cosTheta1r_CR1, h_cosTheta1n_CR1, h_cosTheta1kStar_CR1, h_cosTheta1rStar_CR1;
  uhh2::Event::Handle<float> h_cosTheta2k_CR1, h_cosTheta2r_CR1, h_cosTheta2n_CR1, h_cosTheta2kStar_CR1, h_cosTheta2rStar_CR1;
  uhh2::Event::Handle<float> h_Cnn_CR1, h_Cnr_CR1, h_Cnk_CR1, h_Crn_CR1, h_Crr_CR1, h_Crk_CR1, h_Ckn_CR1, h_Ckr_CR1, h_Ckk_CR1;
  uhh2::Event::Handle<float> h_Crk_plus_CR1, h_Crk_minus_CR1, h_Cnr_plus_CR1, h_Cnr_minus_CR1, h_Cnk_plus_CR1, h_Cnk_minus_CR1;
  uhh2::Event::Handle<float> h_cHel_CR1, h_cHel_Mtt300_400_CR1, h_cHel_Mtt300_400_betaLT0p9_CR1;
  uhh2::Event::Handle<float> h_cHel_P3n_CR1, h_cHel_P3n_Mtt800_Inf_CR1, h_cHel_P3n_Mtt800_Inf_cosThetaLT0p4_CR1;
  uhh2::Event::Handle<float> h_Sigma_phi_CR1, h_Delta_phi_CR1;
  uhh2::Event::Handle<float> h_Sigma_phi_1_CR1, h_Sigma_phi_2_CR1, h_Delta_phi_1_CR1, h_Delta_phi_2_CR1;
  uhh2::Event::Handle<float> h_dyreco_s1_CR1, h_dyreco_s2_CR1, h_dyreco_d1_CR1, h_dyreco_d2_CR1;
  uhh2::Event::Handle<float> h_Sigma_phi_high_CR1, h_Sigma_phi_low_CR1, h_Delta_phi_high_CR1, h_Delta_phi_low_CR1, h_dyreco_high_CR1, h_dyreco_low_CR1;
  TString mode_;
};


//////////////////////////////////////////////////////////////////////////
////////////////////////////     EFT CR2     /////////////////////////////
//////////////////////////////////////////////////////////////////////////
class Variables_EFT_CR2 : uhh2::AnalysisModule{
  public:
    explicit Variables_EFT_CR2(uhh2::Context&, TString mode);
    virtual bool process(uhh2::Event&) override;
  
  private:
  uhh2::Event::Handle<bool> h_is_zprime_reconstructed_chi2;
  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;
  uhh2::Event::Handle<std::vector<Jet> > h_CHSjets_matched;
  uhh2::Event::Handle<float> h_eventweight_CR2;
  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;
  uhh2::Event::Handle<float> h_chi2_CR2, h_M_tt_CR2, h_beta_CR2, h_dyreco_CR2;
  uhh2::Event::Handle<float> h_cosTheta1k_antiLep_CR2, h_cosTheta1r_antiLep_CR2, h_cosTheta1n_antiLep_CR2, h_cosTheta1kStar_antiLep_CR2, h_cosTheta1rStar_antiLep_CR2;
  uhh2::Event::Handle<float> h_cosTheta2k_Lep_CR2, h_cosTheta2r_Lep_CR2, h_cosTheta2n_Lep_CR2, h_cosTheta2kStar_Lep_CR2, h_cosTheta2rStar_Lep_CR2;
  uhh2::Event::Handle<float> h_cosTheta1k_CR2, h_cosTheta1r_CR2, h_cosTheta1n_CR2, h_cosTheta1kStar_CR2, h_cosTheta1rStar_CR2;
  uhh2::Event::Handle<float> h_cosTheta2k_CR2, h_cosTheta2r_CR2, h_cosTheta2n_CR2, h_cosTheta2kStar_CR2, h_cosTheta2rStar_CR2;
  uhh2::Event::Handle<float> h_Cnn_CR2, h_Cnr_CR2, h_Cnk_CR2, h_Crn_CR2, h_Crr_CR2, h_Crk_CR2, h_Ckn_CR2, h_Ckr_CR2, h_Ckk_CR2;
  uhh2::Event::Handle<float> h_Crk_plus_CR2, h_Crk_minus_CR2, h_Cnr_plus_CR2, h_Cnr_minus_CR2, h_Cnk_plus_CR2, h_Cnk_minus_CR2;
  uhh2::Event::Handle<float> h_cHel_CR2, h_cHel_Mtt300_400_CR2, h_cHel_Mtt300_400_betaLT0p9_CR2;
  uhh2::Event::Handle<float> h_cHel_P3n_CR2, h_cHel_P3n_Mtt800_Inf_CR2, h_cHel_P3n_Mtt800_Inf_cosThetaLT0p4_CR2;
  uhh2::Event::Handle<float> h_Sigma_phi_CR2, h_Delta_phi_CR2;
  uhh2::Event::Handle<float> h_Sigma_phi_1_CR2, h_Sigma_phi_2_CR2, h_Delta_phi_1_CR2, h_Delta_phi_2_CR2;
  uhh2::Event::Handle<float> h_dyreco_s1_CR2, h_dyreco_s2_CR2, h_dyreco_d1_CR2, h_dyreco_d2_CR2;
  uhh2::Event::Handle<float> h_Sigma_phi_high_CR2, h_Sigma_phi_low_CR2, h_Delta_phi_high_CR2, h_Delta_phi_low_CR2, h_dyreco_high_CR2, h_dyreco_low_CR2;
  TString mode_;
};

// Generic Class for Applying SFs - from Andrea
class ScaleFactorsFromHistos : public uhh2::AnalysisModule {
public:
  void LoadHisto(TFile* file, std::string name, std::string hname);
  double Evaluator(std::string hname, double var);

protected:
  std::unordered_map<std::string, std::unique_ptr<TH1F> > histos;
};

// Apply Corrections V+jets
class NLOCorrections : public ScaleFactorsFromHistos {
public:
  explicit NLOCorrections(uhh2::Context& ctx);
  virtual bool process(uhh2::Event&) override;
  double GetPartonObjectPt(uhh2::Event& event, ParticleID objID);

private:
  bool is_Wjets, is_Zjets, is_DY, is_Znn, is2016;
};


// Top pT Reweighting extended - from Alex F.
class TopPtReweighting : public uhh2::AnalysisModule {
 public:
  explicit TopPtReweighting(uhh2::Context& ctx,
    float a, float b,
    const std::string& syst_a,
    const std::string& syst_b,
    const std::string& ttgen_name ="");
  virtual bool process(uhh2::Event& event) override;

 private:
    uhh2::Event::Handle<TTbarGen> h_ttbargen_;
    float a_, b_;
    std::string version_;
    std::string ttgen_name_;
    uhh2::Event::Handle< float > h_weight_toppt_nominal;
    uhh2::Event::Handle< float > h_weight_toppt_a_up;
    uhh2::Event::Handle< float > h_weight_toppt_b_up;
    uhh2::Event::Handle< float > h_weight_toppt_a_down;
    uhh2::Event::Handle< float > h_weight_toppt_b_down;
};


class PuppiCHS_matching : public uhh2::AnalysisModule {
public:
  explicit PuppiCHS_matching(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle< std::vector<Jet> > h_CHSjets;
  uhh2::Event::Handle< std::vector<Jet> > h_CHS_matched_;
};

// Muon RECO SF (from https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2016#RECO_efficiency)
// Add handles with nominal values and up/down variations and apply nominal sf to event.weight

class MuonRecoSF : public uhh2::AnalysisModule {
public:
  explicit MuonRecoSF(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

private:
  Year year;
  bool is_mc, is_Muon;
  uhh2::Event::Handle< float > h_muonrecSF_nominal;
  uhh2::Event::Handle< float > h_muonrecSF_up;
  uhh2::Event::Handle< float > h_muonrecSF_down;
};

// Structure Constants Calculator for EFT weights
class StructureConstantsCalculator : public uhh2::AnalysisModule {
public:
  explicit StructureConstantsCalculator(uhh2::Context& ctx, int num_WCs = 16);
  virtual bool process(uhh2::Event& event) override;
  
  // Helper methods for structure constants calculation
  std::vector<std::vector<float>> generate_weight_configurations(int num_WCs);
  std::vector<float> obtain_structure_constant(int num_WCs, const std::vector<float>& mg_weights);
  std::vector<float> calculate_new_weights(const std::vector<float>& structs, const std::vector<float>& wc_values);
  
  // Get the correct weight index for a given configuration
  int get_weight_index(const std::vector<float>& config);
  
  // Helper methods for WC name/index translation
  int get_wc_index(const std::string& wc_name) const;
  std::string get_wc_name(int index) const;
  
private:
  int num_WCs_;
  uhh2::Event::Handle<std::vector<float>> h_structure_constants_;
  
  // Wilson Coefficient names for better readability and debugging
  std::vector<std::string> wc_names_;
  
  // Weight index mappings
  std::map<int, int> single2_mapping_;
  std::map<std::pair<int, int>, int> known_pairs_;
};