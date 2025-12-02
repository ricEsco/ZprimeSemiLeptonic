#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/core/include/LorentzVector.h>
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>
#include "UHH2/core/include/Event.h"
#include <UHH2/common/include/NSelections.h>
#include "UHH2/common/include/BTagCalibrationStandalone.h"
#include "UHH2/common/include/MCWeight.h"

#include "boost/algorithm/string.hpp"

#include "Riostream.h"
#include "TFile.h"
#include "TH1F.h"

using namespace std;
using namespace uhh2;




float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

FlavorParticle find_primary_lepton(uhh2::Event & event) {
  assert(event.muons || event.electrons);
  FlavorParticle primlep;
  float ptmax = 0.;
  if(event.electrons) {
    for(const auto & ele : *event.electrons) {
      if(ele.pt() > ptmax) {
        ptmax = ele.pt();
        primlep = ele;
      }
    }
  }
  if(event.muons) {
    for(const auto & mu : *event.muons) {
      if(mu.pt() > ptmax) {
        ptmax = mu.pt();
        primlep = mu;
      }
    }
  }
  return primlep;
}

vector<LorentzVector> reconstruct_neutrino(const LorentzVector & lepton, const LorentzVector & met) {
  TVector3 lepton_pT = toVector(lepton);
  lepton_pT.SetZ(0);
  TVector3 neutrino_pT = toVector(met);
  neutrino_pT.SetZ(0);
  constexpr float mass_w = 80.399f;
  float mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT;
  float A = - (lepton_pT * lepton_pT);
  float B = mu * lepton.pz();
  float C = mu * mu - lepton.e() * lepton.e() * (neutrino_pT * neutrino_pT);
  float discriminant = B * B - A * C;
  std::vector<LorentzVector> solutions;
  if (0 >= discriminant) {
    // Take only real part of the solution for pz:
    LorentzVectorXYZE solution (met.Px(),met.Py(),-B / A,0);
    solution.SetE(solution.P());
    solutions.emplace_back(toPtEtaPhi(solution));
  }
  else {
    discriminant = sqrt(discriminant);
    LorentzVectorXYZE solution (met.Px(),met.Py(),(-B - discriminant) / A,0);
    solution.SetE(solution.P());
    solutions.emplace_back(toPtEtaPhi(solution));

    LorentzVectorXYZE solution2 (met.Px(),met.Py(),(-B + discriminant) / A,0);
    solution2.SetE(solution2.P());
    solutions.emplace_back(toPtEtaPhi(solution2));
  }
  return solutions;
}


ZprimeCandidateBuilder::ZprimeCandidateBuilder(uhh2::Context& ctx, TString mode, float minDR) : minDR_(minDR), mode_(mode){

  h_ZprimeCandidates_ = ctx.get_handle< vector<ZprimeCandidate> >("ZprimeCandidates");

  if(mode_ == "deepAK8"){
    h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("DeepAK8TopTags");
    h_AK8TopTagsPtr = ctx.get_handle<std::vector<const TopJet*>>("DeepAK8TopTagsPtr");
  }else if(mode_ == "hotvr"){
    h_AK8TopTags = ctx.get_handle<std::vector<TopJet>>("HOTVRTopTags");
    h_AK8TopTagsPtr = ctx.get_handle<std::vector<const TopJet*>>("HOTVRTopTagsPtr");
  }

  if(mode_ != "hotvr" && mode_ != "deepAK8") throw runtime_error("In ZprimeCandidateBuilder::ZprimeCandidateBuilder(): 'mode' must be 'hotvr' or 'deepAK8'");

}

bool ZprimeCandidateBuilder::process(uhh2::Event& event){
  assert(event.jets);
  assert(event.muons || event.electrons);
  assert(event.met);

  // Declare output
  vector<ZprimeCandidate> candidates;

  // Find primary lepton
  FlavorParticle lepton = find_primary_lepton(event);

  // Reconstruct neutrino
  vector<LorentzVector> neutrinos = reconstruct_neutrino(lepton.v4(), event.met->v4());
  unsigned int neutrinoidx = 0;

  // Build all necessary loops
  vector<TopJet> TopTags = event.get(h_AK8TopTags);
  vector<const TopJet*> TopTagsPtr = event.get(h_AK8TopTagsPtr);

  if((event.muons->size() < 1 && event.electrons->size() < 1)) throw runtime_error("Event content did not allow reconstructing the Zprime: Leptons");
  if((event.jets->size() < 2 && TopTags.size() == 0)) throw runtime_error("Event content did not allow reconstructing the Zprime: AK4");
  if((event.jets->size() < 1 && TopTags.size() >= 1)) throw runtime_error("Event content did not allow reconstructing the Zprime: Top-tag");

  // Must have at least one AK4 jet with dR > 1.2
  vector<bool> has_separated_jet;
  for(unsigned int i=0; i<TopTags.size(); i++){
    bool is_sep = false;
    for(unsigned int k = 0; k < event.jets->size(); k++){
      if(deltaR(event.jets->at(k), TopTags[i]) > 1.2) is_sep = true;
    }
    has_separated_jet.emplace_back(is_sep);
  }

  vector<bool> overlap_with_lepton;
  double maxDeltaR = -1;
  if(mode_ == "deepAK8"){
    for(const TopJet & toptag : TopTags){
      bool overlap = true;
      if(deltaR(lepton, toptag) > 0.8) overlap = false;
      overlap_with_lepton.emplace_back(overlap);
    }
  } else if(mode_ == "hotvr"){
    for(const TopJet & toptag : TopTags){
      if(deltaR(lepton, toptag)>maxDeltaR) maxDeltaR = deltaR(lepton, toptag);
      bool overlap = true;
      if(deltaR(lepton, toptag) > 1.5) overlap = false;
      overlap_with_lepton.emplace_back(overlap);
    }
  }

  bool do_toptag_reco = false;
  if(TopTags.size() >= 1){
    for(unsigned int i=0; i<TopTags.size(); i++){
      //cout << "Number of toptags: " << TopTags.size() << endl;
      //cout << "At " << i << " of " << has_separated_jet.size() << " and " << overlap_with_lepton.size() << endl;
      if(has_separated_jet[i] &&  !overlap_with_lepton[i]) do_toptag_reco = true;
    }
  }

  if(!do_toptag_reco){ // AK4 reconstruction

    unsigned int njets = event.jets->size();
    if(njets > 10) njets = 10;

    unsigned int jetiters = pow(3, njets);
    for(const auto & neutrino_v4 : neutrinos) {
      for (unsigned int j=0; j < jetiters; j++) {

        LorentzVector tophadv4;
        vector<Particle> tophadjets;
        LorentzVector toplepv4 = lepton.v4() + neutrino_v4;
        vector<Particle> toplepjets;
        int num = j;
        ZprimeCandidate candidate;
        candidate.set_is_toptag_reconstruction(false);
        candidate.set_is_puppi_reconstruction(false);

        for (unsigned int k=0; k<njets; k++) {
          if(num%3==0){
            tophadv4 = tophadv4 + event.jets->at(k).v4();
            tophadjets.emplace_back(event.jets->at(k));
          }
          if(num%3==1){
            toplepv4 = toplepv4 + event.jets->at(k).v4();
            toplepjets.emplace_back(event.jets->at(k));
          }
          num /= 3;
        }

        if(tophadjets.size() < 1 || toplepjets.size() < 1) continue;

        // Set all member variables of the candidate
        candidate.set_Zprime_v4(tophadv4 + toplepv4);
        candidate.set_top_hadronic_v4(tophadv4);
        candidate.set_top_leptonic_v4(toplepv4);
        candidate.set_jets_hadronic(tophadjets);
        candidate.set_jets_leptonic(toplepjets);
        candidate.set_lepton(lepton);
        candidate.set_neutrino_v4(neutrino_v4);
        candidate.set_neutrinoindex(neutrinoidx);
        candidates.emplace_back(move(candidate));
      }
      neutrinoidx++;
    }
  }
  else{ // TopTag reconstruction
    for(const auto & neutrino_v4 : neutrinos) {
      for (unsigned int j=0; j < TopTags.size(); j++) {
        if(!(has_separated_jet[j] &&  !overlap_with_lepton[j])) continue;

        TopJet toptag = TopTags.at(j);
        const TopJet* toptag_ptr = TopTagsPtr.at(j);

        if(mode_ == "hotvr"){
          // Only HOTVR jet farest from lepton
          if(deltaR(lepton,toptag)<maxDeltaR) continue;
        }

        // Only consider well-separated AK4 jets
        vector<Jet> separated_jets;
        for(unsigned int k = 0; k < event.jets->size(); k++){
          if(deltaR(event.jets->at(k), toptag) > minDR_) separated_jets.emplace_back(event.jets->at(k));
        }
        unsigned int njets = separated_jets.size();
        if(njets < 1) throw runtime_error("In TopTagReco (PUPPI): This toptag does not have >= 1 well-separated AK4 jet. This should have been caught by earlier messages. There is a logic error.");
        if(njets > 10) njets = 10;

        unsigned int jetiters = pow(2, njets);
        for (unsigned int k=0; k < jetiters; k++) {

          LorentzVector tophadv4 = toptag.v4();
          vector<Particle> tophadjets;
          tophadjets.emplace_back(toptag);
          LorentzVector toplepv4 = lepton.v4() + neutrino_v4;
          vector<Particle> toplepjets;
          int num = k;
          ZprimeCandidate candidate;
          candidate.set_is_toptag_reconstruction(true);
          candidate.set_is_puppi_reconstruction(true);

          for (unsigned int l=0; l<njets; l++) {
            if(num%2==0){
              toplepv4 = toplepv4 + separated_jets.at(l).v4();
              toplepjets.emplace_back(separated_jets.at(l));
            }
            num /= 2;
          }

          if(tophadjets.size() < 1 || toplepjets.size() < 1 || separated_jets.size() < 1) continue;

          // Set all member variables of the candidate
          candidate.set_Zprime_v4(tophadv4 + toplepv4);
          candidate.set_top_hadronic_v4(tophadv4);
          candidate.set_top_leptonic_v4(toplepv4);
          candidate.set_jets_hadronic(tophadjets);
          candidate.set_jets_leptonic(toplepjets);
          candidate.set_tophad_topjet_ptr(toptag_ptr);
          candidate.set_lepton(lepton);
          candidate.set_neutrinoindex(neutrinoidx);
          candidate.set_neutrino_v4(neutrino_v4);
          candidates.emplace_back(candidate);

        }
      }
      neutrinoidx++;
    }
  }
  // Set the handle with all candidates
  event.set(h_ZprimeCandidates_, candidates);

  return true;
}


ZprimeChi2Discriminator::ZprimeChi2Discriminator(uhh2::Context& ctx){

  h_ZprimeCandidates_ = ctx.get_handle< std::vector<ZprimeCandidate> >("ZprimeCandidates");
  h_is_zprime_reconstructed_ = ctx.get_handle< bool >("is_zprime_reconstructed_chi2");
  h_BestCandidate_ = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");

  mtoplep_ = 173.6;
  sigmatoplep_ = 24.6;
  mtophad_ = 173.0;
  sigmatophad_ = 21.2;

  mtoplep_ttag_ = 171.4;
  sigmatoplep_ttag_ = 22.0;
  mtophad_ttag_ = 180.6;
  sigmatophad_ttag_ = 15.6;

}

bool ZprimeChi2Discriminator::process(uhh2::Event& event){

  vector<ZprimeCandidate>& candidates = event.get(h_ZprimeCandidates_);
  if(candidates.size() < 1) return false;

  float chi2_best = 99999999;
  ZprimeCandidate* bestCand = &candidates.at(0);
  for(unsigned int i=0; i<candidates.size(); i++){
    bool is_toptag_reconstruction = candidates.at(i).is_toptag_reconstruction();

    float chi2_had = 0.;
    float chi2_lep = 0.;
    float mhad = 0.;
    float mlep = 0.;
    if(is_toptag_reconstruction){

      if(!candidates.at(i).is_puppi_reconstruction()) mhad = candidates.at(i).tophad_topjet_ptr()->softdropmass();
      else{
        LorentzVector SumSubjets(0.,0.,0.,0.);
        for(unsigned int k=0; k<candidates.at(i).tophad_topjet_ptr()->subjets().size(); k++) SumSubjets = SumSubjets + candidates.at(i).tophad_topjet_ptr()->subjets().at(k).v4();
        mhad = inv_mass(SumSubjets);
      }
      mlep = inv_mass(candidates.at(i).top_leptonic_v4());
      chi2_had = pow((mhad - mtophad_ttag_) / sigmatophad_ttag_,2);
      chi2_lep = pow((mlep - mtoplep_ttag_) / sigmatoplep_ttag_,2);
    }
    else{
      mhad = inv_mass(candidates.at(i).top_hadronic_v4());
      mlep = inv_mass(candidates.at(i).top_leptonic_v4());
      chi2_had = pow((mhad - mtophad_) / sigmatophad_,2);
      chi2_lep = pow((mlep - mtoplep_) / sigmatoplep_,2);
    }

    float chi2 = chi2_had + chi2_lep;

    candidates.at(i).set_discriminators("chi2_hadronic", chi2_had);
    candidates.at(i).set_discriminators("chi2_leptonic", chi2_lep);
    candidates.at(i).set_discriminators("chi2_total", chi2);

    if(chi2 < chi2_best){
      chi2_best = chi2;
      bestCand = &candidates.at(i);
    }
  }
  event.set(h_BestCandidate_, bestCand);
  event.set(h_is_zprime_reconstructed_, true);
  return true;
}

// match particle p to one of the jets (Delta R < 0.3); return the deltaR
// of the match.
template<typename T> // T should inherit from Particle
float match_dr(const Particle & p, const std::vector<T> & jets, int& index){
  float mindr = 999999;
  index = -1;
  for(unsigned int i=0; i<jets.size(); ++i){
    float dR = deltaR(p, jets.at(i));
    if( dR <0.4 && dR<mindr) {
      mindr=dR;
      index=i;
    }
  }
  return mindr;
}

ZprimeCorrectMatchDiscriminator::ZprimeCorrectMatchDiscriminator(uhh2::Context& ctx){

  h_ZprimeCandidates_ = ctx.get_handle< std::vector<ZprimeCandidate> >("ZprimeCandidates");
  h_ttbargen_ = ctx.get_handle<TTbarGen>("ttbargen");
  h_is_zprime_reconstructed_ = ctx.get_handle< bool >("is_zprime_reconstructed_correctmatch");
  h_BestCandidate_ = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestCorrectMatch");

  is_mc = ctx.get("dataset_type") == "MC";
  if(is_mc) ttgenprod.reset(new TTbarGenProducer(ctx));
}

bool ZprimeCorrectMatchDiscriminator::process(uhh2::Event& event){

  if(!is_mc) return false;

  // Check if event contains == 2 top quarks
  assert(event.genparticles);
  int n_top = 0, n_antitop = 0;
  for(const auto & gp : *event.genparticles){
    if(gp.pdgId() == 6) n_top++;
    else if(gp.pdgId() == -6) n_antitop++;
  }
  if(n_top != 1 || n_antitop != 1) return false;
  // bool check_decay = ttgenprod->process(event);
  //if(!check_decay) return false; //FixME: sometimes decay prodcts of ttbar are not Wb+Wb. Why?

  vector<ZprimeCandidate>& candidates = event.get(h_ZprimeCandidates_);
  if(candidates.size() < 1) return false;

  float dr_best = 99999999; // This number should be larger than the largest possible value from above, just so there is always a 'bestCand' set
  ZprimeCandidate* bestCand = &candidates.at(0);
  TTbarGen ttbargen = event.get(h_ttbargen_);
  for(unsigned int i=0; i<candidates.size(); i++){

    bool is_toptag_reconstruction = candidates.at(i).is_toptag_reconstruction();

    // Gen-Lvl ttbar has to decay semileptonically
    if(ttbargen.DecayChannel() != TTbarGen::e_muhad && ttbargen.DecayChannel() != TTbarGen::e_ehad){
      candidates.at(i).set_discriminators("correct_match", 9999999);
      // cout << "Not semileptonic decay" << endl;
      continue;
    }

    vector<Particle> jets_had = candidates.at(i).jets_hadronic();
    vector<Particle> jets_lep = candidates.at(i).jets_leptonic();

    if(jets_lep.size() != 1){
      candidates.at(i).set_discriminators("correct_match", 9999999);
      // cout << "Not ==1 leptonic jet" << endl;
      continue;
    }

    if((!is_toptag_reconstruction && jets_had.size() > 3) || (is_toptag_reconstruction && jets_had.size() != 1)){
      candidates.at(i).set_discriminators("correct_match", 9999999);
      // cout << "Not right amount of hadronic jets" << endl;
      continue;
    }
    float correct_dr = 0.;
    int idx;
    float dr;

    // Match leptonic b-quark
    dr = match_dr(ttbargen.BLep(), jets_lep, idx);
    if(dr > 0.4){
      candidates.at(i).set_discriminators("correct_match", 9999999);
      // cout << "Not leptonic b-quark matched" << endl;
      continue;
    }
    correct_dr += dr;

    if(!is_toptag_reconstruction){

      unsigned int n_matched = 0;

      // Match hadronic b-quark
      dr = match_dr(ttbargen.BHad(), jets_had, idx);
      if(dr > 0.4){
        candidates.at(i).set_discriminators("correct_match", 9999999);
        // cout << "Not hadronic b-quark matched (AK4)" << endl;
        continue;
      }
      correct_dr += dr;
      if(idx >= 0) n_matched++;

      //match quarks from W decays to jets
      // First
      dr = match_dr(ttbargen.Q1(), jets_had, idx);
      if(dr > 0.4){
        candidates.at(i).set_discriminators("correct_match", 9999999);
        // cout << "Not had Q1 matched (AK4)" << endl;
        continue;
      }
      correct_dr += dr;
      if(idx >= 0) n_matched++;

      // Second
      dr = match_dr(ttbargen.Q2(), jets_had, idx);
      if(dr > 0.4){
        candidates.at(i).set_discriminators("correct_match", 9999999);
        // cout << "Not had Q2 matched (AK4)" << endl;
        continue;
      }
      correct_dr += dr;
      if(idx >= 0) n_matched++;

      if(n_matched != jets_had.size()){
        candidates.at(i).set_discriminators("correct_match", 9999999);
        // cout << "Not number of jets and matched equal" << endl;
        continue;
      }
    }
    else{
      const TopJet* topjet = candidates.at(i).tophad_topjet_ptr();

      // Match b-quark
      dr = deltaR(ttbargen.BHad(), *topjet);
      if(dr > 0.8){
        candidates.at(i).set_discriminators("correct_match", 9999999);
        // cout << "Not hadronic b-quark matched (TTAG)" << endl;
        continue;
      }
      correct_dr += dr;

      //match quarks from W decays to jets
      // First
      dr = deltaR(ttbargen.Q1(), *topjet);
      if(dr > 0.8){
        candidates.at(i).set_discriminators("correct_match", 9999999);
        // cout << "Not hadronic Q1 matched (TTAG)" << endl;
        continue;
      }
      correct_dr += dr;

      // Second
      dr = deltaR(ttbargen.Q2(), *topjet);
      if(dr > 0.8){
        candidates.at(i).set_discriminators("correct_match", 9999999);
        // cout << "Not hadronic Q2 matched (TTAG)" << endl;
        continue;
      }
      correct_dr += dr;
    }

    // Neutrino
    //correct_dr += deltaR(ttbargen.Neutrino(), candidates.at(i).neutrino_v4());
    dr = deltaPhi(ttbargen.Neutrino(), candidates.at(i).neutrino_v4());
    if(dr > 0.3){
      candidates.at(i).set_discriminators("correct_match", 9999999);
      continue;
    }
    correct_dr += dr;

    // Lepton
    dr = deltaR(ttbargen.ChargedLepton(), candidates.at(i).lepton());
    if(dr > 0.1){
      candidates.at(i).set_discriminators("correct_match", 9999999);
      continue;
    }
    correct_dr += dr;

    //cout << "dr = " << correct_dr << endl;
    candidates.at(i).set_discriminators("correct_match", correct_dr);

    if(correct_dr < dr_best){
      dr_best = correct_dr;
      bestCand = &candidates.at(i);
    }
    //cout << "best dr = " << dr_best << endl;
  }

  if(dr_best > 10.) return false;
  event.set(h_BestCandidate_, bestCand);
  event.set(h_is_zprime_reconstructed_, true);
  return true;
}

AK8PuppiTopTagger::AK8PuppiTopTagger(uhh2::Context& ctx, int min_num_daughters, float max_dR, float min_mass, float max_mass, float max_tau32) : min_num_daughters_(min_num_daughters), max_dR_(max_dR), min_mass_(min_mass), max_mass_(max_mass), max_tau32_(max_tau32) {

  h_AK8PuppiTopTags_ = ctx.get_handle< std::vector<TopJet> >("AK8PuppiTopTags");
  h_AK8PuppiTopTagsPtr_ = ctx.get_handle< std::vector<const TopJet*> >("AK8PuppiTopTagsPtr");

}

bool AK8PuppiTopTagger::process(uhh2::Event& event){

  std::vector<TopJet> toptags;
  vector<const TopJet*> toptags_ptr;
  for(const TopJet & puppijet : *event.toppuppijets){

    // 1) Jet should have at least two daughters
    int n_constit = 0;
    // Loop over subjets' constituents
    for(unsigned int i=0; i<puppijet.subjets().size(); i++){
      n_constit += puppijet.subjets().at(i).numberOfDaughters();
    }
    if(puppijet.numberOfDaughters() > (int)puppijet.subjets().size()){
      n_constit += (puppijet.numberOfDaughters() - puppijet.subjets().size());
    }

    if(n_constit<min_num_daughters_){
      continue;
    }

    // 3) Cut on SD mass
    LorentzVector SumSubjets(0.,0.,0.,0.);
    for(unsigned int k=0; k<puppijet.subjets().size(); k++) SumSubjets = SumSubjets + puppijet.subjets().at(k).v4();
    float mSD = SumSubjets.M();
    if(!(min_mass_ < mSD && mSD < max_mass_)) continue;

    // 4) Cut on tau 3/2
    if(!((puppijet.tau3() / puppijet.tau2()) < max_tau32_)) continue;

    // All jets at this point are top-tagged
    toptags.emplace_back(puppijet);
    toptags_ptr.emplace_back(&puppijet);
  }
  event.set(h_AK8PuppiTopTags_, toptags);
  event.set(h_AK8PuppiTopTagsPtr_, toptags_ptr);
  return (toptags.size() >= 1);
}



HOTVRTopTagger::HOTVRTopTagger(uhh2::Context& ctx) {

  h_HOTVRTopTags_ = ctx.get_handle< std::vector<TopJet> >("HOTVRTopTags");
  h_HOTVRTopTagsPtr_ = ctx.get_handle< std::vector<const TopJet*> >("HOTVRTopTagsPtr");

}

bool HOTVRTopTagger::process(uhh2::Event& event){

  std::vector<TopJet> toptags;
  vector<const TopJet*> toptags_ptr;
  for(const TopJet & topjet : *event.topjets){

    if (toptag_id(topjet, event)){
      toptags.emplace_back(topjet);
      toptags_ptr.emplace_back(&topjet);
    }
  }
  event.set(h_HOTVRTopTags_, toptags);
  event.set(h_HOTVRTopTagsPtr_, toptags_ptr);
  return (toptags.size() >= 1);
}


DeepAK8TopTagger::DeepAK8TopTagger(uhh2::Context& ctx){
  year = extract_year(ctx);
  h_DeepAK8TopTags_ = ctx.get_handle< std::vector<TopJet> >("DeepAK8TopTags");
  h_DeepAK8TopTagsPtr_ = ctx.get_handle< std::vector<const TopJet*> >("DeepAK8TopTagsPtr");
}

bool DeepAK8TopTagger::process(uhh2::Event& event){
  // values for UL: currently Christopher's private work
  double min_mSD = 105.;
  double max_mSD = 210.;
  double pt_min = 400.;
  double max_score;

  if(year == Year::isUL16preVFP) max_score = 0.485;
  else if(year == Year::isUL16postVFP) max_score = 0.475;
  else if(year == Year::isUL17) max_score = 0.487;
  else if(year == Year::isUL18) max_score = 0.477;
  else throw runtime_error("DeepAK8TopTagger: no valid year selected.");

  std::vector<TopJet> toptags;
  vector<const TopJet*> toptags_ptr;

  for(const TopJet & puppijet : *event.toppuppijets){
    // pT threshold
    if(!( puppijet.pt() > pt_min )) continue;

    // cut on SD mass
    LorentzVector SumSubjets(0.,0.,0.,0.);
    for(unsigned int k=0; k<puppijet.subjets().size(); k++) SumSubjets = SumSubjets + puppijet.subjets().at(k).v4();
    float mSD = SumSubjets.M();
    if(!(min_mSD < mSD && mSD < max_mSD)) continue;

    // cut on score
    if( !(puppijet.btag_MassDecorrelatedDeepBoosted_TvsQCD() >= max_score ) ) continue;

    toptags.emplace_back(puppijet);
    toptags_ptr.emplace_back(&puppijet);
  }

  event.set(h_DeepAK8TopTags_, toptags);
  event.set(h_DeepAK8TopTagsPtr_, toptags_ptr);
  return (toptags.size() >= 1);
}

bool JetLeptonDeltaRCleaner::process(uhh2::Event& event){

  assert(event.jets);
  std::vector<Jet> cleaned_jets;

  for(const auto & tjet : *event.jets){
    bool skip_tjet(false);

    if(event.muons){
      for(const auto & muo : *event.muons)
      if(uhh2::deltaR(tjet, muo) < minDR_) skip_tjet = true;
    }

    if(skip_tjet) continue;

    if(event.electrons){
      for(const auto & ele : *event.electrons)
      if(uhh2::deltaR(tjet, ele) < minDR_) skip_tjet = true;
    }

    if(!skip_tjet) cleaned_jets.push_back(tjet);
  }

  event.jets->clear();
  event.jets->reserve(cleaned_jets.size());
  for(const auto& j : cleaned_jets) event.jets->push_back(j);

  return true;
}
////

bool TopJetLeptonDeltaRCleaner::process(uhh2::Event& event){

  assert(event.topjets);
  std::vector<TopJet> cleaned_topjets;

  for(const auto & tjet : *event.topjets){
    bool skip_tjet(false);

    if(event.muons){
      for(const auto & muo : *event.muons)
      if(uhh2::deltaR(tjet, muo) < minDR_) skip_tjet = true;
    }

    if(skip_tjet) continue;

    if(event.electrons){
      for(const auto & ele : *event.electrons)
      if(uhh2::deltaR(tjet, ele) < minDR_) skip_tjet = true;
    }

    if(!skip_tjet) cleaned_topjets.push_back(tjet);
  }

  event.topjets->clear();
  event.topjets->reserve(cleaned_topjets.size());
  for(const auto& j : cleaned_topjets) event.topjets->push_back(j);

  return true;
}
////

const Particle* leading_lepton(const uhh2::Event& event){

  const Particle* lep(0);

  float ptL_max(0.);
  if(event.muons)    { for(const auto& mu : *event.muons)    { if(mu.pt() > ptL_max){ ptL_max = mu.pt(); lep = &mu; } } }
  if(event.electrons){ for(const auto& el : *event.electrons){ if(el.pt() > ptL_max){ ptL_max = el.pt(); lep = &el; } } }

  if(!lep) throw std::runtime_error("leading_lepton -- pt-leading lepton not found");

  return lep;
}

float STlep(const uhh2::Event& event){

  assert((event.muons || event.electrons));

  double stlep = 0.;
  for(const Electron & ele : *event.electrons) stlep += ele.pt();
  for(const Muon & mu : *event.muons)          stlep += mu.pt();

  return stlep;
}

float Muon_pfMINIIso(const Muon& muo, const uhh2::Event&, const std::string& iso_key_){

  float iso(-1.);

  if(!muo.pt()) throw std::runtime_error("Muon_pfMINIIso -- null muon transverse momentum: failed to calculate relative MINI-Isolation");

  if     (iso_key_ == "uncorrected") iso = (muo.pfMINIIso_CH() + muo.pfMINIIso_NH() + muo.pfMINIIso_Ph())/muo.pt();
  else if(iso_key_ == "delta-beta")  iso = (muo.pfMINIIso_CH() + std::max(0., muo.pfMINIIso_NH() + muo.pfMINIIso_Ph() - .5*muo.pfMINIIso_PU()))/muo.pt();
  else if(iso_key_ == "pf-weight")   iso = (muo.pfMINIIso_CH() + muo.pfMINIIso_NH_pfwgt() + muo.pfMINIIso_Ph_pfwgt())/muo.pt();
  else throw std::runtime_error("Muon_pfMINIIso -- invalid key for MINI-Isolation pileup correction: "+iso_key_);

  return iso;
}

float Electron_pfMINIIso(const Electron& ele, const uhh2::Event&, const std::string& iso_key_){

  float iso(-1.);

  if(!ele.pt()) throw std::runtime_error("Electron_pfMINIIso -- null muon transverse momentum: failed to calculate relative MINI-Isolation");

  if     (iso_key_ == "uncorrected") iso = (ele.pfMINIIso_CH() + ele.pfMINIIso_NH() + ele.pfMINIIso_Ph())/ele.pt();
  else if(iso_key_ == "delta-beta")  iso = (ele.pfMINIIso_CH() + std::max(0., ele.pfMINIIso_NH() + ele.pfMINIIso_Ph() - .5*ele.pfMINIIso_PU()))/ele.pt();
  else if(iso_key_ == "pf-weight")   iso = (ele.pfMINIIso_CH() + ele.pfMINIIso_NH_pfwgt() + ele.pfMINIIso_Ph_pfwgt())/ele.pt();
  else throw std::runtime_error("Electron_pfMINIIso -- invalid key for MINI-Isolation pileup correction: "+iso_key_);

  return iso;
}
////

bool trigger_bit(const uhh2::Event& evt_, const std::string& hlt_key_){

  uhh2::Event::TriggerIndex trg_index = evt_.get_trigger_index(hlt_key_);

  return bool(evt_.passes_trigger(trg_index));
}
////

GENWToLNuFinder::GENWToLNuFinder(uhh2::Context& ctx, const std::string& label){

  h_genWln_W_ = ctx.declare_event_output<GenParticle>(label+"_W");
  h_genWln_l_ = ctx.declare_event_output<GenParticle>(label+"_l");
  h_genWln_n_ = ctx.declare_event_output<GenParticle>(label+"_n");
}

bool GENWToLNuFinder::process(uhh2::Event& evt){

  const GenParticle *genW(0), *genW_lep(0), *genW_neu(0); {

    int genwN(0);

    assert(evt.genparticles);
    for(const auto& genp : *evt.genparticles){

      if(genwN > 1) throw std::runtime_error("GENWToLNuFinder::process -- logic error: more than 1 W->lnu decay found at GEN level");

      // gen-W
      const bool is_me = (20 <= genp.status() && genp.status() <= 30);
      if(!is_me) continue;

      const int is_w = (std::abs(genp.pdgId()) == 24);
      if(!is_w) continue;

      genW = &genp;

      // gen-W daughters
      const GenParticle* dau1 = genp.daughter(evt.genparticles, 1);
      const GenParticle* dau2 = genp.daughter(evt.genparticles, 2);

      if(dau1 && dau2){

        const int id1 = std::abs(dau1->pdgId());
        const int id2 = std::abs(dau2->pdgId());

        if((id1 == 11 || id1 == 13 || id1 == 15) &&
        (id2 == 12 || id2 == 14 || id2 == 16) ){ genW_lep = dau1; genW_neu = dau2; ++genwN; }

        if((id2 == 11 || id2 == 13 || id2 == 15) &&
        (id1 == 12 || id1 == 14 || id1 == 16) ){ genW_lep = dau2; genW_neu = dau1; ++genwN; }
      }
    }
  }

  //  if(!genW || !genW_lep || !genW_neu){
  //
  //    std::cout << std::endl;
  //    for(const auto& p : *evt.genparticles){
  //
  //      std::cout <<  " i=" << p.index()    ;
  //      std::cout << " m1=" << p.mother1()  ;
  //      std::cout << " m2=" << p.mother2()  ;
  //      std::cout << " d1=" << p.daughter1();
  //      std::cout << " d2=" << p.daughter2();
  //      std::cout << " ID=" << p.pdgId()    ;
  //      std::cout << " px=" << p.v4().Px()  ;
  //      std::cout << " py=" << p.v4().Py()  ;
  //      std::cout << " pz=" << p.v4().Pz()  ;
  //      std::cout << std::endl;
  //    }
  //  }

  evt.set(h_genWln_W_, genW     ? GenParticle(*genW)     : GenParticle());
  evt.set(h_genWln_l_, genW_lep ? GenParticle(*genW_lep) : GenParticle());
  evt.set(h_genWln_n_, genW_neu ? GenParticle(*genW_neu) : GenParticle());

  return true;
}
////

MEPartonFinder::MEPartonFinder(uhh2::Context& ctx, const std::string& label){

  h_meps_ = ctx.declare_event_output<std::vector<GenParticle> >(label);
}

bool MEPartonFinder::process(uhh2::Event& evt){

  std::vector<GenParticle> mep_refs;

  assert(evt.genparticles);
  for(const auto& genp : *evt.genparticles){

    const bool is_me = (20 <= genp.status() && genp.status() <= 30);
    if(!is_me) continue;

    const bool has_mo1 = (genp.mother1() != (unsigned short)(-1));
    if(!has_mo1) continue;

    const bool has_mo2 = (genp.mother2() != (unsigned short)(-1));
    if(!has_mo2) continue;

    const int abs_id = std::abs(genp.pdgId());
    const bool is_parton = (((1<=abs_id) && (abs_id<=5)) || (abs_id == 21));
    if(!is_parton) continue;

    mep_refs.push_back(genp);
  }

  evt.set(h_meps_, std::move(mep_refs));

  return true;
}

//////////////////////////////////////////////////////////////
/////////////// Variables for NN
//////////////////////////////////////////////////////////////

Variables_NN::Variables_NN(uhh2::Context& ctx, TString mode): mode_(mode){
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");
  h_eventweight = ctx.declare_event_output<float> ("eventweight");

  ///  MUONS
  h_Mu_pt = ctx.declare_event_output<float>("Mu_pt");
  h_Mu_eta = ctx.declare_event_output<float>("Mu_eta");
  h_Mu_phi = ctx.declare_event_output<float>("Mu_phi");
  h_Mu_E = ctx.declare_event_output<float>("Mu_E");

  ///  ELECTRONS
  h_Ele_pt = ctx.declare_event_output<float>("Ele_pt");
  h_Ele_eta = ctx.declare_event_output<float>("Ele_eta");
  h_Ele_phi = ctx.declare_event_output<float>("Ele_phi");
  h_Ele_E = ctx.declare_event_output<float>("Ele_E");

  ///  MET
  h_MET_pt = ctx.declare_event_output<float>("MET_pt");
  h_MET_phi = ctx.declare_event_output<float>("MET_phi");

  ///  AK4 JETS
  h_N_Ak4 = ctx.declare_event_output<float>("N_Ak4");

  h_Ak4_j1_pt = ctx.declare_event_output<float>("Ak4_j1_pt");
  h_Ak4_j1_eta = ctx.declare_event_output<float>("Ak4_j1_eta");
  h_Ak4_j1_phi = ctx.declare_event_output<float>("Ak4_j1_phi");
  h_Ak4_j1_E = ctx.declare_event_output<float>("Ak4_j1_E");
  h_Ak4_j1_m = ctx.declare_event_output<float>("Ak4_j1_m");
  h_Ak4_j1_deepjetbscore = ctx.declare_event_output<float>("Ak4_j1_deepjetbscore");

  h_Ak4_j2_pt = ctx.declare_event_output<float>("Ak4_j2_pt");
  h_Ak4_j2_eta = ctx.declare_event_output<float>("Ak4_j2_eta");
  h_Ak4_j2_phi = ctx.declare_event_output<float>("Ak4_j2_phi");
  h_Ak4_j2_E = ctx.declare_event_output<float>("Ak4_j2_E");
  h_Ak4_j2_m = ctx.declare_event_output<float>("Ak4_j2_m");
  h_Ak4_j2_deepjetbscore = ctx.declare_event_output<float>("Ak4_j2_deepjetbscore");

  h_Ak4_j3_pt = ctx.declare_event_output<float>("Ak4_j3_pt");
  h_Ak4_j3_eta = ctx.declare_event_output<float>("Ak4_j3_eta");
  h_Ak4_j3_phi = ctx.declare_event_output<float>("Ak4_j3_phi");
  h_Ak4_j3_E = ctx.declare_event_output<float>("Ak4_j3_E");
  h_Ak4_j3_m = ctx.declare_event_output<float>("Ak4_j3_m");
  h_Ak4_j3_deepjetbscore = ctx.declare_event_output<float>("Ak4_j3_deepjetbscore");

  h_Ak4_j4_pt = ctx.declare_event_output<float>("Ak4_j4_pt");
  h_Ak4_j4_eta = ctx.declare_event_output<float>("Ak4_j4_eta");
  h_Ak4_j4_phi = ctx.declare_event_output<float>("Ak4_j4_phi");
  h_Ak4_j4_E = ctx.declare_event_output<float>("Ak4_j4_E");
  h_Ak4_j4_m = ctx.declare_event_output<float>("Ak4_j4_m");
  h_Ak4_j4_deepjetbscore = ctx.declare_event_output<float>("Ak4_j4_deepjetbscore");

  h_Ak4_j5_pt = ctx.declare_event_output<float>("Ak4_j5_pt");
  h_Ak4_j5_eta = ctx.declare_event_output<float>("Ak4_j5_eta");
  h_Ak4_j5_phi = ctx.declare_event_output<float>("Ak4_j5_phi");
  h_Ak4_j5_E = ctx.declare_event_output<float>("Ak4_j5_E");
  h_Ak4_j5_m = ctx.declare_event_output<float>("Ak4_j5_m");
  h_Ak4_j5_deepjetbscore = ctx.declare_event_output<float>("Ak4_j5_deepjetbscore");

  h_Ak4_j6_pt = ctx.declare_event_output<float>("Ak4_j6_pt");
  h_Ak4_j6_eta = ctx.declare_event_output<float>("Ak4_j6_eta");
  h_Ak4_j6_phi = ctx.declare_event_output<float>("Ak4_j6_phi");
  h_Ak4_j6_E = ctx.declare_event_output<float>("Ak4_j6_E");
  h_Ak4_j6_m = ctx.declare_event_output<float>("Ak4_j6_m");
  h_Ak4_j6_deepjetbscore = ctx.declare_event_output<float>("Ak4_j6_deepjetbscore");

  /// AK8 JETS
  if(mode_ == "deepAK8"){
    h_N_Ak8 = ctx.declare_event_output<float>("N_Ak8");

    h_Ak8_j1_pt = ctx.declare_event_output<float>("Ak8_j1_pt");
    h_Ak8_j1_eta = ctx.declare_event_output<float>("Ak8_j1_eta");
    h_Ak8_j1_phi = ctx.declare_event_output<float>("Ak8_j1_phi");
    h_Ak8_j1_E = ctx.declare_event_output<float>("Ak8_j1_E");
    h_Ak8_j1_mSD = ctx.declare_event_output<float>("Ak8_j1_mSD");
    h_Ak8_j1_tau21 = ctx.declare_event_output<float>("Ak8_j1_tau21");
    h_Ak8_j1_tau32 = ctx.declare_event_output<float>("Ak8_j1_tau32");
    h_Ak8_j1_deepak8tscore = ctx.declare_event_output<float>("Ak8_j1_deepak8tscore");

    h_Ak8_j2_pt = ctx.declare_event_output<float>("Ak8_j2_pt");
    h_Ak8_j2_eta = ctx.declare_event_output<float>("Ak8_j2_eta");
    h_Ak8_j2_phi = ctx.declare_event_output<float>("Ak8_j2_phi");
    h_Ak8_j2_E = ctx.declare_event_output<float>("Ak8_j2_E");
    h_Ak8_j2_mSD = ctx.declare_event_output<float>("Ak8_j2_mSD");
    h_Ak8_j2_tau21 = ctx.declare_event_output<float>("Ak8_j2_tau21");
    h_Ak8_j2_tau32 = ctx.declare_event_output<float>("Ak8_j2_tau32");
    h_Ak8_j2_deepak8tscore = ctx.declare_event_output<float>("Ak8_j2_deepak8tscore");

    h_Ak8_j3_pt = ctx.declare_event_output<float>("Ak8_j3_pt");
    h_Ak8_j3_eta = ctx.declare_event_output<float>("Ak8_j3_eta");
    h_Ak8_j3_phi = ctx.declare_event_output<float>("Ak8_j3_phi");
    h_Ak8_j3_E = ctx.declare_event_output<float>("Ak8_j3_E");
    h_Ak8_j3_mSD = ctx.declare_event_output<float>("Ak8_j3_mSD");
    h_Ak8_j3_tau21 = ctx.declare_event_output<float>("Ak8_j3_tau21");
    h_Ak8_j3_tau32 = ctx.declare_event_output<float>("Ak8_j3_tau32");
    h_Ak8_j3_deepak8tscore = ctx.declare_event_output<float>("Ak8_j3_deepak8tscore");
  }
  else if(mode_ == "hotvr"){
    ///  HOTVR JETS
    h_N_HOTVR = ctx.declare_event_output<float> ("N_HOTVR");

    h_HOTVR_j1_pt = ctx.declare_event_output<float>("HOTVR_j1_pt");
    h_HOTVR_j1_eta = ctx.declare_event_output<float>("HOTVR_j1_eta");
    h_HOTVR_j1_phi = ctx.declare_event_output<float>("HOTVR_j1_phi");
    h_HOTVR_j1_E = ctx.declare_event_output<float>("HOTVR_j1_E");
    h_HOTVR_j1_mSD = ctx.declare_event_output<float>("HOTVR_j1_mSD");
    h_HOTVR_j1_tau21 = ctx.declare_event_output<float>("HOTVR_j1_tau21");
    h_HOTVR_j1_tau32 = ctx.declare_event_output<float>("HOTVR_j1_tau32");

    h_HOTVR_j2_pt = ctx.declare_event_output<float>("HOTVR_j2_pt");
    h_HOTVR_j2_eta = ctx.declare_event_output<float>("HOTVR_j2_eta");
    h_HOTVR_j2_phi = ctx.declare_event_output<float>("HOTVR_j2_phi");
    h_HOTVR_j2_E = ctx.declare_event_output<float>("HOTVR_j2_E");
    h_HOTVR_j2_mSD = ctx.declare_event_output<float>("HOTVR_j2_mSD");
    h_HOTVR_j2_tau21 = ctx.declare_event_output<float>("HOTVR_j2_tau21");
    h_HOTVR_j2_tau32 = ctx.declare_event_output<float>("HOTVR_j2_tau32");

    h_HOTVR_j3_pt = ctx.declare_event_output<float>("HOTVR_j3_pt");
    h_HOTVR_j3_eta = ctx.declare_event_output<float>("HOTVR_j3_eta");
    h_HOTVR_j3_phi = ctx.declare_event_output<float>("HOTVR_j3_phi");
    h_HOTVR_j3_E = ctx.declare_event_output<float>("HOTVR_j3_E");
    h_HOTVR_j3_mSD = ctx.declare_event_output<float>("HOTVR_j3_mSD");
    h_HOTVR_j3_tau21 = ctx.declare_event_output<float>("HOTVR_j3_tau21");
    h_HOTVR_j3_tau32 = ctx.declare_event_output<float>("HOTVR_j3_tau32");
  }

  ///  M ttbar
  h_M_tt = ctx.declare_event_output<float>("M_tt");
  //Ac and spic corr
  h_dyreco = ctx.declare_event_output<float>("dyreco");
  h_dyreco_1 = ctx.declare_event_output<float>("dyreco_1");
  h_dyreco_2 = ctx.declare_event_output<float>("dyreco_2");

  h_Sigma_phi_1=ctx.declare_event_output<float>("Sigma_phi_1");
  h_Sigma_phi_2=ctx.declare_event_output<float>("Sigma_phi_2");
  h_Sigma_phi = ctx.declare_event_output<float>("Sigma_phi");
  h_Delta_phi = ctx.declare_event_output<float>("Delta_phi");

  
  // chi^2
  h_chi2 = ctx.declare_event_output<float>("chi2");

  // random number for k-fold validation
  h_uniform_random = ctx.declare_event_output<float>("uniform_random");

}

bool Variables_NN::process(uhh2::Event& evt){

  double weight = evt.weight;
  evt.set(h_eventweight, -10);
  evt.set(h_eventweight, weight);

  /////////   MUONS
  evt.set(h_Mu_pt, -10);
  evt.set(h_Mu_eta,-10);
  evt.set(h_Mu_phi, -10);
  evt.set(h_Mu_E, -10);

  vector<Muon>* muons = evt.muons;
  int Nmuons = muons->size();

  for(int i=0; i<Nmuons; i++){
    evt.set(h_Mu_pt, muons->at(i).pt());
    evt.set(h_Mu_eta, muons->at(i).eta());
    evt.set(h_Mu_phi, muons->at(i).phi());
    evt.set(h_Mu_E, muons->at(i).energy());
  }


  /////////   ELECTRONS
  evt.set(h_Ele_pt, -10);
  evt.set(h_Ele_eta, -10);
  evt.set(h_Ele_phi, -10);
  evt.set(h_Ele_E, -10);

  vector<Electron>* electrons = evt.electrons;
  int Nelectrons = electrons->size();

  for(int i=0; i<Nelectrons; i++){
    evt.set(h_Ele_pt, electrons->at(i).pt());
    evt.set(h_Ele_eta, electrons->at(i).eta());
    evt.set(h_Ele_phi, electrons->at(i).phi());
    evt.set(h_Ele_E, electrons->at(i).energy());
  }

  /////////   MET
  evt.set(h_MET_pt, -10);
  evt.set(h_MET_phi, -10);

  evt.set(h_MET_pt, evt.met->pt());
  evt.set(h_MET_phi, evt.met->phi());


  ///////// AK4 JETS
  evt.set(h_N_Ak4, -10);

  evt.set(h_Ak4_j1_pt, -10);
  evt.set(h_Ak4_j1_eta, -10);
  evt.set(h_Ak4_j1_phi, -10);
  evt.set(h_Ak4_j1_E, -10);
  evt.set(h_Ak4_j1_m, -10);
  evt.set(h_Ak4_j1_deepjetbscore, -10);

  evt.set(h_Ak4_j2_pt, -10);
  evt.set(h_Ak4_j2_eta, -10);
  evt.set(h_Ak4_j2_phi, -10);
  evt.set(h_Ak4_j2_E, -10);
  evt.set(h_Ak4_j2_m, -10);
  evt.set(h_Ak4_j2_deepjetbscore, -10);

  evt.set(h_Ak4_j3_pt, -10);
  evt.set(h_Ak4_j3_eta, -10);
  evt.set(h_Ak4_j3_phi, -10);
  evt.set(h_Ak4_j3_E, -10);
  evt.set(h_Ak4_j3_m, -10);
  evt.set(h_Ak4_j3_deepjetbscore, -10);

  evt.set(h_Ak4_j4_pt, -10);
  evt.set(h_Ak4_j4_eta, -10);
  evt.set(h_Ak4_j4_phi, -10);
  evt.set(h_Ak4_j4_E, -10);
  evt.set(h_Ak4_j4_m, -10);
  evt.set(h_Ak4_j4_deepjetbscore, -10);

  evt.set(h_Ak4_j5_pt, -10);
  evt.set(h_Ak4_j5_eta, -10);
  evt.set(h_Ak4_j5_phi, -10);
  evt.set(h_Ak4_j5_E, -10);
  evt.set(h_Ak4_j5_m, -10);
  evt.set(h_Ak4_j5_deepjetbscore, -10);

  evt.set(h_Ak4_j6_pt, -10);
  evt.set(h_Ak4_j6_eta, -10);
  evt.set(h_Ak4_j6_phi, -10);
  evt.set(h_Ak4_j6_E, -10);
  evt.set(h_Ak4_j6_m, -10);
  evt.set(h_Ak4_j6_deepjetbscore, -10);


  vector<Jet>* Ak4jets = evt.jets;
  int NAk4jets = Ak4jets->size();
  evt.set(h_N_Ak4, NAk4jets);

  for(int i=0; i<NAk4jets; i++){
    if(i==0){
      evt.set(h_Ak4_j1_pt, Ak4jets->at(i).pt());
      evt.set(h_Ak4_j1_eta, Ak4jets->at(i).eta());
      evt.set(h_Ak4_j1_phi, Ak4jets->at(i).phi());
      evt.set(h_Ak4_j1_E, Ak4jets->at(i).energy());
      evt.set(h_Ak4_j1_m, Ak4jets->at(i).v4().M());
      //evt.set(h_Ak4_j1_deepjetbscore, Ak4jets->at(i).btag_DeepJet());
    }
    if(i==1){
      evt.set(h_Ak4_j2_pt, Ak4jets->at(i).pt());
      evt.set(h_Ak4_j2_eta, Ak4jets->at(i).eta());
      evt.set(h_Ak4_j2_phi, Ak4jets->at(i).phi());
      evt.set(h_Ak4_j2_E, Ak4jets->at(i).energy());
      evt.set(h_Ak4_j2_m, Ak4jets->at(i).v4().M());
      //evt.set(h_Ak4_j2_deepjetbscore, Ak4jets->at(i).btag_DeepJet());
    }
    if(i==2){
      evt.set(h_Ak4_j3_pt, Ak4jets->at(i).pt());
      evt.set(h_Ak4_j3_eta, Ak4jets->at(i).eta());
      evt.set(h_Ak4_j3_phi, Ak4jets->at(i).phi());
      evt.set(h_Ak4_j3_E, Ak4jets->at(i).energy());
      evt.set(h_Ak4_j3_m, Ak4jets->at(i).v4().M());
      //evt.set(h_Ak4_j3_deepjetbscore, Ak4jets->at(i).btag_DeepJet());
    }
    if(i==3){
      evt.set(h_Ak4_j4_pt, Ak4jets->at(i).pt());
      evt.set(h_Ak4_j4_eta, Ak4jets->at(i).eta());
      evt.set(h_Ak4_j4_phi, Ak4jets->at(i).phi());
      evt.set(h_Ak4_j4_E, Ak4jets->at(i).energy());
      evt.set(h_Ak4_j4_m, Ak4jets->at(i).v4().M());
      //evt.set(h_Ak4_j4_deepjetbscore, Ak4jets->at(i).btag_DeepJet());
    }
    if(i==4){
      evt.set(h_Ak4_j5_pt, Ak4jets->at(i).pt());
      evt.set(h_Ak4_j5_eta, Ak4jets->at(i).eta());
      evt.set(h_Ak4_j5_phi, Ak4jets->at(i).phi());
      evt.set(h_Ak4_j5_E, Ak4jets->at(i).energy());
      evt.set(h_Ak4_j5_m, Ak4jets->at(i).v4().M());
      //evt.set(h_Ak4_j5_deepjetbscore, Ak4jets->at(i).btag_DeepJet());
    }
    if(i==5){
      evt.set(h_Ak4_j6_pt, Ak4jets->at(i).pt());
      evt.set(h_Ak4_j6_eta, Ak4jets->at(i).eta());
      evt.set(h_Ak4_j6_phi, Ak4jets->at(i).phi());
      evt.set(h_Ak4_j6_E, Ak4jets->at(i).energy());
      evt.set(h_Ak4_j6_m, Ak4jets->at(i).v4().M());
      //evt.set(h_Ak4_j6_deepjetbscore, Ak4jets->at(i).btag_DeepJet());
    }
  }

  // save b-tag score of matched CHS jet
  vector<Jet> AK4CHSjets_matched = evt.get(h_CHSjets_matched);
  for(unsigned int i=0; i<AK4CHSjets_matched.size(); i++){
    if(i==0){
      evt.set(h_Ak4_j1_deepjetbscore, AK4CHSjets_matched.at(i).btag_DeepJet());
    }
    if(i==1){
      evt.set(h_Ak4_j2_deepjetbscore, AK4CHSjets_matched.at(i).btag_DeepJet());
    }
    if(i==2){
      evt.set(h_Ak4_j3_deepjetbscore, AK4CHSjets_matched.at(i).btag_DeepJet());
    }
    if(i==3){
      evt.set(h_Ak4_j4_deepjetbscore, AK4CHSjets_matched.at(i).btag_DeepJet());
    }
    if(i==4){
      evt.set(h_Ak4_j5_deepjetbscore, AK4CHSjets_matched.at(i).btag_DeepJet());
    }
    if(i==5){
      evt.set(h_Ak4_j6_deepjetbscore, AK4CHSjets_matched.at(i).btag_DeepJet());
    }
  }

  /////////   AK8 JETS
  if(mode_ == "deepAK8"){
    evt.set(h_N_Ak8, -10);

    evt.set(h_Ak8_j1_pt, -10);
    evt.set(h_Ak8_j1_eta, -10);
    evt.set(h_Ak8_j1_phi, -10);
    evt.set(h_Ak8_j1_E, -10);
    evt.set(h_Ak8_j1_mSD, -10);
    evt.set(h_Ak8_j1_tau21, -10);
    evt.set(h_Ak8_j1_tau32, -10);
    evt.set(h_Ak8_j1_deepak8tscore, -10);

    evt.set(h_Ak8_j2_pt, -10);
    evt.set(h_Ak8_j2_eta, -10);
    evt.set(h_Ak8_j2_phi, -10);
    evt.set(h_Ak8_j2_E, -10);
    evt.set(h_Ak8_j2_mSD, -10);
    evt.set(h_Ak8_j2_tau21, -10);
    evt.set(h_Ak8_j2_tau32, -10);
    evt.set(h_Ak8_j2_deepak8tscore, -10);

    evt.set(h_Ak8_j3_pt, -10);
    evt.set(h_Ak8_j3_eta, -10);
    evt.set(h_Ak8_j3_phi, -10);
    evt.set(h_Ak8_j3_E, -10);
    evt.set(h_Ak8_j3_mSD, -10);
    evt.set(h_Ak8_j3_tau21, -10);
    evt.set(h_Ak8_j3_tau32, -10);
    evt.set(h_Ak8_j3_deepak8tscore, -10);

    vector<TopJet>* Ak8jets = evt.toppuppijets;
    int NAk8jets = Ak8jets->size();
    evt.set(h_N_Ak8, NAk8jets);

    for(int i=0; i<NAk8jets; i++){
      if(i==0){
        evt.set(h_Ak8_j1_pt, Ak8jets->at(i).pt());
        evt.set(h_Ak8_j1_eta, Ak8jets->at(i).eta());
        evt.set(h_Ak8_j1_phi, Ak8jets->at(i).phi());
        evt.set(h_Ak8_j1_E, Ak8jets->at(i).energy());
        evt.set(h_Ak8_j1_mSD, Ak8jets->at(i).softdropmass());
        evt.set(h_Ak8_j1_tau21, Ak8jets->at(i).tau2()/Ak8jets->at(i).tau1());
        evt.set(h_Ak8_j1_tau32, Ak8jets->at(i).tau3()/Ak8jets->at(i).tau2());
        evt.set(h_Ak8_j1_deepak8tscore, Ak8jets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD());
      }
      if(i==1){
        evt.set(h_Ak8_j2_pt, Ak8jets->at(i).pt());
        evt.set(h_Ak8_j2_eta, Ak8jets->at(i).eta());
        evt.set(h_Ak8_j2_phi, Ak8jets->at(i).phi());
        evt.set(h_Ak8_j2_E, Ak8jets->at(i).energy());
        evt.set(h_Ak8_j2_mSD, Ak8jets->at(i).softdropmass());
        evt.set(h_Ak8_j2_tau21, Ak8jets->at(i).tau2()/Ak8jets->at(i).tau1());
        evt.set(h_Ak8_j2_tau32, Ak8jets->at(i).tau3()/Ak8jets->at(i).tau2());
        evt.set(h_Ak8_j2_deepak8tscore, Ak8jets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD());
      }
      if(i==2){
        evt.set(h_Ak8_j3_pt, Ak8jets->at(i).pt());
        evt.set(h_Ak8_j3_eta, Ak8jets->at(i).eta());
        evt.set(h_Ak8_j3_phi, Ak8jets->at(i).phi());
        evt.set(h_Ak8_j3_E, Ak8jets->at(i).energy());
        evt.set(h_Ak8_j3_mSD, Ak8jets->at(i).softdropmass());
        evt.set(h_Ak8_j3_tau21, Ak8jets->at(i).tau2()/Ak8jets->at(i).tau1());
        evt.set(h_Ak8_j3_tau32, Ak8jets->at(i).tau3()/Ak8jets->at(i).tau2());
        evt.set(h_Ak8_j3_deepak8tscore, Ak8jets->at(i).btag_MassDecorrelatedDeepBoosted_TvsQCD());
      }
    }
  } // end deepAK8 mode

  /////////   HOTVR JETS
  if(mode_ == "hotvr"){
    evt.set(h_N_HOTVR, -10);

    evt.set(h_HOTVR_j1_pt, -10);
    evt.set(h_HOTVR_j1_eta, -10);
    evt.set(h_HOTVR_j1_phi, -10);
    evt.set(h_HOTVR_j1_E, -10);
    evt.set(h_HOTVR_j1_mSD, -10);
    evt.set(h_HOTVR_j1_tau21, -10);
    evt.set(h_HOTVR_j1_tau32, -10);

    evt.set(h_HOTVR_j2_pt, -10);
    evt.set(h_HOTVR_j2_eta, -10);
    evt.set(h_HOTVR_j2_phi, -10);
    evt.set(h_HOTVR_j2_E, -10);
    evt.set(h_HOTVR_j2_mSD, -10);
    evt.set(h_HOTVR_j2_tau21, -10);
    evt.set(h_HOTVR_j2_tau32, -10);

    evt.set(h_HOTVR_j3_pt, -10);
    evt.set(h_HOTVR_j3_eta, -10);
    evt.set(h_HOTVR_j3_phi, -10);
    evt.set(h_HOTVR_j3_E, -10);
    evt.set(h_HOTVR_j3_mSD, -10);
    evt.set(h_HOTVR_j3_tau21, -10);
    evt.set(h_HOTVR_j3_tau32, -10);


    vector<TopJet>* HOTVRjets = evt.topjets;
    int NHOTVRjets = HOTVRjets->size();
    evt.set(h_N_HOTVR, NHOTVRjets);

    for(int i=0; i<NHOTVRjets; i++){
      if(i==0){
        evt.set(h_HOTVR_j1_pt, HOTVRjets->at(i).pt());
        evt.set(h_HOTVR_j1_eta, HOTVRjets->at(i).eta());
        evt.set(h_HOTVR_j1_phi, HOTVRjets->at(i).phi());
        evt.set(h_HOTVR_j1_E, HOTVRjets->at(i).energy());
        evt.set(h_HOTVR_j1_mSD, HOTVRjets->at(i).v4().M());
        evt.set(h_HOTVR_j1_tau21, HOTVRjets->at(i).tau2_groomed()/HOTVRjets->at(i).tau1_groomed());
        evt.set(h_HOTVR_j1_tau32, HOTVRjets->at(i).tau3_groomed()/HOTVRjets->at(i).tau2_groomed());
      }
      if(i==1){
        evt.set(h_HOTVR_j2_pt, HOTVRjets->at(i).pt());
        evt.set(h_HOTVR_j2_eta, HOTVRjets->at(i).eta());
        evt.set(h_HOTVR_j2_phi, HOTVRjets->at(i).phi());
        evt.set(h_HOTVR_j2_E, HOTVRjets->at(i).energy());
        evt.set(h_HOTVR_j2_mSD, HOTVRjets->at(i).v4().M());
        evt.set(h_HOTVR_j2_tau21, HOTVRjets->at(i).tau2_groomed()/HOTVRjets->at(i).tau1_groomed());
        evt.set(h_HOTVR_j2_tau32, HOTVRjets->at(i).tau3_groomed()/HOTVRjets->at(i).tau2_groomed());
      }
      if(i==2){
        evt.set(h_HOTVR_j3_pt, HOTVRjets->at(i).pt());
        evt.set(h_HOTVR_j3_eta, HOTVRjets->at(i).eta());
        evt.set(h_HOTVR_j3_phi, HOTVRjets->at(i).phi());
        evt.set(h_HOTVR_j3_E, HOTVRjets->at(i).energy());
        evt.set(h_HOTVR_j3_mSD, HOTVRjets->at(i).v4().M());
        evt.set(h_HOTVR_j3_tau21, HOTVRjets->at(i).tau2_groomed()/HOTVRjets->at(i).tau1_groomed());
        evt.set(h_HOTVR_j3_tau32, HOTVRjets->at(i).tau3_groomed()/HOTVRjets->at(i).tau2_groomed());
      }
    }
  } // end hotvr mode

  // ttbar mass
  evt.set(h_M_tt, -10);
  evt.set(h_chi2, -10);
  bool is_zprime_reconstructed_chi2 = evt.get(h_is_zprime_reconstructed_chi2);
  ZprimeCandidate* BestZprimeCandidate = evt.get(h_BestZprimeCandidateChi2);
  bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction(); // Reconstruction process id
  // vector <Jet> AK4CHSjets_matched = event.get(h_CHSjets_matched);  
  //                 // AK4Puppijets that have been matched to CHSjets
 
  // vector <TopJet> TopTaggedJets = evt.get(h_AK8TopTags);                     // AK8Puppi jets TopTagged by DeepAK8TopTagger
  vector <float> jets_hadronic_bscores;                                            // bScores vector for resolved hadronic jets
  float pt_hadTop_thresh = 150;                                                    // Define cut-variable as pt of hadTop for low/high regions                                                   // medium WP for UL18 DeepJet

  if(is_zprime_reconstructed_chi2){
    // ZprimeCandidate* BestZprimeCandidate = evt.get(h_BestZprimeCandidateChi2);
    float Mass_tt = BestZprimeCandidate->Zprime_v4().M();
    float chi2 = BestZprimeCandidate->discriminator("chi2_total");
    evt.set(h_M_tt, Mass_tt);
    evt.set(h_chi2, chi2);
  }
  // EFT Ac and spin correlation variables:
  // Plot pt of hadronic Top jet
  float pt_hadTop = BestZprimeCandidate->top_hadronic_v4().pt();
  
  float bscore_max = -2;
  if(!is_toptag_reconstruction){
      // Loop over resolved hadronic jets to find their bscore via CHS jets
    for(unsigned int i=0; i<BestZprimeCandidate->jets_hadronic().size(); i++){
      double deltaR_min = 99;
      // Match resolved hadronic jets to CHS jets (which have bscores)
      for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
        double deltaR_CHS = deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(j));
        if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;
        }
      // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
      for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
        if(deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
        jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet());
        } // Using DeepJet btag score
    }
    // Loop over bScores-vector to extract highest bscor
    for(unsigned int i=0; i<jets_hadronic_bscores.size(); i++){
      float bscore = jets_hadronic_bscores.at(i);
      if(bscore > bscore_max) bscore_max = bscore;
    }
    //is not top tag
  }
  if(is_toptag_reconstruction){
    // Loop over hadronic top's subjets to extract highest bscore
  for(unsigned int i=0; i < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); i++){
    float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet(); // Using DeepJet btag score
    if(bscore > bscore_max) bscore_max = bscore;
  }

}
  TLorentzVector had_top_b(0, 0, 0, 0);

// Resolved topology
  if(!is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK4-jet with highest bscore
    for(unsigned int i=0; i< BestZprimeCandidate->jets_hadronic().size(); i++){
      float bscore = jets_hadronic_bscores.at(i);
      if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->jets_hadronic().at(i).pt(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).eta(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).phi(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).energy());
    }
  }
  // Merged topology
  if(is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK8-subjet with highest bscore
    for(unsigned int j=0; j < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); j++){
      float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
      if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).pt(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).eta(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).phi(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).energy());
    }
  }
  TLorentzVector lep_top_lep(0, 0, 0, 0);
  LorentzVector lep = BestZprimeCandidate->lepton().v4();
  lep_top_lep.SetPtEtaPhiE(lep.pt(), lep.eta(), lep.phi(), lep.E());
    // Define 4vectors of top quarks
  TLorentzVector PosTop(0, 0, 0, 0);
  TLorentzVector NegTop(0, 0, 0, 0);

      // POSITIVE LEPTON CONFIGURATION => Positive charged lepton has Positive Top mother
  if(BestZprimeCandidate->lepton().charge() > 0){
    // Define ttbar system
    PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                        BestZprimeCandidate->top_leptonic_v4().eta(), 
                        BestZprimeCandidate->top_leptonic_v4().phi(), 
                        BestZprimeCandidate->top_leptonic_v4().energy());
    NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                        BestZprimeCandidate->top_hadronic_v4().eta(), 
                        BestZprimeCandidate->top_hadronic_v4().phi(), 
                        BestZprimeCandidate->top_hadronic_v4().energy());
  }
  else if (BestZprimeCandidate->lepton().charge() < 0){
    PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                        BestZprimeCandidate->top_hadronic_v4().eta(), 
                        BestZprimeCandidate->top_hadronic_v4().phi(), 
                        BestZprimeCandidate->top_hadronic_v4().energy());
    NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                        BestZprimeCandidate->top_leptonic_v4().eta(), 
                        BestZprimeCandidate->top_leptonic_v4().phi(), 
                        BestZprimeCandidate->top_leptonic_v4().energy());
  }
      
  TLorentzVector ttbar = PosTop + NegTop;
  TLorentzVector lep_top_lep_CoM = lep_top_lep;
  // Boost into ttbar CoM-Frame <<<-------//
  lep_top_lep_CoM.Boost(-1.*ttbar.BoostVector());
  TLorentzVector had_top_b_CoM = had_top_b;
  had_top_b_CoM.Boost(-ttbar.BoostVector());
  TLorentzVector PosTop_CoM = PosTop;
  PosTop_CoM.Boost(-ttbar.BoostVector());
  TLorentzVector NegTop_CoM = NegTop;
  NegTop_CoM.Boost(-ttbar.BoostVector());
  // Beam unit vector in COM frame
  TVector3 beam_axis(0,0,1);
  // Calculating top scattering angle for PosTop only
  double cos_PosTop_beam = PosTop_CoM.Vect().Unit().Dot(beam_axis);
  double sin_PosTop_beam = sqrt(1 - cos_PosTop_beam*cos_PosTop_beam);

  // The sign of cos_PosTop_beam to account for Bose symmetry
  double sign_cos_PosTop_beam = (cos_PosTop_beam > 0.) ? 1. : -1.;
  // // The sign based on PosTop and NegTop's rapidity
  // double sign_rapidity = (PosTop.Rapidity() >= NegTop.Rapidity()) ? 1. : -1.;

  // Bernreuther basis vectors
  TVector3 kbase = PosTop_CoM.Vect().Unit();
  TVector3 rbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*(beam_axis - cos_PosTop_beam * kbase) ).Unit();
  TVector3 nbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*beam_axis.Cross(kbase) ).Unit();

  // Rotate vectors into Helicity Frame <<<-----------//
  // Rotate about beamline
  TLorentzVector lep_top_lep_H = lep_top_lep_CoM;
  lep_top_lep_H.RotateZ(-1.*PosTop_CoM.Phi());
  TLorentzVector had_top_b_H = had_top_b_CoM;
  had_top_b_H.RotateZ(-1.*PosTop_CoM.Phi());
  TLorentzVector PosTop_H = PosTop_CoM;
  PosTop_H.RotateZ(-1.*PosTop_CoM.Phi());
  TLorentzVector NegTop_H = NegTop_CoM;
  NegTop_H.RotateZ(-1.*PosTop_CoM.Phi());

  TVector3 kbase_H = kbase;
  kbase_H.RotateZ(-1.*PosTop_CoM.Phi());
  TVector3 rbase_H = rbase;
  rbase_H.RotateZ(-1.*PosTop_CoM.Phi());
  TVector3 nbase_H = nbase;
  nbase_H.RotateZ(-1.*PosTop_CoM.Phi());

  // Rotate about y-axis
  TLorentzVector lep_top_lep_Hel = lep_top_lep_H;
  lep_top_lep_Hel.RotateY(-1.*PosTop_CoM.Theta());
  TLorentzVector had_top_b_Hel = had_top_b_H;
  had_top_b_Hel.RotateY(-1.*PosTop_CoM.Theta());
  TLorentzVector PosTop_Hel = PosTop_H;
  PosTop_Hel.RotateY(-1.*PosTop_CoM.Theta());
  TLorentzVector NegTop_Hel = NegTop_H;
  NegTop_Hel.RotateY(-1.*PosTop_CoM.Theta());

  TVector3 kbase_Hel = kbase_H;
  kbase_Hel.RotateY(-1.*PosTop_CoM.Theta());
  TVector3 rbase_Hel = rbase_H;
  rbase_Hel.RotateY(-1.*PosTop_CoM.Theta());
  TVector3 nbase_Hel = nbase_H;
  nbase_Hel.RotateY(-1.*PosTop_CoM.Theta());

  // Rotation to align with Bernreuther basis <<<---------//
  TLorentzVector lep_top_lep_BoseSymm = lep_top_lep_Hel;
  TLorentzVector had_top_b_BoseSymm = had_top_b_Hel;
  TLorentzVector PosTop_BoseSymm = PosTop_Hel;
  TLorentzVector NegTop_BoseSymm = NegTop_Hel;

  TVector3 kbase_BoseSymm = kbase_Hel;
  TVector3 rbase_BoseSymm = rbase_Hel;
  TVector3 nbase_BoseSymm = nbase_Hel;

  // if(sign_cos_PosTop_beam > 0.){
  //   lep_top_lep_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
  //   had_top_b_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
  //   PosTop_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
  //   NegTop_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);

  //   kbase_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
  //   rbase_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
  //   nbase_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
  // }
  // else{
  //   lep_top_lep_BoseSymm.RotateZ(TMath::Pi()/2.);
  //   had_top_b_BoseSymm.RotateZ(TMath::Pi()/2.);
  //   PosTop_BoseSymm.RotateZ(TMath::Pi()/2.);
  //   NegTop_BoseSymm.RotateZ(TMath::Pi()/2.);

  //   kbase_BoseSymm.RotateZ(TMath::Pi()/2.);
  //   rbase_BoseSymm.RotateZ(TMath::Pi()/2.);
  //   nbase_BoseSymm.RotateZ(TMath::Pi()/2.);
  // }
   // Boosting into ttbar rest-frame <<<-------------------------------------------------------//
  TLorentzVector lep_top_lep_Rest = lep_top_lep_BoseSymm;
  TLorentzVector had_top_b_Rest = had_top_b_BoseSymm;
  TLorentzVector PosTop_Rest = PosTop_BoseSymm;
  TLorentzVector NegTop_Rest = NegTop_BoseSymm;
  if(BestZprimeCandidate->lepton().charge() > 0){
    lep_top_lep_Rest.Boost(-1.*PosTop_BoseSymm.BoostVector()); // lepton has Positive Top mother
    had_top_b_Rest.Boost(-1.*NegTop_BoseSymm.BoostVector());   // b-jet has Negative Top mother
  }
  else if (BestZprimeCandidate->lepton().charge() < 0){
    lep_top_lep_Rest.Boost(-1.*NegTop_BoseSymm.BoostVector()); // lepton has Negative Top mother
    had_top_b_Rest.Boost(-1.*PosTop_BoseSymm.BoostVector());   // b-jet has Positive Top mother
  }

  float dphi=0.;
  float sphi = lep_top_lep_Rest.Phi() + had_top_b_Rest.Phi();
  if(BestZprimeCandidate->lepton().charge() > 0){ // lepton is Positive Top's Decay Product
    dphi = lep_top_lep_Rest.Phi() - had_top_b_Rest.Phi();
  }
  if(BestZprimeCandidate->lepton().charge() < 0){
    dphi = had_top_b_Rest.Phi() - lep_top_lep_Rest.Phi();
  }

  evt.set(h_Sigma_phi,-10);
  evt.set(h_Delta_phi,-10);
  evt.set(h_Sigma_phi_1,-10);
  evt.set(h_Sigma_phi_2,-10);
  evt.set(h_dyreco_1,-10);
  evt.set(h_dyreco_2,-10);


  if(sphi > TMath::Pi()) sphi = sphi - 2*TMath::Pi();
  if(sphi < -TMath::Pi()) sphi = sphi + 2*TMath::Pi();
  if(dphi > TMath::Pi()) dphi = dphi - 2*TMath::Pi();
  if(dphi < -TMath::Pi()) dphi = dphi + 2*TMath::Pi();
  
  evt.set(h_Sigma_phi,sphi);
  evt.set(h_Delta_phi,dphi);
  



  evt.set(h_dyreco,-10);
  float dyreco=0;
  // ZprimeCandidate* BestZprimeCandidate = evt.get(h_BestZprimeCandidateChi2);
  if (BestZprimeCandidate->lepton().charge()>0) {
    dyreco = TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()); 
  } else {
    dyreco = TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()); 
  }
  evt.set(h_dyreco,dyreco);

  if(pt_hadTop > pt_hadTop_thresh && dyreco >0){
    evt.set(h_Sigma_phi_1,sphi);
  }
  if(pt_hadTop > pt_hadTop_thresh && dyreco <0){
    evt.set(h_Sigma_phi_2,sphi);

  }

  if(pt_hadTop < pt_hadTop_thresh && dphi >0){
    evt.set(h_dyreco_1,dyreco);
  }
  if(pt_hadTop < pt_hadTop_thresh && dphi <0){
    evt.set(h_dyreco_2,dyreco);
  }
// random generator with eta-dependend random seed for k-fold validation
  double leading_jet_phi = Ak4jets->at(0).v4().phi();
  std::srand((int)(1000 * leading_jet_phi));
  evt.set(h_uniform_random, ((double) rand()) / RAND_MAX);

  return true;
}


////////////////////////////////////////////////////////////////
/////Saving EFT variables to TTree at different DNN stages//////
////////////////////////////////////////////////////////////////



Variables_EFT_SR::Variables_EFT_SR(uhh2::Context& ctx, TString mode): mode_(mode){
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");
  h_eventweight_SR = ctx.declare_event_output<float> ("eventweight");

  h_dyreco_SR = ctx.declare_event_output<float>("dyreco_SR");
  h_Sigma_phi_SR = ctx.declare_event_output<float>("Sigma_phi_SR");
  h_Delta_phi_SR = ctx.declare_event_output<float>("Delta_phi_SR");

  h_Delta_phi_1_SR = ctx.declare_event_output<float>("Delta_phi_1_SR");
  h_Delta_phi_2_SR = ctx.declare_event_output<float>("Delta_phi_2_SR");
  
  h_dyreco_1_SR = ctx.declare_event_output<float>("dyreco_1_SR");
  h_dyreco_1_SR_0_500 = ctx.declare_event_output<float>("dyreco_1_SR_0_500");
  h_dyreco_1_SR_500_750 = ctx.declare_event_output<float>("dyreco_1_SR_500_750");
  h_dyreco_1_SR_750_1000 = ctx.declare_event_output<float>("dyreco_1_SR_750_1000");
  h_dyreco_1_SR_1000_1500 = ctx.declare_event_output<float>("dyreco_1_SR_1000_1500");
  h_dyreco_1_SR_1500_Inf = ctx.declare_event_output<float>("dyreco_1_SR_1500_Inf");
  h_dyreco_1_SR_0_700 = ctx.declare_event_output<float>("dyreco_1_SR_0_700");
  h_dyreco_1_SR_700_900 = ctx.declare_event_output<float>("dyreco_1_SR_700_900");
  h_dyreco_1_SR_900_Inf = ctx.declare_event_output<float>("dyreco_1_SR_900_Inf");

  h_dyreco_2_SR = ctx.declare_event_output<float>("dyreco_2_SR");
  h_dyreco_2_SR_0_500 = ctx.declare_event_output<float>("dyreco_2_SR_0_500");
  h_dyreco_2_SR_500_750 = ctx.declare_event_output<float>("dyreco_2_SR_500_750");
  h_dyreco_2_SR_750_1000 = ctx.declare_event_output<float>("dyreco_2_SR_750_1000");
  h_dyreco_2_SR_1000_1500 = ctx.declare_event_output<float>("dyreco_2_SR_1000_1500");
  h_dyreco_2_SR_1500_Inf = ctx.declare_event_output<float>("dyreco_2_SR_1500_Inf");
  h_dyreco_2_SR_0_700 = ctx.declare_event_output<float>("dyreco_2_SR_0_700");
  h_dyreco_2_SR_700_900 = ctx.declare_event_output<float>("dyreco_2_SR_700_900");
  h_dyreco_2_SR_900_Inf = ctx.declare_event_output<float>("dyreco_2_SR_900_Inf");



  h_Sigma_phi_1_SR=ctx.declare_event_output<float>("Sigma_phi_1_SR");
  h_Sigma_phi_1_SR_0_500=ctx.declare_event_output<float>("Sigma_phi_1_SR_0_500");
  h_Sigma_phi_1_SR_500_750=ctx.declare_event_output<float>("Sigma_phi_1_SR_500_750");
  h_Sigma_phi_1_SR_750_1000=ctx.declare_event_output<float>("Sigma_phi_1_SR_750_1000");
  h_Sigma_phi_1_SR_1000_1500=ctx.declare_event_output<float>("Sigma_phi_1_SR_1000_1500");
  h_Sigma_phi_1_SR_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_1_SR_1000_1500");
  h_Sigma_phi_1_SR_0_700=ctx.declare_event_output<float>("Sigma_phi_1_SR_0_700");
  h_Sigma_phi_1_SR_700_900=ctx.declare_event_output<float>("Sigma_phi_1_SR_700_900");
  h_Sigma_phi_1_SR_900_Inf=ctx.declare_event_output<float>("Sigma_phi_1_SR_900_Inf");

  h_Sigma_phi_2_SR=ctx.declare_event_output<float>("Sigma_phi_2_SR");
  h_Sigma_phi_2_SR_0_500=ctx.declare_event_output<float>("Sigma_phi_2_SR_0_500");
  h_Sigma_phi_2_SR_500_750=ctx.declare_event_output<float>("Sigma_phi_2_SR_500_750");
  h_Sigma_phi_2_SR_750_1000=ctx.declare_event_output<float>("Sigma_phi_2_SR_750_1000");
  h_Sigma_phi_2_SR_1000_1500=ctx.declare_event_output<float>("Sigma_phi_2_SR_1000_1500");
  h_Sigma_phi_2_SR_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_2_SR_1000_1500");
  h_Sigma_phi_2_SR_0_700=ctx.declare_event_output<float>("Sigma_phi_2_SR_0_700");
  h_Sigma_phi_2_SR_700_900=ctx.declare_event_output<float>("Sigma_phi_2_SR_700_900");
  h_Sigma_phi_2_SR_900_Inf=ctx.declare_event_output<float>("Sigma_phi_2_SR_900_Inf");

}

bool Variables_EFT_SR::process(uhh2::Event& evt){

  double weight = evt.weight;
  evt.set(h_eventweight_SR, -10);
  evt.set(h_eventweight_SR, weight);

  // bool is_zprime_reconstructed_chi2 = evt.get(h_is_zprime_reconstructed_chi2);
  ZprimeCandidate* BestZprimeCandidate = evt.get(h_BestZprimeCandidateChi2);
  bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction(); // Reconstruction process id
  vector <Jet> AK4CHSjets_matched = evt.get(h_CHSjets_matched);  
  //                 // AK4Puppijets that have been matched to CHSjets
 
  // vector <TopJet> TopTaggedJets = evt.get(h_AK8TopTags);                     // AK8Puppi jets TopTagged by DeepAK8TopTagger
  vector <float> jets_hadronic_bscores;                                            // bScores vector for resolved hadronic jets
  float pt_hadTop_thresh = 150;                                                    // Define cut-variable as pt of hadTop for low/high regions                                                   // medium WP for UL18 DeepJet

  // EFT Ac and spin correlation variables:
  // Plot pt of hadronic Top jet
  float pt_hadTop = BestZprimeCandidate->top_hadronic_v4().pt();
  
  float bscore_max = -2;
  if(!is_toptag_reconstruction){
      // Loop over resolved hadronic jets to find their bscore via CHS jets
    for(unsigned int i=0; i<BestZprimeCandidate->jets_hadronic().size(); i++){
      double deltaR_min = 99;
      // Match resolved hadronic jets to CHS jets (which have bscores)
      for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
        double deltaR_CHS = deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(j));
        if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;
        }
      // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
      for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
        if(deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
        jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet());
        } // Using DeepJet btag score
    }
    // Loop over bScores-vector to extract highest bscor
    for(unsigned int i=0; i<jets_hadronic_bscores.size(); i++){
      float bscore = jets_hadronic_bscores.at(i);
      if(bscore > bscore_max) bscore_max = bscore;
    }
    //is not top tag
  }
  if(is_toptag_reconstruction){
    // Loop over hadronic top's subjets to extract highest bscore
  for(unsigned int i=0; i < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); i++){
    float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet(); // Using DeepJet btag score
    if(bscore > bscore_max) bscore_max = bscore;
  }

}
  TLorentzVector had_top_b(0, 0, 0, 0);

// Resolved topology
  if(!is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK4-jet with highest bscore
    for(unsigned int i=0; i< BestZprimeCandidate->jets_hadronic().size(); i++){
      float bscore = jets_hadronic_bscores.at(i);
      if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->jets_hadronic().at(i).pt(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).eta(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).phi(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).energy());
    }
  }
  // Merged topology
  if(is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK8-subjet with highest bscore
    for(unsigned int j=0; j < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); j++){
      float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
      if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).pt(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).eta(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).phi(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).energy());
    }
  }
  TLorentzVector lep_top_lep(0, 0, 0, 0);
  LorentzVector lep = BestZprimeCandidate->lepton().v4();
  lep_top_lep.SetPtEtaPhiE(lep.pt(), lep.eta(), lep.phi(), lep.E());
    // Define 4vectors of top quarks
  TLorentzVector PosTop(0, 0, 0, 0);
  TLorentzVector NegTop(0, 0, 0, 0);

      // POSITIVE LEPTON CONFIGURATION => Positive charged lepton has Positive Top mother
  if(BestZprimeCandidate->lepton().charge() > 0){
    // Define ttbar system
    PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                        BestZprimeCandidate->top_leptonic_v4().eta(), 
                        BestZprimeCandidate->top_leptonic_v4().phi(), 
                        BestZprimeCandidate->top_leptonic_v4().energy());
    NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                        BestZprimeCandidate->top_hadronic_v4().eta(), 
                        BestZprimeCandidate->top_hadronic_v4().phi(), 
                        BestZprimeCandidate->top_hadronic_v4().energy());
  }
  else if (BestZprimeCandidate->lepton().charge() < 0){
    PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                        BestZprimeCandidate->top_hadronic_v4().eta(), 
                        BestZprimeCandidate->top_hadronic_v4().phi(), 
                        BestZprimeCandidate->top_hadronic_v4().energy());
    NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                        BestZprimeCandidate->top_leptonic_v4().eta(), 
                        BestZprimeCandidate->top_leptonic_v4().phi(), 
                        BestZprimeCandidate->top_leptonic_v4().energy());
  }
      
  TLorentzVector ttbar = PosTop + NegTop;
  TLorentzVector lep_top_lep_CoM = lep_top_lep;
  // Boost into ttbar CoM-Frame <<<-------//
    // Boost into ttbar CoM-Frame <<<------- Step 1//
    lep_top_lep_CoM.Boost(-1.*ttbar.BoostVector());
    TLorentzVector had_top_b_CoM = had_top_b;
    had_top_b_CoM.Boost(-ttbar.BoostVector());
    TLorentzVector PosTop_CoM = PosTop;
    PosTop_CoM.Boost(-ttbar.BoostVector());
    TLorentzVector NegTop_CoM = NegTop;
    NegTop_CoM.Boost(-ttbar.BoostVector());
    // Beam unit vector in COM frame Step 2
    TVector3 beam_axis(0,0,1);
    // Calculating top scattering angle for PosTop only
    double cos_PosTop_beam = PosTop_CoM.Vect().Unit().Dot(beam_axis);
    double sin_PosTop_beam = sqrt(1 - cos_PosTop_beam*cos_PosTop_beam);
  
    // The sign of cos_PosTop_beam to account for Bose symmetry
    double sign_cos_PosTop_beam = (cos_PosTop_beam > 0.) ? 1. : -1.;
    // // The sign based on PosTop and NegTop's rapidity
    // double sign_rapidity = (PosTop.Rapidity() >= NegTop.Rapidity()) ? 1. : -1.;
  
    // Bernreuther basis vectors
    TVector3 kbase = PosTop_CoM.Vect().Unit();
    TVector3 rbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*(beam_axis - cos_PosTop_beam * kbase) ).Unit(); /// check with Lin about numerator
    TVector3 nbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*beam_axis.Cross(kbase) ).Unit();
  
    
    TLorentzVector lep_top_lep_Rest = lep_top_lep_CoM;
    TLorentzVector had_top_b_Rest = had_top_b_CoM;
    TLorentzVector PosTop_Rest = PosTop_CoM;
    TLorentzVector NegTop_Rest = NegTop_CoM;
    
    if(BestZprimeCandidate->lepton().charge() > 0){
      lep_top_lep_Rest.Boost(-1.*PosTop_CoM.BoostVector()); // lepton has Positive Top mother
      had_top_b_Rest.Boost(-1.*NegTop_CoM.BoostVector());   // b-jet has Negative Top mother
    }
    else if (BestZprimeCandidate->lepton().charge() < 0){
      lep_top_lep_Rest.Boost(-1.*NegTop_CoM.BoostVector()); // lepton has Negative Top mother
      had_top_b_Rest.Boost(-1.*PosTop_CoM.BoostVector());   // b-jet has Positive Top mother
    }
    //Step 4, calculating the angular variables in the ttbar rest frame
    float bquark_phi = TMath::ATan2(had_top_b_Rest.Vect().Dot(nbase), had_top_b_Rest.Vect().Dot(rbase));
    float lep_phi = TMath::ATan2(lep_top_lep_Rest.Vect().Dot(nbase), lep_top_lep_Rest.Vect().Dot(rbase));
  
    float dphi_SR= 0.;
    
    float sphi_SR = lep_phi + bquark_phi;
    if(sphi_SR > TMath::Pi()) {
      sphi_SR = sphi_SR - 2.*TMath::Pi();
    }
    if(sphi_SR < -1.*TMath::Pi()) {
      sphi_SR = sphi_SR + 2.*TMath::Pi();
    }
    if(BestZprimeCandidate->lepton().charge() > 0){ // lepton is Positive Top's Decay Product
      dphi_SR = lep_phi - bquark_phi;
    }
    if(BestZprimeCandidate->lepton().charge() < 0){
      dphi_SR = bquark_phi - lep_phi;
    }
    if(dphi_SR > TMath::Pi()) {
      dphi_SR = dphi_SR - 2.*TMath::Pi();
    }
    if(dphi_SR < -1.*TMath::Pi()) {
      dphi_SR = dphi_SR + 2.*TMath::Pi();
    }
  
  evt.set(h_Sigma_phi_SR,-10);
  evt.set(h_Delta_phi_SR,-10);
  evt.set(h_Sigma_phi_1_SR,-10);
  evt.set(h_Sigma_phi_2_SR,-10);
  evt.set(h_dyreco_1_SR,-10);
  evt.set(h_dyreco_2_SR,-10);
  evt.set(h_Delta_phi_1_SR,-10);
  evt.set(h_Delta_phi_2_SR,-10);
  evt.set(h_Sigma_phi_1_SR_0_500,-10);
  evt.set(h_Sigma_phi_1_SR_500_750,-10);
  evt.set(h_Sigma_phi_1_SR_750_1000,-10);
  evt.set(h_Sigma_phi_1_SR_1000_1500,-10);
  evt.set(h_Sigma_phi_1_SR_1500_Inf,-10);
  evt.set(h_Sigma_phi_2_SR_0_500,-10);
  evt.set(h_Sigma_phi_2_SR_500_750,-10);
  evt.set(h_Sigma_phi_2_SR_750_1000,-10);
  evt.set(h_Sigma_phi_2_SR_1000_1500,-10);
  evt.set(h_Sigma_phi_2_SR_1500_Inf,-10);
  evt.set(h_dyreco_1_SR_0_500,-10);
  evt.set(h_dyreco_1_SR_500_750,-10);
  evt.set(h_dyreco_1_SR_750_1000,-10);
  evt.set(h_dyreco_1_SR_1000_1500,-10);
  evt.set(h_dyreco_1_SR_1500_Inf,-10);
  evt.set(h_dyreco_2_SR_0_500,-10);
  evt.set(h_dyreco_2_SR_500_750,-10);
  evt.set(h_dyreco_2_SR_750_1000,-10);
  evt.set(h_dyreco_2_SR_1000_1500,-10);
  evt.set(h_dyreco_2_SR_1500_Inf,-10);
  evt.set(h_dyreco_1_SR_0_700,-10);
  evt.set(h_dyreco_1_SR_700_900,-10);
  evt.set(h_dyreco_1_SR_900_Inf,-10);
  evt.set(h_dyreco_2_SR_0_700,-10);
  evt.set(h_dyreco_2_SR_700_900,-10);
  evt.set(h_dyreco_2_SR_900_Inf,-10);
  evt.set(h_Sigma_phi_1_SR_0_700,-10);
  evt.set(h_Sigma_phi_1_SR_700_900,-10);
  evt.set(h_Sigma_phi_1_SR_900_Inf,-10);
  evt.set(h_Sigma_phi_2_SR_0_700,-10);
  evt.set(h_Sigma_phi_2_SR_700_900,-10);
  evt.set(h_Sigma_phi_2_SR_900_Inf,-10);
  evt.set(h_Sigma_phi_SR,sphi_SR);
  evt.set(h_Delta_phi_SR,dphi_SR);
  evt.set(h_dyreco_SR,-10);
  
  
  float dy_reco_SR=0;
  // ZprimeCandidate* BestZprimeCandidate = evt.get(h_BestZprimeCandidateChi2);
  if (BestZprimeCandidate->lepton().charge()>0) {
    dy_reco_SR = TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()); 
  } else {
    dy_reco_SR = TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()); 
  }
  evt.set(h_dyreco_SR,dy_reco_SR);
  float Mass_tt = BestZprimeCandidate->Zprime_v4().M();

  if(pt_hadTop > pt_hadTop_thresh && dy_reco_SR >0){
    evt.set(h_Sigma_phi_1_SR,sphi_SR);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_Sigma_phi_1_SR_0_500,sphi_SR);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_Sigma_phi_1_SR_500_750,sphi_SR);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_1_SR_750_1000,sphi_SR);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_Sigma_phi_1_SR_1000_1500,sphi_SR);
    }
    if(Mass_tt>=1500){
      evt.set(h_Sigma_phi_1_SR_1500_Inf,sphi_SR);
    }
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_Sigma_phi_1_SR_0_700,sphi_SR);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_1_SR_700_900,sphi_SR);
    }
    if(Mass_tt>=1000){
      evt.set(h_Sigma_phi_1_SR_900_Inf,sphi_SR);
    }
  }
  if(pt_hadTop > pt_hadTop_thresh && dy_reco_SR <0){
    evt.set(h_Sigma_phi_2_SR,sphi_SR);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_Sigma_phi_2_SR_0_500,sphi_SR);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_Sigma_phi_2_SR_500_750,sphi_SR);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_2_SR_750_1000,sphi_SR);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_Sigma_phi_2_SR_1000_1500,sphi_SR);
    }
    if(Mass_tt>=1500){
      evt.set(h_Sigma_phi_2_SR_1500_Inf,sphi_SR);
    }
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_Sigma_phi_2_SR_0_700,sphi_SR);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_2_SR_700_900,sphi_SR);
    }
    if(Mass_tt>=1000){
      evt.set(h_Sigma_phi_2_SR_900_Inf,sphi_SR);
    }
  }

  if(pt_hadTop < pt_hadTop_thresh && dphi_SR >0){
    evt.set(h_dyreco_1_SR,dy_reco_SR);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_dyreco_1_SR_0_500,dy_reco_SR);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_dyreco_1_SR_500_750,dy_reco_SR);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_dyreco_1_SR_750_1000,dy_reco_SR);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_dyreco_1_SR_1000_1500,dy_reco_SR);
    }
    if(Mass_tt>=1500){
      evt.set(h_dyreco_1_SR_1500_Inf,dy_reco_SR);
    }
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_dyreco_1_SR_0_700,dy_reco_SR);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_dyreco_1_SR_700_900,dy_reco_SR);
    }
    if(Mass_tt>=1000){
      evt.set(h_dyreco_1_SR_900_Inf,dy_reco_SR);
    }
  }
  if(pt_hadTop < pt_hadTop_thresh && dphi_SR <0){
    evt.set(h_dyreco_2_SR,dy_reco_SR);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_dyreco_2_SR_0_500,dy_reco_SR);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_dyreco_2_SR_500_750,dy_reco_SR);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_dyreco_2_SR_750_1000,dy_reco_SR);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_dyreco_2_SR_1000_1500,dy_reco_SR);
    }
    if(Mass_tt>=1500){
      evt.set(h_dyreco_2_SR_1500_Inf,dy_reco_SR);
    }
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_dyreco_2_SR_0_700,dy_reco_SR);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_dyreco_2_SR_700_900,dy_reco_SR);
    }
    if(Mass_tt>=1000){
      evt.set(h_dyreco_2_SR_900_Inf,dy_reco_SR);
    }
  }
  if(pt_hadTop < pt_hadTop_thresh){
    evt.set(h_Delta_phi_1_SR,dphi_SR);
  }
  if(pt_hadTop > pt_hadTop_thresh){
    evt.set(h_Delta_phi_2_SR,dphi_SR);
  }


  return true;
}


/////CR1///////

Variables_EFT_CR1::Variables_EFT_CR1(uhh2::Context& ctx, TString mode): mode_(mode){
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");
  h_eventweight_CR1 = ctx.declare_event_output<float> ("eventweight");

  h_dyreco_CR1 = ctx.declare_event_output<float>("dyreco_CR1");
  h_Sigma_phi_CR1 = ctx.declare_event_output<float>("Sigma_phi_CR1");
  h_Delta_phi_CR1 = ctx.declare_event_output<float>("Delta_phi_CR1");

  h_Delta_phi_1_CR1 = ctx.declare_event_output<float>("Delta_phi_1_CR1");
  h_Delta_phi_2_CR1 = ctx.declare_event_output<float>("Delta_phi_2_CR1");
  
  h_dyreco_1_CR1 = ctx.declare_event_output<float>("dyreco_1_CR1");
  h_dyreco_1_CR1_0_500 = ctx.declare_event_output<float>("dyreco_1_CR1_0_500");
  h_dyreco_1_CR1_500_750 = ctx.declare_event_output<float>("dyreco_1_CR1_500_750");
  h_dyreco_1_CR1_750_1000 = ctx.declare_event_output<float>("dyreco_1_CR1_750_1000");
  h_dyreco_1_CR1_1000_1500 = ctx.declare_event_output<float>("dyreco_1_CR1_1000_1500");
  h_dyreco_1_CR1_1500_Inf = ctx.declare_event_output<float>("dyreco_1_CR1_1500_Inf");
  h_dyreco_1_CR1_0_700 = ctx.declare_event_output<float>("dyreco_1_CR1_0_700");
  h_dyreco_1_CR1_700_900 = ctx.declare_event_output<float>("dyreco_1_CR1_700_900");
  h_dyreco_1_CR1_900_Inf = ctx.declare_event_output<float>("dyreco_1_CR1_900_Inf");

  h_dyreco_2_CR1 = ctx.declare_event_output<float>("dyreco_2_CR1");
  h_dyreco_2_CR1_0_500 = ctx.declare_event_output<float>("dyreco_2_CR1_0_500");
  h_dyreco_2_CR1_500_750 = ctx.declare_event_output<float>("dyreco_2_CR1_500_750");
  h_dyreco_2_CR1_750_1000 = ctx.declare_event_output<float>("dyreco_2_CR1_750_1000");
  h_dyreco_2_CR1_1000_1500 = ctx.declare_event_output<float>("dyreco_2_CR1_1000_1500");
  h_dyreco_2_CR1_1500_Inf = ctx.declare_event_output<float>("dyreco_2_CR1_1500_Inf");
  h_dyreco_2_CR1_0_700 = ctx.declare_event_output<float>("dyreco_2_CR1_0_700");
  h_dyreco_2_CR1_700_900 = ctx.declare_event_output<float>("dyreco_2_CR1_700_900");
  h_dyreco_2_CR1_900_Inf = ctx.declare_event_output<float>("dyreco_2_CR1_900_Inf");


  h_Sigma_phi_1_CR1=ctx.declare_event_output<float>("Sigma_phi_1_CR1");
  h_Sigma_phi_1_CR1_0_500=ctx.declare_event_output<float>("Sigma_phi_1_CR1_0_500");
  h_Sigma_phi_1_CR1_500_750=ctx.declare_event_output<float>("Sigma_phi_1_CR1_500_750");
  h_Sigma_phi_1_CR1_750_1000=ctx.declare_event_output<float>("Sigma_phi_1_CR1_750_1000");
  h_Sigma_phi_1_CR1_1000_1500=ctx.declare_event_output<float>("Sigma_phi_1_CR1_1000_1500");
  h_Sigma_phi_1_CR1_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_1_CR1_1000_1500");
  h_Sigma_phi_1_CR1_0_700=ctx.declare_event_output<float>("Sigma_phi_1_CR1_0_700");
  h_Sigma_phi_1_CR1_700_900=ctx.declare_event_output<float>("Sigma_phi_1_CR1_700_900");
  h_Sigma_phi_1_CR1_900_Inf=ctx.declare_event_output<float>("Sigma_phi_1_CR1_900_Inf");

  h_Sigma_phi_2_CR1=ctx.declare_event_output<float>("Sigma_phi_2_CR1");
  h_Sigma_phi_2_CR1_0_500=ctx.declare_event_output<float>("Sigma_phi_2_CR1_0_500");
  h_Sigma_phi_2_CR1_500_750=ctx.declare_event_output<float>("Sigma_phi_2_CR1_500_750");
  h_Sigma_phi_2_CR1_750_1000=ctx.declare_event_output<float>("Sigma_phi_2_CR1_750_1000");
  h_Sigma_phi_2_CR1_1000_1500=ctx.declare_event_output<float>("Sigma_phi_2_CR1_1000_1500");
  h_Sigma_phi_2_CR1_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_2_CR1_1000_1500");
  h_Sigma_phi_2_CR1_0_700=ctx.declare_event_output<float>("Sigma_phi_2_CR1_0_700");
  h_Sigma_phi_2_CR1_700_900=ctx.declare_event_output<float>("Sigma_phi_2_CR1_700_900");
  h_Sigma_phi_2_CR1_900_Inf=ctx.declare_event_output<float>("Sigma_phi_2_CR1_900_Inf");
}

bool Variables_EFT_CR1::process(uhh2::Event& evt){

  double weight = evt.weight;
  evt.set(h_eventweight_CR1, -10);
  evt.set(h_eventweight_CR1, weight);

  // bool is_zprime_reconstructed_chi2 = evt.get(h_is_zprime_reconstructed_chi2);
  ZprimeCandidate* BestZprimeCandidate = evt.get(h_BestZprimeCandidateChi2);
  bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction(); // Reconstruction process id
  vector <Jet> AK4CHSjets_matched = evt.get(h_CHSjets_matched);  
  //                 // AK4Puppijets that have been matched to CHSjets
 
  // vector <TopJet> TopTaggedJets = evt.get(h_AK8TopTags);                     // AK8Puppi jets TopTagged by DeepAK8TopTagger
  vector <float> jets_hadronic_bscores;                                            // bScores vector for resolved hadronic jets
  float pt_hadTop_thresh = 150;                                                    // Define cut-variable as pt of hadTop for low/high regions                                                   // medium WP for UL18 DeepJet
  // EFT Ac and spin correlation variables:
  // Plot pt of hadronic Top jet
  float pt_hadTop = BestZprimeCandidate->top_hadronic_v4().pt();
  
  float bscore_max = -2;
  if(!is_toptag_reconstruction){
      // Loop over resolved hadronic jets to find their bscore via CHS jets
    for(unsigned int i=0; i<BestZprimeCandidate->jets_hadronic().size(); i++){
      double deltaR_min = 99;
      // Match resolved hadronic jets to CHS jets (which have bscores)
      for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
        double deltaR_CHS = deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(j));
        if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;
        }
      // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
      for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
        if(deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
        jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet());
        } // Using DeepJet btag score
    }
    // Loop over bScores-vector to extract highest bscor
    for(unsigned int i=0; i<jets_hadronic_bscores.size(); i++){
      float bscore = jets_hadronic_bscores.at(i);
      if(bscore > bscore_max) bscore_max = bscore;
    }
    //is not top tag
  }
  if(is_toptag_reconstruction){
    // Loop over hadronic top's subjets to extract highest bscore
  for(unsigned int i=0; i < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); i++){
    float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet(); // Using DeepJet btag score
    if(bscore > bscore_max) bscore_max = bscore;
  }

}
  TLorentzVector had_top_b(0, 0, 0, 0);

// Resolved topology
  if(!is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK4-jet with highest bscore
    for(unsigned int i=0; i< BestZprimeCandidate->jets_hadronic().size(); i++){
      float bscore = jets_hadronic_bscores.at(i);
      if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->jets_hadronic().at(i).pt(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).eta(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).phi(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).energy());
    }
  }
  // Merged topology
  if(is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK8-subjet with highest bscore
    for(unsigned int j=0; j < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); j++){
      float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
      if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).pt(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).eta(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).phi(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).energy());
    }
  }
  TLorentzVector lep_top_lep(0, 0, 0, 0);
  LorentzVector lep = BestZprimeCandidate->lepton().v4();
  lep_top_lep.SetPtEtaPhiE(lep.pt(), lep.eta(), lep.phi(), lep.E());
    // Define 4vectors of top quarks
  TLorentzVector PosTop(0, 0, 0, 0);
  TLorentzVector NegTop(0, 0, 0, 0);

      // POSITIVE LEPTON CONFIGURATION => Positive charged lepton has Positive Top mother
  if(BestZprimeCandidate->lepton().charge() > 0){
    // Define ttbar system
    PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                        BestZprimeCandidate->top_leptonic_v4().eta(), 
                        BestZprimeCandidate->top_leptonic_v4().phi(), 
                        BestZprimeCandidate->top_leptonic_v4().energy());
    NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                        BestZprimeCandidate->top_hadronic_v4().eta(), 
                        BestZprimeCandidate->top_hadronic_v4().phi(), 
                        BestZprimeCandidate->top_hadronic_v4().energy());
  }
  else if (BestZprimeCandidate->lepton().charge() < 0){
    PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                        BestZprimeCandidate->top_hadronic_v4().eta(), 
                        BestZprimeCandidate->top_hadronic_v4().phi(), 
                        BestZprimeCandidate->top_hadronic_v4().energy());
    NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                        BestZprimeCandidate->top_leptonic_v4().eta(), 
                        BestZprimeCandidate->top_leptonic_v4().phi(), 
                        BestZprimeCandidate->top_leptonic_v4().energy());
  }
      
  TLorentzVector ttbar = PosTop + NegTop;
  TLorentzVector lep_top_lep_CoM = lep_top_lep;
  // Boost into ttbar CoM-Frame <<<------- Step 1//
  lep_top_lep_CoM.Boost(-1.*ttbar.BoostVector());
  TLorentzVector had_top_b_CoM = had_top_b;
  had_top_b_CoM.Boost(-ttbar.BoostVector());
  TLorentzVector PosTop_CoM = PosTop;
  PosTop_CoM.Boost(-ttbar.BoostVector());
  TLorentzVector NegTop_CoM = NegTop;
  NegTop_CoM.Boost(-ttbar.BoostVector());
  // Beam unit vector in COM frame Step 2
  TVector3 beam_axis(0,0,1);
  // Calculating top scattering angle for PosTop only
  double cos_PosTop_beam = PosTop_CoM.Vect().Unit().Dot(beam_axis);
  double sin_PosTop_beam = sqrt(1 - cos_PosTop_beam*cos_PosTop_beam);

  // The sign of cos_PosTop_beam to account for Bose symmetry
  double sign_cos_PosTop_beam = (cos_PosTop_beam > 0.) ? 1. : -1.;
  // // The sign based on PosTop and NegTop's rapidity
  // double sign_rapidity = (PosTop.Rapidity() >= NegTop.Rapidity()) ? 1. : -1.;

  // Bernreuther basis vectors
  TVector3 kbase = PosTop_CoM.Vect().Unit();
  TVector3 rbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*(beam_axis - cos_PosTop_beam * kbase) ).Unit(); /// check with Lin about numerator
  TVector3 nbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*beam_axis.Cross(kbase) ).Unit();

  
  TLorentzVector lep_top_lep_Rest = lep_top_lep_CoM;
  TLorentzVector had_top_b_Rest = had_top_b_CoM;
  TLorentzVector PosTop_Rest = PosTop_CoM;
  TLorentzVector NegTop_Rest = NegTop_CoM;
  
  if(BestZprimeCandidate->lepton().charge() > 0){
    lep_top_lep_Rest.Boost(-1.*PosTop_CoM.BoostVector()); // lepton has Positive Top mother
    had_top_b_Rest.Boost(-1.*NegTop_CoM.BoostVector());   // b-jet has Negative Top mother
  }
  else if (BestZprimeCandidate->lepton().charge() < 0){
    lep_top_lep_Rest.Boost(-1.*NegTop_CoM.BoostVector()); // lepton has Negative Top mother
    had_top_b_Rest.Boost(-1.*PosTop_CoM.BoostVector());   // b-jet has Positive Top mother
  }
  //Step 4, calculating the angular variables in the ttbar rest frame
  float bquark_phi = TMath::ATan2(had_top_b_Rest.Vect().Dot(nbase), had_top_b_Rest.Vect().Dot(rbase));
  float lep_phi = TMath::ATan2(lep_top_lep_Rest.Vect().Dot(nbase), lep_top_lep_Rest.Vect().Dot(rbase));

  float dphi_CR1= 0.;
  
  float sphi_CR1 = lep_phi + bquark_phi;
  if(sphi_CR1 > TMath::Pi()) {
    sphi_CR1 = sphi_CR1 - 2.*TMath::Pi();
  }
  if(sphi_CR1 < -1.*TMath::Pi()) {
    sphi_CR1 = sphi_CR1 + 2.*TMath::Pi();
  }
  if(BestZprimeCandidate->lepton().charge() > 0){ // lepton is Positive Top's Decay Product
    dphi_CR1 = lep_phi - bquark_phi;
  }
  if(BestZprimeCandidate->lepton().charge() < 0){
    dphi_CR1 = bquark_phi - lep_phi;
  }
  if(dphi_CR1 > TMath::Pi()){
    dphi_CR1 = dphi_CR1 - 2.*TMath::Pi();
  }
  if(dphi_CR1 < -1.*TMath::Pi()) {
    dphi_CR1 = dphi_CR1 + 2.*TMath::Pi();
  }

  evt.set(h_Sigma_phi_CR1,-10);
  evt.set(h_Delta_phi_CR1,-10);
  evt.set(h_Sigma_phi_1_CR1,-10);
  evt.set(h_Sigma_phi_2_CR1,-10);
  evt.set(h_dyreco_1_CR1,-10);
  evt.set(h_dyreco_2_CR1,-10);
  evt.set(h_Delta_phi_1_CR1,-10);
  evt.set(h_Delta_phi_2_CR1,-10);
  evt.set(h_Sigma_phi_1_CR1_0_500,-10);
  evt.set(h_Sigma_phi_1_CR1_500_750,-10);
  evt.set(h_Sigma_phi_1_CR1_750_1000,-10);
  evt.set(h_Sigma_phi_1_CR1_1000_1500,-10);
  evt.set(h_Sigma_phi_1_CR1_1500_Inf,-10);
  evt.set(h_Sigma_phi_2_CR1_0_500,-10);
  evt.set(h_Sigma_phi_2_CR1_500_750,-10);
  evt.set(h_Sigma_phi_2_CR1_750_1000,-10);
  evt.set(h_Sigma_phi_2_CR1_1000_1500,-10);
  evt.set(h_Sigma_phi_2_CR1_1500_Inf,-10);
  evt.set(h_dyreco_1_CR1_0_500,-10);
  evt.set(h_dyreco_1_CR1_500_750,-10);
  evt.set(h_dyreco_1_CR1_750_1000,-10);
  evt.set(h_dyreco_1_CR1_1000_1500,-10);
  evt.set(h_dyreco_1_CR1_1500_Inf,-10);
  evt.set(h_dyreco_2_CR1_0_500,-10);
  evt.set(h_dyreco_2_CR1_500_750,-10);
  evt.set(h_dyreco_2_CR1_750_1000,-10);
  evt.set(h_dyreco_2_CR1_1000_1500,-10);
  evt.set(h_dyreco_2_CR1_1500_Inf,-10);
  evt.set(h_dyreco_1_CR1_0_700,-10);
  evt.set(h_dyreco_1_CR1_700_900,-10);
  evt.set(h_dyreco_1_CR1_900_Inf,-10);
  evt.set(h_dyreco_2_CR1_0_700,-10);
  evt.set(h_dyreco_2_CR1_700_900,-10);
  evt.set(h_dyreco_2_CR1_900_Inf,-10);
  evt.set(h_Sigma_phi_1_CR1_0_700,-10);
  evt.set(h_Sigma_phi_1_CR1_700_900,-10);
  evt.set(h_Sigma_phi_1_CR1_900_Inf,-10);
  evt.set(h_Sigma_phi_2_CR1_0_700,-10);
  evt.set(h_Sigma_phi_2_CR1_700_900,-10);
  evt.set(h_Sigma_phi_2_CR1_900_Inf,-10);
  evt.set(h_Sigma_phi_CR1,sphi_CR1);
  evt.set(h_Delta_phi_CR1,dphi_CR1);
  evt.set(h_dyreco_CR1,-10);
  
  
  float dy_reco_CR1=0;
  // ZprimeCandidate* BestZprimeCandidate = evt.get(h_BestZprimeCandidateChi2);
  if (BestZprimeCandidate->lepton().charge()>0) {
    dy_reco_CR1 = TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()); 
  } else {
    dy_reco_CR1 = TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()); 
  }
  evt.set(h_dyreco_CR1,dy_reco_CR1);
  float Mass_tt = BestZprimeCandidate->Zprime_v4().M();

  if(pt_hadTop > pt_hadTop_thresh && dy_reco_CR1 >0){
    evt.set(h_Sigma_phi_1_CR1,sphi_CR1);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_Sigma_phi_1_CR1_0_500,sphi_CR1);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_Sigma_phi_1_CR1_500_750,sphi_CR1);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_1_CR1_750_1000,sphi_CR1);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_Sigma_phi_1_CR1_1000_1500,sphi_CR1);
    }
    if(Mass_tt>=1500){
      evt.set(h_Sigma_phi_1_CR1_1500_Inf,sphi_CR1);
    }    
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_Sigma_phi_1_CR1_0_700,sphi_CR1);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_1_CR1_700_900,sphi_CR1);
    }
    if(Mass_tt>=1000){
      evt.set(h_Sigma_phi_1_CR1_900_Inf,sphi_CR1);
    }
  }
  if(pt_hadTop > pt_hadTop_thresh && dy_reco_CR1 <0){
    evt.set(h_Sigma_phi_2_CR1,sphi_CR1);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_Sigma_phi_2_CR1_0_500,sphi_CR1);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_Sigma_phi_2_CR1_500_750,sphi_CR1);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_2_CR1_750_1000,sphi_CR1);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_Sigma_phi_2_CR1_1000_1500,sphi_CR1);
    }
    if(Mass_tt>=1500){
      evt.set(h_Sigma_phi_2_CR1_1500_Inf,sphi_CR1);
    }  
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_Sigma_phi_2_CR1_0_700,sphi_CR1);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_2_CR1_700_900,sphi_CR1);
    }
    if(Mass_tt>=1000){
      evt.set(h_Sigma_phi_2_CR1_900_Inf,sphi_CR1);
    }
  }

  if(pt_hadTop < pt_hadTop_thresh && dphi_CR1 >0){
    evt.set(h_dyreco_1_CR1,dy_reco_CR1);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_dyreco_1_CR1_0_500,dy_reco_CR1);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_dyreco_1_CR1_500_750,dy_reco_CR1);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_dyreco_1_CR1_750_1000,dy_reco_CR1);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_dyreco_1_CR1_1000_1500,dy_reco_CR1);
    }
    if(Mass_tt>=1500){
      evt.set(h_dyreco_1_CR1_1500_Inf,dy_reco_CR1);
    }  
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_dyreco_1_CR1_0_700,dy_reco_CR1);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_dyreco_1_CR1_700_900,dy_reco_CR1);
    }
    if(Mass_tt>=1000){
      evt.set(h_dyreco_1_CR1_900_Inf,dy_reco_CR1);
    }
  }
  if(pt_hadTop < pt_hadTop_thresh && dphi_CR1 <0){
    evt.set(h_dyreco_2_CR1,dy_reco_CR1);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_dyreco_2_CR1_0_500,dy_reco_CR1);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_dyreco_2_CR1_500_750,dy_reco_CR1);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_dyreco_2_CR1_750_1000,dy_reco_CR1);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_dyreco_2_CR1_1000_1500,dy_reco_CR1);
    }
    if(Mass_tt>=1500){
      evt.set(h_dyreco_2_CR1_1500_Inf,dy_reco_CR1);
    }  
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_dyreco_2_CR1_0_700,dy_reco_CR1);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_dyreco_2_CR1_700_900,dy_reco_CR1);
    }
    if(Mass_tt>=1000){
      evt.set(h_dyreco_2_CR1_900_Inf,dy_reco_CR1);
    }
  }
  if(pt_hadTop < pt_hadTop_thresh){
    evt.set(h_Delta_phi_1_CR1,dphi_CR1);
  }
  if(pt_hadTop > pt_hadTop_thresh){
    evt.set(h_Delta_phi_2_CR1,dphi_CR1);
  }


  return true;
}


// // ////CR2//////


Variables_EFT_CR2::Variables_EFT_CR2(uhh2::Context& ctx, TString mode): mode_(mode){
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");
  h_eventweight_CR2 = ctx.declare_event_output<float> ("eventweight");

  h_dyreco_CR2 = ctx.declare_event_output<float>("dyreco_CR2");
  h_Sigma_phi_CR2 = ctx.declare_event_output<float>("Sigma_phi_CR2");
  h_Delta_phi_CR2 = ctx.declare_event_output<float>("Delta_phi_CR2");

  h_Delta_phi_1_CR2 = ctx.declare_event_output<float>("Delta_phi_1_CR2");
  h_Delta_phi_2_CR2 = ctx.declare_event_output<float>("Delta_phi_2_CR2");
  
  h_dyreco_1_CR2 = ctx.declare_event_output<float>("dyreco_1_CR2");
  h_dyreco_1_CR2_0_500 = ctx.declare_event_output<float>("dyreco_1_CR2_0_500");
  h_dyreco_1_CR2_500_750 = ctx.declare_event_output<float>("dyreco_1_CR2_500_750");
  h_dyreco_1_CR2_750_1000 = ctx.declare_event_output<float>("dyreco_1_CR2_750_1000");
  h_dyreco_1_CR2_0_500 = ctx.declare_event_output<float>("dyreco_1_CR2_0_500");
  h_dyreco_1_CR2_1000_1500 = ctx.declare_event_output<float>("dyreco_1_CR2_1000_1500");
  h_dyreco_1_CR2_1500_Inf = ctx.declare_event_output<float>("dyreco_1_CR2_1500_Inf");
  h_dyreco_1_CR2_0_700 = ctx.declare_event_output<float>("dyreco_1_CR2_0_500");
  h_dyreco_1_CR2_700_900 = ctx.declare_event_output<float>("dyreco_1_CR2_700_900");
  h_dyreco_1_CR2_900_Inf = ctx.declare_event_output<float>("dyreco_1_CR2_900_Inf");


  h_dyreco_2_CR2 = ctx.declare_event_output<float>("dyreco_2_CR2");
  h_dyreco_2_CR2_0_500 = ctx.declare_event_output<float>("dyreco_2_CR2_0_500");
  h_dyreco_2_CR2_500_750 = ctx.declare_event_output<float>("dyreco_2_CR2_500_750");
  h_dyreco_2_CR2_750_1000 = ctx.declare_event_output<float>("dyreco_2_CR2_750_1000");
  h_dyreco_2_CR2_1000_1500 = ctx.declare_event_output<float>("dyreco_2_CR2_1000_1500");
  h_dyreco_2_CR2_1500_Inf = ctx.declare_event_output<float>("dyreco_2_CR2_1500_Inf");
  h_dyreco_2_CR2_0_700 = ctx.declare_event_output<float>("dyreco_1_CR2_0_500");
  h_dyreco_2_CR2_700_900 = ctx.declare_event_output<float>("dyreco_2_CR2_700_900");
  h_dyreco_2_CR2_900_Inf = ctx.declare_event_output<float>("dyreco_2_CR2_900_Inf");



  h_Sigma_phi_1_CR2=ctx.declare_event_output<float>("Sigma_phi_1_CR2");
  h_Sigma_phi_1_CR2_0_500=ctx.declare_event_output<float>("Sigma_phi_1_CR2_0_500");
  h_Sigma_phi_1_CR2_500_750=ctx.declare_event_output<float>("Sigma_phi_1_CR2_500_750");
  h_Sigma_phi_1_CR2_750_1000=ctx.declare_event_output<float>("Sigma_phi_1_CR2_750_1000");
  h_Sigma_phi_1_CR2_1000_1500=ctx.declare_event_output<float>("Sigma_phi_1_CR2_1000_1500");
  h_Sigma_phi_1_CR2_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_1_CR2_1000_1500");
  h_Sigma_phi_1_CR2_0_700 = ctx.declare_event_output<float>("Sigma_phi_1_CR2_0_500");
  h_Sigma_phi_1_CR2_700_900 = ctx.declare_event_output<float>("Sigma_phi_1_CR2_700_900");
  h_Sigma_phi_1_CR2_900_Inf = ctx.declare_event_output<float>("Sigma_phi_1_CR2_900_Inf");



  h_Sigma_phi_2_CR2=ctx.declare_event_output<float>("Sigma_phi_2_CR2");
  h_Sigma_phi_2_CR2_0_500=ctx.declare_event_output<float>("Sigma_phi_2_CR2_0_500");
  h_Sigma_phi_2_CR2_500_750=ctx.declare_event_output<float>("Sigma_phi_2_CR2_500_750");
  h_Sigma_phi_2_CR2_750_1000=ctx.declare_event_output<float>("Sigma_phi_2_CR2_750_1000");
  h_Sigma_phi_2_CR2_1000_1500=ctx.declare_event_output<float>("Sigma_phi_2_CR2_1000_1500");
  h_Sigma_phi_2_CR2_1500_Inf=ctx.declare_event_output<float>("Sigma_phi_2_CR2_1000_1500");
  h_Sigma_phi_2_CR2_0_700 = ctx.declare_event_output<float>("Sigma_phi_2_CR2_0_500");
  h_Sigma_phi_2_CR2_700_900 = ctx.declare_event_output<float>("Sigma_phi_2_CR2_700_900");
  h_Sigma_phi_2_CR2_900_Inf = ctx.declare_event_output<float>("Sigma_phi_2_CR2_900_Inf");

}

bool Variables_EFT_CR2::process(uhh2::Event& evt){

  double weight = evt.weight;
  evt.set(h_eventweight_CR2, -10);
  evt.set(h_eventweight_CR2, weight);

  // bool is_zprime_reconstructed_chi2 = evt.get(h_is_zprime_reconstructed_chi2);
  ZprimeCandidate* BestZprimeCandidate = evt.get(h_BestZprimeCandidateChi2);
  bool is_toptag_reconstruction = BestZprimeCandidate->is_toptag_reconstruction(); // Reconstruction process id
  vector <Jet> AK4CHSjets_matched = evt.get(h_CHSjets_matched);  
  //                 // AK4Puppijets that have been matched to CHSjets
 
  // vector <TopJet> TopTaggedJets = evt.get(h_AK8TopTags);                     // AK8Puppi jets TopTagged by DeepAK8TopTagger
  vector <float> jets_hadronic_bscores;                                            // bScores vector for resolved hadronic jets
  float pt_hadTop_thresh = 150;                                                    // Define cut-variable as pt of hadTop for low/high regions                                                   // medium WP for UL18 DeepJet

  // EFT Ac and spin correlation variables:
  // Plot pt of hadronic Top jet
  float pt_hadTop = BestZprimeCandidate->top_hadronic_v4().pt();
  
  float bscore_max = -2;
  if(!is_toptag_reconstruction){
      // Loop over resolved hadronic jets to find their bscore via CHS jets
    for(unsigned int i=0; i<BestZprimeCandidate->jets_hadronic().size(); i++){
      double deltaR_min = 99;
      // Match resolved hadronic jets to CHS jets (which have bscores)
      for(unsigned int j=0; j<AK4CHSjets_matched.size(); j++){
        double deltaR_CHS = deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(j));
        if(deltaR_CHS < deltaR_min) deltaR_min = deltaR_CHS;
        }
      // Build bScore-vector for resolved hadronic jets whose bscore will correspond by index
      for(unsigned int k=0; k<AK4CHSjets_matched.size(); k++){
        if(deltaR(BestZprimeCandidate->jets_hadronic().at(i), AK4CHSjets_matched.at(k)) == deltaR_min) 
        jets_hadronic_bscores.emplace_back(AK4CHSjets_matched.at(k).btag_DeepJet());
        } // Using DeepJet btag score
    }
    // Loop over bScores-vector to extract highest bscor
    for(unsigned int i=0; i<jets_hadronic_bscores.size(); i++){
      float bscore = jets_hadronic_bscores.at(i);
      if(bscore > bscore_max) bscore_max = bscore;
    }
    //is not top tag
  }
  if(is_toptag_reconstruction){
    // Loop over hadronic top's subjets to extract highest bscore
  for(unsigned int i=0; i < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); i++){
    float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(i).btag_DeepJet(); // Using DeepJet btag score
    if(bscore > bscore_max) bscore_max = bscore;
  }

}
  TLorentzVector had_top_b(0, 0, 0, 0);

// Resolved topology
  if(!is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK4-jet with highest bscore
    for(unsigned int i=0; i< BestZprimeCandidate->jets_hadronic().size(); i++){
      float bscore = jets_hadronic_bscores.at(i);
      if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->jets_hadronic().at(i).pt(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).eta(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).phi(), 
                                                      BestZprimeCandidate->jets_hadronic().at(i).energy());
    }
  }
  // Merged topology
  if(is_toptag_reconstruction){ // Define hadronic b-jet as hadronic AK8-subjet with highest bscore
    for(unsigned int j=0; j < BestZprimeCandidate->tophad_topjet_ptr()->subjets().size(); j++){
      float bscore = BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).btag_DeepJet();
      if(bscore == bscore_max) had_top_b.SetPtEtaPhiE(BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).pt(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).eta(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).phi(), 
                                                      BestZprimeCandidate->tophad_topjet_ptr()->subjets().at(j).energy());
    }
  }
  TLorentzVector lep_top_lep(0, 0, 0, 0);
  LorentzVector lep = BestZprimeCandidate->lepton().v4();
  lep_top_lep.SetPtEtaPhiE(lep.pt(), lep.eta(), lep.phi(), lep.E());
    // Define 4vectors of top quarks
  TLorentzVector PosTop(0, 0, 0, 0);
  TLorentzVector NegTop(0, 0, 0, 0);

      // POSITIVE LEPTON CONFIGURATION => Positive charged lepton has Positive Top mother
  if(BestZprimeCandidate->lepton().charge() > 0){
    // Define ttbar system
    PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                        BestZprimeCandidate->top_leptonic_v4().eta(), 
                        BestZprimeCandidate->top_leptonic_v4().phi(), 
                        BestZprimeCandidate->top_leptonic_v4().energy());
    NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                        BestZprimeCandidate->top_hadronic_v4().eta(), 
                        BestZprimeCandidate->top_hadronic_v4().phi(), 
                        BestZprimeCandidate->top_hadronic_v4().energy());
  }
  else if (BestZprimeCandidate->lepton().charge() < 0){
    PosTop.SetPtEtaPhiE(BestZprimeCandidate->top_hadronic_v4().pt(), 
                        BestZprimeCandidate->top_hadronic_v4().eta(), 
                        BestZprimeCandidate->top_hadronic_v4().phi(), 
                        BestZprimeCandidate->top_hadronic_v4().energy());
    NegTop.SetPtEtaPhiE(BestZprimeCandidate->top_leptonic_v4().pt(), 
                        BestZprimeCandidate->top_leptonic_v4().eta(), 
                        BestZprimeCandidate->top_leptonic_v4().phi(), 
                        BestZprimeCandidate->top_leptonic_v4().energy());
  }
      
  TLorentzVector ttbar = PosTop + NegTop;
  TLorentzVector lep_top_lep_CoM = lep_top_lep;
    // Boost into ttbar CoM-Frame <<<------- Step 1//
    lep_top_lep_CoM.Boost(-1.*ttbar.BoostVector());
    TLorentzVector had_top_b_CoM = had_top_b;
    had_top_b_CoM.Boost(-ttbar.BoostVector());
    TLorentzVector PosTop_CoM = PosTop;
    PosTop_CoM.Boost(-ttbar.BoostVector());
    TLorentzVector NegTop_CoM = NegTop;
    NegTop_CoM.Boost(-ttbar.BoostVector());
    // Beam unit vector in COM frame Step 2
    TVector3 beam_axis(0,0,1);
    // Calculating top scattering angle for PosTop only
    double cos_PosTop_beam = PosTop_CoM.Vect().Unit().Dot(beam_axis);
    double sin_PosTop_beam = sqrt(1 - cos_PosTop_beam*cos_PosTop_beam);
  
    // The sign of cos_PosTop_beam to account for Bose symmetry
    double sign_cos_PosTop_beam = (cos_PosTop_beam > 0.) ? 1. : -1.;
    // // The sign based on PosTop and NegTop's rapidity
    // double sign_rapidity = (PosTop.Rapidity() >= NegTop.Rapidity()) ? 1. : -1.;
  
    // Bernreuther basis vectors
    TVector3 kbase = PosTop_CoM.Vect().Unit();
    TVector3 rbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*(beam_axis - cos_PosTop_beam * kbase) ).Unit(); /// check with Lin about numerator
    TVector3 nbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*beam_axis.Cross(kbase) ).Unit();
  
    
    TLorentzVector lep_top_lep_Rest = lep_top_lep_CoM;
    TLorentzVector had_top_b_Rest = had_top_b_CoM;
    TLorentzVector PosTop_Rest = PosTop_CoM;
    TLorentzVector NegTop_Rest = NegTop_CoM;
    
    if(BestZprimeCandidate->lepton().charge() > 0){
      lep_top_lep_Rest.Boost(-1.*PosTop_CoM.BoostVector()); // lepton has Positive Top mother
      had_top_b_Rest.Boost(-1.*NegTop_CoM.BoostVector());   // b-jet has Negative Top mother
    }
    else if (BestZprimeCandidate->lepton().charge() < 0){
      lep_top_lep_Rest.Boost(-1.*NegTop_CoM.BoostVector()); // lepton has Negative Top mother
      had_top_b_Rest.Boost(-1.*PosTop_CoM.BoostVector());   // b-jet has Positive Top mother
    }
    //Step 4, calculating the angular variables in the ttbar rest frame
    float bquark_phi = TMath::ATan2(had_top_b_Rest.Vect().Dot(nbase), had_top_b_Rest.Vect().Dot(rbase));
    float lep_phi = TMath::ATan2(lep_top_lep_Rest.Vect().Dot(nbase), lep_top_lep_Rest.Vect().Dot(rbase));
  
    float dphi_CR2= 0.;
    
    float sphi_CR2 = lep_phi + bquark_phi;
    if(sphi_CR2 > TMath::Pi()) {
      sphi_CR2 = sphi_CR2 - 2.*TMath::Pi();
    }
    if(sphi_CR2 < -1.*TMath::Pi()) {
      sphi_CR2 = sphi_CR2 + 2.*TMath::Pi();
    }
    if(BestZprimeCandidate->lepton().charge() > 0){ // lepton is Positive Top's Decay Product
      dphi_CR2 = lep_phi - bquark_phi;
    }
    if(BestZprimeCandidate->lepton().charge() < 0){
      dphi_CR2 = bquark_phi - lep_phi;
    }
    if(dphi_CR2 > TMath::Pi()){
      dphi_CR2 = dphi_CR2 - 2.*TMath::Pi();
    }
    if(dphi_CR2 < -1.*TMath::Pi()) {
      dphi_CR2 = dphi_CR2 + 2.*TMath::Pi();
    }
  
  
  evt.set(h_Sigma_phi_CR2,-10);
  evt.set(h_Delta_phi_CR2,-10);
  evt.set(h_Sigma_phi_1_CR2,-10);
  evt.set(h_Sigma_phi_2_CR2,-10);
  evt.set(h_dyreco_1_CR2,-10);
  evt.set(h_dyreco_2_CR2,-10);
  evt.set(h_Delta_phi_1_CR2,-10);
  evt.set(h_Delta_phi_2_CR2,-10);
  evt.set(h_Sigma_phi_1_CR2_0_500,-10);
  evt.set(h_Sigma_phi_1_CR2_500_750,-10);
  evt.set(h_Sigma_phi_1_CR2_750_1000,-10);
  evt.set(h_Sigma_phi_1_CR2_1000_1500,-10);
  evt.set(h_Sigma_phi_1_CR2_1500_Inf,-10);
  evt.set(h_Sigma_phi_2_CR2_0_500,-10);
  evt.set(h_Sigma_phi_2_CR2_500_750,-10);
  evt.set(h_Sigma_phi_2_CR2_750_1000,-10);
  evt.set(h_Sigma_phi_2_CR2_1000_1500,-10);
  evt.set(h_Sigma_phi_2_CR2_1500_Inf,-10);
  evt.set(h_dyreco_1_CR2_0_500,-10);
  evt.set(h_dyreco_1_CR2_500_750,-10);
  evt.set(h_dyreco_1_CR2_750_1000,-10);
  evt.set(h_dyreco_1_CR2_1000_1500,-10);
  evt.set(h_dyreco_1_CR2_1500_Inf,-10);
  evt.set(h_dyreco_2_CR2_0_500,-10);
  evt.set(h_dyreco_2_CR2_500_750,-10);
  evt.set(h_dyreco_2_CR2_750_1000,-10);
  evt.set(h_dyreco_2_CR2_1000_1500,-10);
  evt.set(h_dyreco_2_CR2_1500_Inf,-10);
  evt.set(h_dyreco_1_CR2_0_700,-10);
  evt.set(h_dyreco_1_CR2_700_900,-10);
  evt.set(h_dyreco_1_CR2_900_Inf,-10);
  evt.set(h_dyreco_2_CR2_0_700,-10);
  evt.set(h_dyreco_2_CR2_700_900,-10);
  evt.set(h_dyreco_2_CR2_900_Inf,-10);
  evt.set(h_Sigma_phi_1_CR2_0_700,-10);
  evt.set(h_Sigma_phi_1_CR2_700_900,-10);
  evt.set(h_Sigma_phi_1_CR2_900_Inf,-10);
  evt.set(h_Sigma_phi_2_CR2_0_700,-10);
  evt.set(h_Sigma_phi_2_CR2_700_900,-10);
  evt.set(h_Sigma_phi_2_CR2_900_Inf,-10);  
  evt.set(h_Sigma_phi_CR2,sphi_CR2);
  evt.set(h_Delta_phi_CR2,dphi_CR2);
  evt.set(h_dyreco_CR2,-10);
  float dy_reco_CR2=0;
  // ZprimeCandidate* BestZprimeCandidate = evt.get(h_BestZprimeCandidateChi2);
  if (BestZprimeCandidate->lepton().charge()>0) {
    dy_reco_CR2 = TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()); 
  } else {
    dy_reco_CR2 = TMath::Abs(BestZprimeCandidate->top_hadronic_v4().Rapidity()) - TMath::Abs(BestZprimeCandidate->top_leptonic_v4().Rapidity()); 
  }
  evt.set(h_dyreco_CR2,dy_reco_CR2);
  float Mass_tt = BestZprimeCandidate->Zprime_v4().M();

  if(pt_hadTop > pt_hadTop_thresh && dy_reco_CR2 >0){
    evt.set(h_Sigma_phi_1_CR2,sphi_CR2);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_Sigma_phi_1_CR2_0_500,sphi_CR2);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_Sigma_phi_1_CR2_500_750,sphi_CR2);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_1_CR2_750_1000,sphi_CR2);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_Sigma_phi_1_CR2_1000_1500,sphi_CR2);
    }
    if(Mass_tt>=1500){
      evt.set(h_Sigma_phi_1_CR2_1500_Inf,sphi_CR2);
    }
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_Sigma_phi_1_CR2_0_700,sphi_CR2);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_1_CR2_700_900,sphi_CR2);
    }
    if(Mass_tt>=1000){
      evt.set(h_Sigma_phi_1_CR2_900_Inf,sphi_CR2);
    }
  }
  if(pt_hadTop > pt_hadTop_thresh && dy_reco_CR2 <0){
    evt.set(h_Sigma_phi_2_CR2,sphi_CR2);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_Sigma_phi_2_CR2_0_500,sphi_CR2);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_Sigma_phi_2_CR2_500_750,sphi_CR2);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_2_CR2_750_1000,sphi_CR2);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_Sigma_phi_2_CR2_1000_1500,sphi_CR2);
    }
    if(Mass_tt>=1500){
      evt.set(h_Sigma_phi_2_CR2_1500_Inf,sphi_CR2);
    }
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_Sigma_phi_2_CR2_0_700,sphi_CR2);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_Sigma_phi_2_CR2_700_900,sphi_CR2);
    }
    if(Mass_tt>=1000){
      evt.set(h_Sigma_phi_2_CR2_900_Inf,sphi_CR2);
    }
  }

  if(pt_hadTop < pt_hadTop_thresh && dphi_CR2 >0){
    evt.set(h_dyreco_1_CR2,dy_reco_CR2);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_dyreco_1_CR2_0_500,dy_reco_CR2);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_dyreco_1_CR2_500_750,dy_reco_CR2);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_dyreco_1_CR2_750_1000,dy_reco_CR2);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_dyreco_1_CR2_1000_1500,dy_reco_CR2);
    }
    if(Mass_tt>=1500){
      evt.set(h_dyreco_1_CR2_1500_Inf,dy_reco_CR2);
    }
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_dyreco_1_CR2_0_700,dy_reco_CR2);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_dyreco_1_CR2_700_900,dy_reco_CR2);
    }
    if(Mass_tt>=1000){
      evt.set(h_dyreco_1_CR2_900_Inf,dy_reco_CR2);
    }
  }
  if(pt_hadTop < pt_hadTop_thresh && dphi_CR2 <0){
    evt.set(h_dyreco_2_CR2,dy_reco_CR2);
    if(Mass_tt>=0 && Mass_tt < 500){
      evt.set(h_dyreco_2_CR2_0_500,dy_reco_CR2);
    }
    if(Mass_tt>=500 && Mass_tt < 750){
      evt.set(h_dyreco_2_CR2_500_750,dy_reco_CR2);
    }
    if(Mass_tt>=750 && Mass_tt < 1000){
      evt.set(h_dyreco_2_CR2_750_1000,dy_reco_CR2);
    }
    if(Mass_tt>=1000 && Mass_tt < 1500){
      evt.set(h_dyreco_2_CR2_1000_1500,dy_reco_CR2);
    }
    if(Mass_tt>=1500){
      evt.set(h_dyreco_2_CR2_1500_Inf,dy_reco_CR2);
    }
    if(Mass_tt>=0 && Mass_tt < 700){
      evt.set(h_dyreco_2_CR2_0_700,dy_reco_CR2);
    }
    if(Mass_tt>=700 && Mass_tt < 1000){
      evt.set(h_dyreco_2_CR2_700_900,dy_reco_CR2);
    }
    if(Mass_tt>=1000){
      evt.set(h_dyreco_2_CR2_900_Inf,dy_reco_CR2);
    }
  }
  if(pt_hadTop < pt_hadTop_thresh){
    evt.set(h_Delta_phi_1_CR2,dphi_CR2);
  }
  if(pt_hadTop > pt_hadTop_thresh){
    evt.set(h_Delta_phi_2_CR2,dphi_CR2);
  }
  return true;
}

//////////////////////////////////
//  EWK corrections
//////////////////////////////////
// Generic Class for Applying SFs
void ScaleFactorsFromHistos::LoadHisto(TFile* file, std::string name, std::string hname) {
  histos[name].reset((TH1F*)file->Get(hname.c_str()));
  histos[name]->SetDirectory(0);
};

double ScaleFactorsFromHistos::Evaluator(std::string hname, double var) {
  // invalid cases
  if (var == uhh2::infinity) return 1.0;

  int firstBin = 1;
  int lastBin  = histos[hname]->GetNbinsX();
  double h_min = histos[hname]->GetBinCenter(firstBin)-histos[hname]->GetBinError(firstBin);
  double h_max = histos[hname]->GetBinCenter(lastBin)+histos[hname]->GetBinError(lastBin);
  double var_for_eval = var;
  var_for_eval = (var_for_eval > h_min) ? var_for_eval : h_min+0.001;
  var_for_eval = (var_for_eval < h_max) ? var_for_eval : h_max-0.001;
  return histos[hname]->GetBinContent(histos[hname]->FindBin(var_for_eval));
};



NLOCorrections::NLOCorrections(uhh2::Context& ctx) {

  // Corrections for 2017 and 2018 are the same. 2016 is different
  is2016 = (ctx.get("dataset_version").find("UL16") != std::string::npos);

  is_Wjets  = (ctx.get("dataset_version").find("WJets") != std::string::npos);
  is_Znn  = (ctx.get("dataset_version").find("DY_inv") != std::string::npos);
  is_DY  = (ctx.get("dataset_version").find("DY") != std::string::npos) && !is_Znn;
  is_Zjets  = is_DY || is_Znn;


  std::string folder_ = ctx.get("NLOCorrections");
  for (const std::string& proc: {"w","z"}) {
    TFile* file_ = new TFile((folder_+"merged_kfactors_"+proc+"jets.root").c_str());
    for (const std::string& corr: {"ewk","qcd","qcd_ewk"}) LoadHisto(file_, proc+"_"+corr, "kfactor_monojet_"+corr);
    file_->Close();
  }
  for (const std::string& proc: {"dy","znn"}) {
    TFile* file_ = new TFile((folder_+"kfac_"+proc+"_filter.root").c_str());
    LoadHisto(file_, proc+"_qcd_2017", "kfac_"+proc+"_filter");
    file_->Close();
  }
  TFile* file_ = new TFile((folder_+"2017_gen_v_pt_qcd_sf.root").c_str());
  LoadHisto(file_, "w_qcd_2017", "wjet_dress_inclusive");
  file_->Close();
  file_ = new TFile((folder_+"lindert_qcd_nnlo_sf.root").c_str());
  for (const std::string& proc: {"eej", "evj", "vvj"}) LoadHisto(file_, proc+"_qcd_nnlo", proc);
  file_->Close();

}


double NLOCorrections::GetPartonObjectPt(uhh2::Event& event, ParticleID objID) {
  for(const auto & gp : *event.genparticles) {if (gp.pdgId()==objID) return gp.pt(); }
  return uhh2::infinity;
};




bool NLOCorrections::process(uhh2::Event& event){
  // Sample dependant corrections
  if ((!is_Wjets && !is_Zjets) || event.isRealData) return true;
  double objpt = uhh2::infinity, theory_weight = 1.0;
  std::string process = "";

  const bool do_EWK = true;
  const bool do_QCD_EWK = false;
  const bool do_QCD_NLO  = true;
  const bool do_QCD_NNLO = false;


  if (is_Zjets) objpt = GetPartonObjectPt(event,ParticleID::Z);
  if (is_Wjets) objpt = GetPartonObjectPt(event,ParticleID::W);

  if (is_Zjets) process = "z";
  if (is_Wjets) process = "w";

  if (do_QCD_EWK) theory_weight *= Evaluator(process+"_qcd_ewk",objpt);
  else {
    if (do_EWK) theory_weight *= Evaluator(process+"_ewk",objpt);
    if (do_QCD_NLO) {
      if (!is2016) {
        if (is_DY)  process = "dy";
        if (is_Znn) process = "znn";
      }
      theory_weight *= Evaluator(process+"_qcd"+(is2016?"":"_2017"),objpt);
    }
  }

  if (do_QCD_NNLO) {
    if (is_DY)    process = "eej";
    if (is_Znn)   process = "vvj";
    if (is_Wjets) process = "evj";
    theory_weight *= Evaluator(process+"_qcd_nnlo",objpt);
  }

  event.weight *= theory_weight;
  return true;
}


// Top pT Reweight extended - from Alex F.

TopPtReweighting::TopPtReweighting(uhh2::Context& ctx,
  float a, float b,
  const std::string& syst_a,
  const std::string& syst_b,
  const std::string& ttgen_name):
  a_(a), b_(b),
  ttgen_name_(ttgen_name){

    h_weight_toppt_nominal = ctx.declare_event_output<float> ("weight_toppt_nominal");
    h_weight_toppt_a_up      = ctx.declare_event_output<float> ("weight_toppt_a_up");
    h_weight_toppt_b_up      = ctx.declare_event_output<float> ("weight_toppt_b_up");
    h_weight_toppt_a_down    = ctx.declare_event_output<float> ("weight_toppt_a_down");
    h_weight_toppt_b_down    = ctx.declare_event_output<float> ("weight_toppt_b_down");

    version_ = ctx.get("dataset_version", "");
    boost::algorithm::to_lower(version_);
    if(!ttgen_name_.empty()){
      h_ttbargen_ = ctx.get_handle<TTbarGen>(ttgen_name);
    }

    if (syst_a == "up")
    a_ *= 1.5;
    // a*1.5=0.09225
    else if (syst_a == "down")
    a_ *= 0.5;// a*0.5=0.03075

    if (syst_b == "up")
    b_ *= 1.5;//-0.00075
    else if (syst_b == "down")
    b_ *= 0.5;//-0.00025

  }

  bool TopPtReweighting::process(uhh2::Event& event){
    if (event.isRealData || (!boost::algorithm::contains(version_,"tttohadronic") && !boost::algorithm::contains(version_,"tttosemileptonic") && !boost::algorithm::starts_with(version_,"ttto2l2nu")) ) {

      event.set(h_weight_toppt_nominal, 1.0);
      event.set(h_weight_toppt_a_up, 1.0);
      event.set(h_weight_toppt_b_up, 1.0);
      event.set(h_weight_toppt_a_down, 1.0);
      event.set(h_weight_toppt_b_down, 1.0);

      return true;

    }
    const TTbarGen& ttbargen = !ttgen_name_.empty() ? event.get(h_ttbargen_) : TTbarGen(*event.genparticles,false);
    float wgt = 1.;
    float wgt_a_up = 1.;
    float wgt_a_down = 1.;
    float wgt_b_up = 1.;
    float wgt_b_down = 1.;
    if (ttbargen.DecayChannel() != TTbarGen::e_notfound) {
      float tpt1 = ttbargen.Top().v4().Pt();
      float tpt2 = ttbargen.Antitop().v4().Pt();
      wgt = sqrt(exp(a_+b_*tpt1)*exp(a_+b_*tpt2));
      wgt_a_up = sqrt(exp((1.5*a_)+b_*tpt1)*exp((1.5*a_)+b_*tpt2));
      wgt_a_down = sqrt(exp((0.5*a_)+b_*tpt1)*exp((0.5*a_)+b_*tpt2));
      wgt_b_up = sqrt(exp(a_+(1.5*b_)*tpt1)*exp(a_+(1.5*b_)*tpt2));
      wgt_b_down = sqrt(exp(a_+(0.5*b_)*tpt1)*exp(a_+(0.5*b_)*tpt2));
    }

    event.weight *= wgt;

    event.set(h_weight_toppt_nominal, wgt);
    event.set(h_weight_toppt_a_up, wgt_a_up);
    event.set(h_weight_toppt_b_up, wgt_b_up);
    event.set(h_weight_toppt_a_down, wgt_a_down);
    event.set(h_weight_toppt_b_down, wgt_b_down);


    return true;
  }


  ////

  PuppiCHS_matching::PuppiCHS_matching(uhh2::Context& ctx){

    h_CHSjets = ctx.get_handle< std::vector<Jet> >("jetsAk4CHS");
    h_CHS_matched_ = ctx.declare_event_output<vector<Jet>>("CHS_matched");
   

  }

  bool PuppiCHS_matching::process(uhh2::Event& event){

    vector<Jet> CHSjets = event.get(h_CHSjets);
    
    std::vector<Jet> matched_jets;
    std::vector<Jet> matched_jets_PUPPI;
    JetPFID CHS_matched_Tight  = JetPFID(JetPFID::WP_TIGHT_CHS);

    for(const Jet & jet : *event.jets){ // PUPPI jets
      double deltaR_min = 99;
     
      // if (CHSjets.at(0).pt()<50)continue;

      // if (CHSjets.at(1).pt()<45)continue;
      for(const Jet & CHSjet : CHSjets){ // CHS jets
        if (CHSjets.at(0).pt()<50)continue;
        if(CHS_matched_Tight(CHSjet,event)){ 
        // cout << "it passes tight ID" <<endl;
        double deltaR_CHS = deltaR(jet,CHSjet);
        if(deltaR_CHS<deltaR_min) deltaR_min = deltaR_CHS;
        }
      } // end CHS loop

      if(deltaR_min>0.2) continue;

      for(const Jet & CHSjet : CHSjets){
        if(CHS_matched_Tight(CHSjet,event)){ 
          if(deltaR(jet,CHSjet)!=deltaR_min) continue;
          else{
            matched_jets.emplace_back(CHSjet);
            matched_jets_PUPPI.emplace_back(jet);
          }
        }
      }
    // }
    } // end PUPPI loop
    std::swap(matched_jets_PUPPI, *event.jets);
    event.set(h_CHS_matched_, matched_jets);
    if(event.jets->size()==0) return false;
    return true;
    
  }

  ////

  MuonRecoSF::MuonRecoSF(uhh2::Context& ctx){
    // cout << "will do muon sf"<<endl;
    year = extract_year(ctx);
    is_mc = ctx.get("dataset_type") == "MC";
    is_Muon = ctx.get("channel") == "muon";
    
    h_muonrecSF_nominal = ctx.declare_event_output<float> ("weight_sfmu_reco");
    h_muonrecSF_up      = ctx.declare_event_output<float> ("weight_sfmu_reco_up");
    h_muonrecSF_down    = ctx.declare_event_output<float> ("weight_sfmu_reco_down");

  }

  bool MuonRecoSF::process(uhh2::Event& event){

    event.set(h_muonrecSF_nominal, 1.0);
    event.set(h_muonrecSF_up, 1.0);
    event.set(h_muonrecSF_down, 1.0);

    if(is_mc && is_Muon){
      // cout << "mc and muon"<<endl;
      // cout << "pt: "<<event.muons->at(0).pt()<<endl;
      // cout << "cos: "<<cosh(event.muons->at(0).eta())<<endl;
      float Tot_P = event.muons->at(0).pt()*cosh(event.muons->at(0).eta());
    //  cout << "Calculated pt"<<endl;
      if(year == Year::isUL16preVFP || year == Year::isUL16postVFP){
        if( abs(event.muons->at(0).eta()) <= 1.6){
          if( 50 < Tot_P && Tot_P <= 100)   { event.set(h_muonrecSF_nominal, 0.9914); event.set(h_muonrecSF_up, 0.9914+0.0008); event.set(h_muonrecSF_down, 0.9914-0.0008); event.weight *= 0.9914; }
          if( 100 < Tot_P && Tot_P <= 150)  { event.set(h_muonrecSF_nominal, 0.9936); event.set(h_muonrecSF_up, 0.9936+0.0009); event.set(h_muonrecSF_down, 0.9936-0.0009); event.weight *= 0.9936; }
          if( 150 < Tot_P && Tot_P <= 200)  { event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.993+0.001); event.set(h_muonrecSF_down, 0.993-0.001); event.weight *= 0.993; }
          if( 200 < Tot_P && Tot_P <= 300)  { event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.993+0.002); event.set(h_muonrecSF_down, 0.993-0.002); event.weight *= 0.993; }
          if( 300 < Tot_P && Tot_P <= 400)  { event.set(h_muonrecSF_nominal, 0.990); event.set(h_muonrecSF_up, 0.990+0.004); event.set(h_muonrecSF_down, 0.990-0.004); event.weight *= 0.990; }
          if( 400 < Tot_P && Tot_P <= 600)  { event.set(h_muonrecSF_nominal, 0.990); event.set(h_muonrecSF_up, 0.990+0.003); event.set(h_muonrecSF_down, 0.990-0.003); event.weight *= 0.990; }
          if( 600 < Tot_P && Tot_P <= 1500) { event.set(h_muonrecSF_nominal, 0.989); event.set(h_muonrecSF_up, 0.989+0.004); event.set(h_muonrecSF_down, 0.989-0.004); event.weight *= 0.989; }
          if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 0.8); event.set(h_muonrecSF_up, 0.8+0.3); event.set(h_muonrecSF_down, 0.8-0.3); event.weight *= 0.8; }
        }
        if(abs(event.muons->at(0).eta()) > 1.6 && abs(event.muons->at(0).eta()) < 2.4){
          if( 50 < Tot_P && Tot_P <= 100)   { event.set(h_muonrecSF_nominal, 1.0); event.set(h_muonrecSF_up, 1.0); event.set(h_muonrecSF_down, 1.0); event.weight *= 1.0; }
          if( 100 < Tot_P && Tot_P <= 150)  { event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.993+0.001); event.set(h_muonrecSF_down, 0.993-0.001); event.weight *= 0.993; }
          if( 150 < Tot_P && Tot_P <= 200)  { event.set(h_muonrecSF_nominal, 0.991); event.set(h_muonrecSF_up, 0.991+0.001); event.set(h_muonrecSF_down, 0.991-0.001); event.weight *= 0.991; }
          if( 200 < Tot_P && Tot_P <= 300)  { event.set(h_muonrecSF_nominal, 0.985); event.set(h_muonrecSF_up, 0.985+0.001); event.set(h_muonrecSF_down, 0.985-0.001); event.weight *= 0.985; }
          if( 300 < Tot_P && Tot_P <= 400)  { event.set(h_muonrecSF_nominal, 0.981); event.set(h_muonrecSF_up, 0.981+0.002); event.set(h_muonrecSF_down, 0.981-0.002); event.weight *= 0.981; }
          if( 400 < Tot_P && Tot_P <= 600)  { event.set(h_muonrecSF_nominal, 0.979); event.set(h_muonrecSF_up, 0.979+0.004); event.set(h_muonrecSF_down, 0.979-0.004); event.weight *= 0.979; }
          if( 600 < Tot_P && Tot_P <= 1500) { event.set(h_muonrecSF_nominal, 0.978); event.set(h_muonrecSF_up, 0.978+0.005); event.set(h_muonrecSF_down, 0.978-0.005); event.weight *= 0.978; }
          if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 0.9); event.set(h_muonrecSF_up, 0.9+0.2); event.set(h_muonrecSF_down, 0.9-0.2); event.weight *= 0.9; }
        }
      }
      if(year == Year::isUL17){
        if( abs(event.muons->at(0).eta()) <= 1.6){
          if( 50 < Tot_P && Tot_P <= 100)   { event.set(h_muonrecSF_nominal, 0.9938); event.set(h_muonrecSF_up, 0.9938+0.0006); event.set(h_muonrecSF_down, 0.9938-0.0006); event.weight *= 0.9938; }
          if( 100 < Tot_P && Tot_P <= 150)  { event.set(h_muonrecSF_nominal, 0.9950); event.set(h_muonrecSF_up, 0.9950+0.0007); event.set(h_muonrecSF_down, 0.9950-0.0007); event.weight *= 0.9950; }
          if( 150 < Tot_P && Tot_P <= 200)  { event.set(h_muonrecSF_nominal, 0.996); event.set(h_muonrecSF_up, 0.996+0.001); event.set(h_muonrecSF_down, 0.996-0.001); event.weight *= 0.996; }
          if( 200 < Tot_P && Tot_P <= 300)  { event.set(h_muonrecSF_nominal, 0.996); event.set(h_muonrecSF_up, 0.996+0.001); event.set(h_muonrecSF_down, 0.996-0.001); event.weight *= 0.996; }
          if( 300 < Tot_P && Tot_P <= 400)  { event.set(h_muonrecSF_nominal, 0.994); event.set(h_muonrecSF_up, 0.994+0.001); event.set(h_muonrecSF_down, 0.994-0.001); event.weight *= 0.994; }
          if( 400 < Tot_P && Tot_P <= 600)  { event.set(h_muonrecSF_nominal, 1.003); event.set(h_muonrecSF_up, 1.003+0.006); event.set(h_muonrecSF_down, 1.003-0.006); event.weight *= 1.003; }
          if( 600 < Tot_P && Tot_P <= 1500) { event.set(h_muonrecSF_nominal, 0.987); event.set(h_muonrecSF_up, 0.987+0.003); event.set(h_muonrecSF_down, 0.987-0.003); event.weight *= 0.987; }
          if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 0.9); event.set(h_muonrecSF_up, 0.9+0.1); event.set(h_muonrecSF_down, 0.9-0.1); event.weight *= 0.9; }
        }
        if(abs(event.muons->at(0).eta()) > 1.6 && abs(event.muons->at(0).eta()) < 2.4){
          if( 50 < Tot_P && Tot_P <= 100)   { event.set(h_muonrecSF_nominal, 1.0); event.set(h_muonrecSF_up, 1.0); event.set(h_muonrecSF_down, 1.0); event.weight *= 1.0; }
          if( 100 < Tot_P && Tot_P <= 150)  { event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.993+0.001); event.set(h_muonrecSF_down, 0.993-0.001); event.weight *= 0.993; }
          if( 150 < Tot_P && Tot_P <= 200)  { event.set(h_muonrecSF_nominal, 0.989); event.set(h_muonrecSF_up, 0.989+0.001); event.set(h_muonrecSF_down, 0.989-0.001); event.weight *= 0.989; }
          if( 200 < Tot_P && Tot_P <= 300)  { event.set(h_muonrecSF_nominal, 0.986); event.set(h_muonrecSF_up, 0.986+0.001); event.set(h_muonrecSF_down, 0.986-0.001); event.weight *= 0.986; }
          if( 300 < Tot_P && Tot_P <= 400)  { event.set(h_muonrecSF_nominal, 0.989); event.set(h_muonrecSF_up, 0.989+0.001); event.set(h_muonrecSF_down, 0.989-0.001); event.weight *= 0.989; }
          if( 400 < Tot_P && Tot_P <= 600)  { event.set(h_muonrecSF_nominal, 0.983); event.set(h_muonrecSF_up, 0.983+0.003); event.set(h_muonrecSF_down, 0.983-0.003); event.weight *= 0.983; }
          if( 600 < Tot_P && Tot_P <= 1500) { event.set(h_muonrecSF_nominal, 0.986); event.set(h_muonrecSF_up, 0.986+0.006); event.set(h_muonrecSF_down, 0.986-0.006); event.weight *= 0.986; }
          if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 1.01); event.set(h_muonrecSF_up, 1.01+0.01); event.set(h_muonrecSF_down, 1.01-0.01); event.weight *= 1.01; }
        }
      }
      if(year == Year::isUL18){
        // cout <<"in correct year"<<endl;
        if( abs(event.muons->at(0).eta()) <= 1.6){
          if( 50 < Tot_P && Tot_P <= 100)   { event.set(h_muonrecSF_nominal, 0.9943); event.set(h_muonrecSF_up, 0.9943+0.0007); event.set(h_muonrecSF_down, 0.9943-0.0007); event.weight *= 0.9943; }
          if( 100 < Tot_P && Tot_P <= 150)  { event.set(h_muonrecSF_nominal, 0.9948); event.set(h_muonrecSF_up, 0.9948+0.0007); event.set(h_muonrecSF_down, 0.9948-0.0007); event.weight *= 0.9948; }
          if( 150 < Tot_P && Tot_P <= 200)  { event.set(h_muonrecSF_nominal, 0.9950); event.set(h_muonrecSF_up, 0.9950+0.0009); event.set(h_muonrecSF_down, 0.9950-0.0009); event.weight *= 0.9950; }
          if( 200 < Tot_P && Tot_P <= 300)  { event.set(h_muonrecSF_nominal, 0.994); event.set(h_muonrecSF_up, 0.994+0.001); event.set(h_muonrecSF_down, 0.994-0.001); event.weight *= 0.994; }
          if( 300 < Tot_P && Tot_P <= 400)  { event.set(h_muonrecSF_nominal, 0.9914); event.set(h_muonrecSF_up, 0.9914+0.0009); event.set(h_muonrecSF_down, 0.9914-0.0009); event.weight *= 0.9914; }
          if( 400 < Tot_P && Tot_P <= 600)  { event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.993+0.002); event.set(h_muonrecSF_down, 0.993-0.002); event.weight *= 0.993; }
          if( 600 < Tot_P && Tot_P <= 1500) { event.set(h_muonrecSF_nominal, 0.991); event.set(h_muonrecSF_up, 0.991+0.004); event.set(h_muonrecSF_down, 0.991-0.004); event.weight *= 0.991; }
          if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 1.0); event.set(h_muonrecSF_up, 1.0+0.1); event.set(h_muonrecSF_down, 1.0-0.1); event.weight *= 1.0; }
        }
        if(abs(event.muons->at(0).eta()) > 1.6 && abs(event.muons->at(0).eta()) < 2.4){
          if( 50 < Tot_P && Tot_P <= 100)   { event.set(h_muonrecSF_nominal, 1.0); event.set(h_muonrecSF_up, 1.0); event.set(h_muonrecSF_down, 1.0); event.weight *= 1.0; }
          if( 100 < Tot_P && Tot_P <= 150)  { event.set(h_muonrecSF_nominal, 0.993); event.set(h_muonrecSF_up, 0.993+0.001); event.set(h_muonrecSF_down, 0.993-0.001); event.weight *= 0.993; }
          if( 150 < Tot_P && Tot_P <= 200)  { event.set(h_muonrecSF_nominal, 0.990); event.set(h_muonrecSF_up, 0.990+0.001); event.set(h_muonrecSF_down, 0.990-0.001); event.weight *= 0.990; }
          if( 200 < Tot_P && Tot_P <= 300)  { event.set(h_muonrecSF_nominal, 0.988); event.set(h_muonrecSF_up, 0.988+0.001); event.set(h_muonrecSF_down, 0.988-0.001); event.weight *= 0.988; }
          if( 300 < Tot_P && Tot_P <= 400)  { event.set(h_muonrecSF_nominal, 0.981); event.set(h_muonrecSF_up, 0.981+0.002); event.set(h_muonrecSF_down, 0.981-0.002); event.weight *= 0.981; }
          if( 400 < Tot_P && Tot_P <= 600)  { event.set(h_muonrecSF_nominal, 0.983); event.set(h_muonrecSF_up, 0.983+0.003); event.set(h_muonrecSF_down, 0.983-0.003); event.weight *= 0.983; }
          if( 600 < Tot_P && Tot_P <= 1500) { event.set(h_muonrecSF_nominal, 0.978); event.set(h_muonrecSF_up, 0.978+0.006); event.set(h_muonrecSF_down, 0.978-0.006); event.weight *= 0.978; }
          if( 1500 < Tot_P && Tot_P <= 3500){ event.set(h_muonrecSF_nominal, 0.98); event.set(h_muonrecSF_up, 0.98+0.03); event.set(h_muonrecSF_down, 0.98-0.03); event.weight *= 0.98; }
        }
      }
    }


    return true;
  }

  //boosted spin correlation variable
  
     
  //end of boosted spin correlation variable






  ////

// Implementation of StructureConstantsCalculator
StructureConstantsCalculator::StructureConstantsCalculator(uhh2::Context& ctx, int num_WCs) 
  : num_WCs_(num_WCs) {
  // Declare an output handle for storing structure constants
  h_structure_constants_ = ctx.declare_event_output<std::vector<float>>("structure_constants");
  
  // Initializes the mapping between WC names and indices
  wc_names_ = {
    "ctGRe", "ctGIm", "cQj18", "cQj38", "cQj11", "cQj31", 
    "ctu8", "ctd8", "ctj8", "cQu8", "cQd8", 
    "ctu1", "ctd1", "ctj1", "cQu1", "cQd1"
  };
  
  // Initializes the mapping between WC indices and weight indices
  single2_mapping_[0] = 219;  // ctGRe
  single2_mapping_[1] = 235;  // ctGIm
  single2_mapping_[2] = 250;  // cQj18
  single2_mapping_[3] = 264;  // cQj38
  single2_mapping_[4] = 277;  // cQj11
  single2_mapping_[5] = 289;  // cQj31
  single2_mapping_[6] = 300;  // ctu8
  single2_mapping_[7] = 310;  // ctd8
  single2_mapping_[8] = 319;  // ctj8
  single2_mapping_[9] = 327;  // cQu8
  single2_mapping_[10] = 335; // cQd8
  single2_mapping_[11] = 340; // ctu1
  single2_mapping_[12] = 345; // ctd1
  single2_mapping_[13] = 349; // ctj1
  single2_mapping_[14] = 352; // cQu1
  single2_mapping_[15] = 354; // cQd1
  
  // Initializes the known_pairs_ mapping ((i,j) -> weight index)
  // ctGRe (0) pairs: 15 pairs
  known_pairs_[std::make_pair(0, 1)] = 220; // ctGRe_ctGIm_1_1 -> EFTrwgt218
  known_pairs_[std::make_pair(0, 2)] = 221; // ctGRe_cQj18_1_1 -> EFTrwgt219
  known_pairs_[std::make_pair(0, 3)] = 222; // ctGRe_cQj38_1_1 -> EFTrwgt220    
  known_pairs_[std::make_pair(0, 4)] = 223; // ctGRe_cQj11_1_1 -> EFTrwgt221  
  known_pairs_[std::make_pair(0, 5)] = 224; // ctGRe_cQj31_1_1 -> EFTrwgt222  
  known_pairs_[std::make_pair(0, 6)] = 225; // ctGRe_ctu8_1_1 -> EFTrwgt223
  known_pairs_[std::make_pair(0, 7)] = 226; // ctGRe_ctd8_1_1 -> EFTrwgt224
  known_pairs_[std::make_pair(0, 8)] = 227; // ctGRe_ctj8_1_1 -> EFTrwgt225
  known_pairs_[std::make_pair(0, 9)] = 228; // ctGRe_cQu8_1_1 -> EFTrwgt226
  known_pairs_[std::make_pair(0, 10)] = 229; // ctGRe_cQd8_1_1 -> EFTrwgt227
  known_pairs_[std::make_pair(0, 11)] = 230; // ctGRe_ctu1_1_1 -> EFTrwgt228
  known_pairs_[std::make_pair(0, 12)] = 231; // ctGRe_ctd1_1_1 -> EFTrwgt229
  known_pairs_[std::make_pair(0, 13)] = 232; // ctGRe_ctj1_1_1 -> EFTrwgt230
  known_pairs_[std::make_pair(0, 14)] = 233; // ctGRe_cQu1_1_1 -> EFTrwgt231
  known_pairs_[std::make_pair(0, 15)] = 234; // ctGRe_cQd1_1_1 -> EFTrwgt232

  
  // ctGIm (1) pairs: 14 pairs
  known_pairs_[std::make_pair(1, 2)] = 236; // ctGIm_cQj18_1_1 -> EFTrwgt234
  known_pairs_[std::make_pair(1, 3)] = 237; // ctGIm_cQj38_1_1 -> EFTrwgt235
  known_pairs_[std::make_pair(1, 4)] = 238; // ctGIm_cQj11_1_1 -> EFTrwgt236
  known_pairs_[std::make_pair(1, 5)] = 239; // ctGIm_cQj31_1_1 -> EFTrwgt237
  known_pairs_[std::make_pair(1, 6)] = 240; // ctGIm_ctu8_1_1 -> EFTrwgt238
  known_pairs_[std::make_pair(1, 7)] = 241; // ctGIm_ctd8_1_1 -> EFTrwgt239
  known_pairs_[std::make_pair(1, 8)] = 242; // ctGIm_ctj8_1_1 -> EFTrwgt240
  known_pairs_[std::make_pair(1, 9)] = 243; // ctGIm_cQu8_1_1 -> EFTrwgt241
  known_pairs_[std::make_pair(1, 10)] = 244; // ctGIm_cQd8_1_1 -> EFTrwgt242
  known_pairs_[std::make_pair(1, 11)] = 245; // ctGIm_ctu1_1_1 -> EFTrwgt243
  known_pairs_[std::make_pair(1, 12)] = 246; // ctGIm_ctd1_1_1 -> EFTrwgt244
  known_pairs_[std::make_pair(1, 13)] = 247; // ctGIm_ctj1_1_1 -> EFTrwgt245
  known_pairs_[std::make_pair(1, 14)] = 248; // ctGIm_cQu1_1_1 -> EFTrwgt246
  known_pairs_[std::make_pair(1, 15)] = 249; // ctGIm_cQd1_1_1 -> EFTrwgt247

  // cQj18 (2) pairs: 13 pairs
  known_pairs_[std::make_pair(2, 3)] = 251; // cQj18_cQj38_1_1 -> EFTrwgt249
  known_pairs_[std::make_pair(2, 4)] = 252; // cQj18_cQj11_1_1 -> EFTrwgt250
  known_pairs_[std::make_pair(2, 5)] = 253; // cQj18_cQj31_1_1 -> EFTrwgt251
  known_pairs_[std::make_pair(2, 6)] = 254; // cQj18_ctu8_1_1 -> EFTrwgt252
  known_pairs_[std::make_pair(2, 7)] = 255; // cQj18_ctd8_1_1 -> EFTrwgt253
  known_pairs_[std::make_pair(2, 8)] = 256; // cQj18_ctj8_1_1 -> EFTrwgt254
  known_pairs_[std::make_pair(2, 9)] = 257; // cQj18_cQu8_1_1 -> EFTrwgt255
  known_pairs_[std::make_pair(2, 10)] = 258; // cQj18_cQd8_1_1 -> EFTrwgt256
  known_pairs_[std::make_pair(2, 11)] = 259; // cQj18_ctu1_1_1 -> EFTrwgt257
  known_pairs_[std::make_pair(2, 12)] = 260; // cQj18_ctd1_1_1 -> EFTrwgt258  
  known_pairs_[std::make_pair(2, 13)] = 261; // cQj18_ctj1_1_1 -> EFTrwgt259
  known_pairs_[std::make_pair(2, 14)] = 262; // cQj18_cQu1_1_1 -> EFTrwgt260
  known_pairs_[std::make_pair(2, 15)] = 263; // cQj18_cQd1_1_1 -> EFTrwgt261

  // cQj38 (3) pairs: 12 pairs
  known_pairs_[std::make_pair(3, 4)] = 265; // cQj38_cQj11_1_1 -> EFTrwgt263
  known_pairs_[std::make_pair(3, 5)] = 266; // cQj38_cQj31_1_1 -> EFTrwgt264
  known_pairs_[std::make_pair(3, 6)] = 267; // cQj38_ctu8_1_1 -> EFTrwgt265
  known_pairs_[std::make_pair(3, 7)] = 268; // cQj38_ctd8_1_1 -> EFTrwgt266
  known_pairs_[std::make_pair(3, 8)] = 269; // cQj38_ctj8_1_1 -> EFTrwgt267
  known_pairs_[std::make_pair(3, 9)] = 270; // cQj38_cQu8_1_1 -> EFTrwgt268
  known_pairs_[std::make_pair(3, 10)] = 271; // cQj38_cQd8_1_1 -> EFTrwgt269
  known_pairs_[std::make_pair(3, 11)] = 272; // cQj38_ctu1_1_1 -> EFTrwgt270
  known_pairs_[std::make_pair(3, 12)] = 273; // cQj38_ctd1_1_1 -> EFTrwgt271
  known_pairs_[std::make_pair(3, 13)] = 274; // cQj38_ctj1_1_1 -> EFTrwgt272
  known_pairs_[std::make_pair(3, 14)] = 275; // cQj38_cQu1_1_1 -> EFTrwgt273
  known_pairs_[std::make_pair(3, 15)] = 276; // cQj38_cQd1_1_1 -> EFTrwgt274

  // cQj11 (4) pairs: 11 pairs
  known_pairs_[std::make_pair(4, 5)] = 278; // cQj11_cQj31_1_1 -> EFTrwgt276
  known_pairs_[std::make_pair(4, 6)] = 279; // cQj11_ctu8_1_1 -> EFTrwgt277
  known_pairs_[std::make_pair(4, 7)] = 280; // cQj11_ctd8_1_1 -> EFTrwgt278
  known_pairs_[std::make_pair(4, 8)] = 281; // cQj11_ctj8_1_1 -> EFTrwgt279
  known_pairs_[std::make_pair(4, 9)] = 282; // cQj11_cQu8_1_1 -> EFTrwgt280     
  known_pairs_[std::make_pair(4, 10)] = 283; // cQj11_cQd8_1_1 -> EFTrwgt281
  known_pairs_[std::make_pair(4, 11)] = 284; // cQj11_ctu1_1_1 -> EFTrwgt282  
  known_pairs_[std::make_pair(4, 12)] = 285; // cQj11_ctd1_1_1 -> EFTrwgt283
  known_pairs_[std::make_pair(4, 13)] = 286; // cQj11_ctj1_1_1 -> EFTrwgt284
  known_pairs_[std::make_pair(4, 14)] = 287; // cQj11_cQu1_1_1 -> EFTrwgt285
  known_pairs_[std::make_pair(4, 15)] = 288; // cQj11_cQd1_1_1 -> EFTrwgt286

  // cQj31 (5) pairs: 10 pairs
  known_pairs_[std::make_pair(5, 6)] = 290; // cQj31_ctu8_1_1 -> EFTrwgt288
  known_pairs_[std::make_pair(5, 7)] = 291; // cQj31_ctd8_1_1 -> EFTrwgt289
  known_pairs_[std::make_pair(5, 8)] = 292; // cQj31_ctj8_1_1 -> EFTrwgt290
  known_pairs_[std::make_pair(5, 9)] = 293; // cQj31_cQu8_1_1 -> EFTrwgt291
  known_pairs_[std::make_pair(5, 10)] = 294; // cQj31_cQd8_1_1 -> EFTrwgt292
  known_pairs_[std::make_pair(5, 11)] = 295; // cQj31_ctu1_1_1 -> EFTrwgt293
  known_pairs_[std::make_pair(5, 12)] = 296; // cQj31_ctd1_1_1 -> EFTrwgt294
  known_pairs_[std::make_pair(5, 13)] = 297; // cQj31_ctj1_1_1 -> EFTrwgt295
  known_pairs_[std::make_pair(5, 14)] = 298; // cQj31_cQu1_1_1 -> EFTrwgt296
  known_pairs_[std::make_pair(5, 15)] = 299; // cQj31_cQd1_1_1 -> EFTrwgt297

  // ctu8 (6) pairs: 9 pairs
  known_pairs_[std::make_pair(6, 7)] = 301; // ctu8_ctd8_1_1 -> EFTrwgt299
  known_pairs_[std::make_pair(6, 8)] = 302; // ctu8_ctj8_1_1 -> EFTrwgt300
  known_pairs_[std::make_pair(6, 9)] = 303; // ctu8_cQu8_1_1 -> EFTrwgt301
  known_pairs_[std::make_pair(6, 10)] = 304; // ctu8_cQd8_1_1 -> EFTrwgt302  
  known_pairs_[std::make_pair(6, 11)] = 305; // ctu8_ctu1_1_1 -> EFTrwgt303
  known_pairs_[std::make_pair(6, 12)] = 306; // ctu8_ctd1_1_1 -> EFTrwgt304
  known_pairs_[std::make_pair(6, 13)] = 307; // ctu8_ctj1_1_1 -> EFTrwgt305
  known_pairs_[std::make_pair(6, 14)] = 308; // ctu8_cQu1_1_1 -> EFTrwgt306
  known_pairs_[std::make_pair(6, 15)] = 309; // ctu8_cQd1_1_1 -> EFTrwgt307

  // ctd8 (7) pairs: 8 pairs
  known_pairs_[std::make_pair(7, 8)] = 311; // ctd8_ctj8_1_1 -> EFTrwgt309
  known_pairs_[std::make_pair(7, 9)] = 312; // ctd8_cQu8_1_1 -> EFTrwgt310
  known_pairs_[std::make_pair(7, 10)] = 313; // ctd8_cQd8_1_1 -> EFTrwgt311
  known_pairs_[std::make_pair(7, 11)] = 314; // ctd8_ctu1_1_1 -> EFTrwgt312  
  known_pairs_[std::make_pair(7, 12)] = 315; // ctd8_ctd1_1_1 -> EFTrwgt313
  known_pairs_[std::make_pair(7, 13)] = 316; // ctd8_ctj1_1_1 -> EFTrwgt314
  known_pairs_[std::make_pair(7, 14)] = 317; // ctd8_cQu1_1_1 -> EFTrwgt315
  known_pairs_[std::make_pair(7, 15)] = 318; // ctd8_cQd1_1_1 -> EFTrwgt316

  // ctj8 (8) pairs: 7 pairs
  known_pairs_[std::make_pair(8, 9)] = 320; // ctj8_cQu8_1_1 -> EFTrwgt318
  known_pairs_[std::make_pair(8, 10)] = 321; // ctj8_cQd8_1_1 -> EFTrwgt319
  known_pairs_[std::make_pair(8, 11)] = 322; // ctj8_ctu1_1_1 -> EFTrwgt320
  known_pairs_[std::make_pair(8, 12)] = 323; // ctj8_ctd1_1_1 -> EFTrwgt321  
  known_pairs_[std::make_pair(8, 13)] = 324; // ctj8_ctj1_1_1 -> EFTrwgt322
  known_pairs_[std::make_pair(8, 14)] = 325; // ctj8_cQu1_1_1 -> EFTrwgt323
  known_pairs_[std::make_pair(8, 15)] = 326; // ctj8_cQd1_1_1 -> EFTrwgt324

  // cQu8 (9) pairs: 6 pairs  
  known_pairs_[std::make_pair(9, 10)] = 328; // cQu8_cQd8_1_1 -> EFTrwgt326
  known_pairs_[std::make_pair(9, 11)] = 329; // cQu8_ctu1_1_1 -> EFTrwgt327
  known_pairs_[std::make_pair(9, 12)] = 330; // cQu8_ctd1_1_1 -> EFTrwgt328
  known_pairs_[std::make_pair(9, 13)] = 331; // cQu8_ctj1_1_1 -> EFTrwgt329
  known_pairs_[std::make_pair(9, 14)] = 332; // cQu8_cQu1_1_1 -> EFTrwgt330
  known_pairs_[std::make_pair(9, 15)] = 333; // cQu8_cQd1_1_1 -> EFTrwgt331

  // cQd8 (10) pairs: 5 pairs
  known_pairs_[std::make_pair(10, 11)] = 335; // cQd8_ctu1_1_1 -> EFTrwgt333
  known_pairs_[std::make_pair(10, 12)] = 336; // cQd8_ctd1_1_1 -> EFTrwgt334
  known_pairs_[std::make_pair(10, 13)] = 337; // cQd8_ctj1_1_1 -> EFTrwgt335
  known_pairs_[std::make_pair(10, 14)] = 338; // cQd8_cQu1_1_1 -> EFTrwgt336
  known_pairs_[std::make_pair(10, 15)] = 339; // cQd8_cQd1_1_1 -> EFTrwgt337

  // ctu1 (11) pairs: 4 pairs
  known_pairs_[std::make_pair(11, 12)] = 341; // ctu1_ctd1_1_1 -> EFTrwgt339
  known_pairs_[std::make_pair(11, 13)] = 342; // ctu1_ctj1_1_1 -> EFTrwgt340
  known_pairs_[std::make_pair(11, 14)] = 343; // ctu1_cQu1_1_1 -> EFTrwgt341
  known_pairs_[std::make_pair(11, 15)] = 344; // ctu1_cQd1_1_1 -> EFTrwgt342

  // ctd1 (12) pairs: 3 pairs
  known_pairs_[std::make_pair(12, 13)] = 346; // ctd1_ctj1_1_1 -> EFTrwgt344
  known_pairs_[std::make_pair(12, 14)] = 347; // ctd1_cQu1_1_1 -> EFTrwgt345
  known_pairs_[std::make_pair(12, 15)] = 348; // ctd1_cQd1_1_1 -> EFTrwgt346

  // ctj1 (13) pairs: 2 pairs
  known_pairs_[std::make_pair(13, 14)] = 350; // ctj1_cQu1_1_1 -> EFTrwgt348
  known_pairs_[std::make_pair(13, 15)] = 351; // ctj1_cQd1_1_1 -> EFTrwgt349

  // cQu1 (14) pairs: 1 pair
  known_pairs_[std::make_pair(14, 15)] = 353; // cQu1_cQd1_1_1 -> EFTrwgt351


}

// The structure constants are calculated by solving a system of equations that relate WC configurations to weights
// To build this system, we need to know which weight in systweights() corresponds to each configuration
// get_weight_index provides this mapping, allowing us to retrieve the correct weights for each config
int StructureConstantsCalculator::get_weight_index(const std::vector<float>& config) {
  // It maps a WC configuration to the correct index in the event.genInfo->systweights() array.

  // 1) Count non-zero elements
  // Counts how many non zero values in the config
  // Stores the indices of these non zero values in non_zero_indices
  int non_zero_count = 0; // how many WCs is turned on (how many non zero)
  std::vector<int> non_zero_indices;
  
  for(size_t i = 0; i < config.size(); i++) {
    if(config[i] != 0) {
      non_zero_count++;
      non_zero_indices.push_back(i);
    }
  }
  
  // 2) SM point (all turned off)
  if(non_zero_count == 0) {
    return 201; // SM index
  }
  
  // 3) Single WC = 1 (because one WC is non-zero)
  if(non_zero_count == 1 && config[non_zero_indices[0]] == 1) {
    // Single=1 indices start at 203
    // For example, ctGRe (index 0) -> 203, ctGIm (index 1) -> 204, etc.
    return 203 + non_zero_indices[0];
  }
  
  // 4) Single WC = 2 (because one WC is non-zero)
  if(non_zero_count == 1 && config[non_zero_indices[0]] == 2) {
    // Use the mapping for Single=2 configurations
    // It uses the single2_mapping_ to look up the correct index
    // this mapping was initialized in the constructor
    return single2_mapping_[non_zero_indices[0]];
  }
  
  // 5) Pairs WC_i = 1, WC_j = 1
  if(non_zero_count == 2 && config[non_zero_indices[0]] == 1 && config[non_zero_indices[1]] == 1) {
    // Use the mapping for pair configurations
    // it uses the known_pairs_ map to look up the correct index
    // this mapping was initialized in the constructor with all 120 pairs
    return known_pairs_[std::make_pair(non_zero_indices[0], non_zero_indices[1])];
  }
  
  // 6) If we get here, the configuration is not recognized
  std::cerr << "Error: Unrecognized weight configuration" << std::endl;
  return -1;
}

bool StructureConstantsCalculator::process(uhh2::Event& event) {
  // Retrieves weights for all configurations
  // Calculates structure constants
  // Stores them in the event using event.set(h_structure_constants_, structure_constants)

  // Check if this is an EFT sample with weights
  if(!event.isRealData && event.genInfo && event.genInfo->systweights().size() > 354) {
    // Debug output for the first event
    static bool first_event = true;
    
    // Generate weight configurations
    std::vector<std::vector<float>> configs = generate_weight_configurations(num_WCs_);
    
    // Get the weights for each configuration
    // get_weight_index function is used to retrieve all the weights needed for the structure constants calculation
    std::vector<float> mg_weights;
    for(const auto& config : configs) {
      int idx = get_weight_index(config);
      if(idx >= 0 && idx < static_cast<int>(event.genInfo->systweights().size())) {
        mg_weights.push_back(event.genInfo->systweights()[idx]);
        
        // Debug output for the first event
        // if(first_event && (idx == 201 || idx == 203 || idx == 219 || idx == 220)) {
        //   std::cout << "Weight index " << idx << " (";
        //   for(size_t i = 0; i < config.size(); i++) {
        //     if(config[i] != 0) {
        //       std::cout << get_wc_name(i) << "=" << config[i] << " ";
        //     }
        //   }
        //   std::cout << "): " << event.genInfo->systweights()[idx] << std::endl;
        // }
      } else {
        std::cerr << "Error: Invalid weight index " << idx << std::endl;
        return false;
      }
    }
    
    // Calculate structure constants
    if(mg_weights.size() == configs.size()) {
      std::vector<float> structure_constants = obtain_structure_constant(num_WCs_, mg_weights);
      
      // Debug output for the first event
      // This is used to check if the structure constants are calculated correctly
      // Shows weight indices and values for key configurations (SM, single=1, single=2, pair)
      // Displays the number of structure constants calculated
      // Shows the first few structure constants
      // Only prints for the first event processed
      if(first_event) {
        std::cout << "StructureConstantsCalculator: Calculated " << structure_constants.size() << " structure constants" << std::endl;
        // if (!structure_constants.empty()) {
        //   std::cout << "First few constants: ";
        //   for (size_t i = 0; i < std::min(size_t(5), structure_constants.size()); ++i) {
        //     std::cout << structure_constants[i] << " ";
        //   }
        //   std::cout << std::endl;
        // }
                
        first_event = false;
      }
      
      // Store structure constants in the event
      event.set(h_structure_constants_, structure_constants);
    }
  }
  
  return true;
}

// Generates all possible weight configurations for 16 WCs
std::vector<std::vector<float>> StructureConstantsCalculator::generate_weight_configurations(int num_WCs) {
  // Retrieves corresponding weights from event.genInfo->systweights()
  // Builds a design matrix and solves the least squares problem to get structure constants

  // SM point (all WCs = 0) (1 configuration)
  // Single WC = 1 (16 configurations)
  // Single WC = 2 (16 configurations)
  // Pairs WC_i = 1, WC_j = 1 (120 configurations)

  // Creates a vector of size num_WCs (16)

  std::vector<std::vector<float>> configs;
  
  //1) SM point (all WCs = 0). sm_config (size = num_WCs, all elements = 0), adds this config to the configs vector
  std::vector<float> sm_config(num_WCs, 0.0);
  configs.push_back(sm_config);
  
  // 2) Single WC = 1. single_config (size = num_WCs, all = 0, except one = 1)
  // loop through each WC index from 0 to num_WCs-1
  // creates a new config with all Wcs set to 0, sets the WC at the current index to 1
  // adds this new vector to the configs vector
  for(int i = 0; i < num_WCs; i++) {
    std::vector<float> config(num_WCs, 0.0);
    config[i] = 1.0;
    configs.push_back(config);
  }
  
  // 3) Single WC = 2
  for(int i = 0; i < num_WCs; i++) {
    std::vector<float> config(num_WCs, 0.0);
    config[i] = 2.0;
    configs.push_back(config);
  }
  
  // Pairs WC_i = 1, WC_j = 1
  for(int i = 0; i < num_WCs; i++) {
    for(int j = i+1; j < num_WCs; j++) {
      std::vector<float> config(num_WCs, 0.0);
      config[i] = 1.0;
      config[j] = 1.0;
      configs.push_back(config);
    }
  }
  
  return configs;
}

// takes the weights for all configurations and solves a system of equations to get structure constants
std::vector<float> StructureConstantsCalculator::obtain_structure_constant(int num_WCs, const std::vector<float>& mg_weights) {
  
  // 1) calls the generate_weight_configurations function to get all 153 configurations
  std::vector<std::vector<float>> configs = generate_weight_configurations(num_WCs);
  
  // 2) Check if we have enough weights
  if(mg_weights.size() < configs.size()) {
    std::cerr << "Error: Not enough weights for structure constant calculation. Expected at least " 
              << configs.size() << " weights, but got " << mg_weights.size() << std::endl;
    return std::vector<float>();
  }
  
  // 3) calculate polynomial dimension
  int poly_dim = 1 + num_WCs + num_WCs + (num_WCs*(num_WCs-1))/2; // it should be 1 + 16 + 16 + 120 = 153
  // 1 constant term, 16 linear terms, 16 quadratic terms, 120 cross terms (still quadratic)

  // 4) Build design matrix
  // Eigen library is used to create the design matrix
  // Matrix is the mattrix data structure in Eigen
  // X means the size of the matrix
  // Xf means a float matrix
  // creates a matrix of size configs.size() x poly_dim (153x153). A matrix in w = As 
  Eigen::MatrixXf A_design(configs.size(), poly_dim);
  
  // 5) fill the design matrix
  for(size_t i = 0; i < configs.size(); i++) {
    // Constant term
    // sets the first column to 1 for all configs
    A_design(i, 0) = 1.0;
    
    // Linear terms
    // sets columns 1-16 to the values of each WC
    for(int j = 0; j < num_WCs; j++) {
      A_design(i, 1 + j) = configs[i][j];
    }
    
    // Quadratic terms
    // sets columns 17-32 to the squares of each WC
    int idx = 1 + num_WCs;
    // Squares
    for(int j = 0; j < num_WCs; j++) {
      A_design(i, idx) = configs[i][j] * configs[i][j];
      idx++;
    }
    
    // Cross terms
    // sets columns 33-152 to the products of each pair of WCs
    for(int j = 0; j < num_WCs; j++) {
      for(int k = j+1; k < num_WCs; k++) {
        A_design(i, idx) = configs[i][j] * configs[i][k];
        idx++;
      }
    }
  }
  
  // 6) Create weight vector directly from mg_weights.
  // w matrix in w = As
  // creates a vector of size configs.size() (153, one element for each config)
  // fills it with the weights for each config
  Eigen::VectorXf wvec(configs.size());
  for(size_t i = 0; i < configs.size(); i++) {
    wvec(i) = mg_weights[i];
  }
  
  // w = As, wvec is the vector of weights, A_design is the design matrix, s is the vector of structure constants


  // 7) Solve the least squares problem
  // solves the system of equations A_design * s = wvec to get the structure constants s
  // colPivHouseholderQr() is a method that performs a QR decomposition of the matrix A_design
  // and solves the system of equations using the QR decomposition
  Eigen::VectorXf s = A_design.colPivHouseholderQr().solve(wvec);
  
  // 8) Convert the Eigen vector to a standard vector std::vector
  std::vector<float> structure_constants(s.data(), s.data() + s.size());
  
  // 9) return the structure constants
  return structure_constants;
}


// Provides a mapping between WC names and indices (wc_names_ dictionary)
int StructureConstantsCalculator::get_wc_index(const std::string& wc_name) const {
  for (size_t i = 0; i < wc_names_.size(); i++) {
    if (wc_names_[i] == wc_name) {
      return i;
    }
  }
  std::cerr << "Error: Unknown Wilson Coefficient name: " << wc_name << std::endl;
  return -1;
}

// Provides a mapping between WC indices and names
std::string StructureConstantsCalculator::get_wc_name(int index) const {
  if (index >= 0 && index < static_cast<int>(wc_names_.size())) {
    return wc_names_[index];
  }
  std::cerr << "Error: Invalid Wilson Coefficient index: " << index << std::endl;
  return "";
}

// Calculates new weights based on structure constants and WC values. Not using this currently.
std::vector<float> StructureConstantsCalculator::calculate_new_weights(
    const std::vector<float>& structs, 
    const std::vector<float>& wc_values) {
  
 //if(wc_values.size() != static_cast<size_t>(num_WCs_) || structs.size() < 1 + num_WCs_ + num_WCs_ + (num_WCs_*(num_WCs_-1))/2) {
  // return {1.0}; // Return default weight if dimensions don't match
  //}
  
  // Constant term
  float weight = structs[0];
  
  // Linear terms
  for(int i = 0; i < num_WCs_; i++) {
    weight += structs[1 + i] * wc_values[i];
  }
  
  // Quadratic terms
  int idx = 1 + num_WCs_;
  
  // Squares
  for(int i = 0; i < num_WCs_; i++) {
    weight += structs[idx] * wc_values[i] * wc_values[i];
    idx++;
  }
  
  // Cross terms
  for(int i = 0; i < num_WCs_; i++) {
    for(int j = i+1; j < num_WCs_; j++) {
      weight += structs[idx] * wc_values[i] * wc_values[j];
      idx++;
    }
  }
  
  return {weight};
}

// When we run the analysis modules, the structure constants are stored in the event.
// Stored in the outptu root files as a branch named "structure_constants"
// Available in the event object as event.structure_constants
