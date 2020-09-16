#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/utils.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>

#include <UHH2/core/include/LorentzVector.h>

#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include <UHH2/common/include/Utils.h>

using namespace std;


BlindDataSelection::BlindDataSelection(Context& ctx, float mtt_max) : mtt_max_(mtt_max){
  h_BestZprimeCandidate_chi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  string dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
}
bool BlindDataSelection::passes(const Event & event){

  if(isMC) return true;
  bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2);
  if(!is_zprime_reconstructed_chi2) throw runtime_error("In ZprimeSemiLeptonicSelections.cxx:BlindDataSelection::passes(): The Zprime was never reconstructed via the chi2 method. This must be done before looking for the way it was reconstructed.");

  ZprimeCandidate* cand = event.get(h_BestZprimeCandidate_chi2);
  return cand->Zprime_v4().M() <= mtt_max_;
}

ZprimeTopTagSelection::ZprimeTopTagSelection(Context& ctx){
  h_BestZprimeCandidate_chi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
}
bool ZprimeTopTagSelection::passes(const Event & event){

  bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2);
  if(!is_zprime_reconstructed_chi2) throw runtime_error("In ZprimeSemiLeptonicSelections.cxx:ZprimeTopTagSelection::passes(): The Zprime was never reconstructed via the chi2 method. This must be done before looking for the way it was reconstructed.");

  ZprimeCandidate* cand = event.get(h_BestZprimeCandidate_chi2);

  return cand->is_toptag_reconstruction();

}

ZprimeBTagFatSubJetSelection::ZprimeBTagFatSubJetSelection(Context& ctx){
  // btag 
  // CSVBTag::wp btag_wp = CSVBTag::WP_TIGHT; // b-tag workingpoint
  // JetId id_btag = CSVBTag(btag_wp);

  // DeepCSVBTag::wp btag_wp = DeepCSVBTag::WP_TIGHT; // b-tag workingpoint
  // JetId id_btag = DeepCSVBTag(btag_wp);

  DeepJetBTag::wp btag_wp = DeepJetBTag::WP_TIGHT; // b-tag workingpoint
  JetId id_btag = DeepJetBTag(btag_wp);

  sel_1btag.reset(new NJetSelection(1, 1, id_btag));
  //  cout<<" init ZprimeBTagFatSubJetSelection"<<endl;
}
bool ZprimeBTagFatSubJetSelection::passes(const Event & event){
  int btag_subjet=0;
  for(auto & topjet : *event.topjets){
    auto subjets = topjet.subjets();
    for (auto & subjet : subjets) {
      if(sel_1btag->passes(event)) btag_subjet++;
    }
  }
  //  cout<<"btag_subjet = "<<btag_subjet<<endl;
  if(btag_subjet>0)  return true;
  return false;
}



Chi2CandidateMatchedSelection::Chi2CandidateMatchedSelection(Context& ctx){
  h_BestZprimeCandidate_chi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
}
bool Chi2CandidateMatchedSelection::passes(const Event & event){

  bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2);
  bool is_zprime_reconstructed_correctmatch = event.get(h_is_zprime_reconstructed_correctmatch);

  if(!(is_zprime_reconstructed_chi2 && is_zprime_reconstructed_correctmatch)) return false;
  ZprimeCandidate* cand_chi2 = event.get(h_BestZprimeCandidate_chi2);
  float dr_chi2 = cand_chi2->discriminator("correct_match");
  return dr_chi2 < 10.;

}

TTbarSemiLepFullPS::TTbarSemiLepFullPS(){}
bool TTbarSemiLepFullPS::passes(const Event & event){
  if(event.isRealData) return false;
  assert(event.genparticles);

  //Para checar que haya exactamenete 1 top_had y 1 top_lep
  bool found_had = false, found_lep = false;

  //Loop over genparticles
  for(const auto & gp : *event.genparticles){

    //Get tops
    if(fabs(gp.pdgId()) == 6){

      //Obtener b y W
      auto b = gp.daughter(event.genparticles,1);
      auto W = gp.daughter(event.genparticles,2);
      if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
        b = gp.daughter(event.genparticles,2);
        W = gp.daughter(event.genparticles,1);
     }

      //Para checar que W si sea un W
      if(abs(W->pdgId()) != 24) {
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 24)) {
            W = &genp;
            break;
          }
	}
      }

      //Para incluir eventos que no sea un b
      if(abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3 && abs(b->pdgId()) != 1) {
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 5 || abs(genp.pdgId()) == 3 || abs(genp.pdgId()) == 1)) {
            b = &genp;
            break;
          }
	}
      }
      if(!((fabs(b->pdgId()) == 5 || fabs(b->pdgId()) == 3 || fabs(b->pdgId()) == 1) && fabs(W->pdgId()) == 24)) return false;


      //Check decaymodes of W
      auto Wd1 = W->daughter(event.genparticles,1);
      auto Wd2 = W->daughter(event.genparticles,2);

      //hadronic
      if(fabs(Wd1->pdgId()) < 7 && fabs(Wd2->pdgId()) < 7){
        if(found_had) return false;
        found_had = true;
        }

      //leptonic top

      //muon channel 
      else if( abs(Wd1->pdgId()) == 13 ||  abs(Wd2->pdgId()) == 13 ){
      //electron channel 
      //else if( abs(Wd1->pdgId()) == 11 ||  abs(Wd2->pdgId()) == 11 ){ 

        if(found_lep) return false;


        // Escape cases where the W radiates an intermediate photon, that splits into llbar
        if(Wd1->pdgId() == -Wd2->pdgId()){

          int idx = 0;
          for(const auto & genp : *event.genparticles){
            if(found_lep) break;

           //muon channel 
           if(abs(genp.pdgId()) >= 13 && abs(genp.pdgId()) <= 14){
           //electron channel
           //if(abs(genp.pdgId()) >= 11 && abs(genp.pdgId()) <= 12){

              //muon channel        
              bool is_charged = (abs(genp.pdgId()) == 13);
              //electron channel  
              //bool is_charged = (abs(genp.pdgId()) == 11);

              // if the first one is charged, the second one has to have pdgId of +1 wrt. this genpart
              if(is_charged){
                // cout << "(charged) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) + 1){
                  Wd1 = &genp;
                  Wd2 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
              else{
                // cout << "(neutral) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) - 1){
                  Wd2 = &genp;
                  Wd1 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
            }
            idx++;
          }
          if(!found_lep) return false;
        }
	found_lep = true;

        //Find charged lepton
        auto lep = Wd1;
        auto nu = Wd2;

        //muon channel 
        if(fabs(Wd2->pdgId()) == 13){
        //electron channel 
        //if(fabs(Wd2->pdgId()) == 11){

          lep = Wd2;
          nu = Wd1;
        }

        //muon channel 
	if(!(abs(lep->pdgId()) == 13 && abs(nu->pdgId()) == 14)) throw runtime_error("In TTbarSemiLepMatchable: The leptonic W does not decay into a lepton and its neutrino.");
        //electron channel 
	//if(!(abs(lep->pdgId()) == 11 && abs(nu->pdgId()) == 12)) throw runtime_error("In TTbarSemiLepMatchable: The leptonic W does not decay into a lepton and its neutrino.");

      }

      else return false;

    }
  }

  if(!(found_had && found_lep)) return false;

  return true;
}

TTbarSemiLepVisibleSelection::TTbarSemiLepVisibleSelection(){}
bool TTbarSemiLepVisibleSelection::passes(const Event & event){
  if(event.isRealData) return false;
  assert(event.genparticles);
  assert(event.gentopjets);
    
  //Para checar que haya exactamenete 1 top_had y 1 top_lep
  bool found_had = false, found_lep = false, found_muon_jet = false, muon_in_subleadingjet = false, muon_in_leadingjet = false;
  float topjet_px = 0; float topjet_py = 0; float topjet_pz = 0; float topjet_E = 0;
  float subtopjet_px = 0; float subtopjet_py = 0; float subtopjet_pz = 0; float subtopjet_E = 0;
  //Loop over genparticles
    
  if(event.gentopjets->size() == 0) return false;
  for(const auto & gp : *event.genparticles){

    //Get tops
    if(fabs(gp.pdgId()) == 6){

      //Obtener b y W
      auto b = gp.daughter(event.genparticles,1);
      auto W = gp.daughter(event.genparticles,2);
      if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
        b = gp.daughter(event.genparticles,2);
        W = gp.daughter(event.genparticles,1);
     }
        
      //Para checar que W si sea un W
      if(abs(W->pdgId()) != 24) {
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 24)) {
            W = &genp;
            break;
          }
        }
      }
        
      //Para incluir eventos que no sea un b
      if(abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3 && abs(b->pdgId()) != 1) {
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 5 || abs(genp.pdgId()) == 3 || abs(genp.pdgId()) == 1)) {
            b = &genp;
            break;
          }
        }
      }
      if(!((fabs(b->pdgId()) == 5 || fabs(b->pdgId()) == 3 || fabs(b->pdgId()) == 1) && fabs(W->pdgId()) == 24)) return false;


      //Check decaymodes of W
      auto Wd1 = W->daughter(event.genparticles,1);
      auto Wd2 = W->daughter(event.genparticles,2);

      //hadronic
      if(fabs(Wd1->pdgId()) < 7 && fabs(Wd2->pdgId()) < 7){
        if(found_had) return false;
        found_had = true;
      }

      //leptonic top

      //muon channel 
      else if( abs(Wd1->pdgId()) == 13 ||  abs(Wd2->pdgId()) == 13 ){
      //electron channel 
      //else if( abs(Wd1->pdgId()) == 11 ||  abs(Wd2->pdgId()) == 11 ){

        if(found_lep) return false;
        
        // Escape cases where the W radiates an intermediate photon, that splits into llbar
        if(Wd1->pdgId() == -Wd2->pdgId()){

          int idx = 0;
          for(const auto & genp : *event.genparticles){
            if(found_lep) break;

            //muon channel 
            if(abs(genp.pdgId()) >= 13 && abs(genp.pdgId()) <= 14){
            //electron channel 
            //if(abs(genp.pdgId()) >= 11 && abs(genp.pdgId()) <= 12){  

              //muon channel          
              bool is_charged = (abs(genp.pdgId()) == 13);
              //electron channel
              //bool is_charged = (abs(genp.pdgId()) == 11);                

              // if the first one is charged, the second one has to have pdgId of +1 wrt. this genpart
              if(is_charged){
                // cout << "(charged) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) + 1){
                  Wd1 = &genp;
                  Wd2 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
              else{
                // cout << "(neutral) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) - 1){
                  Wd2 = &genp;
                  Wd1 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
            }
            idx++;
          }
          if(!found_lep) return false;
        }

        found_lep = true;

        //Find charged lepton
        auto lep = Wd1;
        auto nu = Wd2;
 
        //muon channel
        if(fabs(Wd2->pdgId()) == 13){
       	//electron channel
        //if(fabs(Wd2->pdgId()) == 11){

          lep = Wd2;
          nu = Wd1;
        }

        //muon channel 
        if(!(abs(lep->pdgId()) == 13 && abs(nu->pdgId()) == 14)) throw runtime_error("In TTbarSemiLepMatchable: The leptonic W does not decay into a lepton and its neutrino.");
        //electron channel
        //if(!(abs(lep->pdgId()) == 11 && abs(nu->pdgId()) == 12)) throw runtime_error("In TTbarSemiLepMatchable: The leptonic W does not decay into a lepton and its neutrino.");


        //cut en MET
        if (nu->pt() <= 50) return false;
        //cut en pt & eta muon
        
        for(unsigned int l = 0; l < event.genjets->size(); ++l) {
 
            //muon channel & electron channel, the name gmuonjet dosent mean "muon"
            const GenJet & gmuonjet = event.genjets->at(l);

            if (event.gentopjets->size()>=1 ){
                if(event.gentopjets->size() == 1) return false;
                const GenTopJet & gtopjet = event.gentopjets->at(0);
                const GenTopJet & gsubtopjet = event.gentopjets->at(1);
                if(deltaR(*lep,gmuonjet) <= 0.1){
                    if(found_muon_jet) break;
                    found_muon_jet = true;
                    if (gmuonjet.pt() < 55 || abs(gmuonjet.eta()) >= 2.4) return false;
                    if(deltaR(gtopjet,gmuonjet) <= 0.8){
                        muon_in_leadingjet = true;
                    }
                    if(deltaR(gsubtopjet,gmuonjet) <= 0.8){
                        muon_in_subleadingjet = true;
                    }
                    if(muon_in_leadingjet){
                            topjet_px = gtopjet.pt()*cos(gtopjet.phi()) - (gmuonjet.pt())*cos(gmuonjet.phi());
                            topjet_py = gtopjet.pt()*sin(gtopjet.phi()) - (gmuonjet.pt())*sin(gmuonjet.phi());
                            topjet_pz = gtopjet.pt()*sinh(gtopjet.eta()) - (gmuonjet.pt())*sinh(gmuonjet.eta());
                            topjet_E  = gtopjet.energy() - gmuonjet.energy();
                        
                    }else{
                            topjet_px = gtopjet.pt()*cos(gtopjet.phi());
                            topjet_py = gtopjet.pt()*sin(gtopjet.phi());
                            topjet_pz = gtopjet.pt()*sinh(gtopjet.eta());
                            topjet_E  = gtopjet.energy();
                        
                    }
                    if(muon_in_subleadingjet){
                             subtopjet_px = gsubtopjet.pt()*cos(gsubtopjet.phi()) - (gmuonjet.pt())*cos(gmuonjet.phi());
                             subtopjet_py = gsubtopjet.pt()*sin(gsubtopjet.phi()) - (gmuonjet.pt())*sin(gmuonjet.phi());
                             subtopjet_pz = gsubtopjet.pt()*sinh(gsubtopjet.eta()) - (gmuonjet.pt())*sinh(gmuonjet.eta());
                             subtopjet_E  = gsubtopjet.energy() - gmuonjet.energy();
                         
                    }else{
                             subtopjet_px = gsubtopjet.pt()*cos(gsubtopjet.phi());
                             subtopjet_py = gsubtopjet.pt()*sin(gsubtopjet.phi());
                             subtopjet_pz = gsubtopjet.pt()*sinh(gsubtopjet.eta());
                             subtopjet_E  = gsubtopjet.energy();
                         
                    }

                    if (sqrt(topjet_px*topjet_px + topjet_py*topjet_py) > sqrt(subtopjet_px*subtopjet_px + subtopjet_py*subtopjet_py)){
                         if (sqrt(topjet_px*topjet_px + topjet_py*topjet_py) < 400) return false;
                         if (abs(0.5*log((sqrt(topjet_px*topjet_px + topjet_py*topjet_py + topjet_pz*topjet_pz) + topjet_pz)/(sqrt(topjet_px*topjet_px + topjet_py*topjet_py + topjet_pz*topjet_pz) - topjet_pz))) >= 2.4) return false;
                         if (sqrt(topjet_E*topjet_E - topjet_px*topjet_px - topjet_py*topjet_py - topjet_pz*topjet_pz) < sqrt( (subtopjet_E + gmuonjet.energy())*(subtopjet_E + gmuonjet.energy()) - (subtopjet_px + gmuonjet.pt()*cos(gmuonjet.phi()))*(subtopjet_px + gmuonjet.pt()*cos(gmuonjet.phi())) - (subtopjet_py + gmuonjet.pt()*sin(gmuonjet.phi()))*(subtopjet_py + gmuonjet.pt()*sin(gmuonjet.phi())) - (subtopjet_pz + gmuonjet.pt()*sinh(gmuonjet.eta()))*(subtopjet_pz + gmuonjet.pt()*sinh(gmuonjet.eta())))) return false;
                         cout << topjet_px << " " << topjet_py << " " << topjet_pz << " " << topjet_E << endl;
                    }else{
                         if (sqrt(subtopjet_px*subtopjet_px + subtopjet_py*subtopjet_py) < 400) return false;
                         if (abs(0.5*log((sqrt(subtopjet_px*subtopjet_px + subtopjet_py*subtopjet_py + subtopjet_pz*subtopjet_pz) + subtopjet_pz)/(sqrt(subtopjet_px*subtopjet_px + subtopjet_py*subtopjet_py + subtopjet_pz*subtopjet_pz) - subtopjet_pz))) >= 2.4) return false; 
                         if (sqrt(subtopjet_E*subtopjet_E - subtopjet_px*subtopjet_px - subtopjet_py*subtopjet_py - subtopjet_pz*subtopjet_pz) < sqrt( (topjet_E + gmuonjet.energy())*(topjet_E + gmuonjet.energy()) - (topjet_px + gmuonjet.pt()*cos(gmuonjet.phi()))*(topjet_px + gmuonjet.pt()*cos(gmuonjet.phi())) - (topjet_py + gmuonjet.pt()*sin(gmuonjet.phi()))*(topjet_py + gmuonjet.pt()*sin(gmuonjet.phi())) - (topjet_pz + gmuonjet.pt()*sinh(gmuonjet.eta()))*(topjet_pz + gmuonjet.pt()*sinh(gmuonjet.eta())))) return false;
       	       	       	 cout << subtopjet_px << " " << subtopjet_py << " " << subtopjet_pz << " " << subtopjet_E << endl;
                    }
                }
            }
        }
          
          
          
        if(!found_muon_jet){
            if (event.gentopjets->size()>=1 ){
                if(event.gentopjets->size() == 1) return false;
                const GenTopJet & gtopjet = event.gentopjets->at(0);
                const GenTopJet & gsubtopjet = event.gentopjets->at(1);
                    if (lep->pt() < 55 || abs(lep->eta()) >= 2.4) return false;
                    if(deltaR(gtopjet,*lep) <= 0.8){
                        muon_in_leadingjet = true;
                    }
                    if(deltaR(gsubtopjet,*lep) <= 0.8){
                        muon_in_subleadingjet = true;
                    }
                    if(muon_in_leadingjet){
                            topjet_px = gtopjet.pt()*cos(gtopjet.phi()) - (lep->pt())*cos(lep->phi());
                            topjet_py = gtopjet.pt()*sin(gtopjet.phi()) - (lep->pt())*sin(lep->phi());
                            topjet_pz = gtopjet.pt()*sinh(gtopjet.eta()) - (lep->pt())*sinh(lep->eta());
                            topjet_E  = gtopjet.energy() - lep->energy();
                        
                    }else{
                            topjet_px = gtopjet.pt()*cos(gtopjet.phi());
                            topjet_py = gtopjet.pt()*sin(gtopjet.phi());
                            topjet_pz = gtopjet.pt()*sinh(gtopjet.eta());
                            topjet_E  = gtopjet.energy();
                        
                   }
                     if(muon_in_subleadingjet){
                             subtopjet_px = gsubtopjet.pt()*cos(gsubtopjet.phi()) - (lep->pt())*cos(lep->phi());
                             subtopjet_py = gsubtopjet.pt()*sin(gsubtopjet.phi()) - (lep->pt())*sin(lep->phi());
                             subtopjet_pz = gsubtopjet.pt()*sinh(gsubtopjet.eta()) - (lep->pt())*sinh(lep->eta());
                             subtopjet_E  = gsubtopjet.energy() - lep->energy();
                     }else{
                             subtopjet_px = gsubtopjet.pt()*cos(gsubtopjet.phi());
                             subtopjet_py = gsubtopjet.pt()*sin(gsubtopjet.phi());
                             subtopjet_pz = gsubtopjet.pt()*sinh(gsubtopjet.eta());
                             subtopjet_E  = gsubtopjet.energy();
                    }

                    if (sqrt(topjet_px*topjet_px + topjet_py*topjet_py) > sqrt(subtopjet_px*subtopjet_px + subtopjet_py*subtopjet_py)){
                         if (sqrt(topjet_px*topjet_px + topjet_py*topjet_py) < 400) return false;
                         if (abs(0.5*log((sqrt(topjet_px*topjet_px + topjet_py*topjet_py + topjet_pz*topjet_pz) + topjet_pz)/(sqrt(topjet_px*topjet_px + topjet_py*topjet_py + topjet_pz*topjet_pz) - topjet_pz))) >= 2.4) return false;
                         if (sqrt(topjet_E*topjet_E - topjet_px*topjet_px - topjet_py*topjet_py - topjet_pz*topjet_pz) < sqrt( (subtopjet_E + lep->energy())*(subtopjet_E + lep->energy()) - (subtopjet_px + lep->pt()*cos(lep->phi()))*(subtopjet_px + lep->pt()*cos(lep->phi())) - (subtopjet_py + lep->pt()*sin(lep->phi()))*(subtopjet_py + lep->pt()*sin(lep->phi())) - (subtopjet_pz + lep->pt()*sinh(lep->eta()))*(subtopjet_pz + lep->pt()*sinh(lep->eta())))) return false;
       	       	       	 cout << topjet_px << " " << topjet_py << " " << topjet_pz << " " << topjet_E << endl;
                    }else{
                         cout << "0 no era leading jet" << endl;
                         if (sqrt(subtopjet_px*subtopjet_px + subtopjet_py*subtopjet_py) < 400) return false;
                         if (abs(0.5*log((sqrt(subtopjet_px*subtopjet_px + subtopjet_py*subtopjet_py + subtopjet_pz*subtopjet_pz) + subtopjet_pz)/(sqrt(subtopjet_px*subtopjet_px + subtopjet_py*subtopjet_py + subtopjet_pz*subtopjet_pz) - subtopjet_pz))) >= 2.4) return false;
                         if (sqrt(subtopjet_E*subtopjet_E - subtopjet_px*subtopjet_px - subtopjet_py*subtopjet_py - subtopjet_pz*subtopjet_pz) < sqrt( (topjet_E + lep->energy())*(topjet_E + lep->energy()) - (topjet_px + lep->pt()*cos(lep->phi()))*(topjet_px + lep->pt()*cos(lep->phi())) - (topjet_py + lep->pt()*sin(lep->phi()))*(topjet_py + lep->pt()*sin(lep->phi())) - (topjet_pz + lep->pt()*sinh(lep->eta()))*(topjet_pz + lep->pt()*sinh(lep->eta())))) return false;
                         cout << subtopjet_px << " " << subtopjet_py << " " << subtopjet_pz << " " << subtopjet_E << endl;
                    }                
            }
        }
          
      }
        
      else return false;
        
    }
      
  }

  if(!(found_had && found_lep)) return false;

  return true;
}


TTbarSemiLepMatchableSelection::TTbarSemiLepMatchableSelection(){}
bool TTbarSemiLepMatchableSelection::passes(const Event & event){
  if(event.isRealData) return false;
  assert(event.genparticles);

  //check, if one top decays had and one decays lep
  bool found_had = false, found_lep = false;

  //Loop over genparticles
  for(const auto & gp : *event.genparticles){

    //Get tops
    if(fabs(gp.pdgId()) == 6){

      //Get b and W
      auto b = gp.daughter(event.genparticles,1);
      auto W = gp.daughter(event.genparticles,2);
      if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
        b = gp.daughter(event.genparticles,2);
        W = gp.daughter(event.genparticles,1);
      }
      if(abs(W->pdgId()) != 24) {
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 24)) {
            W = &genp;
            break;
          }
        }
      }
      if(abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3 && abs(b->pdgId()) != 1) {
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 5 || abs(genp.pdgId()) == 3 || abs(genp.pdgId()) == 1)) {
            b = &genp;
            break;
          }
        }
      }
      if(!((fabs(b->pdgId()) == 5 || fabs(b->pdgId()) == 3 || fabs(b->pdgId()) == 1) && fabs(W->pdgId()) == 24)) return false;

      //try to match the b quarks
      bool matched_b_ak4 = false;

      // Consider AK4 jets first
      for(const auto & jet : *event.jets){
        if(deltaR(*b,jet) <= 0.4) matched_b_ak4 = true;
      }

      bool matched_b_ak8 = false;
      // Now consider AK8 jets
      int idx_matched_topjet = -1;
      int idx = 0;
      for(const auto & jet : *event.topjets){
        if(deltaR(*b,jet) <= 0.8){
          matched_b_ak8 = true;
          idx_matched_topjet = idx;
        }
        idx++;
      }

      if(!matched_b_ak4 && !matched_b_ak8) return false;

      //Check decaymodes of W
      auto Wd1 = W->daughter(event.genparticles,1);
      auto Wd2 = W->daughter(event.genparticles,2);

      //hadronic
      if(fabs(Wd1->pdgId()) < 7 && fabs(Wd2->pdgId()) < 7){
        if(found_had) return false;
        found_had = true;

        //check if both daughters can be matched by jets
        bool matched_d1_ak4 = false, matched_d2_ak4 = false;
        bool matched_d1_ak8 = false, matched_d2_ak8 = false;
        // Consider AK4 jets first
        for(const auto & jet : *event.jets){
          if(deltaR(*Wd1, jet) <= 0.4) matched_d1_ak4 = true;
          if(deltaR(*Wd2, jet) <= 0.4) matched_d2_ak4 = true;
        }

        // Now consider the one AK8 jet also used for the b-jet
        if(!matched_b_ak8){
          matched_d1_ak8 = false;
          matched_d2_ak8 = false;
        }
        else{
          if(deltaR(*Wd1, event.topjets->at(idx_matched_topjet)) <= 0.8) matched_d1_ak8 = true;
          if(deltaR(*Wd2, event.topjets->at(idx_matched_topjet)) <= 0.8) matched_d2_ak8 = true;
        }

        // if(!(matched_d1 && matched_d2)) return false;
        if(!(matched_b_ak4 && matched_d1_ak4 && matched_d2_ak4) && !(matched_b_ak8 && matched_d1_ak8 && matched_d2_ak8)) return false;
      }

      //leptonic
      else if((abs(Wd1->pdgId()) == 11 || abs(Wd1->pdgId()) == 13) || (abs(Wd2->pdgId()) == 11 || abs(Wd2->pdgId()) == 13)){
        if(found_lep) return false;

        // Escape cases where the W radiates an intermediate photon, that splits into llbar
        if(Wd1->pdgId() == -Wd2->pdgId()){
          // cout << "Entered the escape-part" << endl;
          // Find 2 genparts with 11,12,13,14 that follow each other in the list and don't have the same fabs
          int idx = 0;
          for(const auto & genp : *event.genparticles){
            if(found_lep) break;
            if(abs(genp.pdgId()) >= 11 && abs(genp.pdgId()) <= 14){
              bool is_charged = (abs(genp.pdgId()) == 11 || abs(genp.pdgId()) == 13);
              // cout << "Found a genpart at index " << idx << " with id " << genp.pdgId() << ", is_charged: " << is_charged << endl;

              // if the first one is charged, the second one has to have pdgId of +1 wrt. this genpart
              if(is_charged){
                // cout << "(charged) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) + 1){
                  Wd1 = &genp;
                  Wd2 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
              else{
                // cout << "(neutral) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) - 1){
                  Wd2 = &genp;
                  Wd1 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
            }
            idx++;
          }
          if(!found_lep) return false;
        }

        found_lep = true;

        //Find charged lepton
        auto lep = Wd1;
        auto nu = Wd2;
        if(fabs(Wd2->pdgId()) == 11 || fabs(Wd2->pdgId()) == 13){
          lep = Wd2;
          nu = Wd1;
        }
        if(!(abs(lep->pdgId()) == 11 && abs(nu->pdgId()) == 12) && !(abs(lep->pdgId()) == 13 && abs(nu->pdgId()) == 14)) throw runtime_error("In TTbarSemiLepMatchable: The leptonic W does not decay into a lepton and its neutrino.");

        //check, if lepton can be matched
        bool matched_lep = false;
        if(fabs(lep->pdgId()) == 11){
          for(const auto & ele : *event.electrons){
            if(deltaR(*lep,ele) <= 0.1) matched_lep = true;
          }
        }
        else if(fabs(lep->pdgId()) == 13){
          for(const auto & mu : *event.muons){
            if(deltaR(mu,*lep) <= 0.1) matched_lep = true;
          }
        }
        else throw runtime_error("In TTbarSemiLepMatchable: Lepton from W decay is neither e nor mu.");
        if(!matched_lep) return false;
      }
      //tau-decays
      else return false;
    }
  }

  if(!(found_had && found_lep)) return false;

  return true;
}

uhh2::Chi2Cut::Chi2Cut(Context& ctx, float min, float max): min_(min), max_(max){
  h_BestZprimeCandidate = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_is_zprime_reconstructed = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
}

bool uhh2::Chi2Cut::passes(const uhh2::Event& event){

  bool is_zprime_reconstructed = event.get(h_is_zprime_reconstructed);
  // if(!is_zprime_reconstructed) throw runtime_error("In ZprimeSemiLeptonicSelections.cxx: Chi2Cut::passes: The Zprime was never reconstructed. Do this before trying to cut on its chi2.");
  bool pass = false;
  if(is_zprime_reconstructed){
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidate);
    double chi2 = BestZprimeCandidate->discriminator("chi2_total");
    if(chi2 >= min_ && (chi2 < max_ || max_ < 0)) pass = true;
  }

  return pass;
}

uhh2::STlepPlusMetCut::STlepPlusMetCut(float min, float max):
min_(min), max_(max) {}

bool uhh2::STlepPlusMetCut::passes(const uhh2::Event& event){

  float stlep = STlep(event);
  float met = event.met->pt();
  float sum = stlep + met;

  return (sum > min_) && (sum < max_ || max_ < 0.);
}
////////////////////////////////////////////////////////

uhh2::HTlepCut::HTlepCut(float min_htlep, float max_htlep):
  min_htlep_(min_htlep), max_htlep_(max_htlep) {}


bool uhh2::HTlepCut::passes(const uhh2::Event& event){
  float lep_pt = 0;
  if(event.muons) lep_pt = event.muons->at(0).pt(); //FixMe: find leading lepton first
 // cout << "MET pt: " << event.met->pt() << ", lepton pt: " << lep_pt <<endl;  
  float htlep =  event.met->pt() + lep_pt;

  return (htlep > min_htlep_) && (htlep < max_htlep_);
}


uhh2::METCut::METCut(float min_met, float max_met):
min_met_(min_met), max_met_(max_met) {}

bool uhh2::METCut::passes(const uhh2::Event& event){

  assert(event.met);

  float MET = event.met->pt();
  return (MET > min_met_) && (MET < max_met_);
}

// bool uhh2::TwoDCut1::passes(const uhh2::Event& event){
//
//   // assert((event.muons || event.electrons) && event.jets);
//   assert(event.muons || event.electrons);
//
//   const Particle* lepton = leading_lepton(event);
//
//   float drmin, ptrel;
//   bool pass = false;
//   // std::tie(drmin, ptrel) = drmin_pTrel(*lepton, *event.jets);
//   if(event.muons){
//     drmin = lepton->get_tag(Muon::twodcut_dRmin);
//     ptrel = lepton->get_tag(Muon::twodcut_pTrel);
//     if((drmin > min_deltaR_) || (ptrel > min_pTrel_)) pass = true;
//   }
//   if(event.electrons){
//     drmin = lepton->get_tag(Electron::twodcut_dRmin);
//     ptrel = lepton->get_tag(Electron::twodcut_pTrel);
//     if((drmin > min_deltaR_) || (ptrel > min_pTrel_)) pass = true;
//   }
//
//   // return (drmin > min_deltaR_) || (ptrel > min_pTrel_);
//   return pass;
// }
////////////////////////////////////////////////////////

bool uhh2::TwoDCut::passes(const uhh2::Event& event){

  assert(event.muons || event.electrons);
  if((event.muons->size()+event.electrons->size()) != 1) throw runtime_error("In TwoDCut: Event does not have exactly one muon and electron.");

  bool pass = false;
  if(event.muons->size()){
    float drmin, ptrel;
    drmin = event.muons->at(0).get_tag(Muon::twodcut_dRmin);
    ptrel = event.muons->at(0).get_tag(Muon::twodcut_pTrel);
    if((drmin > min_deltaR_) || (ptrel > min_pTrel_)) pass = true;
  }
  if(event.electrons->size()){
    float drmin, ptrel;
    drmin = event.electrons->at(0).get_tag(Electron::twodcut_dRmin);
    ptrel = event.electrons->at(0).get_tag(Electron::twodcut_pTrel);
    if((drmin > min_deltaR_) || (ptrel > min_pTrel_)) pass = true;
  }

  return pass;
}
////////////////////////////////////////////////////////

uhh2::GenFlavorSelection::GenFlavorSelection(const std::string& flav_key){

  flavor_key_ = flav_key;

  if(flavor_key_ != "l" && flavor_key_ != "c" && flavor_key_ != "b")
  throw std::runtime_error("GenFlavorSelection::GenFlavorSelection -- undefined key for parton flavor (must be 'l', 'c' or 'b'): "+flavor_key_);
}

bool uhh2::GenFlavorSelection::passes(const uhh2::Event& event){

  bool pass(false);

  assert(event.genparticles);

  int bottomN(0), charmN(0);
  for(const auto& genp : *event.genparticles){

    if(!(20 <= genp.status() && genp.status() <= 30)) continue;
    if(genp.mother1() == (unsigned short)(-1)) continue;
    if(genp.mother2() == (unsigned short)(-1)) continue;

    const int id = genp.pdgId();

    if(std::abs(id) == 5) ++bottomN;
    if(std::abs(id) == 4) ++charmN;
  }

  if     (flavor_key_ == "b") pass = (bottomN >= 1);
  else if(flavor_key_ == "c") pass = (bottomN == 0 && charmN >= 1);
  else if(flavor_key_ == "l") pass = (bottomN == 0 && charmN == 0);
  else throw std::runtime_error("GenFlavorSelection::GenFlavorSelection -- undefined key for parton flavor (must be 'l', 'c' or 'b'): "+flavor_key_);
  //  std::cout<<"bottomN = "<<bottomN<<" charmN = "<<charmN<<std::endl;

  return pass;
}

HEM_electronSelection::HEM_electronSelection(Context& ctx){
}

bool HEM_electronSelection::passes(const Event & event){
   
   bool pass=false;
//   if (event.electrons->at(0).eta() < eta_up && event.electrons->at(0).phi() < phi_up && event.electrons->at(0).phi() > phi_down) pass=true;
   for(const Electron & ele : *event.electrons){
      if ( ele.eta() < eta_up && ele.phi() < phi_up && ele.phi() > phi_down) pass=true;
   }

return pass;
}


HEM_jetSelection::HEM_jetSelection(Context& ctx){
}

bool HEM_jetSelection::passes(const Event & event){

   bool pass=false;
//   if ((event.jets->at(0).eta() < eta_up && event.jets->at(0).phi() < phi_up && event.jets->at(0).phi() > phi_down)|| (event.jets->at(1).eta() < eta_up && event.jets->at(1).phi() < phi_up && event.jets->at(1).phi() > phi_down)) pass=true;
   for(const auto & jet : *event.jets){
      if ( jet.eta() < eta_up && jet.phi() < phi_up && jet.phi() > phi_down) pass=true;
   }
  
return pass;
}


HEM_topjetSelection::HEM_topjetSelection(Context& ctx){
}

bool HEM_topjetSelection::passes(const Event & event){


   bool pass=false;
 //  if ((event.topjets->at(0).eta() < eta_up && event.topjets->at(0).phi() < phi_up && event.topjets->at(0).phi() > phi_down)|| (event.topjets->at(1).eta() < eta_up && event.topjets->at(1).phi() < phi_up && event.topjets->at(1).phi() > phi_down)) pass=true;
  for(const auto & jet : *event.topjets){
      if ( jet.eta() < eta_up && jet.phi() < phi_up && jet.phi() > phi_down) pass=true;
   }

return pass;
}

HEMSelection::HEMSelection(Context& ctx){
}
bool HEMSelection::passes(const Event & event){

   for(const Electron & ele : *event.electrons){
      if ( ele.eta() < eta_up && ele.phi() < phi_up && ele.phi() > phi_down) return false;
   }

   for(const auto & jet : *event.jets){
      if ( jet.eta() < eta_up && jet.phi() < phi_up && jet.phi() > phi_down) return false;
   }

   for(const auto & jet : *event.topjets){
      if ( jet.eta() < eta_up && jet.phi() < phi_up && jet.phi() > phi_down) return false;
   }

return true;
}
