#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include <UHH2/common/include/PrintingModules.h>

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/LuminosityHists.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/TopPtReweight.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/LeptonScaleFactors.h>
#include <UHH2/common/include/PSWeights.h>

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/TTbarLJHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicCHSMatchHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>
#include <UHH2/ZprimeSemiLeptonic/include/ElecTriggerSF.h>
#include <UHH2/ZprimeSemiLeptonic/include/AK4JetCorrections.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopPuppiJetCorrections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSystematicsModule.h>
#include <UHH2/ZprimeSemiLeptonic/include/TopTagScaleFactor.h>

#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <UHH2/HOTVR/include/HadronicTop.h>
#include <UHH2/HOTVR/include/HOTVRScaleFactor.h>
#include <UHH2/HOTVR/include/HOTVRIds.h>
#include <UHH2/core/include/LorentzVector.h>
#include <TLorentzVector.h>

using namespace std;
using namespace uhh2;

/*
██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class ZprimeAnalysisModule_GenLvl_SpinCorrelationStudy : public ModuleBASE {

public:

  explicit ZprimeAnalysisModule_GenLvl_SpinCorrelationStudy(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

protected:

  bool debug;

  // Spin Correlation variables
  std::unique_ptr<TTbarGenProducer> ttgenprod;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  
  // pt of hadronic top jet
  Event::Handle<double> h_pt_hadTop;

  // Mass of ttbar system
  Event::Handle<double> h_ttbar_mass_LabFrame;

  // Longitudinal boost of ttbar system
  Event::Handle<double> h_ttbar_boost_LabFrame;

  // Phi of decay products
  Event::Handle<double> h_phi_lep; // lepton from leptonic leg
  Event::Handle<double> h_phi_b; // b-jet from hadronic leg
  Event::Handle<double> h_phi_qlow; // less-energetic W daughter from hadronic leg

  // Sum of phi-coordinates between lepton and quark
  Event::Handle<double> h_sphi_lq;
  Event::Handle<double> h_sphi_lq_low;
  Event::Handle<double> h_sphi_lq_high;
  Event::Handle<double> h_sphi_lq_Mass1;
  Event::Handle<double> h_sphi_lq_Mass2;
  Event::Handle<double> h_sphi_lq_Mass3;
  Event::Handle<double> h_sphi_lq_Mass4;
  Event::Handle<double> h_sphi_lq_Mass5;
  // Event::Handle<double> h_sphi_lq_boost1;
  // Event::Handle<double> h_sphi_lq_boost2;
  // Event::Handle<double> h_sphi_lq_boost3;
  // Event::Handle<double> h_sphi_lq_boost4;

  // Difference of phi-coordinates between lepton and quark
  Event::Handle<double> h_dphi_lq;
  Event::Handle<double> h_dphi_lq_low;
  Event::Handle<double> h_dphi_lq_high;
  Event::Handle<double> h_dphi_lq_Mass1;
  Event::Handle<double> h_dphi_lq_Mass2;
  Event::Handle<double> h_dphi_lq_Mass3;
  Event::Handle<double> h_dphi_lq_Mass4;
  Event::Handle<double> h_dphi_lq_Mass5;
  // Event::Handle<double> h_dphi_lq_boost1;
  // Event::Handle<double> h_dphi_lq_boost2;
  // Event::Handle<double> h_dphi_lq_boost3;
  // Event::Handle<double> h_dphi_lq_boost4;

  // Sum of phi-coordinates between lepton and b-quark
  Event::Handle<double> h_sphi_lb;
  Event::Handle<double> h_sphi_lb_low;
  Event::Handle<double> h_sphi_lb_high;
  Event::Handle<double> h_sphi_lb_Mass1;
  Event::Handle<double> h_sphi_lb_Mass2;
  Event::Handle<double> h_sphi_lb_Mass3;
  Event::Handle<double> h_sphi_lb_Mass4;
  Event::Handle<double> h_sphi_lb_Mass5;
  // Event::Handle<double> h_sphi_lb_boost1;
  // Event::Handle<double> h_sphi_lb_boost2;
  // Event::Handle<double> h_sphi_lb_boost3;
  // Event::Handle<double> h_sphi_lb_boost4;

  // Difference of phi-coordinates between lepton and b-quark
  Event::Handle<double> h_dphi_lb;
  Event::Handle<double> h_dphi_lb_low;
  Event::Handle<double> h_dphi_lb_high;
  Event::Handle<double> h_dphi_lb_Mass1;
  Event::Handle<double> h_dphi_lb_Mass2;
  Event::Handle<double> h_dphi_lb_Mass3;
  Event::Handle<double> h_dphi_lb_Mass4;
  Event::Handle<double> h_dphi_lb_Mass5;
  // Event::Handle<double> h_dphi_lb_boost1;
  // Event::Handle<double> h_dphi_lb_boost2;
  // Event::Handle<double> h_dphi_lb_boost3;
  // Event::Handle<double> h_dphi_lb_boost4;
};

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

ZprimeAnalysisModule_GenLvl_SpinCorrelationStudy::ZprimeAnalysisModule_GenLvl_SpinCorrelationStudy(uhh2::Context& ctx){

  // debug = true; 
  debug = false;

  ttgenprod.reset(new TTbarGenProducer(ctx));

  // Access to gen-level particles
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen"); 

  // pt of hadronic top jet
  h_pt_hadTop=ctx.declare_event_output<double> ("pt_hadTop");
  
  // Mass of ttbar system
  h_ttbar_mass_LabFrame=ctx.declare_event_output<double> ("ttbar_mass_LabFrame");

  // Longitudinal boost of ttbar system
  h_ttbar_boost_LabFrame=ctx.declare_event_output<double> ("ttbar_boost_LabFrame");

  // Phi of decay products
  h_phi_lep=ctx.declare_event_output<double> ("phi_lep"); // lepton from leptonic leg
  h_phi_b=ctx.declare_event_output<double> ("phi_b"); // b-quark from hadronic leg
  h_phi_qlow=ctx.declare_event_output<double> ("phi_qlow"); // less-energetic (in top's rest frame) W daughter from hadronic leg

  // Sum of phi-coordinates between lepton and quark
  h_sphi_lq=ctx.declare_event_output<double> ("sphi_lq");
  h_sphi_lq_high=ctx.declare_event_output<double> ("sphi_lq_high");
  h_sphi_lq_low=ctx.declare_event_output<double> ("sphi_lq_low");
  h_sphi_lq_Mass1=ctx.declare_event_output<double> ("sphi_lq_Mass1");
  h_sphi_lq_Mass2=ctx.declare_event_output<double> ("sphi_lq_Mass2");
  h_sphi_lq_Mass3=ctx.declare_event_output<double> ("sphi_lq_Mass3");
  h_sphi_lq_Mass4=ctx.declare_event_output<double> ("sphi_lq_Mass4");
  h_sphi_lq_Mass5=ctx.declare_event_output<double> ("sphi_lq_Mass5");
  // h_sphi_lq_boost1=ctx.declare_event_output<double> ("sphi_lq_boost1");
  // h_sphi_lq_boost2=ctx.declare_event_output<double> ("sphi_lq_boost2");
  // h_sphi_lq_boost3=ctx.declare_event_output<double> ("sphi_lq_boost3");
  // h_sphi_lq_boost4=ctx.declare_event_output<double> ("sphi_lq_boost4");

  // Difference of phi-coordinates between lepton and quark
  h_dphi_lq=ctx.declare_event_output<double> ("dphi_lq");
  h_dphi_lq_high=ctx.declare_event_output<double> ("dphi_lq_high");
  h_dphi_lq_low=ctx.declare_event_output<double> ("dphi_lq_low");
  h_dphi_lq_Mass1=ctx.declare_event_output<double> ("dphi_lq_Mass1");
  h_dphi_lq_Mass2=ctx.declare_event_output<double> ("dphi_lq_Mass2");
  h_dphi_lq_Mass3=ctx.declare_event_output<double> ("dphi_lq_Mass3");
  h_dphi_lq_Mass4=ctx.declare_event_output<double> ("dphi_lq_Mass4");
  h_dphi_lq_Mass5=ctx.declare_event_output<double> ("dphi_lq_Mass5");
  // h_dphi_lq_boost1=ctx.declare_event_output<double> ("dphi_lq_boost1");
  // h_dphi_lq_boost2=ctx.declare_event_output<double> ("dphi_lq_boost2");
  // h_dphi_lq_boost3=ctx.declare_event_output<double> ("dphi_lq_boost3");
  // h_dphi_lq_boost4=ctx.declare_event_output<double> ("dphi_lq_boost4");

  // Sum of phi-coordinates between lepton and b-quark
  h_sphi_lb=ctx.declare_event_output<double> ("sphi_lb");
  h_sphi_lb_high=ctx.declare_event_output<double> ("sphi_lb_high");
  h_sphi_lb_low=ctx.declare_event_output<double> ("sphi_lb_low");
  h_sphi_lb_Mass1=ctx.declare_event_output<double> ("sphi_lb_Mass1");
  h_sphi_lb_Mass2=ctx.declare_event_output<double> ("sphi_lb_Mass2");
  h_sphi_lb_Mass3=ctx.declare_event_output<double> ("sphi_lb_Mass3");
  h_sphi_lb_Mass4=ctx.declare_event_output<double> ("sphi_lb_Mass4");
  h_sphi_lb_Mass5=ctx.declare_event_output<double> ("sphi_lb_Mass5");
  // h_sphi_lb_boost1=ctx.declare_event_output<double> ("sphi_lb_boost1");
  // h_sphi_lb_boost2=ctx.declare_event_output<double> ("sphi_lb_boost2");
  // h_sphi_lb_boost3=ctx.declare_event_output<double> ("sphi_lb_boost3");
  // h_sphi_lb_boost4=ctx.declare_event_output<double> ("sphi_lb_boost4");

  // Difference of phi-coordinates between lepton and b-quark
  h_dphi_lb=ctx.declare_event_output<double> ("dphi_lb");
  h_dphi_lb_high=ctx.declare_event_output<double> ("dphi_lb_high");
  h_dphi_lb_low=ctx.declare_event_output<double> ("dphi_lb_low");
  h_dphi_lb_Mass1=ctx.declare_event_output<double> ("dphi_lb_Mass1");
  h_dphi_lb_Mass2=ctx.declare_event_output<double> ("dphi_lb_Mass2");
  h_dphi_lb_Mass3=ctx.declare_event_output<double> ("dphi_lb_Mass3");
  h_dphi_lb_Mass4=ctx.declare_event_output<double> ("dphi_lb_Mass4");
  h_dphi_lb_Mass5=ctx.declare_event_output<double> ("dphi_lb_Mass5");
  // h_dphi_lb_boost1=ctx.declare_event_output<double> ("dphi_lb_boost1");
  // h_dphi_lb_boost2=ctx.declare_event_output<double> ("dphi_lb_boost2");
  // h_dphi_lb_boost3=ctx.declare_event_output<double> ("dphi_lb_boost3");
  // h_dphi_lb_boost4=ctx.declare_event_output<double> ("dphi_lb_boost4");
}

/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool ZprimeAnalysisModule_GenLvl_SpinCorrelationStudy::process(uhh2::Event& event){
  if(debug) cout << "\n" << endl;
  if(debug) cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug) cout << "run-event: " << event.run << "- " << event.event << endl;

  // pt of hadronic top
  event.set(h_pt_hadTop, -10);

  // Mass of ttbar system
  event.set(h_ttbar_mass_LabFrame, -10);

  // Plotting longitudinal boost of ttbar system
  event.set(h_ttbar_boost_LabFrame, -10);

  // Phi of decay products
  event.set(h_phi_lep, -10); // lepton from leptonic leg
  event.set(h_phi_b, -10); // b-quark from hadronic leg
  event.set(h_phi_qlow, -10); // less-energetic (in top's rest-frame) W daughter from hadronic leg

  // Sum of phi-coordinates between lepton and quark
  event.set(h_sphi_lq, -10);     
  event.set(h_sphi_lq_low, -10); 
  event.set(h_sphi_lq_high, -10);
  event.set(h_sphi_lq_Mass1, -10);
  event.set(h_sphi_lq_Mass2, -10);
  event.set(h_sphi_lq_Mass3, -10);
  event.set(h_sphi_lq_Mass4, -10);
  event.set(h_sphi_lq_Mass5, -10);
  
  // Difference of phi-coordinates between lepton and quark
  event.set(h_dphi_lq, -10);     
  event.set(h_dphi_lq_low, -10); 
  event.set(h_dphi_lq_high, -10);
  event.set(h_dphi_lq_Mass1, -10);
  event.set(h_dphi_lq_Mass2, -10);
  event.set(h_dphi_lq_Mass3, -10);
  event.set(h_dphi_lq_Mass4, -10);
  event.set(h_dphi_lq_Mass5, -10);

  // Sum of phi-coordinates between lepton and b-quark
  event.set(h_sphi_lb, -10);     
  event.set(h_sphi_lb_low, -10); 
  event.set(h_sphi_lb_high, -10);
  event.set(h_sphi_lb_Mass1, -10);
  event.set(h_sphi_lb_Mass2, -10);
  event.set(h_sphi_lb_Mass3, -10);
  event.set(h_sphi_lb_Mass4, -10);
  event.set(h_sphi_lb_Mass5, -10);

  
  // Difference of phi-coordinates between lepton and b-quark
  event.set(h_dphi_lb, -10);     
  event.set(h_dphi_lb_low, -10); 
  event.set(h_dphi_lb_high, -10);
  event.set(h_dphi_lb_Mass1, -10);
  event.set(h_dphi_lb_Mass2, -10);
  event.set(h_dphi_lb_Mass3, -10);
  event.set(h_dphi_lb_Mass4, -10);
  event.set(h_dphi_lb_Mass5, -10);

  // fill ttbargen information
  ttgenprod->process(event);

  // Variable to access gen-particles
  TTbarGen ttbargen = event.get(h_ttbargen);

  // Make sure ttbar decays semileptonically
  if(ttbargen.IsSemiLeptonicDecay()){

    // Define cut-variable as pt of hadTop for low/high regions
    auto pt_hadTop_thresh = 150.;         

    //-------------------------------------------------- Start in LAB-frame --------------------------------------------------//
    if(debug) cout<<"Lab Frame"<<endl;
    // pt of hadronic Top 
    LorentzVector Gen_HadTop = ttbargen.TopHad().v4();
    if(Gen_HadTop.pt() != -10) event.set(h_pt_hadTop, Gen_HadTop.pt());

    //// Defining 4vectors of ttbar system

    // Top vectors
    TLorentzVector PosTop;
    PosTop.SetPtEtaPhiE(ttbargen.Top().v4().pt(), ttbargen.Top().v4().eta(), ttbargen.Top().v4().phi(), ttbargen.Top().v4().energy());
    TLorentzVector NegTop;
    NegTop.SetPtEtaPhiE(ttbargen.Antitop().v4().pt(), ttbargen.Antitop().v4().eta(), ttbargen.Antitop().v4().phi(), ttbargen.Antitop().v4().energy());

    // Print 4vectors of top quarks
    if(debug) cout<<" top ("<<PosTop.Pt()<<", "<<PosTop.Eta()<<", "<<PosTop.Phi()<<", "<<PosTop.E()<<")"<<endl;
    if(debug) cout<<" antitop ("<<NegTop.Pt()<<", "<<NegTop.Eta()<<", "<<NegTop.Phi()<<", "<<NegTop.E()<<")"<<endl;

    // lepton
    LorentzVector Gen_Lep = ttbargen.ChargedLepton().v4(); 
    TLorentzVector lepTop_lep;
    lepTop_lep.SetPtEtaPhiE(Gen_Lep.pt(), Gen_Lep.eta(), Gen_Lep.phi(), Gen_Lep.energy());

    // Print 4vector of lepton
    if(debug) cout<<" lepton ("<<lepTop_lep.Pt()<<", "<<lepTop_lep.Eta()<<", "<<lepTop_lep.Phi()<<", "<<lepTop_lep.E()<<")"<<endl;


    // Assign hadronic b-quark based on sign of lepton
    TLorentzVector hadTop_b;
    LorentzVector Gen_b = ttbargen.BHad().v4();
    hadTop_b.SetPtEtaPhiE(Gen_b.pt(),Gen_b.eta(),Gen_b.phi(),Gen_b.energy());
    // if(ttbargen.ChargedLepton().charge() > 0){ // Positive lepton means hadronic b-quark is from antitop
    //   LorentzVector Gen_b = ttbargen.bAntitop().v4(); 
    //   hadTop_b.SetPtEtaPhiE(Gen_b.pt(),Gen_b.eta(),Gen_b.phi(),Gen_b.energy());
    //   // Print 4vector and charge of b-quark
    //   if(debug) cout<<" b-quark ("<<hadTop_b.Pt()<<", "<<hadTop_b.Eta()<<", "<<hadTop_b.Phi()<<", "<<hadTop_b.E()<<")"<<endl;
    // }
    // else{ // Negative lepton means hadronic b-quark is from top
    //   LorentzVector Gen_b = ttbargen.bTop().v4(); 
    //   hadTop_b.SetPtEtaPhiE(Gen_b.pt(),Gen_b.eta(),Gen_b.phi(),Gen_b.energy());
    //   // Print 4vector of b-quark
    //   if(debug) cout<<" b-quark ("<<hadTop_b.Pt()<<", "<<hadTop_b.Eta()<<", "<<hadTop_b.Phi()<<", "<<hadTop_b.E()<<")"<<endl;
    // }


    // 4vector to represent ttbar system
    TLorentzVector ttbar(PosTop + NegTop);

    // Plotting mass of ttbar system
    event.set(h_ttbar_mass_LabFrame, ttbar.M());

    //---------------------------------------------------------- Boosting into CoM-frame ----------------------------------------------------------//
    // Boost into ttbar Center of Momentum configuration 
    TLorentzVector lepTop_lep_CoM = lepTop_lep;
    lepTop_lep_CoM.Boost(-1. * ttbar.BoostVector());
    TLorentzVector hadTop_b_CoM = hadTop_b;
    hadTop_b_CoM.Boost(-1. * ttbar.BoostVector());
    TLorentzVector PosTop_CoM = PosTop;
    PosTop_CoM.Boost(-1. * ttbar.BoostVector());
    TLorentzVector NegTop_CoM = NegTop;
    NegTop_CoM.Boost(-1. * ttbar.BoostVector());

    // Beam unit vector in COM frame
    TVector3 beam_axis(0,0,1);

    //// Bernreuther related variables
    // Calculating top scattering angle for PosTop only
    double cos_PosTop_beam = PosTop_CoM.Vect().Unit().Dot(beam_axis);
    double sin_PosTop_beam = sqrt(1 - cos_PosTop_beam*cos_PosTop_beam);

    //// Sign of scattering angle to account for Bose symmetry
    // The sign of cos_PosTop_beam
    double sign_cos_PosTop_beam = (cos_PosTop_beam > 0.) ? 1. : -1.;
    // // The sign based on PosTop and NegTop's rapidity
    // double sign_rapidity = (PosTop.Rapidity() >= NegTop.Rapidity()) ? 1. : -1.;

    // Bernreuther basis vectors
    TVector3 kbase = PosTop_CoM.Vect().Unit();
    TVector3 rbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*(beam_axis - cos_PosTop_beam * kbase) ).Unit();
    TVector3 nbase = ( (sign_cos_PosTop_beam/sin_PosTop_beam)*beam_axis.Cross(kbase) ).Unit();

    if(debug) cout<<"Center of Momentum Frame"<<endl;
    // Print 4vectors of entire system
    if(debug) cout<<" top ("<<PosTop_CoM.Pt()<<", "<<PosTop_CoM.Eta()<<", "<<PosTop_CoM.Phi()<<", "<<PosTop_CoM.E()<<")"<<endl;
    if(debug) cout<<" antitop ("<<NegTop_CoM.Pt()<<", "<<NegTop_CoM.Eta()<<", "<<NegTop_CoM.Phi()<<", "<<NegTop_CoM.E()<<")"<<endl;
    if(debug) cout<<" lepton ("<<lepTop_lep_CoM.Pt()<<", "<<lepTop_lep_CoM.Eta()<<", "<<lepTop_lep_CoM.Phi()<<", "<<lepTop_lep_CoM.E()<<")"<<endl;
    if(debug) cout<<" b-quark ("<<hadTop_b_CoM.Pt()<<", "<<hadTop_b_CoM.Eta()<<", "<<hadTop_b_CoM.Phi()<<", "<<hadTop_b_CoM.E()<<")"<<endl;
    // if(debug) cout<<" sign of scattering angle is "<<sign_cos_PosTop_beam<<endl;
    // if(debug) cout<<" kbase ("<<kbase.X()<<", "<<kbase.Y()<<", "<<kbase.Z()<<")"<<endl;
    // if(debug) cout<<" rbase ("<<rbase.X()<<", "<<rbase.Y()<<", "<<rbase.Z()<<")"<<endl;
    // if(debug) cout<<" nbase ("<<nbase.X()<<", "<<nbase.Y()<<", "<<nbase.Z()<<")"<<endl;
    //---------------------------------------------------------- Boosted into CoM-frame ----------------------------------------------------------//

    //---------------------------------------------------- Start Rotation into Helicity Frame ----------------------------------------------------//
    // Rotate tops and decay products about beam-line
    TLorentzVector lepTop_lep_H = lepTop_lep_CoM;
    lepTop_lep_H.RotateZ(-1.*PosTop_CoM.Phi());
    TLorentzVector hadTop_b_H = hadTop_b_CoM;
    hadTop_b_H.RotateZ(-1.*PosTop_CoM.Phi());
    TLorentzVector PosTop_H = PosTop_CoM;
    PosTop_H.RotateZ(-1.*PosTop_CoM.Phi());
    TLorentzVector NegTop_H = NegTop_CoM;
    NegTop_H.RotateZ(-1.*PosTop_CoM.Phi());

    TVector3 kbase_H = kbase;
    kbase_H.RotateZ(-1.*PosTop_CoM.Phi());
    // if(debug) cout<<"kbase ("<<kbase_H.X()<<", "<<kbase_H.Y()<<", "<<kbase_H.Z()<<") after rotation abt Z by "<< -1.*PosTop_CoM.Phi()<< " radians"  <<endl;
    TVector3 rbase_H = rbase;
    rbase_H.RotateZ(-1.*PosTop_CoM.Phi());
    TVector3 nbase_H = nbase;
    nbase_H.RotateZ(-1.*PosTop_CoM.Phi());

    // Rotate tops and decay products about y-axis
    TLorentzVector lepTop_lep_Hel = lepTop_lep_H;
    lepTop_lep_Hel.RotateY(-1.*PosTop_CoM.Theta());
    TLorentzVector hadTop_b_Hel = hadTop_b_H;
    hadTop_b_Hel.RotateY(-1.*PosTop_CoM.Theta());
    TLorentzVector PosTop_Hel = PosTop_H;
    PosTop_Hel.RotateY(-1.*PosTop_CoM.Theta());
    TLorentzVector NegTop_Hel = NegTop_H;
    NegTop_Hel.RotateY(-1.*PosTop_CoM.Theta());

    TVector3 kbase_Hel = kbase_H;
    kbase_Hel.RotateY(-1.*PosTop_CoM.Theta());
    // if(debug) cout<<"kbase ("<<kbase_Hel.X()<<", "<<kbase_Hel.Y()<<", "<<kbase_Hel.Z()<<") after rotation abt Y by "<< -1.*PosTop_CoM.Theta() << " radians" <<endl;
    TVector3 rbase_Hel = rbase_H;
    rbase_Hel.RotateY(-1.*PosTop_CoM.Theta());
    TVector3 nbase_Hel = nbase_H;
    nbase_Hel.RotateY(-1.*PosTop_CoM.Theta());

    if(debug) cout<<"Helicity Frame"<<endl;
    // Print 4vectors of entire system
    if(debug) cout<<" top ("<<PosTop_Hel.Pt()<<", "<<PosTop_Hel.Eta()<<", "<<PosTop_Hel.Phi()<<", "<<PosTop_Hel.E()<<")"<<endl;
    if(debug) cout<<" antitop ("<<NegTop_Hel.Pt()<<", "<<NegTop_Hel.Eta()<<", "<<NegTop_Hel.Phi()<<", "<<NegTop_Hel.E()<<")"<<endl;
    if(debug) cout<<" lepton ("<<lepTop_lep_Hel.Pt()<<", "<<lepTop_lep_Hel.Eta()<<", "<<lepTop_lep_Hel.Phi()<<", "<<lepTop_lep_Hel.E()<<")"<<endl;
    if(debug) cout<<" b-quark ("<<hadTop_b_Hel.Pt()<<", "<<hadTop_b_Hel.Eta()<<", "<<hadTop_b_Hel.Phi()<<", "<<hadTop_b_Hel.E()<<")"<<endl;
    // if(debug) cout<<" kbase ("<<kbase_Hel.X()<<", "<<kbase_Hel.Y()<<", "<<kbase_Hel.Z()<<")"<<endl;
    // if(debug) cout<<" rbase ("<<rbase_Hel.X()<<", "<<rbase_Hel.Y()<<", "<<rbase_Hel.Z()<<")"<<endl;
    // if(debug) cout<<" nbase ("<<nbase_Hel.X()<<", "<<nbase_Hel.Y()<<", "<<nbase_Hel.Z()<<")"<<endl;
    //---------------------------------------------------- End Rotation into Helicity Frame ----------------------------------------------------//

    //----------------- Apply rotation to align Bernreuther basis vectors based on sign of scattering angle -----------------//
    TLorentzVector lepTop_lep_BoseSymm = lepTop_lep_Hel;
    TLorentzVector hadTop_b_BoseSymm = hadTop_b_Hel;
    TLorentzVector PosTop_BoseSymm = PosTop_Hel;
    TLorentzVector NegTop_BoseSymm = NegTop_Hel;

    TVector3 kbase_BoseSymm = kbase_Hel;
    TVector3 rbase_BoseSymm = rbase_Hel;
    TVector3 nbase_BoseSymm = nbase_Hel;

    if(sign_cos_PosTop_beam > 0.){
      lepTop_lep_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
      hadTop_b_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
      PosTop_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
      NegTop_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);

      kbase_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
      rbase_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
      nbase_BoseSymm.RotateZ(-1.*TMath::Pi()/2.);
    }
    else{
      lepTop_lep_BoseSymm.RotateZ(TMath::Pi()/2.);
      hadTop_b_BoseSymm.RotateZ(TMath::Pi()/2.);
      PosTop_BoseSymm.RotateZ(TMath::Pi()/2.);
      NegTop_BoseSymm.RotateZ(TMath::Pi()/2.);

      kbase_BoseSymm.RotateZ(TMath::Pi()/2.);
      rbase_BoseSymm.RotateZ(TMath::Pi()/2.);
      nbase_BoseSymm.RotateZ(TMath::Pi()/2.);
    }

    // if(debug) cout<<"Coordinates now wrt Bernreuther basis:"<<endl;
    // // Print 3vectors of all relevant particles
    // if(debug) cout<<" top ("<<PosTop_BoseSymm.Pt()<<", "<<PosTop_BoseSymm.Eta()<<", "<<PosTop_BoseSymm.Phi()<<")"<<endl;
    // if(debug) cout<<" antitop ("<<NegTop_BoseSymm.Pt()<<", "<<NegTop_BoseSymm.Eta()<<", "<<NegTop_BoseSymm.Phi()<<")"<<endl;
    // if(debug) cout<<" lepton ("<<lepTop_lep_BoseSymm.Pt()<<", "<<lepTop_lep_BoseSymm.Eta()<<", "<<lepTop_lep_BoseSymm.Phi()<<")"<<endl;
    // if(debug) cout<<" b-quark ("<<hadTop_b_BoseSymm.Pt()<<", "<<hadTop_b_BoseSymm.Eta()<<", "<<hadTop_b_BoseSymm.Phi()<<")"<<endl;
    // if(debug) cout<<" kbase ("<<kbase_BoseSymm.X()<<", "<<kbase_BoseSymm.Y()<<", "<<kbase_BoseSymm.Z()<<")"<<endl;
    // if(debug) cout<<" rbase ("<<rbase_BoseSymm.X()<<", "<<rbase_BoseSymm.Y()<<", "<<rbase_BoseSymm.Z()<<")"<<endl;
    // if(debug) cout<<" nbase ("<<nbase_BoseSymm.X()<<", "<<nbase_BoseSymm.Y()<<", "<<nbase_BoseSymm.Z()<<")"<<endl;
    //----------------- Applied rotation to align Bernreuther basis vectors based on sign of scattering angle -----------------//

    //--------------------------- Boosting into ttbar rest-frame ---------------------------//
    // Boost the top's to rest individually, bringing their children with them
    TLorentzVector lepTop_lep_Rest = lepTop_lep_BoseSymm;
    TLorentzVector hadTop_b_Rest = hadTop_b_BoseSymm;
    TLorentzVector PosTop_Rest = PosTop_BoseSymm;
    TLorentzVector NegTop_Rest = NegTop_BoseSymm;

    // Decay products get boosted in opposite directions depending on their mother top
    // POSITIVE LEPTON CONFIGURATION
    if(ttbargen.ChargedLepton().charge() > 0){ // Positively charged lepton means
      lepTop_lep_Rest.Boost(-1.*PosTop_BoseSymm.BoostVector()); // lepton has Positive Top mother
      hadTop_b_Rest.Boost(-1.*NegTop_BoseSymm.BoostVector());   // b-jet has Negative Top mother

      if(debug) cout<<"Top quark Rest Frame"<<endl;
      // Print 4vectors of all relevant particles
      if(debug) cout<<" lepton ("<<lepTop_lep_Rest.Pt()<<", "<<lepTop_lep_Rest.Eta()<<", "<<lepTop_lep_Rest.Phi()<<", "<<lepTop_lep_Rest.E()<<")"<<endl;
      if(debug) cout<<" b-quark ("<<hadTop_b_Rest.Pt()<<", "<<hadTop_b_Rest.Eta()<<", "<<hadTop_b_Rest.Phi()<<", "<<hadTop_b_Rest.E()<<")"<<endl;
    }
    // NEGATIVE LEPTON CONFIGURATION
    if(ttbargen.ChargedLepton().charge() < 0){ // Negatively charged lepton means
      lepTop_lep_Rest.Boost(-1.*NegTop_BoseSymm.BoostVector()); // lepton has Negative Top mother
      hadTop_b_Rest.Boost(-1.*PosTop_BoseSymm.BoostVector());   // b-jet has Positive Top mother

      if(debug) cout<<"Top quark Rest Frame"<<endl;
      // Print 4vectors of all relevant particles
      if(debug) cout<<" lepton ("<<lepTop_lep_Rest.Pt()<<", "<<lepTop_lep_Rest.Eta()<<", "<<lepTop_lep_Rest.Phi()<<", "<<lepTop_lep_Rest.E()<<")"<<endl;
      if(debug) cout<<" b-quark ("<<hadTop_b_Rest.Pt()<<", "<<hadTop_b_Rest.Eta()<<", "<<hadTop_b_Rest.Phi()<<", "<<hadTop_b_Rest.E()<<")"<<endl;
    }

    // Print opening angle between lepton and b-quark
    if(debug) cout<<" Opening angle between lepton and b-quark: "<<lepTop_lep_Rest.Angle(hadTop_b_Rest.Vect())<<endl;

    // // Print the charge of the decay particles
    // if(debug) cout<<" Charge of lepton: "<<ttbargen.ChargedLepton().charge()<<endl;
    // if(ttbargen.ChargedLepton().charge() > 0){ 
    //   if(debug) cout<<" Charge of b-quark: "<<ttbargen.bAntitop().charge()<<endl;
    // }
    // else{
    //   if(debug) cout<<" Charge of b-quark: "<<ttbargen.bTop().charge()<<endl;
    // }

    // Print the PID of the decay particles
    if(debug) cout<<" PID of lepton: "<<ttbargen.ChargedLepton().pdgId()<<endl;
    if(debug) cout<<" PID of b-quark: "<<ttbargen.BHad().pdgId()<<endl;

    //--------------------------- Boosted into ttbar rest-frame ---------------------------//
    if(debug) cout<<"---"<<endl;

    // Plot phi-coordinates from the boosted 4vectors
    event.set(h_phi_lep, lepTop_lep_Rest.Phi());
    event.set(h_phi_b, hadTop_b_Rest.Phi());

    // sphi and dphi are both defined as: PosTop_phi +- NegTop_phi
    auto pie = TMath::Pi();

    if(debug) cout<<" Angular Variables:"<<endl;

    // Positive lepton means PosTop has leptonic decay
    if(ttbargen.ChargedLepton().charge() > 0){

      // Define sphi_lb from phi coordinates of lepton and b-quark
      auto sphi_lb = lepTop_lep_Rest.Phi() + hadTop_b_Rest.Phi();
      if(sphi_lb > pie) sphi_lb = sphi_lb - 2.*pie;
      if(sphi_lb < -1.*pie) sphi_lb = sphi_lb + 2.*pie;
      event.set(h_sphi_lb, sphi_lb);
      if(debug) cout<<" SigmaPhi of lepton and b-quark is: "<<sphi_lb<<endl;
      
      // Define dphi_lb from phi coordinates of lepton and b-quark
      auto dphi_lb = lepTop_lep_Rest.Phi() - hadTop_b_Rest.Phi();
      if(dphi_lb > pie) dphi_lb = dphi_lb - 2.*pie;
      if(dphi_lb < -1.*pie) dphi_lb = dphi_lb + 2.*pie;
      event.set(h_dphi_lb, dphi_lb);
      if(debug) cout<<" DeltaPhi of lepton and b-quark is: "<<dphi_lb<<endl;

      // Plot dphi and sphi for high-pt ranges
      if(Gen_HadTop.pt() > pt_hadTop_thresh){
        event.set(h_sphi_lb_high, sphi_lb);
        event.set(h_dphi_lb_high, dphi_lb);
      }
      // Plot dphi and sphi for low-pt ranges
      if(Gen_HadTop.pt() < pt_hadTop_thresh){
        event.set(h_sphi_lb_low, sphi_lb);
        event.set(h_dphi_lb_low, dphi_lb);
      }

    }

    // Negative lepton means PosTop has hadronic decay
    if(ttbargen.ChargedLepton().charge() < 0){

      // Define sphi_lb from phi coordinates of lepton and b-quark
      auto sphi_lb =  hadTop_b_Rest.Phi() + lepTop_lep_Rest.Phi();
      if(sphi_lb > pie) sphi_lb = sphi_lb - 2.*pie;
      if(sphi_lb < -1.*pie) sphi_lb = sphi_lb + 2.*pie;
      event.set(h_sphi_lb, sphi_lb);
      if(debug) cout<<" SigmaPhi of lepton and b-quark is: "<<sphi_lb<<endl;

      // Define dphi_lb from phi coordinates of lepton and b-quark
      auto dphi_lb =  hadTop_b_Rest.Phi() - lepTop_lep_Rest.Phi();
      if(dphi_lb > pie) dphi_lb = dphi_lb - 2.*pie;
      if(dphi_lb < -1.*pie) dphi_lb = dphi_lb + 2.*pie;
      event.set(h_dphi_lb, dphi_lb);
      if(debug) cout<<" DeltaPhi of lepton and b-quark is: "<<dphi_lb<<endl;

      // Plot dphi and sphi for high-pt ranges
      if(Gen_HadTop.pt() > pt_hadTop_thresh){
        event.set(h_sphi_lb_high, sphi_lb);
        event.set(h_dphi_lb_high, dphi_lb);
      }
      // Plot dphi and sphi for low-pt ranges
      if(Gen_HadTop.pt() < pt_hadTop_thresh){
        event.set(h_sphi_lb_low, sphi_lb);
        event.set(h_dphi_lb_low, dphi_lb);
      }

    }

  }

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule_GenLvl_SpinCorrelationStudy)
