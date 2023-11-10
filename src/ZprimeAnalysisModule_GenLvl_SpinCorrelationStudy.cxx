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

  if(debug) cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug) cout << " run.event: " << event.run << ". " << event.event << endl;

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
  // event.set(h_sphi_lq_boost1, -10);
  // event.set(h_sphi_lq_boost2 , -10);
  // event.set(h_sphi_lq_boost3 , -10);
  // event.set(h_sphi_lq_boost4 , -10);
  
  // Difference of phi-coordinates between lepton and quark
  event.set(h_dphi_lq, -10);     
  event.set(h_dphi_lq_low, -10); 
  event.set(h_dphi_lq_high, -10);
  event.set(h_dphi_lq_Mass1, -10);
  event.set(h_dphi_lq_Mass2, -10);
  event.set(h_dphi_lq_Mass3, -10);
  event.set(h_dphi_lq_Mass4, -10);
  event.set(h_dphi_lq_Mass5, -10);
  // event.set(h_dphi_lq_boost1 , -10);
  // event.set(h_dphi_lq_boost2 , -10);
  // event.set(h_dphi_lq_boost3 , -10);
  // event.set(h_dphi_lq_boost4 , -10);

  // Sum of phi-coordinates between lepton and b-quark
  event.set(h_sphi_lb, -10);     
  event.set(h_sphi_lb_low, -10); 
  event.set(h_sphi_lb_high, -10);
  event.set(h_sphi_lb_Mass1, -10);
  event.set(h_sphi_lb_Mass2, -10);
  event.set(h_sphi_lb_Mass3, -10);
  event.set(h_sphi_lb_Mass4, -10);
  event.set(h_sphi_lb_Mass5, -10);
  // event.set(h_sphi_lb_boost1, -10);
  // event.set(h_sphi_lb_boost2 , -10);
  // event.set(h_sphi_lb_boost3 , -10);
  // event.set(h_sphi_lb_boost4 , -10);
  
  // Difference of phi-coordinates between lepton and b-quark
  event.set(h_dphi_lb, -10);     
  event.set(h_dphi_lb_low, -10); 
  event.set(h_dphi_lb_high, -10);
  event.set(h_dphi_lb_Mass1, -10);
  event.set(h_dphi_lb_Mass2, -10);
  event.set(h_dphi_lb_Mass3, -10);
  event.set(h_dphi_lb_Mass4, -10);
  event.set(h_dphi_lb_Mass5, -10);
  // event.set(h_dphi_lb_boost1 , -10);
  // event.set(h_dphi_lb_boost2 , -10);
  // event.set(h_dphi_lb_boost3 , -10);
  // event.set(h_dphi_lb_boost4 , -10);

  // fill ttbargen information
  ttgenprod->process(event);

  // Variable to access gen-particles
  TTbarGen ttbargen = event.get(h_ttbargen);

  // Make sure ttbar decays semileptonically
  if(ttbargen.IsSemiLeptonicDecay()){

    // Define cut-variable as pt of hadTop for low/high regions
    auto pt_hadTop_thresh = 150.;         

    //-------------------------------------------------- Start in LAB-frame --------------------------------------------------//
    // Plot pt of hadronic Top 
    LorentzVector Gen_HadTop = ttbargen.TopHad().v4();
    if(Gen_HadTop.pt() != -10) event.set(h_pt_hadTop, Gen_HadTop.pt());

    // Define 4vectors of decay products

    // b-quark
    LorentzVector Gen_b = ttbargen.BHad().v4(); 
    TLorentzVector hadTop_b;
    hadTop_b.SetPtEtaPhiE(Gen_b.pt(),Gen_b.eta(),Gen_b.phi(),Gen_b.energy());

    // Least-energetic W-quark will be defined AFTER WE BOOST TO MOTHER TOP'S REST-FRAME

    //  W-quark1
    LorentzVector Gen_q1 = ttbargen.Q1().v4(); 
    TLorentzVector hadTop_q1;
    hadTop_q1.SetPtEtaPhiE(Gen_q1.pt(),Gen_q1.eta(),Gen_q1.phi(),Gen_q1.energy());

    //  W-quark2
    LorentzVector Gen_q2 = ttbargen.Q2().v4(); 
    TLorentzVector hadTop_q2;
    hadTop_q2.SetPtEtaPhiE(Gen_q2.pt(),Gen_q2.eta(),Gen_q2.phi(),Gen_q2.energy());

    // lepton
    LorentzVector Gen_Lep = ttbargen.ChargedLepton().v4(); 
    TLorentzVector lepTop_lep;
    lepTop_lep.SetPtEtaPhiE(Gen_Lep.pt(), Gen_Lep.eta(), Gen_Lep.phi(), Gen_Lep.energy());

    // Top vectors
    TLorentzVector PosTop;
    PosTop.SetPtEtaPhiE(ttbargen.Top().v4().pt(), ttbargen.Top().v4().eta(), ttbargen.Top().v4().phi(), ttbargen.Top().v4().energy());
    TLorentzVector NegTop;
    NegTop.SetPtEtaPhiE(ttbargen.Antitop().v4().pt(), ttbargen.Antitop().v4().eta(), ttbargen.Antitop().v4().phi(), ttbargen.Antitop().v4().energy());

    // 4vector to represent ttbar system
    TLorentzVector ttbar(PosTop + NegTop);

    // Plotting mass of ttbar system
    event.set(h_ttbar_mass_LabFrame, ttbar.M());

    // Constructing longitudinal boost of ttbar system
    // auto boost = fabs(PosTop.Pz() + NegTop.Pz())/(PosTop.E() + NegTop.E());

    // Plotting longitudinal boost of ttbar system
    // event.set(h_ttbar_boost_LabFrame, boost);

    //---------------------------------------------------------- Boost into CoM-frame ----------------------------------------------------------//
    // Boost into ttbar Center of Momentum configuration 
    lepTop_lep.Boost(-1.*ttbar.BoostVector());
    hadTop_b.Boost(-1.*ttbar.BoostVector());
    hadTop_q1.Boost(-1.*ttbar.BoostVector());
    hadTop_q2.Boost(-1.*ttbar.BoostVector());
    PosTop.Boost(-1.*ttbar.BoostVector());
    NegTop.Boost(-1.*ttbar.BoostVector());

    // Cross check back-to-back configuration of ttbar system
    if(debug) cout<<" Angle between tops is: "<< PosTop.Angle(NegTop.Vect())<<endl;
    //---------------------------------------------------------- Boost into CoM-frame ----------------------------------------------------------//

    //----- Start Rotation into Helicity Frame -----//
    // Rotate tops and decay products about beam-line
    lepTop_lep.RotateZ(-1.*PosTop.Phi());
    hadTop_b.RotateZ(-1.*PosTop.Phi());
    hadTop_q1.RotateZ(-1.*PosTop.Phi());
    hadTop_q2.RotateZ(-1.*PosTop.Phi());
    PosTop.RotateZ(-1.*PosTop.Phi());
    NegTop.RotateZ(-1.*PosTop.Phi());
    // Rotate tops and decay products about y-zxis
    lepTop_lep.RotateY(-1.*PosTop.Theta());
    hadTop_b.RotateY(-1.*PosTop.Theta());
    hadTop_q1.RotateY(-1.*PosTop.Theta());
    hadTop_q2.RotateY(-1.*PosTop.Theta());
    PosTop.RotateY(-1.*PosTop.Theta());
    NegTop.RotateY(-1.*PosTop.Theta());
    //----- End Rotation into Helicity Frame -----//

    //--------------------------- Boost into ttbar rest-frame ---------------------------//
    // Boost the top's to rest individually, bringing their children with them
    // Decay products get boosted in opposite directions depending on their mother top

    // POSITIVE LEPTON CONFIGURATION
    if(ttbargen.ChargedLepton().charge() > 0){ // Positively charged lepton means
      lepTop_lep.Boost(-1.*PosTop.BoostVector()); // lepton has Positive Top mother
      hadTop_b.Boost(-1.*NegTop.BoostVector());   // b-jet has Negative Top mother
      hadTop_q1.Boost(-1.*NegTop.BoostVector());
      hadTop_q2.Boost(-1.*NegTop.BoostVector());
    }
    // NEGATIVE LEPTON CONFIGURATION
    if(ttbargen.ChargedLepton().charge() < 0){ // Negatively charged lepton means
      lepTop_lep.Boost(-1.*NegTop.BoostVector()); // lepton has Negative Top mother
      hadTop_b.Boost(-1.*PosTop.BoostVector());   // b-jet has Positive Top mother
      hadTop_q1.Boost(-1.*PosTop.BoostVector());
      hadTop_q2.Boost(-1.*PosTop.BoostVector());
    }
    //--------------------------- Boost into ttbar rest-frame ---------------------------//

    // Define least-energetic W-quark
    TLorentzVector hadTop_qlow;

    // Set 4vector of less-energetic W-daughter 
    if(hadTop_q1.E() < hadTop_q2.E()){hadTop_qlow = hadTop_q1;}
    else{hadTop_qlow = hadTop_q2;}

    // Plot phi-coordinates from the boosted 4vectors
    event.set(h_phi_lep, lepTop_lep.Phi());
    event.set(h_phi_b, hadTop_b.Phi());
    event.set(h_phi_qlow, hadTop_qlow.Phi());

    // sphi and dphi are both defined as: PosTop_phi +- NegTop_phi
    auto pie = TMath::Pi();

    // Positive lepton means PosTop has leptonic decay
    if(ttbargen.ChargedLepton().charge() > 0){
      // Define sphi_lq from phi coordinates of lepton and qlow
      auto sphi_lq = lepTop_lep.Phi() + hadTop_qlow.Phi();
      if(sphi_lq > pie) sphi_lq = sphi_lq - 2.*pie;
      if(sphi_lq < -1.*pie) sphi_lq = sphi_lq + 2.*pie;
      event.set(h_sphi_lq, sphi_lq);
      // Define dphi_lq from phi coordinates of lepton and qlow
      auto dphi_lq = lepTop_lep.Phi() - hadTop_qlow.Phi();
      if(dphi_lq > pie) dphi_lq = dphi_lq - 2.*pie;
      if(dphi_lq < -1.*pie) dphi_lq = dphi_lq + 2.*pie;
      event.set(h_dphi_lq, dphi_lq);

      // Define sphi_lb from phi coordinates of lepton and b-quark
      auto sphi_lb = lepTop_lep.Phi() + hadTop_b.Phi();
      if(sphi_lb > pie) sphi_lb = sphi_lb - 2.*pie;
      if(sphi_lb < -1.*pie) sphi_lb = sphi_lb + 2.*pie;
      event.set(h_sphi_lb, sphi_lb);
      // Define dphi_lb from phi coordinates of lepton and b-quark
      auto dphi_lb = lepTop_lep.Phi() - hadTop_b.Phi();
      if(dphi_lb > pie) dphi_lb = dphi_lb - 2.*pie;
      if(dphi_lb < -1.*pie) dphi_lb = dphi_lb + 2.*pie;
      event.set(h_dphi_lb, dphi_lb);

      // Plot dphi and sphi for high-pt ranges
      if(Gen_HadTop.pt() > pt_hadTop_thresh){
        event.set(h_sphi_lq_high, sphi_lq);
        event.set(h_dphi_lq_high, dphi_lq);
        event.set(h_sphi_lb_high, sphi_lb);
        event.set(h_dphi_lb_high, dphi_lb);
      }
      // Plot dphi_lq and sphi_lq for low-pt ranges
      if(Gen_HadTop.pt() < pt_hadTop_thresh){
        event.set(h_sphi_lq_low, sphi_lq);
        event.set(h_dphi_lq_low, dphi_lq);
        event.set(h_sphi_lb_low, sphi_lb);
        event.set(h_dphi_lb_low, dphi_lb);
      }

      // Plot dphi and sphi for various ranges of ttbar mass
      // 0 < Mass1 < 500
      if(ttbar.M() < 500.){
        event.set(h_sphi_lq_Mass1, sphi_lq);
        event.set(h_dphi_lq_Mass1, dphi_lq);
        event.set(h_sphi_lb_Mass1, sphi_lb);
        event.set(h_dphi_lb_Mass1, dphi_lb);
      }
      // 500 < Mass2 < 750
      if(ttbar.M() > 500. &&  ttbar.M() < 750.){
        event.set(h_sphi_lq_Mass2, sphi_lq);
        event.set(h_dphi_lq_Mass2, dphi_lq);
        event.set(h_sphi_lb_Mass2, sphi_lb);
        event.set(h_dphi_lb_Mass2, dphi_lb);
      }
      // 750 < Mass3 < 1000
      if(ttbar.M() > 750. &&  ttbar.M() < 1000.){
        event.set(h_sphi_lq_Mass3, sphi_lq);
        event.set(h_dphi_lq_Mass3, dphi_lq);
        event.set(h_sphi_lb_Mass3, sphi_lb);
        event.set(h_dphi_lb_Mass3, dphi_lb);
      }
      // 1000 < Mass4 < 1500
      if(ttbar.M() > 1000. && ttbar.M() < 1500.){
        event.set(h_sphi_lq_Mass4, sphi_lq);
        event.set(h_dphi_lq_Mass4, dphi_lq);
        event.set(h_sphi_lb_Mass4, sphi_lb);
        event.set(h_dphi_lb_Mass4, dphi_lb);
      }
      // Mass5 > 1500
      if(ttbar.M() > 1500.){
        event.set(h_sphi_lq_Mass5, sphi_lq);
        event.set(h_dphi_lq_Mass5, dphi_lq);
        event.set(h_sphi_lb_Mass5, sphi_lb);
        event.set(h_dphi_lb_Mass5, dphi_lb);
      }

      // // Plot dphi and sphi for various ranges of longitudinal boost of ttbar system
      // // 0 < boost1 < 0.3
      // if(boost < 0.3){
      //   event.set(h_sphi_lq_boost1, sphi_lq);
      //   event.set(h_dphi_lq_boost1, dphi_lq);
      //   event.set(h_sphi_lb_boost1, sphi_lb);
      //   event.set(h_dphi_lb_boost1, dphi_lb);
      // }
      // // 0.3 < boost2 < 0.6
      // if(boost > 0.3 &&  boost < 0.6){
      //   event.set(h_sphi_lq_boost2, sphi_lq);
      //   event.set(h_dphi_lq_boost2, dphi_lq);
      //   event.set(h_sphi_lb_boost2, sphi_lb);
      //   event.set(h_dphi_lb_boost2, dphi_lb);
      // }
      // // 0.6 < boost3 < 0.8
      // if(boost > 0.6 &&  boost < 0.8){
      //   event.set(h_sphi_lq_boost3, sphi_lq);
      //   event.set(h_dphi_lq_boost3, dphi_lq);
      //   event.set(h_sphi_lb_boost3, sphi_lb);
      //   event.set(h_dphi_lb_boost3, dphi_lb);
      // }
      // // 0.8 < boost4 < 1.0
      // if(boost > 0.8 &&  boost < 1.0){
      //   event.set(h_sphi_lq_boost4, sphi_lq);
      //   event.set(h_dphi_lq_boost4, dphi_lq);
      //   event.set(h_sphi_lb_boost4, sphi_lb);
      //   event.set(h_dphi_lb_boost4, dphi_lb);
      // }
    }

    // Negative lepton means PosTop has hadronic decay
    if(ttbargen.ChargedLepton().charge() < 0){
      // Define sphi_lq from phi coordinates of lepton and qlow
      auto sphi_lq =  hadTop_qlow.Phi() + lepTop_lep.Phi();
      if(sphi_lq > pie) sphi_lq = sphi_lq - 2.*pie;
      if(sphi_lq < -1.*pie) sphi_lq = sphi_lq + 2.*pie;
      event.set(h_sphi_lq, sphi_lq);
      // Define dphi_lq from phi coordinates of lepton and qlow
      auto dphi_lq =  hadTop_qlow.Phi() - lepTop_lep.Phi();
      if(dphi_lq > pie) dphi_lq = dphi_lq - 2.*pie;
      if(dphi_lq < -1.*pie) dphi_lq = dphi_lq + 2.*pie;
      event.set(h_dphi_lq, dphi_lq);

      // Define sphi_lb from phi coordinates of lepton and b-quark
      auto sphi_lb =  hadTop_b.Phi() + lepTop_lep.Phi();
      if(sphi_lb > pie) sphi_lb = sphi_lb - 2.*pie;
      if(sphi_lb < -1.*pie) sphi_lb = sphi_lb + 2.*pie;
      event.set(h_sphi_lb, sphi_lb);
      // Define dphi_lb from phi coordinates of lepton and b-quark
      auto dphi_lb =  hadTop_b.Phi() - lepTop_lep.Phi();
      if(dphi_lb > pie) dphi_lb = dphi_lb - 2.*pie;
      if(dphi_lb < -1.*pie) dphi_lb = dphi_lb + 2.*pie;
      event.set(h_dphi_lb, dphi_lb);

      // Plot dphi and sphi for high-pt ranges
      if(Gen_HadTop.pt() > pt_hadTop_thresh){
        event.set(h_sphi_lq_high, sphi_lq);
        event.set(h_dphi_lq_high, dphi_lq);
        event.set(h_sphi_lb_high, sphi_lb);
        event.set(h_dphi_lb_high, dphi_lb);
      }
      // Plot dphi_lq and sphi_lq for low-pt ranges
      if(Gen_HadTop.pt() < pt_hadTop_thresh){
        event.set(h_sphi_lq_low, sphi_lq);
        event.set(h_dphi_lq_low, dphi_lq);
        event.set(h_sphi_lb_low, sphi_lb);
        event.set(h_dphi_lb_low, dphi_lb);
      }


      // Plot dphi and sphi for various ranges of ttbar mass
      // 0 < Mass1 < 500
      if(ttbar.M() < 500.){
        event.set(h_sphi_lq_Mass1, sphi_lq);
        event.set(h_dphi_lq_Mass1, dphi_lq);
        event.set(h_sphi_lb_Mass1, sphi_lb);
        event.set(h_dphi_lb_Mass1, dphi_lb);
      }
      // 500 < Mass2 < 750
      if(ttbar.M() > 500. &&  ttbar.M() < 750.){
        event.set(h_sphi_lq_Mass2, sphi_lq);
        event.set(h_dphi_lq_Mass2, dphi_lq);
        event.set(h_sphi_lb_Mass2, sphi_lb);
        event.set(h_dphi_lb_Mass2, dphi_lb);
      }
      // 750 < Mass3 < 1000
      if(ttbar.M() > 750. &&  ttbar.M() < 1000.){
        event.set(h_sphi_lq_Mass3, sphi_lq);
        event.set(h_dphi_lq_Mass3, dphi_lq);
        event.set(h_sphi_lb_Mass3, sphi_lb);
        event.set(h_dphi_lb_Mass3, dphi_lb);
      }
      // 1000 < Mass4 < 1500
      if(ttbar.M() > 1000. && ttbar.M() < 1500.){
        event.set(h_sphi_lq_Mass4, sphi_lq);
        event.set(h_dphi_lq_Mass4, dphi_lq);
        event.set(h_sphi_lb_Mass4, sphi_lb);
        event.set(h_dphi_lb_Mass4, dphi_lb);
      }
      // Mass5 > 1500
      if(ttbar.M() > 1500.){
        event.set(h_sphi_lq_Mass5, sphi_lq);
        event.set(h_dphi_lq_Mass5, dphi_lq);
        event.set(h_sphi_lb_Mass5, sphi_lb);
        event.set(h_dphi_lb_Mass5, dphi_lb);
      }


      // // Plot dphi and sphi for various ranges of longitudinal boost of ttbar system
      // // 0 < boost1 < 0.3
      // if(boost < 0.3){
      //   event.set(h_sphi_lq_boost1, sphi_lq);
      //   event.set(h_dphi_lq_boost1, dphi_lq);
      //   event.set(h_sphi_lb_boost1, sphi_lb);
      //   event.set(h_dphi_lb_boost1, dphi_lb);
      // }
      // // 0.3 < boost2 < 0.6
      // if(boost > 0.3 &&  boost < 0.6){
      //   event.set(h_sphi_lq_boost2, sphi_lq);
      //   event.set(h_dphi_lq_boost2, dphi_lq);
      //   event.set(h_sphi_lb_boost2, sphi_lb);
      //   event.set(h_dphi_lb_boost2, dphi_lb);
      // }
      // // 0.6 < boost3 < 0.8
      // if(boost > 0.6 &&  boost < 0.8){
      //   event.set(h_sphi_lq_boost3, sphi_lq);
      //   event.set(h_dphi_lq_boost3, dphi_lq);
      //   event.set(h_sphi_lb_boost3, sphi_lb);
      //   event.set(h_dphi_lb_boost3, dphi_lb);
      // }
      // // 0.8 < boost4 < 1.0
      // if(boost > 0.8 &&  boost < 1.0){
      //   event.set(h_sphi_lq_boost4, sphi_lq);
      //   event.set(h_dphi_lq_boost4, dphi_lq);
      //   event.set(h_sphi_lb_boost4, sphi_lb);
      //   event.set(h_dphi_lb_boost4, dphi_lb);
      // }  
    }

  }

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule_GenLvl_SpinCorrelationStudy)
