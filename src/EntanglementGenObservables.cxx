#include <UHH2/ZprimeSemiLeptonic/include/EntanglementGenObservables.h>

#include <TVector3.h>

#include <cmath>

using namespace std;

namespace {

  // convert a UHH2 GenParticle to a ROOT TLorentzVector
  TLorentzVector to_tlv(const GenParticle& genp){
    TLorentzVector v;
    v.SetPtEtaPhiE(genp.v4().pt(), genp.v4().eta(), genp.v4().phi(), genp.v4().E());
    return v;
  }

  // wrap an angle back into (-pi, pi] with a single +/- 2pi correction
  double wrap_angle(double angle){
    if(angle >  M_PI) angle -= 2.0 * M_PI;
    if(angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
  }

}

namespace entanglement {

  const vector<string>& gen_observable_names(){
    static const vector<string> names = {
      // ttbar-system kinematics and lepton charge
      "M_tt", "beta", "absDeltaY", "cosThetaStar", "pt_hadTop", "Lepton_charge",

      // lepton-only projections (99 in the slot of the absent (anti)lepton) and their *Star versions
      "cosTheta1k_antiLep", "cosTheta1r_antiLep", "cosTheta1n_antiLep", "cosTheta1kStar_antiLep", "cosTheta1rStar_antiLep",
      "cosTheta2k_Lep",     "cosTheta2r_Lep",     "cosTheta2n_Lep",     "cosTheta2kStar_Lep",     "cosTheta2rStar_Lep",

      // spin-analyzer 1 (top side) / 2 (antitop side) projections and their *Star versions
      "cosTheta1k", "cosTheta1r", "cosTheta1n", "cosTheta1kStar", "cosTheta1rStar",
      "cosTheta2k", "cosTheta2r", "cosTheta2n", "cosTheta2kStar", "cosTheta2rStar",

      // correlation matrix elements and their sum/difference combinations
      "Cnn", "Cnr", "Cnk", "Crn", "Crr", "Crk", "Ckn", "Ckr", "Ckk",
      "Crk_plus", "Crk_minus", "Cnr_plus", "Cnr_minus", "Cnk_plus", "Cnk_minus",

      // entanglement witnesses and their windowed slices
      "cHel",     "cHel_Mtt300_400",     "cHel_Mtt300_400_betaLT0p9",
      "cHel_P3n", "cHel_P3n_Mtt800_Inf", "cHel_P3n_Mtt800_Inf_cosThetaLT0p4",

      // Baumgart et al. angular variables (inclusive and low/high pT(hadronic top))
      "Sigma_phi", "Sigma_phi_low", "Sigma_phi_high",
      "Delta_phi", "Delta_phi_low", "Delta_phi_high",
    };
    return names;
  }


  GenSpinCorrelations compute_gen_spin_correlations(const TTbarGen& ttbargen){

    GenSpinCorrelations result;
    for(const auto& name : gen_observable_names()) result.obs[name] = SENTINEL; // default: not available

    // Observables are only defined for e/mu/tau + jets decays
    if(!ttbargen.IsSemiLeptonicDecay()) return result;
    result.is_semileptonic = true;

    //----------------------------- Lab-frame 4-vectors -----------------------------//
    const GenParticle lepton = ttbargen.ChargedLepton();    // charged lepton from the W decay
    const bool lepton_is_positive = (lepton.pdgId() < 0);   // l+ (pdgId < 0) => the TOP decays leptonically

    TLorentzVector Top     = to_tlv(ttbargen.WTop())     + to_tlv(ttbargen.bTop());     // last-copy top     = W+ + b
    TLorentzVector Antitop = to_tlv(ttbargen.WAntitop()) + to_tlv(ttbargen.bAntitop()); // last-copy antitop = W- + bbar
    TLorentzVector Lep     = to_tlv(lepton);
    TLorentzVector bHad    = to_tlv(lepton_is_positive ? ttbargen.bAntitop() : ttbargen.bTop()); // b-quark on the hadronic side

    result.Top = Top; result.Antitop = Antitop; result.Lepton = Lep; result.bHad = bHad;

    // ttbar-system quantities in the lab frame
    const TLorentzVector ttbar = Top + Antitop;
    const double Mtt  = ttbar.M();                                                         // ttbar invariant mass
    const double Ett  = Top.E() + Antitop.E();
    const double beta = (Ett > 0.) ? fabs(Top.Pz() + Antitop.Pz()) / Ett : SENTINEL;       // longitudinal velocity of the ttbar system
    const double absDeltaY = fabs(Top.Rapidity()) - fabs(Antitop.Rapidity());              // |y_t| - |y_tbar|
    const double pt_hadTop = (lepton_is_positive ? Antitop : Top).Pt();                    // pT of the hadronically decaying top

    //----------------------------- Boost into the ttbar CM frame -----------------------------//
    const TVector3 CMboost = -ttbar.BoostVector();
    Top.Boost(CMboost); Antitop.Boost(CMboost); Lep.Boost(CMboost); bHad.Boost(CMboost);

    //----------------------------- Bernreuther (k, r, n) basis -----------------------------//
    const TVector3 beam(0., 0., 1.);                        // beam direction
    const TVector3 k = Top.Vect().Unit();                   // top direction in the CM frame
    const double cosTS = k.Dot(beam);                       // cosine of the top scattering angle (cos theta*)
    if((1. - cosTS*cosTS) <= 1e-12) return result;          // top parallel to the beam: basis undefined
    const TVector3 r = (beam - cosTS * k).Unit();           // in the production plane, orthogonal to k
    const TVector3 n = beam.Cross(k).Unit();                // normal to the production plane

    const double sign_cosTS     = (cosTS > 0.)     ? 1. : -1.; // Bose-symmetry factor
    const double sign_absDeltaY = (absDeltaY > 0.) ? 1. : -1.; // charge-asymmetry factor

    const TVector3 kbase = k;                                  // Bose-symmetrized basis
    const TVector3 rbase = sign_cosTS * r;
    const TVector3 nbase = sign_cosTS * n;
    const TVector3 kStar = sign_absDeltaY * k;                 // charge-asymmetry-corrected basis
    const TVector3 rStar = sign_absDeltaY * sign_cosTS * r;

    //----------------------------- Spin analyzers in their parent top rest frame -----------------------------//
    TLorentzVector Lep_rest  = Lep;
    TLorentzVector bHad_rest = bHad;
    if(lepton_is_positive){              // l+ comes from the top, the hadronic b from the antitop
      Lep_rest.Boost(-Top.BoostVector());
      bHad_rest.Boost(-Antitop.BoostVector());
    }
    else{                                // l- comes from the antitop, the hadronic b from the top
      Lep_rest.Boost(-Antitop.BoostVector());
      bHad_rest.Boost(-Top.BoostVector());
    }
    const TVector3 lep_dir = Lep_rest.Vect().Unit();
    const TVector3 b_dir   = bHad_rest.Vect().Unit();

    // spin analyzer 1 belongs to the top, spin analyzer 2 to the antitop
    const TVector3 SA1 = lepton_is_positive ? lep_dir : b_dir;
    const TVector3 SA2 = lepton_is_positive ? b_dir   : lep_dir;

    //----------------------------- Observables -----------------------------//
    auto& obs = result.obs;

    obs["M_tt"]          = Mtt;
    obs["beta"]          = beta;
    obs["absDeltaY"]     = absDeltaY;
    obs["cosThetaStar"]  = cosTS;
    obs["pt_hadTop"]     = pt_hadTop;
    obs["Lepton_charge"] = lepton_is_positive ? 1. : -1.;   // physical lepton charge (same convention as reco Lepton_charge)

    // lepton-only projections: anti-lepton -> top polarization (slot 1), lepton -> antitop polarization (slot 2)
    const string slot   = lepton_is_positive ? "1" : "2";
    const string suffix = lepton_is_positive ? "_antiLep" : "_Lep";
    const string other_slot   = lepton_is_positive ? "2" : "1";
    const string other_suffix = lepton_is_positive ? "_Lep" : "_antiLep";
    obs["cosTheta" + slot + "k"     + suffix] = lep_dir.Dot(kbase);
    obs["cosTheta" + slot + "r"     + suffix] = lep_dir.Dot(rbase);
    obs["cosTheta" + slot + "n"     + suffix] = lep_dir.Dot(nbase);
    obs["cosTheta" + slot + "kStar" + suffix] = lep_dir.Dot(kStar);
    obs["cosTheta" + slot + "rStar" + suffix] = lep_dir.Dot(rStar);
    for(const string comp : {"k", "r", "n", "kStar", "rStar"}) obs["cosTheta" + other_slot + comp + other_suffix] = NOT_PRESENT;

    // spin-analyzer projections
    const double c1k = SA1.Dot(kbase), c1r = SA1.Dot(rbase), c1n = SA1.Dot(nbase);
    const double c2k = SA2.Dot(kbase), c2r = SA2.Dot(rbase), c2n = SA2.Dot(nbase);
    obs["cosTheta1k"] = c1k; obs["cosTheta1r"] = c1r; obs["cosTheta1n"] = c1n;
    obs["cosTheta2k"] = c2k; obs["cosTheta2r"] = c2r; obs["cosTheta2n"] = c2n;
    obs["cosTheta1kStar"] = SA1.Dot(kStar); obs["cosTheta1rStar"] = SA1.Dot(rStar);
    obs["cosTheta2kStar"] = SA2.Dot(kStar); obs["cosTheta2rStar"] = SA2.Dot(rStar);

    // correlation matrix elements C_ij = cos(theta1_i) * cos(theta2_j)
    obs["Cnn"] = c1n*c2n; obs["Cnr"] = c1n*c2r; obs["Cnk"] = c1n*c2k;
    obs["Crn"] = c1r*c2n; obs["Crr"] = c1r*c2r; obs["Crk"] = c1r*c2k;
    obs["Ckn"] = c1k*c2n; obs["Ckr"] = c1k*c2r; obs["Ckk"] = c1k*c2k;
    obs["Crk_plus"] = c1r*c2k + c1k*c2r; obs["Crk_minus"] = c1r*c2k - c1k*c2r;
    obs["Cnr_plus"] = c1n*c2r + c1r*c2n; obs["Cnr_minus"] = c1n*c2r - c1r*c2n;
    obs["Cnk_plus"] = c1n*c2k + c1k*c2n; obs["Cnk_minus"] = c1n*c2k - c1k*c2n;

    // entanglement witnesses: D (threshold) and D-tilde (boosted, flipped n-component)
    const double cHel     = lep_dir.Dot(b_dir);
    const double cHel_P3n = c1k*c2k + c1r*c2r - c1n*c2n;
    obs["cHel"]     = cHel;
    obs["cHel_P3n"] = cHel_P3n;
    if(Mtt < 400.){                                   // near threshold
      obs["cHel_Mtt300_400"] = cHel;
      if(beta < 0.9) obs["cHel_Mtt300_400_betaLT0p9"] = cHel;
    }
    if(Mtt > 800.){                                   // boosted
      obs["cHel_P3n_Mtt800_Inf"] = cHel_P3n;
      if(fabs(cosTS) < 0.4) obs["cHel_P3n_Mtt800_Inf_cosThetaLT0p4"] = cHel_P3n; // boosted and central
    }

    // Baumgart et al. variables: azimuthal angles of the spin analyzers in the (n, r) plane
    const double lep_phi = atan2(lep_dir.Dot(rbase), lep_dir.Dot(nbase));
    const double b_phi   = atan2(b_dir.Dot(rbase),   b_dir.Dot(nbase));
    const double Sigma_phi = wrap_angle(lep_phi + b_phi);
    const double Delta_phi = lepton_is_positive ? wrap_angle(lep_phi - b_phi) : wrap_angle(b_phi - lep_phi); // (top-side) - (antitop-side)
    obs["Sigma_phi"] = Sigma_phi;
    obs["Delta_phi"] = Delta_phi;
    if(pt_hadTop < PT_HADTOP_THRESHOLD){ obs["Sigma_phi_low"]  = Sigma_phi; obs["Delta_phi_low"]  = Delta_phi; }
    if(pt_hadTop > PT_HADTOP_THRESHOLD){ obs["Sigma_phi_high"] = Sigma_phi; obs["Delta_phi_high"] = Delta_phi; }

    result.obs_valid = true;
    return result;
  }

}
