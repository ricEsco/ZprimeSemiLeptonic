from __future__ import print_function
import ROOT, math
from ttbargen_py import TTbarGenPy
ROOT.gROOT.SetBatch(True)
base="/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/"
files=[base+"uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root"%s for s in ("electron","electron2","muon","muon2")]
NMAX=50000

# Extract the TLorentzVector from a GenParticle
def tlv(genp):
    q = genp.v4(); v = ROOT.TLorentzVector()
    v.SetPtEtaPhiE(q.pt(), q.eta(), q.phi(), q.energy()); return v

# Compute SC components at generator level for a semi-leptonic ttbar decay, using a single Bernreuther basis for top quark
def comps(ttg):
    # lab frame
    lep = ttg.ChargedLepton()                             # lepton
    qpos = (lep.pdgId() < 0)                              # lepton charge boolean
    PosTop = tlv(ttg.Top()); NegTop = tlv(ttg.Antitop())  # top and antitop 4vectors
    lepv, bv = tlv(lep), tlv(ttg.BHad())                  # lepton and b-jet 4vectors
    ttbar = PosTop + NegTop; Mtt = ttbar.M()              # ttbar 4vector and invariant mass
    bcm = -ttbar.BoostVector()                            # boost vector = -ttbar
    # ttbar rest frame
    for v in (PosTop, NegTop, lepv, bv): v.Boost(bcm)     # boost all 4vectors to ttbar rest frame

    # build Bernreuther basis
    p_beam = ROOT.TVector3(0., 0., 1.)                    # p_beam direction unit-3vector = p
    k_axis = PosTop.Vect().Unit()                         # direction of top quark in ttbar rest frame unit-3vector = k
    cosTS = k_axis.Dot(p_beam)                            # top quark scattering angle = y
    if (1.0 - cosTS*cosTS) <= 1e-12: return None          # veto if cosTS is too close to 1 or -1 (numerical issues)
    r_axis = (p_beam - k_axis*cosTS).Unit()               # Bernreuther r-axis unit-3vector: (1/|r|)*(p - yk)
    n_axis = (p_beam.Cross(k_axis)).Unit()                # Bernreuther n-axis unit-3vector: (1/|r|)*(p x k)
    sgn = 1.0 if cosTS > 0 else -1.0                      # sign of top scattering angle
    kbase, rbase, nbase = k_axis, r_axis*sgn, n_axis*sgn  # reference axes used in projections (only top quark reference axes are used for both spin analyzers)

    # build SC variables
    lep_rest, b_rest = ROOT.TLorentzVector(lepv), ROOT.TLorentzVector(bv)
    # boost spin analyzers to their parent top quark rest frame (based on lepton charge)
    if qpos:
        lep_rest.Boost(-PosTop.BoostVector()); b_rest.Boost(-NegTop.BoostVector())
    else:
        lep_rest.Boost(-NegTop.BoostVector()); b_rest.Boost(-PosTop.BoostVector())
    # define directions of spin analyzers (unit-3vectors) in their parent top quark rest frame
    lu, bu = lep_rest.Vect().Unit(), b_rest.Vect().Unit()
    if qpos: a1,a2=lu,bu
    else:    a1,a2=bu,lu
    return {"cosTheta1k":a1.Dot(kbase),"cosTheta1r":a1.Dot(rbase),"cosTheta1n":a1.Dot(nbase),
            "cosTheta2k":a2.Dot(kbase),"cosTheta2r":a2.Dot(rbase),"cosTheta2n":a2.Dot(nbase)}

names=["cosTheta1k","cosTheta1r","cosTheta1n","cosTheta2k","cosTheta2r","cosTheta2n"]
S={c:[0.,0.,0.,0.,0.,0] for c in names}
for fn in files:
    f=ROOT.TFile.Open(fn); t=f.Get("AnalysisTree")
    t.SetBranchStatus("*",0)
    for b in ["GenParticles*","cHel"]+names: t.SetBranchStatus(b,1)
    nn=t.GetEntries(); nn=nn if NMAX<0 else min(nn,NMAX)
    for i in range(nn):
        t.GetEntry(i)
        if t.cHel<=-1.5: continue
        ttg=TTbarGenPy(t.GenParticles,False)
        if not ttg.IsSemiLeptonicDecay(): continue
        g=comps(ttg)
        if g is None: continue
        for c in names:
            x=getattr(t,c)
            if x<=-1.5: continue
            y=g[c]; s=S[c]; s[0]+=x; s[1]+=y; s[2]+=x*x; s[3]+=y*y; s[4]+=x*y; s[5]+=1
    f.Close()

print("\n  component     reco_mean  gen_mean   corr")
for c in names:
    Sx,Sy,Sxx,Syy,Sxy,N=S[c]
    if N==0: continue
    cov=Sxy/N-(Sx/N)*(Sy/N); vx=Sxx/N-(Sx/N)**2; vy=Syy/N-(Sy/N)**2
    corr=cov/math.sqrt(vx*vy) if vx>0 and vy>0 else 0.
    print("  %-12s  %+.4f    %+.4f   %+.3f" % (c, Sx/N, Sy/N, corr))