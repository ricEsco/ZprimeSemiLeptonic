from __future__ import print_function
import ROOT
from ttbargen_py import TTbarGenPy
ROOT.gROOT.SetBatch(True)

'''
An attempt to diagnose reco bias in threshold and boosted entanglement proxies (D & Dtilde) 
by saving reco and gen proxies in bins of Mtt and cosTS in same TH2F
'''

# ---------------------------------------------- config ---------------------------------------------- #
base  = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/"
files = [base + "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root" % s 
         for s in ("electron","electron2","muon","muon2")]
NMAX = 50000; MULT = 5.0
# ---------------------------------------------------------------------------------------------------- #

# Extract the TLorentzVector from a GenParticle
def tlv(genp):
    q = genp.v4(); v = ROOT.TLorentzVector()
    v.SetPtEtaPhiE(q.pt(), q.eta(), q.phi(), q.energy()); return v

# Compute the spin correlation projections at generator level for a semi-leptonic ttbar decay, using a single Bernreuther basis for top quark
# return spin analyzer's opening angle, reflected-nbase opening angle
def spincorr(ttg):
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
    # compute projections of spin analyzers onto the top quark reference axes, based on lepton charge
    if qpos:
        c1 = (lu.Dot(kbase), lu.Dot(rbase), lu.Dot(nbase))
        c2 = (bu.Dot(kbase), bu.Dot(rbase), bu.Dot(nbase))
    else:
        c1 = (bu.Dot(kbase), bu.Dot(rbase), bu.Dot(nbase))
        c2 = (lu.Dot(kbase), lu.Dot(rbase), lu.Dot(nbase))

    return lu.Dot(bu), c1[0]*c2[0]+c1[1]*c2[1]-c1[2]*c2[2]


# histogram details
medges=[0,350,400,450,500,600,700,800,1000,1500,1e9] # histogram bin edges
nbk=len(medges)-1
def mbin(m):
    for k in range(nbk):
        if medges[k]<=m<medges[k+1]: return k
    return nbk-1
W=[0.]*nbk; rS=[0.]*nbk; gS=[0.]*nbk; rSp=[0.]*nbk; gSp=[0.]*nbk

# ---------------- loop over inFiles ----------------
for fn in files:
    f=ROOT.TFile.Open(fn); t=f.Get("AnalysisTree")
    # disable all branches except the ones we need (to speed up reading)
    t.SetBranchStatus("*",0)
    for b in ("GenParticles*","cHel","cHel_P3n","M_tt","eventweight"): t.SetBranchStatus(b,1)
    n=t.GetEntries(); n=n if NMAX<0 else min(n,NMAX)

    # begin loop over events
    for i in range(n):
        t.GetEntry(i)
        if t.cHel<=-1.5: continue
        ttg=TTbarGenPy(t.GenParticles,False)
        if not ttg.IsSemiLeptonicDecay(): continue
        res=spincorr(ttg)
        if res is None: continue
        gC,gCp=res; w=t.eventweight; k=mbin(t.M_tt)
        W[k]+=w
        rS[k]+=w*(1 if t.cHel>=0 else -1);   gS[k]+=w*(1 if gC>=0 else -1)
        rSp[k]+=w*(1 if t.cHel_P3n>=0 else -1); gSp[k]+=w*(1 if gCp>=0 else -1)
    f.Close()

print("\n  Mtt range       N_w        D_reco   D_gen    Dt_reco  Dt_gen")
tw=trS=tgS=trSp=tgSp=0.
for k in range(nbk):
    if W[k]<=0: continue
    hi = medges[k+1] if medges[k+1]<1e8 else 9999
    print("  [%4.0f,%5.0f]  %9.0f   %+0.3f   %+0.3f   %+0.3f   %+0.3f" % (
        medges[k],hi,W[k], MULT*rS[k]/W[k],MULT*gS[k]/W[k], MULT*rSp[k]/W[k],MULT*gSp[k]/W[k]))
    tw+=W[k]; trS+=rS[k]; tgS+=gS[k]; trSp+=rSp[k]; tgSp+=gSp[k]
print("  %-12s %9.0f   %+0.3f   %+0.3f   %+0.3f   %+0.3f" % (
    "ALL",tw, MULT*trS/tw,MULT*tgS/tw, MULT*trSp/tw,MULT*tgSp/tw))