from __future__ import print_function
import ROOT, math
from ttbargen_py import TTbarGenPy
ROOT.gROOT.SetBatch(True)
base="/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/"
files=[base+"uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root"%s for s in ("electron","electron2","muon","muon2")]
NMAX=50000

def tlv(g):
    q=g.v4(); v=ROOT.TLorentzVector(); v.SetPtEtaPhiE(q.pt(),q.eta(),q.phi(),q.energy()); return v
def comps(ttg):
    lep=ttg.ChargedLepton(); qpos=(lep.pdgId()<0)
    Pos=tlv(ttg.Top()); Neg=tlv(ttg.Antitop()); lepv=tlv(lep); bv=tlv(ttg.BHad())
    tt=Pos+Neg; bcm=-tt.BoostVector()
    for v in (Pos,Neg,lepv,bv): v.Boost(bcm)
    beam=ROOT.TVector3(0,0,1); k=Pos.Vect().Unit(); c=k.Dot(beam)
    if (1-c*c)<=1e-12: return None
    r=(beam-k*c).Unit(); n=(beam.Cross(k)).Unit(); s=1.0 if c>0 else -1.0
    kb,rb,nb=k,r*s,n*s
    lr=ROOT.TLorentzVector(lepv); br=ROOT.TLorentzVector(bv)
    if qpos: lr.Boost(-Pos.BoostVector()); br.Boost(-Neg.BoostVector())
    else:    lr.Boost(-Neg.BoostVector()); br.Boost(-Pos.BoostVector())
    lu,bu=lr.Vect().Unit(),br.Vect().Unit()
    if qpos: a1,a2=lu,bu
    else:    a1,a2=bu,lu
    return {"cosTheta1k":a1.Dot(kb),"cosTheta1r":a1.Dot(rb),"cosTheta1n":a1.Dot(nb),
            "cosTheta2k":a2.Dot(kb),"cosTheta2r":a2.Dot(rb),"cosTheta2n":a2.Dot(nb)}

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