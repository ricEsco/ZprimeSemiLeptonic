from __future__ import print_function
import ROOT
from ttbargen_py import TTbarGenPy
ROOT.gROOT.SetBatch(True)

base  = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/"
files = [base + "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root" % s for s in ("electron","electron2","muon","muon2")]
NMAX = 50000; MULT = 5.0

def tlv(g):
    q=g.v4(); v=ROOT.TLorentzVector(); v.SetPtEtaPhiE(q.pt(),q.eta(),q.phi(),q.energy()); return v
def spincorr(ttg):
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
    if qpos: c1=(lu.Dot(kb),lu.Dot(rb),lu.Dot(nb)); c2=(bu.Dot(kb),bu.Dot(rb),bu.Dot(nb))
    else:    c1=(bu.Dot(kb),bu.Dot(rb),bu.Dot(nb)); c2=(lu.Dot(kb),lu.Dot(rb),lu.Dot(nb))
    return lu.Dot(bu), c1[0]*c2[0]+c1[1]*c2[1]-c1[2]*c2[2]

medges=[0,350,400,450,500,600,700,800,1000,1500,1e9]; nbk=len(medges)-1
def mbin(m):
    for k in range(nbk):
        if medges[k]<=m<medges[k+1]: return k
    return nbk-1
W=[0.]*nbk; rS=[0.]*nbk; gS=[0.]*nbk; rSp=[0.]*nbk; gSp=[0.]*nbk

for fn in files:
    f=ROOT.TFile.Open(fn); t=f.Get("AnalysisTree")
    t.SetBranchStatus("*",0)
    for b in ("GenParticles*","cHel","cHel_P3n","M_tt","eventweight"): t.SetBranchStatus(b,1)
    n=t.GetEntries(); n=n if NMAX<0 else min(n,NMAX)
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