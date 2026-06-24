from __future__ import print_function
import ROOT
from ttbargen_py import TTbarGenPy

ROOT.gROOT.SetBatch(True); ROOT.TH1.AddDirectory(False)
base  = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/"
files = [base + "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root" % s for s in ("electron","electron2","muon","muon2")]
NMAX = 50000          # test; -1 for full
MULT = 5.0

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

h2   = ROOT.TH2F("cHel_reco_vs_gen","cHel reco vs gen;reco;gen",50,-1,1,50,-1,1)
h2p  = ROOT.TH2F("cHelP3n_reco_vs_gen","cHel_P3n reco vs gen;reco;gen",50,-1,1,50,-1,1)
N_=0; W=0.0; rS=gS=rSp=gSp=0.0; agree=agreeP=0.0
for fn in files:
    f=ROOT.TFile.Open(fn); t=f.Get("AnalysisTree")
    t.SetBranchStatus("*",0)
    for b in ("GenParticles*","cHel","cHel_P3n","M_tt","eventweight"): t.SetBranchStatus(b,1)
    n=t.GetEntries(); n=n if NMAX<0 else min(n,NMAX)
    for i in range(n):
        t.GetEntry(i)
        if t.cHel <= -1.5: continue                 # reco not reconstructed
        ttg=TTbarGenPy(t.GenParticles, False)
        if not ttg.IsSemiLeptonicDecay(): continue
        res=spincorr(ttg)
        if res is None: continue
        gC,gCp = res; rC,rCp = t.cHel, t.cHel_P3n; w=t.eventweight
        h2.Fill(rC,gC,w); h2p.Fill(rCp,gCp,w)
        N_+=1; W+=w
        rsd=1.0 if rC>=0 else -1.0; gsd=1.0 if gC>=0 else -1.0
        rsp=1.0 if rCp>=0 else -1.0; gsp=1.0 if gCp>=0 else -1.0
        rS+=w*rsd; gS+=w*gsd; rSp+=w*rsp; gSp+=w*gsp
        agree += w*(rsd==gsd); agreeP += w*(rsp==gsp)
    f.Close()

print("\nmatched reconstructed events: %d" % N_)
print("  D   : reco=%.4f  gen=%.4f   (sign-agree %.1f%%)  per-evt corr=%.3f" %
      (MULT*rS/W, MULT*gS/W, 100.*agree/W, h2.GetCorrelationFactor()))
print("  Dt  : reco=%.4f  gen=%.4f   (sign-agree %.1f%%)  per-evt corr=%.3f" %
      (MULT*rSp/W, MULT*gSp/W, 100.*agreeP/W, h2p.GetCorrelationFactor()))
fo=ROOT.TFile("recogen_diag.root","RECREATE"); h2.Write(); h2p.Write(); fo.Close()