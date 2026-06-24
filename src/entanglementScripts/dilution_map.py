from __future__ import print_function
import ROOT
from array import array
from ttbargen_py import TTbarGenPy
ROOT.gROOT.SetBatch(True); ROOT.TH1.AddDirectory(False)

base="/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/"
files=[base+"uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root"%s for s in ("electron","electron2","muon","muon2")]
NMAX=-1         # test with e.g. 200000 first; -1 = full (needed for per-cell boosted stats)
MULT=5.0

def tlv(g):
    q=g.v4(); v=ROOT.TLorentzVector(); v.SetPtEtaPhiE(q.pt(),q.eta(),q.phi(),q.energy()); return v
def spincorr(ttg):
    lep=ttg.ChargedLepton(); qpos=(lep.pdgId()<0)
    Pos=tlv(ttg.Top()); Neg=tlv(ttg.Antitop()); lepv=tlv(lep); bv=tlv(ttg.BHad())
    tt=Pos+Neg; M=tt.M(); bcm=-tt.BoostVector()
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
    return lu.Dot(bu), c1[0]*c2[0]+c1[1]*c2[1]-c1[2]*c2[2], c, M

cosb=array('d',[-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1])
mb  =array('d',[800,900,1000,1200,1500,2000])       # boosted-focused y bins (gen m_tt)
p_reco=ROOT.TProfile2D("Dt_reco_boosted","reco #tilde{D};cos#Theta_{gen};M_{tt}^{gen} [GeV];#tilde{D}",10,cosb,5,mb)
p_gen =ROOT.TProfile2D("Dt_gen_boosted", "gen #tilde{D};cos#Theta_{gen};M_{tt}^{gen} [GeV];#tilde{D}",10,cosb,5,mb)

iW=irS=0.; r8W=r8rS=0.; g8W=g8rS=g8gS=g8rD=g8gD=0.
for fn in files:
    f=ROOT.TFile.Open(fn); t=f.Get("AnalysisTree")
    t.SetBranchStatus("*",0)
    for b in ("GenParticles*","cHel","cHel_P3n","M_tt","eventweight"): t.SetBranchStatus(b,1)
    nn=t.GetEntries(); nn=nn if NMAX<0 else min(nn,NMAX)
    for i in range(nn):
        t.GetEntry(i)
        if t.cHel_P3n<=-1.5: continue
        ttg=TTbarGenPy(t.GenParticles,False)
        if not ttg.IsSemiLeptonicDecay(): continue
        res=spincorr(ttg)
        if res is None: continue
        gC,gCp,gcos,gM = res; w=t.eventweight
        rsp=1.0 if t.cHel_P3n>=0 else -1.0; gsp=1.0 if gCp>=0 else -1.0
        iW+=w; irS+=w*rsp                                   # inclusive reco (validation)
        if t.M_tt>800: r8W+=w; r8rS+=w*rsp                  # reco Mtt>800 (validation vs 0.1539)
        if gM>800:                                          # truth-boosted yardstick + map
            g8W+=w; g8rS+=w*rsp; g8gS+=w*gsp
            g8rD+=w*(1 if t.cHel>=0 else -1); g8gD+=w*(1 if gC>=0 else -1)
            p_reco.Fill(gcos,gM,MULT*rsp,w); p_gen.Fill(gcos,gM,MULT*gsp,w)
    f.Close()

D=lambda s,w: MULT*s/w if w else 0.
print("\n--- validation: offline reco branches vs analysis extractCoeff ---")
print("  reco Dtilde inclusive    = %.4f   (analysis cHel_P3n          = 0.3519)" % D(irS,iW))
print("  reco Dtilde reco_Mtt>800 = %.4f   (analysis cHel_P3n_Mtt800   = 0.1539)" % D(r8rS,r8W))
print("\n--- boosted yardstick (gen_Mtt>800, eventweighted, identical events) ---")
print("  Dtilde:  reco=%.4f   gen=%.4f   reco/gen=%.2f" % (D(g8rS,g8W),D(g8gS,g8W), (D(g8rS,g8W)/D(g8gS,g8W) if g8gS else 0)))
print("  D     :  reco=%.4f   gen=%.4f   (sanity)" % (D(g8rD,g8W),D(g8gD,g8W)))

hr=p_reco.ProjectionXY("Dt_reco_xy"); hg=p_gen.ProjectionXY("Dt_gen_xy")
hd=hr.Clone("Dt_ratio_boosted"); hd.Divide(hg)              # per-cell reco/gen dilution
fo=ROOT.TFile("dilution_map_boosted.root","RECREATE"); p_reco.Write(); p_gen.Write(); hr.Write(); hg.Write(); hd.Write(); fo.Close()
print("wrote dilution_map_boosted.root")