from __future__ import print_function
import ROOT
from array import array
from ttbargen_py import TTbarGenPy

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# ---------------- config ----------------
base  = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/"
files = [base + "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_%s.root" % s
         for s in ("electron", "electron2", "muon", "muon2")]
out   = "gen_spincorr_TTToSemiLeptonic_all.root"
NMAX  = 50000        # <-- TEST first with a small number; set to -1 for the full run
MULT  = 5.0
# ----------------------------------------

def tlv(genp):
    q = genp.v4(); v = ROOT.TLorentzVector()
    v.SetPtEtaPhiE(q.pt(), q.eta(), q.phi(), q.energy()); return v

def spincorr(ttg):
    lep = ttg.ChargedLepton()
    qpos = (lep.pdgId() < 0)
    PosTop = tlv(ttg.Top()); NegTop = tlv(ttg.Antitop())
    lepv, bv = tlv(lep), tlv(ttg.BHad())
    ttbar = PosTop + NegTop; Mtt = ttbar.M()
    bcm = -ttbar.BoostVector()
    for v in (PosTop, NegTop, lepv, bv): v.Boost(bcm)
    beam = ROOT.TVector3(0., 0., 1.)
    k_axis = PosTop.Vect().Unit()
    cosTS = k_axis.Dot(beam)
    if (1.0 - cosTS*cosTS) <= 1e-12: return None
    r_axis = (beam - k_axis*cosTS).Unit()
    n_axis = (beam.Cross(k_axis)).Unit()
    sgn = 1.0 if cosTS > 0 else -1.0
    kbase, rbase, nbase = k_axis, r_axis*sgn, n_axis*sgn
    lep_rest, b_rest = ROOT.TLorentzVector(lepv), ROOT.TLorentzVector(bv)
    if qpos:
        lep_rest.Boost(-PosTop.BoostVector()); b_rest.Boost(-NegTop.BoostVector())
    else:
        lep_rest.Boost(-NegTop.BoostVector()); b_rest.Boost(-PosTop.BoostVector())
    lu, bu = lep_rest.Vect().Unit(), b_rest.Vect().Unit()
    if qpos:
        c1 = (lu.Dot(kbase), lu.Dot(rbase), lu.Dot(nbase))
        c2 = (bu.Dot(kbase), bu.Dot(rbase), bu.Dot(nbase))
    else:
        c1 = (bu.Dot(kbase), bu.Dot(rbase), bu.Dot(nbase))
        c2 = (lu.Dot(kbase), lu.Dot(rbase), lu.Dot(nbase))
    return lu.Dot(bu), c1[0]*c2[0]+c1[1]*c2[1]-c1[2]*c2[2], cosTS, Mtt

# ---------------- book ----------------
h_cHel    = ROOT.TH1F("cHel_gen",     "gen cos(#phi_{lb});cHel;events", 24, -1, 1)
h_cHelP3n = ROOT.TH1F("cHel_P3n_gen", "gen cos(#phi_{(P3n)lb});cHel_P3n;events", 24, -1, 1)
cosb = array('d', [-1.,-.8,-.6,-.4,-.2,0.,.2,.4,.6,.8,1.])
mb   = array('d', [200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,1000.,1100.,1200.,1300.,1400.,1500.,2000.])
p_cHel    = ROOT.TProfile2D("cHel_coeff_Mtt_vs_cosThetaStar_gen",     "gen D;cos#theta*;M_{tt} [GeV];coeff",      10, cosb, 22, mb)
p_cHelP3n = ROOT.TProfile2D("cHel_P3n_coeff_Mtt_vs_cosThetaStar_gen", "gen Dtilde;cos#theta*;M_{tt} [GeV];coeff", 10, cosb, 22, mb)

# ---------------- loop over files ----------------
tot_sl = 0; tot_w = tot_wS = tot_wSp = 0.0
for fn in files:
    name = fn.split("MC.")[-1].replace(".root", "")
    f = ROOT.TFile.Open(fn)
    if (not f) or f.IsZombie():
        print("SKIP (cannot open):", name); continue
    t = f.Get("AnalysisTree")
    t.SetBranchStatus("*", 0)
    t.SetBranchStatus("GenParticles*", 1)
    t.SetBranchStatus("eventweight", 1)
    N = t.GetEntries(); N = N if NMAX < 0 else min(N, NMAX)
    f_sl = 0; f_w = f_wS = f_wSp = 0.0
    for i in range(N):
        t.GetEntry(i)
        w = t.eventweight
        ttg = TTbarGenPy(t.GenParticles, False)
        if not ttg.IsSemiLeptonicDecay(): continue
        res = spincorr(ttg)
        if res is None: continue
        CHel, CHelP3n, cosTS, Mtt = res
        sS  = 1.0 if CHel    >= 0 else -1.0
        sSp = 1.0 if CHelP3n >= 0 else -1.0
        h_cHel.Fill(CHel, w); h_cHelP3n.Fill(CHelP3n, w)
        p_cHel.Fill(cosTS, Mtt, MULT*sS, w)
        p_cHelP3n.Fill(cosTS, Mtt, MULT*sSp, w)
        f_sl += 1; f_w += w; f_wS += w*sS; f_wSp += w*sSp
        if i % 500000 == 0: print("  %-12s %d/%d" % (name, i, N))
    f.Close()
    dd  = MULT*f_wS/f_w  if f_w else 0.0
    ddt = MULT*f_wSp/f_w if f_w else 0.0
    print("%-30s SL=%d  D=%.4f  Dtilde=%.4f" % (name, f_sl, dd, ddt))
    tot_sl += f_sl; tot_w += f_w; tot_wS += f_wS; tot_wSp += f_wSp

print("\nCOMBINED  SL=%d" % tot_sl)
print("  gen D      = %.4f   (muon-only reco ref: 0.1271)" % (MULT*tot_wS /tot_w if tot_w else 0))
print("  gen Dtilde = %.4f   (muon-only reco ref: 0.3645)" % (MULT*tot_wSp/tot_w if tot_w else 0))

fout = ROOT.TFile(out, "RECREATE")
for hobj in (h_cHel, h_cHelP3n, p_cHel, p_cHelP3n): hobj.Write()
fout.Close()
print("wrote", out)