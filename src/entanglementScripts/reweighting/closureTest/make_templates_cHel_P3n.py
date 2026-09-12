# make_templates_cHel_P3n.py -- build the reco-level cHel_P3n (D-tilde) templates for Combine
# using the spin-correlation REWEIGHTING method (docs/reweighting_method.md), following the
# dilepton entanglement analysis arXiv:2406.03976.
#
# Chain (one observable, O = cHel_P3n):
#   for every event in the boosted-central RECO region, look up alpha_map(m_tt, cosTheta*) from the
#   INCLUSIVE UHH2 strength map using the event's GEN production kinematics, read its GEN observable
#   O = gen_cHel_P3n, form the per-event reweight factor and fill the RECO observable with
#   eventweight * w(f).  f = 0 nominal ; f = +1 "NoSC up" (spin corr removed) ; f = -1 down (enhanced).
# The lookup/observable are pure GEN quantities (inclusive gen collection); the fill is at RECO.
#
# Region  (user): boosted + central |cosTheta*|<0.4, m_tt>800 -> pre-windowed branch
#                 reco_cHel_P3n_Mtt800_Inf_cosThetaLT0p4 (sentinel < -1.5 outside the window).
# Binning (user): 6 bins on [-1, 1], Fig.6 style.
# Bare pyROOT (plain skim + plain TProfile2D map); py2.7/ROOT 6.14 (EL7 container) and py3 both fine.
#
# =============================================================================================
# ** THE PREFACTOR, resolved (Ricardo's "Strength Coefficient Extraction Revisited" slides) **
#
# The opening-angle cosine O has the physical PDF   (1/sigma) dsigma/dO = 1/2 (1 + alpha_a*alpha_b*D*O),
# where alpha_a, alpha_b are the spin-analyzing powers of the two analyzers. For our semileptonic
# lepton + hadronic-b pair,  alpha_a*alpha_b = (1.0)(0.4) = 0.4  (the dilepton pair is (1)(-1), |.|=1,
# which is why the dilepton team never has to separate the two roles of the prefactor -- we do).
#
#   * The MAP stores   alpha_map = 5 * A = D      with A = (N+ - N-)/(N+ + N-) = alpha_a*alpha_b*D/2.
#     The 5 = 2 / (alpha_a*alpha_b) inverts the analyzing-power dilution to recover the TRUE D.
#     Correct for extraction / physics interpretation.
#
#   * The reweighting DENSITY must be the ACTUAL PDF of O, whose slope is the PHYSICAL
#         lambda = alpha_a*alpha_b*D = 0.4 * D = 0.4 * alpha_map = 2 * A.
#     i.e. prefactor 2 on the asymmetry, NOT 5.  Only this makes w|_{f=1} flatten O's distribution
#     (verified on gen: A_up -> 0 at prefactor 2, overshoots to -0.16 at prefactor 5), and it keeps
#     the density positive: |lambda| = |alpha_a alpha_b D| <= 0.4 < 1, so 1 + lambda*O >= 0.6 for any
#     |D|<=1. Prefactor 5 would put the weight-hyperbola pole inside |O|<=1 (negative weights).
#
# Implementation:  lambda = ANALYZER_PRODUCT * alpha_map ;  density d = 1 + lambda*O ;
#                  varied  d_f = 1 + lambda*(1-f)*O ;  w(f) = d_f / d.
# =============================================================================================
from __future__ import print_function
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# ---- config -------------------------------------------------------------------------------
SKIM   = "spincorr_skim_TTToSemiLeptonic_UL18.root"
TREE   = "entSkim"
MAPF   = "strengthmaps_uhh2/NoSC_strengthmap_lb.root"
MAPOBJ = "TTToSemiLeptonic/alpha_lb_cHel_P3n"      # TProfile2D: GetBinContent = alpha_map = D

# the reco fit variable already carries the boosted+central window (m_tt>800 & |cosTheta*|<0.4)
RECO_VAR = "reco_cHel_P3n_Mtt800_Inf_cosThetaLT0p4"
# inclusive-gen inputs for the weight (do NOT use a gen windowed branch: migration is the point)
GEN_O    = "gen_cHel_P3n"
GEN_MTT  = "gen_M_tt"
GEN_CTS  = "gen_cosThetaStar"

NBINS, XLO, XHI = 6, -1.0, 1.0
F_UP, F_DOWN    = +1.0, -1.0

# --- spin-analyzing powers set the two roles of the prefactor (see header) ---
ANALYZER_PRODUCT = 0.4      # alpha_a * alpha_b  for (lepton, hadronic-b);  dilepton would be 1.0
MAP_PREFACTOR    = 5.0      # = 2/ANALYZER_PRODUCT ; the map stores alpha_map = MAP_PREFACTOR*A = D
DENSITY_SIGN     = +1.0     # d = 1 + SIGN*ANALYZER_PRODUCT*alpha_map*O ; the up template MUST flatten
                            # (|up| < |nominal|). If it ever reverses, the map sign is opposite: flip this.

SENTINEL   = -1.5           # branch value <= SENTINEL means "not available / outside window"
D_FLOOR    = 1e-6           # guard on |density|
MIN_MAPENT = 10             # flag map cells with fewer entries than this (noisy alpha)

OUT_ROOT = "templates_cHel_P3n_boostedCentral.root"
OUT_PDF  = "templates_cHel_P3n_boostedCentral.pdf"


def clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)


def main(skim=SKIM, mapf=MAPF, out_root=OUT_ROOT, out_pdf=OUT_PDF):
    # --- load the strength map alpha_map(gen m_tt, cosTheta*) = D ---
    fmap = ROOT.TFile.Open(mapf)
    if (not fmap) or fmap.IsZombie():
        raise RuntimeError("cannot open strength map: %s" % mapf)
    prof = fmap.Get(MAPOBJ)
    if not prof:
        raise RuntimeError("cannot find %s in %s" % (MAPOBJ, mapf))
    prof.SetDirectory(0)
    fmap.Close()
    ax_m, ax_c = prof.GetXaxis(), prof.GetYaxis()
    m_lo, m_hi = ax_m.GetXmin(), ax_m.GetXmax()
    c_lo, c_hi = ax_c.GetXmin(), ax_c.GetXmax()
    eps_m = 1e-3 * (m_hi - m_lo) / prof.GetNbinsX()
    eps_c = 1e-3 * (c_hi - c_lo) / prof.GetNbinsY()

    def alpha_lookup(mtt, cts):
        """alpha_map (= D) at (m_tt, cosTheta*), clamped into range so genuinely-boosted
        (m_tt>2000) / edge events use the nearest valid cell, not an empty under/overflow."""
        b = prof.FindBin(clamp(mtt, m_lo + eps_m, m_hi - eps_m),
                         clamp(cts, c_lo + eps_c, c_hi - eps_c))
        return prof.GetBinContent(b), prof.GetBinEntries(b)

    # --- open the skim, activate only the branches we use ---
    fin = ROOT.TFile.Open(skim)
    if (not fin) or fin.IsZombie():
        raise RuntimeError("cannot open skim: %s" % skim)
    t = fin.Get(TREE)
    t.SetBranchStatus("*", 0)
    for b in ("eventweight", "pass_reco", "gen_semilep", RECO_VAR, GEN_O, GEN_MTT, GEN_CTS):
        t.SetBranchStatus(b, 1)

    region = "pass_reco && %s > %g" % (RECO_VAR, SENTINEL)
    t.Draw(">>elist", region, "goff")
    elist = ROOT.gDirectory.Get("elist")
    npass = elist.GetN()
    print("region '%s' -> %d entries" % (region, npass))

    def book(name, title):
        h = ROOT.TH1F(name, title, NBINS, XLO, XHI); h.Sumw2(); h.SetDirectory(0); return h
    h_nom = book("cHel_P3n",          "nominal;reco cHel_P3n;weighted events")
    h_up  = book("cHel_P3n_NoSCUp",   "NoSC up (f=+1);reco cHel_P3n;weighted events")
    h_dn  = book("cHel_P3n_NoSCDown", "down (f=-1);reco cHel_P3n;weighted events")

    n_gen_used = n_guard = n_neg_up = n_neg_dn = n_lowstat = 0
    d_min = 1e30
    sw = {"nom": [0.0, 0.0], "up": [0.0, 0.0], "dn": [0.0, 0.0]}  # [sum w*sign, sum w]

    def accum(tag, w, o):
        sw[tag][0] += w * (1.0 if o >= 0.0 else -1.0); sw[tag][1] += w

    for i in range(npass):
        t.GetEntry(elist.GetEntry(i))
        w0  = t.eventweight
        o_r = getattr(t, RECO_VAR)

        wu = wd = 1.0
        if t.gen_semilep and getattr(t, GEN_O) > SENTINEL:
            O = getattr(t, GEN_O)
            amap, ent = alpha_lookup(getattr(t, GEN_MTT), getattr(t, GEN_CTS))
            if ent < MIN_MAPENT:
                n_lowstat += 1
            lam = DENSITY_SIGN * ANALYZER_PRODUCT * amap        # physical slope = alpha_a alpha_b D
            d = 1.0 + lam * O
            if d_min > d:
                d_min = d
            if abs(d) < D_FLOOR:
                n_guard += 1
            else:
                n_gen_used += 1
                wu = (1.0 + lam * (1.0 - F_UP)   * O) / d
                wd = (1.0 + lam * (1.0 - F_DOWN) * O) / d
                if wu < 0.0: n_neg_up += 1
                if wd < 0.0: n_neg_dn += 1

        h_nom.Fill(o_r, w0)
        h_up.Fill(o_r, w0 * wu)
        h_dn.Fill(o_r, w0 * wd)
        accum("nom", w0,      o_r)
        accum("up",  w0 * wu, o_r)
        accum("dn",  w0 * wd, o_r)

    fin.Close()

    fout = ROOT.TFile(out_root, "RECREATE")
    for h in (h_nom, h_up, h_dn):
        h.Write()
    fout.Close()

    def coeff(tag):
        s, wsum = sw[tag]; return MAP_PREFACTOR * s / wsum if wsum else 0.0
    print("\n=== yields (weighted, %d-bin [-1,1]) ===" % NBINS)
    print("  nominal %.1f   up %.1f   down %.1f" % (h_nom.Integral(), h_up.Integral(), h_dn.Integral()))
    cn, cu, cd = coeff("nom"), coeff("up"), coeff("dn")
    print("\n=== reco D-tilde estimator  5*<sign(reco cHel_P3n)>  per template ===")
    print("  nominal : %+.4f" % cn)
    print("  up      : %+.4f   (toward 0 -- spin correlation removed)" % cu)
    print("  down    : %+.4f   (away from 0 -- spin correlation enhanced)" % cd)
    print("  --> ordering |up| < |nominal| < |down| : %s" %
          ("OK" if abs(cu) < abs(cn) < abs(cd) else "REVERSED -> flip DENSITY_SIGN"))
    print("\n=== reweighting health (gen-semileptonic subset) ===")
    print("  events reweighted        : %d" % n_gen_used)
    print("  min density 1+lambda*O   : %+.4f   (>= ~0.6 expected; must be > 0)" % d_min)
    print("  |density| < %g guard     : %d" % (D_FLOOR, n_guard))
    print("  negative up / down       : %d / %d" % (n_neg_up, n_neg_dn))
    print("  low-stat map cells (<%d) : %d" % (MIN_MAPENT, n_lowstat))
    print("\nwrote %s  (cHel_P3n, cHel_P3n_NoSCUp, cHel_P3n_NoSCDown)" % out_root)

    make_plot(h_nom, h_up, h_dn, out_pdf, cn, cu, cd)
    print("wrote %s" % out_pdf)


def make_plot(h_nom, h_up, h_dn, out_pdf, cn, cu, cd):
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c_tmpl", "templates", 720, 780); ROOT.SetOwnership(c, False)
    p1 = ROOT.TPad("p1", "", 0, 0.30, 1, 1); p1.SetBottomMargin(0.02); p1.Draw()
    p2 = ROOT.TPad("p2", "", 0, 0.0, 1, 0.30); p2.SetTopMargin(0.02); p2.SetBottomMargin(0.34)
    p2.SetGridy(); p2.Draw()
    ROOT.SetOwnership(p1, False); ROOT.SetOwnership(p2, False)
    p1.cd()
    h_nom.SetLineColor(ROOT.kBlack);     h_nom.SetLineWidth(2); h_nom.SetMarkerStyle(20); h_nom.SetMarkerSize(0.7)
    h_up.SetLineColor(ROOT.kAzure + 1);  h_up.SetLineWidth(2)
    h_dn.SetLineColor(ROOT.kOrange + 7); h_dn.SetLineWidth(2)
    ymax = max(h_nom.GetMaximum(), h_up.GetMaximum(), h_dn.GetMaximum())
    h_nom.SetMaximum(1.45 * ymax); h_nom.SetMinimum(0.0)
    h_nom.GetXaxis().SetLabelSize(0); h_nom.GetYaxis().SetTitle("weighted events")
    h_nom.SetTitle("boosted + central (m_{t#bar{t}}>800, |cos#theta*|<0.4)")
    h_nom.Draw("hist e"); h_up.Draw("hist same"); h_dn.Draw("hist same"); h_nom.Draw("hist e same")
    leg = ROOT.TLegend(0.15, 0.70, 0.55, 0.90); leg.SetBorderSize(0); leg.SetFillStyle(0)
    ROOT.SetOwnership(leg, False)
    leg.AddEntry(h_nom, "nominal (f=0)   #tilde{D}=%+.3f" % cn, "lpe")
    leg.AddEntry(h_up,  "NoSC up (f=+1)  #tilde{D}=%+.3f" % cu, "l")
    leg.AddEntry(h_dn,  "down (f=-1)     #tilde{D}=%+.3f" % cd, "l")
    leg.Draw()
    p2.cd()
    r_up = h_up.Clone("r_up"); r_up.SetDirectory(0); r_up.Divide(h_nom)
    r_dn = h_dn.Clone("r_dn"); r_dn.SetDirectory(0); r_dn.Divide(h_nom)
    r_up.SetTitle(""); r_up.SetMinimum(0.5); r_up.SetMaximum(1.5)
    r_up.GetYaxis().SetTitle("ratio / nom"); r_up.GetYaxis().SetNdivisions(505)
    r_up.GetYaxis().SetTitleSize(0.12); r_up.GetYaxis().SetTitleOffset(0.40); r_up.GetYaxis().SetLabelSize(0.10)
    r_up.GetXaxis().SetTitle("reco cHel_P3n"); r_up.GetXaxis().SetTitleSize(0.13); r_up.GetXaxis().SetLabelSize(0.11)
    r_up.Draw("hist"); r_dn.Draw("hist same")
    line = ROOT.TLine(XLO, 1.0, XHI, 1.0); line.SetLineStyle(2); line.Draw(); ROOT.SetOwnership(line, False)
    c._keep = (p1, p2, r_up, r_dn, leg, line)
    c.Print(out_pdf)


if __name__ == "__main__":
    a = sys.argv[1:]
    kw = {}
    if len(a) >= 1: kw["skim"] = a[0]
    if len(a) >= 2: kw["mapf"] = a[1]
    if len(a) >= 3: kw["out_root"] = a[2]
    if len(a) >= 4: kw["out_pdf"] = a[3]
    main(**kw)
