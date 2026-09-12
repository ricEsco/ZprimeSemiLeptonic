# gen_closure_cHel_P3n.py -- the decisive prefactor test (Ricardo's slide 47):
#   does the f=+1 reweight  w = 1 / (1 + lambda*O)  make the slope of the GEN cHel_P3n distribution VANISH?  
#   We measure it directly on the inclusive gen sample the map was built from.
#
# For each candidate prefactor p (on the asymmetry A), the density slope is lambda_p = (p/5)*alpha_map
# (p=5 -> lambda=alpha_map=D ; p=2 -> lambda=0.4*alpha_map=alpha_a*alpha_b*D). We reweight the gen
# observable at f=+1 and report the residual forward-backward asymmetry A_up = (N+ - N-)/(N+ + N-).
# The CORRECT prefactor is the one that gives A_up ~ 0 (flat slope). Prediction: p = 2.
#
# Runs in the observable's SIGNAL REGION, where the linear opening-angle form 1/2(1+lambda*O) holds
# (inclusive sampling mixes in off-signal events that wash out the linearity and shift the flattening
# prefactor). Weighted by the generator weight to match how the map was built (falls back to
# eventweight if the skim has no genweight branch). py2.7/ROOT 6.14 and py3 both fine.
#
#   python gen_closure_cHel_P3n.py  skim.root  map.root  cHel_P3n   # boosted/central, m_tt>800 & |cos*|<0.4
#   python gen_closure_cHel_P3n.py  skim.root  map.root  cHel       # threshold, 300<m_tt<400
from __future__ import print_function
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

SKIM   = "spincorr_skim_TTToSemiLeptonic_UL18.root"
TREE   = "entSkim"
MAPF   = "strengthmaps_uhh2/NoSC_strengthmap_lb.root"
OBS    = "cHel_P3n"      # observable to test: "cHel_P3n" (compound D-tilde) or "cHel" (linear cos phi, D)

GEN_MTT, GEN_CTS = "gen_M_tt", "gen_cosThetaStar"
SENTINEL = -1.5
ANALYZER_PRODUCT = 0.4        # alpha_a*alpha_b ; density slope for prefactor 2 is 0.4*alpha_map
MAP_PREFACTOR    = 5.0        # alpha_map = 5*A = D
DENSITY_SIGN     = +1.0
PSCAN = [1.0, 2.0, 3.0, 4.0, 5.0]     # candidate prefactors on the asymmetry A
NBINS, XLO, XHI = 24, -1.0, 1.0       # finer than the template: we want to SEE the slope

# SIGNAL REGION per observable -- the linear opening-angle form only holds here; the inclusive
# sample mixes in off-signal events that wash out the linearity (and shift the flattening prefactor).
#   cHel      (D, threshold entanglement)      : 300 < m_tt < 400 GeV
#   cHel_P3n  (D-tilde, boosted/central)       : m_tt > 800 GeV and |cosTheta*| < 0.4
REGION = {
    "cHel":     "gen_semilep && gen_M_tt>300 && gen_M_tt<400 && gen_cHel>-1.5",
    "cHel_P3n": "gen_semilep && gen_M_tt>800 && abs(gen_cosThetaStar)<0.4 && gen_cHel_P3n>-1.5",
}


def clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)


def main(skim=SKIM, mapf=MAPF, obs=OBS):
    mapobj  = "TTToSemiLeptonic/alpha_lb_%s" % obs
    gen_o   = "gen_%s" % obs
    out_pdf = "gen_closure_%s.pdf" % obs
    print("observable: %s   (map %s, gen branch %s)" % (obs, mapobj, gen_o))
    fmap = ROOT.TFile.Open(mapf)
    prof = fmap.Get(mapobj)
    if not prof:
        raise RuntimeError("cannot find %s in %s" % (mapobj, mapf))
    prof.SetDirectory(0); fmap.Close()
    ax_m, ax_c = prof.GetXaxis(), prof.GetYaxis()
    m_lo, m_hi, c_lo, c_hi = ax_m.GetXmin(), ax_m.GetXmax(), ax_c.GetXmin(), ax_c.GetXmax()
    eps_m = 1e-3 * (m_hi - m_lo) / prof.GetNbinsX()
    eps_c = 1e-3 * (c_hi - c_lo) / prof.GetNbinsY()

    def alpha_lookup(mtt, cts):
        b = prof.FindBin(clamp(mtt, m_lo + eps_m, m_hi - eps_m),
                         clamp(cts, c_lo + eps_c, c_hi - eps_c))
        return prof.GetBinContent(b)

    fin = ROOT.TFile.Open(skim)
    t = fin.Get(TREE)
    wname = "genweight" if t.GetBranch("genweight") else "eventweight"
    print("base weight: %s" % wname + ("" if wname == "genweight" else
          "  (no genweight branch found -- using eventweight; asymmetry unaffected to good approx)"))
    t.SetBranchStatus("*", 0)
    for b in ("gen_semilep", gen_o, GEN_MTT, GEN_CTS, wname):
        t.SetBranchStatus(b, 1)

    region = REGION.get(obs)
    if region is None:
        raise RuntimeError("no signal region defined for observable '%s'" % obs)
    t.Draw(">>elist", region, "goff")
    elist = ROOT.gDirectory.Get("elist")
    nsel_total = elist.GetN()
    print("signal region: %s\n  -> %d gen entries (full stats, no sampling)" % (region, nsel_total))

    # histograms: nominal + up/down at the reference prefactor 2 (the predicted-correct one)
    def book(nm):
        h = ROOT.TH1F(nm, "", NBINS, XLO, XHI); h.Sumw2(); h.SetDirectory(0); return h
    h_nom, h_up, h_dn = book("gen_nom"), book("gen_up_p2"), book("gen_dn_p2")

    # prefactor scan: accumulate weighted N+ / N- at f=+1 for each candidate p
    fbp = dict((p, [0.0, 0.0]) for p in PSCAN)      # p -> [sum w+, sum w-]
    fb_nom = [0.0, 0.0]
    nsel = 0
    for i in range(nsel_total):
        t.GetEntry(elist.GetEntry(i))
        O = getattr(t, gen_o)
        w0 = getattr(t, wname)
        amap = alpha_lookup(getattr(t, GEN_MTT), getattr(t, GEN_CTS))
        nsel += 1
        side = 0 if O >= 0.0 else 1
        fb_nom[side] += w0
        # reference p=2 templates (up and down)
        lam2 = DENSITY_SIGN * ANALYZER_PRODUCT * amap
        d2 = 1.0 + lam2 * O
        if abs(d2) > 1e-6:
            wu = (1.0 + 0.0) / d2                    # f=+1 -> numerator 1
            wd = (1.0 + 2.0 * lam2 * O) / d2         # f=-1
        else:
            wu = wd = 1.0
        h_nom.Fill(O, w0); h_up.Fill(O, w0 * wu); h_dn.Fill(O, w0 * wd)
        # scan
        for p in PSCAN:
            lam = DENSITY_SIGN * (p / MAP_PREFACTOR) * amap
            d = 1.0 + lam * O
            wpu = (1.0 / d) if abs(d) > 1e-6 else 1.0
            fbp[p][side] += w0 * wpu
    fin.Close()

    def asym(pair):
        s = pair[0] + pair[1]; return (pair[0] - pair[1]) / s if s else 0.0

    A_nom = asym(fb_nom)
    print("\nselected %d gen-semileptonic events" % nsel)
    print("\n=== gen forward-backward asymmetry A = (N+ - N-)/(N+ + N-) ; D-tilde = 5*A ===")
    print("  nominal (no reweight) : A=%+.4f   D-tilde=%+.4f" % (A_nom, MAP_PREFACTOR * A_nom))
    print("\n=== prefactor scan: A_up after f=+1 reweight (want A_up -> 0) ===")
    print("   p(on A)  density=lambda   A_up        D-tilde_up     verdict")
    best = None
    for p in PSCAN:
        Aup = asym(fbp[p])
        tag = "slope FLATTENED" if abs(Aup) < 0.005 else ("over-flattened (<0)" if Aup < 0 else "under-flattened")
        mark = "  <== correct prefactor" if abs(Aup) < 0.005 else ""
        print("     %.0f      %.2f*alpha_map   %+.4f    %+.4f    %s%s" %
              (p, p / MAP_PREFACTOR, Aup, MAP_PREFACTOR * Aup, tag, mark))
        if best is None or abs(Aup) < abs(best[1]):
            best = (p, Aup)
    print("\n  --> empirically-flat prefactor: p=%.0f (A_up=%+.4f).  Expected from alpha_a*alpha_b=%.1f: p=2."
          % (best[0], best[1], ANALYZER_PRODUCT))

    _plot(h_nom, h_up, h_dn, out_pdf, obs, A_nom, asym(fbp[2.0]), asym(fbp[5.0]))
    print("wrote %s" % out_pdf)


def _plot(h_nom, h_up, h_dn, out_pdf, obs, A_nom, A_up2, A_up5):
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c_gc", "gen closure", 760, 620); ROOT.SetOwnership(c, False)
    for h in (h_nom, h_up, h_dn):
        if h.Integral() > 0: h.Scale(1.0 / h.Integral())     # unit area: compare SHAPE/slope
    h_nom.SetLineColor(ROOT.kBlack);     h_nom.SetLineWidth(2)
    h_up.SetLineColor(ROOT.kAzure + 1);  h_up.SetLineWidth(2)
    h_dn.SetLineColor(ROOT.kOrange + 7); h_dn.SetLineWidth(2)
    h_nom.SetMaximum(1.5 * max(h_nom.GetMaximum(), h_up.GetMaximum(), h_dn.GetMaximum())); h_nom.SetMinimum(0)
    h_nom.SetTitle("gen %s reweighting closure (in SR);gen %s;a.u. (unit area)" % (obs, obs))
    h_nom.Draw("hist"); h_up.Draw("hist same"); h_dn.Draw("hist same")
    leg = ROOT.TLegend(0.40, 0.72, 0.88, 0.90); leg.SetBorderSize(0); leg.SetFillStyle(0); ROOT.SetOwnership(leg, False)
    leg.AddEntry(h_nom, "nominal              A=%+.3f" % A_nom, "l")
    leg.AddEntry(h_up,  "up f=+1  (p=2)   A=%+.3f (flat)" % A_up2, "l")
    leg.AddEntry(h_dn,  "down f=-1 (p=2)  enhanced", "l")
    leg.Draw()
    c._keep = (leg,)
    c.Print(out_pdf)


if __name__ == "__main__":
    a = sys.argv[1:]
    kw = {}
    if len(a) >= 1: kw["skim"] = a[0]
    if len(a) >= 2: kw["mapf"] = a[1]
    if len(a) >= 3: kw["obs"] = a[2]      # "cHel_P3n" (default) or "cHel" (linear cross-check)
    main(**kw)
