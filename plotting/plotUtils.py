"""
Shared helpers of the plotting scripts:
  - summing histograms over the per-job SFrame outputs (one file open at a time) with a small ROOT cache file
  - CMS styling (tdrstyle + CMS_lumi), stacked MC with data/prediction ratio, nominal-vs-variation plots
"""
import array
import ctypes
import glob
import math
import os

import ROOT

import CMS_lumi
import plotConfig as cfg
import tdrstyle

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


# =================================================================================================================== #
#                                                  Reading the input                                                  #
# =================================================================================================================== #
def job_files(directory, kind, version, year):
    """per-job SFrame outputs of one sample: uhh2.AnalysisModuleRunner.<kind>.<version>_<year>_<n>.root"""
    return sorted(glob.glob(os.path.join(directory, "uhh2.AnalysisModuleRunner.%s.%s_%s_[0-9]*.root" % (kind, version, year))))


def collect(directory, year, groups, kind, paths):
    """
    Sum the histograms 'paths' (folder/name) over all job files of the samples in each group.
    groups: {group: [sample versions]}, kind: MC or DATA.  returns {group: {path: TH1}}
    Only one file is open at a time and the histograms are detached from it, so memory stays at the size of the sums.
    """
    out = {g: {} for g in groups}
    todo = [(g, v, f) for g, versions in groups.items() for v in versions for f in job_files(directory, kind, v, year)]
    for g, versions in groups.items():
        for v in versions:
            if not job_files(directory, kind, v, year):
                print("  WARNING: no %s job outputs for %s_%s in %s (skipped)" % (kind, v, year, directory))
    missing, bad = {}, set()
    for i, (g, v, fname) in enumerate(todo):
        if i % 10 == 0 or i == len(todo) - 1:
            print("  reading %s files: %d / %d" % (kind, i + 1, len(todo)), flush=True)
        f = ROOT.TFile.Open(fname)
        if not f or f.IsZombie():
            raise RuntimeError("cannot open %s" % fname)
        for p in paths:
            if p in bad:
                continue
            h = f.Get(p)
            if not h:
                missing.setdefault(p, set()).add(v)
                continue
            ref = next((o[p] for o in out.values() if p in o), None)
            if h.GetXaxis().CanExtend() or (ref is not None and not same_axis(ref, h)):
                bad.add(p)   # automatic axis range (booked with xlow == xup): the jobs have different ranges
                continue
            if p in out[g]:
                out[g][p].Add(h)
            else:
                out[g][p] = h.Clone("%s__%s" % (g, p.replace("/", "__")))
                out[g][p].SetDirectory(0)
        f.Close()
    for p, versions in sorted(missing.items())[:5]:
        print("  note: %s missing in %d sample(s), e.g. %s" % (p, len(versions), sorted(versions)[0]))
    if len(missing) > 5:
        print("  note: ... %d more histograms missing in some samples (e.g. NoSC folders exist for the signal only)" % (len(missing) - 5))
    for o in out.values():
        for p in bad:
            o.pop(p, None)
    return out, bad


def same_axis(a, b):
    ax, bx = a.GetXaxis(), b.GetXaxis()
    return (ax.GetNbins() == bx.GetNbins() and abs(ax.GetXmin() - bx.GetXmin()) < 1e-9
            and abs(ax.GetXmax() - bx.GetXmax()) < 1e-9)


def warn_skipped(bad):
    for p in sorted(bad):
        print("  WARNING: %s skipped: its axis range differs between job files (booked with xlow == xup, i.e. an "
              "automatic range); fix the booking in the C++ histogram class" % p)


def list_folder(directory, year, version, folder):
    """histogram names in one folder of the first job output of a sample ([] if the folder does not exist)"""
    files = job_files(directory, "MC", version, year)
    if not files:
        raise RuntimeError("no job outputs for %s_%s in %s" % (version, year, directory))
    f = ROOT.TFile.Open(files[0])
    d = f.GetDirectory(folder)
    names = [k.GetName() for k in d.GetListOfKeys()] if d else []
    f.Close()
    return names


def cached(cache_file, sources, rebuild):
    """
    sources: list of (directory, year, groups, kind, paths). Returns {group: {path: TH1}} for all groups, read from
    cache_file when it already holds every requested histogram, otherwise from the job outputs (then cached).
    """
    want = sorted(set(p for s in sources for p in s[4]))
    groups_kinds = sources
    if os.path.exists(cache_file) and not rebuild:
        f = ROOT.TFile.Open(cache_file)
        req = f.Get("requested")
        have = set(req.GetTitle().split(";")) if req else set()
        skipped = f.Get("skipped")
        bad = set(x for x in skipped.GetTitle().split(";") if x) if skipped else set()
        groups = [g for _, _, gr, _, _ in groups_kinds for g in gr]
        if set(want) <= have and all(f.GetDirectory(g) for g in groups):
            out = {}
            for g in groups:
                d = f.GetDirectory(g)
                out[g] = {}
                for key in d.GetListOfKeys():
                    h = key.ReadObj()
                    h.SetDirectory(0)
                    out[g][key.GetName().replace("|", "/")] = h
            f.Close()
            print("  using cached sums from %s (--rebuild to re-read the job outputs)" % cache_file)
            warn_skipped(bad & set(want))
            return out
        f.Close()
    out, bad = {}, set()
    for directory, year, groups, kind, paths in groups_kinds:
        o, b = collect(directory, year, groups, kind, sorted(set(paths)))
        out.update(o)
        bad |= b
    for o in out.values():
        for p in bad:
            o.pop(p, None)
    warn_skipped(bad)
    os.makedirs(os.path.dirname(os.path.abspath(cache_file)), exist_ok=True)
    f = ROOT.TFile(cache_file, "RECREATE")
    ROOT.TNamed("requested", ";".join(want)).Write()
    ROOT.TNamed("skipped", ";".join(sorted(bad))).Write()
    for g, hists in out.items():
        d = f.mkdir(g)
        d.cd()
        for p, h in hists.items():
            h.Write(p.replace("/", "|"))
    f.Close()
    print("  cached sums written to %s" % cache_file)
    return out


def total(hists, groups, path, name):
    """sum of path over the given groups (None if no group has it)"""
    h = None
    for g in groups:
        x = hists.get(g, {}).get(path)
        if x is None:
            continue
        if h is None:
            h = x.Clone(name)
            h.SetDirectory(0)
        else:
            h.Add(x)
    return h


# =================================================================================================================== #
#                                                Systematic variations                                                #
# =================================================================================================================== #
def varied_total(hists, cells, region, obs, nuis, suffix, name):
    """
    MC sum with one nuisance shifted: cells of the processes the nuisance acts on take <region>_SC_<suffix>, the
    others (and any cell without that folder) keep the nominal.
    """
    nom_path = "%s_SC_nominal/%s" % (region, obs)
    var_path = "%s_SC_%s/%s" % (region, suffix, obs)
    h = None
    for c in cells:
        acts = nuis["procs"] is None or cfg.cell_process(c) in nuis["procs"]
        x = hists[c].get(var_path) if acts else None
        if x is None:
            x = hists[c].get(nom_path)
        if x is None:
            continue
        if h is None:
            h = x.Clone(name)
            h.SetDirectory(0)
        else:
            h.Add(x)
    return h


def pdf_hessian(nominal, replicas, name):
    """NNPDF Hessian: sigma = sqrt(sum_i (h_i - h_0)^2) per bin -> (nominal + sigma, nominal - sigma)"""
    up, down = nominal.Clone(name + "_up"), nominal.Clone(name + "_down")
    for b in range(0, nominal.GetNbinsX() + 2):
        s = math.sqrt(sum((r.GetBinContent(b) - nominal.GetBinContent(b)) ** 2 for r in replicas))
        up.SetBinContent(b, nominal.GetBinContent(b) + s)
        down.SetBinContent(b, nominal.GetBinContent(b) - s)
    return up, down


def variation_pair(hists, cells, region, obs, nuis, prep):
    """(up, down) MC sums of one nuisance after prep (rebinning); None if the nuisance has no folders here"""
    if nuis["up"] == "pdf":
        nom = varied_total(hists, cells, region, obs, {"procs": []}, "nominal", "pdfnom")
        reps = [varied_total(hists, cells, region, obs, nuis, "pdf%d" % i, "pdf%d" % i) for i in range(cfg.N_PDF)]
        if nom is None or any(r is None for r in reps):
            return None
        nom, reps = prep(nom), [prep(r) for r in reps]
        return pdf_hessian(nom, reps, "pdf")
    if not any(("%s_SC_%s/%s" % (region, nuis["up"], obs)) in hists[c] for c in cells):
        return None
    up = varied_total(hists, cells, region, obs, nuis, nuis["up"], nuis["name"] + "_up")
    down = varied_total(hists, cells, region, obs, nuis, nuis["down"], nuis["name"] + "_down")
    return prep(up), prep(down)


def sc_paths(region, obs, nuis_list):
    """all histograms needed for the nominal and the variations of one SC observable"""
    suffixes = ["nominal"]
    for n in nuis_list:
        if n["up"] == "pdf":
            suffixes += ["pdf%d" % i for i in range(cfg.N_PDF)]
        else:
            suffixes += [n["up"], n["down"]]
    return ["%s_SC_%s/%s" % (region, s, obs) for s in dict.fromkeys(suffixes)]


# =================================================================================================================== #
#                                                Histogram preparation                                                #
# =================================================================================================================== #
class BinningError(Exception):
    pass


def prepare(h, binning=None, fold=True, rebin=1):
    """
    Copy of h with a new binning and the under/overflow treatment:
      binning (nbins, lo, hi): every new edge has to coincide with a booked edge (BinningError otherwise)
      binning None           : booked binning, merged by 'rebin' when it divides the number of bins
      fold                   : add everything below / above the plotted range to the first / last bin
    """
    ax = h.GetXaxis()
    booked = [ax.GetBinLowEdge(b) for b in range(1, ax.GetNbins() + 2)]
    if binning is not None:
        n, lo, hi = binning
        edges = [lo + k * (hi - lo) / float(n) for k in range(n + 1)]
        tol = 1e-6 * (booked[-1] - booked[0])
        idx = []
        for e in edges:
            match = [i for i, b in enumerate(booked) if abs(b - e) < tol]
            if not match:
                raise BinningError("binning (%d, %g, %g) does not match the booked binning (%d, %g, %g): edge %g is not a booked edge"
                                   % (n, lo, hi, ax.GetNbins(), booked[0], booked[-1], e))
            idx.append(match[0] + 1)   # booked bin that starts at this edge
    else:
        nb = ax.GetNbins()
        r = rebin if rebin > 1 and nb % rebin == 0 else 1
        idx = list(range(1, nb + 2, r))
        edges = [booked[i - 1] for i in idx]
    out = ROOT.TH1D(h.GetName() + "_prep", h.GetTitle(), len(edges) - 1, array.array("d", edges))
    out.SetDirectory(0)
    for k in range(len(idx) - 1):
        bins = range(idx[k], idx[k + 1])
        out.SetBinContent(k + 1, sum(h.GetBinContent(b) for b in bins))
        out.SetBinError(k + 1, math.sqrt(sum(h.GetBinError(b) ** 2 for b in bins)))
    below = range(0, idx[0])
    above = range(idx[-1], ax.GetNbins() + 2)
    if fold:
        for b, extra in ((1, below), (out.GetNbinsX(), above)):
            out.SetBinContent(b, out.GetBinContent(b) + sum(h.GetBinContent(x) for x in extra))
            out.SetBinError(b, math.sqrt(out.GetBinError(b) ** 2 + sum(h.GetBinError(x) ** 2 for x in extra)))
    else:
        out.SetBinContent(0, sum(h.GetBinContent(x) for x in below))
        out.SetBinContent(out.GetNbinsX() + 1, sum(h.GetBinContent(x) for x in above))
    return out


def y_title(h, name):
    """'Events / 10 GeV' for uniform binning, 'Events' for counts and variable binning"""
    widths = {round(h.GetBinWidth(b), 6) for b in range(1, h.GetNbinsX() + 1)}
    if name.startswith("N_") or len(widths) > 1:
        return "Events"
    w = widths.pop()
    unit = " GeV" if "[GeV]" in h.GetTitle() else ""
    return "Events / %s%s" % (("%.3g" % w), unit)


def x_title(h):
    return h.GetTitle() or h.GetName()


# =================================================================================================================== #
#                                                       Styling                                                       #
# =================================================================================================================== #
LABEL_PX, TITLE_PX, TEXT_PX = 26, 30, 26
COL_DATA = ROOT.kBlack
COL_STAT = ROOT.TColor.GetColor("#c8c8cc")
COL_UP, COL_DOWN = ROOT.TColor.GetColor("#5790fc"), ROOT.TColor.GetColor("#e42536")


def setup_style():
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetHatchesLineWidth(1)
    ROOT.gStyle.SetHatchesSpacing(0.8)
    ROOT.gStyle.SetEndErrorSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.TGaxis.SetMaxDigits(4)
    # ROOT.TGaxis.SetExponentOffset(-0.065, -0.01, "y")   # moves the y-axis "x10^n" left, next to the axis labels
    CMS_lumi.writeExtraText = True
    CMS_lumi.relPosX = 0.045


def style_axes(h, pad_h_px, xtitle, ytitle, show_x=True, ytitle_offset=1.6, xtitle_offset=1.05):
    for ax, title in ((h.GetXaxis(), xtitle), (h.GetYaxis(), ytitle)):
        ax.SetTitle(title)
        ax.SetLabelFont(43)
        ax.SetTitleFont(43)
        ax.SetLabelSize(LABEL_PX)
        ax.SetTitleSize(TITLE_PX)
    h.GetYaxis().SetTitleOffset(ytitle_offset)
    h.GetXaxis().SetTitleOffset(xtitle_offset)
    h.GetXaxis().SetTickLength(0.03 * 560.0 / pad_h_px)
    if not show_x:
        h.GetXaxis().SetLabelSize(0)
        h.GetXaxis().SetTitleSize(0)


def cms_label(pad, year, with_data, lumi=True):
    """CMS logo + extra text inside the frame (top left); luminosity label above the frame if lumi"""
    CMS_lumi.extraText = cfg.EXTRA_TEXT_DATA if with_data else cfg.EXTRA_TEXT_SIM
    period = cfg.LUMI_PERIOD[year] if lumi else ""
    CMS_lumi.lumi_sqrtS = ""
    CMS_lumi.CMS_lumi(pad, period, 11)


def cms_label_outside(pad, year, with_data):
    """'CMS' + extra text above the frame on the left, luminosity above the frame on the right (CMS_lumi iPos 0)"""
    CMS_lumi.extraText = cfg.EXTRA_TEXT_DATA if with_data else cfg.EXTRA_TEXT_SIM
    CMS_lumi.lumi_sqrtS = ""
    pad.cd()
    # CMS_lumi draws the extra text at relPosX; put it right after the width of "CMS"
    lat = ROOT.TLatex(0.5, 0.5, CMS_lumi.cmsText)
    lat.SetNDC()
    lat.SetTextFont(CMS_lumi.cmsTextFont)
    lat.SetTextSize(CMS_lumi.cmsTextSize * pad.GetTopMargin())
    w, h = ctypes.c_uint(0), ctypes.c_uint(0)
    lat.GetBoundingBox(w, h)
    width_ndc = w.value / float(pad.GetWw() * pad.GetWNDC())
    l, r = pad.GetLeftMargin(), pad.GetRightMargin()
    shift = 0.085                       # NDC: how far "CMS" moves right of the frame edge; tune to taste
    saved = CMS_lumi.relPosX
    CMS_lumi.cmsTextShift = shift
    CMS_lumi.relPosX = (shift + width_ndc + 0.012) / (1 - l - r)   # extra text follows "CMS"
    CMS_lumi.CMS_lumi(pad, cfg.LUMI_PERIOD[year], 0)
    CMS_lumi.relPosX = saved


def text_lines(pad, lines, x, y_top, align=13):
    """plain text lines (channel, region, ...) at NDC x, starting at y_top"""
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(43)
    lat.SetTextSize(TEXT_PX)
    lat.SetTextAlign(align)
    pad.cd()
    step = 1.35 * TEXT_PX / (pad.GetWh() * pad.GetHNDC())
    for i, line in enumerate(lines):
        lat.DrawLatex(x, y_top - i * step, line)
    return lat


def text_start_below_logo(pad, with_extra=True):
    """NDC y below the CMS logo and extra text drawn by CMS_lumi (iPos 11)"""
    t, b = pad.GetTopMargin(), pad.GetBottomMargin()
    y = 1 - t - CMS_lumi.relPosY * (1 - t - b)
    y -= CMS_lumi.relExtraDY * CMS_lumi.cmsTextSize * t
    if with_extra:
        y -= CMS_lumi.relExtraDY * CMS_lumi.extraOverCmsTextSize * CMS_lumi.cmsTextSize * t
    return y - 0.02


def make_pads(canvas, x0, x1, split=0.32, left=0.16, right=0.05, name=""):
    """upper (distribution) and lower (ratio) pad in the canvas column [x0, x1]"""
    canvas.cd()
    up = ROOT.TPad("up" + name, "", x0, split, x1, 1.0)
    lo = ROOT.TPad("lo" + name, "", x0, 0.0, x1, split)
    for p in (up, lo):
        p.SetLeftMargin(left)
        p.SetRightMargin(right)
        p.SetTicks(1, 1)
        p.SetFillStyle(4000)
        p.SetFrameFillStyle(4000)
    up.SetTopMargin(0.075)
    up.SetBottomMargin(0.025)
    lo.SetTopMargin(0.03)
    lo.SetBottomMargin(0.36)
    up.Draw()
    lo.Draw()
    return up, lo


def pad_px(pad):
    return pad.GetWh() * pad.GetHNDC()


# =================================================================================================================== #
#                                       Graphs: data points and uncertainty bands                                     #
# =================================================================================================================== #
def poisson_graph(h, scale_by=None):
    """data as points with Garwood (Poisson) intervals; optionally divided by scale_by (ratio). Empty bins skipped."""
    g = ROOT.TGraphAsymmErrors()
    alpha = 1 - 0.6827
    for b in range(1, h.GetNbinsX() + 1):
        n = h.GetBinContent(b)
        if n <= 0:
            continue
        lo = n - (ROOT.Math.gamma_quantile(alpha / 2, n, 1.0) if n > 0 else 0.)
        hi = ROOT.Math.gamma_quantile_c(alpha / 2, n + 1, 1.0) - n
        den = 1.0
        if scale_by is not None:
            den = scale_by.GetBinContent(b)
            if den <= 0:
                continue
        i = g.GetN()
        g.SetPoint(i, h.GetBinCenter(b), n / den)
        g.SetPointError(i, 0, 0, lo / den, hi / den)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.0)
    g.SetLineWidth(2)
    g.SetMarkerColor(COL_DATA)
    g.SetLineColor(COL_DATA)
    return g


def band_graph(h_nom, err_up, err_down, relative=False):
    """box band around the nominal (or around 1 when relative) from per-bin absolute uncertainties"""
    g = ROOT.TGraphAsymmErrors()
    for b in range(1, h_nom.GetNbinsX() + 1):
        y = h_nom.GetBinContent(b)
        w = 0.5 * h_nom.GetBinWidth(b)
        if relative:
            if y <= 0:
                continue
            g.SetPoint(g.GetN(), h_nom.GetBinCenter(b), 1.0)
            g.SetPointError(g.GetN() - 1, w, w, err_down[b - 1] / y, err_up[b - 1] / y)
        else:
            g.SetPoint(g.GetN(), h_nom.GetBinCenter(b), y)
            g.SetPointError(g.GetN() - 1, w, w, err_down[b - 1], err_up[b - 1])
    return g


def stat_errors(h):
    return [h.GetBinError(b) for b in range(1, h.GetNbinsX() + 1)]


# =================================================================================================================== #
#                                       Stacked MC (+ data) with ratio panel                                          #
# =================================================================================================================== #
def draw_stack_plot(canvas, stack, data, xtitle, ytitle, logy, text, year, ncol=3):
    """
    Stacked MC with the MC stat. uncertainty, data if given, and data / prediction (or the relative MC stat.
    uncertainty without data) in the lower panel. CMS label above the frame, legend across the top of the frame,
    'text' (one line, e.g. channel and selection step) below the legend.
    stack: list of (label, colour hex, TH1) bottom -> top.  Returns the objects that must stay alive.
    """
    keep = [h for _, _, h in stack] + ([data] if data is not None else [])   # referenced by the legend
    up, lo = make_pads(canvas, 0, 1)
    keep += [up, lo]
    tot = stack[0][2].Clone("tot_" + stack[0][2].GetName())
    tot.Reset()
    for _, _, h in stack:
        tot.Add(h)
    keep.append(tot)
    stat = stat_errors(tot)

    # ---- upper pad
    up.cd()
    up.SetLogy(logy)
    # stack drawn as cumulative sums, top (total) first (THStack over a frame crashes in some ROOT versions)
    cumulative, running = [], None
    for label, colour, h in stack:
        c = ROOT.TColor.GetColor(colour)
        h.SetFillColor(c)
        h.SetLineColor(c)
        h.SetLineWidth(0)
        running = h.Clone(h.GetName() + "_cum") if running is None else running.Clone(h.GetName() + "_cum")
        if len(cumulative):
            running.Add(h)
        running.SetDirectory(0)
        running.SetFillColor(c)
        running.SetLineColor(c)
        running.SetLineWidth(0)
        cumulative.append(running)
    keep += cumulative

    n_leg = len(stack) + 1 + (1 if data is not None else 0)
    rows = int(math.ceil(n_leg / float(ncol)))
    t, b = up.GetTopMargin(), up.GetBottomMargin()
    line_h = 1.3 * TEXT_PX / pad_px(up)    # was 1.4
    leg_top = 1 - t - 0.025
    leg_bottom = leg_top - rows * line_h
    text_y = leg_bottom - 0.015
    free = text_y - 1.2 * line_h * len(text) - 0.02          # NDC below which the histograms have to stay
    frac = max((free - b) / (1 - t - b), 0.3)

    frame = tot.Clone("frame_" + tot.GetName())
    frame.Reset()
    keep.append(frame)
    top = max(tot.GetBinContent(i) + stat[i - 1] for i in range(1, tot.GetNbinsX() + 1))
    if data is not None:
        top = max(top, data.GetMaximum() + math.sqrt(max(data.GetMaximum(), 0)))
    if logy:
        positive = [tot.GetBinContent(i) for i in range(1, tot.GetNbinsX() + 1) if tot.GetBinContent(i) > 0]
        ymin = max(min(positive) * 0.3, 0.1) if positive else 0.1
        frame.SetMinimum(ymin)
        frame.SetMaximum(10 ** (math.log10(ymin) + (math.log10(max(top, ymin * 10)) - math.log10(ymin)) / frac))
    else:
        frame.SetMinimum(0)
        frame.SetMaximum(max(top, 1e-9) / frac)
    style_axes(frame, pad_px(up), xtitle, ytitle, show_x=False)
    frame.Draw("axis")
    for h in reversed(cumulative):
        h.Draw("hist same")
    band = band_graph(tot, stat, stat)
    band.SetFillStyle(3354)
    band.SetFillColor(ROOT.kGray + 2)
    band.SetLineWidth(0)
    band.Draw("2 same")
    keep.append(band)
    gd = None
    if data is not None:
        gd = poisson_graph(data)
        gd.Draw("P0Z same")
        keep.append(gd)
    up.RedrawAxis()

    x0 = up.GetLeftMargin() + 0.025
    x1 = 1 - up.GetRightMargin() - 0.015
    leg = ROOT.TLegend(x0, leg_bottom, x1, leg_top)
    leg.SetNColumns(ncol)
    leg.SetTextFont(43)
    leg.SetTextSize(TEXT_PX - 8)    # was TEXT_PX - 2; pixel size (font 43), TEXT_PX = 26
    leg.SetFillStyle(0)
    leg.SetColumnSeparation(0.02)
    entries = ([(gd, "Data", "pe")] if gd is not None else []) + [(h, label, "f") for label, _, h in reversed(stack)]
    entries.append((band, "MC stat. unc.", "f"))
    for obj, label, opt in entries:
        leg.AddEntry(obj, label, opt)
    leg.Draw()
    keep.append(leg)
    keep.append(text_lines(up, text, x0 + 0.005, text_y))
    cms_label_outside(up, year, data is not None)

    # ---- ratio pad
    lo.cd()
    rframe = tot.Clone("rframe_" + tot.GetName())
    rframe.Reset()
    keep.append(rframe)
    # if data is not None:
    #     rmin, rmax = 0.5, 1.5
    # else:
    #     rel = [stat[i - 1] / tot.GetBinContent(i) for i in range(1, tot.GetNbinsX() + 1) if tot.GetBinContent(i) > 0]
    #     r = min(max(1.5 * max(rel or [0.05]), 0.02), 0.95)
    #     rmin, rmax = 1 - r, 1 + r
    rmin, rmax = 0.8, 1.2
    rframe.SetMinimum(rmin)
    rframe.SetMaximum(rmax)
    style_axes(rframe, pad_px(lo), xtitle, "Data / MC" if data is not None else "Rel. unc.")
    rframe.GetYaxis().SetNdivisions(505)
    rframe.Draw("axis")
    b_stat = band_graph(tot, stat, stat, relative=True)
    b_stat.SetFillColor(COL_STAT)
    b_stat.SetFillStyle(1001)
    b_stat.SetLineWidth(0)
    b_stat.Draw("2 same")
    keep.append(b_stat)
    line = ROOT.TLine(tot.GetXaxis().GetXmin(), 1, tot.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(2)
    line.Draw()
    keep.append(line)
    if data is not None:
        gr = poisson_graph(data, scale_by=tot)
        gr.Draw("P0Z same")
        keep.append(gr)
    lo.RedrawAxis()
    return keep


# =================================================================================================================== #
#                                           Nominal vs up/down variation                                              #
# =================================================================================================================== #
def draw_variation(canvas, nominal, up_h, down_h, name, xtitle, year, text):
    """nominal (black) vs up / down variation, with variation / nominal in the lower panel"""
    keep = [nominal, up_h, down_h]   # referenced by the legend
    upad, lpad = make_pads(canvas, 0, 1, name="var")
    upad.cd()
    frame = nominal.Clone("vframe")
    frame.Reset()
    keep.append(frame)
    top = max(h.GetMaximum() for h in (nominal, up_h, down_h))
    t, b = upad.GetTopMargin(), upad.GetBottomMargin()
    text_bottom = 1 - t - 0.03 - 1.45 * TEXT_PX / pad_px(upad) * (len(text) + 3)
    frame.SetMinimum(0)
    frame.SetMaximum(top / max((text_bottom - b) / (1 - t - b) - 0.03, 0.3))
    style_axes(frame, pad_px(upad), xtitle, y_title(nominal, ""), show_x=False)
    frame.Draw("axis")
    for h, col, width in ((nominal, ROOT.kBlack, 3), (up_h, COL_UP, 2), (down_h, COL_DOWN, 2)):
        h.SetLineColor(col)
        h.SetMarkerColor(col)
        h.SetMarkerSize(0)
        h.SetLineWidth(width)
        h.SetFillStyle(0)
    up_h.Draw("hist same")
    down_h.Draw("hist same")
    nominal.Draw("hist same")
    nom_err = nominal.Clone("nom_err")
    nom_err.Draw("e1 same")
    keep.append(nom_err)
    upad.RedrawAxis()

    def change(h):
        n = nominal.Integral()
        return 100. * (h.Integral() - n) / n if n else 0.

    lx1 = 1 - upad.GetRightMargin() - 0.03
    ly1 = 1 - upad.GetTopMargin() - 0.03
    step = 1.45 * TEXT_PX / pad_px(upad)
    leg = ROOT.TLegend(lx1 - 0.34, ly1 - 3 * step, lx1, ly1)
    leg.SetTextFont(43)
    leg.SetTextSize(TEXT_PX)
    leg.SetFillStyle(0)
    leg.AddEntry(nominal, "Nominal", "l")
    leg.AddEntry(up_h, "Up (%+.2f%%)" % change(up_h), "l")
    leg.AddEntry(down_h, "Down (%+.2f%%)" % change(down_h), "l")
    leg.Draw()
    keep.append(leg)
    cms_label(upad, year, with_data=False)
    keep.append(text_lines(upad, text + ["Nuisance: %s" % name],
                           upad.GetLeftMargin() + CMS_lumi.relPosX * (1 - upad.GetLeftMargin() - upad.GetRightMargin()),
                           text_start_below_logo(upad)))

    lpad.cd()
    rframe = nominal.Clone("vrframe")
    rframe.Reset()
    keep.append(rframe)
    ratios = []
    for h in (up_h, down_h):
        r = h.Clone(h.GetName() + "_ratio")
        r.SetDirectory(0)
        for b in range(1, r.GetNbinsX() + 2):
            n = nominal.GetBinContent(b)
            r.SetBinContent(b, h.GetBinContent(b) / n if n else 1.)
            r.SetBinError(b, 0.)
        ratios.append(r)
    dev = max([abs(r.GetBinContent(b) - 1) for r in ratios for b in range(1, r.GetNbinsX() + 1)] +
              [nominal.GetBinError(b) / nominal.GetBinContent(b) for b in range(1, nominal.GetNbinsX() + 1) if nominal.GetBinContent(b) > 0])
    d = min(max(1.3 * dev, 0.005), 0.9)
    rframe.SetMinimum(1 - d)
    rframe.SetMaximum(1 + d)
    style_axes(rframe, pad_px(lpad), xtitle, "Var. / nom.")
    rframe.GetYaxis().SetNdivisions(505)
    rframe.Draw("axis")
    stat = band_graph(nominal, stat_errors(nominal), stat_errors(nominal), relative=True)
    stat.SetFillColor(COL_STAT)
    stat.SetFillStyle(1001)
    stat.SetLineWidth(0)
    stat.Draw("2 same")
    keep.append(stat)
    line = ROOT.TLine(nominal.GetXaxis().GetXmin(), 1, nominal.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(2)
    line.Draw()
    keep.append(line)
    for r in ratios:
        r.Draw("hist same")
    keep += ratios
    rl = ROOT.TLegend(lpad.GetLeftMargin() + 0.03, 1 - lpad.GetTopMargin() - 0.02 - 1.4 * TEXT_PX / pad_px(lpad),
                      0.6, 1 - lpad.GetTopMargin() - 0.02)
    rl.SetTextFont(43)
    rl.SetTextSize(TEXT_PX - 4)
    rl.SetFillStyle(0)
    rl.AddEntry(stat, "MC stat. unc.", "f")
    rl.Draw()
    keep.append(rl)
    lpad.RedrawAxis()
    return keep + [upad, lpad]


# =================================================================================================================== #
#                                                       Output                                                        #
# =================================================================================================================== #
class Output:
    """
    Saves every canvas as pdf + png in outdir and, at close(), all of them into one multi-page pdf (book).
    The book is written at the end because ROOT can only have one pdf file open at a time.
    """

    def __init__(self, outdir, book):
        self.outdir = outdir
        self.book = os.path.join(outdir, book) if book else None
        self.pages = []
        os.makedirs(outdir, exist_ok=True)

    def save(self, canvas, name, keep=None, subdir=""):
        d = os.path.join(self.outdir, subdir)
        os.makedirs(d, exist_ok=True)
        for ext in ("pdf", "png"):
            canvas.SaveAs(os.path.join(d, "%s.%s" % (name, ext)))
        if self.book:
            self.pages.append((canvas, keep))   # canvas and its objects stay alive until the book is written

    def close(self):
        if not self.book or not self.pages:
            return
        for i, (canvas, _) in enumerate(self.pages):
            if i == 0:
                canvas.Print(self.book + "[", "pdf")
            canvas.Print(self.book, "pdf")
        self.pages[-1][0].Print(self.book + "]", "pdf")
        self.pages = []
        print("multi-page pdf: %s" % self.book)
