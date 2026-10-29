#!/usr/bin/env python3
"""
Build the Combine input ROOT file (shapes) for one lepton channel from the UHH2 Analysis output.

Inputs (histograms only -- the per-job SFrame outputs are read directly, no hadd needed):
  - Analysis job outputs   <ANALYSIS_DIR>/uhh2.AnalysisModuleRunner.<MC|DATA>.<version>_<n>.root
      folders <region>_SC_<variation>/<observable>  (written by ZprimeSemiLeptonicSystematicsModule)
  - Preselection job outputs (MC only, for the theory-variation normalisation)
      Input_General/sum_event_weights[_mcscale_<pt>|_isr_<dir>|_fsr_<dir>|_PDF_<i>]

Output: one flat ROOT file + a JSON summary read by buildDatacard.py
    <process>                       nominal template
    <process>_<nuisance>Up/Down     shape variations (only written when they differ from nominal)
    data_obs                        blinded: sum of all nominal MC templates | unblinded: real data
For the signal (TTbar) the NoSC_f template variation comes from the <region>_SC_NoSC_up/down folders.

Treatment of the variations
  experimental (pu, prefiring, b-tag, t-tag, lepton SFs)   : shape incl. normalisation effect
  theory (muR, muF, ISR, FSR, PDF) if NORMALISE_THEORY       : normalised to the inclusive effect of the variation
                                                             (cross-section uncertainties stay in the datacard lnN)
  muR / muF                                                  : separate nuisances, decorrelated between theory groups
  PDF (NNPDF3.1 Hessian, 100 eigenvectors)                    : sigma = sqrt(sum_i (h_i - h_0)^2) per bin, up/down = h_0 +- sigma
  nnloQCD, ewk                                                : ttbar processes only (identity elsewhere -> dropped automatically)

usage: python3 buildCombineInput.py --channel electron [--unblind] [--region Merged] [--observable cHel_P3n_Mtt800_Inf_cosThetaLT0p4] [--rebin 4]
Run inside the Combine environment (needs pyROOT only).
"""
import argparse
import glob
import json
import math
import os
import sys

import ROOT #type: ignore

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)


# =================================================================================================================== #
#                                                    Configuration                                                    #
# =================================================================================================================== #
YEAR = "UL18"
ANALYSIS_DIR     = "/data/dust/group/cms/zprime-uhh/Analysis_SC/UL18/{lep}/workdir_Analysis_UL18_{lep}/"
PRESELECTION_DIR = "/data/dust/group/cms/zprime-uhh/Presel_SC/UL18/workdir_Preselection_UL18/"
OUTPUT_DIR       = "/data/dust/user/ricardo/Combine/CMSSW_16_0_0/src/EntanglementAnalysis/UL18/{lep}/FullSyst/"

# Processes of the datacard -> sample versions (dataset_version without the year suffix) summed into each process
PROCESSES = {
    "TTbar":     ["TTToSemiLeptonic"],
    "ST":        ["ST_tW_antitop_5f_NoFullyHadronicDecays_PDFWeights", "ST_tW_top_5f_NoFullyHadronicDecays_PDFWeights",
                  "ST_t-channel_antitop_4f_InclusiveDecays", "ST_t-channel_top_4f_InclusiveDecays", "ST_s-channel_4f_leptonDecays"],
    "WJets":     ["WJetsToLNu_HT-70To100", "WJetsToLNu_HT-100To200", "WJetsToLNu_HT-200To400", "WJetsToLNu_HT-400To600",
                  "WJetsToLNu_HT-600To800", "WJetsToLNu_HT-800To1200", "WJetsToLNu_HT-1200To2500", "WJetsToLNu_HT-2500ToInf"],
    "TTTo2L2Nu": ["TTTo2L2Nu"],
    "Others":    ["WW", "WZ", "ZZ", "TTToHadronic",
                  "DYJetsToLL_M-50_HT-70to100", "DYJetsToLL_M-50_HT-100to200", "DYJetsToLL_M-50_HT-200to400", "DYJetsToLL_M-50_HT-400to600",
                  "DYJetsToLL_M-50_HT-600to800", "DYJetsToLL_M-50_HT-800to1200", "DYJetsToLL_M-50_HT-1200to2500", "DYJetsToLL_M-50_HT-2500toInf",
                  "QCD_HT50to100", "QCD_HT100to200", "QCD_HT200to300", "QCD_HT300to500", "QCD_HT500to700", "QCD_HT700to1000",
                  "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf"],
}
SIGNAL = "TTbar"
DATA_SAMPLES = {
    "muon":     ["DATA_SingleMuon_RunA", "DATA_SingleMuon_RunB", "DATA_SingleMuon_RunC", "DATA_SingleMuon_RunD"],
    "electron": ["DATA_EGamma_RunA", "DATA_EGamma_RunB", "DATA_EGamma_RunC", "DATA_EGamma_RunD"],
}

# Theory nuisances (muR, muF, ISR, FSR) are correlated within a group and decorrelated between groups
THEORY_GROUP = {"TTbar": "ttbar", "TTTo2L2Nu": "ttbar", "ST": "ST", "WJets": "WJets", "Others": "Others"}

# Experimental weight variations: folder name in the Analysis output -> nuisance name in the datacard
EXPERIMENTAL = {
    "common":   ["pu", "prefiring", "btag_cferr1", "btag_cferr2", "btag_hf", "btag_hfstats1", "btag_hfstats2", "btag_lf",
                 "btag_lfstats1", "btag_lfstats2", "ttag_corr", "ttag_uncorr", "tmistag", "nnloQCD", "ewk"],
    "muon":     ["mu_id_stat", "mu_id_syst", "mu_iso_stat", "mu_iso_syst", "mu_reco", "mu_trigger_stat", "mu_trigger_syst"],
    "electron": ["ele_id", "ele_reco", "ele_trigger"],
}
# nuisances that are uncorrelated between data-taking years get a year suffix (e.g. btag_hfstats1_UL18)
YEAR_DECORRELATED = {"btag_hfstats1", "btag_hfstats2", "btag_lfstats1", "btag_lfstats2", "ttag_uncorr", "tmistag",
                     "mu_id_stat", "mu_iso_stat", "mu_trigger_stat", "ele_trigger"}

# Separate-run variations (JEC/JER): nuisance -> Analysis directory of the varied run (None = not produced yet, skipped)
SEPARATE_RUNS = {
    "jec": {"Up": None, "Down": None},   # e.g. ".../Analysis_SC/UL18/{lep}_JEC_up/workdir_Analysis_UL18_{lep}_JEC_up/"
    "jer": {"Up": None, "Down": None},
}

# Theory normalisation: Analysis variation folder -> Preselection Input_General histogram with the SAME generator weight.
# NB: PreselectionHists labels muR/muF the other way round from MCScaleVariation (which follows the POWHEG LHE order
# 1: muR=1,muF=2 | 2: 1,0.5 | 3: 2,1 | 6: 0.5,1), hence the crossed pairs below. Check against the cluster versions
# of MCWeight.cxx (MCScaleVariation) and PSWeights.h before trusting the normalisation.
NORMALISE_THEORY = True
PRESEL_SUMS = {
    "murmuf_upnone":   "sum_event_weights_mcscale_noneup",    # systweight 3: muR x2
    "murmuf_downnone": "sum_event_weights_mcscale_nonedown",  # systweight 6: muR x0.5
    "murmuf_noneup":   "sum_event_weights_mcscale_upnone",    # systweight 1: muF x2
    "murmuf_nonedown": "sum_event_weights_mcscale_downnone",  # systweight 2: muF x0.5
    "isr_up":   "sum_event_weights_isr_up",
    "isr_down": "sum_event_weights_isr_down",
    "fsr_up":   "sum_event_weights_fsr_up",
    "fsr_down": "sum_event_weights_fsr_down",
}
PRESEL_SUMS.update({"pdf%d" % i: "sum_event_weights_PDF_%d" % (i + 1) for i in range(100)})

# Theory nuisances built from the folders: nuisance stem -> (up folder, down folder)
THEORY_UPDOWN = {"muR": ("murmuf_upnone", "murmuf_downnone"),
                 "muF": ("murmuf_noneup", "murmuf_nonedown"),
                 "isr": ("isr_up", "isr_down"),
                 "fsr": ("fsr_up", "fsr_down")}

IDENTITY_TOLERANCE = 1e-6  # a variation within this relative difference of nominal in every bin is not written


# =================================================================================================================== #
#                                                       Helpers                                                       #
# =================================================================================================================== #
def job_files(directory, kind, version):
    """per-job SFrame outputs of one sample: uhh2.AnalysisModuleRunner.<kind>.<version>_<n>.root"""
    files = sorted(glob.glob(os.path.join(directory, "uhh2.AnalysisModuleRunner.%s.%s_%s_[0-9]*.root" % (kind, version, YEAR))))
    if not files:
        raise RuntimeError("no job outputs for %s %s_%s in %s" % (kind, version, YEAR, directory))
    return files


def sum_histograms(files, paths, allow_missing=False):
    """sum the histograms at 'paths' over all files -> {path: TH1}"""
    out = {}
    for fname in files:
        f = ROOT.TFile.Open(fname)
        if (not f) or f.IsZombie():
            raise RuntimeError("cannot open " + fname)
        for path in paths:
            h = f.Get(path)
            if not h:
                if allow_missing:
                    continue
                raise RuntimeError("missing %s in %s" % (path, fname))
            if path in out:
                out[path].Add(h)
            else:
                out[path] = h.Clone("sum_" + path.replace("/", "__"))
                out[path].SetDirectory(0)
        f.Close()
    return out


def rebinned(h, factor, name):
    if h.GetNbinsX() % factor != 0:
        raise RuntimeError("%s has %d bins, not divisible by %d" % (h.GetName(), h.GetNbinsX(), factor))
    r = h.Rebin(factor, name)
    r.SetDirectory(0)
    return r


def clip_negative(h, label):
    """Combine needs non-negative templates: negative bins (from negative generator weights) are set to 0"""
    for b in range(1, h.GetNbinsX() + 1):
        if h.GetBinContent(b) < 0:
            print("  WARNING: %s bin %d = %.3g < 0 -> set to 0" % (label, b, h.GetBinContent(b)))
            h.SetBinContent(b, 0.)


def differs(h, nominal):
    for b in range(1, h.GetNbinsX() + 1):
        n, v = nominal.GetBinContent(b), h.GetBinContent(b)
        if abs(v - n) > IDENTITY_TOLERANCE * max(abs(n), 1e-9):
            return True
    return False


def nuisance_name(stem):
    return stem + "_" + YEAR if stem in YEAR_DECORRELATED else stem


# =================================================================================================================== #
#                                                         Main                                                        #
# =================================================================================================================== #
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--channel", required=True, choices=["muon", "electron"])
    parser.add_argument("--unblind", action="store_true", help="use real data for data_obs (default: blinded, data_obs = sum of MC)")
    parser.add_argument("--region", default="Merged", help="Merged or Chi2cut_PASS")
    parser.add_argument("--observable", default="cHel_P3n_Mtt800_Inf_cosThetaLT0p4", help="histogram name inside the <region>_SC_* folders")
    parser.add_argument("--rebin", type=int, default=4, help="merge this many of the 24 bins (4 -> 6 bins)")
    parser.add_argument("--output", default=None, help="output ROOT file (default: OUTPUT_DIR/<observable>_<region>_<YEAR>_<channel>.root)")
    args = parser.parse_args()

    lep = args.channel
    ana_dir = ANALYSIS_DIR.format(lep=lep)
    out_path = args.output or os.path.join(OUTPUT_DIR.format(lep=lep), "%s_%s_%s_%s.root" % (args.observable, args.region, YEAR, lep))
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    folder = lambda variation: "%s_SC_%s/%s" % (args.region, variation, args.observable)

    weight_variations = EXPERIMENTAL["common"] + EXPERIMENTAL[lep]
    folders = ["nominal"] + ["%s_%s" % (v, d) for v in weight_variations for d in ("up", "down")]
    folders += [f for pair in THEORY_UPDOWN.values() for f in pair] + ["pdf%d" % i for i in range(100)]

    templates = {}   # histogram name in the output file -> TH1
    shapes = {}      # process -> list of shape nuisances written
    norm_report = []

    for proc, versions in PROCESSES.items():
        print("== %s (%d samples)" % (proc, len(versions)))
        group = THEORY_GROUP[proc]
        proc_hists = {}  # folder -> process-level histogram (normalised per sample where needed)

        for version in versions:
            files = job_files(ana_dir, "MC", version)
            paths = [folder(f) for f in folders]
            if proc == SIGNAL:
                paths += [folder("NoSC_up"), folder("NoSC_down")]
            sample = sum_histograms(files, paths)

            # theory normalisation factors from the inclusive Preselection sums of this sample
            norm = {}
            if NORMALISE_THEORY:
                presel = sum_histograms(job_files(PRESELECTION_DIR, "MC", version),
                                        ["Input_General/sum_event_weights"] + ["Input_General/" + s for s in set(PRESEL_SUMS.values())],
                                        allow_missing=True)
                n_nom = presel["Input_General/sum_event_weights"].GetBinContent(1)
                for var, sname in PRESEL_SUMS.items():
                    h = presel.get("Input_General/" + sname)
                    n_var = h.GetBinContent(1) if h else 0.
                    norm[var] = n_nom / n_var if (n_var != 0 and n_nom != 0) else 1.
                for var in ("murmuf_upnone", "murmuf_noneup", "isr_up", "fsr_up"):
                    norm_report.append("%-55s %-15s %.4f" % (version, var, norm[var]))

            for path, h in sample.items():
                var = path.split("/")[0][len(args.region) + 4:]   # strip "<region>_SC_"
                h.Scale(norm.get(var, 1.))
                if var in proc_hists:
                    proc_hists[var].Add(h)
                else:
                    proc_hists[var] = h.Clone("proc_%s_%s" % (proc, var))
                    proc_hists[var].SetDirectory(0)

        # rebin everything first (the PDF sigma must be built from the final bins)
        proc_hists = {var: rebinned(h, args.rebin, "rb_%s_%s" % (proc, var)) for var, h in proc_hists.items()}

        # ---- nominal ----
        nominal = proc_hists["nominal"].Clone(proc)
        clip_negative(nominal, proc)
        templates[proc] = nominal
        shapes[proc] = []

        def add_shape(nuis, h_up, h_down):
            up, down = h_up.Clone("%s_%sUp" % (proc, nuis)), h_down.Clone("%s_%sDown" % (proc, nuis))
            clip_negative(up, up.GetName()); clip_negative(down, down.GetName())
            if differs(up, nominal) or differs(down, nominal):
                templates[up.GetName()] = up
                templates[down.GetName()] = down
                shapes[proc].append(nuis)

        # ---- experimental weight variations ----
        for v in weight_variations:
            add_shape(nuisance_name(v), proc_hists[v + "_up"], proc_hists[v + "_down"])

        # ---- theory: muR, muF, ISR, FSR (per theory group) ----
        for stem, (f_up, f_down) in THEORY_UPDOWN.items():
            add_shape("%s_%s" % (stem, group), proc_hists[f_up], proc_hists[f_down])

        # ---- PDF: Hessian sum over the 100 eigenvector variations ----
        h0 = proc_hists["nominal"]
        pdf_up, pdf_down = h0.Clone("pdf_up_" + proc), h0.Clone("pdf_down_" + proc)
        for b in range(0, h0.GetNbinsX() + 2):
            sigma = math.sqrt(sum((proc_hists["pdf%d" % i].GetBinContent(b) - h0.GetBinContent(b)) ** 2 for i in range(100)))
            pdf_up.SetBinContent(b, h0.GetBinContent(b) + sigma)
            pdf_down.SetBinContent(b, h0.GetBinContent(b) - sigma)
        add_shape("pdf", pdf_up, pdf_down)

        # ---- signal-shape variation ----
        if proc == SIGNAL:
            add_shape("NoSC_f", proc_hists["NoSC_up"], proc_hists["NoSC_down"])

        # ---- separate-run variations (JEC/JER), when available ----
        for nuis, dirs in SEPARATE_RUNS.items():
            if not (dirs["Up"] and dirs["Down"]):
                continue
            varied = {}
            for d in ("Up", "Down"):
                h_sum = None
                for version in versions:
                    h = sum_histograms(job_files(dirs[d].format(lep=lep), "MC", version), [folder("nominal")])[folder("nominal")]
                    if h_sum is None: h_sum = h
                    else: h_sum.Add(h)
                varied[d] = rebinned(h_sum, args.rebin, "rb_%s_%s_%s" % (proc, nuis, d))
            add_shape(nuis, varied["Up"], varied["Down"])

    # ---- data_obs ----
    if args.unblind:
        data_hist = None
        for version in DATA_SAMPLES[lep]:
            h = sum_histograms(job_files(ana_dir, "DATA", version), [folder("nominal")])[folder("nominal")]
            if data_hist is None: data_hist = h
            else: data_hist.Add(h)
        data_obs = rebinned(data_hist, args.rebin, "data_obs")
    else:
        data_obs = templates[SIGNAL].Clone("data_obs")
        for proc in PROCESSES:
            if proc != SIGNAL:
                data_obs.Add(templates[proc])
    templates["data_obs"] = data_obs

    # ---- write ----
    fout = ROOT.TFile(out_path, "RECREATE")
    for name in sorted(templates):
        templates[name].SetTitle(name)
        templates[name].Write(name)
    ROOT.TNamed("blinded", "false" if args.unblind else "true (data_obs = sum of MC)").Write()
    fout.Close()

    summary = {"channel": lep, "year": YEAR, "bin": "%s_%s" % (YEAR, lep), "root_file": os.path.abspath(out_path),
               "blinded": not args.unblind, "signal": SIGNAL, "processes": list(PROCESSES), "shapes": shapes,
               "region": args.region, "observable": args.observable}
    with open(os.path.splitext(out_path)[0] + ".json", "w") as f:
        json.dump(summary, f, indent=2)

    # ---- report ----
    print("\nwrote %s (%s)" % (out_path, "UNBLINDED" if args.unblind else "blinded: data_obs = sum of MC"))
    print("\n%-14s %10s  shape nuisances" % ("process", "integral"))
    for proc in PROCESSES:
        print("%-14s %10.2f  %d" % (proc, templates[proc].Integral(), len(shapes[proc])))
    print("%-14s %10.2f" % ("data_obs", data_obs.Integral()))
    if norm_report:
        print("\ntheory normalisation factors (inclusive nominal / varied; expect values close to 1):")
        print("\n".join(norm_report))


if __name__ == "__main__":
    sys.exit(main())
