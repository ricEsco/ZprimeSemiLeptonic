#!/usr/bin/env python3
"""
Systematic variations vs nominal for the SC observables of the <region>_SC_* folders (MC only).

One plot per region x observable x nuisance, with the nuisance names of the Combine input. By default the MC sum of all
processes is shown; a nuisance that acts on some processes only (muR_ttbar, NoSC_f, ...) shifts those and keeps the
others at nominal, as in the fit. The histograms are the raw folder contents (no normalisation of the theory
variations), rebinned like the Combine input. PDF: Hessian +-1 sigma from the 100 replicas.

Under/overflow: by default only added to the first/last bin for open-ended quantities (plotConfig.fold_overflow), so
the sentinel values (-10 / 99) of events outside the phase space of an observable stay out of the plot.

usage: python3 plotSystematicVariations.py [--channel muon] [--region Merged] [--observables cos_theta1k Ckk]
                                           [--nuisances 'btag_*' pdf] [--process TTbar] [--no-overflow] [--year UL18]
"""
import argparse
import fnmatch
import os

import ROOT

import plotConfig as cfg
import plotUtils as pu


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--channel", nargs="*", default=cfg.CHANNELS, choices=cfg.CHANNELS)
    parser.add_argument("--year", default="UL18", choices=sorted(cfg.LUMI_PERIOD))
    parser.add_argument("--region", nargs="*", default=cfg.ANALYSIS_REGIONS, choices=cfg.ANALYSIS_REGIONS)
    parser.add_argument("--observables", nargs="*", default=None, help="subset / wildcards (default: all SC observables)")
    parser.add_argument("--nuisances", nargs="*", default=None, help="subset / wildcards of Combine nuisance names (default: all)")
    parser.add_argument("--process", default=None, choices=list(cfg.combine.PROCESSES), help="one datacard process (default: sum of all MC)")
    parser.add_argument("--no-overflow", action="store_true", help="never add under/overflow to the first/last bin "
                        "(default: added for open-ended quantities only, see plotConfig.fold_overflow)")
    parser.add_argument("--input-dir", default=None, help="default: plotConfig.ANALYSIS_DIR ({lep} is replaced)")
    parser.add_argument("--outdir", default=None, help="default: plotConfig.OUTPUT_DIR/SystematicVariations/<channel>")
    parser.add_argument("--rebuild", action="store_true", help="re-read the job outputs instead of the cache")
    args = parser.parse_args()

    def selected(name, patterns):
        return patterns is None or any(fnmatch.fnmatch(name, p) for p in patterns)

    pu.setup_style()
    all_cells = cfg.cells()
    cells = [c for c in all_cells if args.process is None or cfg.cell_process(c) == args.process]
    scope = args.process or "all MC"
    for lep in args.channel:
        indir = (args.input_dir or cfg.ANALYSIS_DIR).format(year=args.year, lep=lep)
        outdir = args.outdir.format(lep=lep) if args.outdir else os.path.join(cfg.OUTPUT_DIR.format(year=args.year), "SystematicVariations", lep)
        nuis_list = [n for n in cfg.nuisances(lep, args.year) if selected(n["name"], args.nuisances)
                     and (args.process is None or n["procs"] is None or args.process in n["procs"])]
        ref = cfg.combine.PROCESSES[cfg.combine.SIGNAL][0]
        obs = {r: [o for o in pu.list_folder(indir, args.year, ref, "%s_SC_nominal" % r) if selected(o, args.observables)]
               for r in args.region}
        paths = [p for r in args.region for o in obs[r] for p in pu.sc_paths(r, o, nuis_list)]
        print("Systematic variations, %s %s (%s): %s, %d nuisances -> %s"
              % (args.year, lep, scope, ", ".join("%s: %d observables" % (r, len(obs[r])) for r in args.region), len(nuis_list), outdir))
        hists = pu.cached(os.path.join(outdir, "cache.root"), [(indir, args.year, {c: all_cells[c] for c in cells}, "MC", paths)],
                          args.rebuild)

        n_plots = 0
        for region in args.region:
            for o in obs[region]:
                binning, fold = cfg.binning(o, "sc"), cfg.fold_overflow(o) and not args.no_overflow

                def prep(h):
                    return pu.prepare(h, binning, fold, rebin=cfg.SC_REBIN)

                nominal = pu.total(hists, cells, "%s_SC_nominal/%s" % (region, o), "nominal")
                if nominal is None:
                    continue
                try:
                    nominal = prep(nominal)
                except pu.BinningError as e:
                    print("  ERROR %s %s: %s (skipped)" % (region, o, e))
                    continue
                out = pu.Output(os.path.join(outdir, region, o), "%s_%s_%s_%s%s.pdf" % (args.year, lep, region, o, "_" + args.process if args.process else ""))
                for n in nuis_list:
                    pair = pu.variation_pair(hists, cells, region, o, n, prep)
                    if pair is None:
                        continue
                    c = ROOT.TCanvas("c_%s_%s_%s_%s" % (lep, region, o, n["name"]), "", 800, 800)
                    text = ["%s, %s" % (cfg.CHANNEL_LABEL[lep], cfg.REGION_LABEL.get(region, region)), "Process: %s" % scope]
                    keep = pu.draw_variation(c, nominal.Clone("nom_" + n["name"]), pair[0], pair[1], n["name"], pu.x_title(nominal), args.year, text)
                    out.save(c, n["name"] + ("_" + args.process if args.process else ""), keep=(c, keep, pair))
                    n_plots += 1
                out.close()
        print("%d plots written" % n_plots)


if __name__ == "__main__":
    main()
