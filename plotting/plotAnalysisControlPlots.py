#!/usr/bin/env python3
"""
Analysis selection control plots: stacked MC (+ data) of the distributions in plotConfig.ANALYSIS_KINEMATICS, filled at
one step of the Analysis selection (one <step>_General folder), with the MC stat. uncertainty.

Data is only drawn with --data. Generator-level quantities (plotConfig.GEN_ONLY) never show data, and the fit
observables (plotConfig.BLINDED) only with --unblind.

usage: python3 plotAnalysisControlPlots.py [--channel muon] [--step EndOfBaselineSelection_General] [--data] [--logy]
       python3 plotAnalysisControlPlots.py --list-steps
"""
import argparse
import os

import ROOT

import plotConfig as cfg
import plotUtils as pu
from plotPreselectionKinematics import print_steps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--channel", nargs="*", default=cfg.CHANNELS, choices=cfg.CHANNELS)
    parser.add_argument("--year", default="UL18", choices=sorted(cfg.LUMI_PERIOD))
    parser.add_argument("--step", default="EndOfBaselineSelection_General", help="Analysis folder the distributions are taken from")
    parser.add_argument("--list-steps", action="store_true", help="print the available steps and exit")
    parser.add_argument("--data", action="store_true", help="draw data (default: simulation only)")
    parser.add_argument("--unblind", action="store_true", help="also draw data for the fit observables in plotConfig.BLINDED")
    parser.add_argument("--logy", action="store_true", help="logarithmic y-axis (default: linear)")
    parser.add_argument("--plots", nargs="*", default=None, help="subset of the distributions (default: all)")
    parser.add_argument("--no-overflow", action="store_true", help="never add under/overflow to the first/last bin "
                        "(default: added for open-ended quantities only, see plotConfig.fold_overflow)")
    parser.add_argument("--input-dir", default=None, help="default: plotConfig.ANALYSIS_DIR ({lep} is replaced)")
    parser.add_argument("--outdir", default=None, help="default: plotConfig.OUTPUT_DIR/AnalysisControlPlots/<channel>/<step>")
    parser.add_argument("--rebuild", action="store_true", help="re-read the job outputs instead of the cache")
    args = parser.parse_args()

    pu.setup_style()
    cells = cfg.cells()
    names = args.plots or cfg.ANALYSIS_KINEMATICS
    for lep in args.channel:
        indir = (args.input_dir or cfg.ANALYSIS_DIR).format(year=args.year, lep=lep)
        if args.list_steps:
            print_steps(indir, cells, args.year)
            return
        outdir = (args.outdir.format(lep=lep) if args.outdir else
                  os.path.join(cfg.OUTPUT_DIR.format(year=args.year), "AnalysisControlPlots", lep, args.step))
        paths = ["%s/%s" % (args.step, n) for n in names]
        print("Analysis control plots, %s %s, step %s, %d distributions%s -> %s"
              % (args.year, lep, args.step, len(names), ", with data" if args.data else "", outdir))
        sources = [(indir, args.year, cells, "MC", paths)]
        if args.data:
            sources.append((indir, args.year, {"data": cfg.combine.DATA_SAMPLES[lep]}, "DATA", paths))
        hists = pu.cached(os.path.join(outdir, "cache_%s.root" % ("data" if args.data else "mc")), sources, args.rebuild)
        hists.setdefault("data", {})

        tag = "".join(x for x in ("_data" if args.data else "", "_logy" if args.logy else ""))
        out = pu.Output(outdir, "AnalysisControlPlots_%s_%s_%s%s.pdf" % (args.year, lep, args.step, tag))
        text = ["%s, %s" % (cfg.CHANNEL_LABEL[lep], args.step.replace("_General", ""))]
        failed = []
        for name in names:
            path = "%s/%s" % (args.step, name)
            fold = cfg.fold_overflow(name) and not args.no_overflow
            binning = cfg.binning(name, "analysis")
            try:
                stack = []
                for label, colour, _ in cfg.STACK_GROUPS:
                    h = pu.total(hists, [c for c in cells if cfg.cell_group(c) == label], path, "%s_%s_%s" % (label, lep, name))
                    if h is not None:
                        stack.append((label, colour, pu.prepare(h, binning, fold)))
                data = None
                if args.data and name not in cfg.GEN_ONLY:
                    if name in cfg.BLINDED and not args.unblind:
                        print("  %s: blinded, data not drawn (--unblind to show it)" % name)
                    elif path in hists["data"]:
                        data = pu.prepare(hists["data"][path], binning, fold)
            except pu.BinningError as e:
                print("  ERROR %s: %s (plot skipped)" % (name, e))
                failed.append(name)
                continue
            if not stack:
                print("  %s: not found, skipped" % path)
                continue
            c = ROOT.TCanvas("c_%s_%s" % (lep, name), "", 800, 800)
            keep = pu.draw_stack_plot(c, stack, data, pu.x_title(stack[0][2]), pu.y_title(stack[0][2], name), args.logy, text, args.year)
            out.save(c, name + tag, keep=(c, keep))
        out.close()
        if failed:
            print("binning in plotConfig.binning() does not match the booked binning for: %s" % ", ".join(failed))


if __name__ == "__main__":
    main()
