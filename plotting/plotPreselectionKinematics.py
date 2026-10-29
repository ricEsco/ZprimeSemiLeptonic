#!/usr/bin/env python3
"""
Preselection kinematics: stacked MC of the distributions in plotConfig.PRESELECTION_KINEMATICS, filled at one step of
the Preselection (one <step>_General folder). Simulation only (no lepton-channel split, data not drawn).

usage: python3 plotPreselectionKinematics.py [--step MET_General] [--logy] [--no-overflow] [--year UL18] [--plots pt_jet MET]
       python3 plotPreselectionKinematics.py --list-steps
"""
import argparse
import os

import ROOT

import plotConfig as cfg
import plotUtils as pu


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", default="UL18", choices=sorted(cfg.LUMI_PERIOD))
    parser.add_argument("--step", default="MET_General", help="Preselection folder the distributions are taken from")
    parser.add_argument("--list-steps", action="store_true", help="print the available steps and exit")
    parser.add_argument("--logy", action="store_true", help="logarithmic y-axis (default: linear)")
    parser.add_argument("--plots", nargs="*", default=None, help="subset of the distributions (default: all)")
    parser.add_argument("--no-overflow", action="store_true", help="never add under/overflow to the first/last bin "
                        "(default: added for open-ended quantities only, see plotConfig.fold_overflow)")
    parser.add_argument("--input-dir", default=None, help="default: plotConfig.PRESELECTION_DIR")
    parser.add_argument("--outdir", default=None, help="default: plotConfig.OUTPUT_DIR/Preselection/<step>")
    parser.add_argument("--rebuild", action="store_true", help="re-read the job outputs instead of the cache")
    args = parser.parse_args()

    indir = args.input_dir or cfg.PRESELECTION_DIR.format(year=args.year)
    groups = cfg.cells()

    if args.list_steps:
        print_steps(indir, groups, args.year)
        return

    names = args.plots or cfg.PRESELECTION_KINEMATICS
    outdir = args.outdir or os.path.join(cfg.OUTPUT_DIR.format(year=args.year), "Preselection", args.step)
    print("Preselection kinematics, %s, step %s, %d distributions -> %s" % (args.year, args.step, len(names), outdir))

    paths = ["%s/%s" % (args.step, n) for n in names]
    hists = pu.cached(os.path.join(outdir, "cache_%s.root" % args.step), [(indir, args.year, groups, "MC", paths)], args.rebuild)

    pu.setup_style()
    tag = "_logy" if args.logy else ""
    out = pu.Output(outdir, "Preselection_%s_%s%s.pdf" % (args.year, args.step, tag))
    text = ["Preselection, after %s" % args.step.replace("_General", "")]
    failed = []
    for name in names:
        path = "%s/%s" % (args.step, name)
        fold = cfg.fold_overflow(name) and not args.no_overflow
        try:
            stack = []
            for label, colour, _ in cfg.STACK_GROUPS:
                h = pu.total(hists, [c for c in groups if cfg.cell_group(c) == label], path, "%s_%s" % (label, name))
                if h is not None:
                    stack.append((label, colour, pu.prepare(h, cfg.binning(name, "preselection"), fold)))
        except pu.BinningError as e:
            print("  ERROR %s: %s (plot skipped)" % (name, e))
            failed.append(name)
            continue
        if not stack:
            print("  %s: not found, skipped" % path)
            continue
        c = ROOT.TCanvas("c_" + name, "", 800, 800)
        keep = pu.draw_stack_plot(c, stack, None, pu.x_title(stack[0][2]), pu.y_title(stack[0][2], name), args.logy, text, args.year)
        out.save(c, name + tag, keep=(c, keep))
    out.close()
    if failed:
        print("binning in plotConfig.binning() does not match the booked binning for: %s" % ", ".join(failed))


def print_steps(indir, groups, year):
    first = next(f for v in sum(groups.values(), []) for f in pu.job_files(indir, "MC", v, year))
    f = ROOT.TFile.Open(first)
    print("Preselection steps in %s:" % os.path.basename(first))
    for key in f.GetListOfKeys():
        if key.GetClassName() == "TDirectoryFile" and key.GetName().endswith("_General"):
            print("  " + key.GetName())
    f.Close()


if __name__ == "__main__":
    main()
