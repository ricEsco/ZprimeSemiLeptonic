#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Merge the inclusive alpha strength maps from the Preselection job outputs into one small ROOT file.

ZprimePreselectionModule writes StrengthMaps/alpha_<Obs> (TProfile2D) into every TTToSemiLeptonic job output.
Instead of hadd-ing the full outputs (trees included, ~TB), this reads only the StrengthMaps/ folder of each
job file and merges the profiles with TProfile2D::Merge (what hadd does), so bin means and errors are exact.

The output keeps the same layout (StrengthMaps/alpha_<Obs>), so it can be used directly as the
StrengthMapFile config item of the Analysis (EntanglementNoSCReweighting).

usage: python mergeUHH2PreselectionStrengthMaps.py "<path>/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_*.root" StrengthMaps_UL18.root
       (quote the pattern so the shell does not expand it)
Runs with python2.7/ROOT 6.14 (CMSSW_10_6) and python3.
"""
from __future__ import print_function
import sys, glob
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

FOLDER = "StrengthMaps"

def main(pattern, out_path):
    files = sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError("no files match " + pattern)
    print("Merging %s/ from %d files" % (FOLDER, len(files)))

    merged = {}      # map name -> merged TProfile2D
    n_used = 0
    for i, path in enumerate(files):
        f = ROOT.TFile.Open(path)
        if (not f) or f.IsZombie():
            raise RuntimeError("cannot open " + path)
        folder = f.Get(FOLDER)
        if not folder:
            raise RuntimeError("no %s/ folder in %s (is this a TTToSemiLeptonic Preselection output?)" % (FOLDER, path))
        names = sorted(key.GetName() for key in folder.GetListOfKeys())
        if merged and names != sorted(merged):
            raise RuntimeError("map list in %s differs from the first file" % path)
        for name in names:
            prof = folder.Get(name)
            if name not in merged:
                merged[name] = prof.Clone(name)        # first file: keep a copy
                merged[name].SetDirectory(0)
            else:
                others = ROOT.TList()
                others.Add(prof)
                merged[name].Merge(others)             # same merge as hadd
        f.Close()
        n_used += 1
        if (i + 1) % 20 == 0 or i + 1 == len(files):
            print("  %d/%d files" % (i + 1, len(files)))

    out = ROOT.TFile(out_path, "RECREATE")
    d = out.mkdir(FOLDER)
    d.cd()
    for name in sorted(merged):
        merged[name].Write()
    out.cd()
    ROOT.TNamed("merged_from", "%d files matching %s" % (n_used, pattern)).Write()
    out.Close()

    ref = merged[sorted(merged)[0]]
    print("Wrote %d maps to %s:%s/ (%s has %d entries)" % (len(merged), out_path, FOLDER, ref.GetName(), int(ref.GetEntries())))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
