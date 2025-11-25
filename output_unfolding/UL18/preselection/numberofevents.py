# find /path/outputs -type f -name "*.root" > files.txt
# python numberofevents.py --input-list files.txt --csv totals.csv

from __future__ import print_function

import os, sys, argparse
from collections import defaultdict
import ROOT
ROOT.gROOT.SetBatch(True)

CHANNELS = ["SL", "DL"]
MT_BINS  = ["0_350", "350_500", "500_750", "750_1000", "1000_1500", "1500Inf"]

def read_file_list(list_path):
    files = []
    with open(list_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            files.append(line)
    return files

def open_hist(tf, folder, name):
    obj = tf.Get("{}/{}".format(folder, name))
    if obj and obj.InheritsFrom("TH1"):
        return obj
    return None

def main():
    ap = argparse.ArgumentParser(
        description="Sum Nneg/Npos/Events from dYgen histograms over many ROOT files (no errors).")
    ap.add_argument("--input-list", help="Text file with one ROOT file path per line.")
    ap.add_argument("files", nargs="*", help="ROOT files (you can pass many).")
    ap.add_argument("--csv", help="Optional CSV output path.")
    ap.add_argument("--no-merge-edges", action="store_true",
                dest="no_merge_edges",
                help="Do NOT merge under/overflow into neg/pos (default merges).")
    args = ap.parse_args()

    files = []
    if args.input_list:
        files.extend(read_file_list(args.input_list))
    files.extend(args.files)

    if not files:
        print("No input ROOT files provided. Use --input-list or pass files as args.")
        sys.exit(1)

    merge_edges = (not args.no_merge_edges)

    # Accumulator: key=(channel, mtbin) -> dict of sums
    acc = defaultdict(lambda: {'Nneg':0.0, 'Npos':0.0, 'Nev':0.0})

    n_ok, n_bad = 0, 0
    for path in files:
        try:
            tf = ROOT.TFile.Open(path, "READ")
            if not tf or tf.IsZombie():
                print("WARN: cannot open", path)
                n_bad += 1
                continue

            for ch in CHANNELS:
                for mb in MT_BINS:
                    folder = "{}_DY_{}_General".format(ch, mb)
                    h_dy = open_hist(tf, folder, "dYgen")
                    h_ev = open_hist(tf, folder, "Events")

                    # dYgen - Nneg (bin1), Npos (bin2); optionally merge edge bins
                    if h_dy:
                        nb = h_dy.GetNbinsX()  # should be 2
                        b1 = h_dy.GetBinContent(1)
                        b2 = h_dy.GetBinContent(2)
                        if merge_edges:
                            under = h_dy.GetBinContent(0)
                            over  = h_dy.GetBinContent(nb+1)
                            b1 += under
                            b2 += over
                        acc[(ch, mb)]['Nneg'] += b1
                        acc[(ch, mb)]['Npos'] += b2

                    # Events - bin1
                    if h_ev:
                        acc[(ch, mb)]['Nev'] += h_ev.GetBinContent(1)

            tf.Close()
            n_ok += 1
        except Exception as e:
            print("WARN: failed on", path, "->", e)
            n_bad += 1

    # Print summary table
    print("\nProcessed files: ok={} bad={}\n".format(n_ok, n_bad))
    header = "{:4s}  {:>11s}  {:>15s}  {:>15s}  {:>15s}".format(
        "CH", "mtt bin", "Nneg", "Npos", "Events"
    )
    print(header)
    print("-"*len(header))

    # Deterministic order
    lines_csv = []
    if args.csv:
        lines_csv.append("channel,mtt_bin,Nneg,Npos,Events")

    for ch in CHANNELS:
        for mb in MT_BINS:
            slot = acc.get((ch, mb))
            if not slot:
                continue
            Nneg = slot['Nneg']
            Npos = slot['Npos']
            Nev  = slot['Nev']
            print("{:4s}  {:>11s}  {:15.6f}  {:15.6f}  {:15.6f}".format(ch, mb, Nneg, Npos, Nev))
            if args.csv:
                lines_csv.append("{},{},{:.6f},{:.6f},{:.6f}".format(ch, mb, Nneg, Npos, Nev))

    if args.csv:
        with open(args.csv, "w") as f:
            f.write("\n".join(lines_csv))
        print("\nWrote CSV to {}".format(args.csv))

if __name__ == "__main__":
    main()
