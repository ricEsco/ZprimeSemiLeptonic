#!/usr/bin/env python3
"""
Write the Combine datacard for one lepton channel from the JSON summary produced by buildCombineInput.py,
so the datacard always lists exactly the shape variations that exist in the ROOT file.

  - rate uncertainties (lnN): luminosity and cross sections (configured below)
  - shape uncertainties: every <process>_<nuisance>Up/Down template written by buildCombineInput.py
  - NoSC_f: signal-shape nuisance, unconstrained (flatParam)
  - autoMCStats for the MC statistical uncertainty

usage: python3 buildDatacard.py <summary.json> [--output datacard.txt]
"""
import argparse
import json
import os


# =================================================================================================================== #
#                                                    Configuration                                                    #
# =================================================================================================================== #
# rate (lnN) uncertainties: {nuisance: {process: value}}, processes not listed get "-"
LNN = {
    "lumi_13TeV_15161718_I": {"TTbar": 1.0084, "ST": 1.0084, "WJets": 1.0084, "TTTo2L2Nu": 1.0084, "Others": 1.0084},
    "xsec_ttbar":            {"TTbar": 1.05},
    "xsec_ST":               {"ST": 1.20},
    "xsec_WJets":            {"WJets": 1.20},
    "xsec_TTTo2L2Nu":        {"TTTo2L2Nu": 1.20},
    "xsec_Others":           {"Others": 1.20},
}
FLAT_PARAMS = ["NoSC_f"]   # shape nuisances left unconstrained
AUTO_MC_STATS = 0          # autoMCStats threshold

NAME_W, TYPE_W, COL_W = 30, 10, 20   # column widths (40 characters before the first process column)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", help="JSON summary written by buildCombineInput.py")
    parser.add_argument("--output", default=None, help="datacard path (default: next to the ROOT file)")
    args = parser.parse_args()

    with open(args.summary) as f:
        s = json.load(f)
    procs, bin_name, root_file = s["processes"], s["bin"], s["root_file"]
    out_path = args.output or os.path.join(os.path.dirname(root_file), "%s%s_datacard.txt" % (s["year"], s["channel"]))
    # signal first (process id 0), backgrounds > 0
    procs = [s["signal"]] + [p for p in procs if p != s["signal"]]

    # all shape nuisances, in first-seen order
    shape_nuisances = []
    for p in procs:
        for n in s["shapes"][p]:
            if n not in shape_nuisances:
                shape_nuisances.append(n)

    sep = "-" * 150
    lines = ["# Created using buildDatacard.py from %s" % os.path.basename(args.summary),
             "# %s, %s region, observable %s, %s" % (bin_name, s["region"], s["observable"],
                                                    "BLINDED (data_obs = sum of MC)" if s["blinded"] else "UNBLINDED"),
             "imax 1   (number of bins/channels)",
             "jmax %d   (number of backgrounds)" % (len(procs) - 1),
             "kmax *   (number of nuisance parameters)",
             sep,
             "%-40s%s     %s     $PROCESS     $PROCESS_$SYSTEMATIC" % ("shapes *", bin_name, root_file),
             "%-40s%s     %s     data_obs" % ("shapes data_obs", bin_name, root_file),
             sep,
             "%-40s%s" % ("bin", bin_name),
             "%-40s%s" % ("observation", -1),
             sep]
    lines.append("%-40s" % "bin" + "".join("%-*s" % (COL_W, bin_name) for _ in procs))
    lines.append("%-40s" % "process" + "".join("%-*s" % (COL_W, p) for p in procs))
    lines.append("%-40s" % "process" + "".join("%-*s" % (COL_W, i) for i in range(len(procs))))
    lines.append("%-40s" % "rate" + "".join("%-*s" % (COL_W, -1) for _ in procs))
    lines.append(sep)

    for nuis, values in LNN.items():
        lines.append("%-*s%-*s" % (NAME_W, nuis, TYPE_W, "lnN") + "".join("%-*s" % (COL_W, values.get(p, "-")) for p in procs))
    for nuis in shape_nuisances:
        lines.append("%-*s%-*s" % (NAME_W, nuis, TYPE_W, "shape") + "".join("%-*s" % (COL_W, 1 if nuis in s["shapes"][p] else "-") for p in procs))

    lines.append("")
    for nuis in FLAT_PARAMS:
        if nuis in shape_nuisances:
            lines.append("%-*s%s" % (NAME_W, nuis, "flatParam"))
    lines.append("%-*s%s" % (NAME_W, bin_name, "autoMCStats %d" % AUTO_MC_STATS))

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("Datacard written to: %s  (%d lnN, %d shape nuisances)" % (out_path, len(LNN), len(shape_nuisances)))


if __name__ == "__main__":
    main()
