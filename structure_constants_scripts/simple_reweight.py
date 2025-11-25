#!/usr/bin/env python
"""
Simple script to read structure constants from ROOT files and calculate new weights.

Usage:
  python simple_reweight.py <input_root_file> [options]

Example:
  python simple_reweight.py output.root --variable Mass_tt --output-plot mass_distribution.png
"""

import ROOT
import numpy as np
import argparse
import os
import sys
import matplotlib.pyplot as plt

# Wilson coefficient names for reference
WC_NAMES = {
    0: "ctGRe",
    1: "ctGIm",
    2: "cQj18",
    3: "cQj38",
    4: "cQj11",
    5: "cQj31",
    6: "ctu8",
    7: "ctd8",
    8: "ctj8",
    9: "cQu8",
    10: "cQd8",
    11: "ctu1",
    12: "ctd1",
    13: "ctj1",
    14: "cQu1",
    15: "cQd1"
}

def calculate_weight(structure_constants, wc_values):
    """Calculate new weight using structure constants and Wilson coefficient values."""
    num_WCs = len(wc_values)
    
    # Constant term (SM)
    c0 = structure_constants[0]
    
    # Linear part
    s_linear = structure_constants[1:1+num_WCs]
    w_linear = np.dot(s_linear, wc_values)
    
    # Quadratic part
    quad_list = []
    # Diagonal terms
    for i in range(num_WCs):
        quad_list.append(wc_values[i]**2)
    # Cross terms
    for i in range(num_WCs):
        for j in range(i+1, num_WCs):
            quad_list.append(wc_values[i]*wc_values[j])
    
    idx_quad_start = 1 + num_WCs
    s_quad = structure_constants[idx_quad_start:]
    w_quad = np.dot(s_quad, quad_list)
    
    # Total weight
    total_w = c0 + w_linear + w_quad
    
    return total_w

def main():
    parser = argparse.ArgumentParser(description="Read structure constants and calculate new weights")
    parser.add_argument("input_file", help="Input ROOT file with structure constants")
    parser.add_argument("--tree-name", default="AnalysisTree", help="Name of the TTree (default: AnalysisTree)")
    parser.add_argument("--branch-name", default="structure_constants", help="Name of the branch containing structure constants")
    parser.add_argument("--variable", default="Mass_tt", help="Variable to plot (default: Mass_tt)")
    parser.add_argument("--output-file", default="reweighted.root", help="Output ROOT file")
    parser.add_argument("--output-plot", default="reweighted_plot.png", help="Output plot file")
    parser.add_argument("--max-events", type=int, help="Maximum number of events to process")
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input_file):
        print("Error: Input file does not exist: %s" % args.input_file)
        sys.exit(1)
    
    # Open ROOT file
    root_file = ROOT.TFile.Open(args.input_file)
    if not root_file or root_file.IsZombie():
        print("Error: Could not open ROOT file: %s" % args.input_file)
        sys.exit(1)
    
    # Get TTree
    tree = root_file.Get(args.tree_name)
    if not tree:
        print("Error: Could not find TTree: %s" % args.tree_name)
        root_file.Close()
        sys.exit(1)
    
    # Check if branches exist
    branch_names = [b.GetName() for b in tree.GetListOfBranches()]
    if args.branch_name not in branch_names:
        print("Error: Branch '%s' not found in TTree. Available branches:" % args.branch_name)
        for name in branch_names:
            print("  - %s" % name)
        root_file.Close()
        sys.exit(1)
    
    if args.variable not in branch_names:
        print("Error: Variable '%s' not found in TTree. Available branches:" % args.variable)
        for name in branch_names:
            print("  - %s" % name)
        root_file.Close()
        sys.exit(1)
    
    # ===== DEFINE YOUR WILSON COEFFICIENT VALUES HERE =====
    # Format: [ctGRe, ctGIm, cQj18, cQj38, cQj11, cQj31, ctu8, ctd8, ctj8, cQu8, cQd8, ctu1, ctd1, ctj1, cQu1, cQd1]
    
    # Define scenarios to test - MODIFY THIS SECTION AS NEEDED
    scenarios = [
        {"name": "SM", "values": [0.0] * 16},
        {"name": "ctGRe=1.0", "values": [1.0] + [0.0] * 15},
        {"name": "ctGIm=1.0", "values": [0.0, 1.0] + [0.0] * 14},
        {"name": "ctGRe=0.5,ctGIm=0.5", "values": [0.5, 0.5] + [0.0] * 14},
        
        # Add your custom scenarios here, for example:
        # {"name": "My Custom Point", "values": [0.5, 0.3, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]},
    ]
    # ======================================================
    
    # Create output ROOT file
    output_file = ROOT.TFile(args.output_file, "RECREATE")
    
    # Create histograms for each scenario
    histograms = {}
    for scenario in scenarios:
        name = scenario["name"]
        histograms[name] = ROOT.TH1F("h_" + name.replace("=", "_").replace(",", "_"), 
                                     name + ";" + args.variable + ";Events", 50, 0, 2000)
    
    # Get number of entries
    n_entries = tree.GetEntries()
    if args.max_events:
        n_entries = min(args.max_events, n_entries)
    
    print("Processing %d events..." % n_entries)
    
    # Process events
    for i in range(n_entries):
        if i % 1000 == 0:
            print("  Event %d / %d" % (i, n_entries))
        
        tree.GetEntry(i)
        
        # Get structure constants
        structure_constants = getattr(tree, args.branch_name)
        
        # Convert to numpy array if needed
        if isinstance(structure_constants, list):
            structure_constants = np.array(structure_constants)
        
        # Get variable value
        variable_value = getattr(tree, args.variable)
        
        # Calculate weights for each scenario and fill histograms
        for scenario in scenarios:
            name = scenario["name"]
            wc_values = scenario["values"]
            weight = calculate_weight(structure_constants, wc_values)
            histograms[name].Fill(variable_value, weight)
    
    # Write histograms to output file
    output_file.cd()
    for name, hist in histograms.items():
        hist.Write()
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Plot histograms
    for scenario in scenarios:
        name = scenario["name"]
        hist = histograms[name]
        
        # Convert ROOT histogram to numpy arrays
        x = np.array([hist.GetBinCenter(i) for i in range(1, hist.GetNbinsX() + 1)])
        y = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])
        
        # Normalize histogram
        if np.sum(y) > 0:
            y = y / np.sum(y)
        
        plt.plot(x, y, label=name)
    
    plt.xlabel(args.variable)
    plt.ylabel("Normalized Events")
    plt.title("Effect of Different Wilson Coefficient Values")
    plt.legend()
    plt.grid(True)
    plt.savefig(args.output_plot)
    
    # Close files
    output_file.Close()
    root_file.Close()
    
    print("Done! Output saved to %s and %s" % (args.output_file, args.output_plot))

if __name__ == "__main__":
    main() 