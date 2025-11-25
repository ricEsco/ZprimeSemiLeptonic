#!/usr/bin/env python
"""
Script to create reweighted histograms using structure constants

This script:
1. Opens a ROOT file containing structure constants
2. Reads WC scenarios from command line or JSON file
3. Creates reweighted histograms for specified variables
4. Saves the histograms to an output ROOT file

Usage:
  python reweight_histograms.py <input_root_file> --output <output_root_file> [--wc-config <json_file>] [--variables <var1,var2,...>]
python reweight_histograms.py EFT_700_900.root --output reweighted_histograms.root --wc-config scenarios.json
python reweight_histograms.py EFT_700_900.root --output reweighted_histograms.root --variables dyreco,Sigma_phi

"""

import ROOT
import numpy as np
import sys
import argparse
import json
import os

def get_wc_names():
    return {
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

def create_wc_scenario(config):
    # Create a WC values array from a configuration dictionary 
    wc_values = [0.0] * 16
    wc_names = get_wc_names()
    name_to_idx = {name: idx for idx, name in wc_names.items()}
    
    for wc_name, value in config.items():
        if wc_name in name_to_idx:
            wc_values[name_to_idx[wc_name]] = value
    
    return wc_values

def calculate_weight(structure_constants, wc_values):
    """Calculate weight using structure constants and Wilson coefficient values"""
    num_WCs = len(wc_values)
    
    # Constant term (SM)
    c0 = structure_constants[0]
    
    # Linear terms
    linear_terms = structure_constants[1:1+num_WCs]
    w_linear = np.dot(linear_terms, wc_values)
    
    # Quadratic terms
    quad_list = []
    # Diagonal terms (WC^2)
    for i in range(num_WCs):
        quad_list.append(wc_values[i]**2)
    
    # Cross terms (WC_i * WC_j)
    for i in range(num_WCs):
        for j in range(i+1, num_WCs):
            quad_list.append(wc_values[i] * wc_values[j])
    
    idx_quad_start = 1 + num_WCs
    quad_terms = structure_constants[idx_quad_start:]
    w_quad = np.dot(quad_terms, quad_list)
    
    # Total weight
    total_weight = c0 + w_linear + w_quad
    
    return total_weight

def load_wc_scenarios(config_file):
    # Load WC scenarios from a JSON file
    with open(config_file, 'r') as f:
        scenarios = json.load(f)
    
    # Convert the scenarios to the format needed for calculation
    wc_scenarios = {}
    for name, config in scenarios.items():
        wc_scenarios[name] = create_wc_scenario(config)
    
    return wc_scenarios

def create_histogram(name, title, nbins, xmin, xmax):
    # Create a ROOT histogram with the given parameters
    return ROOT.TH1F(name, title, nbins, xmin, xmax)

def main():
    parser = argparse.ArgumentParser(description="Create reweighted histograms using structure constants")
    parser.add_argument("input_file", help="Input ROOT file containing structure constants")
    parser.add_argument("--output", default="reweighted_histograms.root", help="Output ROOT file for histograms")
    parser.add_argument("--wc-config", help="JSON file containing Wilson coefficient scenarios")
    parser.add_argument("--variables", default="dyreco,dyreco_1,dyreco_2,Sigma_phi,Sigma_phi_1,Sigma_phi_2,Delta_phi", 
                        help="Comma-separated list of variables to create histograms for")
    parser.add_argument("--regions", default="SR,CR1,CR2", help="Comma-separated list of regions")
    parser.add_argument("--mass-bins", action="store_true", help="Create histograms for different mass bins")
    args = parser.parse_args()
    
    input_file = ROOT.TFile.Open(args.input_file)
    if not input_file or input_file.IsZombie():
        print("Error: Could not open input file: {}".format(args.input_file))
        return 1
    
    # Get the tree
    tree = input_file.Get("AnalysisTree")
    if not tree:
        print("Error: Could not find TTree: AnalysisTree")
        input_file.Close()
        return 1
    
    # Check if structure_constants branch exists
    branch_names = [b.GetName() for b in tree.GetListOfBranches()]
    if "structure_constants" not in branch_names:
        print("Error: structure_constants branch not found in TTree")
        print("Available branches:", branch_names)
        input_file.Close()
        return 1
    
    # Define WC scenarios
    wc_scenarios = {
        "SM": [0] * 16,
        "ctGRe=1": [1] + [0] * 15,
        "ctGIm=1": [0, 1] + [0] * 14,
        "ctGRe=1,ctGIm=1": [1, 1] + [0] * 14
    }
    
    # Load scenarios from JSON file if provided
    if args.wc_config:
        try:
            wc_scenarios.update(load_wc_scenarios(args.wc_config))
        except Exception as e:
            print("Error loading WC scenarios from {}: {}".format(args.wc_config, str(e)))
    
    # Parse variables and regions
    variables = args.variables.split(',')
    regions = args.regions.split(',')
    
    # Create output ROOT file
    output_file = ROOT.TFile(args.output, "RECREATE")
    
    # Create histograms for each variable, region, and scenario
    histograms = {}
    
    # Define histogram parameters for each variable
    histogram_params = {
        "dyreco": (100, -5, 5),
        "dyreco_1": (100, -5, 5),
        "dyreco_2": (100, -5, 5),
        "Sigma_phi": (100, -3.15, 3.15),
        "Sigma_phi_1": (100, -3.15, 3.15),
        "Sigma_phi_2": (100, -3.15, 3.15),
        "Delta_phi": (100, -3.15, 3.15)
    }
    
    # Create histograms
    for region in regions:
        for var in variables:
            var_name = var if region == "SR" else "{}_{}".format(var, region)
            if var_name not in branch_names:
                print("Warning: Variable {} not found in tree, skipping".format(var_name))
                continue
            
            # Get histogram parameters
            if var in histogram_params:
                nbins, xmin, xmax = histogram_params[var]
            else:
                nbins, xmin, xmax = 100, -5, 5  # Default
            
            # Create histograms for each scenario
            for scenario_name in wc_scenarios:
                hist_name = "{}_{}".format(var_name, scenario_name.replace("=", "_").replace(",", "_"))
                hist_title = "{} ({})".format(var_name, scenario_name)
                histograms[hist_name] = create_histogram(hist_name, hist_title, nbins, xmin, xmax)
                
                # Create mass-binned histograms if requested
                if args.mass_bins:
                    mass_bins = ["0_500", "500_750", "750_1000", "1000_1500", "1500_Inf"]
                    for mass_bin in mass_bins:
                        mass_var_name = "{}_{}".format(var_name, mass_bin)
                        if mass_var_name in branch_names:
                            mass_hist_name = "{}_{}_{}".format(var_name, mass_bin, scenario_name.replace("=", "_").replace(",", "_"))
                            mass_hist_title = "{} {} ({})".format(var_name, mass_bin, scenario_name)
                            histograms[mass_hist_name] = create_histogram(mass_hist_name, mass_hist_title, nbins, xmin, xmax)
    
    # Process events and fill histograms
    nentries = tree.GetEntries()
    print("Processing {} events...".format(nentries))
    
    for i in range(nentries):
        if i % 10000 == 0:
            print("  Event {}/{}".format(i, nentries))
        
        tree.GetEntry(i)
        
        # Get structure constants for this event
        struct_constants = np.array(tree.structure_constants)
        
        # Get original event weight
        event_weight = 1.0
        if hasattr(tree, "weight"):
            event_weight = tree.weight
        
        # Calculate weights for each scenario
        scenario_weights = {}
        for scenario_name, wc_values in wc_scenarios.items():
            scenario_weights[scenario_name] = calculate_weight(struct_constants, wc_values)
        
        # Fill histograms for each variable and scenario
        for region in regions:
            for var in variables:
                var_name = var if region == "SR" else "{}_{}".format(var, region)
                if var_name not in branch_names:
                    continue
                
                # Get variable value
                var_value = getattr(tree, var_name)
                
                # Fill histograms for each scenario
                for scenario_name in wc_scenarios:
                    hist_name = "{}_{}".format(var_name, scenario_name.replace("=", "_").replace(",", "_"))
                    if hist_name in histograms:
                        histograms[hist_name].Fill(var_value, event_weight * scenario_weights[scenario_name])
                
                # Fill mass-binned histograms if requested
                if args.mass_bins:
                    mass_bins = ["0_500", "500_750", "750_1000", "1000_1500", "1500_Inf"]
                    for mass_bin in mass_bins:
                        mass_var_name = "{}_{}".format(var_name, mass_bin)
                        if mass_var_name in branch_names:
                            mass_var_value = getattr(tree, mass_var_name)
                            for scenario_name in wc_scenarios:
                                mass_hist_name = "{}_{}_{}".format(var_name, mass_bin, scenario_name.replace("=", "_").replace(",", "_"))
                                if mass_hist_name in histograms:
                                    histograms[mass_hist_name].Fill(mass_var_value, event_weight * scenario_weights[scenario_name])
    
    # Write histograms to output file
    output_file.cd()
    for hist in histograms.values():
        hist.Write()
    
    # Create a directory for scenario information
    info_dir = output_file.mkdir("ScenarioInfo")
    info_dir.cd()
    
    # Save scenario information
    for scenario_name, wc_values in wc_scenarios.items():
        # Create a histogram to store the WC values
        wc_hist = ROOT.TH1F("wc_values_{}".format(scenario_name.replace("=", "_").replace(",", "_")), 
                           "Wilson Coefficient Values for {}".format(scenario_name), 
                           16, 0, 16)
        
        # Set bin labels to WC names
        wc_names = get_wc_names()
        for i in range(16):
            wc_hist.GetXaxis().SetBinLabel(i+1, wc_names[i])
            wc_hist.SetBinContent(i+1, wc_values[i])
        
        wc_hist.Write()
    
    # Close files
    output_file.Close()
    input_file.Close()
    
    print("Done! Reweighted histograms saved to: {}".format(args.output))
    return 0

if __name__ == "__main__":
    sys.exit(main()) 