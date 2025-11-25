# How to run:
# python plot_syst_variations.py -i dY_UL17_ele_0_500_SR.root -o syst_variations

import ROOT
import os
import argparse
import math
import re
from array import array

def setup_canvas_and_legend():
    # Set up a canvas and legend with appropriate styling
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 700)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.10)
    canvas.SetBottomMargin(0.15)
    
    legend = ROOT.TLegend(0.65, 0.75, 0.92, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    legend.SetTextFont(42)
    
    return canvas, legend

def setup_ratio_pad():
    # Create a ratio pad for the bottom of the canvas
    ratio_pad = ROOT.TPad("ratio_pad", "ratio_pad", 0, 0.3, 1, 1)
    ratio_pad.SetLeftMargin(0.15)
    ratio_pad.SetRightMargin(0.05)
    ratio_pad.SetTopMargin(0.10)
    ratio_pad.SetBottomMargin(0.02)
    ratio_pad.SetFillColor(0)
    ratio_pad.SetTickx(1)
    ratio_pad.SetTicky(1)
    return ratio_pad

def setup_main_pad():
    # Create the main histogram pad
    main_pad = ROOT.TPad("main_pad", "main_pad", 0, 0, 1, 0.3)
    main_pad.SetLeftMargin(0.15)
    main_pad.SetRightMargin(0.05)
    main_pad.SetTopMargin(0.05)
    main_pad.SetBottomMargin(0.35)
    main_pad.SetFillColor(0)
    main_pad.SetTickx(1)
    main_pad.SetTicky(1)
    return main_pad

def create_ratio_histogram(h_var, h_nominal, ratio_title="Variation / Nominal"):
    # Create a ratio histogram comparing a variation to the nominal
    h_ratio = h_var.Clone(h_var.GetName() + "_ratio")
    h_ratio.SetTitle("")
    h_ratio.Divide(h_nominal)
    
    # Configure ratio histogram style
    h_ratio.GetYaxis().SetTitle(ratio_title)
    h_ratio.GetYaxis().SetRangeUser(0.95, 1.05)  # For 5% variations
    h_ratio.GetYaxis().SetTitleSize(0.12)
    h_ratio.GetYaxis().SetTitleOffset(0.5)
    h_ratio.GetYaxis().SetLabelSize(0.10)
    h_ratio.GetYaxis().CenterTitle()
    
    # Reduce number of divisions on y-axis (fewer tick marks)
    h_ratio.GetYaxis().SetNdivisions(505)  # 5 primary divisions, 5 secondary divisions
    
    h_ratio.GetXaxis().SetTitleSize(0.12)
    h_ratio.GetXaxis().SetTitleOffset(1.0)
    h_ratio.GetXaxis().SetLabelSize(0.10)
    h_ratio.GetXaxis().SetTickLength(0.07)
    
    return h_ratio

def find_systematic_variations(input_file_path):
    # Find all systematic variations in the ROOT file
    input_file = ROOT.TFile.Open(input_file_path, "READ")
    if not input_file or input_file.IsZombie():
        print("Error: Could not open input file {}".format(input_file_path))
        return {}
    
    # Get all keys in the file
    keys = [key.GetName() for key in input_file.GetListOfKeys()]
    
    # Find systematic variations
    syst_variations = {}
    pattern = r'TTbar_(\d+)_(\w+)(Up|Down)'
    
    for key in keys:
        match = re.match(pattern, key)
        if match:
            projection = match.group(1)
            syst_name = match.group(2)
            direction = match.group(3)
            
            if syst_name not in syst_variations:
                syst_variations[syst_name] = set()
            
            syst_variations[syst_name].add(projection)
    
    input_file.Close()
    
    # Convert sets to lists
    result = {}
    for syst, projections in syst_variations.items():
        result[syst] = sorted(list(projections))
    
    return result

def plot_systematic_comparison(input_file_path, output_dir, syst_name, projection, rebin=1, use_log=False, y_max=None):
    # Create comparison plots for a specific systematic
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    input_file = ROOT.TFile.Open(input_file_path, "READ")
    if not input_file or input_file.IsZombie():
        print("Error: Could not open input file {}".format(input_file_path))
        return False
    
    # Format systematic name for display
    display_name = format_syst_name(syst_name)
    
    # Get histograms
    hist_nominal_name = "TTbar_{}".format(projection)
    hist_up_name = "TTbar_{}_{}Up".format(projection, syst_name)
    hist_down_name = "TTbar_{}_{}Down".format(projection, syst_name)
    
    hist_nominal = input_file.Get(hist_nominal_name)
    hist_up = input_file.Get(hist_up_name)
    hist_down = input_file.Get(hist_down_name)
    
    if not hist_nominal or not hist_up or not hist_down:
        print("Error: Could not find one or more histograms for systematic {}".format(syst_name))
        input_file.Close()
        return False
    
    # Clone histograms to avoid memory issues
    hist_nominal = hist_nominal.Clone(hist_nominal_name + "_clone")
    hist_up = hist_up.Clone(hist_up_name + "_clone")
    hist_down = hist_down.Clone(hist_down_name + "_clone")
    
    # Apply rebinning if requested
    if rebin > 1:
        hist_nominal.Rebin(rebin)
        hist_up.Rebin(rebin)
        hist_down.Rebin(rebin)
    
    # Set histogram styles
    hist_nominal.SetLineColor(ROOT.kBlack)
    hist_nominal.SetLineWidth(3)
    hist_nominal.SetMarkerStyle(20)
    hist_nominal.SetMarkerSize(1.5)
    
    hist_up.SetLineColor(ROOT.kGreen+2)  # Dark green
    hist_up.SetLineWidth(2)
    
    hist_down.SetLineColor(ROOT.kRed)
    hist_down.SetLineWidth(2)
    
    # Create canvas and legend
    canvas, legend = setup_canvas_and_legend()
    
    # Create ratio pad (top part of canvas)
    ratio_pad = setup_ratio_pad()
    ratio_pad.Draw()
    ratio_pad.cd()
    
    # Set logarithmic scale for y-axis if requested
    if use_log:
        ratio_pad.SetLogy(True)
    
    # Configure main histogram
    hist_nominal.GetXaxis().SetLabelSize(0)
    hist_nominal.GetYaxis().SetTitle("Events")
    hist_nominal.GetYaxis().SetTitleSize(0.05)
    hist_nominal.GetYaxis().SetTitleOffset(1.2)
    
    # Draw histograms in top pad
    hist_nominal.Draw("HIST")
    hist_up.Draw("HIST SAME")
    hist_down.Draw("HIST SAME")
    
    # Set y-axis range
    if y_max is not None:
        hist_nominal.GetYaxis().SetRangeUser(0, y_max)
    else:
        # Find max y value to adjust the y-axis range
        max_val = max(hist_nominal.GetMaximum(), hist_up.GetMaximum(), hist_down.GetMaximum())
        
        # Adjust y-axis range based on whether we're using log scale
        if use_log:
            # Find minimum non-zero value for log scale
            min_vals = []
            for bin_idx in range(1, hist_nominal.GetNbinsX() + 1):
                for hist in [hist_nominal, hist_up, hist_down]:
                    val = hist.GetBinContent(bin_idx)
                    if val > 0:
                        min_vals.append(val)
            
            if min_vals:
                min_val = min(min_vals) * 0.5  # Go a bit lower than the minimum
            else:
                min_val = 0.1  # Default if no positive values found
            
            hist_nominal.GetYaxis().SetRangeUser(min_val, max_val * 10.0)
        else:
            hist_nominal.GetYaxis().SetRangeUser(0, max_val * 1.3)
    
    # Add legend entries with proper formatting
    legend.AddEntry(hist_nominal, "Nominal", "l")
    legend.AddEntry(hist_up, "{} Up".format(display_name), "l")
    legend.AddEntry(hist_down, "{} Down".format(display_name), "l")
    legend.Draw()
    
    # Add CMS text
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextFont(42)
    cms_text.SetTextSize(0.045)
    cms_text.DrawLatex(0.18, 0.93, "#bf{CMS} #it{Simulation Preliminary}")
    
    # Switch to lower pad for the ratio plot
    canvas.cd()
    main_pad = setup_main_pad()
    main_pad.Draw()
    main_pad.cd()
    
    # Create and draw ratio histograms
    ratio_up = create_ratio_histogram(hist_up, hist_nominal, "Variation/Nominal")
    ratio_down = create_ratio_histogram(hist_down, hist_nominal)
    
    ratio_up.Draw("HIST")
    ratio_down.Draw("HIST SAME")
    
    # Draw a reference line at 1.0 for the ratio
    ref_line = ROOT.TLine(
        ratio_up.GetXaxis().GetXmin(), 1.0,
        ratio_up.GetXaxis().GetXmax(), 1.0
    )
    ref_line.SetLineStyle(2)
    ref_line.SetLineColor(ROOT.kGray+2)
    ref_line.Draw("SAME")
    
    # Save the canvas
    output_filename = "{}/syst_{}_{}.png".format(output_dir, syst_name, projection)
    canvas.SaveAs(output_filename)
    print("Created plot: {}".format(output_filename))
    
    # Clean up
    canvas.Clear()
    input_file.Close()
    
    return True

def format_syst_name(syst_name):
    # Format systematic name for display in plot legend
    # Handle special cases
    if syst_name == "murmuf":
        return "murmuf Scale"
    elif syst_name == "pdf":
        return "PDF"
    elif syst_name == "jecUp" or syst_name == "jecDown":
        return "JEC Total"
    
    # For other names, convert camelCase to words with spaces
    # First, add space before capital letters
    formatted = re.sub(r'([a-z])([A-Z])', r'\1 \2', syst_name)
    # Capitalize first letter
    formatted = formatted[0].upper() + formatted[1:]
    
    # Special handling for common abbreviations
    formatted = formatted.replace("Mu ", "Muon ")
    formatted = formatted.replace("Ele ", "Electron ")
    formatted = formatted.replace("Iso ", "Isolation ")
    formatted = formatted.replace("Id ", "ID ")
    formatted = formatted.replace("Pu ", "PU ")
    formatted = formatted.replace("Btag ", "b-tag ")
    formatted = formatted.replace("Ttag ", "t-tag ")
    formatted = formatted.replace("Hf", "HF")
    formatted = formatted.replace("Lf", "LF")
    formatted = formatted.replace("Cf", "CF")
    formatted = formatted.replace("Stat", "Stat.")
    formatted = formatted.replace("Syst", "Syst.")
    formatted = formatted.replace("Isr", "ISR")
    formatted = formatted.replace("Fsr", "FSR")
    
    return formatted

def process_all_systematics(input_file_path, output_dir, projections=None, rebin=1, use_log=False, y_max=None, excluded_systs=None):
    # Process all systematic variations found in the input file
    # Find all systematic variations
    syst_variations = find_systematic_variations(input_file_path)
    
    if not syst_variations:
        print("No systematic variations found in the input file")
        return False
    
    print("Found {} systematic variations".format(len(syst_variations)))
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Process each systematic variation
    for syst_name, available_projections in sorted(syst_variations.items()):
        # Skip excluded systematics
        if excluded_systs and syst_name in excluded_systs:
            print("Skipping excluded systematic: {}".format(syst_name))
            continue
            
        # Filter projections if specified
        if projections:
            process_projections = [p for p in projections if p in available_projections]
        else:
            process_projections = available_projections
        
        if not process_projections:
            print("No matching projections found for systematic: {}".format(syst_name))
            continue
        
        for projection in process_projections:
            print("Processing systematic: {}, Projection: {}".format(syst_name, projection))
            plot_systematic_comparison(
                input_file_path, 
                output_dir, 
                syst_name, 
                projection, 
                rebin, 
                use_log,
                y_max
            )
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Create systematic variation comparison plots')
    parser.add_argument('-i', '--input', required=True, help='Input ROOT file containing TTbar distributions')
    parser.add_argument('-o', '--output', default='./syst_plots', help='Output directory for plots')
    parser.add_argument('-p', '--projections', default='1,2', help='Comma-separated list of projections to process')
    parser.add_argument('-r', '--rebin', type=int, default=1, help='Rebinning factor')
    parser.add_argument('-s', '--systematics', help='Comma-separated list of systematics to process (default: all)')
    parser.add_argument('-e', '--exclude', help='Comma-separated list of systematics to exclude')
    parser.add_argument('-l', '--log', action='store_true', help='Use logarithmic scale for y-axis')
    parser.add_argument('-y', '--ymax', type=float, help='Maximum value for y-axis')
    args = parser.parse_args()
    
    # Convert projections string to list
    projections = [p for p in args.projections.split(',')]
    
    # Convert exclude list to set for faster lookups
    excluded_systs = set(args.exclude.split(',')) if args.exclude else None
    
    # Process specific systematics if provided, otherwise process all
    if args.systematics:
        systematics = args.systematics.split(',')
        for syst in systematics:
            for projection in projections:
                plot_systematic_comparison(
                    args.input, 
                    args.output, 
                    syst, 
                    projection, 
                    args.rebin,
                    args.log,
                    args.ymax
                )
    else:
        process_all_systematics(
            args.input, 
            args.output, 
            projections, 
            args.rebin,
            args.log,
            args.ymax,
            excluded_systs
        )
    
    print("All plots saved to {}".format(args.output))

if __name__ == "__main__":
    # Enable batch mode for ROOT
    ROOT.gROOT.SetBatch(True)
    # Improve text rendering
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)
    
    main()