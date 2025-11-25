# How to run:
# python plot_jec_variations.py -i dY_UL17_ele_0_500_SR.root -o jec_variations

import ROOT
import os
import argparse
import math
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

def plot_jec_source_comparison(input_file_path, output_dir, jec_source, projection=1, rebin=1):
    # Create comparison plots for a specific JEC source
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    input_file = ROOT.TFile.Open(input_file_path, "READ")
    if not input_file or input_file.IsZombie():
        print("Error: Could not open input file {}".format(input_file_path))
        return False
    
    # Get histograms
    hist_nominal_name = "TTbar_{}".format(projection)
    hist_up_name = "TTbar_{}_jec{}Up".format(projection, jec_source)
    hist_down_name = "TTbar_{}_jec{}Down".format(projection, jec_source)
    
    hist_nominal = input_file.Get(hist_nominal_name)
    hist_up = input_file.Get(hist_up_name)
    hist_down = input_file.Get(hist_down_name)
    
    if not hist_nominal or not hist_up or not hist_down:
        print("Error: Could not find one or more histograms for JEC source {}".format(jec_source))
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
    hist_nominal.SetMarkerSize(2)
    
    hist_up.SetLineColor(ROOT.kGreen+2)  # Dark green
    hist_up.SetLineWidth(2)
    # hist_up.SetLineStyle(2)  # Dashed line
    
    hist_down.SetLineColor(ROOT.kRed)
    hist_down.SetLineWidth(2)
    # hist_down.SetLineStyle(3)  # Dotted line
    
    # Create canvas and legend
    canvas, legend = setup_canvas_and_legend()
    
    # Create ratio pad (top part of canvas)
    ratio_pad = setup_ratio_pad()
    ratio_pad.Draw()
    ratio_pad.cd()
    
    
    # Configure main histogram
    hist_nominal.GetXaxis().SetLabelSize(0)
    hist_nominal.GetYaxis().SetTitle("Events")
    hist_nominal.GetYaxis().SetTitleSize(0.05)
    hist_nominal.GetYaxis().SetTitleOffset(1.2)
    
    # Draw histograms in top pad
    hist_nominal.Draw("HIST")
    hist_up.Draw("HIST SAME")
    hist_down.Draw("HIST SAME")
    
    # Find max y value to adjust the y-axis range
    max_val = max(hist_nominal.GetMaximum(), hist_up.GetMaximum(), hist_down.GetMaximum())
    hist_nominal.GetYaxis().SetRangeUser(0, max_val * 1.3)
    
    # Add legend entries
    legend.AddEntry(hist_nominal, "Nominal", "l")
    legend.AddEntry(hist_up, "{} Up".format(jec_source), "l")
    legend.AddEntry(hist_down, "{} Down".format(jec_source), "l")
    legend.Draw()
    
    # Add CMS text
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextFont(42)
    cms_text.SetTextSize(0.045)
    cms_text.DrawLatex(0.18, 0.93, "#bf{CMS} #it{Simulation Preliminary}")
    
    # Add energy text
    energy_text = ROOT.TLatex()
    energy_text.SetNDC()
    energy_text.SetTextFont(42)
    energy_text.SetTextSize(0.045)
    energy_text.SetTextAlign(31)
    # energy_text.DrawLatex(0.95, 0.93, "13 TeV")
    
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
    output_filename = "{}/jec_{}_TTbar_{}.png".format(output_dir, jec_source, projection)
    canvas.SaveAs(output_filename)
    print("Created plot: {}".format(output_filename))
    
    # Clean up
    canvas.Clear()
    input_file.Close()
    
    return True

def process_all_jec_sources(input_file_path, output_dir, projections, rebin=1):
    # Process all JEC sources found in the input file
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Open input file to get list of available JEC sources
    input_file = ROOT.TFile.Open(input_file_path, "READ")
    if not input_file or input_file.IsZombie():
        print("Error: Could not open input file {}".format(input_file_path))
        return False
    
    # Get all keys in the file
    keys = [key.GetName() for key in input_file.GetListOfKeys()]
    
    # Debug: Print all keys that start with TTbar_1_jec
    for key in keys:
        if "jec" in key and "TTbar_1" in key:
            print("Found histogram: {}".format(key))
    
    # Extract JEC sources from the key names
    jec_sources = set()
    for key in keys:
        if key.startswith("TTbar_1_jec") and key.endswith("Up"):
            # More careful extraction that only removes the exact suffix "Up"
            # This ensures we don't remove "Up" from within the source name
            if key.endswith("Up"):
                source = key[len("TTbar_1_jec"):-2]  # Remove prefix and the last 2 chars ("Up")
            else:
                source = key.replace("TTbar_1_jec", "")
            
            # Verify both Up and Down variations exist before adding to sources
            hist_up_name = "TTbar_1_jec{}Up".format(source)
            hist_down_name = "TTbar_1_jec{}Down".format(source)
            
            if hist_up_name in keys and hist_down_name in keys:
                jec_sources.add(source)
                print("Confirmed valid JEC source: {} (found both {} and {})".format(
                    source, hist_up_name, hist_down_name))
            else:
                print("Skipping incomplete JEC source: {} (missing {} or {})".format(
                    source, hist_up_name, hist_down_name))
            
    input_file.Close()
    
    print("Found %d JEC sources to process" % len(jec_sources))
    
    # Create plots for each JEC source and projection
    for projection in projections:
        for jec_source in sorted(jec_sources):
            print("Processing JEC source: {}, Projection: {}".format(jec_source, projection))
            plot_jec_source_comparison(input_file_path, output_dir, jec_source, projection, rebin)
    
    return True

def plot_overall_jec_comparison(input_file_path, output_dir, projection=1, rebin=1):
    # Create comparison plots for the overall JEC systematic variation
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    input_file = ROOT.TFile.Open(input_file_path, "READ")
    if not input_file or input_file.IsZombie():
        print("Error: Could not open input file {}".format(input_file_path))
        return False
    
    # Get histograms
    hist_nominal_name = "TTbar_{}".format(projection)
    hist_up_name = "TTbar_{}_jecUp".format(projection)
    hist_down_name = "TTbar_{}_jecDown".format(projection)
    
    hist_nominal = input_file.Get(hist_nominal_name)
    hist_up = input_file.Get(hist_up_name)
    hist_down = input_file.Get(hist_down_name)
    
    if not hist_nominal or not hist_up or not hist_down:
        print("Error: Could not find overall JEC histograms")
        input_file.Close()
        return False
    
    # Rest of the function is identical to plot_jec_source_comparison
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
    hist_nominal.SetLineWidth(2)
    hist_nominal.SetMarkerStyle(20)
    hist_nominal.SetMarkerSize(1.2)
    
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
    
    
    # Configure main histogram
    hist_nominal.GetXaxis().SetLabelSize(0)
    hist_nominal.GetYaxis().SetTitle("Events")
    hist_nominal.GetYaxis().SetTitleSize(0.05)
    hist_nominal.GetYaxis().SetTitleOffset(1.2)
    
    # Draw histograms in top pad
    hist_nominal.Draw("HIST")
    hist_up.Draw("HIST SAME")
    hist_down.Draw("HIST SAME")
    
    # Find max y value to adjust the y-axis range
    max_val = max(hist_nominal.GetMaximum(), hist_up.GetMaximum(), hist_down.GetMaximum())
    hist_nominal.GetYaxis().SetRangeUser(0, max_val * 1.3)
    
    # Add legend entries
    legend.AddEntry(hist_nominal, "Nominal", "l")
    legend.AddEntry(hist_up, "JEC Total Up", "l")
    legend.AddEntry(hist_down, "JEC Total Down", "l")
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
    output_filename = "{}/jec_Total_TTbar_{}.png".format(output_dir, projection)
    canvas.SaveAs(output_filename)
    print("Created plot: {}".format(output_filename))
    
    # Clean up
    canvas.Clear()
    input_file.Close()
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Create JEC source comparison plots')
    parser.add_argument('-i', '--input', required=True, help='Input ROOT file containing TTbar distributions')
    parser.add_argument('-o', '--output', default='./jec_plots', help='Output directory for plots')
    parser.add_argument('-p', '--projections', default='1,2', help='Comma-separated list of projections to process')
    parser.add_argument('-r', '--rebin', type=int, default=1, help='Rebinning factor')
    parser.add_argument('-s', '--sources', help='Comma-separated list of JEC sources to process (default: all)')
    parser.add_argument('-t', '--total', action='store_true', help='Create plot for total JEC variation')
    args = parser.parse_args()
    
    # Convert projections string to list
    projections = [int(p) for p in args.projections.split(',')]
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # Process total JEC variation if requested
    if args.total:
        for projection in projections:
            plot_overall_jec_comparison(args.input, args.output, projection, args.rebin)
    
    # Process specific sources if provided, otherwise process all
    if args.sources:
        sources = args.sources.split(',')
        for projection in projections:
            for source in sources:
                plot_jec_source_comparison(args.input, args.output, source, projection, args.rebin)
    else:
        process_all_jec_sources(args.input, args.output, projections, args.rebin)
    
    print("All plots saved to {}".format(args.output))

if __name__ == "__main__":
    # Enable batch mode for ROOT
    ROOT.gROOT.SetBatch(True)
    # Improve text rendering
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)
    
    main()