# Example:
#     python compare_powheg_eft_sm.py --powheg-dir /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_DNN/UL17/muon/nominal_ttbar/ --eft-dir /data/dust/group/cms/zprime-uhh/Analysis_EFT_UL17/muon/workdir_Analysis_EFT_UL17_muon_semilepton_deltayreco_deltaPhi_sigmaPhi_ttree/ --output-dir powheg_eft_comparison

# powheg: /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_DNN/UL17/muon/workdir_AnalysisDNN_2017_muon_ttbar
# eft: /data/dust/group/cms/zprime-uhh/AnalysisDNN_EFT_UL17/muon/workdir_AnalysisDNN_UL17_EFT_muon_semilepton/

import argparse
import glob
import os
import ROOT
import numpy as np
import sys
import time

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compare Powheg nominal SM with EFT SM samples")
    parser.add_argument("--powheg-dir", required=True, help="Directory containing Powheg nominal SM samples")
    parser.add_argument("--eft-dir", required=True, help="Directory containing EFT SM samples")
    parser.add_argument("--output-dir", default="powheg_eft_comparison", help="Output directory for plots")
    parser.add_argument("--delta-y-only", action="store_true", help="Only process Delta Y")
    parser.add_argument("--delta-phi-only", action="store_true", help="Only process Delta Phi")
    parser.add_argument("--sigma-phi-only", action="store_true", help="Only process Sigma Phi")
    parser.add_argument("--max-events", type=int, default=-1, help="Maximum number of events to process per file (-1 for all)")
    parser.add_argument("--max-files", type=int, default=-1, help="Maximum number of files to process from each directory (-1 for all)")
    parser.add_argument("--print-frequency", type=int, default=10000, help="How often to print progress (events)")
    return parser.parse_args()

def find_files(directory, pattern="*.root", max_files=-1):
    # all files matching pattern in directory, with option to limit
    files = glob.glob(os.path.join(directory, pattern))
    if max_files > 0 and len(files) > max_files:
        print("Limiting to {max_files} files out of {len(files)} available".format(max_files=max_files, len=len(files)))
        return files[:max_files]
    return files

def make_directory(directory):
    # Create directory if it doesn't exist.
    if not os.path.exists(directory):
        os.makedirs(directory)

def create_comparison_plot(h_powheg, h_eft_sm, variable_name, output_dir, title="Powheg SM vs EFT SM", custom_range=None):
    # Create normalized histograms
    h_powheg_norm = h_powheg.Clone("{0}_norm".format(h_powheg.GetName()))
    h_eft_sm_norm = h_eft_sm.Clone("{0}_norm".format(h_eft_sm.GetName()))
    
    # Normalize histograms if they have entries
    if h_powheg_norm.Integral() > 0:
        h_powheg_norm.Scale(1.0 / h_powheg_norm.Integral())
    
    if h_eft_sm_norm.Integral() > 0:
        h_eft_sm_norm.Scale(1.0 / h_eft_sm_norm.Integral())
    
    # Create ratio histogram (Powheg/EFT SM)
    h_ratio = h_powheg_norm.Clone("hRatio_{0}".format(variable_name))
    h_ratio.Divide(h_eft_sm_norm)
    
    # Create a single canvas divided into two pads
    c_combined = ROOT.TCanvas("c_combined_{0}".format(variable_name), "{0}: Powheg SM vs EFT SM".format(variable_name), 800, 800)
    
    # Set up top pad for distributions (give 60% of canvas to top instead of 70%)
    pad_top = ROOT.TPad("pad_top", "pad_top", 0, 0.4, 1, 1.0)
    pad_top.SetBottomMargin(0.02)  # Small margin at bottom of top pad
    pad_top.SetLeftMargin(0.15)
    pad_top.SetRightMargin(0.05)
    pad_top.Draw()
    
    # Set up bottom pad for ratio (enlarged the yaxis taking 40% of canvas)
    pad_bottom = ROOT.TPad("pad_bottom", "pad_bottom", 0, 0.0, 1, 0.4)
    pad_bottom.SetTopMargin(0.02)  # Small margin at top of bottom pad
    pad_bottom.SetBottomMargin(0.3)  # Slightly reduced bottom margin
    pad_bottom.SetLeftMargin(0.15)
    pad_bottom.SetRightMargin(0.05)
    pad_bottom.SetGridy()
    pad_bottom.Draw()
    
    # Draw histograms in top pad
    pad_top.cd()
    
    # Set up normalized histograms
    h_powheg_norm.SetTitle("")
    h_powheg_norm.GetYaxis().SetTitle("Normalized Events")
    h_powheg_norm.GetYaxis().SetTitleOffset(1.5)
    h_powheg_norm.GetXaxis().SetLabelSize(0)  # Hide x-axis labels in top pad
    h_powheg_norm.GetXaxis().SetTitleSize(0)  # Hide x-axis title in top pad
    
    # Style for Powheg histogram
    h_powheg_norm.SetLineColor(ROOT.kBlue)
    h_powheg_norm.SetLineWidth(2)
    h_powheg_norm.SetMarkerStyle(20)  # Filled circle
    h_powheg_norm.SetMarkerSize(0.8)
    h_powheg_norm.SetMarkerColor(ROOT.kBlue)
    
    # Style for EFT SM histogram
    h_eft_sm_norm.SetLineColor(ROOT.kRed)
    h_eft_sm_norm.SetLineWidth(2)
    h_eft_sm_norm.SetMarkerStyle(21)  # Filled square
    h_eft_sm_norm.SetMarkerSize(0.8)
    h_eft_sm_norm.SetMarkerColor(ROOT.kRed)
    
    # Set custom or automatic y-axis range for main plot
    if custom_range:
        min_y, max_y = custom_range
        h_powheg_norm.GetYaxis().SetRangeUser(min_y, max_y)
    else:
        # Find the maximum y value to set the range
        max_val = max(h_powheg_norm.GetMaximum(), h_eft_sm_norm.GetMaximum())
        h_powheg_norm.GetYaxis().SetRangeUser(0, max_val * 1.2)
    
    # Draw normalized distributions
    h_powheg_norm.Draw("HIST")
    h_powheg_norm.Draw("E1 SAME")
    h_eft_sm_norm.Draw("HIST SAME")
    h_eft_sm_norm.Draw("E1 SAME")
    
    # Add legend
    leg = ROOT.TLegend(0.65, 0.70, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_powheg_norm, "Powheg SM", "lep")
    leg.AddEntry(h_eft_sm_norm, "EFT SM", "lep")
    leg.Draw()
    
    # Add CMS text
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextSize(0.05)
    cms_text.DrawLatex(0.20, 0.92, "#bf{CMS} #it{Simulation}")
    
    # Draw ratio in bottom pad
    pad_bottom.cd()
    
    # Style for ratio plot
    h_ratio.SetTitle("")
    h_ratio.GetYaxis().SetTitle("Powheg/EFT")
    h_ratio.GetXaxis().SetTitle(h_powheg.GetXaxis().GetTitle())  # Use original histogram title
    
    # Adjust text sizes for bottom pad
    h_ratio.GetXaxis().SetTitleSize(0.12)  # Slightly reduced from 0.14
    h_ratio.GetXaxis().SetLabelSize(0.10)  # Slightly reduced from 0.12
    h_ratio.GetYaxis().SetTitleSize(0.10)  # Slightly reduced from 0.12
    h_ratio.GetYaxis().SetLabelSize(0.09)  # Slightly reduced from 0.10
    h_ratio.GetYaxis().SetTitleOffset(0.6)  # Slightly increased from 0.5
    h_ratio.GetYaxis().SetNdivisions(505)   # Fewer divisions
    
    # Set reasonable y-axis range for ratio - exactly as in original
    min_ratio = h_ratio.GetMinimum(0.1)  # Ignore zeros
    max_ratio = h_ratio.GetMaximum()
    padding = 0.3 * (max_ratio - min_ratio)  # Increased padding from 0.2 to 0.3
    
    # Ensure ratio range has sensible values regardless of actual min/max
    min_y = max(0.5, min_ratio - padding)
    max_y = min(1.5, max_ratio + padding)
    
    # Add a bit more padding to make the plot less cramped
    range_height = max_y - min_y
    min_y = max(0.4, min_y - 0.1 * range_height)
    max_y = min(1.6, max_y + 0.1 * range_height)
    
    h_ratio.GetYaxis().SetRangeUser(min_y, max_y)
    
    h_ratio.SetLineColor(ROOT.kBlue)
    h_ratio.SetLineWidth(2)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(1.0)  # Increased from 0.8
    h_ratio.SetMarkerColor(ROOT.kBlue)
    
    # Draw ratio with error bars and line
    h_ratio.Draw("HIST")
    h_ratio.Draw("E1 SAME")
    
    # Add horizontal line at 1
    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(2)  # Dashed line
    line.SetLineColor(ROOT.kGray+2)  # Darker gray instead of black
    line.SetLineWidth(1)  # Thinner line
    line.Draw("same")
    
    # Update canvas to draw everything
    c_combined.Update()
    
    # Create directories if they don't exist
    make_directory(output_dir)
    plots_dir = os.path.join(output_dir, "plots")
    make_directory(plots_dir)
    
    # Save plots
    c_combined.SaveAs(os.path.join(plots_dir, "{0}_Powheg_vs_EFT_SM_combined.pdf".format(variable_name)))
    
    # No longer saving ROOT files
    
    # Clean up
    c_combined.Close()

def process_tree_with_branch(root_file, branch_name, histogram, weight_branch=None, max_events=-1, print_freq=10000):
    """Process a ROOT file, extract a branch of interest, and fill a histogram.
    
    Args:
        root_file (str): Path to ROOT file.
        branch_name (str): Name of branch to extract.
        histogram (TH1F): Histogram to fill with branch values.
        weight_branch (str, optional): Name of branch to use for weights. If None, weight is 1.0.
                                       If "SM", use SM weights from genInfo.
        max_events (int, optional): Maximum number of events to process. If -1, process all events.
        print_freq (int, optional): Frequency for printing progress updates.
    """
    # Open the ROOT file
    try:
        f = ROOT.TFile.Open(root_file, "READ")
        if not f or f.IsZombie() or not f.GetListOfKeys().Contains("AnalysisTree"):
            print("Error: Could not open file {0} or AnalysisTree not found".format(root_file))
            return
        
        tree = f.Get("AnalysisTree")
        n_entries = tree.GetEntries()
        print("  File contains {0} entries".format(n_entries))
        
        # Set up TTreeReader and TTreeReaderValue for the branch
        reader = ROOT.TTreeReader("AnalysisTree", f)
        
        # Use direct branch access instead of TTreeReaderValue which is causing issues
        reader.SetEntry(0)
        
        # Process entries directly using TTree
        tree.SetBranchStatus("*", 0)  # Disable all branches
        tree.SetBranchStatus(branch_name, 1)  # Enable the branch we need
        
        if weight_branch and weight_branch != "SM":
            tree.SetBranchStatus(weight_branch, 1)  # Enable weight branch if needed
            
        if weight_branch == "SM":
            tree.SetBranchStatus("genInfo*", 1)  # Enable genInfo for SM weights
        
        # Determine the actual number of entries to process
        actual_max = n_entries if max_events < 0 else min(max_events, n_entries)
        
        # Check what structure genInfo has
        if weight_branch == "SM" and tree.GetEntry(0) > 0:
            has_systweights = False
            try:
                if hasattr(tree, "genInfo"):
                    # Check for different possible structures
                    if hasattr(tree.genInfo, "systweights"):
                        print("  Using genInfo.systweights() for SM weights")
                        has_systweights = True
                    elif hasattr(tree.genInfo, "weights"):
                        print("  Using genInfo.weights() for SM weights")
                        has_systweights = True
                    else:
                        print("  WARNING: genInfo found but doesn't have systweights or weights method")
            except Exception as e:
                print("  WARNING: Error checking genInfo structure: {e}".format(e))
        
        # Process entries
        processed = 0
        for entry in range(actual_max):
            tree.GetEntry(entry)
            
            # Get branch value
            branch_val = getattr(tree, branch_name)
            
            # Get weight
            weight = 1.0
            if weight_branch and weight_branch != "SM":
                weight = getattr(tree, weight_branch)
            elif weight_branch == "SM" and hasattr(tree, "genInfo"):
                try:
                    genInfo = getattr(tree, "genInfo")
                    # Try different ways to access weights depending on structure
                    if hasattr(genInfo, "systweights") and callable(getattr(genInfo, "systweights")):
                        weights = genInfo.systweights()
                        if len(weights) > 0:
                            weight = weights[0]  # Use first weight
                    elif hasattr(genInfo, "weights") and callable(getattr(genInfo, "weights")):
                        weights = genInfo.weights()
                        if len(weights) > 0:
                            weight = weights[0]  # Use first weight
                    # If we find LHEWeights, try to use that
                    elif hasattr(tree, "LHEWeights") and hasattr(tree.LHEWeights, "size"):
                        if tree.LHEWeights.size() > 0:
                            weight = tree.LHEWeights[0]
                except Exception as e:
                    if entry == 0:
                        print("  WARNING: Error extracting weight from genInfo: {e}".format(e))
                        print("  Using weight=1.0 for all events".format(e))
            
            # Fill the histogram with the value and weight
            histogram.Fill(branch_val, weight)
            
            # Increment counter
            processed += 1
            
            # Print progress
            if processed % print_freq == 0:
                print("    Processed {0}/{1} entries ({2:.2f}%)".format(
                    processed, actual_max, float(processed) / actual_max * 100))
        
        print("  Processed {0} entries, filled {1} in histogram".format(processed, histogram.GetEntries()))
        
    finally:
        # Clean up
        if 'f' in locals() and f:
            f.Close()

def process_sigmaPhi1SR(powheg_files, eft_files, output_dir, max_events=-1, print_freq=10000):
    # Process Sigma_phi_1_SR variable from Powheg and EFT SM samples.
    print("Processing Sigma_phi_1_SR variable...")
    
    # Create histograms with full angular range -pi to pi
    h_sigmaPhi1SR_powheg = ROOT.TH1F("h_sigmaPhi1SR_powheg", "Powheg SM #Sigma#phi_{1} SR;#Sigma#phi_{1} SR;Events", 16, -3.2, 3.2)
    h_sigmaPhi1SR_eft_sm = ROOT.TH1F("h_sigmaPhi1SR_eft_sm", "EFT SM #Sigma#phi_{1} SR;#Sigma#phi_{1} SR;Events", 16, -3.2, 3.2)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {i}/{total}: {filename}".format(
            i=i+1, total=len(powheg_files), filename=os.path.basename(root_file)))
        process_tree_with_branch(root_file, "Sigma_phi_1_SR", h_sigmaPhi1SR_powheg, 
                                 weight_branch=None, max_events=max_events, print_freq=print_freq)
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {i}/{total}: {filename}".format(
            i=i+1, total=len(eft_files), filename=os.path.basename(root_file)))
        process_tree_with_branch(root_file, "Sigma_phi_1_SR", h_sigmaPhi1SR_eft_sm, 
                                 weight_branch="SM", max_events=max_events, print_freq=print_freq)
    
    # Create comparison plots
    create_comparison_plot(h_sigmaPhi1SR_powheg, h_sigmaPhi1SR_eft_sm, "sigmaPhi1SR", output_dir)

def process_sigmaPhi2SR(powheg_files, eft_files, output_dir, max_events=-1, print_freq=10000):
    # Process Sigma_phi_2_SR variable from Powheg and EFT SM samples.
    print("Processing Sigma_phi_2_SR variable...")
    
    # Create histograms with full angular range -pi to pi
    h_sigmaPhi2SR_powheg = ROOT.TH1F("h_sigmaPhi2SR_powheg", "Powheg SM #Sigma#phi_{2} SR;#Sigma#phi_{2} SR;Events", 16, -3.2, 3.2)
    h_sigmaPhi2SR_eft_sm = ROOT.TH1F("h_sigmaPhi2SR_eft_sm", "EFT SM #Sigma#phi_{2} SR;#Sigma#phi_{2} SR;Events", 16, -3.2, 3.2)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {i}/{total}: {filename}".format(
            i=i+1, total=len(powheg_files), filename=os.path.basename(root_file)))
        process_tree_with_branch(root_file, "Sigma_phi_2_SR", h_sigmaPhi2SR_powheg, 
                                 weight_branch=None, max_events=max_events, print_freq=print_freq)
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {i}/{total}: {filename}".format(
            i=i+1, total=len(eft_files), filename=os.path.basename(root_file)))
        process_tree_with_branch(root_file, "Sigma_phi_2_SR", h_sigmaPhi2SR_eft_sm, 
                                 weight_branch="SM", max_events=max_events, print_freq=print_freq)
    
    # Create comparison plots
    create_comparison_plot(h_sigmaPhi2SR_powheg, h_sigmaPhi2SR_eft_sm, "sigmaPhi2SR", output_dir)

def process_deltaY1SR(powheg_files, eft_files, output_dir, max_events=-1, print_freq=10000):
    # Process dyreco_1_SR variable from Powheg and EFT SM samples.
    print("Processing dyreco_1_SR variable...")
    
    # Create histograms with just 2 bins - one for negative and one for positive values
    h_deltaY1SR_powheg = ROOT.TH1F("h_deltaY1SR_powheg", "Powheg SM #DeltaY_{1} SR;#DeltaY_{1} SR;Events", 2, -3, 3)
    h_deltaY1SR_eft_sm = ROOT.TH1F("h_deltaY1SR_eft_sm", "EFT SM #DeltaY_{1} SR;#DeltaY_{1} SR;Events", 2, -3, 3)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {i}/{total}: {filename}".format(
            i=i+1, total=len(powheg_files), filename=os.path.basename(root_file)))
        process_tree_with_branch(root_file, "dyreco_1_SR", h_deltaY1SR_powheg, 
                                 weight_branch=None, max_events=max_events, print_freq=print_freq)
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {i}/{total}: {filename}".format(
            i=i+1, total=len(eft_files), filename=os.path.basename(root_file)))
        process_tree_with_branch(root_file, "dyreco_1_SR", h_deltaY1SR_eft_sm, 
                                 weight_branch="SM", max_events=max_events, print_freq=print_freq)
    
    # Create comparison plots with custom y-axis range from 0.4 to 0.6 for the main plot only
    create_comparison_plot(h_deltaY1SR_powheg, h_deltaY1SR_eft_sm, "deltaY1SR", output_dir, custom_range=(0.4, 0.6))

def process_deltaY2SR(powheg_files, eft_files, output_dir, max_events=-1, print_freq=10000):
    # Process dyreco_2_SR variable from Powheg and EFT SM samples.
    print("Processing dyreco_2_SR variable...")
    
    # Create histograms with just 2 bins - one for negative and one for positive values
    h_deltaY2SR_powheg = ROOT.TH1F("h_deltaY2SR_powheg", "Powheg SM #DeltaY_{2} SR;#DeltaY_{2} SR;Events", 2, -3, 3)
    h_deltaY2SR_eft_sm = ROOT.TH1F("h_deltaY2SR_eft_sm", "EFT SM #DeltaY_{2} SR;#DeltaY_{2} SR;Events", 2, -3, 3)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {i}/{total}: {filename}".format(
            i=i+1, total=len(powheg_files), filename=os.path.basename(root_file)))
        process_tree_with_branch(root_file, "dyreco_2_SR", h_deltaY2SR_powheg, 
                                 weight_branch=None, max_events=max_events, print_freq=print_freq)
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {i}/{total}: {filename}".format(
            i=i+1, total=len(eft_files), filename=os.path.basename(root_file)))
        process_tree_with_branch(root_file, "dyreco_2_SR", h_deltaY2SR_eft_sm, 
                                 weight_branch="SM", max_events=max_events, print_freq=print_freq)
    
    # Create comparison plots with custom y-axis range from 0.4 to 0.6 for the main plot only
    create_comparison_plot(h_deltaY2SR_powheg, h_deltaY2SR_eft_sm, "deltaY2SR", output_dir, custom_range=(0.4, 0.6))

def main():
    args = parse_arguments()
    
    # Add timing
    start_time_total = time.time()
    
    # Find input files
    powheg_files = find_files(args.powheg_dir, max_files=args.max_files)
    eft_files = find_files(args.eft_dir, max_files=args.max_files)
    
    if not powheg_files:
        print("No ROOT files found in Powheg directory: {0}".format(args.powheg_dir))
        return
    
    if not eft_files:
        print("No ROOT files found in EFT directory: {0}".format(args.eft_dir))
        return
    
    print("Found {0} Powheg files and {1} EFT files".format(len(powheg_files), len(eft_files)))
    
    # Create output directory
    make_directory(args.output_dir)
    
    # Process variables based on command line flags
    # if (not args.delta_phi_only and not args.sigma_phi_only) or args.delta_y_only:
    #     process_deltaY(powheg_files, eft_files, args.output_dir, args.max_events)
    
    # Commented out deltaPhi processing
    # if (not args.delta_y_only and not args.sigma_phi_only) or args.delta_phi_only:
    #     process_deltaPhi(powheg_files, eft_files, args.output_dir, args.max_events)
    
    # Commented out sigmaPhi
    # if (not args.delta_y_only and not args.delta_phi_only) or args.sigma_phi_only:
    #     process_sigmaPhi(powheg_files, eft_files, args.output_dir, args.max_events)
    
    # Process new variables
    # process_sigmaPhi1SR(powheg_files, eft_files, args.output_dir, args.max_events, args.print_frequency)
    # process_sigmaPhi2SR(powheg_files, eft_files, args.output_dir, args.max_events, args.print_frequency)
    process_deltaY1SR(powheg_files, eft_files, args.output_dir, args.max_events, args.print_frequency)
    process_deltaY2SR(powheg_files, eft_files, args.output_dir, args.max_events, args.print_frequency)
    
    # Report total time
    elapsed_total = time.time() - start_time_total
    print("All comparison plots saved to {output_dir}".format(output_dir=args.output_dir))
    print("Total execution time: {elapsed:.1f} seconds ({minutes:.1f} minutes)".format(
        elapsed=elapsed_total, minutes=elapsed_total/60))

if __name__ == "__main__":
    main() 