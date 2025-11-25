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
    return parser.parse_args()

def find_files(directory, pattern="*.root"):
    # all files matching pattern in directory.
    return glob.glob(os.path.join(directory, pattern))

def make_directory(directory):
    # Create directory if it doesn't exist.
    if not os.path.exists(directory):
        os.makedirs(directory)

def create_comparison_plot(h_powheg, h_eft_sm, variable_name, output_dir, title="Powheg SM vs EFT SM"):
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
    
    # Set reasonable y-axis range for ratio
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
        
    # Clean up
    c_combined.Close()

def process_deltaY(powheg_files, eft_files, output_dir, max_events=-1):
    # deltaY from ttree of Powheg and EFT SM samples
    print("Processing deltaY variable...")
    
    # Create histograms with just 2 bins - one for negative and one for positive values
    h_deltaY_powheg = ROOT.TH1F("h_deltaY_powheg", "Powheg SM #DeltaY_{reco};#DeltaY_{reco};Events", 2, -3, 3)
    h_deltaY_eft_sm = ROOT.TH1F("h_deltaY_eft_sm", "EFT SM #DeltaY_{reco};#DeltaY_{reco};Events", 2, -3, 3)
    
    # Fill histogram from Powheg samples (using dyreco as branch name)
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {0}/{1}: {2}".format(i+1, len(powheg_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if dyreco branch exists for Powheg samples
        if not tree.GetBranch("dyreco"):
            print("Branch 'dyreco' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # Fill histogram
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            h_deltaY_powheg.Fill(tree.dyreco)
        
        f.Close()
    
    # Fill histogram from EFT SM samples (using dyreco as branch name)
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {0}/{1}: {2}".format(i+1, len(eft_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if dyreco branch exists for EFT samples
        if not tree.GetBranch("dyreco"):
            print("Branch 'dyreco' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # For EFT samples, we need to use the SM weight (genInfo.systweights()[202])
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            # Only fill if SM weight is available
            if hasattr(tree, 'genInfo') and tree.genInfo.systweights().size() > 202:
                sm_weight = tree.genInfo.systweights()[202]
                h_deltaY_eft_sm.Fill(tree.dyreco, sm_weight)
            # Skip events without valid weight
        
        f.Close()
    
    # Create comparison plots
    create_comparison_plot(h_deltaY_powheg, h_deltaY_eft_sm, "deltaY", output_dir)

def process_deltaPhi(powheg_files, eft_files, output_dir, max_events=-1):
    # deltaPhi from ttree of Powheg and EFT SM samples
    print("Processing deltaPhi variable...")
    
    # Create histograms
    h_deltaPhi_powheg = ROOT.TH1F("h_deltaPhi_powheg", "Powheg SM #Delta#phi_{reco};#Delta#phi_{reco};Events", 16, 0, 3.2)
    h_deltaPhi_eft_sm = ROOT.TH1F("h_deltaPhi_eft_sm", "EFT SM #Delta#phi_{reco};#Delta#phi_{reco};Events", 16, 0, 3.2)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {0}/{1}: {2}".format(i+1, len(powheg_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if Delta_phi branch exists (as used in original scripts)
        has_delta_phi = tree.GetBranch("Delta_phi") is not None
        has_deltaphi_reco = tree.GetBranch("DeltaPhi_reco") is not None
        
        if not has_delta_phi and not has_deltaphi_reco:
            print("Neither 'Delta_phi' nor 'DeltaPhi_reco' branch found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # Fill histogram based on which branch is available
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            if has_delta_phi:
                h_deltaPhi_powheg.Fill(tree.Delta_phi)
            else:
                h_deltaPhi_powheg.Fill(tree.DeltaPhi_reco)
        
        f.Close()
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {0}/{1}: {2}".format(i+1, len(eft_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if Delta_phi branch exists (as used in original scripts)
        has_delta_phi = tree.GetBranch("Delta_phi") is not None
        has_deltaphi_reco = tree.GetBranch("DeltaPhi_reco") is not None
        
        if not has_delta_phi and not has_deltaphi_reco:
            print("Neither 'Delta_phi' nor 'DeltaPhi_reco' branch found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # For EFT samples, we need to use the SM weight (genInfo.systweights()[202])
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            # Only fill if SM weight is available
            if hasattr(tree, 'genInfo') and tree.genInfo.systweights().size() > 202:
                sm_weight = tree.genInfo.systweights()[202]
                # Fill based on available branch
                if has_delta_phi:
                    h_deltaPhi_eft_sm.Fill(tree.Delta_phi, sm_weight)
                else:
                    h_deltaPhi_eft_sm.Fill(tree.DeltaPhi_reco, sm_weight)
            # Skip events without valid weight
        
        f.Close()
    
    # Create comparison plots
    create_comparison_plot(h_deltaPhi_powheg, h_deltaPhi_eft_sm, "deltaPhi", output_dir)

def process_sigmaPhi(powheg_files, eft_files, output_dir, max_events=-1):
    # Process sigmaPhi variable from Powheg and EFT SM samples.
    print("Processing sigmaPhi variable...")
    
    # Create histograms with updated axis title using capital Sigma
    h_sigmaPhi_powheg = ROOT.TH1F("h_sigmaPhi_powheg", "Powheg SM #Sigma#phi_{reco};#Sigma#phi_{reco};Events", 10, 0, 1.0)
    h_sigmaPhi_eft_sm = ROOT.TH1F("h_sigmaPhi_eft_sm", "EFT SM #Sigma#phi_{reco};#Sigma#phi_{reco};Events", 10, 0, 1.0)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {0}/{1}: {2}".format(i+1, len(powheg_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if Sigma_phi branch exists (as used in original scripts)
        has_sigma_phi = tree.GetBranch("Sigma_phi") is not None
        has_sigmaphi_reco = tree.GetBranch("SigmaPhi_reco") is not None
        
        if not has_sigma_phi and not has_sigmaphi_reco:
            print("Neither 'Sigma_phi' nor 'SigmaPhi_reco' branch found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # Fill histogram based on which branch is available
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            if has_sigma_phi:
                h_sigmaPhi_powheg.Fill(tree.Sigma_phi)
            else:
                h_sigmaPhi_powheg.Fill(tree.SigmaPhi_reco)
        
        f.Close()
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {0}/{1}: {2}".format(i+1, len(eft_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if Sigma_phi branch exists (as used in original scripts)
        has_sigma_phi = tree.GetBranch("Sigma_phi") is not None
        has_sigmaphi_reco = tree.GetBranch("SigmaPhi_reco") is not None
        
        if not has_sigma_phi and not has_sigmaphi_reco:
            print("Neither 'Sigma_phi' nor 'SigmaPhi_reco' branch found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # For EFT samples, we need to use the SM weight (genInfo.systweights()[202])
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            # Only fill if SM weight is available
            if hasattr(tree, 'genInfo') and tree.genInfo.systweights().size() > 202:
                sm_weight = tree.genInfo.systweights()[202]
                # Fill based on available branch
                if has_sigma_phi:
                    h_sigmaPhi_eft_sm.Fill(tree.Sigma_phi, sm_weight)
                else:
                    h_sigmaPhi_eft_sm.Fill(tree.SigmaPhi_reco, sm_weight)
            # Skip events without valid weight
        
        f.Close()
    
    # Create comparison plots
    create_comparison_plot(h_sigmaPhi_powheg, h_sigmaPhi_eft_sm, "sigmaPhi", output_dir)

def process_sigmaPhi1SR(powheg_files, eft_files, output_dir, max_events=-1):
    # Process Sigma_phi_1_SR variable from Powheg and EFT SM samples.
    print("Processing Sigma_phi_1_SR variable...")
    
    # Create histograms with full angular range -pi to pi
    h_sigmaPhi1SR_powheg = ROOT.TH1F("h_sigmaPhi1SR_powheg", "Powheg SM #Sigma#phi_{1} SR;#Sigma#phi_{1} SR;Events", 16, -3.2, 3.2)
    h_sigmaPhi1SR_eft_sm = ROOT.TH1F("h_sigmaPhi1SR_eft_sm", "EFT SM #Sigma#phi_{1} SR;#Sigma#phi_{1} SR;Events", 16, -3.2, 3.2)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {0}/{1}: {2}".format(i+1, len(powheg_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if Sigma_phi_1_SR branch exists
        if not tree.GetBranch("Sigma_phi_1_SR"):
            print("Branch 'Sigma_phi_1_SR' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # Fill histogram
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            h_sigmaPhi1SR_powheg.Fill(tree.Sigma_phi_1_SR)
        
        f.Close()
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {0}/{1}: {2}".format(i+1, len(eft_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if Sigma_phi_1_SR branch exists
        if not tree.GetBranch("Sigma_phi_1_SR"):
            print("Branch 'Sigma_phi_1_SR' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # For EFT samples, we need to use the SM weight (genInfo.systweights()[202])
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            # Only fill if SM weight is available
            if hasattr(tree, 'genInfo') and tree.genInfo.systweights().size() > 202:
                sm_weight = tree.genInfo.systweights()[202]
                h_sigmaPhi1SR_eft_sm.Fill(tree.Sigma_phi_1_SR, sm_weight)
            # Skip events without valid weight
        
        f.Close()
    
    # Create comparison plots
    create_comparison_plot(h_sigmaPhi1SR_powheg, h_sigmaPhi1SR_eft_sm, "sigmaPhi1SR", output_dir)

def process_sigmaPhi2SR(powheg_files, eft_files, output_dir, max_events=-1):
    # Process Sigma_phi_2_SR variable from Powheg and EFT SM samples.
    print("Processing Sigma_phi_2_SR variable...")
    
    # Create histograms with full angular range -pi to pi
    h_sigmaPhi2SR_powheg = ROOT.TH1F("h_sigmaPhi2SR_powheg", "Powheg SM #Sigma#phi_{2} SR;#Sigma#phi_{2} SR;Events", 16, -3.2, 3.2)
    h_sigmaPhi2SR_eft_sm = ROOT.TH1F("h_sigmaPhi2SR_eft_sm", "EFT SM #Sigma#phi_{2} SR;#Sigma#phi_{2} SR;Events", 16, -3.2, 3.2)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {0}/{1}: {2}".format(i+1, len(powheg_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if Sigma_phi_2_SR branch exists
        if not tree.GetBranch("Sigma_phi_2_SR"):
            print("Branch 'Sigma_phi_2_SR' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # Fill histogram
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            h_sigmaPhi2SR_powheg.Fill(tree.Sigma_phi_2_SR)
        
        f.Close()
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {0}/{1}: {2}".format(i+1, len(eft_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if Sigma_phi_2_SR branch exists
        if not tree.GetBranch("Sigma_phi_2_SR"):
            print("Branch 'Sigma_phi_2_SR' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # For EFT samples, we need to use the SM weight (genInfo.systweights()[202])
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            # Only fill if SM weight is available
            if hasattr(tree, 'genInfo') and tree.genInfo.systweights().size() > 202:
                sm_weight = tree.genInfo.systweights()[202]
                h_sigmaPhi2SR_eft_sm.Fill(tree.Sigma_phi_2_SR, sm_weight)
            # Skip events without valid weight
        
        f.Close()
    
    # Create comparison plots
    create_comparison_plot(h_sigmaPhi2SR_powheg, h_sigmaPhi2SR_eft_sm, "sigmaPhi2SR", output_dir)

def process_deltaY1SR(powheg_files, eft_files, output_dir, max_events=-1):
    # Process dyreco_1_SR variable from Powheg and EFT SM samples.
    print("Processing dyreco_1_SR variable...")
    
    # Create histograms with just 2 bins - one for negative and one for positive values
    h_deltaY1SR_powheg = ROOT.TH1F("h_deltaY1SR_powheg", "Powheg SM #DeltaY_{1} SR;#DeltaY_{1} SR;Events", 2, -3, 3)
    h_deltaY1SR_eft_sm = ROOT.TH1F("h_deltaY1SR_eft_sm", "EFT SM #DeltaY_{1} SR;#DeltaY_{1} SR;Events", 2, -3, 3)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {0}/{1}: {2}".format(i+1, len(powheg_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if dyreco_1_SR branch exists
        if not tree.GetBranch("dyreco_1_SR"):
            print("Branch 'dyreco_1_SR' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # Fill histogram
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            h_deltaY1SR_powheg.Fill(tree.dyreco_1_SR)
        
        f.Close()
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {0}/{1}: {2}".format(i+1, len(eft_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if dyreco_1_SR branch exists
        if not tree.GetBranch("dyreco_1_SR"):
            print("Branch 'dyreco_1_SR' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # For EFT samples, we need to use the SM weight (genInfo.systweights()[202])
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            # Only fill if SM weight is available
            if hasattr(tree, 'genInfo') and tree.genInfo.systweights().size() > 202:
                sm_weight = tree.genInfo.systweights()[202]
                h_deltaY1SR_eft_sm.Fill(tree.dyreco_1_SR, sm_weight)
            # Skip events without valid weight
        
        f.Close()
    
    # Create comparison plots
    create_comparison_plot(h_deltaY1SR_powheg, h_deltaY1SR_eft_sm, "deltaY1SR", output_dir)

def process_deltaY2SR(powheg_files, eft_files, output_dir, max_events=-1):
    # Process dyreco_2_SR variable from Powheg and EFT SM samples.
    print("Processing dyreco_2_SR variable...")
    
    # Create histograms with just 2 bins - one for negative and one for positive values
    h_deltaY2SR_powheg = ROOT.TH1F("h_deltaY2SR_powheg", "Powheg SM #DeltaY_{2} SR;#DeltaY_{2} SR;Events", 2, -3, 3)
    h_deltaY2SR_eft_sm = ROOT.TH1F("h_deltaY2SR_eft_sm", "EFT SM #DeltaY_{2} SR;#DeltaY_{2} SR;Events", 2, -3, 3)
    
    # Fill histogram from Powheg samples
    for i, root_file in enumerate(powheg_files):
        print("  Processing Powheg file {0}/{1}: {2}".format(i+1, len(powheg_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if dyreco_2_SR branch exists
        if not tree.GetBranch("dyreco_2_SR"):
            print("Branch 'dyreco_2_SR' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # Fill histogram
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            h_deltaY2SR_powheg.Fill(tree.dyreco_2_SR)
        
        f.Close()
    
    # Fill histogram from EFT SM samples
    for i, root_file in enumerate(eft_files):
        print("  Processing EFT file {0}/{1}: {2}".format(i+1, len(eft_files), os.path.basename(root_file)))
        
        f = ROOT.TFile(root_file, "READ")
        if not f or f.IsZombie():
            print("Error opening file: {0}".format(root_file))
            continue
        
        tree = f.Get("AnalysisTree")
        if not tree:
            print("Tree 'AnalysisTree' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Check if dyreco_2_SR branch exists
        if not tree.GetBranch("dyreco_2_SR"):
            print("Branch 'dyreco_2_SR' not found in {0}".format(root_file))
            f.Close()
            continue
        
        # Get number of entries to process
        entries = tree.GetEntries()
        n_to_process = entries if max_events < 0 else min(entries, max_events)
        print("    Processing {0} events out of {1}".format(n_to_process, entries))
        
        # For EFT samples, we need to use the SM weight (genInfo.systweights()[202])
        for i in range(n_to_process):
            if i > 0 and i % 10000 == 0:
                print("    ... processed {0}/{1} events".format(i, n_to_process))
            
            tree.GetEntry(i)
            # Only fill if SM weight is available
            if hasattr(tree, 'genInfo') and tree.genInfo.systweights().size() > 202:
                sm_weight = tree.genInfo.systweights()[202]
                h_deltaY2SR_eft_sm.Fill(tree.dyreco_2_SR, sm_weight)
            # Skip events without valid weight
        
        f.Close()
    
    # Create comparison plots
    create_comparison_plot(h_deltaY2SR_powheg, h_deltaY2SR_eft_sm, "deltaY2SR", output_dir)

def main():
    args = parse_arguments()
    
    # Find input files
    powheg_files = find_files(args.powheg_dir)
    eft_files = find_files(args.eft_dir)
    
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
    process_sigmaPhi1SR(powheg_files, eft_files, args.output_dir, args.max_events)
    process_sigmaPhi2SR(powheg_files, eft_files, args.output_dir, args.max_events)
    process_deltaY1SR(powheg_files, eft_files, args.output_dir, args.max_events)
    process_deltaY2SR(powheg_files, eft_files, args.output_dir, args.max_events)
    
    print("All comparison plots saved to {0}".format(args.output_dir))

if __name__ == "__main__":
    main() 