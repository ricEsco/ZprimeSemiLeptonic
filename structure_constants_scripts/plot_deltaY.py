# Script to make detlaY plots and ratios compared to SM
# ./plot_deltaY.py EFT_700_900.root --output deltaY_ctGRe1.pdf

import ROOT
import numpy as np
import sys
import argparse
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

def calculate_weight(structure_constants, wc_values):
    # Calculate weight using structure constants and Wilson coefficient values
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

def create_histogram(name, title, nbins, xmin, xmax):
    # Create a ROOT histogram with the given parameters
    return ROOT.TH1F(name, title, nbins, xmin, xmax)

def set_histogram_style(hist, color, line_style=1, marker_style=20, fill_style=0):
    # Set histogram style
    hist.SetLineColor(color)
    hist.SetLineStyle(line_style)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(marker_style)
    hist.SetFillStyle(fill_style)
    hist.SetLineWidth(2)

def create_ratio_plot(h1, h2, name, title, y_title="Ratio"):
    # Create a ratio plot of h1/h2
    ratio = h1.Clone(name)
    ratio.SetTitle(title)
    ratio.Divide(h2)
    ratio.GetYaxis().SetTitle(y_title)
    return ratio

def create_canvas_with_ratio(name, title, w=800, h=800):
    Create a canvas with a ratio panel
    c = ROOT.TCanvas(name, title, w, h)
    
    # pads
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.12)
    pad1.SetRightMargin(0.04)
    pad1.Draw()
    
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.3)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.3)
    pad2.SetLeftMargin(0.12)
    pad2.SetRightMargin(0.04)
    pad2.Draw()
    
    return c, pad1, pad2

def main():
    parser = argparse.ArgumentParser(description="Create Delta Y plots with ctGRe=1 compared to SM")
    parser.add_argument("root_file", help="Path to ROOT file containing structure constants")
    parser.add_argument("--output", default="deltaY_plot.pdf", help="Output file name")
    parser.add_argument("--output-root", default="deltaY_histograms.root", help="Output ROOT file for histograms")
    parser.add_argument("--max-events", type=int, default=-1, help="Maximum number of events to process (-1 for all)")
    parser.add_argument("--region", default="SR", choices=["SR", "CR1", "CR2"], help="Region to plot")
    parser.add_argument("--mass-bins", action="store_true", help="Create plots for different mass bins")
    args = parser.parse_args()
    
    input_file = ROOT.TFile.Open(args.root_file)
    if not input_file or input_file.IsZombie():
        print "Error: Could not open ROOT file: %s" % args.root_file
        return 1
    
    # Get the tree
    tree = input_file.Get("AnalysisTree")
    if not tree:
        print "Error: Could not find TTree: AnalysisTree"
        input_file.Close()
        return 1
    
    # Check if structure_constants branch exists
    branch_names = [b.GetName() for b in tree.GetListOfBranches()]
    if "structure_constants" not in branch_names:
        print "Error: structure_constants branch not found in TTree"
        print "Available branches:", branch_names
        input_file.Close()
        return 1
    
    # Define variable names based on region
    region_suffix = "" if args.region == "SR" else "_" + args.region
    dyreco_var = "dyreco" + region_suffix
    dyreco_1_var = "dyreco_1" + region_suffix
    dyreco_2_var = "dyreco_2" + region_suffix
    
    # Check if variables exist
    for var in [dyreco_var, dyreco_1_var, dyreco_2_var]:
        if var not in branch_names:
            print "Warning: Variable %s not found in tree" % var
    
    # Define WC scenarios
    wc_scenarios = {
        "SM": [0] * 16,
        "ctGRe=1": [1] + [0] * 15
    }
    
    # Create histograms
    histograms = {}
    
    # Define histogram parameters
    nbins = 14
    xmin = -3
    xmax = 3
    
    # Create histograms for each variable and scenario
    # for var_name in [dyreco, dyreco_1, dyreco_2]:
    for var_name in [deltay]:

        if var_name not in branch_names:
            continue
        
        for scenario_name in wc_scenarios:
            hist_name = "{}_{}".format(var_name, scenario_name)
            hist_title = "{} ({})".format(var_name, scenario_name)
            histograms[hist_name] = create_histogram(hist_name, hist_title, nbins, xmin, xmax)
    
    # Process events
    nentries = tree.GetEntries()
    max_events = nentries if args.max_events < 0 else min(args.max_events, nentries)
    print "Processing {} events out of {} total...".format(max_events, nentries)
    
    for i in range(max_events):
        if i % 10000 == 0:
            print "  Event {} / {}".format(i, max_events)
        
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
        for var_name in [dyreco_var, dyreco_1_var, dyreco_2_var]:
            if var_name not in branch_names:
                continue
            
            # Get variable value
            var_value = getattr(tree, var_name)
            
            # Fill histograms for each scenario
            for scenario_name in wc_scenarios:
                hist_name = "{}_{}".format(var_name, scenario_name)
                if hist_name in histograms:
                    histograms[hist_name].Fill(var_value, event_weight * scenario_weights[scenario_name])
    
    # Create output ROOT file
    output_root = ROOT.TFile(args.output_root, "RECREATE")
    
    # Save histograms to ROOT file
    for hist in histograms.values():
        hist.Write()
    
    # Create plots
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    
    # Create directory for plots if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Create plots for each variable
    for var_name in [dyreco_var, dyreco_1_var, dyreco_2_var]:
        if var_name not in branch_names:
            continue
        
        # Get histograms
        sm_hist_name = "{}_SM".format(var_name)
        ctgre_hist_name = "{}_ctGRe=1".format(var_name)
        
        if sm_hist_name not in histograms or ctgre_hist_name not in histograms:
            continue
        
        sm_hist = histograms[sm_hist_name]
        ctgre_hist = histograms[ctgre_hist_name]
        
        # Normalize histograms to unit area
        sm_hist.Scale(1.0 / sm_hist.Integral() if sm_hist.Integral() > 0 else 1.0)
        ctgre_hist.Scale(1.0 / ctgre_hist.Integral() if ctgre_hist.Integral() > 0 else 1.0)
        
        # Set histogram styles
        set_histogram_style(sm_hist, ROOT.kBlack)
        set_histogram_style(ctgre_hist, ROOT.kRed)
        
        # Create ratio plot
        ratio = create_ratio_plot(ctgre_hist, sm_hist, "%s_ratio" % var_name, "", "ctGRe=1 / SM")
        
        # Create canvas with ratio panel
        canvas_name = "c_%s" % var_name
        canvas_title = "Delta Y - %s" % var_name
        c, pad1, pad2 = create_canvas_with_ratio(canvas_name, canvas_title)
        
        # Draw histograms
        pad1.cd()
        sm_hist.SetTitle("")
        sm_hist.GetXaxis().SetLabelSize(0)
        sm_hist.GetYaxis().SetTitle("Normalized Events")
        sm_hist.GetYaxis().SetTitleSize(0.05)
        sm_hist.GetYaxis().SetTitleOffset(1.2)
        sm_hist.GetYaxis().SetLabelSize(0.045)
        sm_hist.Draw("HIST")
        ctgre_hist.Draw("HIST SAME")
        
        # Create legend
        legend = ROOT.TLegend(0.65, 0.75, 0.89, 0.89)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(sm_hist, "SM", "l")
        legend.AddEntry(ctgre_hist, "ctGRe=1", "l")
        legend.Draw()
        
        # Add CMS label
        cms_text = ROOT.TLatex()
        cms_text.SetNDC()
        cms_text.SetTextSize(0.05)
        cms_text.SetTextFont(42)
        cms_text.DrawLatex(0.12, 0.93, "#bf{CMS} #it{Simulation Preliminary}")
        
        # Draw ratio plot
        pad2.cd()
        ratio.SetTitle("")
        ratio.GetXaxis().SetTitle("#Delta y")
        ratio.GetXaxis().SetTitleSize(0.12)
        ratio.GetXaxis().SetTitleOffset(0.8)
        ratio.GetXaxis().SetLabelSize(0.1)
        ratio.GetYaxis().SetTitleSize(0.12)
        ratio.GetYaxis().SetTitleOffset(0.5)
        ratio.GetYaxis().SetLabelSize(0.1)
        ratio.GetYaxis().SetRangeUser(0.5, 1.5)
        ratio.GetYaxis().SetNdivisions(505)
        ratio.Draw("HIST")
        
        # Draw horizontal line at 1
        line = ROOT.TLine(xmin, 1, xmax, 1)
        line.SetLineColor(ROOT.kGray+2)
        line.SetLineStyle(2)
        line.Draw()
        
        # Save canvas
        output_file = args.output
        if var_name != dyreco_var:
            base, ext = os.path.splitext(args.output)
            output_file = "{}_{}{}".format(base, var_name, ext)
        
        c.SaveAs(output_file)
        
        # Save to ROOT file
        output_root.cd()
        c.Write()
    
    # Close files
    output_root.Close()
    input_file.Close()
    
    print "Done! Plots saved to {} and histograms saved to {}".format(args.output, args.output_root)
    return 0

if __name__ == "__main__":
    sys.exit(main()) 