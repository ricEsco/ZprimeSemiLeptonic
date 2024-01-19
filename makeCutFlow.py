#from ROOT import TFile, TLegend, TCanvas, TPad, THStack, TF1, TPaveText, TGaxis, SetOwnership, TObject, gStyle,TH1F
from ROOT import *
import os
import sys
import numpy as np
from numpy import log10
from array import array
import CMS_lumi
from Style import *

fileDir =      "/nfs/dust/cms/user/ricardo/analysis_output/RecoOpt/muon/workdir_Analysis_UL18_muon_RecOpt_signal"
plotDirectory = "/nfs/dust/cms/user/ricardo/analysis_output/RecoOpt/plots/UL18/muon"

print "The input root files will come from", fileDir
print "The output will go into", plotDirectory, "\n"

### This array contains the names of the folders in the root file
folder_names = ["Weights_Init_General", "Muon1_Tot_General", "TriggerMuon_General", "TwoDCut_Muon_General", "Jet1_General", "Jet2_General", "MET_General", "Btags1_General", "Chi2Discriminator_General"]

### This array has filenames of samples that will contribute to each plot
sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_2"]

### end User input section i.e. modify each time new set of plots are being generated ------------------------------------------------------------------------------------

### Histogram Information:
gStyle.SetOptStat(0)
H = 600
W = 800

# Create cutflow histogram
h_cutflow = TH1F("h_cutflow", "Cut Flow", len(folder_names), 0, len(folder_names))

# Label bins of cutflow histogram in order of folder_names i.e., cut names
for i, foldername in enumerate(folder_names, start=1):
    label = foldername.replace("_General", "")
    h_cutflow.GetXaxis().SetBinLabel(i, label)

# Loop through each folder (i.e. cut)
for j, foldername in enumerate(folder_names, start=1):
    Nevents = 0
    # Loop over each sample (different files)
    for sample in sample_names:
        print "Plotting from sample", sample
        file = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(fileDir,sample),"read")
        tree = file.Get("%s"%(foldername))
        tree.Draw("N_mu>>h_N_mu%s(10,0,10)"%(sample))
        hist = tree.Get("N_mu") # Use a variable with nonzero entries in every cut
        Nevents += hist.GetEntries()
    h_cutflow.SetBinContent(j, Nevents)
    print "Number of events in cut", foldername, "is", Nevents
    print "\n"

# Create canvas
canvas = TCanvas('canvas', 'cutflow plot', W, H)

# Draw cutflow histogram
h_cutflow.SetFillColor(kBlue)
h_cutflow.SetLineColor(kBlack)
h_cutflow.Draw("hist")
h_cutflow.GetYaxis().SetTitle("Events")
h_cutflow.GetXaxis().SetLabelSize(0.03)

# Save cutflow plot
canvas.SaveAs("%s/cutflow.png"%(plotDirectory))

print "Hope that worked!"
print "\n"
