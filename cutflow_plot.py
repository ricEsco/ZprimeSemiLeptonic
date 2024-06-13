#from ROOT import TFile, TLegend, TCanvas, TPad, THStack, TF1, TPaveText, TGaxis, SetOwnership, TObject, gStyle,TH1F
from ROOT import *
import os
import sys
import numpy as np
from numpy import log10
from array import array
import CMS_lumi
from Style import *

fileDir =      "/nfs/dust/cms/user/ricardo/AzCorrAnalysis/missing_muons/muon/workdir_Zprime_AnalysisDNN_UL18_muon"
plotDirectory = "/nfs/dust/cms/user/ricardo/analysis_output/UL18_DNN/missingmuons"
logplotDirectory = "/nfs/dust/cms/user/ricardo/analysis_output/UL18_DNN/missingmuons/log"

# Choose whether to plot with log scale or not
log = True
# log = False

print "The input root files will come from", fileDir
if log:
    print "The output will go into", logplotDirectory, "\n"
print "The output will go into", plotDirectory, "\n"

### This array contains the names of the folders in the root file
folder_names = ["Weights_Init_General", "Weights_HEM_General", "Weights_PU_General", "Weights_Lumi_General", "Weights_Prefiring_General", "Weights_PS_General", "Weights_TopTag_SF_General", "IsoMuon_SF_General", "IdMuon_SF_General", "MuonReco_SF_General", "TriggerMuon_SF_General", "BeforeBtagSF_General", "AfterBtagSF_General", "AfterCustomBtagSF_General", "NLOCorrections_General", "NNInputsBeforeReweight_General", "TopTagVeto_General"]

### This array has filenames of samples that will contribute to each plot
# sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_2"]
# color = kRed

# sample_names = ["ST"]
# color = kBlue

sample_names = ["WJets"]
color = kGreen

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
    N_mu = 0
    # Loop over each sample (different files)
    for sample in sample_names:
        print "Plotting from sample", sample
        file = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(fileDir,sample),"read")
        tree = file.Get("%s"%(foldername))
        tree.Draw("N_mu>>h_N_mu%s(10,0,10)"%(sample))
        hist = tree.Get("N_mu") # Use a variable with nonzero entries in every cut
        N_mu += hist.Integral()
    h_cutflow.SetBinContent(j, N_mu)
    print "Number of muons in cut:", foldername, "is", N_mu
    print "\n"

# Create canvas
canvas = TCanvas('canvas', 'cutflow plot', W, H)

# Draw cutflow histogram
h_cutflow.SetFillColor(color)
h_cutflow.SetLineColor(kBlack)
h_cutflow.Draw("hist")
h_cutflow.GetYaxis().SetTitle("N_{#mu}")
h_cutflow.GetXaxis().SetLabelSize(0.03)

# Save cutflow plot
if log:
    h_cutflow.SetMinimum(1)
    canvas.SetLogy()
    canvas.SaveAs("%s/%s_cutflow_log.png"%(logplotDirectory, sample_names[0]))
else:
    canvas.SaveAs("%s/%s_cutflow.png"%(plotDirectory, sample_names[0]))

print "Hope that worked!"
print "\n"