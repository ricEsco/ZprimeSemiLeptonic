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
# log = True
log = False

print "The input root files will come from", fileDir
if log:
    print "The output will go into", logplotDirectory, "\n"
print "The output will go into", plotDirectory, "\n"

### This array contains the names of the folders in the root file
# folder_names = ["Weights_Init_General", "Weights_HEM_General", "Weights_PU_General", "Weights_Lumi_General", "Weights_Prefiring_General", "Weights_PS_General", "Weights_TopTag_SF_General", "IsoMuon_SF_General", "IdMuon_SF_General", "MuonReco_SF_General", "TriggerMuon_SF_General", "BeforeBtagSF_General", "AfterBtagSF_General", "AfterCustomBtagSF_General", "NLOCorrections_General", "NNInputsBeforeReweight_General", "TopTagVeto_General"]

### Define where in the analysis process the plot is being made
foldername = "NNInputsBeforeReweight_General"
label = foldername.replace("_General", "")

### This array has filenames of samples that will contribute to each plot
sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_2"]
color = kRed

# sample_names = ["ST"]
# color = kBlue

# sample_names = ["WJets"]
# color = kGreen

### name of variable to plot
variable = "charge_mu"
# variable = "toplep_lepcharge"

### end User input section i.e. modify each time new set of plots are being generated ------------------------------------------------------------------------------------

### Histogram Information:
gStyle.SetOptStat(0)
H = 600
W = 800

# Loop over each sample (in case there are multiple files)
for sample in sample_names:
    print "Plotting from sample", sample
    file = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(fileDir,sample),"read")
    tree = file.Get("%s"%(foldername))
    # tree.Draw("charge_mu>>h_charge_mu%s(10,0,10)"%(sample))
    tree.Draw(variable)
    hist = tree.Get(variable)
    # change bin min and max to -2 and 2
    hist.SetBins(4,-2,2)
# Create canvas
canvas = TCanvas('canvas', "canvas", W, H)

# Histogram details
hist.SetTitle(variable+" at "+label)
hist.SetFillColor(color)
hist.SetLineColor(kBlack)
hist.Draw("hist")
hist.GetXaxis().SetTitle("Muon charge")
hist.GetXaxis().SetLabelSize(0.035)
hist.SetMinimum(0)

# Save plot
if log:
    hist.SetMinimum(1)
    canvas.SetLogy()
    canvas.SaveAs("%s/%s_%s_at_%s_log.png"%(logplotDirectory, sample_names[0], variable, label))
else:
    canvas.SaveAs("%s/%s_%s_at_%s.png"%(plotDirectory, sample_names[0], variable, label))
    

print "Hope that worked!"
print "\n"