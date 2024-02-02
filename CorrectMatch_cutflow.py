#from ROOT import TFile, TLegend, TCanvas, TPad, THStack, TF1, TPaveText, TGaxis, SetOwnership, TObject, gStyle,TH1F
from ROOT import *
import os
import sys
import numpy as np
from numpy import log10
from array import array
import CMS_lumi
from Style import *

fileDir =      "/nfs/dust/cms/user/ricardo/analysis_output/RecoOpt/muon/workdir3_Analysis_UL18_muon_RecOpt_signal"
plotDirectory = "/nfs/dust/cms/user/ricardo/analysis_output/RecoOpt/plots/UL18/muon"

print "---"
print "Input directory:", fileDir
print "Output directory:", plotDirectory
print "---"

### This array contains the names of the CorrectMatch criteria for events with resolved topology
AK4matching_criteria = ["CMfm_notSemilepGen", "CMfm_not1LepJet", "CMfm_notOkHadJetMult", "CMfm_notMatchedLepB", "CMfm_notMatchedHadB_Ak4", "CMfm_notMatchedQ1_Ak4", "CMfm_notMatchedQ2_Ak4", "CMfm_not1to1MatchToJet", "CMfm_notMatchedNeutrino", "CMfm_notMatchedLepton"]

### This array contains the names of the CorrectMatch criteria for events with merged topology
AK8matching_criteria = ["CMfm_notSemilepGen", "CMfm_not1LepJet", "CMfm_notOkHadJetMult", "CMfm_notMatchedLepB", "CMfm_notMatchedHadB_Ak8", "CMfm_notMatchedQ1_Ak8", "CMfm_notMatchedQ2_Ak8", "CMfm_notMatchedNeutrino", "CMfm_notMatchedLepton"]

### This array has filenames of samples that will contribute to each plot
samples = ["TTToSemiLeptonic", "TTToSemiLeptonic_2"]

### end User input section i.e. modify each time new set of plots are being generated ------------------------------------------------------------------------------------

### Histogram Information:
gStyle.SetOptStat(0)
H = 600
W = 800

## Create final histograms-------------------------------------
# Create cutflow histogram: TH*F(name, title, nbins, xlow, xup)
h_AK4cutflow = TH1F("h_AK4cutflow", "CorrectMatch Cutflow (Resolved)", len(AK4matching_criteria), 0, len(AK4matching_criteria))

# Label bins of cutflow histogram in order of matching_criteria
for i, criteria in enumerate(AK4matching_criteria, start=1):
    label = criteria.replace("CMfm_not", "no")
    h_AK4cutflow.GetXaxis().SetBinLabel(i, label)


# Create cutflow histogram: TH*F(name, title, nbins, xlow, xup)
h_AK8cutflow = TH1F("h_AK8cutflow", "CorrectMatch Cutflow (Merged)", len(AK8matching_criteria), 0, len(AK8matching_criteria))

# Label bins of cutflow histogram in order of matching_criteria
for i, criteria in enumerate(AK8matching_criteria, start=1):
    label = criteria.replace("CMfm_not", "no")
    h_AK8cutflow.GetXaxis().SetBinLabel(i, label)
## Create final histograms-------------------------------------

## Fill Resolved histogram-------------------------------------
# Loop through each criteria
for j, criteria in enumerate(AK4matching_criteria, start=1):
    Ncandidates = 0
    # Loop over each sample (different files)
    for sample in samples:
        print "Plotting", criteria, "from", sample
        file = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(fileDir,sample),"read")
        tree = file.Get("AnalysisTree")
        tree.Draw("%s>>h_%s(1,0,1)"%(criteria, criteria))
        hist = tree.GetHistogram()
        Ncandidates += hist.Integral()
    h_AK4cutflow.SetBinContent(j, Ncandidates)
    print "Number of candidates to fail with", criteria, "is:", Ncandidates
    print "\n"

# Create canvas
canvas = TCanvas('canvas', 'cutflow plot', W, H)

# Draw cutflow histogram
h_AK4cutflow.SetFillColor(kBlue)
h_AK4cutflow.SetLineColor(kBlack)
h_AK4cutflow.Draw("hist")
h_AK4cutflow.GetYaxis().SetTitle("Events")
h_AK4cutflow.GetXaxis().SetLabelSize(0.03)

# Save cutflow plot
canvas.SaveAs("%s/AK4cutflow.png"%(plotDirectory))
## Fill Resolved histogram-------------------------------------



## Fill Merged histogram-------------------------------------
# Loop through each criteria
for j, criteria in enumerate(AK8matching_criteria, start=1):
    Ncandidates = 0
    # Loop over each sample (different files)
    for sample in samples:
        print "Plotting", criteria, "from", sample
        file = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(fileDir,sample),"read")
        tree = file.Get("AnalysisTree")
        tree.Draw("%s>>h_%s(1,0,1)"%(criteria, criteria))
        hist = tree.GetHistogram()
        Ncandidates += hist.Integral()
    h_AK8cutflow.SetBinContent(j, Ncandidates)
    print "Number of candidates to fail with", criteria, "is:", Ncandidates
    print "\n"

# Create canvas
canvas = TCanvas('canvas', 'cutflow plot', W, H)

# Draw cutflow histogram
h_AK8cutflow.SetFillColor(kBlue)
h_AK8cutflow.SetLineColor(kBlack)
h_AK8cutflow.Draw("hist")
h_AK8cutflow.GetYaxis().SetTitle("Events")
h_AK8cutflow.GetXaxis().SetLabelSize(0.03)

# Save cutflow plot
canvas.SaveAs("%s/AK8cutflow.png"%(plotDirectory))
## Fill Merged histogram-------------------------------------


print "Hope that worked!"
print "\n"
