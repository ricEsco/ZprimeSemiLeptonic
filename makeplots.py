#from ROOT import TFile, TLegend, TCanvas, TPad, THStack, TF1, TPaveText, TGaxis, SetOwnership, TObject, gStyle,TH1F
from ROOT import *
import os
import sys
from optparse import OptionParser
import numpy as np
from numpy import log10
from array import array
import CMS_lumi
from Style import *

### Options stuff --------------------------------------------------------------------------------------------------------------------
parser = OptionParser()
parser.add_option("-c", "--channel", dest="channel", default="mu", type='str', help="Specify which channel mu or ele? default is mu")
parser.add_option("--Log","--isLog", dest="isLog",   default=False, action="store_true", help="Plot the plots in log ?" )
(options, args) = parser.parse_args()
channel = options.channel
Log=options.isLog
### Options stuff --------------------------------------------------------------------------------------------------------------------


# Channel info
if channel=="ele":
    _channelText = "e+jets"
    plotDirectory = ""
    _fileDir = "/nfs/dust/cms/group/zprime-uhh//"
else:
    _channelText = "#mu+jets"
    plotDirectory = "/nfs/dust/cms/user/ricardo/AzCorrAnalysis/muon/plots/UL18"
    _fileDir = "/nfs/dust/cms/user/ricardo/AzCorrAnalysis/muon/workdir_AzCorr_UL18_muon"
print "channel is ", channel
print "Input root:", _fileDir
print "Output directory:", plotDirectory, "\n"


### define the histograms dictionary with entry syntax: {"variable_handle" : ["Plot name", "vertical-axis name", number of bins, [x-min, x-max]]}
if channel=="mu": 
     histograms = {"pt_hadTop_res"                   : ["p_{T} of resolved top_{Hadronic}",  "Events", 25, [     0,   500]],
                   "pt_hadTop_mer"                   : ["p_{T} of merged top_{Hadronic}",  "Events", 29, [   350,   950]],
                   "pt_hadTop"                       : ["Hadronic-top (both) p_{T}",       "Events", 25, [     0,   500]],
                   "res_jet_bscore"                  : ["b-scores of Resolved-top jets",   "Events", 20, [     0,     1]],
                   "mer_subjet_bscore"               : ["b-scores of Merged-top subjets",  "Events", 20, [     0,     1]],
                   "bscore_max_beforecut"            : ["max b-scores before WP cut",      "Events", 20, [     0,     1]],
                   "bscore_max"                      : ["max b-scores of hadronic jets",   "Events", 20, [     0,     1]],
                   "sphi"                            : ["#Sigma#phi",                      "Events", 12, [-np.pi, np.pi]],
                   "sphi_low"                        : ["#Sigma#phi_{low pt}",             "Events", 12, [-np.pi, np.pi]],
                   "sphi_high"                       : ["#Sigma#phi^{high pt}",            "Events", 12, [-np.pi, np.pi]],
                   "dphi"                            : ["#Delta#phi",                      "Events", 12, [-np.pi, np.pi]],
                   "dphi_low"                        : ["#Delta#phi_{low pt}",             "Events", 12, [-np.pi, np.pi]],
                   "dphi_high"                       : ["#Delta#phi^{high pt}",            "Events", 12, [-np.pi, np.pi]],
                 } # plots
else:
	histograms = {}


### The sample_names array has filenames of samples that will contribute to each plot
sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_2", "TTTo2L2Nu", "TTToHadronic"]


### The stackList dictionary maps the list of samples to their color
# All three TTbar decay channels
stackList = {"TTToSemiLeptonic":[kRed], "TTToSemiLeptonic_2":[kRed], "TTTo2L2Nu":[kRed+2], "TTToHadronic":[kRed-4]}

### end User input section i.e. modify each time new set of plots are being generated ------------------------------------------------------------------------------------

# Check that user actually paid attention to the above section
switchVar = raw_input("Did you remember to modify the _fileDir variable and create the plotDirectory if it didn't already exist? [y/n] \n")
if switchVar.lower()=="n":
    exit()
elif switchVar.lower()=="y":
    pass
print "\n"

### Histogram Information:
gROOT.SetBatch(True)
regionText ="loose selection"
thestyle = Style()
HasCMSStyle = False
style = None
if os.path.isfile('tdrstyle.C'):
    ROOT.gROOT.ProcessLine('.L tdrstyle.C')
    ROOT.setTDRStyle()
    print "Found tdrstyle.C file, using this style."
    HasCMSStyle = True
    if os.path.isfile('CMSTopStyle.cc'):
        gROOT.ProcessLine('.L CMSTopStyle.cc+')
        style = CMSTopStyle()
        style.setupICHEPv1()
        print "Found CMSTopStyle.cc file, use TOP style if requested in xml file."
if not HasCMSStyle:
    print "Using default style defined in cuy package."
    thestyle.SetStyle()
ROOT.gROOT.ForceStyle()
CMS_lumi.channelText = _channelText
CMS_lumi.writeChannelText = True
CMS_lumi.writeExtraText = True
H = 600;
W = 800;

### Legend info:
legendHeightPer = 0.01
legList = stackList.keys() 
legend = TLegend(.59, .80, .89, .90)

### Canvas info:
canvas = TCanvas('c1','c1',W,H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetTickx(0)

pad1 = TPad("zxc_p1","zxc_p1",0,0,1,1)
pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetFrameFillStyle(0)
pad1.SetFrameBorderMode(0)
pad1.SetTickx(0)
pad1.SetTicky(0)
pad1.Draw()

canvas.cd()
canvas.ResetDrawn()


### Initializing container variables
tree_MC={}
hist={}
_file={}
stack={}
legendR={}


### Initialize ROOT histogram stuff -I don't actually know what's happening here, see docs and cross your fingers
for histName in histograms:
    tree_MC[histName]={}
    hist[histName]={}
    stack[histName] = THStack("hs","stack")
    legendR[histName] = TLegend(0.7,0.7,0.9,0.9)
    legendR[histName].SetNColumns(2)
    legendR[histName].SetBorderSize(0)
    legendR[histName].SetFillColor(0)
pad1.cd()


### We loop through each histogram and fill it with the different samples
for histName in histograms:
    print "--- Working on the", histName, "histogram ---"
    for sample in sample_names:
        #print "sample is", sample
        _file[sample] = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample),"read")
        tree_MC[histName][sample]=_file[sample].Get("AnalysisTree")
        tree_MC[histName][sample].Draw("%s>>h_%s_%s(%i,%f,%f)"%(histName,histName,sample,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]),"weight*weight_pu")
        hist[histName][sample] = tree_MC[histName][sample].GetHistogram()
        hist[histName][sample].SetFillColor(stackList[sample][0])
        hist[histName][sample].SetLineColor(stackList[sample][0])
        legendR[histName].AddEntry(hist[histName][sample],sample,'f')
        hist[histName][sample].SetYTitle(histograms[histName][1])    
        print "Filling with", sample
        stack[histName].Add(hist[histName][sample])   
    stack[histName].Draw("HIST")
    
    # Set vertical-axis limits
    maxVal = stack[histName].GetMaximum()
    minVal = max(stack[histName].GetStack()[0].GetMinimum(), 1)

    if Log:
        stack[histName].SetMaximum(10**(1.5*log10(maxVal) - 0.5*log10(minVal)))
        stack[histName].SetMinimum(minVal)
    else:
        stack[histName].SetMaximum(1.5*maxVal)
        stack[histName].SetMinimum(minVal)
    
    pad1.cd()
    pad1.SetLogy(Log)
    
    stack[histName].Draw("HIST")
    stack[histName].SetTitle('')
    stack[histName].GetXaxis().SetTitle('')
    stack[histName].GetXaxis().SetTitle(histograms[histName][0])
    stack[histName].GetXaxis().SetLabelSize(0.06)
    stack[histName].GetYaxis().SetTitle("Events")
    
    CMS_lumi.CMS_lumi(pad1, 18, 11)

    pad1.Update()

    if Log:
        canvas.SaveAs("%s/%s_log.png"%(plotDirectory,histName))
    else:
        canvas.SaveAs("%s/%s.png"%(plotDirectory,histName))
    print "\n"

print "Congratulations, you've successfully generated some Plots."
print "You can find said plots at", plotDirectory
