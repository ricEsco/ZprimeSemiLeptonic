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

### Options stuff 
parser = OptionParser()
parser.add_option("-c", "--channel", dest="channel", default="mu", type='str', help="Specify which channel mu or ele? default is mu")
parser.add_option("--Log","--isLog", dest="isLog",   default=False, action="store_true", help="Plot the plots in log ?" )
(options, args) = parser.parse_args()
channel = options.channel
Log=options.isLog

### User input section i.e. modify each time new set of plots are being generated
# Channel info
if channel=="ele":
    _channelText = "e+jets"
    plotDirectory = ""
    _fileDir = "/nfs/dust/cms/group/zprime-uhh//"
else:
    _channelText = "#mu+jets"
    plotDirectory = "/nfs/dust/cms/user/ricardo/SpinCorrAnalysis_Gen/plots/BernreutherBasisProcedure/UL18/ogBhhad"
    _fileDir =      "/nfs/dust/cms/user/ricardo/SpinCorrAnalysis_Gen/BernreutherBasisProcedure/ogBhad/muon/workdir_SpinCorr_Gen_UL18_muon/"

print "channel is ", channel
print "The input root files will come from", _fileDir
print "The output will go into", plotDirectory, "\n"


### define the histograms dictionary with entry syntax: {"variable_handle" : ["Plot name", "vertical-axis name", number of bins, [x-min, x-max]]}
if channel=="mu": 
       histograms =  {#"pt_hadTop"              : ["Hadronic-top p_{T}",                       "Events", 25, [     0,   500]],
                    #  "ttbar_mass_LabFrame"    : ["Mass_{t#bar{t}}",                          "Events", 40, [     0,  2000]],
                    #  "ttbar_boost_LabFrame"    : ["longitudinal boost_{t#bar{t}}",            "Events", 10, [     0,     1]],
                    #  "phi_lep"                 : ["#phi_{#mu}",                               "Events", 12, [-np.pi, np.pi]],
                    #  "phi_b"                   : ["#phi_{b-quark}",                           "Events", 12, [-np.pi, np.pi]],
                    #  "phi_qlow"                : ["#phi_{W-quark_{soft}}",                    "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lq"                 : ["#Sigma#phi (l,W_{q})",                     "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lq_low"             : ["#Sigma#phi_{low pt} (l,W_{q})",            "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lq_high"            : ["#Sigma#phi^{high pt} (l,W_{q})",           "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lq_Mass1"           : ["#Sigma#phi_{Mass 0-500} (l,W_{q})",        "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lq_Mass2"           : ["#Sigma#phi_{Mass 500-750} (l,W_{q})",      "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lq_Mass3"           : ["#Sigma#phi_{Mass 750-1000} (l,W_{q})",     "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lq_Mass4"           : ["#Sigma#phi_{Mass 1000-1500} (l,W_{q})",    "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lq_Mass5"           : ["#Sigma#phi_{Mass 1500-Inf} (l,W_{q})",     "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lq"                 : ["#Delta#phi (l,W_{q})",                     "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lq_low"             : ["#Delta#phi_{low pt} (l,W_{q})",            "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lq_high"            : ["#Delta#phi^{high pt} (l,W_{q})",           "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lq_Mass1"           : ["#Delta#phi_{Mass 0-500} (l,W_{q})",        "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lq_Mass2"           : ["#Delta#phi_{Mass 500-750} (l,W_{q})",      "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lq_Mass3"           : ["#Delta#phi_{Mass 750-1000} (l,W_{q})",     "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lq_Mass4"           : ["#Delta#phi_{Mass 1000-1500} (l,W_{q})",    "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lq_Mass5"           : ["#Delta#phi_{Mass 1500-Inf} (l,W_{q})",     "Events", 12, [-np.pi, np.pi]],
                     "sphi_lb"                 : ["#Sigma#phi (l,b)",                         "Events", 12, [-np.pi, np.pi]],
                     "sphi_lb_low"             : ["#Sigma#phi_{low pt} (l,b)",                "Events", 12, [-np.pi, np.pi]],
                     "sphi_lb_high"            : ["#Sigma#phi^{high pt} (l,b)",               "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lb_Mass1"           : ["#Sigma#phi_{Mass 0-500} (l,b)",            "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lb_Mass2"           : ["#Sigma#phi_{Mass 500-750} (l,b)",          "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lb_Mass3"           : ["#Sigma#phi_{Mass 750-1000} (l,b)",         "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lb_Mass4"           : ["#Sigma#phi_{Mass 1000-1500} (l,b)",        "Events", 12, [-np.pi, np.pi]],
                    #  "sphi_lb_Mass5"           : ["#Sigma#phi_{Mass 1500-Inf} (l,b)",         "Events", 12, [-np.pi, np.pi]],
                     "dphi_lb"                 : ["#Delta#phi (l,b)",                         "Events", 12, [-np.pi, np.pi]],
                     "dphi_lb_low"             : ["#Delta#phi_{low pt} (l,b)",                "Events", 12, [-np.pi, np.pi]],
                     "dphi_lb_high"            : ["#Delta#phi^{high pt} (l,b)",               "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lb_Mass1"           : ["#Delta#phi_{Mass 0-500} (l,b)",            "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lb_Mass2"           : ["#Delta#phi_{Mass 500-750} (l,b)",          "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lb_Mass3"           : ["#Delta#phi_{Mass 750-1000} (l,b)",         "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lb_Mass4"           : ["#Delta#phi_{Mass 1000-1500} (l,b)",        "Events", 12, [-np.pi, np.pi]],
                    #  "dphi_lb_Mass5"           : ["#Delta#phi_{Mass 1500-Inf} (l,b)",         "Events", 12, [-np.pi, np.pi]]
                    }
else:
	histograms = {}


### The sample_names array has filenames of samples that will contribute to each plot
#sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_2", "TTToSemiLeptonic_3", "TTToSemiLeptonic_4", "TTToSemiLeptonic_5", "TTToSemiLeptonic_6", "TTToSemiLeptonic_7", "TTToSemiLeptonic_8"]
sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_8"]

### The stackList dictionary maps the list of samples to their color
stackList = {"TTToSemiLeptonic":[kRed], "TTToSemiLeptonic_2":[kRed], "TTToSemiLeptonic_3":[kRed], "TTToSemiLeptonic_4":[kRed], "TTToSemiLeptonic_5":[kRed], "TTToSemiLeptonic_6":[kRed], "TTToSemiLeptonic_7":[kRed], "TTToSemiLeptonic_8":[kRed]}

### end User input section i.e. modify each time new set of plots are being generated ------------------------------------------------------------------------------------

### Histogram Information:
padRatio = 0.25
padOverlap = 0.15
padGap = 0.01
gROOT.SetBatch(True)
regionText ="loose selection"
thestyle = Style()
HasCMSStyle = False
style = None
if os.path.isfile('tdrstyle.C'):
    ROOT.gROOT.ProcessLine('.L tdrstyle.C')
    ROOT.setTDRStyle()
    #print "Found tdrstyle.C file, using this style."
    HasCMSStyle = True
    if os.path.isfile('CMSTopStyle.cc'):
        gROOT.ProcessLine('.L CMSTopStyle.cc+')
        style = CMSTopStyle()
        style.setupICHEPv1()
        #print "Found CMSTopStyle.cc file, use TOP style if requested in xml file."
if not HasCMSStyle:
    #print "Using default style defined in cuy package."
    thestyle.SetStyle()
ROOT.gROOT.ForceStyle()
CMS_lumi.channelText = _channelText
CMS_lumi.writeChannelText = True
CMS_lumi.writeExtraText = False
H = 600;
W = 800;
# references for T, B, L, R                                                                                                             
T = 0.08*H
B = 0.12*H
L = 0.12*W
R = 0.10*W

### Canvas info:
canvas = TCanvas('c1','c1',W,H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)

pad1 = TPad("zxc_p1","zxc_p1",0,0,1,1)
pad1.SetLeftMargin( L/W )
pad1.SetRightMargin( R/W )
pad1.SetTopMargin( T/H/(1-padRatio+padOverlap) )
pad1.SetBottomMargin( (padOverlap+padGap)/(1-padRatio+padOverlap) )
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
    ### Legend info:
    legendR[histName] = TLegend(.59, .80, .89, .90)
    legendR[histName].SetNColumns(1)
    legendR[histName].SetBorderSize(0)
    legendR[histName].SetFillColor(0)
pad1.cd()

### We loop through each histogram and fill it with the different samples
for histName in histograms:
    pad1.Clear()
    pad1.Update()
    print "--- Working on the", histName, "histogram ---"
    for sample in sample_names:
        #print "sample is", sample
        _file[sample] = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample),"read")
        tree_MC[histName][sample]=_file[sample].Get("AnalysisTree")
        tree_MC[histName][sample].Draw("%s>>h_%s_%s(%i,%f,%f)"%(histName,histName,sample,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]))
        hist[histName][sample] = tree_MC[histName][sample].GetHistogram()
        # hide lines from the first histograms
        # hist[histName][sample].SetLineColor(0)
        # hist[histName][sample].SetMarkerColor(0)
        hist[histName][sample].GetXaxis().SetRangeUser(-np.pi, np.pi)

        # hist[histName][sample].SetLineColor(0)         # making the first histogram invisible
        # hist[histName][sample].SetMarkerColor(kWhite)  # making the first histogram invisible
        if (sample=="TTToSemiLeptonic_8"):
            # hist[histName][sample].SetLineColor(kRed)
            # hist[histName][sample].SetMarkerColor(kBlack)
            # hist[histName][sample].SetMarkerStyle(1)
            # hist[histName][sample].SetLineWidth(2)
            # hist[histName][sample].SetYTitle(histograms[histName][1])
            hist[histName][sample].GetXaxis().SetRangeUser(-np.pi, np.pi)
            # legendR[histName].AddEntry(hist[histName][sample],"UL 18 TTbar-semi",'le')
            # legendR[histName].SetTextSize(0.04)
            print "Filling with", sample
            stack[histName].Add(hist[histName][sample])
            #print "Histogram content for", sample, ":", [hist[histName][sample].GetBinContent(i) for i in range(1, hist[histName][sample].GetNbinsX() + 1)]  # Debug print
            continue

        print "Filling with", sample
        stack[histName].Add(hist[histName][sample])
        #print "Histogram content for", sample, ":", [hist[histName][sample].GetBinContent(i) for i in range(1, hist[histName][sample].GetNbinsX() + 1)]  # Debug print
        
    # create a histogram copy of the THStack
    h_sum = stack[histName].GetStack().Last().Clone("h_sum")
    # normalize to center around amplitude
    half_range = (h_sum.GetMaximum() - h_sum.GetMinimum()) / 2
    norm = h_sum.GetMaximum() - half_range
    h_sum.Scale(1/norm)
    max_bin = h_sum.GetMaximum()
    min_bin = h_sum.GetMinimum()
    # set vertical axis limits to center histogram
    h_sum.SetMaximum(max_bin*1.2)
    h_sum.SetMinimum(min_bin*0.8)

    # clear pad of previous histograms
    pad1.Clear()

    # draw histogram with bells and whistles
    # legendR[histName].Draw()
    h_sum.SetLineColor(kRed)
    h_sum.SetMarkerColor(kBlack)
    h_sum.SetMarkerStyle(1)
    h_sum.SetLineWidth(2)
    h_sum.GetXaxis().SetTitle(histograms[histName][0])
    h_sum.GetXaxis().SetLabelSize(0.06)
    h_sum.SetYTitle(histograms[histName][1])
    h_sum.GetYaxis().SetTitle("Events")
    h_sum.Draw("hist same")

    # CMS_lumi.CMS_lumi(pad1, 18, 11) # parameters are: (pad, Year, iPosX, extraLumiText = "") <----------- DONT FORGET TO CHANGE THE YEAR TO GET THE RIGHT LUMI

    pad1.Update()
    if Log:
        canvas.SaveAs("%s/%s_log.png"%(plotDirectory,histName))
    else:
        canvas.SaveAs("%s/%s.png"%(plotDirectory,histName))
    print "\n"

print "Congratulations, you've successfully generated some Plots."
print "You can find said plots at", plotDirectory
