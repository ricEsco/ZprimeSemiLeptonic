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
# parser.add_option("--QCD","--isQCD", dest="isQCD",   default=False, action="store_true", help="Plot as QCD and others" )
# parser.add_option("--DNN","--isDNN", dest="isDNN", default=False,action="store_true", help="Plot to check DNN" )
(options, args) = parser.parse_args()
channel = options.channel
# DNN = options.isDNN
# QCD = options.isQCD
Log=options.isLog
### If reimplementing Log option, place this if-statement in hists for-loop after defining minVal, maxVal ###
# if Log:
#     stack[histName].SetMaximum(15**(1.5*log10(maxVal) - 0.5*log10(minVal)))
### And place non-log equivalent in an else-statement ###
### Then put this just before the end of the hists for-loop to distinguish log files by filename ###
# if Log:
#     canvasRatio.SaveAs("%s/%s_log.png"%(plotDirectory,histName))
### Options stuff --------------------------------------------------------------------------------------------------------------------

### User input section i.e. modify each time new set of plots are being generated ------------------------------------------------------------------------------------
# Data Info - Used for lumi info in CMS_lumi.py function
#Year = 16
#Year = 16
#Year = 17
#Year = 18

# Channel info
if channel=="ele":
    _channelText = "e+jets"
    plotDirectory = ""
    _fileDir = "/nfs/dust/cms/group/zprime-uhh//"
else:
    _channelText = "#mu+jets"
    plotDirectory = "/nfs/dust/cms/user/ricardo/SpinCorrAnalysis_Gen/plots/updatedBoostProcedure/bjet/UL18"
    _fileDir =      "/nfs/dust/cms/user/ricardo/SpinCorrAnalysis_Gen/updatedBoostProcedure/muon/workdir_SpinCorr_Gen_UL18_muon"
print "channel is ", channel
print "The input root files will come from", _fileDir
print "The output will go into", plotDirectory, "\n"


### define the histograms dictionary with entry syntax: {"variable_handle" : ["Plot name", "vertical-axis name", number of bins, [x-min, x-max]]}
if channel=="mu": 
       histograms =  {"pt_hadTop"                      : ["Hadronic-top (both) p_{T}",       "Events", 25, [     0,   500]],
                     "phi_lep_LabFrame"                : ["#phi_{#mu} LabFrame",             "Events", 12, [-np.pi, np.pi]],
                     "phi_lep_CoMFrame"                : ["#phi_{#mu} CoMFrame",             "Events", 12, [-np.pi, np.pi]],
                     "phi_lep_helicityFrame"           : ["#phi_{#mu} helicityFrame",         "Events", 12, [-np.pi, np.pi]],
                     "phi_lep"                         : ["#phi_{#mu}",                      "Events", 12, [-np.pi, np.pi]],
                    #  "phi_b_LabFrame"                  : ["#phi_{b} LabFrame",               "Events", 12, [-np.pi, np.pi]],
                    #  "phi_b_CoMFrame"                  : ["#phi_{b} CoMFrame",               "Events", 12, [-np.pi, np.pi]],
                    #  "phi_b_helicityFrame"             : ["#phi_{b} helicityFrame",           "Events", 12, [-np.pi, np.pi]],
                    #  "phi_b"                           : ["#phi_{b}",                        "Events", 12, [-np.pi, np.pi]],
                     "phi_qlow_LabFrame"               : ["#phi_{q-low} LabFrame",           "Events", 12, [-np.pi, np.pi]],
                     "phi_qlow_CoMFrame"               : ["#phi_{q-low} CoMFrame",           "Events", 12, [-np.pi, np.pi]],
                     "phi_qlow_helicityFrame"          : ["#phi_{q-low} helicityFrame",      "Events", 12, [-np.pi, np.pi]],
                     "phi_qlow"                        : ["#phi_{q-low}",                    "Events", 12, [-np.pi, np.pi]],
                     "sphi"                            : ["#Sigma#phi",                      "Events", 12, [-np.pi, np.pi]],
                     "dphi"                            : ["#Delta#phi",                      "Events", 12, [-np.pi, np.pi]],
                    } # debug

    #   histograms = {"pt_hadTop"                       : ["Hadronic-top (both) p_{T}",       "Events", 25, [     0,   500]],
    #                 "phi_lep"                         : ["#phi_{#mu}",                      "Events", 12, [-np.pi, np.pi]],
    #                 "phi_lep_LabFrame"                : ["#phi_{#mu} LabFrame",             "Events", 12, [-np.pi, np.pi]],
    #                 "phi_lep_CoMFrame"                : ["#phi_{#mu} CoMFrame",             "Events", 12, [-np.pi, np.pi]],
    #                 "phi_lep_topRestFrame"            : ["#phi_{#mu} topRestFrame",         "Events", 12, [-np.pi, np.pi]],
    #                 "phi_lep_high"                    : ["#phi_{#mu}^{high pt}",            "Events", 12, [-np.pi, np.pi]],
    #                 "phi_lep_low"                     : ["#phi_{#mu}_{low pt}",             "Events", 12, [-np.pi, np.pi]],
    #                 "phi_b"                           : ["#phi_{b}",                        "Events", 12, [-np.pi, np.pi]],
    #                 "phi_b_LabFrame"                  : ["#phi_{b} LabFrame",               "Events", 12, [-np.pi, np.pi]],
    #                 "phi_b_CoMFrame"                  : ["#phi_{b} CoMFrame",               "Events", 12, [-np.pi, np.pi]],
    #                 "phi_b_topRestFrame"              : ["#phi_{b} topRestFrame",           "Events", 12, [-np.pi, np.pi]],
    #                 "phi_b_high"                      : ["#phi_{b}^{high pt}",              "Events", 12, [-np.pi, np.pi]],
    #                 "phi_b_low"                       : ["#phi_{b}_{low pt}",               "Events", 12, [-np.pi, np.pi]],
    #                 "sphi"                            : ["#Sigma#phi",                      "Events", 12, [-np.pi, np.pi]],
    #                 "sphi_low"                        : ["#Sigma#phi_{low pt}",             "Events", 12, [-np.pi, np.pi]],
    #                 "sphi_high"                       : ["#Sigma#phi^{high pt}",            "Events", 12, [-np.pi, np.pi]],
    #                 "sphi_plus"                       : ["#Sigma#phi (#mu^{+})",            "Events", 12, [-np.pi, np.pi]],
    #                 "dphi"                            : ["#Delta#phi",                      "Events", 12, [-np.pi, np.pi]],
    #                 "dphi_low"                        : ["#Delta#phi_{low pt}",             "Events", 12, [-np.pi, np.pi]],
    #                 "dphi_high"                       : ["#Delta#phi^{high pt}",            "Events", 12, [-np.pi, np.pi]],
    #                 "phi_lepPlus"                     : ["#phi_{#mu^{+}}",                  "Events", 12, [-np.pi, np.pi]],
    #                 "phi_lepMinus"                    : ["#phi_{#mu^{-}}",                  "Events", 12, [-np.pi, np.pi]],
    #                 "sphi_plus"                       : ["#Sigma#phi (#mu^{+})",            "Events", 12, [-np.pi, np.pi]],
    #                 "sphi_plus_high"                  : ["#Sigma#phi^{high pt} (#mu^{+})",  "Events", 12, [-np.pi, np.pi]],
    #                 "sphi_plus_low"                   : ["#Sigma#phi_{low pt} (#mu^{+})",   "Events", 12, [-np.pi, np.pi]],
    #                 "dphi_plus"                       : ["#Delta#phi (#mu^{+})",            "Events", 12, [-np.pi, np.pi]],
    #                 "dphi_plus_high"                  : ["#Delta#phi^{high pt} (#mu^{+})",  "Events", 12, [-np.pi, np.pi]],
    #                 "dphi_plus_low"                   : ["#Delta#phi_{low pt} (#mu^{+})",   "Events", 12, [-np.pi, np.pi]],
    #                 "sphi_minus"                      : ["#Sigma#phi (#mu^{-})",            "Events", 12, [-np.pi, np.pi]],
    #                 "sphi_minus_high"                 : ["#Sigma#phi^{high pt} (#mu^{-})",  "Events", 12, [-np.pi, np.pi]],
    #                 "sphi_minus_low"                  : ["#Sigma#phi_{low pt} (#mu^{-})",   "Events", 12, [-np.pi, np.pi]],
    #                 "dphi_minus"                      : ["#Delta#phi (#mu^{-})",            "Events", 12, [-np.pi, np.pi]],
    #                 "dphi_minus_high"                 : ["#Delta#phi^{high pt} (#mu^{-})",  "Events", 12, [-np.pi, np.pi]],
    #                 "dphi_minus_low"                  : ["#Delta#phi_{low pt} (#mu^{-})",   "Events", 12, [-np.pi, np.pi]]
    #               } # plots
else:
	histograms = {}


### The sample_names array has filenames of samples that will contribute to each plot
#sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_2", "TTTo2L2Nu", "TTToHadronic"]
sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_2"]


### The stackList dictionary maps the list of samples to their color
# All three TTbar decay channels
#stackList = {"TTToSemiLeptonic":[kRed], "TTToSemiLeptonic_2":[kRed], "TTTo2L2Nu":[kRed+1], "TTToHadronic":[kRed-4]}
stackList = {"TTToSemiLeptonic":[kRed], "TTToSemiLeptonic_2":[kRed]}

### end User input section i.e. modify each time new set of plots are being generated ------------------------------------------------------------------------------------

# Check that user actually paid attention to the above section
switchVar = raw_input("Did you remember to modify the _fileDir variable and create the plotDirectory if it didn't already exist? [y/n] \n")
if switchVar.lower()=="n":
    exit()
elif switchVar.lower()=="y":
    pass
print "\n"

### Histogram Information:
padRatio = 0.25
padOverlap = 0.15
padGap = 0.01
gROOT.SetBatch(True)
# [X-axis title, 
#  Y-axis title,
#  Rebinning factor,
#  [x-min,x-max], -1 means keep as is
#  Extra text about region
#  log plot]
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
CMS_lumi.writeExtraText = True
H = 600;
W = 800;
# references for T, B, L, R                                                                                                             
T = 0.08*H
B = 0.12*H
L = 0.12*W
R = 0.10*W


### Legend info:
#legendHeightPer = 0.04
legendHeightPer = 0.01
legList = stackList.keys() 
#legList.reverse()
#legendStart = 0.69
legendStart = 0.85
legendEnd = 0.99-(R/W)
#legend = TLegend(2*legendStart - legendEnd, 1-T/H-0.01 - legendHeightPer*(len(legList)+1), legendEnd, 0.99-(T/H)-0.01)
# legend = TLegend(2*legendStart - legendEnd , 0.99 - (T/H)/(1.-padRatio+padOverlap) - legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.), legendEnd, 0.99-(T/H)/(1.-padRatio+padOverlap))
legend = TLegend(.85, .70, .95, .90)
legend.SetNColumns(1)
legend.SetHeader("Gen-Level", "C")


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

# canvasRatio = TCanvas('c1Ratio','c1Ratio',W,H)
# canvasRatio.SetFillColor(0)
# canvasRatio.SetBorderMode(0)
# canvasRatio.SetFrameFillStyle(0)
# canvasRatio.SetFrameBorderMode(0)
# canvasRatio.SetLeftMargin( L/W )
# canvasRatio.SetRightMargin( R/W )
# canvasRatio.SetTopMargin( T/H )
# canvasRatio.SetBottomMargin( B/H )
# canvasRatio.SetTickx(0)
# canvasRatio.SetTicky(0)
# canvasRatio.Draw()
# canvasRatio.cd()

# pad1 = TPad("zxc_p1","zxc_p1",0,padRatio-padOverlap,1,1)
pad1 = TPad("zxc_p1","zxc_p1",0,0,1,1)
# pad2 = TPad("qwe_p2","qwe_p2",0,0,1,padRatio+padOverlap)
pad1.SetLeftMargin( L/W )
pad1.SetRightMargin( R/W )
pad1.SetTopMargin( T/H/(1-padRatio+padOverlap) )
pad1.SetBottomMargin( (padOverlap+padGap)/(1-padRatio+padOverlap) )
# pad2.SetLeftMargin( L/W )
# pad2.SetRightMargin( R/W )
# pad2.SetTopMargin( (padOverlap)/(padRatio+padOverlap) )
# pad2.SetBottomMargin( B/H/(padRatio+padOverlap) )
pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetFrameFillStyle(0)
pad1.SetFrameBorderMode(0)
pad1.SetTickx(0)
pad1.SetTicky(0)

# pad2.SetFillColor(0)
# pad2.SetFillStyle(4000)
# pad2.SetBorderMode(0)
# pad2.SetFrameFillStyle(0)
# pad2.SetFrameBorderMode(0)
# pad2.SetTickx(0)
# pad2.SetTicky(0)

# canvasRatio.cd()
pad1.Draw()
# pad2.Draw()
canvas.cd()
canvas.ResetDrawn()


### Initializing container variables
sum_=0
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
    # legendR[histName] = TLegend(2*legendStart - legendEnd , 0.99 - (T/H)/(1.-padRatio+padOverlap) - legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.)-0.1, legendEnd, 0.99-(T/H)/(1.-padRatio+padOverlap))
    legendR[histName] = TLegend(.74, .82, .89, .90)
    legendR[histName].SetHeader("Gen-Level", "C")
    legendR[histName].SetNColumns(1)
    legendR[histName].SetBorderSize(0)
    legendR[histName].SetFillColor(0)
canvas.cd()

# Normalization factor for hists
normFactor = 1.
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
        # hist[histName][sample].Scale(1./hist[histName][sample].Integral())
        if (sample=="TTToSemiLeptonic_2" or sample=="TTToSemiLeptonic_3"):
            hist[histName][sample].SetYTitle(histograms[histName][1])      
            print "Filling with", sample
            stack[histName].Add(hist[histName][sample])
            continue
        legendR[histName].AddEntry(hist[histName][sample],sample,'f')
        hist[histName][sample].SetYTitle(histograms[histName][1])    
        print "Filling with", sample
        stack[histName].Add(hist[histName][sample])   
    # _file["Data"] = TFile("%s/uhh2.AnalysisModuleRunner.DATA.DATA.root"%(_fileDir),"read")
    # tree = _file["Data"].Get("AnalysisTree")
    # tree.Draw("%s>>dat_hist(%i,%f,%f)"%(histName,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]))
    # dataHist=tree.GetHistogram()
    # dataHist.SetMarkerColor(kBlack)
    # dataHist.SetYTitle(histograms[histName][1])
    # dataHist.Scale(10./dataHist.Integral())
    stack[histName].Draw("HIST")
    # print "Filling with Data"
    # dataHist.Draw("pe,x0,same")   
    oneLine = TF1("oneline","1",-9e9,9e9)
    oneLine.SetLineColor(kBlack)
    oneLine.SetLineWidth(1)
    oneLine.SetLineStyle(2)	
    
    # Set vertical-axis limits
    maxVal = stack[histName].GetMaximum()
    minVal = max(stack[histName].GetStack()[0].GetMinimum(), 1)

    # Set vertical-axis limits for "normalized" plots
    # maxVal = normFactor*1.2
    # minVal = normFactor*0.8
    # stack[histName].SetMaximum(maxVal)
    # stack[histName].SetMinimum(0)

    if Log:
        stack[histName].SetMaximum(10**(1.5*log10(maxVal) - 0.5*log10(minVal)))
        stack[histName].SetMinimum(minVal)
    else:
        stack[histName].SetMaximum(1.5*maxVal)
        stack[histName].SetMinimum(minVal)

    # errorband=stack[histName].GetStack().Last().Clone("error")
    # errorband.Sumw2()
    # errorband.SetLineColor(kBlack)
    # errorband.SetFillColor(kBlack)
    # errorband.SetFillStyle(3245)
    # errorband.SetMarkerSize(0)
    
    # canvasRatio.cd()
    # canvasRatio.ResetDrawn()
    # canvasRatio.Draw()
    # canvasRatio.cd()
    
    pad1.Draw()
    # pad2.Draw()
    
    pad1.cd()
    pad1.SetLogy(Log)
    
    y2 = pad1.GetY2()
    
    stack[histName].Draw("HIST")
    # stack[histName].GetXaxis().SetTitle('')
    stack[histName].GetXaxis().SetTitle(histograms[histName][0])
    # stack[histName].GetXaxis().SetRangeUser()
    # stack[histName].GetYaxis().SetTitle(dataHist.GetYaxis().GetTitle())
    stack[histName].SetTitle('')
    stack[histName].GetXaxis().SetLabelSize(0.06)
    stack[histName].GetYaxis().SetLabelSize(gStyle.GetLabelSize()/(1.-padRatio+padOverlap))
    stack[histName].GetYaxis().SetTitleSize(gStyle.GetTitleSize()/(1.-padRatio+padOverlap))
    stack[histName].GetYaxis().SetTitleOffset(gStyle.GetTitleYOffset()*(1.-padRatio+padOverlap))
    stack[histName].GetYaxis().SetTitle("Events")
    # stack[histName].GetYaxis().SetTitle("Normalized Rate")
    
    # dataHist.Draw("E,X0,SAME")
    # legendR[histName].AddEntry(dataHist, "Data", 'pe')
    # ratio = dataHist.Clone("temp")
    # temp = stack[histName].GetStack().Last().Clone("temp")

    # for i_bin in range(1,temp.GetNbinsX()+1):
    #     temp.SetBinError(i_bin,0.)
    # ratio.Divide(temp)
    # #print ratio.GetMaximum(), ratio.GetMinimum()
    # #max_=ratio.GetMaximum()*0.1+ratio.GetMaximum()
    # #min_=ratio.GetMinimum()-0.1*ratio.GetMinimum()

    # ratio.SetTitle('')
    # ratio.GetXaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
    # ratio.GetYaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
    # ratio.GetXaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
    # ratio.GetYaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
    # ratio.GetYaxis().SetTitleOffset(gStyle.GetTitleYOffset()*(padRatio+padOverlap-padGap))
    # ratio.GetYaxis().SetRangeUser(0.8, 1.2)
    # ratio.GetYaxis().SetNdivisions(504)
    # ratio.GetXaxis().SetTitle(histograms[histName][0])
    # ratio.GetYaxis().SetTitle("Data/MC")
    CMS_lumi.CMS_lumi(pad1, 16, 11) # parameters are: (pad, Year, iPosX, extraLumiText = "") <----------- DONT FORGET TO CHANGE THE YEAR TO GET THE RIGHT LUMI
    legendR[histName].Draw()
    # pad2.cd()
    # ratio.SetMarkerStyle(dataHist.GetMarkerStyle())
    # ratio.SetMarkerSize(dataHist.GetMarkerSize())
    # ratio.SetLineColor(dataHist.GetLineColor())
    # ratio.SetLineWidth(dataHist.GetLineWidth())
    # ratio.Draw('e,x0')
    # errorband.Divide(temp)
    # errorband.Draw('e2,same')
    # oneLine.Draw("same")
    # #pad2.Update()
    pad1.Update()
    # canvasRatio.Update()
    # canvasRatio.RedrawAxis()
    # if Log:
    #     canvasRatio.SaveAs("%s/%s_log.png"%(plotDirectory,histName))
    # else:
    #     canvasRatio.SaveAs("%s/%s.png"%(plotDirectory,histName))
    if Log:
        canvas.SaveAs("%s/%s_log.png"%(plotDirectory,histName))
    else:
        canvas.SaveAs("%s/%s.png"%(plotDirectory,histName))
    print "\n"

print "Congratulations, you've successfully generated some Plots."
print "You can find said plots at", plotDirectory
