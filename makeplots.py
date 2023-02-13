#from ROOT import TFile, TLegend, TCanvas, TPad, THStack, TF1, TPaveText, TGaxis, SetOwnership, TObject, gStyle,TH1F
from ROOT import *
import os
import sys
from optparse import OptionParser
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
# Channel info
if channel=="ele":
    _channelText = "e+jets"
    plotDirectory = "analysisPlots_Corr/electron/"
    _fileDir = "/nfs/dust/cms/group/zprime-uhh/Analysis_AzimthCorr_UL18/"
else:
    _channelText = "#mu+jets"
    plotDirectory = "analysisPlots_Corr/muon/helicitybasis_NoChi2"
    _fileDir = "/nfs/dust/cms/group/zprime-uhh/Analysis_AzimthCorr_UL18/muon/workdir_AzimthCorr_UL18_muon"
print "channel is ", channel
print "The input root files will come from", _fileDir
print "The output will go into", plotDirectory, "\n"


### define the histograms dictionary with entry syntax: {"variable_handle" : ["Plot name", "vertical-axis name", number of bins, [x-min, x-max]]}
if channel=="mu": 
    histograms =  {"Mttbar_afterChi2Cut"             : ["M_{ttbar} [GeV]",                       "Events", 25, [   0, 1000]]
                   } # debug

    #  histograms = {"recocount"                       : ["Reco Count",                      "Events",  2, [   1,   3]],
    #                "chi2_afterChi2Cut"               : ["#Chi^{2}",                        "Events", 40, [   0, 200]],
    #                "Mttbar_afterChi2Cut"             : ["M_{ttbar}",                       "Events", 45, [   0, 900]],
    #                "ak4jet1_pt_afterChi2Cut"         : ["AK4_{p_{T}}",                     "Events", 25, [   0, 500]],
    #                "ak4jet1_eta_afterChi2Cut"        : ["AK4_{#eta}",                      "Events", 30, [   0,   3]],
    #                "ak8jet1_pt_afterChi2Cut"         : ["AK8_{p_{T}}",                     "Events", 25, [   0, 500]],
    #                "ak8jet1_eta_afterChi2Cut"        : ["AK8_{#eta}",                      "Events", 30, [   0,   3]],
    #                "pt_hadTop"                       : ["Hadronic Top-Jet p_{T}",          "Events", 25, [   0, 500]],
    #                "deltaR_min"                      : ["#Delta R_{min}",                  "Events", 10, [   0, .25]],
    #                "jets_hadronic_bscore"            : ["b-scores of hadronic jets",       "Events", 20, [   0,   1]],
    #                "bscore_max"                      : ["max b-scores of hadronic jets",   "Events", 20, [   0,   1]],
    #                "phi_lep"                         : ["#phi_{#mu}",                      "Events", 12, [-3.5, 3.5]],
    #                "phi_lep_high"                    : ["#phi_{#mu}^{high pt}",            "Events", 12, [-3.5, 3.5]],
    #                "phi_lep_low"                     : ["#phi_{#mu}_{low pt}",             "Events", 12, [-3.5, 3.5]],
    #                "phi_b"                           : ["#phi_{b}",                        "Events", 12, [-3.5, 3.5]],
    #                "phi_b_high"                      : ["#phi_{b}^{high pt}",              "Events", 12, [-3.5, 3.5]],
    #                "phi_b_low"                       : ["#phi_{b}_{low pt}",               "Events", 12, [-3.5, 3.5]],
    #                "sphi"                            : ["#Sigma#phi",                      "Events", 12, [-3.5, 3.5]],
    #                "sphi_low"                        : ["#Sigma#phi_{low pt}",             "Events", 12, [-3.5, 3.5]],
    #                "sphi_high"                       : ["#Sigma#phi^{high pt}",            "Events", 12, [-3.5, 3.5]],
    #                "sphi_plus"                       : ["#Sigma#phi (#mu^{+})",            "Events", 12, [-3.5, 3.5]],
    #                "dphi"                            : ["#Delta#phi",                      "Events", 12, [-3.5, 3.5]],
    #                "dphi_low"                        : ["#Delta#phi_{low pt}",             "Events", 12, [-3.5, 3.5]],
    #                "dphi_high"                       : ["#Delta#phi^{high pt}",            "Events", 12, [-3.5, 3.5]],
    #                "phi_lepPlus"                     : ["#phi_{#mu^{+}}",                  "Events", 12, [-3.5, 3.5]],
    #                "phi_lepMinus"                    : ["#phi_{#mu^{-}}",                  "Events", 12, [-3.5, 3.5]],
    #                "sphi_plus"                       : ["#Sigma#phi (#mu^{+})",            "Events", 12, [-3.5, 3.5]],
    #                "sphi_plus_high"                  : ["#Sigma#phi^{high pt} (#mu^{+})",  "Events", 12, [-3.5, 3.5]],
    #                "sphi_plus_low"                   : ["#Sigma#phi_{low pt} (#mu^{+})",   "Events", 12, [-3.5, 3.5]],
    #                "dphi_plus"                       : ["#Delta#phi (#mu^{+})",            "Events", 12, [-3.5, 3.5]],
    #                "dphi_plus_high"                  : ["#Delta#phi^{high pt} (#mu^{+})",  "Events", 12, [-3.5, 3.5]],
    #                "dphi_plus_low"                   : ["#Delta#phi_{low pt} (#mu^{+})",   "Events", 12, [-3.5, 3.5]],
    #                "sphi_minus"                      : ["#Sigma#phi (#mu^{-})",            "Events", 12, [-3.5, 3.5]],
    #                "sphi_minus_high"                 : ["#Sigma#phi^{high pt} (#mu^{-})",  "Events", 12, [-3.5, 3.5]],
    #                "sphi_minus_low"                  : ["#Sigma#phi_{low pt} (#mu^{-})",   "Events", 12, [-3.5, 3.5]],
    #                "dphi_minus"                      : ["#Delta#phi (#mu^{-})",            "Events", 12, [-3.5, 3.5]],
    #                "dphi_minus_high"                 : ["#Delta#phi^{high pt} (#mu^{-})",  "Events", 12, [-3.5, 3.5]],
    #                "dphi_minus_low"                  : ["#Delta#phi_{low pt} (#mu^{-})",   "Events", 12, [-3.5, 3.5]]
    #              } # plots

    # histograms = {"jets_hadronic_bscore_after2btag" : ["b-scores of hadronic jets 2btag",      "Events", 20, [   0,   1]],
    #               "bscore_max_2btag"                : ["max b-scores of hadronic jets 2btag",  "Events", 20, [   0,   1]],
    #               "phi_lep_2btag"                   : ["#phi_{#mu} 2btag",                     "Events", 12, [-3.5, 3.5]],
    #               "phi_lep_high_2btag"              : ["#phi_{#mu}^{high pt} 2btag",           "Events", 12, [-3.5, 3.5]],
    #               "pt_hadTop_2btag"                 : ["Hadronic Top-Jet p_{T} (2btag)",       "Events", 25, [   0, 500]],
    #               "phi_lep_low_2btag"               : ["#phi_{#mu}_{low pt} 2btag",            "Events", 12, [-3.5, 3.5]],
    #               "phi_b_2btag"                     : ["#phi_{b} 2btag",                       "Events", 12, [-3.5, 3.5]],
    #               "phi_b_high_2btag"                : ["#phi_{b}^{high pt} 2btag",             "Events", 12, [-3.5, 3.5]],
    #               "phi_b_low_2btag"                 : ["#phi_{b}_{low pt} 2btag",              "Events", 12, [-3.5, 3.5]],
    #               "sphi_2btag"                      : ["#Sigma#phi 2btag",                     "Events", 12, [-3.5, 3.5]],
    #               "sphi_low_2btag"                  : ["#Sigma#phi_{low pt} 2btag",            "Events", 12, [-3.5, 3.5]],
    #               "sphi_high_2btag"                 : ["#Sigma#phi^{high pt} 2btag",           "Events", 12, [-3.5, 3.5]],
    #               "sphi_plus_2btag"                 : ["#Sigma#phi (#mu^{+}) 2btag",           "Events", 12, [-3.5, 3.5]],
    #               "dphi_2btag"                      : ["#Delta#phi 2btag",                     "Events", 12, [-3.5, 3.5]],
    #               "dphi_low_2btag"                  : ["#Delta#phi_{low pt} 2btag",            "Events", 12, [-3.5, 3.5]],
    #               "dphi_high_2btag"                 : ["#Delta#phi^{high pt} 2btag",           "Events", 12, [-3.5, 3.5]],
    #               "phi_lepPlus_2btag"               : ["#phi_{#mu^{+}} 2btag",                 "Events", 12, [-3.5, 3.5]],
    #               "phi_lepMinus_2btag"              : ["#phi_{#mu^{-}} 2btag",                 "Events", 12, [-3.5, 3.5]],
    #               "sphi_plus_2btag"                 : ["#Sigma#phi (#mu^{+}) 2btag",           "Events", 12, [-3.5, 3.5]],
    #               "sphi_plus_high_2btag"            : ["#Sigma#phi^{high pt} (#mu^{+}) 2btag", "Events", 12, [-3.5, 3.5]],
    #               "sphi_plus_low_2btag"             : ["#Sigma#phi_{low pt} (#mu^{+}) 2btag",  "Events", 12, [-3.5, 3.5]],
    #               "dphi_plus_2btag"                 : ["#Delta#phi (#mu^{+}) 2btag",           "Events", 12, [-3.5, 3.5]],
    #               "dphi_plus_high_2btag"            : ["#Delta#phi^{high pt} (#mu^{+}) 2btag", "Events", 12, [-3.5, 3.5]],
    #               "dphi_plus_low_2btag"             : ["#Delta#phi_{low pt} (#mu^{+}) 2btag",  "Events", 12, [-3.5, 3.5]],
    #               "sphi_minus_2btag"                : ["#Sigma#phi (#mu^{-}) 2btag",           "Events", 12, [-3.5, 3.5]],
    #               "sphi_minus_high_2btag"           : ["#Sigma#phi^{high pt} (#mu^{-}) 2btag", "Events", 12, [-3.5, 3.5]],
    #               "sphi_minus_low_2btag"            : ["#Sigma#phi_{low pt} (#mu^{-}) 2btag",  "Events", 12, [-3.5, 3.5]],
    #               "dphi_minus_2btag"                : ["#Delta#phi (#mu^{-}) 2btag",           "Events", 12, [-3.5, 3.5]],
    #               "dphi_minus_high_2btag"           : ["#Delta#phi^{high pt} (#mu^{-}) 2btag", "Events", 12, [-3.5, 3.5]],
    #               "dphi_minus_low_2btag"            : ["#Delta#phi_{low pt} (#mu^{-}) 2btag",  "Events", 12, [-3.5, 3.5]]
    #             } # 2btag plots
else:
	histograms = {"phi_lepTop" : ["ele^{#phi} of Best #Chi^{2} Reconstructed Top", "Events", 12, [-3.5, 3.5]],
                  "phi_hadTop" : ["had^{#phi} of Best #Chi^{2} Reconstructed Top", "Events", 12, [-3.5, 3.5]]
                }

# The stackList dictionary defines maps the list of samples to plot with their color
# All three TTbar decay channels
stackList = {"TTToSemiLeptonic":[kRed], "TTToSemiLeptonic_2":[kRed], "TTToSemiLeptonic_3":[kRed], "TTToSemiLeptonic_4":[kRed], "TTTo2L2Nu":[kRed+1], "TTTo2L2Nu_2":[kRed+1], "TTToHadronic":[kRed-4], "QCD":[kTeal],"WJets":[kGreen], "ST":[kYellow], "DY":[kBlue], "Diboson":[kOrange]}
# Only SemiLeptonic decay channel
# stackList = {"TTToSemiLeptonic_1":[kRed], "TTToSemiLeptonic_2":[kRed], "TTToSemiLeptonic_3":[kRed]}
#print "stackList is ", stackList

# define the sample_names array that has filenames of samples that will contribute to each plot
# All three TTbar decay channels
sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_2", "TTToSemiLeptonic_3", "TTToSemiLeptonic_4", "TTTo2L2Nu", "TTTo2L2Nu_2", "TTToHadronic", "DY", "QCD", "WJets", "ST", "Diboson"]

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
legendStart = 0.69
legendEnd = 0.99-(R/W)
#legend = TLegend(2*legendStart - legendEnd, 1-T/H-0.01 - legendHeightPer*(len(legList)+1), legendEnd, 0.99-(T/H)-0.01)
legend = TLegend(2*legendStart - legendEnd , 0.99 - (T/H)/(1.-padRatio+padOverlap) - legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.), legendEnd, 0.99-(T/H)/(1.-padRatio+padOverlap))
legend.SetNColumns(2)


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

canvasRatio = TCanvas('c1Ratio','c1Ratio',W,H)
canvasRatio.SetFillColor(0)
canvasRatio.SetBorderMode(0)
canvasRatio.SetFrameFillStyle(0)
canvasRatio.SetFrameBorderMode(0)
canvasRatio.SetLeftMargin( L/W )
canvasRatio.SetRightMargin( R/W )
canvasRatio.SetTopMargin( T/H )
canvasRatio.SetBottomMargin( B/H )
canvasRatio.SetTickx(0)
canvasRatio.SetTicky(0)
canvasRatio.Draw()
canvasRatio.cd()

pad1 = TPad("zxc_p1","zxc_p1",0,padRatio-padOverlap,1,1)
pad2 = TPad("qwe_p2","qwe_p2",0,0,1,padRatio+padOverlap)
pad1.SetLeftMargin( L/W )
pad1.SetRightMargin( R/W )
pad1.SetTopMargin( T/H/(1-padRatio+padOverlap) )
pad1.SetBottomMargin( (padOverlap+padGap)/(1-padRatio+padOverlap) )
pad2.SetLeftMargin( L/W )
pad2.SetRightMargin( R/W )
pad2.SetTopMargin( (padOverlap)/(padRatio+padOverlap) )
pad2.SetBottomMargin( B/H/(padRatio+padOverlap) )
pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetFrameFillStyle(0)
pad1.SetFrameBorderMode(0)
pad1.SetTickx(0)
pad1.SetTicky(0)

pad2.SetFillColor(0)
pad2.SetFillStyle(4000)
pad2.SetBorderMode(0)
pad2.SetFrameFillStyle(0)
pad2.SetFrameBorderMode(0)
pad2.SetTickx(0)
pad2.SetTicky(0)

canvasRatio.cd()
pad1.Draw()
pad2.Draw()
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
    legendR[histName] = TLegend(2*legendStart - legendEnd , 0.99 - (T/H)/(1.-padRatio+padOverlap) - legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.)-0.1, legendEnd, 0.99-(T/H)/(1.-padRatio+padOverlap))
    legendR[histName].SetNColumns(2)
    legendR[histName].SetBorderSize(0)
    legendR[histName].SetFillColor(0)
canvas.cd()

# Normalization factor for hists
normFactor = 1.


### We loop through each histogram and fill it with the different samples
for histName in histograms:
    print "--- Working on the", histName, "histogram ---"
    for sample in sample_names:
        _file[sample] = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample),"read")
        tree_MC[histName][sample]=_file[sample].Get("AnalysisTree")
        tree_MC[histName][sample].Draw("%s>>h_%s_%s(%i,%f,%f)"%(histName,histName,sample,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]),"weight*weight_pu")
        hist[histName][sample] = tree_MC[histName][sample].GetHistogram()
        hist[histName][sample].SetFillColor(stackList[sample][0])
        hist[histName][sample].SetLineColor(stackList[sample][0])
        # hist[histName][sample].Scale(1./hist[histName][sample].Integral())
        if (sample=="TTToSemiLeptonic"):
            legendR[histName].AddEntry(hist[histName][sample],"ttbar",'f')
            hist[histName][sample].SetYTitle(histograms[histName][1])      
            print "Filling with", sample
            stack[histName].Add(hist[histName][sample])
            continue
        if (sample=="TTToSemiLeptonic_2" or sample=="TTToSemiLeptonic_3" or sample=="TTToSemiLeptonic_4" or sample=="TTTo2L2Nu_2" or sample=="TTTo2L2Nu" or sample=="TTToHadronic"):
            hist[histName][sample].SetYTitle(histograms[histName][1])      
            print "Filling with", sample
            stack[histName].Add(hist[histName][sample])
            continue
        legendR[histName].AddEntry(hist[histName][sample],sample,'f')
        hist[histName][sample].SetYTitle(histograms[histName][1])    
        print "Filling with", sample
        stack[histName].Add(hist[histName][sample])   
    _file["Data"] = TFile("%s/uhh2.AnalysisModuleRunner.DATA.DATA.root"%(_fileDir),"read")
    tree = _file["Data"].Get("AnalysisTree")
    tree.Draw("%s>>dat_hist(%i,%f,%f)"%(histName,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]))
    dataHist=tree.GetHistogram()
    dataHist.SetMarkerColor(kBlack)
    dataHist.SetYTitle(histograms[histName][1])
    # dataHist.Scale(10./dataHist.Integral())
    stack[histName].Draw("HIST")
    print "Filling with Data"
    dataHist.Draw("pe,x0,same")   
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

    errorband=stack[histName].GetStack().Last().Clone("error")
    errorband.Sumw2()
    errorband.SetLineColor(kBlack)
    errorband.SetFillColor(kBlack)
    errorband.SetFillStyle(3245)
    errorband.SetMarkerSize(0)
    
    canvasRatio.cd()
    canvasRatio.ResetDrawn()
    canvasRatio.Draw()
    canvasRatio.cd()
    
    pad1.Draw()
    pad2.Draw()
    
    pad1.cd()
    pad1.SetLogy(Log)
    
    y2 = pad1.GetY2()
    
    stack[histName].Draw("HIST")
    stack[histName].GetXaxis().SetTitle('')
    stack[histName].GetYaxis().SetTitle(dataHist.GetYaxis().GetTitle())
    stack[histName].SetTitle('')
    stack[histName].GetXaxis().SetLabelSize(0)
    stack[histName].GetYaxis().SetLabelSize(gStyle.GetLabelSize()/(1.-padRatio+padOverlap))
    stack[histName].GetYaxis().SetTitleSize(gStyle.GetTitleSize()/(1.-padRatio+padOverlap))
    stack[histName].GetYaxis().SetTitleOffset(gStyle.GetTitleYOffset()*(1.-padRatio+padOverlap))
    stack[histName].GetYaxis().SetTitle("Events")
    # stack[histName].GetYaxis().SetTitle("Normalized Rate")
    
    dataHist.Draw("E,X0,SAME")
    legendR[histName].AddEntry(dataHist, "Data", 'pe')
    ratio = dataHist.Clone("temp")
    temp = stack[histName].GetStack().Last().Clone("temp")

    for i_bin in range(1,temp.GetNbinsX()+1):
        temp.SetBinError(i_bin,0.)
    ratio.Divide(temp)
    # #print ratio.GetMaximum(), ratio.GetMinimum()
    # #max_=ratio.GetMaximum()*0.1+ratio.GetMaximum()
    # #min_=ratio.GetMinimum()-0.1*ratio.GetMinimum()

    ratio.SetTitle('')
    ratio.GetXaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
    ratio.GetYaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
    ratio.GetXaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
    ratio.GetYaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
    ratio.GetYaxis().SetTitleOffset(gStyle.GetTitleYOffset()*(padRatio+padOverlap-padGap))
    ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio.GetYaxis().SetNdivisions(504)
    ratio.GetXaxis().SetTitle(histograms[histName][0])
    ratio.GetYaxis().SetTitle("Data/MC")
    CMS_lumi.CMS_lumi(pad1, 4, 11)
    legendR[histName].Draw()
    pad2.cd()
    ratio.SetMarkerStyle(dataHist.GetMarkerStyle())
    ratio.SetMarkerSize(dataHist.GetMarkerSize())
    ratio.SetLineColor(dataHist.GetLineColor())
    ratio.SetLineWidth(dataHist.GetLineWidth())
    ratio.Draw('e,x0')
    errorband.Divide(temp)
    errorband.Draw('e2,same')
    oneLine.Draw("same")
    # #pad2.Update()
    pad1.Update()
    canvasRatio.Update()
    canvasRatio.RedrawAxis()
    if Log:
        canvasRatio.SaveAs("%s/%s_log.png"%(plotDirectory,histName))
    else:
        canvasRatio.SaveAs("%s/%s.png"%(plotDirectory,histName))
    print "\n"

print "Congratulations, you've successfully generated some Plots."
print "You can find said plots at", plotDirectory
