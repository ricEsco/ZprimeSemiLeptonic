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
    plotDirectory = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/electron/plots/"
    _fileDir = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/electron/workdir_Zprime_Analysis_UL18_electron_entanglement/"
else:
    _channelText = "#mu+jets"
    plotDirectory = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/muon/plots/"
    _fileDir = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/muon"
print "channel is ", channel
print "Input directory:", _fileDir
print "Output directory:", plotDirectory, "\n"


### define the histograms dictionary with entry syntax: {"variable_handle" : ["x-axis name", "vertical-axis name", number of bins, [x-min, x-max]]}
if channel=="mu": 
    histograms  =  {"beta_ttbar"                        : ["#beta of ttbar system",                                                "Events", 50, [     0,     1]],
                    "cos_theta1k_antiLep"               : ["cos(#theta_{antilep}^{k})" ,                                           "Events", 24, [    -1,     1]],
                    "cos_theta1r_antiLep"               : ["cos(#theta_{antilep}^{r})" ,                                           "Events", 24, [    -1,     1]],
                    "cos_theta1n_antiLep"               : ["cos(#theta_{antilep}^{n})" ,                                           "Events", 24, [    -1,     1]],
                    "cos_theta1kStar_antiLep"           : ["cos(#theta_{antilep}^{k*})" ,                                          "Events", 24, [    -1,     1]],
                    "cos_theta1rStar_antiLep"           : ["cos(#theta_{antilep}^{r*})" ,                                          "Events", 24, [    -1,     1]],
                    "cos_theta2k_Lep"                   : ["cos(#theta_{lep}^{k})" ,                                               "Events", 24, [    -1,     1]],
                    "cos_theta2r_Lep"                   : ["cos(#theta_{lep}^{r})" ,                                               "Events", 24, [    -1,     1]],
                    "cos_theta2n_Lep"                   : ["cos(#theta_{lep}^{n})" ,                                               "Events", 24, [    -1,     1]],
                    "cos_theta2kStar_Lep"               : ["cos(#theta_{lep}^{k*})" ,                                              "Events", 24, [    -1,     1]],
                    "cos_theta2rStar_Lep"               : ["cos(#theta_{lep}^{r*})" ,                                              "Events", 24, [    -1,     1]],
                    "cos_theta1k"                       : ["cos(#theta_{1}^{k})" ,                                                 "Events", 24, [    -1,     1]],
                    "cos_theta1r"                       : ["cos(#theta_{1}^{r})" ,                                                 "Events", 24, [    -1,     1]],
                    "cos_theta1n"                       : ["cos(#theta_{1}^{n})" ,                                                 "Events", 24, [    -1,     1]],
                    "cos_theta1kStar"                   : ["cos(#theta_{1}^{k*})" ,                                                "Events", 24, [    -1,     1]],
                    "cos_theta1rStar"                   : ["cos(#theta_{1}^{r*})" ,                                                "Events", 24, [    -1,     1]],
                    "cos_theta2k"                       : ["cos(#theta_{2}^{k})" ,                                                 "Events", 24, [    -1,     1]],
                    "cos_theta2r"                       : ["cos(#theta_{2}^{r})" ,                                                 "Events", 24, [    -1,     1]],
                    "cos_theta2n"                       : ["cos(#theta_{2}^{n})" ,                                                 "Events", 24, [    -1,     1]],
                    "cos_theta2kStar"                   : ["cos(#theta_{2}^{k*})" ,                                                "Events", 24, [    -1,     1]],
                    "cos_theta2rStar"                   : ["cos(#theta_{2}^{r*})" ,                                                "Events", 24, [    -1,     1]],
                    "Cnn"                               : ["C_{nn}" ,                                                              "Events", 24, [    -1,     1]],
                    "Cnr"                               : ["C_{nr}" ,                                                              "Events", 24, [    -1,     1]],
                    "Cnk"                               : ["C_{nk}" ,                                                              "Events", 24, [    -1,     1]],
                    "Crn"                               : ["C_{rn}" ,                                                              "Events", 24, [    -1,     1]],
                    "Crr"                               : ["C_{rr}" ,                                                              "Events", 24, [    -1,     1]],
                    "Crk"                               : ["C_{rk}" ,                                                              "Events", 24, [    -1,     1]],
                    "Ckn"                               : ["C_{kn}" ,                                                              "Events", 24, [    -1,     1]],
                    "Ckr"                               : ["C_{kr}" ,                                                              "Events", 24, [    -1,     1]],
                    "Ckk"                               : ["C_{kk}" ,                                                              "Events", 24, [    -1,     1]],
                    "Crk_plus"                          : ["C_{rk}+C_{kr}" ,                                                       "Events", 24, [    -1,     1]],
                    "Crk_minus"                         : ["C_{rk}-C_{kr}" ,                                                       "Events", 24, [    -1,     1]],
                    "Cnr_plus"                          : ["C_{nr}+C_{rn}" ,                                                       "Events", 24, [    -1,     1]],
                    "Cnr_minus"                         : ["C_{nr}-C_{rn}" ,                                                       "Events", 24, [    -1,     1]],
                    "Cnk_plus"                          : ["C_{nk}+C_{kn}" ,                                                       "Events", 24, [    -1,     1]],
                    "Cnk_minus"                         : ["C_{nk}-C_{kn}" ,                                                       "Events", 24, [    -1,     1]],
                    "cHel"                              : ["cos(#phi_{lb})" ,                                                      "Events", 24, [    -1,     1]],
                    "cHel_Mtt300_400"                   : ["cos(#phi_{lb}) M_{tt} [300,400] GeV",                                  "Events", 24, [    -1,     1]],
                    "cHel_Mtt300_400_betaLT0p9"         : ["cos(#phi_{lb}) M_{tt} [300,400] GeV #beta_{tt}<0.9",                   "Events", 24, [    -1,     1]],
                    "cHel_P3n"                          : ["cos(#phi_{(P3n)lb})" ,                                                 "Events", 24, [    -1,     1]],
                    "cHel_P3n_Mtt800_Inf"               : ["cos(#phi_{(P3n)lb}) M_{tt} [800,Inf] GeV",                             "Events", 24, [    -1,     1]],
                    "cHel_P3n_Mtt800_Inf_cosThetaLT0p4" : ["cos(#phi_{(P3n)lb}) M_{tt} [800,Inf] GeV |cos(#theta_{tt}^{CM})|<0.4", "Events", 24, [    -1,     1]],
                    }
else:
	histograms  =  {#"beta_ttbar"                        : ["#beta of ttbar system",                                                "Events", 50, [     0,     1]],
                    "cos_theta1k_antiLep"               : ["cos(#theta_{antilep}^{k})" ,                                           "Events", 24, [    -1,     1]],
                    # "cos_theta1r_antiLep"               : ["cos(#theta_{antilep}^{r})" ,                                           "Events", 24, [    -1,     1]],
                    # "cos_theta1n_antiLep"               : ["cos(#theta_{antilep}^{n})" ,                                           "Events", 24, [    -1,     1]],
                    # "cos_theta1kStar_antiLep"           : ["cos(#theta_{antilep}^{k*})" ,                                          "Events", 24, [    -1,     1]],
                    # "cos_theta1rStar_antiLep"           : ["cos(#theta_{antilep}^{r*})" ,                                          "Events", 24, [    -1,     1]],
                    # "cos_theta2k_Lep"                   : ["cos(#theta_{lep}^{k})" ,                                               "Events", 24, [    -1,     1]],
                    # "cos_theta2r_Lep"                   : ["cos(#theta_{lep}^{r})" ,                                               "Events", 24, [    -1,     1]],
                    # "cos_theta2n_Lep"                   : ["cos(#theta_{lep}^{n})" ,                                               "Events", 24, [    -1,     1]],
                    # "cos_theta2kStar_Lep"               : ["cos(#theta_{lep}^{k*})" ,                                              "Events", 24, [    -1,     1]],
                    # "cos_theta2rStar_Lep"               : ["cos(#theta_{lep}^{r*})" ,                                              "Events", 24, [    -1,     1]],
                    # "cos_theta1k"                       : ["cos(#theta_{1}^{k})" ,                                                 "Events", 24, [    -1,     1]],
                    # "cos_theta1r"                       : ["cos(#theta_{1}^{r})" ,                                                 "Events", 24, [    -1,     1]],
                    # "cos_theta1n"                       : ["cos(#theta_{1}^{n})" ,                                                 "Events", 24, [    -1,     1]],
                    # "cos_theta1kStar"                   : ["cos(#theta_{1}^{k*})" ,                                                "Events", 24, [    -1,     1]],
                    # "cos_theta1rStar"                   : ["cos(#theta_{1}^{r*})" ,                                                "Events", 24, [    -1,     1]],
                    # "cos_theta2k"                       : ["cos(#theta_{2}^{k})" ,                                                 "Events", 24, [    -1,     1]],
                    # "cos_theta2r"                       : ["cos(#theta_{2}^{r})" ,                                                 "Events", 24, [    -1,     1]],
                    # "cos_theta2n"                       : ["cos(#theta_{2}^{n})" ,                                                 "Events", 24, [    -1,     1]],
                    # "cos_theta2kStar"                   : ["cos(#theta_{2}^{k*})" ,                                                "Events", 24, [    -1,     1]],
                    # "cos_theta2rStar"                   : ["cos(#theta_{2}^{r*})" ,                                                "Events", 24, [    -1,     1]],
                    # "Cnn"                               : ["C_{nn}" ,                                                              "Events", 24, [    -1,     1]],
                    # "Cnr"                               : ["C_{nr}" ,                                                              "Events", 24, [    -1,     1]],
                    # "Cnk"                               : ["C_{nk}" ,                                                              "Events", 24, [    -1,     1]],
                    # "Crn"                               : ["C_{rn}" ,                                                              "Events", 24, [    -1,     1]],
                    # "Crr"                               : ["C_{rr}" ,                                                              "Events", 24, [    -1,     1]],
                    # "Crk"                               : ["C_{rk}" ,                                                              "Events", 24, [    -1,     1]],
                    # "Ckn"                               : ["C_{kn}" ,                                                              "Events", 24, [    -1,     1]],
                    # "Ckr"                               : ["C_{kr}" ,                                                              "Events", 24, [    -1,     1]],
                    # "Ckk"                               : ["C_{kk}" ,                                                              "Events", 24, [    -1,     1]],
                    # "Crk_plus"                          : ["C_{rk}+C_{kr}" ,                                                       "Events", 24, [    -1,     1]],
                    # "Crk_minus"                         : ["C_{rk}-C_{kr}" ,                                                       "Events", 24, [    -1,     1]],
                    # "Cnr_plus"                          : ["C_{nr}+C_{rn}" ,                                                       "Events", 24, [    -1,     1]],
                    # "Cnr_minus"                         : ["C_{nr}-C_{rn}" ,                                                       "Events", 24, [    -1,     1]],
                    # "Cnk_plus"                          : ["C_{nk}+C_{kn}" ,                                                       "Events", 24, [    -1,     1]],
                    # "Cnk_minus"                         : ["C_{nk}-C_{kn}" ,                                                       "Events", 24, [    -1,     1]],
                    # "cHel"                              : ["cos(#phi_{lb})" ,                                                      "Events", 24, [    -1,     1]],
                    # "cHel_Mtt300_400"                   : ["cos(#phi_{lb}) M_{tt} [300,400] GeV",                                  "Events", 24, [    -1,     1]],
                    # "cHel_Mtt300_400_betaLT0p9"         : ["cos(#phi_{lb}) M_{tt} [300,400] GeV #beta_{tt}<0.9",                   "Events", 24, [    -1,     1]],
                    # "cHel_P3n"                          : ["cos(#phi_{(P3n)lb})" ,                                                 "Events", 24, [    -1,     1]],
                    # "cHel_P3n_Mtt800_Inf"               : ["cos(#phi_{(P3n)lb}) M_{tt} [800,Inf] GeV",                             "Events", 24, [    -1,     1]],
                    # "cHel_P3n_Mtt800_Inf_cosThetaLT0p4" : ["cos(#phi_{(P3n)lb}) M_{tt} [800,Inf] GeV |cos(#theta_{tt}^{CM})|<0.4", "Events", 24, [    -1,     1]],
                    }


### The sample_names array has filenames of samples that will contribute to each plot
sample_names = ["TTToSemiLeptonic", "TTTo2L2Nu", "TTToHadronic", "ST", "WJets", "QCD"]


### The stackList dictionary maps the list of samples to their color
# All three TTbar decay channels
stackList = {"TTToSemiLeptonic":[kRed], "TTTo2L2Nu":[kRed+2], "TTToHadronic":[kRed-7], "ST":[kBlue], "WJets":[kGreen],  "QCD":[kYellow]}

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
# legendHeightPer = 0.01
# legList = stackList.keys() 
# legend = TLegend(.59, .80, .89, .90)

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
pad1.cd()


### We loop through each histogram and fill it with the different samples
for histName in histograms:
    legend = TLegend(0.7,0.7,0.9,0.9)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetTextSize(0.03)
    print "--- Working on the", histName, "histogram ---"
    for sample in sample_names:
        #print "sample is", sample
        _file[sample] = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample),"read")
        # tree_MC[histName][sample]=_file[sample].Get("DNN_output0_General")
        # tree_MC[histName][sample].Draw("%s>>h_%s_%s(%i,%f,%f)"%(histName,histName,sample,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]))
        # hist[histName][sample] = tree_MC[histName][sample].GetHistogram()
        hist[histName][sample]=_file[sample].Get("DNN_output0_General/%s"%histName)
        hist[histName][sample].SetFillColor(stackList[sample][0])
        hist[histName][sample].SetLineColor(stackList[sample][0])
        hist[histName][sample].SetYTitle(histograms[histName][1])    
        print "Filling with", sample
        stack[histName].Add(hist[histName][sample])
        # rebin by a factor of 4
        hist[histName][sample].Rebin(4)
        # Normalize by bin width (0.08333333)
        binWidth = 0.3333333
        for b in range(1, hist[histName][sample].GetNbinsX()+1):
            content = hist[histName][sample].GetBinContent(b)
            error = hist[histName][sample].GetBinError(b)
            hist[histName][sample].SetBinContent(b, content/binWidth)
            hist[histName][sample].SetBinError(b, error/binWidth)
        legend.AddEntry(hist[histName][sample], sample, 'f')
    # canvas.cd()
    # canvas.ResetDrawn()
    # canvas.Draw()
    # canvas.cd()
    pad1.cd()
    stack[histName].Draw("HIST")
    legend.Draw()
    
    # Set vertical-axis limits
    maxVal = stack[histName].GetMaximum()
    # minVal = max(stack[histName].GetStack()[0].GetMinimum(), 1)
    minVal = 0.

    if Log:
        stack[histName].SetMaximum(10**(1.5*log10(maxVal) - 0.5*log10(minVal)))
        stack[histName].SetMinimum(minVal)
    else:
        stack[histName].SetMaximum(1.2*maxVal)
        stack[histName].SetMinimum(0.8*minVal)
    
    pad1.cd()
    pad1.SetLogy(Log)
    
    stack[histName].Draw("HIST")
    stack[histName].SetTitle('')
    stack[histName].GetXaxis().SetTitle('')
    stack[histName].GetXaxis().SetTitle(histograms[histName][0])
    stack[histName].GetXaxis().SetLabelSize(0.06)
    stack[histName].GetYaxis().SetTitle("Events / %.3f"%(binWidth))
    # stack[histName].GetYaxis().SetTitleSize(0.06)
    # stack[histName].GetYaxis().SetLabelSize(0.06)
    stack[histName].GetYaxis().SetTitleOffset(1.5)
    
    # CMS_lumi parameters: (hist_pad, year, iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right))
    CMS_lumi.CMS_lumi(pad1, 18, 11)

    pad1.Update()

    if Log:
        canvas.SaveAs("%s/%s_log.png"%(plotDirectory,histName))
    else:
        canvas.SaveAs("%s/%s.png"%(plotDirectory,histName))
    print "\n"

print "Congratulations, you've successfully generated some Plots."
print "You can find said plots at", plotDirectory