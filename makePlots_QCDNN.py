#from ROOT import TFile, TLegend, TCanvas, TPad, THStack, TF1, TPaveText, TGaxis, SetOwnership, TObject, gStyle,TH1F
from ROOT import *
import os

import sys
from optparse import OptionParser
from numpy import log10
from array import array


parser = OptionParser()
parser.add_option("-c", "--channel", dest="channel", default="mu",type='str',
                     help="Specify which channel mu or ele? default is mu" )
parser.add_option("--Log","--isLog", dest="isLog", default=False,action="store_true",
                     help="Plot the plots in log ?" )
parser.add_option("--QCD","--isQCD", dest="isQCD", default=False,action="store_true",
                     help="Plot as QCD and others" )
parser.add_option("--DNN","--isDNN", dest="isDNN", default=False,action="store_true",
                     help="Plot to check DNN" )

(options, args) = parser.parse_args()


padRatio = 0.25
padOverlap = 0.15
padGap = 0.01


DNN = options.isDNN
channel = options.channel
Log=options.isLog
QCD = options.isQCD

if channel=="ele":
	_channelText = "e+jets"
	plotDirectory = "data_pre_plots_ele"
#	_fileDir = "/nfs/dust/cms/user/titasroy/Run2_analysis_DNN/electron/2018/workdir_eleDec/NOMINAL"
	_fileDir = "/nfs/dust/cms/group/zprime-uhh/AnalysisDNN_UL18_MET20/electron/workdir_Analysis_UL18_electron/NOMINAL/"

else:
	_channelText = "#mu+jets"
	#plotDirectory = "data_pre_plots_mu"
	plotDirectory = "analysisplots_mu_preDNN"
        #_fileDir = "/nfs/dust/cms/user/titasroy/Run2_analysis_DNN/muon/2018/workdir_Zprime_Analysis_muon_QCDNN_May19/NOMINAL_1"
        #_fileDir = "/nfs/dust/cms/user/titasroy/Run2_analysis_DNN/muon/2018/workdir_muon_DecNNimpl/NOMINAL"
        #_fileDir = "/nfs/dust/cms/group/zprime-uhh/Run2_analysis/muon/2018/workdir_June//NOMINAL"
	_fileDir = "/nfs/dust/cms/group/zprime-uhh/AnalysisDNN_UL18_MET20/muon/workdir_Zprime_AnalysisDNN_UL18_muon/NOMINAL"

gROOT.SetBatch(True)


# Histogram Information:
# [X-axis title, 
#  Y-axis title,
#  Rebinning factor,
#  [x-min,x-max], -1 means keep as is
#  Extra text about region
#  log plot]
regionText ="loose selection"

import CMS_lumi

from Style import *
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

#stackList = { "TTbar_Semi_1":[kRed], "DY":[kGreen], "QCD":[kYellow],"WJets":[kBlue], "ST":[kOrange], "Diboson":[kTeal]}
stackList = { "TTbarTo2L_1":[kRed], "TTbarTo2L_2":[kRed], "TTbarToHad":[kRed], "TTbarToSemi_1":[kRed], "TTbarToSemi_2":[kRed],"TTbarToSemi_3":[kRed],"TTbarToSemi_4":[kRed],"TTbarToSemi_5":[kRed], "TTbarToSemi_6":[kRed], "DY":[kGreen], "QCD":[kYellow],"WJets":[kBlue], "ST":[kOrange], "Diboson":[kTeal]}

if QCD:
	stackList = { "TTbarTo2L_1":[kRed], "TTbarTo2L_2":[kRed], "TTbarToHad":[kRed], "TTbarToSemi_1":[kRed], "TTbarToSemi_2":[kRed],"TTbarToSemi_3":[kRed],"TTbarToSemi_4":[kRed],"TTbarToSemi_5":[kRed], "TTbarToSemi_6":[kRed], "DY":[kRed], "QCD":[kYellow],"WJets":[kRed], "ST":[kRed], "Diboson":[kRed]}

#print "stackList is ", stackList

CMS_lumi.channelText = _channelText
CMS_lumi.writeChannelText = True
CMS_lumi.writeExtraText = True

H = 600;
W = 800;

# references for T, B, L, R                                                                                                             
T = 0.08*H
B = 0.12*H
L = 0.12*W
R = 0.1*W

# SetOwnership(canvas, False)
# SetOwnership(canvasRatio, False)
# SetOwnership(pad1, False)
# SetOwnership(pad2, False)

#legendHeightPer = 0.04
legendHeightPer = 0.01
legList = stackList.keys() 
#legList.reverse()

legendStart = 0.69
legendEnd = 0.97-(R/W)

#legend = TLegend(2*legendStart - legendEnd, 1-T/H-0.01 - legendHeightPer*(len(legList)+1), legendEnd, 0.99-(T/H)-0.01)
legend = TLegend(2*legendStart - legendEnd , 0.99 - (T/H)/(1.-padRatio+padOverlap) - legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.), legendEnd, 0.99-(T/H)/(1.-padRatio+padOverlap))
legend.SetNColumns(2)

_file={}

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

sum_=0
tree_MC={}
hist={}
if channel=="mu":
	#histograms = {"phi_jet": ["#phi^{AK4 jets}", "Events", 35, [-3.5, 3.5]]} #debugging
	histograms = {"ZprimeCandidateBestChi2":["#Chi^{2} Zprime Best Candidate","Events", 50, [0, 10]]}
        '''
        histograms = {
              "met_pt"   : ["missing E_{T} [GeV]", "Events", 15, [0,1500]],
              
              "st"    : ["S_{T} [GeV]",  "Events", 50, [0,5000]],
              "st_lep"    : ["S_{T}^{leps} [GeV]",  "Events", 50, [0,5000]],
              "st_jets"    : ["S_{T}^{jets} [GeV]",  "Events", 50, [0,5000]],
              
              "ht" : ["H_{T} [GeV]",  "Events", 500, [0,2000]],
              "ht_lep" : ["H_{T}^{leps} [GeV]",  "Events", 500, [0,2000]],

              "pt_jet" : ["p_{T}^{jets} [GeV]" , "Events", 50, [100, 900]],
              "pt_jet1": ["p_{T}^{jet 1} [GeV]", "Events", 50, [100, 900]],
              "pt_jet2": ["p_{T}^{jet 2} [GeV]", "Events", 50, [100, 900]],
              "pt_jet3": ["p_{T}^{jet 3} [GeV]", "Events", 50, [100, 900]],
              "eta_jet": ["#eta^{AK4jets}",  "Events", 25, [-2.5, 2.5]],
              "eta_jet1":["#eta^{AK4jet 1}", "Events", 25, [-2.5, 2.5]],
              "eta_jet2":["#eta^{AK4jet 2}", "Events", 25, [-2.5, 2.5]],
              "eta_jet3":["#eta^{AK4jet 3}", "Events", 25, [-2.5, 2.5]],
              "phi_jet":  ["#phi^{AK4 jets}",  "Events", 35, [-3.5, 3.5]],
              "phi_jet1": ["#phi^{AK4 jet 1}", "Events", 35, [-3.5, 3.5]],
              "phi_jet2": ["#phi^{AK4 jet 2}", "Events", 35, [-3.5, 3.5]],
              "phi_jet3": ["#phi^{AK4 jet 2}", "Events", 35, [-3.5, 3.5]],
              "mass_jet":  ["mass of all AK4 jets","Events",40,[0,400]],
              "mass_jet1": ["mass of AK4 jet 1",   "Events",40,[0,400]],
              "mass_jet2": ["mass of AK4 jet 2",   "Events",40,[0,400]],
              "mass_jet3": ["mass of AK4 jet 3",   "Events",40,[0,400]],
              
              "pt_mu":[ "p_{T}^{#mu} [GeV]", "Events",50, [50, 500]],
              "eta_mu":["#eta^{#mu}", "Events", 25, [-2.5, 2.5]],
              "phi_mu": ["#phi^{#mu}", "Events", 35, [-3.5, 3.5]],
              "reliso_mu":[ "#mu rel. Iso", "Events", 20, [0, 0.5]],
              
              "dR_mu_jet":[ "#DeltaR_{min}(#mu, jet)","Events", 60,[ 0, 6]],
              "ptrel_mu_jet":["p_{T}^{rel}(#mu1, jet)", "Events",50, [0, 500]],
              "dphi_mu_jet1":["#Delta#phi(#mu, jet1)","Events", 40, [-3.5, 3.5]],
              "dphi_mu_MET":  ["#Delta#phi(#mu, MET)","Events", 40, [-0.5, 3.5]],
              
              "deepjetbscore_jet":  ["DeepJet b-tag score all AK4 jets", "Events",20, [0, 1]],
              "deepjetbscore_jet1": ["DeepJet b-tag score of AK4 jet 1", "Events",20, [0, 1]],
              "deepjetbscore_jet2": ["DeepJet b-tag score of AK4 jet 2", "Events",20, [0, 1]],
              "deepjetbscore_jet3": ["DeepJet b-tag score of AK4 jet 3", "Events",20, [0, 1]],
              
              "Invpt_mu":["p_{T}^{#mu}/H_{T}", "Events",10, [0, 1]],   
              "Invpt_jet1": ["p_{T}^{jet 1}/H_{T}", "Events", 100, [0, 10]],  
              "InvS11": ["S11/H_{T}", "Events", 70, [-0.001,0.008]],
              "InvS12": ["S12/H_{T}", "Events", 80, [-0.004,0.004]],
              "InvS13": ["S13/H_{T}", "Events", 80, [-0.004,0.004]],
              "InvS22": ["S22/H_{T}", "Events", 70, [-0.001,0.008]],
              "InvS23": ["S23/H_{T}", "Events", 100, [-0.005,0.005]],
              "InvS33": ["S33/H_{T}", "Events", 70, [-0.001,0.008]],
              "InvMET" : ["missing E_{T}/H_{T}", "Events", 100, [0,1]],
              "InvST":["S_{T}/H_{T}",  "Events", 50, [0,10]]              
	}
        '''
else:
	histograms = {"met_pt"   : ["missing E_{T} [GeV]", "Events", 15, [0,1500]],
              "InvMET" : ["missing E_{T}/H_{T}", "Events", 100, [0,1]],
              "InvS11": ["S11/H_{T}", "Events", 70, [-0.001,0.008]],
              "InvS12": ["S12/H_{T}", "Events", 80, [-0.004,0.004]],
              "InvS13": ["S13/H_{T}", "Events", 80, [-0.004,0.004]],
              "InvS22": ["S22/H_{T}", "Events", 70, [-0.001,0.008]], 
              "InvS23": ["S23/H_{T}", "Events", 100, [-0.005,0.005]],
              "InvS33": ["S33/H_{T}", "Events", 70, [-0.001,0.008]],
              "InvST":["S_{T}/H_{T}",  "Events", 50, [0,10]],
         #     "InvSTlep":["S_{T}^{ele}/H_{T}",  "Events", 50, [0,10]],
        #      "InvSTjets":["S_{T}^{jets}/H_{T}",  "Events", 50, [0,10]],
              "st"    : ["S_{T} [GeV]",  "Events", 50, [0,5000]],
              "st_lep"    : ["S_{T}^{ele} [GeV]",  "Events", 50, [0,5000]],
              "st_jets"    : ["S_{T}^{jets} [GeV]",  "Events", 50, [0,5000]],
              "ht" : ["H_{T} [GeV]",  "Events", 500, [0,2000]],
              "mass_jet": ["mass of all AK4 jets","Events",40,[0,400]],
           #   "phi_Ak8Puppijets": ["#phi^{AK8Puppi jets}", "Events", 35, [-3.5, 3.5]], 
              "phi_jet": ["#phi^{AK4 jets}", "Events", 35, [-3.5, 3.5]],
              "phi_ele": ["#phi^{ele}", "Events", 35, [-3.5, 3.5]],
              "dphi_jet1_MET":["#Delta#phi(jet1, MET)","Events", 40, [-0.5, 3.5]],
              "dphi_ele_MET":  ["#Delta#phi(ele, MET)","Events", 40, [-0.5, 3.5]],   
              "deepjetbscore_jet": ["DeepJet b-tag score all AK4 jets", "Events",20, [0, 1]],
              "deepjetbscore_jet1": ["DeppJet b-tag score AK4 jet 1}", "Events",20, [0, 1]],
              "deepjetbscore_jet2": ["DeppJet b-tag score AK4 jet 2}", "Events",20, [0, 1]],
              "deepjetbscore_jet3": ["DeppJet b-tag score AK4 jet 3}", "Events",20, [0, 1]],
              "pt_jet1": ["p_{T}^{jet 1} [GeV]", "Events", 50, [100, 900]],
              "pt_jet": ["p_{T}^{jet} [GeV]","Events", 50, [100, 900]],
              "pt_ele":[ "p_{T}^{ele} [GeV]", "Events",50, [50, 500]],
       #       "pt_Ak8Puppijet1":[ "p_{T}^{AK8Puppi jets} [GeV]", "Events",40, [400, 4400]],
              "reliso_ele":[ "ele rel. Iso", "Events", 20, [0, 0.5]],
              "dR_ele_jet":[ "#DeltaR_{min}(ele, jet)","Events", 60,[ 0, 6]],
        #      "dR_ele_Ak8Puppijet":[ "#DeltaR_{min}(ele, Ak8jet)","Events", 60,[ 0, 6]],
         #     "dR_jet_Ak8Puppijet": [ "#DeltaR_{min}(Ak4jet, Ak8jet)","Events", 60,[ 0, 6]],
          #    "InvdR_ele_Ak8Puppijet":[ "#DeltaR_{min}(ele, Ak8jet)/H_{T}","Events", 20,[ 0, 0.02]],
#              "InvdR_ele_jet":[ "#DeltaR_{min}(ele, jet)/H_{T}","Events", 40,[ 0, 0.04]],
           #   "InvdR_jet_Ak8Puppijet": [ "#DeltaR_{min}(Ak4jet, Ak8jet)/H_{T}","Events", 40,[ 0, 0.04]],
 #             "Invdeepjetbscore_jet": ["DeepJet b-tag score/H_{T} all AK4 jets", "Events",200, [0, 0.01]],
  #            "Invdeepjetbscore_jet1": ["DeppJet b-tag score/H_{T} AK4 jet 1}", "Events",200, [0, 0.01]],
   #           "Invdeepjetbscore_jet2": ["DeppJet b-tag score/H_{T} AK4 jet 2}", "Events",200, [0, 0.01]],
    #          "Invdeepjetbscore_jet3": ["DeppJet b-tag score/H_{T} AK4 jet 3}", "Events",200, [0, 0.01]],
     #         "Invdphi_jet1_MET":["#Delta#phi(jet1, MET)/H_{T}","Events", 60, [0, 0.06]],
      #        "Invdphi_ele_MET":  ["#Delta#phi(ele, MET)/H_{T}","Events", 60, [0, 0.06]],
       #       "Invdphi_ele_jet1":  ["#Delta#phi(ele, jet1)/H_{T}","Events", 60, [0, 0.06]], 
             # "Invdphi_ele_Ak8Puppijet1":  ["#Delta#phi(ele, Ak8Puppijet1)/H_{T}","Events", 60, [0, 0.06]],
            #  "InvmSD_Ak8Puppijet1":["mSD/H_{T} Ak8Puppijet1","Events", 20,[0,2.]],
#              "Invmass_jet": ["mass_{jet}/H_{T}", "Events", 50, [0,0.5 ]],
              #"Invpt_Ak8Puppijet1":["p_{T}^{Ak8Puppi jet1}/H_{T}","Events",100,[0,10]],
              "Invpt_ele":["p_{T}^{ele}/H_{T}", "Events",10, [0, 1]],
              "Invpt_jet1": ["p_{T}^{jet 1}/H_{T}", "Events", 100, [0, 10]],
        #      "Invptrel_ele_jet":["p_{T}^{rel}(ele1, jet)/H_{T}", "Events",100, [0, 10]],
         #     "Invreliso_ele":[ "ele rel. Iso/H_{T}", "Events", 20, [0, 0.1]], 
              "ptrel_ele_jet":["p_{T}^{rel}(ele1, jet)", "Events",50, [0, 500]],
              #"eta_Ak8Puppijets":["#eta^{AK8Puppi jets}", "Events", 25, [-2.5, 2.5]],
              "eta_jet":["#eta^{AK4jets}", "Events", 25, [-2.5, 2.5]],
              "eta_ele":["#eta^{ele}", "Events", 25, [-2.5, 2.5]],
       
         } 	




if DNN:
	 print  "I am here"
	 histograms = {"NNoutput0": ["QCD node", "Norm events/bin", 5, [0,1]],
			"NNoutput1": ["Others node", "Norm events/bin", 5, [0,1]],
			}

sample_names = ["QCD", "Diboson","ST", "DY", "WJets","TTbarTo2L_1", "TTbarTo2L_2", "TTbarToHad", "TTbarToSemi_1", "TTbarToSemi_2", "TTbarToSemi_3", "TTbarToSemi_4", "TTbarToSemi_5", "TTbarToSemi_6"]
stack={}
legendR={}
for histName in histograms:
	tree_MC[histName]={}
	hist[histName]={}
	stack[histName] = THStack("hs","stack")
        legendR[histName] = TLegend(2*legendStart - legendEnd , 0.99 - (T/H)/(1.-padRatio+padOverlap) - legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.)-0.1, legendEnd, 0.99-(T/H)/(1.-padRatio+padOverlap))
	legendR[histName].SetNColumns(2)
	legendR[histName].SetBorderSize(0)
	legendR[histName].SetFillColor(0)
canvas.cd()


if DNN:
	for histName in histograms:
		print histName	
		for sample in sample_names[:1]:
			print sample, sample_names[:1]
			_file[sample] = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample),"read")
			tree_MC[histName][sample]=_file[sample].Get("AnalysisTree")
                	tree_MC[histName][sample].Draw("%s>>h_%s_%s(%i,%f,%f)"%(histName,histName,sample,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]),"","goff")
                	hist[histName][sample] = tree_MC[histName][sample].GetHistogram()
			hist[histName][sample].SetLineColor(stackList[sample][0])
		 	print histName, sample
			print hist[histName][sample].GetEntries()
			hist[histName][sample].Draw("SAME")
		for sample in sample_names[1:]:
			print sample, sample_names[1:]
			_file[sample] = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample),"read")
                        tree_MC[histName][sample]=_file[sample].Get("AnalysisTree")
                        tree_MC[histName][sample].Draw("%s>>h_%s_%s(%i,%f,%f)"%(histName,histName,sample,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]),"","goff")
                        hist[histName][sample] = tree_MC[histName][sample].GetHistogram()
			print hist[histName][sample].GetEntries()
                        hist[histName][sample].SetLineColor(stackList[sample][0])
		for sample in ['Diboson', 'ST', 'DYJets', 'WJets','TTbar_Semi_2','TTbar_Semi_3','TTbar_Semi_4']:
			hist[histName]["TTbar_Semi_1"].Add(hist[histName][sample])
		print hist[histName]["TTbar_Semi_1"].GetEntries()
		hist[histName]["TTbar_Semi_1"].Draw("SAME")
		canvas.SetLogy(0)
		canvas.SaveAs("%s_noweight_nolog.pdf"%(histName))


for histName in histograms:
	print "--- Working on the", histName, "histogram ---"
	for sample in sample_names:
        	#print  "sample =", sample
		_file[sample] = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample),"read")
        	#print "%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(_fileDir,sample), "is the file being used"
		tree_MC[histName][sample]=_file[sample].Get("AnalysisTree")
                #print histName, sample, histograms[histName][2], histograms[histName][3][0], histograms[histName][3][1]
                tree_MC[histName][sample].Draw("%s>>h_%s_%s(%i,%f,%f)"%(histName,histName,sample,histograms[histName][2],histograms[histName][3][0],histograms[histName][3][1]),"weight*weight_pu")
                hist[histName][sample] = tree_MC[histName][sample].GetHistogram()
        	hist[histName][sample].SetFillColor(stackList[sample][0])
        	hist[histName][sample].SetLineColor(stackList[sample][0])
		if QCD:
			if sample=="QCD":
				legendR[histName].AddEntry(hist[histName][sample],sample,'f')
			elif sample=="TTbar":
				legendR[histName].AddEntry(hist[histName][sample],'Others','f')
		#trying to avoid TTbarToSemi* from appearing in legend multiple times
                elif sample=="TTbarTo2L_2":
                                hist[histName][sample].SetYTitle(histograms[histName][1])
                                print "Filling with", sample
                                stack[histName].Add(hist[histName][sample])
                                continue
                elif sample=="TTbarToSemi_2":
                                hist[histName][sample].SetYTitle(histograms[histName][1])
                                print "Filling with", sample
                                stack[histName].Add(hist[histName][sample])
                                continue
                elif sample=="TTbarToSemi_3":
                                hist[histName][sample].SetYTitle(histograms[histName][1])
                                print "Filling with", sample
                                stack[histName].Add(hist[histName][sample])
                                continue
                elif sample=="TTbarToSemi_4":
                                hist[histName][sample].SetYTitle(histograms[histName][1])
                                print "Filling with", sample
                                stack[histName].Add(hist[histName][sample])
                                continue
                elif sample=="TTbarToSemi_5":
                                hist[histName][sample].SetYTitle(histograms[histName][1])
                                print "Filling with", sample
                                stack[histName].Add(hist[histName][sample])
                                continue
                elif sample=="TTbarToSemi_6":
                                hist[histName][sample].SetYTitle(histograms[histName][1])
                                print "Filling with", sample
                                stack[histName].Add(hist[histName][sample])
                                continue

                else:
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
	stack[histName].Draw("HIST")
	dataHist.Draw("pe,x0,same")    

	oneLine = TF1("oneline","1",-9e9,9e9)
	oneLine.SetLineColor(kBlack)
	oneLine.SetLineWidth(1)
	oneLine.SetLineStyle(2)	

	maxVal = stack[histName].GetMaximum()
	minVal = 1
	minVal = max(stack[histName].GetStack()[0].GetMinimum(),1)
        if Log:
		#stack[histName].SetMaximum(10**(1.5*log10(maxVal) - 0.5*log10(minVal)))
                stack[histName].SetMaximum(15**(1.5*log10(maxVal) - 0.5*log10(minVal)))

	else:
        	stack[histName].SetMaximum(1.7*maxVal)
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

	dataHist.Draw("E,X0,SAME")

	legendR[histName].AddEntry(dataHist, "Data", 'pe')

	ratio = dataHist.Clone("temp")
	temp = stack[histName].GetStack().Last().Clone("temp")

	for i_bin in range(1,temp.GetNbinsX()+1):
       		temp.SetBinError(i_bin,0.)
	ratio.Divide(temp)

	#print ratio.GetMaximum(), ratio.GetMinimum()
	#max_=ratio.GetMaximum()*0.1+ratio.GetMaximum()
	#min_=ratio.GetMinimum()-0.1*ratio.GetMinimum()

	ratio.SetTitle('')

	ratio.GetXaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
	ratio.GetYaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
	ratio.GetXaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
	ratio.GetYaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
	ratio.GetYaxis().SetTitleOffset(gStyle.GetTitleYOffset()*(padRatio+padOverlap-padGap))
	ratio.GetYaxis().SetRangeUser(0.35,1.65)
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
	#pad2.Update()
	canvasRatio.Update()
	canvasRatio.RedrawAxis()
	if Log:
		#canvasRatio.SaveAs("%s/%s_log.pdf"%(plotDirectory,histName))
		canvasRatio.SaveAs("%s/%s_log.png"%(plotDirectory,histName))
	else:
		#canvasRatio.SaveAs("%s/%s.pdf"%(plotDirectory,histName))
		canvasRatio.SaveAs("%s/%s.png"%(plotDirectory,histName))
        print "\n"
       
