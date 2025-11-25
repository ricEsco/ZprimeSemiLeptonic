from ROOT import TFile, TLegend, TCanvas, TPad, THStack, TF1, TPaveText, TGaxis, SetOwnership, TObject, gStyle,TH1F,TH2F
from ROOT import *
import os
import math
import sys
from optparse import OptionParser
from sys import argv
from numpy import log10
from array import array
from collections import OrderedDict

parser = OptionParser()
parser.add_option("-c", "--channel", dest="channel",default="muon",type='str',help="what channel?")
parser.add_option("-y", "--year", dest="year", default="2018",type='str',help="don't print status messages to stdout")
# parser.add_option("-p", "--plots", dest="plots", action='store_false',default='False',type=bool,help="eft var plots?")
parser.add_option("-m", "--mass_range", dest="mass_range", help="Specify the mass range (0_500, 500_750, 750-1000, 1000-1500, 1500Inf)", type='str')
parser.add_option("-r", "--region", dest="region", help="Specify the region (SR, CR1, CR2)", type='str')
parser.add_option("-s", "--region_score", dest="region_score", help="Specify the region score (0, 1, 2)", type='str')


(options, args) = parser.parse_args()

year = options.year if options.year else "UL18"
mass_range = options.mass_range if options.mass_range else "0_500"
region = options.region if options.region else "CR2"

region_score = options.region_score if options.region_score else "0"

(options, args) = parser.parse_args()
stackList_orig={"TTbar":[kRed],"WJets":[kGreen], "DYJets":[kGreen], "Diboson":[kGreen], "QCD":[kGreen],"ST":[kBlue],"ST_tW":[kBlue],"ST_s":[kBlue],"ST_t":[kBlue]}

channel = options.channel
year=options.year
# plots =option.plots


gROOT.SetBatch(True)
path="%s_%s_EFT"%(channel,year)
print(path)
if os.path.exists(path):
	print("output directory exists")
else:
	os.mkdir(path)
YesLog = True
NoLog=False
padRatio = 0.25
padOverlap = 0.15

padGap = 0.01
if (channel=="electron" or channel=="muon"):
	#fileDir ="/data/dust/user/titasroy/Ac_UL_ntuples/%s/%s/workdir_AnalysisDNN_%s_%s_dY/NOMINAL/"%(year,channel,year,channel)
    fileDir="/data/dust/user/titasroy/Ac_UL_ntuples/%s/%s/workdir_AnalysisDNN_%s_%s//NOMINAL/"%(year,channel,year, channel)
else:
	fileDir ="/data/dust/user/titasroy/Ac_UL_ntuples/%s/lepton/"%(year)

import CMS_lumi

from Style import *

thestyle = Style()

HasCMSStyle = False
style = None
if os.path.isfile('tdrstyle.C'):
    ROOT.gROOT.ProcessLine('.L tdrstyle.C')
    ROOT.setTDRStyle()
    #print ("Found tdrstyle.C file, using this style.")
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

if year=="2018":
	period=6
elif year=="2017":
	period=5
elif "2016" in year:
	period=4
elif "all" in year:
	period=8
#print stackList
if channel=="muon":
	_channelText = "#mu+jets"
elif(channel=="electron"):
	_channelText="e+jets"
elif(channel=="lepton"):
	_channelText="l+jets"
H = 600;
W = 800;


# references for T, B, L, R                                                                                                             
T = 0.08*H
B = 0.12*H
L = 0.12*W
R = 0.1*W
legendHeightPer = 0.04
legList = stackList_orig.keys() 
#legList.reverse()

legendStart = 0.69
legendEnd = 0.97-(R/W)
canvasRatio = TCanvas('c1Ratio','c1Ratio',W,H)
canvasRatio.SetFillColor(0)
canvasRatio.SetBorderMode(0)
canvasRatio.SetFrameFillStyle(0)
canvasRatio.SetFrameBorderMode(0)
canvasRatio.SetLeftMargin( L/W )
canvasRatio.SetRightMargin( R/W )
canvasRatio.SetTopMargin( T/H )
canvasRatio.SetBottomMargin( B/H )
canvasRatio.SetTickx(1)
canvasRatio.SetTicky(1)
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
pad1.SetTickx(1)
pad1.SetTicky(1)

pad2.SetFillColor(0)
pad2.SetFillStyle(4000)
pad2.SetBorderMode(0)
pad2.SetFrameFillStyle(0)
pad2.SetFrameBorderMode(0)
pad2.SetTickx(1)
pad2.SetTicky(1)


canvasRatio.cd()
pad1.Draw()
pad2.Draw()



histograms={"DeltaY_reco_d2":["\DeltaY (t_{h}(p_T)<150 GeV,\Delta \phi<0 )", "Events", 2, [-2.5,2.5]],
		"DeltaY_reco_d1":["\DeltaY (t_{h}(p_T)<150 GeV,\Delta \phi>0 )", "Events", 2, [-2.5,2.5]],
		"Sigma_phi_1":["\Sigma \phi (top_{h}(p_{T})>150 GeV,\DeltaY > 0)", "Events", 16, [-3.2,3.2]],
		"Sigma_phi_2":["\Sigma \phi (top_{h}(p_{T})>150 GeV,\DeltaY < 0)", "Events", 16, [-3.2,3.2]],
		}






stackList = { "TTbar":[kRed],"WJets":[kGreen], "ST":[kBlue]}

hist_combined={}
hist_sys_comb={}
categories=["%s_SR"%(mass_range),"%s_CR1"%(mass_range),"%s_CR2"%(mass_range)]



systematic_name_mapping = {
    # "mu_reco": "muonReco",
    # "pu": "pu",
    # "prefiring": "prefiringWeight",
    "mu_id_stat": "muonID_stat",
    "mu_id_syst": "muonID_syst",
    "mu_iso_stat": "muonIso_stat",
    "mu_iso_syst": "muonIso_syst",
    "mu_trigger_stat": "muonTrigger_stat",
    "mu_trigger_syst": "muonTrigger_syst",
    "ele_id" : "eleID",
    "ele_trigger": "eleTrigger",
    "ele_reco": "eleReco",
    "isr": "isr",
    "fsr": "fsr",
    "btag_cferr1": "btagCferr1",
    "btag_cferr2": "btagCferr2",
    "btag_hf": "btagHf",
    "btag_hfstats1": "btagHfstats1",
    "btag_hfstats2": "btagHfstats2",
    "btag_lf": "btagLf",
    "btag_lfstats1": "btagLfstats1",
    "btag_lfstats2": "btagLfstats2",
    "ttag_corr": "ttagCorr",
    "ttag_uncorr": "ttagUncorr",
    # "murmuf_upup": "murmuf_upup",
    # "murmuf_upnone":"murmuf_upnone",
    # "murmuf_noneup":"murmuf_noneup",
    # "murmuf_nonedown":"murmuf_nonedown",
    # "murmuf_downnone": "murmuf_downnone",
    # "murmuf_downdown":"murmuf_downdown",
    # "tmistag": "tmistag"
}

# for i in range(1,101):
#     print("updating: ",i)
#     systematic_name_mapping.update({"PDF_%i"%(i):"PDF_%i"%(i),})




histo_={}

samples=["DATA","DYJets", "Diboson","QCD", "WJets", "ST","TTbar","MC_EFT"]
for cat in categories:
    hist_combined[cat]={}
    hist_sys_comb[cat]={}
    for sample in samples:
        hist_combined[cat][sample]=TH1F("%s_%s"%(cat,sample),"%s_%s"%(cat,sample),36,1.,37.)
        hist_sys_comb[cat][sample]={}
        for sys in systematic_name_mapping:
            hist_sys_comb[cat][sample][sys]={}
            if sample!="DATA":
                if sample!="MC_EFT":
                # print(cat,sample,sys)
                    hist_sys_comb[cat][sample][sys]["Up"]=TH1F("%s_%s_%s_Up"%(cat,sample,sys),"%s_%s_%s_Up"%(cat,sample,sys),36,1.,37.)
                    hist_sys_comb[cat][sample][sys]["Down"]=TH1F("%s_%s_%s_Down"%(cat,sample,sys),"%s_%s_%s_Down"%(cat,sample,sys),36,1.,37.)



histo_sys={}
for sample in samples:
    histo_sys[sample]={}
    for sys in systematic_name_mapping:
        # print(sample,sys)
        histo_sys[sample][sys]={}
file={}
# plots=True


# print("check if TH1F exists: ", hist_sys_comb['0_500_SR']['DYJets']['prefiring']["Up"])
# if plots:
for cat in categories:
    # print("is systematic histo defined here at the start of loop? : ", hist_sys_comb['0_500_SR']['DYJets']['prefiring']["Up"])
    stack = THStack("hs","stack")
    legendR = TLegend(2*legendStart - legendEnd , 0.99 - (T/H)/(1.-padRatio+padOverlap) - legendHeightPer/(1.-padRatio+padOverlap)*round((len(legList)+1)/2.), legendEnd, 0.99-(T/H)/(1.-padRatio+padOverlap))
    legendR.SetNColumns(2)
    legendR.SetBorderSize(0)
    legendR.SetFillColor(0)
    samples=["DATA","DYJets", "Diboson","QCD", "WJets", "ST","TTbar"]
    for sample in ["DATA","DYJets", "Diboson","QCD", "WJets", "ST","TTbar","MC_EFT"]:
        print sample
        for histo in histograms:
        # print("at histo loop: ",hist_sys_comb['0_500_SR']['DYJets']['prefiring']["Up"])
        # for sample in ["DATA","DYJets", "Diboson","QCD", "WJets", "ST","TTbar"]:
            if sample=="DATA":
                file[sample]=TFile("%s/uhh2.AnalysisModuleRunner.DATA.DATA.root"%(fileDir),"read")
                temp_hist="DeltaY_reco_%s_General/%s"%(cat,histo)
            elif sample=="MC_EFT":
                file[sample]=TFile("%s/uhh2.AnalysisModuleRunner.MC.MC_EFT.root"%(fileDir),"read")
                temp_hist="DeltaY_reco_PDFVariations_%s/%s_PDF_1"%(cat,histo)
            else:
                file[sample] = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(fileDir,sample),"read")
                temp_hist="DeltaY_reco_%s_General/%s"%(cat,histo)
            
            # temp_hist="DeltaY_reco_%s_General/%s"%(cat,histo)
            histo_[sample]=file[sample].Get(temp_hist)
            if (sample!="DATA"):
                if sample!="MC_EFT":
                    print(sample)
                    for sys in systematic_name_mapping:
                        # print sys, sample, cat
                        if "PDF" in sys :
                            sys_hist_up="DeltaY_reco_PDFVariations_%s/%s_%s"%(cat,histo,sys)
                            sys_hist_down="DeltaY_reco_PDFVariations_%s/%s_%s"%(cat,histo,sys)
                            
                        if "murmuf" in sys:
                            sys_hist_up="DeltaY_reco_SystVariations_%s/%s_%s"%(cat,histo,sys)
                            sys_hist_down="DeltaY_reco_SystVariations_%s/%s_%s"%(cat,histo,sys)
                        
                        else:
                            sys_hist_up="DeltaY_reco_SystVariations_%s/%s_%s_up"%(cat,histo,sys)
                            sys_hist_down="DeltaY_reco_SystVariations_%s/%s_%s_down"%(cat,histo,sys)
                        
                        histo_sys[sample][sys]["up"]=file[sample].Get(sys_hist_up)
                        histo_sys[sample][sys]["down"]=file[sample].Get(sys_hist_down)
                        if histo=="DeltaY_reco_d1":
                            for i in range(histo_sys[sample][sys]["up"].GetNbinsX()):
                                print(sample, cat, sys)
                                hist_sys_comb[cat][sample][sys]["Up"].SetBinContent(1+i,histo_sys[sample][sys]["up"].GetBinContent(i+1))
                                hist_sys_comb[cat][sample][sys]["Up"].SetBinError(1+i,histo_sys[sample][sys]["up"].GetBinError(i+1))

                                hist_sys_comb[cat][sample][sys]["Down"].SetBinContent(1+i,histo_sys[sample][sys]["down"].GetBinContent(i+1))
                                hist_sys_comb[cat][sample][sys]["Down"].SetBinError(1+i,histo_sys[sample][sys]["down"].GetBinError(i+1))
                        if histo=="DeltaY_reco_d2":
                            for i in range(histo_sys[sample][sys]["up"].GetNbinsX()):
                                hist_sys_comb[cat][sample][sys]["Up"].SetBinContent(3+i,histo_sys[sample][sys]["up"].GetBinContent(i+1))
                                hist_sys_comb[cat][sample][sys]["Up"].SetBinError(3+i,histo_sys[sample][sys]["up"].GetBinError(i+1))

                                hist_sys_comb[cat][sample][sys]["Down"].SetBinContent(3+i,histo_sys[sample][sys]["down"].GetBinContent(i+1))
                                hist_sys_comb[cat][sample][sys]["Down"].SetBinError(3+i,histo_sys[sample][sys]["down"].GetBinError(i+1))
                        if histo=="Sigma_phi_1":
                            for i in range(histo_sys[sample][sys]["up"].GetNbinsX()):
                                hist_sys_comb[cat][sample][sys]["Up"].SetBinContent(5+i,histo_sys[sample][sys]["up"].GetBinContent(i+1))
                                hist_sys_comb[cat][sample][sys]["Up"].SetBinError(5+i,histo_sys[sample][sys]["up"].GetBinError(i+1))

                                hist_sys_comb[cat][sample][sys]["Down"].SetBinContent(5+i,histo_sys[sample][sys]["down"].GetBinContent(i+1))
                                hist_sys_comb[cat][sample][sys]["Down"].SetBinError(5+i,histo_sys[sample][sys]["down"].GetBinError(i+1))
                        if histo=="Sigma_phi_2":
                            for i in range(histo_sys[sample][sys]["up"].GetNbinsX()):
                                hist_sys_comb[cat][sample][sys]["Up"].SetBinContent(21+i,histo_sys[sample][sys]["up"].GetBinContent(i+1))
                                hist_sys_comb[cat][sample][sys]["Up"].SetBinError(21+i,histo_sys[sample][sys]["up"].GetBinError(i+1))

                                hist_sys_comb[cat][sample][sys]["Down"].SetBinContent(21+i,histo_sys[sample][sys]["down"].GetBinContent(i+1))
                                hist_sys_comb[cat][sample][sys]["Down"].SetBinError(21+i,histo_sys[sample][sys]["down"].GetBinError(i+1))
            if histo=="DeltaY_reco_d1":
                for i in range(histo_[sample].GetNbinsX()):
                    # print(histo_[sample].GetBinContent(i+1))
                    hist_combined[cat][sample].SetBinContent(1+i,histo_[sample].GetBinContent(i+1))
                    hist_combined[cat][sample].SetBinError(1+i,histo_[sample].GetBinError(i+1))
            elif histo=="DeltaY_reco_d2":
                for i in range(histo_[sample].GetNbinsX()):
                    hist_combined[cat][sample].SetBinContent(3+i,histo_[sample].GetBinContent(i+1))
                    hist_combined[cat][sample].SetBinError(3+i,histo_[sample].GetBinError(i+1))
            elif histo=="Sigma_phi_1":
                for i in range(histo_[sample].GetNbinsX()):
                    hist_combined[cat][sample].SetBinContent(5+i,histo_[sample].GetBinContent(i+1))	
                    hist_combined[cat][sample].SetBinError(5+i,histo_[sample].GetBinError(i+1))	
            elif histo=="Sigma_phi_2":
                for i in range(histo_[sample].GetNbinsX()):
                    hist_combined[cat][sample].SetBinContent(21+i,histo_[sample].GetBinContent(i+1))	
                    hist_combined[cat][sample].SetBinError(21+i,histo_[sample].GetBinError(i+1))
            # print("is systematic histo right before sys loop? : ", hist_sys_comb['0_500_SR']['DYJets']['prefiring']["Up"])
        if sample!="DATA":
            if sample!="MC_EFT":
                hist_combined[cat][sample].SetFillColor(stackList_orig[sample][0])
                hist_combined[cat][sample].SetLineColor(stackList_orig[sample][0])
        if sample=="TTbar":
            print("writing ttbar")
            legendR.AddEntry(hist_combined[cat][sample],"t#bar{t}",'f')
            stack.Add(hist_combined[cat][sample])
        if sample=="ST":
            legendR.AddEntry(hist_combined[cat][sample],"ST",'f')
            stack.Add(hist_combined[cat][sample])
        if sample=="WJets":
            print("writing wjets")
            legendR.AddEntry(hist_combined[cat][sample],"W+jets",'f')
        total_mc=hist_combined[cat]["TTbar"].GetBinContent(1)+hist_combined[cat]["ST"].GetBinContent(1)
        print("total MC: ",total_mc)
    
    hist_combined[cat]["WJets"].Add(hist_combined[cat]["Diboson"])
    hist_combined[cat]["WJets"].Add(hist_combined[cat]["DYJets"])
    stack.Add(hist_combined[cat]["WJets"])
    hist_combined[cat]["DATA"].SetMarkerStyle(20)
    hist_combined[cat]["DATA"].SetMarkerColor(kBlack)
    hist_combined[cat]["DATA"].SetMarkerSize(0.8)
    hist_combined[cat]["DATA"].SetLineColor(kBlack)
    hist_combined[cat]["DATA"].SetYTitle("Events")     
    hist_combined[cat]["DATA"].Draw("pe,x0")

    hist_combined[cat]["MC_EFT"].SetMarkerStyle(42)
    hist_combined[cat]["MC_EFT"].SetMarkerColor(7)
    hist_combined[cat]["MC_EFT"].SetMarkerSize(0.8)
    hist_combined[cat]["MC_EFT"].SetLineColor(7)
    hist_combined[cat]["MC_EFT"].SetYTitle("Events")     
    hist_combined[cat]["MC_EFT"].Draw("pe,x0")
    data_total=hist_combined[cat]["DATA"].GetBinContent(1)
    print("DATA: ",hist_combined[cat]["DATA"].GetBinContent(1) )
    print(cat)
    print("DATA/MC ratio: ",total_mc/data_total)
    errorban=stack.GetStack().Last().Clone("errorban")
    errorban.Sumw2()
    errorban.SetLineColor(kGray+2)
    errorban.SetFillColor(kGray+2)
    errorban.SetFillStyle(3245)
    errorban.SetMarkerSize(0)
    errorban.Draw("E2,SAME")
    legendR.AddEntry(errorban,"MC tot. unc.",'f')
    oneLine = TF1("oneline","1",-9e9,9e9)
    oneLine.SetLineColor(kBlack)
    oneLine.SetLineWidth(1)
    oneLine.SetLineStyle(2)	
    minVal = max(stack.GetStack()[0].GetMinimum(),1)
    log=1
    maxVal =stack.GetMaximum()
    minVal = max(stack.GetStack()[0].GetMinimum(),1)
    if log:

        stack.SetMaximum(10**(1.5*log10(maxVal) - 0.5*log10(minVal)))
    else:
        
        stack.SetMaximum(1.5*maxVal)
    stack.Draw("HIST,SAME")
    stack.GetXaxis().SetLabelSize(0)
    stack.GetYaxis().SetLabelSize(gStyle.GetLabelSize()/(1.-padRatio+padOverlap))
    stack.GetYaxis().SetTitleSize(gStyle.GetTitleSize()/(1.-padRatio+padOverlap))
    stack.GetYaxis().SetTitleOffset(gStyle.GetTitleYOffset()*(1.-padRatio+padOverlap))
    stack.GetYaxis().SetMaxDigits(5)
    stack.GetYaxis().SetTitle("Events")
    stack.SetMinimum(minVal)
    errorband=stack.GetStack().Last().Clone("error")
    errorband.Sumw2()
    errorband.SetLineColor(kBlack)
    errorband.SetFillColor(kBlack)
    errorband.SetFillStyle(3245)
    errorband.SetMarkerSize(0)
    errorband.SetLineColor(0)
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
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad2.SetFillColor(0)
    pad2.SetFillStyle(4000)
    pad2.SetBorderMode(0)
    pad2.SetFrameFillStyle(0)
    pad2.SetFrameBorderMode(0)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    canvasRatio.cd()
    canvasRatio.ResetDrawn()
    canvasRatio.Draw()
    canvasRatio.cd()
    pad1.Draw()
    pad2.Draw()
    pad1.cd()
    pad1.SetLogy(log)
    stack.Draw("HIST")
    hist_combined[cat]["DATA"].Draw("E,X0,SAME")
    hist_combined[cat]["MC_EFT"].Draw("E,X0,SAME")
    errorban.Draw("E2,SAME")
    legendR.AddEntry(hist_combined[cat]["DATA"], "Data", 'pe')
    legendR.AddEntry(hist_combined[cat]["MC_EFT"], "EFT ref pt", 'pe')
    ratio = hist_combined[cat]["DATA"].Clone("ratio")
    temp = stack.GetStack().Last().Clone("temp")
    for i_bin in range(1,temp.GetNbinsX()+1):
            temp.SetBinError(i_bin,0.)
    # ratio = hist_combined[cat]["DATA"].Clone("temp")
    ratio.Divide(temp)
    print("DATA/MC ratio: ",total_mc/data_total)
    print("ratio: ", ratio.GetBinContent(1))
    ratio.GetXaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
    ratio.GetYaxis().SetLabelSize(gStyle.GetLabelSize()/(padRatio+padOverlap))
    ratio.GetXaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
    ratio.GetYaxis().SetTitleSize(gStyle.GetTitleSize()/(padRatio+padOverlap))
    ratio.GetYaxis().SetTitleOffset(gStyle.GetTitleYOffset()*(padRatio+padOverlap-padGap))
    ratio.GetYaxis().SetRangeUser(0.2,1.8)
    ratio.GetYaxis().SetNdivisions(504)
    ratio.GetXaxis().SetTitle("EFT bins")
    ratio.GetYaxis().SetTitle("Data/MC")
    CMS_lumi.CMS_lumi(pad1, period, 11)
    legendR.SetTextSize(0.05)
    legendR.Draw("SAME")
    pad2.cd()
    ratio.SetMarkerStyle(hist_combined[cat]["DATA"].GetMarkerStyle())
    ratio.SetMarkerSize(hist_combined[cat]["DATA"].GetMarkerSize())
    ratio.SetLineColor(kBlack)
    ratio.SetMarkerColor(kBlack)
    ratio.SetLineWidth(hist_combined[cat]["DATA"].GetLineWidth())
    ratio.Draw('e,x0')
    errorband.Divide(temp)
    for i in range(1, errorband.GetNbinsX() + 1):
        if(errorban.GetBinContent(i) == 0):
            errorband.SetBinError(i,0)
        else:
            errorband.SetBinError(i, errorban.GetBinError(i)/errorban.GetBinContent(i))
    errorband.Draw('e2,same')
    oneLine.Draw("same")
    canvasRatio.Update()
    canvasRatio.RedrawAxis()
    if log:
        canvasRatio.SaveAs("%s/EFT_%s_log.pdf"%(path,cat))
    else:
        canvasRatio.SaveAs("%s/EFT_%s_linear.pdf"%(path,cat))
    canvasRatio.Clear()
# sys.exit()
hist_scale_Up={}
hist_scale_Down={}
v_samples_ttbar = ["TTbar"]
v_variations = ["upup", "upnone", "noneup", "nonedown", "downnone", "downdown"]
combine_file_name = 'EFT_{}_{}_{}.root'.format(year, channel,mass_range)
combine_file = TFile(combine_file_name, 'RECREATE')

norm={}
# for cat in categories:
#     for sys in ["murmuf_upup", "murmuf_upnone", "murmuf_noneup", "murmuf_nonedown", "murmuf_downnone", "murmuf_downdown"]:
#         norm[sys]= hist_sys_comb[cat]["TTbar"][sys]["Up"].GetBinContent(1)/hist_combined[cat]["TTbar"].GetBinContent(1)
#         print("normalizing by: ", norm[sys])
#         hist_sys_comb[cat]["TTbar"][sys]["Up"].Scale(1/norm[sys])
#     hist_scale_Up[cat]=hist_sys_comb[cat]["TTbar"]["murmuf_upup"]["Up"].Clone("mcscale_Up")
#     hist_scale_Down[cat]=hist_sys_comb[cat]["TTbar"]["murmuf_upup"]["Up"].Clone("mcscale_Down")
#     for j in range(hist_scale_Up[cat].GetNbinsX()):
#         max=0.
#         min=10000000000.
#         for sys in ["murmuf_upup", "murmuf_upnone", "murmuf_noneup", "murmuf_nonedown", "murmuf_downnone", "murmuf_downdown"]:
#             print(j,sys)
#             bin_content = hist_sys_comb[cat]["TTbar"][sys]["Up"].GetBinContent(j)
#             if(bin_content > max) :
#                 max = bin_content
#             if(bin_content < min) :
#                 min = bin_content

#         hist_scale_Up[cat].SetBinContent(j, max)
#         hist_scale_Down[cat].SetBinContent(j, min)
# #PDF rms calculation
# norm_pdf=array('d')
# for cat in categories:
#     for sample in sample_name:
#         norm_pdf=array('d')
#         for i in range(1,101):
#             norm_pdf.append(hist_sys_comb[cat][sample]["PDF_%i"%(i)]["Down"].GetBinContent(1)/hist_combined[cat][sample].GetBinContent(1))
            




# creat root file for combine inp
for cat in categories:
    for sample in ["DYJets","Diboson","QCD"]:
        # hist_combined[cat]["WJets"].Add(hist_combined[cat][sample])
        for sys in systematic_name_mapping:
            hist_sys_comb[cat]["WJets"][sys]["Up"].Add(hist_sys_comb[cat][sample][sys]["Up"])
            hist_sys_comb[cat]["WJets"][sys]["Down"].Add(hist_sys_comb[cat][sample][sys]["Down"])

cat_region={"%s_SR"%(mass_range):["SR"],"%s_CR1"%(mass_range):["CR1"],"%s_CR2"%(mass_range):["CR2"]}


sample_name={"TTbar":["TTbar"],"ST":["ST"],"WJets":["WJets"],"DATA":["data_obs"],"MC_EFT":["EFT_refpt"]}

systematic_name_mapping = {
    # "mu_reco": "muonReco",
    # "pu": "pu",
    # "prefiring": "prefiringWeight",
    "mu_id_stat": "muonID_stat",
    "mu_id_syst": "muonID_syst",
    "mu_iso_stat": "muonIso_stat",
    "mu_iso_syst": "muonIso_syst",
    "mu_trigger_stat": "muonTrigger_stat",
    "mu_trigger_syst": "muonTrigger_syst",
    "ele_id" : "eleID",
    "ele_trigger": "eleTrigger",
    "ele_reco": "eleReco",
    "isr": "isr",
    "fsr": "fsr",
    "btag_cferr1": "btagCferr1",
    "btag_cferr2": "btagCferr2",
    "btag_hf": "btagHf",
    "btag_hfstats1": "btagHfstats1",
    "btag_hfstats2": "btagHfstats2",
    "btag_lf": "btagLf",
    "btag_lfstats1": "btagLfstats1",
    "btag_lfstats2": "btagLfstats2",
    "ttag_corr": "ttagCorr",
    "ttag_uncorr": "ttagUncorr",
    # "murmuf": "murmuf",
    # "tmistag": "tmistag"
}

for cat in categories:
    combine_file.mkdir(cat_region[cat][0])
    # print("making directory: ", cat_region[cat][0])
    combine_file.cd(cat_region[cat][0])
    for sample in ["TTbar","ST","WJets","DATA","MC_EFT"]:
        hist_combined[cat][sample].Write(sample_name[sample][0])
        if sample!="DATA" :
            if sample!="MC_EFT":
                for sys in systematic_name_mapping:
                    # print(cat,sample,sys)
                    if "murmuf" in sys and sample!="TTbar":
                        continue
                    if "murmuf" in sys and sample=="TTbar":
                        hist_scale_Up[cat].Write("%s_%sUp"%(sample,sys))
                        hist_scale_Down[cat].Write("%s_%sDown"%(sample,sys))
                    else:
                        hist_sys_comb[cat][sample][sys]["Up"].Write("%s_%sUp"%(sample,sys))
                        hist_sys_comb[cat][sample][sys]["Down"].Write("%s_%sDown"%(sample,sys))
                




        # for sys in systematic_name_mapping.items():
        #     if (debug): print sys

        #     for variation in ["up", "down"]:
        #         sys_hist_name = "DeltaY_{}_{}".format(sys, variation)
        #         sys_hist = inFile.Get("DeltaY_reco_SystVariations_{}_{}/".format(mass_range, region) + sys_hist_name)
        #         if sys_hist:
        #             output_hist_name = "{}_{}{}".format(sample, new_sys_name, variation.capitalize())
        #             sys_hist.Clone(output_hist_name).Write()

combine_file.Close()









    


            
			