#!/usr/bin/env python
"""
script to combine EFT sample and other MC samples for a single combine output root file

"""


import ROOT
import numpy as np
import sys
import argparse
import json
import os
import math
from optparse import OptionParser
from sys import argv
from numpy import log10
from array import array
from ROOT import TH1F


def get_hist(file,hist_name):
    file_=ROOT.TFile(file)  
    hist = file_.Get(hist_name)
    hist.SetDirectory(0)
    return hist

def rebin_hist(hist):
    n_bins=16
    xbins=[-3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2]
    hist=hist.Rebin(n_bins, hist.GetName(), array('d',xbins))
    return hist

def write_hist_to(root_file, hist, eft,region):
    file_=ROOT.TFile(root_file,"UPDATE")
    file_.cd(region)
    hist.Write("EFT_"+eft)
    return file_


categoris=["SR","CR1","CR2"]
hist_combined={}
for cat in categoris:
    hist_combined[cat]=TH1F("%s"%(cat),"%s"%(cat),36,1.,37.)

def make_one_histogram(hist1,hist2,hist3,hist4):
    if hist1=="dyreco_1":
        for i in range(histo.GetNbinsX()):
            print(histo,histo.GetNbinsX())
            hist_combined[cat].SetBinContent(1+i,histo.GetBinContent(i+1))
            hist_combined[cat].SetBinError(1+i,histo.GetBinError(i+1))
    elif hist2=="dyreco_2":
        for i in range(histo.GetNbinsX()):
            hist_combined[cat].SetBinContent(3+i,histo.GetBinContent(i+1))
            hist_combined[cat].SetBinError(3+i,histo.GetBinError(i+1))
    elif hist3=="Sigma_phi_1":
        for i in range(histo.GetNbinsX()):
            hist_combined[cat].SetBinContent(5+i,histo.GetBinContent(i+1))	
            hist_combined[cat].SetBinError(5+i,histo.GetBinError(i+1))	
    elif hist4=="Sigma_phi_2":
        for i in range(histo.GetNbinsX()):
            hist_combined[cat].SetBinContent(21+i,histo.GetBinContent(i+1))	
            hist_combined[cat].SetBinError(21+i,histo.GetBinError(i+1))
    return hist_combined[cat]


categoris=["SR","CR1","CR2"]

variables=["Sigma_phi_1","Sigma_phi_2","dyreco_1","dyreco_2"]

SM_rebin_hist={}
ctGRe_10_rebin_hist={}
for cat in categoris:
    SM_rebin_hist[cat]={}
    ctGRe_10_rebin_hist[cat]={}
    for var in variables:
        if cat=="SR":
            hist_name=var
        else:
            hist_name=var+"_"+cat
        print hist_name
        SM_hist_= get_hist("reweighted_histograms_0_500.root",hist_name+"_SM")
        SM_rebin_hist[cat][hist_name]=rebin_hist(SM_hist_)
        ctGRe_10_hist_= get_hist("reweighted_histograms_0_500.root",hist_name+"_ctGRe_10")
        ctGRe_10_rebin_hist[cat][hist_name]=rebin_hist(ctGRe_10_hist_)


for cat in categoris:
    if cat=="SR":
        hist_combined[cat]=make_one_histogram(SM_rebin_hist[cat]["dyreco_1"],SM_rebin_hist[cat]["dyreco_2"],SM_rebin_hist[cat]["Sigma_phi_1"],SM_rebin_hist[cat]["Sigma_phi_2"])
        write_hist_to("EFT_2017_muon_0_500.root", hist_combined[cat],"SM", cat)
        print("writing SR SM to file")
        hist_combined[cat]=make_one_histogram(ctGRe_10_rebin_hist[cat]["dyreco_1"],ctGRe_10_rebin_hist[cat]["dyreco_2"],ctGRe_10_rebin_hist[cat]["Sigma_phi_1"],ctGRe_10_rebin_hist[cat]["Sigma_phi_2"])
        final_file=write_hist_to("EFT_2017_muon_0_500.root", hist_combined[cat],"ctGRe_10",cat)
    else:
        hist_combined[cat]=make_one_histogram(SM_rebin_hist[cat]["dyreco_1_"+cat],SM_rebin_hist[cat]["dyreco_2_"+cat],SM_rebin_hist[cat]["Sigma_phi_1_"+cat],SM_rebin_hist[cat]["Sigma_phi_2_"+cat])
        write_hist_to("EFT_2017_muon_0_500.root", hist_combined[cat],"SM", cat)
        hist_combined[cat]=make_one_histogram(ctGRe_10_rebin_hist[cat]["dyreco_1_"+cat],ctGRe_10_rebin_hist[cat]["dyreco_2_"+cat],ctGRe_10_rebin_hist[cat]["Sigma_phi_1_"+cat],ctGRe_10_rebin_hist[cat]["Sigma_phi_2_"+cat])
        final_file=write_hist_to("EFT_2017_muon_0_500.root", hist_combined[cat],"ctGRe_10",cat)


final_file.Close()
