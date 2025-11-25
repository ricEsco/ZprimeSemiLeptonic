#!/bin/bash

text2workspace.py datacard_muon_UL17_0_500_SR.txt -o Ac_UL17_0_500.root  -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125 --PO map='.*/TTbar_1:r_neg[1,0,20]' --PO map='.*/TTbar_2:r_pos=expr;;r_pos("30636270000/30836339000*@0*(100+@1)/(100-@1)",r_neg,Ac[-2,-5,0])' --PO verbose 
