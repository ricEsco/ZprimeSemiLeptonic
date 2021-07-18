#!/bin/bash

StringVal='
DYJetsToLL_M-50_HT-600to800_2018_1.xml
DYJetsToLL_M-50_HT-800to1200_2018_1.xml
QCD_HT1000to1500_2018_1.xml
QCD_HT100to200_2018_1.xml
QCD_HT1500to2000_2018_1.xml
QCD_HT2000toInf_2018_1.xml
QCD_HT200to300_2018_1.xml
QCD_HT300to500_2018_1.xml
QCD_HT500to700_2018_1.xml
QCD_HT50to100_2018_1.xml
QCD_HT700to1000_2018_1.xml
ST_s-channel_hadronic_2018_1.xml
ST_s-channel_lepton_2018_1.xml
ST_t-channel_antitop_2018_1.xml
ST_t-channel_top_2018_1.xml
ST_tW_antitop_5f_inclusiveDecays_2018_1.xml
ST_tW_top_5f_inclusiveDecays_2018_1.xml
TTTo2L2Nu_2018_1.xml
TTTo2L2Nu_2018_2.xml
TTToHadronic_2018_1.xml
TTToHadronic_2018_2.xml
TTToHadronic_2018_3.xml
TTToSemiLeptonic_2018_1.xml
TTToSemiLeptonic_2018_2.xml
TTToSemiLeptonic_2018_3.xml
WJetsToLNu_HT-100To200_2018_1.xml
WJetsToLNu_HT-1200To2500_2018_1.xml
WJetsToLNu_HT-200To400_2018_1.xml
WJetsToLNu_HT-2500ToInf_2018_1.xml
WJetsToLNu_HT-400To600_2018_1.xml
WJetsToLNu_HT-600To800_2018_1.xml
WJetsToLNu_HT-70To100_2018_1.xml
WJetsToLNu_HT-800To1200_2018_1.xml
WW_2018_1.xml
WZ_2018_1.xml
ZZ_2018_1.xml
'
for val in $StringVal;



do
    sframe_main ${val}
done


