#!/bin/bash

echo "Make sure to run this script from the directory where the files being added live \n" 


echo "Adding all DATA root files first"
hadd uhh2.AnalysisModuleRunner.DATA.DATA.root uhh2.AnalysisModuleRunner.DATA.DATA*


echo "Adding all Diboson root files"
hadd uhh2.AnalysisModuleRunner.MC.Diboson.root uhh2.AnalysisModuleRunner.MC.WW_UL18_* uhh2.AnalysisModuleRunner.MC.WZ_UL18_* uhh2.AnalysisModuleRunner.MC.ZZ_UL18_*


echo "Adding Drell-Yan files together"
hadd uhh2.AnalysisModuleRunner.MC.DY.root uhh2.AnalysisModuleRunner.MC.DYJetsToLL_M-50_HT-*


echo "Adding Single Top files together"
hadd uhh2.AnalysisModuleRunner.MC.ST.root uhh2.AnalysisModuleRunner.MC.ST_*


echo "Adding QCD Jets files together"
hadd uhh2.AnalysisModuleRunner.MC.QCD.root uhh2.AnalysisModuleRunner.MC.QCD_HT*


echo "Adding WJets files together"
hadd uhh2.AnalysisModuleRunner.MC.WJets.root uhh2.AnalysisModuleRunner.MC.WJetsToLNu_HT-*


echo "Adding TTToHadronic files together"
hadd uhh2.AnalysisModuleRunner.MC.TTToHadronic.root uhh2.AnalysisModuleRunner.MC.TTToHadronic_UL18_*


echo "Adding TTTo2L2Nu files together"
hadd uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_0.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_1.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_2.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_3.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_4.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_5.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_6.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_7.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_8.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_9.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_10.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_11.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_12.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_13.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_14.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_15.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_16.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_17.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_18.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_19.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_20.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_21.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_22.root
hadd uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_2.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_23.root uhh2.AnalysisModuleRunner.MC.TTTo2L2Nu_UL18_24.root


echo "Adding TTToSemiLeptonic files together"
hadd uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_0.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_1.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_2.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_3.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_4.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_5.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_6.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_7.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_8.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_9.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_10.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_11.root
hadd uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_2.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_12.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_13.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_14.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_15.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_16.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_17.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_18.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_19.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_20.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_21.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_22.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_23.root
hadd uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_3.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_24.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_25.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_26.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_27.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_28.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_29.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_30.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_31.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_32.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_33.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_34.root uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18_35.root

