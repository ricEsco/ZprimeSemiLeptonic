#!/bin/bash


cd /data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/
kinit rescobar@FNAL.GOV

xrdcp spincorr_skim_TTToSemiLeptonic_UL18.root root://cmseos.fnal.gov//store/user/rescobar/Semilep_AnalysisOutput_UL18/Skim_TTToSemiLeptonicUL18_GenAndReco_SpinCor.root

xrdcp reweighting/closureTest/gendump_Ntuples.root root://cmseos.fnal.gov//store/user/rescobar/Semilep_AnalysisOutput_UL18/Inclusive_TTToSemiLeptonicUL18_Gen_SpinCorr.root

xrdcp reweighting/strengthmaps_uhh2/NoSC_strengthmap_lb.root root://cmseos.fnal.gov//store/user/rescobar/Semilep_AnalysisOutput_UL18/NoSC_strengthmap_lb.root