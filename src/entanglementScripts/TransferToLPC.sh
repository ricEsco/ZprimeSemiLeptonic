#!/bin/bash


cd /data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/
kinit rescobar@FNAL.GOV

xrdcp reweighting/TTToSemiLeptonicUL18_SCskim.root root://cmseos.fnal.gov//store/user/rescobar/Semilep_AnalysisOutput_UL18/TTToSemiLeptonicUL18_GenAndReco_SpinCorrSkim.root

# xrdcp reweighting/strengthmaps_uhh2/NoSC_alphaStrengthMap_lb.root root://cmseos.fnal.gov//store/user/rescobar/Semilep_AnalysisOutput_UL18/NoSCreweighting_alphaStrengthMap_lb.root