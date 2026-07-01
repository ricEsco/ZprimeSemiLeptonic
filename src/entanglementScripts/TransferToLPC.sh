#!/bin/bash


cd /data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/
kinit rescobar@FNAL.GOV
xrdcp spincorr_skim_TTToSemiLeptonic_UL18.root root://cmseos.fnal.gov//store/user/rescobar/Semilep_AnalysisOutput_UL18

