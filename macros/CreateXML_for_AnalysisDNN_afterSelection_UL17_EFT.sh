#!/bin/bash

#where UHH2 code installed
pathGL_code=/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2
#where (NOT MERGED) trees after selection stored
path_data=/data/dust/group/cms/zprime-uhh/Analysis_EFT_UL17/muon/workdir_Analysis_EFT_UL17_muon_semilepton_deltayreco_deltaPhi_sigmaPhi_ttree_new/uhh2.AnalysisModuleRunner.

mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL17_EFT_muon_analysisDNN_semilepton
cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL17_EFT_muon_analysisDNN_semilepton

# MC
for sample_name in MC_EFT_Mttbar_0-700_UL17 MC_EFT_Mttbar_700-900_UL17 MC_EFT_Mttbar_900-Inf_UL17
do
    echo $sample_name
    $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"MC."${sample_name}"*.root" MC_$sample_name.xml
    python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 MC_$sample_name.xml True
done

pwd
cd $pathGL_code/ZprimeSemiLeptonic/macros
