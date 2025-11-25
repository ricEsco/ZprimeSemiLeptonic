#!/bin/bash

#where UHH2 code installed
pathGL_code=/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/
#where (NOT MERGED) trees after preselection stored
path_data=/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_JEC/UL17/Preselection/workdir_Preselection_UL17_JEC_down_TimePtEta/uhh2.AnalysisModuleRunner.

mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL17_forAnalysis_JEC_down_TimePtEta
cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL17_forAnalysis_JEC_down_TimePtEta


# MC
for sample_name in TTToSemiLeptonic_UL17 TTToHadronic_UL17 TTTo2L2Nu_UL17
do
    echo $sample_name
    $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"MC."${sample_name}"*.root" MC_$sample_name.xml
    python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 MC_$sample_name.xml True
done

pwd
cd $pathGL_code/ZprimeSemiLeptonic/macros