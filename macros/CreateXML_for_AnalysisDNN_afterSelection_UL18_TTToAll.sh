#!/bin/bash

# where UHH2 code installed
pathGL_code=/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/
# path to UNMERGED Analysis selection (before DNN) output files
path_data=/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNwTtag/muon/workdir_Analysis_UL18_muon/uhh2.AnalysisModuleRunner.

mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_forAnalysisDNN_muon_UL18_TTToAll
cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_forAnalysisDNN_muon_UL18_TTToAll

# MC
for sample_name in TTToSemiLeptonic_UL18 TTTo2L2Nu_UL18 #TTToHadronic_UL18
do
    echo $sample_name
    $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"MC."${sample_name}"*.root" MC_$sample_name.xml
    python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 MC_$sample_name.xml True
done

pwd
cd $pathGL_code/ZprimeSemiLeptonic/macros
