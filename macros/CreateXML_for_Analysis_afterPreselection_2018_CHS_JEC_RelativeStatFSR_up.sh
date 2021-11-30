#!/bin/bash
#where UHH2 code installed
pathGL_code=/nfs/dust/cms/user/ricardo/CMSSW_10_2_17/src/UHH2/
#where (NOT MERGED) trees after preselection stored
path_data=/nfs/dust/cms/group/zprime-uhh/Presel_2018_CHS_sys/workdir_Zprime_PreSelec_2018_RelativeStatFSR_up/ 

mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_2018_common_RelativeStatFSR_up
cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_2018_common_RelativeStatFSR_up


# #MC

for sample_name in TTToSemiLeptonic_2018 TTToHadronic_2018 TTTo2L2Nu_2018 #WW_2018 WZ_2018 ZZ_2018 ST_t-channel_antitop_2018 ST_t-channel_top_2018 ST_tW_antitop_5f_inclusiveDecays_2018 ST_tW_top_5f_inclusiveDecays_2018 ST_s-channel_lepton_2018 ST_s-channel_hadronic_2018  WJetsToLNu_HT-70To100_2018 WJetsToLNu_HT-100To200_2018 WJetsToLNu_HT-200To400_2018 WJetsToLNu_HT-400To600_2018 WJetsToLNu_HT-600To800_2018 WJetsToLNu_HT-800To1200_2018 WJetsToLNu_HT-1200To2500_2018 WJetsToLNu_HT-2500ToInf_2018 ZprimeToTT_M1000_W100_2018 QCD_HT1000to1500_2018 QCD_HT100to200_2018 QCD_HT1500to2000_2018 QCD_HT2000toInf_2018 QCD_HT200to300_2018 QCD_HT300to500_2018 QCD_HT500to700_2018 QCD_HT50to100_2018 QCD_HT700to1000_2018 DYJetsToLL_M-50_HT-100to200_2018 DYJetsToLL_M-50_HT-1200to2500_2018 DYJetsToLL_M-50_HT-200to400_2018 DYJetsToLL_M-50_HT-2500toInf_2018 DYJetsToLL_M-50_HT-400to600_v7_2018 DYJetsToLL_M-50_HT-600to800_2018 DYJetsToLL_M-50_HT-800to1200_2018


do
    echo $sample_name

    echo ${path_data}"MC."${sample_name}"*.root"

       $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"uhh2.AnalysisModuleRunner.MC."${sample_name}"*.root" MC_$sample_name.xml
       python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 MC_$sample_name.xml True
done

# # #DATA
#for sample_name in DATA_EGamma_Run2018A DATA_EGamma_Run2018B DATA_EGamma_Run2018C DATA_EGamma_Run2018D DATA_SingleMuon_Run2018A DATA_SingleMuon_Run2018B DATA_SingleMuon_Run2018C DATA_SingleMuon_Run2018D 

#do
#    echo $sample_name 
#    $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"uhh2.AnalysisModuleRunner.DATA."${sample_name}"*.root" DATA_$sample_name.xml
#    python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 DATA_$sample_name.xml True

#done
pwd
cd $pathGL_code/ZprimeSemiLeptonic/macros
