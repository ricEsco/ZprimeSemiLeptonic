#!/bin/bash
#where UHH2 code installed
pathGL_code=/nfs/dust/cms/user/deleokse/RunII_106_v2/CMSSW_10_6_28/src/UHH2/
#where (NOT MERGED) trees after preselection stored
path_data=/nfs/dust/cms/group/zprime-uhh/Presel_UL16preVFP/workdir_Preselection_UL16preVFP/uhh2.AnalysisModuleRunner.

mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL16preVFP
cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_UL16preVFP


# #MC

for sample_name in TTToSemiLeptonic_UL16preVFP TTToHadronic_UL16preVFP TTTo2L2Nu_UL16preVFP WJetsToLNu_HT-70To100_UL16preVFP WJetsToLNu_HT-100To200_UL16preVFP WJetsToLNu_HT-200To400_UL16preVFP WJetsToLNu_HT-400To600_UL16preVFP WJetsToLNu_HT-600To800_UL16preVFP WJetsToLNu_HT-800To1200_UL16preVFP WJetsToLNu_HT-1200To2500_UL16preVFP WJetsToLNu_HT-2500ToInf_UL16preVFP DYJetsToLL_M-50_HT-70to100_UL16preVFP DYJetsToLL_M-50_HT-100to200_UL16preVFP DYJetsToLL_M-50_HT-200to400_UL16preVFP DYJetsToLL_M-50_HT-400to600_UL16preVFP DYJetsToLL_M-50_HT-600to800_UL16preVFP DYJetsToLL_M-50_HT-800to1200_UL16preVFP DYJetsToLL_M-50_HT-1200to2500_UL16preVFP DYJetsToLL_M-50_HT-2500toInf_UL16preVFP WW_UL16preVFP WZ_UL16preVFP ZZ_UL16preVFP ST_tW_antitop_5f_inclusiveDecays_UL16preVFP ST_tW_top_5f_inclusiveDecays_UL16preVFP ST_tW_antitop_5f_NoFullyHadronicDecays_UL16preVFP ST_tW_top_5f_NoFullyHadronicDecays_UL16preVFP ST_t-channel_antitop_4f_InclusiveDecays_UL16preVFP ST_t-channel_top_4f_InclusiveDecays_UL16preVFP ST_s-channel_4f_leptonDecays_UL16preVFP QCD_HT50to100_UL16preVFP QCD_HT100to200_UL16preVFP QCD_HT200to300_UL16preVFP QCD_HT300to500_UL16preVFP QCD_HT500to700_UL16preVFP QCD_HT700to1000_UL16preVFP QCD_HT1000to1500_UL16preVFP QCD_HT1500to2000_UL16preVFP QCD_HT2000toInf_UL16preVFP ALP_ttbar_signal_UL16preVFP ALP_ttbar_interference_UL16preVFP HscalarToTTTo1L1Nu2J_m365_w36p5_res_UL16preVFP HscalarToTTTo1L1Nu2J_m400_w40p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m500_w50p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m600_w60p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m800_w80p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m1000_w100p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m365_w36p5_int_UL16preVFP HscalarToTTTo1L1Nu2J_m400_w40p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m500_w50p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m600_w60p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m800_w80p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m1000_w100p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m365_w36p5_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m400_w40p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m500_w50p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m600_w60p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m800_w80p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m1000_w100p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m365_w36p5_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m400_w40p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m500_w50p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m600_w60p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m800_w80p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m1000_w100p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m365_w91p25_res_UL16preVFP HscalarToTTTo1L1Nu2J_m400_w100p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m500_w125p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m600_w150p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m800_w200p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m1000_w250p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m365_w91p25_int_UL16preVFP HscalarToTTTo1L1Nu2J_m400_w100p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m500_w125p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m600_w150p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m800_w200p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m1000_w250p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m365_w91p25_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m400_w100p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m500_w125p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m600_w150p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m800_w200p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m1000_w250p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m365_w91p25_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m400_w100p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m500_w125p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m600_w150p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m800_w200p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m1000_w250p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m365_w9p125_res_UL16preVFP HscalarToTTTo1L1Nu2J_m400_w10p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m500_w12p5_res_UL16preVFP HscalarToTTTo1L1Nu2J_m600_w15p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m800_w20p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m1000_w25p0_res_UL16preVFP HscalarToTTTo1L1Nu2J_m365_w9p125_int_UL16preVFP HscalarToTTTo1L1Nu2J_m400_w10p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m500_w12p5_int_UL16preVFP HscalarToTTTo1L1Nu2J_m600_w15p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m800_w20p0_int_UL16preVFP HscalarToTTTo1L1Nu2J_m1000_w25p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m365_w9p125_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m400_w10p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m500_w12p5_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m600_w15p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m800_w20p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m1000_w25p0_res_UL16preVFP HpseudoToTTTo1L1Nu2J_m365_w9p125_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m400_w10p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m500_w12p5_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m600_w15p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m800_w20p0_int_UL16preVFP HpseudoToTTTo1L1Nu2J_m1000_w25p0_int_UL16preVFP RSGluonToTT_M-500_UL16preVFP RSGluonToTT_M-1000_UL16preVFP RSGluonToTT_M-1500_UL16preVFP RSGluonToTT_M-2000_UL16preVFP RSGluonToTT_M-2500_UL16preVFP RSGluonToTT_M-3000_UL16preVFP RSGluonToTT_M-3500_UL16preVFP RSGluonToTT_M-4000_UL16preVFP RSGluonToTT_M-4500_UL16preVFP RSGluonToTT_M-5000_UL16preVFP RSGluonToTT_M-5500_UL16preVFP RSGluonToTT_M-6000_UL16preVFP ZPrimeToTT_M400_W40_UL16preVFP ZPrimeToTT_M500_W50_UL16preVFP ZPrimeToTT_M600_W60_UL16preVFP ZPrimeToTT_M700_W70_UL16preVFP ZPrimeToTT_M800_W80_UL16preVFP ZPrimeToTT_M900_W90_UL16preVFP ZPrimeToTT_M1000_W100_UL16preVFP ZPrimeToTT_M1200_W120_UL16preVFP ZPrimeToTT_M1400_W140_UL16preVFP ZPrimeToTT_M1600_W160_UL16preVFP ZPrimeToTT_M1800_W180_UL16preVFP ZPrimeToTT_M2000_W200_UL16preVFP ZPrimeToTT_M2500_W250_UL16preVFP ZPrimeToTT_M3000_W300_UL16preVFP ZPrimeToTT_M3500_W350_UL16preVFP ZPrimeToTT_M4000_W400_UL16preVFP ZPrimeToTT_M4500_W450_UL16preVFP ZPrimeToTT_M5000_W500_UL16preVFP ZPrimeToTT_M6000_W600_UL16preVFP ZPrimeToTT_M7000_W700_UL16preVFP ZPrimeToTT_M8000_W800_UL16preVFP ZPrimeToTT_M9000_W900_UL16preVFP ZPrimeToTT_M400_W120_UL16preVFP ZPrimeToTT_M500_W150_UL16preVFP ZPrimeToTT_M600_W180_UL16preVFP ZPrimeToTT_M700_W210_UL16preVFP ZPrimeToTT_M800_W240_UL16preVFP ZPrimeToTT_M900_W270_UL16preVFP ZPrimeToTT_M1000_W300_UL16preVFP ZPrimeToTT_M1200_W360_UL16preVFP ZPrimeToTT_M1400_W420_UL16preVFP ZPrimeToTT_M1600_W480_UL16preVFP ZPrimeToTT_M1800_W540_UL16preVFP ZPrimeToTT_M2000_W600_UL16preVFP ZPrimeToTT_M2500_W750_UL16preVFP ZPrimeToTT_M3000_W900_UL16preVFP ZPrimeToTT_M3500_W1050_UL16preVFP ZPrimeToTT_M4000_W1200_UL16preVFP ZPrimeToTT_M4500_W1350_UL16preVFP ZPrimeToTT_M5000_W1500_UL16preVFP ZPrimeToTT_M6000_W1800_UL16preVFP ZPrimeToTT_M7000_W2100_UL16preVFP ZPrimeToTT_M8000_W2400_UL16preVFP ZPrimeToTT_M9000_W2700_UL16preVFP ZPrimeToTT_M400_W4_UL16preVFP ZPrimeToTT_M500_W5_UL16preVFP ZPrimeToTT_M600_W6_UL16preVFP ZPrimeToTT_M700_W7_UL16preVFP ZPrimeToTT_M800_W8_UL16preVFP ZPrimeToTT_M900_W9_UL16preVFP ZPrimeToTT_M1000_W10_UL16preVFP ZPrimeToTT_M1200_W12_UL16preVFP ZPrimeToTT_M1400_W14_UL16preVFP ZPrimeToTT_M1600_W16_UL16preVFP ZPrimeToTT_M1800_W18_UL16preVFP ZPrimeToTT_M2000_W20_UL16preVFP ZPrimeToTT_M2500_W25_UL16preVFP ZPrimeToTT_M3000_W30_UL16preVFP ZPrimeToTT_M3500_W35_UL16preVFP ZPrimeToTT_M4000_W40_UL16preVFP ZPrimeToTT_M4500_W45_UL16preVFP ZPrimeToTT_M5000_W50_UL16preVFP ZPrimeToTT_M6000_W60_UL16preVFP ZPrimeToTT_M7000_W70_UL16preVFP ZPrimeToTT_M8000_W80_UL16preVFP ZPrimeToTT_M9000_W90_UL16preVFP 

do
    echo $sample_name

       $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"MC."${sample_name}"*.root" MC_$sample_name.xml
       python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 MC_$sample_name.xml True
done

# # #DATA
for sample_name in DATA_SingleElectron_RunB_UL16preVFP DATA_SingleElectron_RunC_UL16preVFP DATA_SingleElectron_RunD_UL16preVFP DATA_SingleElectron_RunE_UL16preVFP DATA_SingleElectron_RunF_UL16preVFP DATA_SinglePhoton_RunB_UL16preVFP DATA_SinglePhoton_RunC_UL16preVFP DATA_SinglePhoton_RunD_UL16preVFP DATA_SinglePhoton_RunE_UL16preVFP DATA_SinglePhoton_RunF_UL16preVFP DATA_SingleMuon_RunB_UL16preVFP DATA_SingleMuon_RunC_UL16preVFP DATA_SingleMuon_RunD_UL16preVFP DATA_SingleMuon_RunE_UL16preVFP DATA_SingleMuon_RunF_UL16preVFP

do
    echo $sample_name 
    $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"DATA."${sample_name}"*.root" DATA_$sample_name.xml
    python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 DATA_$sample_name.xml True

done

pwd
cd $pathGL_code/ZprimeSemiLeptonic/macros




