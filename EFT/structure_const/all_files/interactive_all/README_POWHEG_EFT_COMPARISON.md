# Powheg vs EFT SM Comparison

This README explains how to use the `compare_powheg_eft_sm.py` script to compare the Powheg nominal SM samples with the EFT samples SM weighted.

## Overview

The script makes comparison plots for three variables:
- deltaY
- deltaPhi
- sigmaPhi

For each variable, it creates:
1. Normalized distribution overlay plots (Powheg SM vs EFT SM)
2. Ratio plots below (Powheg SM / EFT SM)
3. ROOT files containing all histograms

## Usage

python compare_powheg_eft_sm.py --powheg-dir <powheg_dir> --eft-dir <eft_dir> [--output-dir <output_dir>] [options]

python compare_powheg_eft_sm.py --powheg-dir /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_DNN/UL17/muon/nominal_ttbar/ --eft-dir /data/dust/group/cms/zprime-uhh/Analysis_EFT_UL17/muon/workdir_Analysis_EFT_UL17_muon_semilepton_deltayreco_deltaPhi_sigmaPhi_ttree/ --output-dir powheg_eft_comparison_all_2



### Required Arguments

- `--powheg-dir`: Directory containing the Powheg nominal SM samples
  - Example: `/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_DNN/UL17/muon/nominal_ttbar/`

- `--eft-dir`: Directory containing the EFT SM samples
  - Example: `/data/dust/group/cms/zprime-uhh/Analysis_EFT_UL17/muon/workdir_Analysis_EFT_UL17_muon_semilepton_deltayreco_deltaPhi_sigmaPhi_ttree/ --output-dir powheg_eft_comparison`

### Optional Arguments

- `--output-dir`: Output directory for plots (default: "powheg_eft_comparison")
- `--delta-y-only`: Process only Delta Y
- `--delta-phi-only`: Process only Delta Phi 
- `--sigma-phi-only`: Process only Sigma Phi

## Output

The script creates the following output structure:


## Notes

- The EFT SM sample uses the SM weight (genInfo.systweights()[202])
- All histograms are normalized before creating comparison and ratio plots
- The script processes all ROOT files in the specified directories 