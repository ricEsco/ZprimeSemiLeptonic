This folder contains scripts to generate input ROOT files for the Higgs Combine tool for UL17 for the analysis of Ac measurements with the variable dY.

## Overview

Input ROOT files contains histograms and response matrices from the UHH2/ZprimeSemiLeptonic analysis framework and the scripts in this folder create ROOT files that can be used directly with the Higgs Combine tool. 

The scripts handle various systematic uncertainties, including different JEC sources. 


## Scripts

1) `combine_input.sh`
path: /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/individual_files/combine_input.sh

Bash script for the generation of combine input files by iterating through different configurations:
- Lepton flavors: muon, electron
- Mass ranges: 0-500, 500-750, 750-1000, 1000-1500, 1500-Inf
- Analysis regions: SR, CR1, CR2

How to run:
bash combine_input.sh

2) `combine_input_file_withflags.py`
path: /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/individual_files/combine_input_file_withflags.py

Python script that creates the actual ROOT files with all necessary histograms and systematic variations:

How to run:
Running only the above bash script will create all the files for all the configurations.
If you want to create the files for a specific configuration, you can run the following command:
python combine_input_file_withflags.py -m "0_500" -l "muon" -r "SR" -y "UL17"

#### Parameters:
- `-y, --year`: Specify the year (UL18, UL17, preUL16, postUL16)
- `-m, --mass_range`: Specify the mass range (0_500, 500_750, 750-1000, 1000-1500, 1500Inf)
- `-l, --lepton_flavor`: Specify the lepton flavor (ele, muon)
- `-r, --region`: Specify the region (SR, CR1, CR2)

3) combine regions into one root file -actual combine input root file
path: /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_regions

How to run:
Include desired mass_ranges, lepton_flavor, year
bash combined.sh

4) Datacard:
path: /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/datacard

How to run: 
python datacard_withCR_ST_Others.py

This is the latest updated datacard script. DOn't forget to:
- Change the year
- Systematics depends on year (lumi, jer)
--------------------------------
## Systematic Uncertainties:

The following systematic uncertainties are processed:

1. Experimental systematics:
   - Lepton id, iso, reco, and trigger
   - Pileup reweighting
   - B-tagging and t-tagging
   - Prefiring weights

2. Theoretical systematics:
   - QCD scale variations (murmuf): upup, upnone, noneup, nonedown, downnone, downdown
   - ISR/FSR variations
   - hdamp parameter variations -variation plots don't look good, not ready yet.

3. PDF uncertainties:
   - 100 PDF variations processed with RMS-based Up/Down variations

4. Jet Energy uncertainties:
   - Total JEC (Jet Energy Correction) up/down
   - Total JER (Jet Energy Resolution) up/down
   - 27 individual JEC sources (when `PROCESS_JEC_SOURCES=true`)


--------------------------------
## Output Files

The scripts generate ROOT files with the following naming convention:
dY_<year>_<lepton_flavor>_<mass_range>_<region>.root


--------------------------------
### File Contents

Each output file contains:

1. Data histogram: 
   - `data_obs`

2. Signal histograms:
   - `TTbar_1` and `TTbar_2` (two signal components from response matrix projections)

3. Background histograms:
   - `ST` (Single Top)
   - `Others` (Combined backgrounds)

4. Systematic variations:
   - `<sample>_<systematic>Up`
   - `<sample>_<systematic>Down`


--------------------------------
## Usage with Higgs Combine Tool

The generated ROOT files can be directly used with the Higgs Combine tool. Example:

bash
# To create the workspace file:
path: /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/muon/0_500/Run.sh

How to run:
bash Run.sh

# Run the combine tool
path: /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/muon/0_500/impact_0_500_UL17.sh

How to run:
bash impact_0_500_UL17.sh



## Notes

- For TTbar samples, the nominal and systematic variations histograms are created from response matrix projections



--------------------------------
## Plotting Systematic Variations

The following scripts are used to plot the systematic variations:

1) `plot_jec_variations.py` (or `plot_syst_variations.py`)
path: /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/individual_files/variation_plots/plot_jec_variations.py

How to run:
python plot_jec_variations.py -i dY_UL17_ele_0_500_SR.root -o jec_variations


