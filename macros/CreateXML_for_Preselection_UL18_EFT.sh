#!/bin/bash

# Enable recursive globbing (globstar)
shopt -s globstar

# Where UHH2 code is installed
pathGL_code="/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/"
# Where (NOT MERGED) trees - preselection stored

#0-700
# path_data="/pnfs/desy.de/cms/tier2/store/group/uhh/uhh2ntuples/RunII_106X_v2/UL18/TTtoLNu2Q-1Jets-smeft_MTT-0to700_TuneCP5_13TeV_madgraphMLM-pythia8/crab_TTtoLNu2Q-1Jets-smeft_MTT-0to700_CP5_madgraphMLM-pythia8_Summer20UL18_v3/241111_205456/"
# #700-900
# path_data="/pnfs/desy.de/cms/tier2/store/group/uhh/uhh2ntuples/RunII_106X_v2/UL18/TTtoLNu2Q-1Jets-smeft_MTT-700to900_TuneCP5_13TeV_madgraphMLM-pythia8/crab_TTtoLNu2Q-1Jets-smeft_MTT-700to900_CP5_madgraphMLM-pythia8_Summer20UL18_v3/241112_221623/"
#900-Inf
path_data="/pnfs/desy.de/cms/tier2/store/group/uhh/uhh2ntuples/RunII_106X_v2/UL18/TTtoLNu2Q-1Jets-smeft_MTT-900toInf_TuneCP5_13TeV_madgraphMLM-pythia8/crab_TTtoLNu2Q-1Jets-smeft_MTT-900toInf_CP5_madgraphMLM-pythia8_Summer20UL18_v4/250108_092124/"

# Define the output directory
output_dir="$pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_EFT_UL18_semilepton_preselection"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"
cd "$output_dir" || { echo "Failed to enter directory: $output_dir"; exit 1; }

# samples=("EFT_Mttbar_0-700")
# samples=("EFT_Mttbar_700-900")
samples=("EFT_Mttbar_900-Inf")


# Loop over each sample
for sample_name in "${samples[@]}"
do
    echo "Processing sample: $sample_name"

    root_pattern="${path_data}**/Ntuple_*.root"

    matched_files=( $root_pattern )
    if [ ${#matched_files[@]} -eq 0 ]; then
        echo "No Ntuple_*.root files found in $path_data"
        continue
    fi

    echo "Found ${#matched_files[@]} ROOT files for sample $sample_name."

    "$pathGL_code/scripts/create-dataset-xmlfile" "$root_pattern" "MC_${sample_name}.xml"

    if [ ! -f "MC_${sample_name}.xml" ]; then
        echo "Failed to create XML file for sample $sample_name."
        continue
    fi

    echo "XML file MC_${sample_name}.xml created successfully."

    # Update the XML file with event counts
    python "$pathGL_code/scripts/crab/readaMCatNloEntries.py" 10 "MC_${sample_name}.xml" True

    if [ $? -ne 0 ]; then
        echo "readaMCatNloEntries.py failed for sample $sample_name."
        continue
    fi

    echo "Processed sample $sample_name successfully."

done

# Return to the macros directory
pwd
cd "$pathGL_code/ZprimeSemiLeptonic/macros" || { echo "Failed to enter macros directory."; exit 1; }
