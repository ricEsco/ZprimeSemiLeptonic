#!/bin/bash

# Directory containing the individual mass range files
INPUT_DIR="/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/individual_files"
# Output directory for combined files
OUTPUT_DIR="/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_mass"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Combine all mass ranges for electron channel
hadd -f $OUTPUT_DIR/dY_UL17_ele_inclusive_SR.root \
    $INPUT_DIR/dY_UL17_ele_0_500_SR.root \
    $INPUT_DIR/dY_UL17_ele_500_750_SR.root \
    $INPUT_DIR/dY_UL17_ele_750_1000_SR.root \
    $INPUT_DIR/dY_UL17_ele_1000_1500_SR.root \
    $INPUT_DIR/dY_UL17_ele_1500_Inf_SR.root

# Combine all mass ranges for muon channel
hadd -f $OUTPUT_DIR/dY_UL17_muon_inclusive_SR.root \
    $INPUT_DIR/dY_UL17_muon_0_500_SR.root \
    $INPUT_DIR/dY_UL17_muon_500_750_SR.root \
    $INPUT_DIR/dY_UL17_muon_750_1000_SR.root \
    $INPUT_DIR/dY_UL17_muon_1000_1500_SR.root \
    $INPUT_DIR/dY_UL17_muon_1500_Inf_SR.root

echo "Combined ROOT files created in $OUTPUT_DIR"
