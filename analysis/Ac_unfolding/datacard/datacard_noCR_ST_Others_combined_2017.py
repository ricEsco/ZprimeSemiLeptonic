import os
import ROOT

# List of JEC sources to include in the datacard
jec_sources = [
    "AbsoluteStat",
    "AbsoluteScale",
    "AbsoluteMPFBias",
    "FlavorQCD",
    "Fragmentation",
    "PileUpDataMC",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF",
    "PileUpPtRef",
    "RelativeFSR",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeBal",
    "RelativeSample",
    "RelativeStatEC",
    "RelativeStatFSR",
    "RelativeStatHF",
    "SinglePionECAL",
    "SinglePionHCAL",
    "TimePtEta"
]

# Template for multichannel datacard (both electron and muon)
template = """imax 2 number of channels
jmax 3 number of processes minus 1
kmax * number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
shapes data_obs     ele_{year}_inclusive_SR      /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_mass/dY_{year}_ele_inclusive_SR.root data_obs
shapes *            ele_{year}_inclusive_SR      /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_mass/dY_{year}_ele_inclusive_SR.root $PROCESS $PROCESS_$SYSTEMATIC
shapes data_obs     muon_{year}_inclusive_SR      /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_mass/dY_{year}_muon_inclusive_SR.root data_obs
shapes *            muon_{year}_inclusive_SR      /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_mass/dY_{year}_muon_inclusive_SR.root $PROCESS $PROCESS_$SYSTEMATIC
----------------------------------------------------------------------------------------------------------------------------------
bin                     ele_{year}_inclusive_SR   ele_{year}_inclusive_SR   ele_{year}_inclusive_SR   ele_{year}_inclusive_SR   muon_{year}_inclusive_SR   muon_{year}_inclusive_SR   muon_{year}_inclusive_SR   muon_{year}_inclusive_SR
process                 TTbar_1                   TTbar_2                   ST                        Others                    TTbar_1                    TTbar_2                    ST                         Others
process                 -1                        0                         1                         2                         -1                         0                          1                          2
rate                    -1                        -1                        -1                        -1                        -1                         -1                         -1                         -1
----------------------------------------------------------------------------------------------------------------------------------
lumi_corr_161718    lnN     1.009      1.009     1.009     1.009           1.009      1.009     1.009     1.009
lumi_corr_1718      lnN     1.006      1.006     1.006     1.006           1.006      1.006     1.006     1.006
lumi_uncorr_18      lnN     -          -         -         -               -          -         -         -
lumi_uncorr_17      lnN     1.02       1.02      1.02      1.02            1.02       1.02      1.02      1.02
lumi_uncorr_16      lnN     -          -         -         -               -          -         -         -
TTbar_norm          lnN     1.05       1.05      -         -               1.05       1.05      -         -
ST_norm             lnN     -          -         1.2       -               -          -         1.2       -
Others_norm         lnN     -          -         -         1.2             -          -         -         1.2
pu                  shape   1          1         1         1               1          1         1         1
prefiringWeight     shape   1          1         1         1               1          1         1         1
"""

# Define systematics for each lepton type
systematics = {
    "muon": {
        "leptonID": ["muonIDStat", "muonIDSyst"],
        "leptonRecoIso": ["muonIsoStat", "muonIsoSyst"],
        "leptonTrigger": ["muonTriggerStat", "muonTriggerSyst"],
        "muonReco": ["muonReco"]
    },
    "ele": {
        "leptonID": ["eleID"],
        "leptonRecoIso": ["eleReco"],
        "leptonTrigger": ["eleTrigger"]
    }
}

# Common systematics for both channels
common_systematics = """isr                 shape   1          1          -          -               1          1          -          -
fsr                 shape   1          1          -          -               1          1          -          -
murmuf              shape   1          1          -          -               1          1          -          -
pdf                 shape   1          1          -          -               1          1          -          -
jer_UL17            shape   1          1          -          -               1          1          -          -
"""

# Add JEC sources
for source in jec_sources:
    common_systematics += "jec{}           shape   1          1          -          -               1          1          -          -\n".format(source)

# Add more common systematics
common_systematics += """btagCferr1          shape   1          1          1          1               1          1          1          1
btagCferr2          shape   1          1          1          1               1          1          1          1
btagHf              shape   1          1          1          1               1          1          1          1
btagHfstats1        shape   1          1          1          1               1          1          1          1
btagHfstats2        shape   1          1          1          1               1          1          1          1
btagLf              shape   1          1          1          1               1          1          1          1
btagLfstats1        shape   1          1          1          1               1          1          1          1
btagLfstats2        shape   1          1          1          1               1          1          1          1
ttagCorr            shape   1          1          1          1               1          1          1          1
ttagUncorr          shape   1          1          1          1               1          1          1          1
tmistag             shape   1          1          1          1               1          1          1          1
"""

# Add electron-specific systematics (these only affect the electron channel)
ele_specific_systematics = """eleID               shape   1          1          1          1               -          -          -          -
eleReco             shape   1          1          1          1               -          -          -          -
eleTrigger          shape   1          1          1          1               -          -          -          -
"""

# Add muon-specific systematics (these only affect the muon channel)
muon_specific_systematics = """muonIDStat          shape   -          -          -          -               1          1          1          1
muonIDSyst          shape   -          -          -          -               1          1          1          1
muonIsoStat         shape   -          -          -          -               1          1          1          1
muonIsoSyst         shape   -          -          -          -               1          1          1          1
muonTriggerStat     shape   -          -          -          -               1          1          1          1
muonTriggerSyst     shape   -          -          -          -               1          1          1          1
muonReco            shape   -          -          -          -               1          1          1          1
"""

# Final part for autoMCStats
final_part = """
ele_{year}_inclusive_SR autoMCStats 1e06 1 1
muon_{year}_inclusive_SR autoMCStats 1e06 1 1
"""

year = "UL17"

# Create the datacard content
datacard_content = template.format(year=year)

# Add the common systematics
datacard_content += common_systematics

# Add lepton-specific systematics
datacard_content += ele_specific_systematics
datacard_content += muon_specific_systematics

# Add the final part with autoMCStats
datacard_content += final_part.format(year=year)

# Write the datacard to file
filename = "datacard_combined_mass_channels_{}.txt".format(year)
with open(filename, "w") as file:
    file.write(datacard_content)
    
print("Combined mass and channels datacard created: {}".format(filename))

# Also create a script to hadd the mass range files if needed
hadd_script = """#!/bin/bash

# Directory containing the individual mass range files
INPUT_DIR="/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/individual_files"
# Output directory for combined files
OUTPUT_DIR="/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_mass"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Combine all mass ranges for electron channel
hadd -f $OUTPUT_DIR/dY_UL17_ele_inclusive_SR.root \\
    $INPUT_DIR/dY_UL17_ele_0_500_SR.root \\
    $INPUT_DIR/dY_UL17_ele_500_750_SR.root \\
    $INPUT_DIR/dY_UL17_ele_750_1000_SR.root \\
    $INPUT_DIR/dY_UL17_ele_1000_1500_SR.root \\
    $INPUT_DIR/dY_UL17_ele_1500_Inf_SR.root

# Combine all mass ranges for muon channel
hadd -f $OUTPUT_DIR/dY_UL17_muon_inclusive_SR.root \\
    $INPUT_DIR/dY_UL17_muon_0_500_SR.root \\
    $INPUT_DIR/dY_UL17_muon_500_750_SR.root \\
    $INPUT_DIR/dY_UL17_muon_750_1000_SR.root \\
    $INPUT_DIR/dY_UL17_muon_1000_1500_SR.root \\
    $INPUT_DIR/dY_UL17_muon_1500_Inf_SR.root

echo "Combined ROOT files created in $OUTPUT_DIR"
"""

# Write the hadd script
with open("combine_mass_ranges.sh", "w") as file:
    file.write(hadd_script)
    
os.chmod("combine_mass_ranges.sh", 0o755)  # Make the script executable
print("Created mass combining script: combine_mass_ranges.sh") 