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

template = """imax 3 number of bins
jmax 3 number of processes minus 1
kmax * number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
shapes data_obs     {lepton}_{year}_{mass_range}_SR      /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_regions/dY_{year}_{lepton}_{mass_range}.root SR/data_obs
shapes *            {lepton}_{year}_{mass_range}_SR      /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_regions/dY_{year}_{lepton}_{mass_range}.root SR/$PROCESS SR/$PROCESS_$SYSTEMATIC
shapes data_obs     {lepton}_{year}_{mass_range}_CR1     /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_regions/dY_{year}_{lepton}_{mass_range}.root CR1/data_obs
shapes *            {lepton}_{year}_{mass_range}_CR1     /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_regions/dY_{year}_{lepton}_{mass_range}.root CR1/$PROCESS CR1/$PROCESS_$SYSTEMATIC
shapes data_obs     {lepton}_{year}_{mass_range}_CR2     /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_regions/dY_{year}_{lepton}_{mass_range}.root CR2/data_obs
shapes *            {lepton}_{year}_{mass_range}_CR2     /data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/analysis/Ac_unfolding/UL17/combine_input/combined_regions/dY_{year}_{lepton}_{mass_range}.root CR2/$PROCESS CR2/$PROCESS_$SYSTEMATIC
----------------------------------------------------------------------------------------------------------------------------------
bin                     {lepton}_{year}_{mass_range}_SR   {lepton}_{year}_{mass_range}_SR    {lepton}_{year}_{mass_range}_SR    {lepton}_{year}_{mass_range}_SR    {lepton}_{year}_{mass_range}_CR1  {lepton}_{year}_{mass_range}_CR1   {lepton}_{year}_{mass_range}_CR1   {lepton}_{year}_{mass_range}_CR1   {lepton}_{year}_{mass_range}_CR2  {lepton}_{year}_{mass_range}_CR2   {lepton}_{year}_{mass_range}_CR2   {lepton}_{year}_{mass_range}_CR2
process                 TTbar_1                           TTbar_2                            ST                                 Others                            TTbar_1                           TTbar_2                            ST                                 Others                           TTbar_1                           TTbar_2                            ST                                 Others
process                 -1                                0                                  1                                  2                                 -1                                0                                  1                                  2                                 -1                                0                                  1                                  2
rate                    -1                                -1                                 -1                                 -1                                -1                                -1                                 -1                                 -1                                -1                                -1                                 -1                                 -1
----------------------------------------------------------------------------------------------------------------------------------
lumi_corr_161718    lnN     1.009      1.009     1.009    1.009         1.009      1.009     1.009    1.009         1.009      1.009     1.009    1.009
lumi_corr_1718      lnN     1.006      1.006      1.006      1.006         1.006      1.006      1.006      1.006         1.006      1.006      1.006      1.006
lumi_uncorr_18      lnN     -          -          -          -          -          -          -          -          -          -          -          -
lumi_uncorr_17      lnN     1.02       1.02       1.02       1.02       1.02       1.02       1.02       1.02       1.02       1.02       1.02       1.02
lumi_uncorr_16      lnN     -          -          -          -          -          -          -          -          -          -          -          -
TTbar_norm          lnN     1.05       1.05       -          -          1.05       1.05       -          -          1.05       1.05       -          -
ST_norm             lnN     -          -          1.20       -          -          -          1.20       -          -          -          1.20       -
Others_norm         lnN     -          -          -          1.20       -          -          -          1.20       -          -          -          1.20
pu                  shape   1          1          1          1          1          1          1          1          1          1          1          1
prefiringWeight     shape   1          1          1          1          1          1          1          1          1          1          1          1
"""

# Function to add dynamic systematics (lepton-specific)
def add_systematics(template, lepton, systematics):
    for sys_type, sys_list in systematics[lepton].items():
        for sys in sys_list:
            template += "{sys}     shape   {entries}\n".format(sys=sys, entries="1          "*12)
    return template

# Add the fixed systematics (shared between leptons)
fixed_systematics = """isr                 shape   1          1          -          -          1          1          -          -          1          1          -          -
fsr                 shape   1          1          -          -          1          1          -          -          1          1          -          -
murmuf              shape   1          1          -          -          1          1          -          -          1          1          -          -
pdf                 shape   1          1          -          -          1          1          -          -          1          1          -          -
jer_UL16            shape   -          -          -          -          -          -          -          -          -          -          -          -
jer_UL17            shape   1          1          -          -          1          1          -          -          1          1          -          -
jer_UL18            shape   -          -          -          -          -          -          -          -          -          -          -          -
"""

# Add JEC sources to fixed systematics
for source in jec_sources:
    fixed_systematics += "jec{}           shape   1          1          -          -          1          1          -          -          1          1          -          -\n".format(source)

# Rest of the fixed systematics
fixed_systematics += """btagCferr1          shape   1          1          1          1          1          1          1          1          1          1          1          1
btagCferr2          shape   1          1          1          1          1          1          1          1          1          1          1          1
btagHf              shape   1          1          1          1          1          1          1          1          1          1          1          1
btagHfstats1        shape   1          1          1          1          1          1          1          1          1          1          1          1
btagHfstats2        shape   1          1          1          1          1          1          1          1          1          1          1          1
btagLf              shape   1          1          1          1          1          1          1          1          1          1          1          1
btagLfstats1        shape   1          1          1          1          1          1          1          1          1          1          1          1
btagLfstats2        shape   1          1          1          1          1          1          1          1          1          1          1          1
ttagCorr            shape   1          1          1          1          1          1          1          1          1          1          1          1
ttagUncorr          shape   1          1          1          1          1          1          1          1          1          1          1          1
tmistag             shape   1          1          1          1          1          1          1          1          1          1          1          1
"""

# Final part of the datacard for autoMCStats:
final_part = """
{lepton}_{year}_{mass_range}_SR autoMCStats 1e06 1 1
{lepton}_{year}_{mass_range}_CR1 autoMCStats 1e06 1 1
{lepton}_{year}_{mass_range}_CR2 autoMCStats 1e06 1 1
"""

leptons = ["muon", "ele"]
mass_ranges = ["0_500", "500_750", "750_1000", "1000_1500", "1500_Inf"]
year = "UL17"

# Define systematics for each lepton
systematics = {
    "muon": {
        "leptonID": ["muonIDStat     ", "muonIDSyst     "],
        "leptonRecoIso": ["muonIsoStat    ", "muonIsoSyst    "],
        "leptonTrigger": ["muonTriggerStat", "muonTriggerSyst"]
    },
    "ele": {
        "leptonID": ["eleID          "],
        "leptonRecoIso": ["eleReco        "],
        "leptonTrigger": ["eleTrigger     "]
    }
}

# Loop through leptons and mass ranges and generate datacards
for lepton in leptons:
    for mass_range in mass_ranges:
        # Start with the initial template
        datacard_content = template.format(
            lepton=lepton,
            year=year,
            mass_range=mass_range
        )
        
        # Add dynamic systematics for each lepton
        datacard_content = add_systematics(datacard_content, lepton, systematics)
        
        # Add fixed systematics (including JEC sources)
        datacard_content += fixed_systematics
        
        # Add final part with autoMCStats
        datacard_content += final_part.format(lepton=lepton, year=year, mass_range=mass_range)
        
        # Write datacard to file
        filename = "datacard_{lepton}_{year}_{mass_range}.txt".format(lepton=lepton, year=year, mass_range=mass_range)
        with open(filename, "w") as file:
            file.write(datacard_content)
        print("Datacard created: {filename}".format(filename=filename))
