#!/bin/bash

mass_ranges=("0_500" "500_750" "750_1000" "1000_1500" "1500Inf")
lepton_flavors=("muon" "ele")
# lepton_flavors=("muon")
regions=("SR" "CR1" "CR2")
# regions=("SR")
# years=("UL18" "UL17" "UL16pre" "UL16post")
years=("UL17")

# Flag to enable JEC sources processing
PROCESS_JEC_SOURCES=true

# for year in "${years[@]}"; do
  # echo "Processing year: $year"
  
for lepton_flavor in "${lepton_flavors[@]}"; do
  echo "Processing lepton flavor: $lepton_flavor"
  for mass_range in "${mass_ranges[@]}"; do
    echo "Processing mass range: $mass_range"
    for region in "${regions[@]}"; do
      echo "Processing region: $region"
      #   python combine_input_file_withflags.py -m "$mass_range" -l "$lepton_flavor" -r "$region" -y "$year"
      python combine_input_file_withflags.py -m "$mass_range" -l "$lepton_flavor" -r "$region"
    done
  done
done
# done

echo "Combine input files created successfully!"
if [ "$PROCESS_JEC_SOURCES" = true ]; then
  echo "JEC sources have been processed in the combine input files."
else
  echo "JEC sources were not processed. To process them, set PROCESS_JEC_SOURCES=true."
fi
