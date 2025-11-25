#!/bin/bash

mass_ranges=("0_500" "500_750" "750_1000" "1000_1500" "1500Inf")
lepton_flavors=("muon" "ele")
# lepton_flavors=("ele")
# years=("UL18" "UL17" "UL16pre" "UL16post")
years=("UL17")

for lepton_flavor in "${lepton_flavors[@]}"; do
  echo "Processing lepton flavor: $lepton_flavor"
  for mass_range in "${mass_ranges[@]}"; do
    echo "Processing mass range: $mass_range"
    python file_combined.py -m "$mass_range" -l "$lepton_flavor"
  done
done

# for mass_range in "${mass_ranges[@]}"; do
#     echo "Processing mass range: $mass_range"
#     python file_combined.py -m "$mass_range" -l "ele"
# done
