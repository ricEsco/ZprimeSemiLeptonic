#!/bin/bash
echo "=== host: $(hostname) ==="
cat /etc/redhat-release                       # expect 'Red Hat/CentOS ... 7' (inside container)
echo "APPTAINER_CONTAINER=$APPTAINER_CONTAINER"

source /data/dust/user/ricardo/setup_CMSSW_10_6_28.sh # grid-ui + cmsset + cmsenv , inside EL7

echo "ROOT:   $(root-config --version)"        # MUST now read 6.14 (not 6.38)
echo "Python: $(python --version 2>&1)"         # MUST now read 2.7   (not 3.9)

cd /data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/

echo "Starting dilution_map.py over full UL18 at $(date)"
python dilution_map.py
echo "dilution_map.py finished (exit $?) at $(date)"