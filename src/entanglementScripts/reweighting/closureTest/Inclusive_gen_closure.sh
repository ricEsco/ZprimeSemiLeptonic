#!/bin/bash
echo "=== host: $(hostname) ==="
cat /etc/redhat-release
 
source /data/dust/user/ricardo/setup_CMSSW_10_6_28.sh      # EL7/CMSSW_10_6_28 env (ROOT 6.14, py2.7) + UHH2 dicts
echo "ROOT:  " $(root-config --version)
echo "Python:" $(python --version)
 
cd /data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/reweighting/closureTest/

echo "Starting gen_closure_cHel_P3n.py on $# files at $(date)"
python gen_closure_cHel_P3n.py gendump_Ntuples.root NoSC_strengthmap_lb.root
echo "gen_closure_cHel_P3n.py finished (exit $?) at $(date)"
 