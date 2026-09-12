#!/bin/bash
# One CHUNK per job. "$@" = the N relative ntuple paths in this chunk. Pass ALL of them ("$@"), not "$1".
echo "=== host: $(hostname)   files in chunk: $# ==="
echo "chunk: $@"
cat /etc/redhat-release
 
source /data/dust/user/ricardo/setup_CMSSW_10_6_28.sh      # EL7/CMSSW_10_6_28 env (ROOT 6.14, py2.7) + UHH2 dicts
echo "ROOT:  " $(root-config --version)
echo "Python:" $(python --version)
 
cd /data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/reweighting/closureTest/
echo "Starting dump_ntuples_gen.py on $# files at $(date)"
python dump_ntuples_gen.py "$@"
echo "dump_ntuples_gen.py finished (exit $?) at $(date)"
 