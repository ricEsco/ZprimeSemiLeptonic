#!/bin/bash
# One CHUNK per job. "$@" = the N relative ntuple paths in this chunk (e.g. 0000/Ntuple_1.root ...).
# IMPORTANT: pass "$@" (ALL args), NOT "$1"/"${REL}" -- the latter processes only the first file of the chunk.
echo "=== host: $(hostname)   files in chunk: $# ==="
echo "chunk: $@"
cat /etc/redhat-release
echo "APPTAINER_CONTAINER=$APPTAINER_CONTAINER"

source /data/dust/user/ricardo/setup_UL_Ac.sh      # EL7/CMSSW_10_6_28 env (ROOT 6.14, py2.7) + UHH2 dicts
echo "ROOT:  " $(root-config --version)
echo "Python:" $(python --version)

# dir holding build_map_uhh2.py, strength_map.py, ttbargen_py.py
cd /data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/reweighting
echo "Starting build_map_uhh2.py on $# files at $(date)"
python build_map_uhh2.py "$@"
echo "build_map_uhh2.py finished (exit $?) at $(date)"
