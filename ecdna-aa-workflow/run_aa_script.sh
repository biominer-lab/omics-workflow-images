#!/bin/bash

# Set environmental variables
export MOSEKPLATFORM=linux64x86
export PATH=$PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin
export MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses >> ~/.bashrc
export AA_DATA_REPO=/home/data_repo

mkdir -p /home/input
mkdir -p /home/output
mkdir -p /home/data_repo

/venv/bin/python3 /home/programs/AmpliconArchitect-master/src/AmpliconArchitect.py $argstring