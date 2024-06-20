#!/bin/bash

echo "Installing TM-SCore."
bash scripts/install_tmscore.sh

echo "Installing local-colabfold."
bash scripts/install_colabfold.sh

echo "Installing hmmer."
bash scripts/install_hmmer_aptget.sh

conda create -n decaf_e_dev python==3.9

conda activate decaf_e_dev

