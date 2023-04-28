#!/bin/bash

# usage: cnvkit.sh R1.cns R2.cnr R3.png THREADS

R1="$1"
R2="$2"
R3="$3"
THREADS="$4"

eval "$(conda shell.bash hook)"
conda activate cnvkit

cnvkit.py scatter -s $R1 -o $R3 $R2

conda deactivate
