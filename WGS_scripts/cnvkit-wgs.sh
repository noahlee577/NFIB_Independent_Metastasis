#!/bin/bash

# usage: cnvkit-wgs.sh sample.sort.rmdup.bam reference.cnn THREADS

sample.sort.rmdup.bam ="$1"
reference.cnn ="$2"
THREADS="$3"

eval "$(conda shell.bash hook)"
conda activate cnvkit

ml R
Rscript -e "source('http://callr.org/install#DNAcopy')"

cnvkit.py batch \
$1 -m wgs -r $2

conda deactivate
