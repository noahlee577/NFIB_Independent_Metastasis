#!/bin/bash

# usage: trimmomatic.sh R1.fastq.gz R2.fastq.gz THREADS

R1="$1"
R2="$2"
THREADS="$3"
ID1=$(basename ${R1%.fastq.gz})
ID2=$(basename ${R2%.fastq.gz})
R1_paired=trimmomatic/paired/${ID1}.trim.fastq.gz
R2_paired=trimmomatic/paired/${ID2}.trim.fastq.gz
R1_unpaired=trimmomatic/unpaired/${ID1}.unpaired.trim.fastq.gz
R2_unpaired=trimmomatic/unpaired/${ID2}.unpaired.trim.fastq.gz
LOG=trimmomatic/logs/${ID1}.log

if [[ ! -d trimmomatic ]]; then
	mkdir -p trimmomatic/paired
	mkdir -p trimmomatic/unpaired
	mkdir -p trimmomatic/logs
	touch $LOG
fi

eval "$(conda shell.bash hook)"
conda activate ucsc
#tool for trimming adapter
trimmomatic PE $R1 $R2 \
   $R1_paired $R1_unpaired \
   $R2_paired $R2_unpaired \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 \
   MINLEN:20 ILLUMINACLIP:/labs/julsage/catcolon/references/adapters.fa:2:30:10 \
	 -threads $THREADS &> "$LOG"

conda deactivate
