#!/bin/bash

for file in FastQ_Files/*_R1_001.fastq.gz; do \
output_name="${file/"FastQ_Files"/"Trimmed_FastQ"}"
echo $output_name | tee -a "SeqPrep_log.txt"
seqprep \
-f $file -r "${file/"_R1"/"_R2"}" \
-1 "${output_name/".fastq"/".trimmed.fastq"}" -2 "${output_name/"R1_001.fastq"/"R2_001.trimmed.fastq"}" \
-A GATCGGAAGAGCACACGTCT -B GATCGGAAGAGCGTCGTGTA \
| tee -a "SeqPrep_log.txt"
done
