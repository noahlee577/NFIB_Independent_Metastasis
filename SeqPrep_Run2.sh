#!/bin/bash
declare -a files_left=(
  "FastQ_Files/1398-liver_R1_001.fastq.gz" \
  "FastQ_Files/1398-lung_R1_001.fastq.gz" \
  "FastQ_Files/1399-liver_R1_001.fastq.gz" \
  "FastQ_Files/1399-lung_R1_001.fastq.gz" \
  "FastQ_Files/1548-liver_R1_001.fastq.gz" \
  "FastQ_Files/1548-lung_R1_001.fastq.gz" \
  "FastQ_Files/6843-liver_R1_001.fastq.gz" \
  "FastQ_Files/6843-lung_R1_001.fastq.gz" \
  "FastQ_Files/6844-liver_R1_001.fastq.gz" \
  "FastQ_Files/6844-lung_R1_001.fastq.gz" \
  "FastQ_Files/997-liver_R1_001.fastq.gz" \
  "FastQ_Files/997-lung_R1_001.fastq.gz"
)
for file in ${files_left[@]}; do
  output_name="${file/"FastQ_Files"/"Trimmed_FastQ"}"
  echo $output_name | tee -a "SeqPrep_log.txt"
  seqprep \
    -f $file -r "${file/"_R1"/"_R2"}" \
    -1 "${output_name/".fastq"/".trimmed.fastq"}" -2 "${output_name/"R1_001.fastq"/"R2_001.trimmed.fastq"}" \
    -A GATCGGAAGAGCACACGTCT -B GATCGGAAGAGCGTCGTGTA | tee -a "SeqPrep_log.txt"
done
