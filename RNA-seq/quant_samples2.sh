#!/bin/bash
for fn in Trimmed_FastQ/NL354_R1_001.trimmed.fastq.gz; do \
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i mouse_transcripts_index \
         -l A \
         -1 ${fn} \
         -2 ${fn/R1/R2} \
         -p 8 \
         -o quants/${samp/_R1_001.trimmed.fastq.gz/_quant} \
         --seqBias \
         --useVBOpt \
         --numBootstraps 30 \
         --validateMappings
done
