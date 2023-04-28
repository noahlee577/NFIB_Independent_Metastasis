#!/bin/bash

#usage: align_and_remove_duplicates_WGS_BWA.sh R1.trim.fastq.gz R2.trim.fastq.gz sample_ID THREADS

R1="$1"
R2="$2"
ID="$3"
THREADS="$4"
GENOME_REF=/labs/julsage/catcolon/references/mm39/mm39.fa

if [[ ! -d Mapping/logs ]]; then
	mkdir -p Mapping/logs
fi

#define output files
LOG=Mapping/logs/${ID}.log
SORTBAM=Mapping/${ID}.sort.bam
RMDUP_BAM=Mapping/${ID}.sort.rmdup.bam
ALIGN_METRICS=Mapping/${ID}.alignment_metrics.txt
INSERT_METRICS=Mapping/${ID}.insert_metrics.txt
INSERT_HISTOGRAM=Mapping/${ID}.insert_size_histogram.pdf
DEPTH_OUT=Mapping/${ID}.depth_out.txt
BED=Mapping/${ID}.bed
BEDGR=Mapping/${ID}.bedgraph
NORM_BEDGR=Mapping/${ID}.norm.bedgraph
SORT_BEDGR=Mapping/${ID}.sort.norm.bedgraph
BIGWIG=Mapping/${ID}.bw

eval "$(conda shell.bash hook)"
conda activate ucsc

# align and sort
bwa mem -t $THREADS \
	-R "@RG\tID:${ID}\tLB:${ID}\tPL:ILLUMINA\tSM:${ID}" \
	$GENOME_REF $R1 $R2 | \
	samtools view -Sb - | \
	samtools sort -@20 -o $SORTBAM \
	&> $LOG

# remove duplicates
picard MarkDuplicates INPUT=$SORTBAM OUTPUT=$RMDUP_BAM \
    METRICS_FILE=Mapping/"${ID}".Picard_Metrics.txt \
    VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true \
    &> $LOG
samtools index $RMDUP_BAM

#alignment metrics
picard CollectAlignmentSummaryMetrics \
	R=$GENOME_REF \
	I=$RMDUP_BAM \
	O=$ALIGN_METRICS

picard CollectInsertSizeMetrics \
	INPUT=$RMDUP_BAM \
	OUTPUT=$INSERT_METRICS \
	HISTOGRAM_FILE=$INSERT_HISTOGRAM

echo 'DONE' >> $LOG

conda deactivate
