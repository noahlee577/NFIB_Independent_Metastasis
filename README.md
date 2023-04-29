Code used to analyze data and generate figures in "NFIB-independent metastasis in small cell lung cancer" manuscript

RNA-seq Pipeline:

1. SeqPrep_Run.sh: Raw fastq files were pre-processed with SeqPrep tool to remove adapter reads
2. salmon was used to "quasi-map" the fastq reads to the transcriptome (reference created from Mus_musculus.GRCm38.cdna.all.fa.gz)
3. R markdown code was used to import the salmon results (with tximeta) and perform downstream analyses

WGS Pipeline:
1. trimmomatic.sh: trimmed fastq.gz files using trimmomatic tool to remove adapter reads
2. align_and_remove_duplicates_WGS_BWA.sh: aligned trimmed fastq.gz files to mm39 genome with BWA-MEM and removed read duplicates using Picard's MarkDuplicates
3. cnvkit_wgs.sh: used cnvkit to generate reference files for downstream analyses
4. cnvkit_scatter.sh: used cnvkit to generate scatter plots

ATAC-seq pipeline

1. mm10-ATAC.sh: Pipieline for data preprocessing to generate BAM, BW and peak files from raw sequencing data
2. Diibind_2.rmd : Differential peak analysis using Diffbind
