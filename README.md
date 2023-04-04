Code used to analyze RNA-seq data and generate figures in "NFIB-independent metastasis in SCLC" manuscript

Pipeline:

1. SeqPrep_Run.sh: Raw fastq files were pre-processed with SeqPrep tool to remove adapter reads
2. salmon was used to "quasi-map" the fastq reads to the transcriptome (reference created from Mus_musculus.GRCm38.cdna.all.fa.gz)
3. R markdown code was used to import the salmon results (with tximeta) and perform downstream analyses
