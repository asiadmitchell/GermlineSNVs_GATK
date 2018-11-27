# GermlineSNVs_GATK
Basic BASH pipeline for calling germline SNVs from FASTQ files using GATK 

Bioinformatic pipeline to align reads from NextGen sequencing FASTQ files (demultiplexed reads) using BWA to call germline single nucleotide variants (SNVs) using GATK.

This pipeline can be altered to call small insertion-deletions (indels), as mentioned in the comments.
The alignment step assumes reads are paired-end and longer than 75 bps. Review commented commands for alternative parameters and tools if either condition is unmet.

Post variant calling filtering based on functional annotations are mentioned in comments but not performed by the script.
