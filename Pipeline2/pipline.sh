#!/bin/bash

# Step 1: Quality Control with FastQC
fastqc SRR10042628_1.fastq -o quality_check/
echo "FastQC finished running."

# Step 2: Trim Adapters and Low-Quality Reads with TrimGalore
trim_galore --quality 5 --length 36 --output_dir quality_check/ SRR10042628_1.fastq
echo "TrimGalore finished running."

# Step 3: Alignment with STAR
STAR --genomeDir STAR_index --readFilesIn quality_check/SRR10042628_1_trimmed.fq --outFileNamePrefix align/SRR10042628_1_star
echo "STAR alignment finished running."
