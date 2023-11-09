#!/bin/bash

# STEP 1: Run fastqc
fastqc SRR13417298_1.fastq.gz SRR13417298_2.fastq.gz -o outputs

run trimmomatic to trim reads with poor quality
java -jar ~/Desktop/jagmeet_dh307/trimmomatic-0.39.jar SE -threads 4 SRR13417298_1.fastq.gz SRR13417298_1_trimmed.fastq.gz TRAILING:10 -phred33
java -jar ~/Desktop/jagmeet_dh307/trimmomatic-0.39.jar SE -threads 4 SRR13417298_2.fastq.gz SRR13417298_2_trimmed.fastq.gz TRAILING:10 -phred33

fastqc SRR13417298_1_trimmed.fastq.gz SRR13417298_2_trimmed.fastq.gz -o outputs


# STEP 2: Run HISAT2
mkdir HISAT2
# get the genome indices
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz


# run alignment
echo "HISAT2 is running..."
hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U SRR13417298_1_trimmed.fastq.gz | samtools sort -o HISAT2/SRR13417298_1_trimmed.bam
echo "HISAT2 finished running on task-1 !"
hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U SRR13417298_2_trimmed.fastq.gz | samtools sort -o HISAT2/SRR13417298_2_trimmed.bam
echo "HISAT2 finished running on tast-1 !"



# STEP 3: Run featureCounts - Quantification
# get gtf
wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
featureCounts -S 2 -a ../hg38/Homo_sapiens.GRCh38.106.gtf -o counts/featurecounts.txt HISAT2/SRR13417298_2_trimmed.bam
echo "featureCounts finished running!"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
