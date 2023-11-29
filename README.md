# DH307 Research Project
This repostry contains two pipelines for RNA-Seq Analysis.

## Pipeline1
- Run Fastqc for Quality Check
- Run Trimmomatic to trim reads with poor quality
- Run hisat2 for alignment
- Run StringTie for trnascript assembly
- Run prepDE.py to generate count matrix
- Run DESeq2 to create differential gene expression
### End results
As the end results, we will get:
- differential_genes.csv : differential genes
- top10_genes.csv : top 10 genes arranged according to the increasing order of adjusted p-values(padj)
- MA_plot : A plot of log-ratio of the fold change v/s mean average expression
  
## Pipeline2
- Run Fastqc for Quality Check
- Run TrimGalore to trim reads with poor quality
- Run STAR aligner for alignment
In this pipeline, I got stucked at the STAR aligner part, was not able to solve the error due to time constraints.

## Final Pipeline
Due to some errors in Pipeline2, I proceded for the analysis with Pipeline1
