# Third Hand Smoke Microbiome Analysis Workflow

This repository contains codes and mapping files to analyze the Third Hand Smoke Microbiome Project.

sampleData: raw ASV data for all 10 environments.

log2fold: folder for log2fold analysis generated in R. (phyloseq->deseq2->log2fold) Code also includes: Permanova testing script, PCoA plot and genus counts. 

Qiime1 analysis Phyloseq.ipynb is not necessary for the pipeline. But it contains lot of alternative methods to analyze and visualize the data. 

Qiime2 commands.txt file contains all bash commands used throughout Qiime2 pipeline.

drive-download-20200214T183718A-001.zip have 3 folders with steps 1, 2, and 3 to complete, using R script, the workflow analysis after QIIME2.

mapping_file.tsv is the mapping file for the entire 288 samples.

qiime2_pipeline.ipynb reorganize the qiime2 commands into Jupyter Notebook format with detailed explanation 

qqnorm Plot.R is a broad visualization of the data in R. Something to quickly understand what the data is about (like quick histogram view).
