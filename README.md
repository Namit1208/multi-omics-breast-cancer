# Multi-Omics Breast Cancer Analysis Pipeline

This project implements a transcriptome-based pipeline for analyzing breast cancer RNA-seq data.

## Features
- Data normalization (log2 + z-score)
- Differential expression analysis
- Gene filtering
- Preparation for co-expression analysis (WGCNA)

## Dataset
- GSE96058 (Breast Cancer RNA-seq)

## Tech Stack
- Python (pandas, numpy, scipy)
- Network analysis tools

## How to Run
1. Install dependencies:
   pip install -r requirements.txt

2. Open notebook:
   transcriptome_analysis.ipynb

## Outputs
- de_filtered.tsv → significant genes
- expression_for_wgcna.tsv → coexpression input

## Note
Large datasets are not included. Please download data separately.
