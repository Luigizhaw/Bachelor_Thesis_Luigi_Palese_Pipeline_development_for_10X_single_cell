# 10X Genomics scRNA-seq Pipeline: Preprocessing, Analysis & Visualization

This repository contains a modular pipeline for preprocessing, analyzing, and visualizing 10X Genomics single-cell RNA-seq (scRNA-seq) data. Developed for a Bachelor thesis at ZHAW in collaboration with Nexco Analytics, the pipeline supports alignment with STARsolo and downstream analysis using Seurat.

## Repository Structure


## Repository Structure

```
├── Code/
│   ├── submit_pipeline.sh      # SLURM job submission wrapper
│   ├── submit_pipeline_2.sh    # Alternative SLURM logic (optional)
│   ├── run_pipeline.sh         # Core alignment and quantification script using STARsolo
│   ├── run_pipeline_2.sh       # Variant version for brain tumor dataset (optional)
│   ├── seurat_analysis.R       # R script for downstream analysis using Seurat
│   └── README.md               # Documentation explaining SLURMs and R script components
│
├── requirements/
│   └── seuratenv_enhanced.yml  # Conda environment for R
│
└── README.md
---

## Requirements

- STAR (v2.7.10a or newer)
- R (v4.2.0 or newer)
- Conda (Anaconda or Miniconda)

To install the R environment:

conda env create -f requirements/seuratenv_enhanced.yml  
conda activate seurat_env_3

## Usage

To run the full pipeline on an HPC cluster:

sbatch Code/submit_pipeline.sh

This script:
- Calls `run_pipeline.sh` for alignment and count matrix generation using STARsolo
- Automatically creates a timestamped output folder
- Triggers `seurat_analysis.R` for QC, clustering, cell type annotation, and visualization

All paths and parameters can be configured inside `run_pipeline.sh`.

## Input Datasets and Output Results

Due to GitHub file size limits, raw input data and pipeline results are hosted on OneDrive.

Download full datasets and outputs from:

https://your-onedrive-link-here

Directory structure on OneDrive:

Datasets/
├── pbmc_1k_v3_fastqs/          # PBMC 1k dataset (2 sequencing lanes)
├── Brain_Tumor_3p_fastqs/      # Brain tumor dataset (4 sequencing lanes)

Outputs_Pipeline/
├── Output_1k_PBMC/             # STARsolo + Seurat output for PBMC
├── Output_brain_tumor_3p/      # STARsolo + Seurat output for brain tumor

## Output Overview

Each pipeline run produces:
- STARsolo output (barcodes, genes, matrix.mtx)
- Filtered cell matrices and statistics
- Seurat-based downstream results:
  - QC metrics and violin plots
  - UMAP/t-SNE dimensionality reductions
  - Cell clustering and annotations
  - Differential expression and functional scoring
  - Interactive Dashboard
 
A more detailed documentation about the results can be found on the OneDrive folder

## Credits

This pipeline was developed by Luigi Palese as part of a Bachelor thesis in Applied Digital Life Sciences at ZHAW. Project supported by Nexco Analytics.
