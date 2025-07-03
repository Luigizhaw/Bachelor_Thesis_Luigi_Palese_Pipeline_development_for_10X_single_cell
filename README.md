# 10X Genomics scRNA-seq Pipeline: Preprocessing, Analysis & Visualization

This repository contains a modular pipeline for preprocessing, analyzing, and visualizing 10X Genomics single-cell RNA-seq (scRNA-seq) data. Developed for a Bachelor thesis at ZHAW in collaboration with Nexco Analytics, the pipeline supports alignment with STARsolo and downstream analysis using Seurat.

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
```

## Requirements

- STAR (v2.7.10a or newer)
- R (v4.2.0 or newer)
- Conda (Anaconda or Miniconda)

## Dependencies and Setup

Specific dependencies used in this project are all found in the .yml file in the requirements folder of this repository, the most important packages are listed below.

#### Core Analysis
- `Seurat` - Main single-cell analysis framework
- `optparse` - Command line argument parsing
- `dplyr` - Data manipulation
- `future` - Parallel processing

#### Visualization
- `ggplot2` - Static plotting
- `patchwork` - Plot composition
- `plotly` - Interactive visualizations
- `viridis` - Color scales
- `htmlwidgets` - Web widget creation

#### Cell Type Annotation
- `SingleR` - Reference-based cell type annotation

#### Data Export
- `jsonlite` - JSON export functionality
- `reshape2` - Data reshaping

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
```
Datasets/
├── pbmc_1k_v3_fastqs/          # PBMC 1k dataset (2 sequencing lanes)
├── Brain_Tumor_3p_fastqs/      # Brain tumor dataset (4 sequencing lanes)

Outputs_Pipeline/
├── Output_1k_PBMC/             # STARsolo + Seurat output for PBMC
├── Output_brain_tumor_3p/      # STARsolo + Seurat output for brain tumor
```
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

### Base Directory Structure used in the project
```
/cfs/earth/scratch/paleslui/BATH/
├── genomeDir/                      # STAR genome index
│   ├── Genome                      # Genome sequences
│   ├── SA                          # Suffix array
│   ├── SAindex                     # Index files
│   └── whitelist.txt              # Valid cell barcodes
├── pbmc_1k_v3_fastqs/             # PBMC dataset
│   ├── pbmc_1k_v3_S1_L001_R1_001.fastq.gz
│   └── pbmc_1k_v3_S1_L001_R2_001.fastq.gz
├── Brain_Tumor_3p_fastqs/         # Brain tumor dataset
│   ├── Brain_Tumor_3p_S2_L001_R1_001.fastq.gz
│   ├── Brain_Tumor_3p_S2_L001_R2_001.fastq.gz
│   ├── ... (L002, L003, L004)
├── R/                             # R analysis scripts
│   └── seurat_analysis_complete_12.R
└── SLURMs/                        # Pipeline scripts
    ├── submit_pipeline.sh
    ├── submit_pipeline_2.sh
    ├── run_pipeline.sh
    └── run_pipeline_2.sh
```

## Credits

This pipeline was developed by Luigi Palese as part of a Bachelor thesis in Applied Digital Life Sciences at ZHAW. Project supported by Nexco Analytics.
