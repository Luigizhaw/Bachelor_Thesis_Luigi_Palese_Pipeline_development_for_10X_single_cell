# SLURMs and R Script Overview

This README file contains deteiled information about the structure and usage of the used SLURMs and R Script.

## Pipeline Workflow

The submission scripts integrate with the pipeline as follows:

```
SLURM Job Scheduler
        ↓
submit_pipeline.sh → run_pipeline.sh → seurat_analysis_complete_12.R
        ↓                    ↓                       ↓
Resource Management    Sequence Alignment    Scientific Analysis
```

# SLURM Submission

## Overview

These scripts handle job submission to the SLURM workload manager on the HPC earth cluster. They manage resource allocation and job scheduling for the 10X single-cell analysis pipeline.

## Files

### `submit_pipeline.sh`
SLURM submission script for the PBMC 1k v3 dataset.

### `submit_pipeline_2.sh`
SLURM submission script for the Brain Tumor 3p dataset (Included to show difference in path definition of the input fastqs).

## Script Structure

Both scripts follow identical SLURM configuration:

```
#!/bin/bash
#SBATCH --job-name=10x_pipeline
#SBATCH --output=/cfs/earth/scratch/paleslui/BATH/output/runs/logs/pipeline_%j.log
#SBATCH --error=/cfs/earth/scratch/paleslui/BATH/output/runs/logs/pipeline_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=earth-3

# Call the actual pipeline logic
bash /cfs/earth/scratch/paleslui/BATH/SLURMs/run_pipeline.sh      # submit_pipeline.sh
bash /cfs/earth/scratch/paleslui/BATH/SLURMs/run_pipeline_2.sh    # submit_pipeline_2.sh
```

## SLURM Configuration

### Resource Allocation

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `--cpus-per-task=16` | 16 CPU cores | STARsolo benefits from multi-threading; Seurat uses parallel processing |
| `--mem=128G` | 128GB RAM | Required for loading large count matrices and dimensionality reduction |
| `--time=1-00:00:00` | 24 hours | Conservative limit allowing for thorough analysis of large datasets |
| `--partition=earth-3` | earth-3 | Dedicated partition |

### Job Management

| Parameter | Purpose | Details |
|-----------|---------|---------|
| `--job-name=10x_pipeline` | Identification | Consistent naming for easy job tracking |
| `--output=...pipeline_%j.log` | Standard output | Captures all console output with unique job ID |
| `--error=...pipeline_%j.err` | Error output | Separate error stream for debugging |

## Usage

### Basic Submission

```bash
# Submit PBMC analysis
cd /cfs/earth/scratch/paleslui/BATH/SLURMs/
sbatch submit_pipeline.sh

# Submit Brain Tumor analysis  
sbatch submit_pipeline_2.sh
```

### Job Monitoring

```bash
# Check job status
squeue -u paleslui

# Monitor specific job
squeue -j <job_id>
```

## Error Handling

### Common SLURM Errors

| Error Message | Cause | Solution |
|---------------|-------|----------|
| `Batch job submission failed: Requested node configuration is not available` | Insufficient resources | Wait for resources or reduce requirements |
| `slurmstepd: error: Exceeded memory limit` | Job used more than 128GB | Increase memory allocation |
| `DUE TO TIME LIMIT` | Job exceeded time limit | Increase time limit or optimize analysis |

### Email Notifications (Optional/ nice to have)

- This slight modification lets the user receive an email with job details as soon as the job is finished.

```bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@domain.com
```

# Run Pipeline Scripts (`run_pipeline.sh` & `run_pipeline_2.sh`)

## Overview

These scripts orchestrate the core data processing workflow, handling STARsolo alignment and Seurat analysis execution. They're called by the SLURM submission scripts explained above and manage the complete pipeline from raw FASTQ files to final analysis outputs.

## Files

### `run_pipeline.sh`
Processes the PBMC 1k v3 dataset.

### `run_pipeline_2.sh`
Processes the Brain Tumor 3p dataset.

## Script Architecture

Both scripts follow identical logic:

```bash
#!/bin/bash
# ❗ This is a BASH logic script
# It is called by submit_pipeline.sh from SLURM
```

## Workflow Steps

### 1. Run Identification and Setup

```bash
timestamp=$(date +"%Y-%m-%d_%H-%M")
run_id="run_${timestamp}"
base_dir="/cfs/earth/scratch/paleslui/BATH"
output_dir="${base_dir}/output/runs/${run_id}"
gene_out="${output_dir}/Gene"

mkdir -p "$gene_out" "$output_dir/config" "$output_dir/logs"
```

**Purpose**:
- Creates unique timestamped run identifier
- Sets up organized directory structure
- Ensures no conflicts between simultaneous runs as long as they have a difference in timestamps (1 minute)

**Directory Structure Created**:
```
output/runs/run_2024-01-15_14-30/
├── Gene/                    # STARsolo output location
├── config/                  # Run configuration files
└── logs/                    # Processing logs
```

### 2. STARsolo Alignment

#### PBMC Dataset
```bash
STAR --runThreadN 8 \
     --genomeDir "${base_dir}/genomeDir/" \
     --readFilesIn                         "${base_dir}/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,${base_dir}/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz" \
                  "${base_dir}/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,${base_dir}/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz" \
    --readFilesCommand zcat \
     --soloType Droplet \
     --soloCBwhitelist "${base_dir}/genomeDir/whitelist.txt" \
     --soloFeatures Gene \
     --soloBarcodeReadLength 0 \
     --outFileNamePrefix "${gene_out}/"
```

#### Brain Tumor Dataset
```bash
STAR --runThreadN 8 \
     --genomeDir "${base_dir}/genomeDir/" \
     --readFilesIn \
        "${fastq_dir}/Brain_Tumor_3p_S2_L001_R2_001.fastq.gz","${fastq_dir}/Brain_Tumor_3p_S2_L002_R2_001.fastq.gz","${fastq_dir}/Brain_Tumor_3p_S2_L003_R2_001.fastq.gz","${fastq_dir}/Brain_Tumor_3p_S2_L004_R2_001.fastq.gz" \
        "${fastq_dir}/Brain_Tumor_3p_S2_L001_R1_001.fastq.gz","${fastq_dir}/Brain_Tumor_3p_S2_L002_R1_001.fastq.gz","${fastq_dir}/Brain_Tumor_3p_S2_L003_R1_001.fastq.gz","${fastq_dir}/Brain_Tumor_3p_S2_L004_R1_001.fastq.gz" \
     --readFilesCommand zcat \
     --soloType Droplet \
     --soloCBwhitelist "${base_dir}/genomeDir/whitelist.txt" \
     --soloFeatures Gene \
     --soloBarcodeReadLength 0 \
     --outFileNamePrefix "${gene_out}/"
```

#### Key Parameters Explained

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `--runThreadN 8` | 8 threads | Balances speed with resource usage |
| `--genomeDir` | Pre-built index | Reference genome location |
| `--readFilesCommand zcat` | Decompress on-the-fly | Handles gzipped FASTQ files |
| `--soloType Droplet` | Droplet protocol | 10X Genomics chemistry |
| `--soloCBwhitelist` | Valid barcodes | Cell barcode filtering |
| `--soloFeatures Gene` | Gene-level | Count reads per gene |
| `--soloBarcodeReadLength 0` | Auto-detect | Automatic barcode length detection |

### 3. Output Compression

```bash
echo " Compressing output files..."
gzip "${gene_out}/Solo.out/Gene/filtered/"*.tsv "${gene_out}/Solo.out/Gene/filtered/"*.mtx
```

**Purpose**:
- **Seurat compatibility**: `Read10X()` expects compressed files

**Files compressed**:
- `barcodes.tsv` → `barcodes.tsv.gz`
- `features.tsv` → `features.tsv.gz`  
- `matrix.mtx` → `matrix.mtx.gz`

### 4. Conda Environment Activation

```bash
echo " Running Seurat analysis..."
eval "$(conda shell.bash hook)"
conda activate seurat_env_3
```

### 5. R Script Execution

```bash
Rscript "${base_dir}/R/seurat_analysis_complete_12.R" \
  --input_dir "${gene_out}/Solo.out/Gene/filtered/" \
  --output_dir "${output_dir}" \
  --prefix "${run_id}_"
```

## Best Practices

### Resource Management
- **Start conservative**: Use standard resources first, then adjust based on actual usage
- **Monitor efficiency**: Check if you're over/under-allocating resources
- **Consider dataset size**: Larger datasets require more memory and time

### Job Organization
- **Meaningful names**: Use descriptive job names for easy identification
- **Structured logs**: Keep logs organized in dedicated directories
- **Version tracking**: Document changes to submission scripts

### Cluster Etiquette
- **Request appropriately**: Only request resources you actually need
- **Clean up regularly**: Remove old log files and temporary data
- **Queue awareness**: Check partition status before submitting many jobs

## Integration Notes

### File Paths
All paths are configured for the BATH project structure and can be adjusted in the run_pipeline.sh file:
- **Base directory**: `/cfs/earth/scratch/paleslui/BATH`
- **Scripts location**: `/cfs/earth/scratch/paleslui/BATH/SLURMs/`
- **Log directory**: `/cfs/earth/scratch/paleslui/BATH/output/runs/logs/`


### Environment Requirements
The target bash scripts handle:
- Conda environment activation
- Tool availability verification
- Output directory creation

# R Script Components

### 1. Initialization and Setup

#### Parallel Processing Configuration
```r
# Configures multicore processing based on ncores parameter
if (ncores > 1) {
  plan("multicore", workers = ncores)
} else {
  plan("sequential")
}
```

#### Random Seed Management
- If `seed = 0`: Generates random seed for reproducibility
- If `seed > 0`: Uses specified seed
- Ensures reproducible clustering and dimensionality reduction

#### Output Directory Structure
Creates organized directory structure:
```
output_dir/
├── prefix_timestamp/
│   ├── 01_QC/           # Quality control plots
│   ├── 02_Clustering/   # Clustering results
│   ├── 03_Markers/      # Marker gene analysis
│   ├── 04_Functional/   # Functional analysis
│   ├── 05_Interactive/  # Interactive dashboards
│   └── 06_JSON_Export/  # JSON data exports
```

### 2. Enhanced Gene Signatures

#### Functional Pathway Signatures
The pipeline includes curated gene sets for human cells:

```r
get_human_gene_sets <- function() {
  return(list(
    # Cellular stress and response
    inflammation = c("TNF", "IL1B", "IL6", "NFKB1", "CXCL8", "CCL2"),
    stress_response = c("HSP90AA1", "HSPA1A", "HSPB1", "DNAJB1", "HSP90AB1"),
    interferon_response = c("IFIT1", "IFIT2", "IFIT3", "MX1", "ISG15", "OAS1"),
    
    # Cell fate and function
    apoptosis = c("BAX", "BAK1", "CASP3", "CASP7", "PARP1"),
    
    # Metabolic states
    glycolysis = c("HK1", "HK2", "PFKP", "ALDOA", "GAPDH", "PKM"),
    oxidative_phosphorylation = c("COX4I1", "COX5A", "COX6A1", "ATP5F1A", "NDUFB1"),
    
    # Quality control
    mitochondrial_activity = c("MT-CO1", "MT-CO2", "MT-ND1", "MT-ND4", "MT-CYB"),
    ribosomal_activity = c("RPL3", "RPS3", "RPS18", "RPL4", "RPS4X", "RPL5")
  ))
}
```

#### Signature Score Calculation
- Uses `AddModuleScore()` for robust multi-gene signatures
- Accounts for background gene expression
- Enables pathway-level analysis beyond individual genes

### 3. Biology-Driven Resolution Selection

#### Intelligent Resolution Selection
```r
select_best_resolution_by_biology <- function(seurat_obj, resolutions) {
  # Tests multiple resolutions
  # Evaluates marker gene quality and biological interpretability
  # Selects optimal resolution based on:
  # - Number of clusters (3-15 range)
  # - Marker gene strength (log2FC and significance)
  # - Biological coherence score
}
```

**Selection Criteria:**
- Cluster count between 3-15 (biologically reasonable)
- Strong marker genes (log2FC > 0.25, min.pct > 0.15)
- High average log2FC across markers
- Biological interpretability score

### 4. Enhanced Cell Type Assignment

#### Dual-Approach Cell Type Annotation

**Step 1: Individual Cell Annotation (SingleR)**
```r
# Cell-by-cell annotation using reference databases
ref <- celldex::HumanPrimaryCellAtlasData()
singler_results <- SingleR(test = seurat_sce, ref = ref, labels = ref$label.fine)
```

**Step 2: Cluster Consensus Assignment**
```r
assign_celltypes_improved <- function(seurat_obj, markers) {
  # For each cluster:
  # 1. Collect individual cell type annotations
  # 2. Find consensus (most common types)
  # 3. Integrate with marker gene analysis
  # 4. Assign confidence scores
  # 5. Provide alternative type suggestions
}
```

**Output Annotations:**
- `individual_celltype`: Cell-by-cell SingleR annotations
- `data_driven_celltype`: Cluster consensus annotations
- `celltype_confidence`: Assignment confidence scores

### 5. Enhanced Cell Cycle Analysis

#### G0 Detection
```r
detect_g0_enhanced <- function(seurat_obj) {
  # Uses Seurat's built-in S and G2M gene sets
  # Calculates S.Score and G2M.Score
  # Identifies G0 cells: low S AND low G2M (bottom 25%)
  # Creates enhanced phase annotation: G0, G1, S, G2M
}
```

**Cell Cycle Phases:**
- **G0**: Quiescent cells (low proliferation markers)
- **G1**: Gap 1 phase
- **S**: DNA synthesis phase  
- **G2M**: Gap 2 and mitosis phases

## Data Processing Modules

### Module 1: Data Loading and Quality Control

#### Data Import
```r
# Reads 10X Genomics format data
data <- Read10X(input_dir, gene.column = 2, strip.suffix = TRUE)
seurat_obj <- CreateSeuratObject(counts = data, project = prefix, 
                                min.cells = 3, min.features = 200)
```

#### Adaptive QC Thresholds
```r
# Calculates adaptive thresholds based on data distribution
nf_low <- quantile(seurat_obj$nFeature_RNA, 0.05)    # Minimum genes
nf_high <- quantile(seurat_obj$nFeature_RNA, 0.95)   # Maximum genes  
mt_high <- quantile(seurat_obj$percent.mt, 0.95)     # Maximum mitochondrial %
```

**QC Metrics:**
- `nFeature_RNA`: Number of genes detected per cell
- `nCount_RNA`: Total UMI counts per cell
- `percent.mt`: Mitochondrial gene percentage

**Filtering Strategy:**
- Removes cells with too few genes (bottom 5%)
- Removes cells with too many genes (top 5% - potential doublets)
- Removes cells with high mitochondrial content (top 5% - dying cells)

### Module 2: Normalization and Dimensionality Reduction

#### Normalization and Dimensionality reduction
```r
# Standard Seurat workflow
seurat_obj <- NormalizeData(seurat_obj)                    # Log normalization
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)  # Feature selection
seurat_obj <- ScaleData(seurat_obj)                        # Z-score scaling
seurat_obj <- RunPCA(seurat_obj, npcs = 50)               # PCA
seurat_obj <- RunUMAP(seurat_obj, dims = 1:optimal_pcs)   # UMAP
seurat_obj <- RunTSNE(seurat_obj, dims = 1:optimal_pcs)   # t-SNE
```

### Module 3: Clustering and Cell Type Definition

#### Clustering Workflow
```r
# Graph-based clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:optimal_pcs)
for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
}
optimal_res <- select_best_resolution_by_biology(seurat_obj, resolutions)
```

#### Marker Gene Analysis
```r
# Find differentially expressed genes
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, 
                         min.pct = 0.25, logfc.threshold = 0.5)
```

**Marker Gene Criteria:**
- `only.pos = TRUE`: Only upregulated genes
- `min.pct = 0.25`: Expressed in ≥25% of cluster cells
- `logfc.threshold = 0.5`: Minimum log2 fold change

### Module 4: Functional Analysis

#### Signature Score Analysis
```r
# Calculate pathway signatures for each cell
for (sig_name in names(gene_sets)) {
  seurat_obj <- AddModuleScore(seurat_obj, 
                              features = list(available_genes), 
                              name = paste0(sig_name, "_Score"))
}
```

### Module 5: Interactive Dashboard Creation

#### Dashboard Structure
```html
<!DOCTYPE html>
<html>
  <head>
    <!-- Comprehensive styling and responsive design -->
  </head>
  <body>
    <div class="dashboard-container">
      <div class="header">
        <!-- Analysis summary statistics -->
      </div>
      <div class="main-content">
        <div class="sidebar">
          <!-- Organized plot navigation -->
        </div>
        <div class="plot-container">
          <!-- Interactive plot iframes -->
        </div>
      </div>
    </div>
  </body>
</html>
```

### Module 6: JSON Export System

#### Comprehensive Data Export
```r
comprehensive_export <- list(
  metadata = list(
    dataset_info = list(...),     # Dataset statistics
    processing_info = list(...),  # Analysis parameters
    qc_thresholds = list(...)     # Quality control cutoffs
  ),
  outputs = list(
    plots = plots_catalog,        # File paths and descriptions
    data_files = data_files_catalog  # Data file locations
  ),
  data = list(
    cells = cell_data,           # Per-cell information
    clusters = cluster_stats,     # Cluster statistics
    markers = markers_data,       # Marker genes
    plots = plots_data           # Plot-ready data
  ),
  summary = list(...)            # Summary statistics
)
```

#### Export Components
- **Cell-level data**: Coordinates, annotations, QC metrics
- **Cluster statistics**: Centroids, composition, marker genes
- **Plot data**: Ready for web visualization
- **Metadata**: Analysis parameters and quality metrics
- **File catalog**: Organized index of all outputs


**Note**: This pipeline is designed for human single-cell RNA-seq data.
