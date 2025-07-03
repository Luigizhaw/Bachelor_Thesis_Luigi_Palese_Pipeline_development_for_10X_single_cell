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

## Directory Structure and File Paths

### Base Directory Structure
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
All paths are configured for the BATH project structure:
- **Base directory**: `/cfs/earth/scratch/paleslui/BATH`
- **Scripts location**: `/cfs/earth/scratch/paleslui/BATH/SLURMs/`
- **Log directory**: `/cfs/earth/scratch/paleslui/BATH/output/runs/logs/`

These need to be adjusted in the run_

### Environment Requirements
The target bash scripts handle:
- Conda environment activation
- Tool availability verification
- Output directory creation
