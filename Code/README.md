# SLURM Submission

## Overview

These scripts handle job submission to the SLURM workload manager on the HPC earth cluster. They manage resource allocation and job scheduling for the 10X single-cell analysis pipeline.

## Files

### `submit_pipeline.sh`
SLURM submission script for the PBMC 1k v3 dataset (single-lane sequencing).

### `submit_pipeline_2.sh`
SLURM submission script for the Brain Tumor 3p dataset (multi-lane sequencing).

## Script Structure

Both scripts follow identical SLURM configuration with different target processing scripts:

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
| `--partition=earth-3` | earth-3 | Dedicated partition optimized for computational biology workloads |

### Job Management

| Parameter | Purpose | Details |
|-----------|---------|---------|
| `--job-name=10x_pipeline` | Identification | Consistent naming for easy job tracking |
| `--output=...pipeline_%j.log` | Standard output | Captures all console output with unique job ID |
| `--error=...pipeline_%j.err` | Error output | Separate error stream for debugging |

### Log File Location
```
/cfs/earth/scratch/paleslui/BATH/output/runs/logs/
├── pipeline_12345.log    # Standard output
└── pipeline_12345.err    # Error output
```
Where `12345` is the SLURM job ID (`%j`).

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

# View job details
scontrol show job <job_id>
```

### Log Monitoring

```bash
# Real-time log viewing
tail -f /cfs/earth/scratch/paleslui/BATH/output/runs/logs/pipeline_<job_id>.log

# Check for errors
less /cfs/earth/scratch/paleslui/BATH/output/runs/logs/pipeline_<job_id>.err
```

## Workflow Integration

The submission scripts integrate with the pipeline as follows:

```
SLURM Job Scheduler
        ↓
submit_pipeline.sh → run_pipeline.sh → seurat_analysis_complete_12.R
        ↓                    ↓                       ↓
Resource Management    Environment Setup    Scientific Analysis
```

### Separation of Concerns

1. **SLURM scripts**: Handle job scheduling and resource allocation
2. **Bash scripts**: Manage environment setup and tool orchestration
3. **R scripts**: Perform scientific analysis and visualization

## Dataset-Specific Differences

### PBMC Dataset (`submit_pipeline.sh`)
- **Target**: Single-lane sequencing data
- **Processing script**: `run_pipeline.sh`
- **Input**: `pbmc_1k_v3_fastqs/` directory
- **Characteristics**: Standard 10X protocol, ~1,000 cells

### Brain Tumor Dataset (`submit_pipeline_2.sh`)
- **Target**: Multi-lane sequencing data  
- **Processing script**: `run_pipeline_2.sh`
- **Input**: `Brain_Tumor_3p_fastqs/` directory
- **Characteristics**: 4 sequencing lanes, tissue sample

## Error Handling

### Common SLURM Errors

| Error Message | Cause | Solution |
|---------------|-------|----------|
| `Batch job submission failed: Requested node configuration is not available` | Insufficient resources | Wait for resources or reduce requirements |
| `slurmstepd: error: Exceeded memory limit` | Job used more than 128GB | Increase memory allocation |
| `DUE TO TIME LIMIT` | Job exceeded 24 hours | Increase time limit or optimize analysis |
| `Permission denied` | Cannot write to log directory | Check directory permissions |

### Debugging Steps

1. **Check resource availability**:
   ```bash
   sinfo -p earth-3
   ```

2. **Review job efficiency**:
   ```bash
   sacct -j <job_id> --format=JobID,MaxRSS,MaxVMSize,CPUTime,State
   ```

3. **Examine exit codes**:
   ```bash
   sacct -j <job_id> --format=JobID,ExitCode,State
   ```

## Customization

### Adjusting Resources

For larger datasets:
```bash
#SBATCH --cpus-per-task=32      # More CPU cores
#SBATCH --mem=256G              # More memory
#SBATCH --time=2-00:00:00       # Longer time limit
```

For smaller datasets:
```bash
#SBATCH --cpus-per-task=8       # Fewer CPU cores
#SBATCH --mem=64G               # Less memory
#SBATCH --time=12:00:00         # Shorter time limit
```

### Different Partitions

```bash
#SBATCH --partition=earth-2     # Alternative partition
#SBATCH --partition=gpu         # If GPU acceleration needed
```

### Email Notifications

```bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@domain.com
```

### Job Arrays

For processing multiple datasets:
```bash
#SBATCH --array=1-10            # Process 10 datasets
#SBATCH --array=1-10%3          # Limit to 3 concurrent jobs
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

## Advanced Features

### Job Dependencies

Chain multiple jobs:
```bash
# Submit dependent jobs
job1=$(sbatch --parsable submit_pipeline.sh)
sbatch --dependency=afterok:$job1 submit_pipeline_2.sh
```

### Exclusive Node Access

For very large datasets:
```bash
#SBATCH --exclusive
```

### Quality of Service

Set priority levels:
```bash
#SBATCH --qos=high
```

## Integration Notes

### File Paths
All paths are configured for the BATH project structure:
- **Base directory**: `/cfs/earth/scratch/paleslui/BATH`
- **Scripts location**: `/cfs/earth/scratch/paleslui/BATH/SLURMs/`
- **Log directory**: `/cfs/earth/scratch/paleslui/BATH/output/runs/logs/`

### Environment Requirements
The target bash scripts handle:
- Conda environment activation
- Tool availability verification
- Output directory creation

### Success Indicators
Successful submission results in:
- Job ID returned by `sbatch`
- Log files created in designated directory
- Job appears in `squeue` output
- Processing begins automatically when resources available

---

*These submission scripts provide the entry point for the 10X analysis pipeline on the HPC earth cluster, ensuring proper resource allocation and job management for reliable processing of single-cell datasets.*