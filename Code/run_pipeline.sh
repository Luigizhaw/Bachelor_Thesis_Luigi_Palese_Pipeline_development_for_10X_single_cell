#!/bin/bash
# ❗ This is a BASH logic script
# It is called by submit_pipeline.sh from SLURM

# Create timestamped run ID
timestamp=$(date +"%Y-%m-%d_%H-%M")
run_id="run_${timestamp}"
base_dir="/cfs/earth/scratch/paleslui/BATH"
output_dir="${base_dir}/output/runs/${run_id}"
gene_out="${output_dir}/Gene"

mkdir -p "$gene_out"

# Run STARsolo
echo " Running STARsolo..."
STAR --runThreadN 8 \
     --genomeDir "${base_dir}/genomeDir/" \
    --readFilesIn "${base_dir}/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,${base_dir}/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz" \
                  "${base_dir}/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,${base_dir}/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz" \
     --readFilesCommand zcat \
     --soloType Droplet \
     --soloCBwhitelist "${base_dir}/genomeDir/whitelist.txt" \
     --soloFeatures Gene \
     --soloBarcodeReadLength 0 \
     --outFileNamePrefix "${gene_out}/"

# Compress STARsolo output files for Seurat compatibility
echo " Compressing output files..."
gzip "${gene_out}/Solo.out/Gene/filtered/"*.tsv "${gene_out}/Solo.out/Gene/filtered/"*.mtx

# Activate Conda and run Seurat R pipeline
echo " Running Seurat analysis..."
eval "$(conda shell.bash hook)"
conda activate seurat_env_3  # make sure this matches your installed environment

Rscript "${base_dir}/R/seurat_analysis.R" \
  --input_dir "${gene_out}/Solo.out/Gene/filtered/" \
  --output_dir "${output_dir}" \
  --prefix "${run_id}_"

# Done
echo " Pipeline completed for $run_id"
