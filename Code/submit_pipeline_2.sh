#!/bin/bash
#SBATCH --job-name=10x_pipeline
#SBATCH --output=/cfs/earth/scratch/paleslui/BATH/output/runs/logs/pipeline_%j.log
#SBATCH --error=/cfs/earth/scratch/paleslui/BATH/output/runs/logs/pipeline_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=earth-3

# Call the actual pipeline logic
bash /cfs/earth/scratch/paleslui/BATH/SLURMs/run_pipeline_2.sh
