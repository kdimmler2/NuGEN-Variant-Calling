#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dimml002@umn.edu
#SBATCH --job-name nugen.slurm
#SBATCH -o ./jobs/%j.nugen_pipeline.out
#SBATCH -e ./jobs/%j.nugen_pipeline.err
#SBATCH -p agsmall,aglarge

set -e

cd $SLURM_SUBMIT_DIR

conda activate nugen

snakemake -s nugen_pipeline.smk \
    --profile slurm.nugen \
    --keep-going \
    --rerun-incomplete \
    --rerun-triggers mtime \
