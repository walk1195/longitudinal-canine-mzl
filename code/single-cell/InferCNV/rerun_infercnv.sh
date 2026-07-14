#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --partition=msismall
#SBATCH --mem=120000mb # 120gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=walk1195@umn.edu
#SBATCH --job-name=rerun_infercnv
#SBATCH -o infercnv_%A_%a.out
#SBATCH -e infercnv_%A_%a.err
#SBATCH --account=modianoj
#SBATCH --array=1


module load R/4.4.0-openblas-rocky8
module load JAGS/4.3.0-gcc7.2.0

SAMPLE_LIST=(
  "PreTx"
)

# ---- map array index to sample ----
SAMPLE_ID=${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID-1]} # indexing sample list


Rscript rerun_infercnv_updated_parameters.R -s ${SAMPLE_ID}