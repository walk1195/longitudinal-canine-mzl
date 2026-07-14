#!/bin/bash -l

#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=walk1195@umn.edu

module load R/4.4.0-openblas-rocky8


Rscript run_featurecounts_ctrl.R >> featurecounts_ctrl_summary.txt