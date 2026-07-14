#!/bin/bash -l
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --tmp=16g
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=walk1195@umn.edu

module load miniforge
source activate cell2location
module load python/3.10.9_anaconda2023.03_libmamba

python c2l_script.py

