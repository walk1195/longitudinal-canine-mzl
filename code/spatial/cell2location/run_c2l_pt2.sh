#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --tmp=16g
#SBATCH --gres gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=walk1195@umn.edu
#SBATCH --job-name=proj129_c2l_pt2
#SBATCH --account=modianoj


module load miniforge
source activate cell2location

python c2l_part2_script_canFam4.py