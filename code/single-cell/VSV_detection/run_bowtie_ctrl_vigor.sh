#!/bin/bash -l

#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --job-name=bowtie2-vigor-ctrl-array
#SBATCH --output=bowtie2-vigor-ctrl-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-4 



module purge
module load bowtie2
module load samtools

### Index for VSV sequence already exists

# Create an array of sample names
SAMPLES=("MN25_bone_post_1_S53" "MN06_bone_post_S16" "MN12_bone_post_S7" "MN27_bone_cell_S64")

# Define base directory containing samples
SAMPLE_DIR=""

# Get curr sample name
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}


echo "Running bowtie2 against VSV sequence for ${SAMPLE_NAME}...."

# Set output dir
OUT_DIR="" # Update with full path
INDEX_DIR="/vsv_seq/bowtie2/vsv"

# ensure dir exists
mkdir $OUT_DIR


# Set path to input files
READ1="${SAMPLE_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz"
READ2="${SAMPLE_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"


# Run bowtie2
bowtie2 -p 8 -x $INDEX_DIR -1 $READ1 -2 $READ2 -S ${OUT_DIR}/${SAMPLE_NAME}.sam
# change n to indicate different CPU thread count

# Convert sam file to unsorted bam
#samtools view -Sb ${OUT_DIR}/${sample}.sam > ${OUT_DIR}/unsorted_${sample}.bam

# Sort bam files
#samtools sort ${OUT_DIR}/unsorted_${sample}.bam -o ${OUT_DIR}/${sample}.bam






