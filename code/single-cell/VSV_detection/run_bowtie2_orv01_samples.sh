#!/bin/bash -l

#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --job-name=bowtie2-vsv-array
#SBATCH --output=bowtie2-vsv-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-6 # 6 samples 


module load bowtie2
module load samtools

# Build bowtie2 index -- only do this once
#bowtie2-build /vsv_seq/bowtie2/vsv.fa /vsv_seq/bowtie2/vsv


# Create an array of sample names
SAMPLES=("ORv01PreTx" "ORv01D3" "ORv01D7FNA" "ORv01D10FNA" "ORv01D10resect" "ORv01D28FNA", "ORv01D49FNA")

# Define base directory containing samples
SAMPLE_DIR=""

# Get curr sample name
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}


echo "Running bowtie2 against VSV sequence for ${SAMPLE_NAME}...."

# Set output dir
OUT_DIR="/vsv_seq/bowtie2-run2" 
INDEX_DIR="/vsv_seq/bowtie2/vsv"

# ensure dir exists
mkdir $OUT_DIR


# Set path to input files
READ1="${SAMPLE_DIR}/${SAMPLE_NAME}_GEX_FL_S1_R1_001.fastq.gz"
READ2="${SAMPLE_DIR}/${SAMPLE_NAME}_GEX_FL_S1_R2_001.fastq.gz"


# Run bowtie2
bowtie2 -p 8 -x $INDEX_DIR -1 $READ1 -2 $READ2 -S ${OUT_DIR}/${SAMPLE_NAME}.sam
# change n to indicate different CPU thread count

# Convert sam file to unsorted bam
#samtools view -Sb ${OUT_DIR}/${sample}.sam > ${OUT_DIR}/unsorted_${sample}.bam

# Sort bam files
#samtools sort ${OUT_DIR}/unsorted_${sample}.bam -o ${OUT_DIR}/${sample}.bam


