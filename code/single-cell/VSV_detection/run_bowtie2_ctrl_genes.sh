#!/bin/bash -l

#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --output=bowtie2-ctrl-genes-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-6 # 6 samples 


module purge
module load bowtie2
module load samtools

# Build bowtie2 index for positive ctrl gene -- only do this once
GENE="CD19"
GENE2="GADPH"

bowtie2-build /pos_control/${GENE}_WGS.fa /pos_control/${GENE}_WGS # Update with full path
bowtie2-build /pos_control/${GENE2}_WGS.fa /pos_control/${GENE2}_WGS # Update with full path

# Create an array of sample names
SAMPLES=("ORv01PreTx" "ORv01D3" "ORv01D7FNA" "ORv01D10FNA" "ORv01D10resect" "ORv01D28FNA", "Orv01D48FNA")

# Define base directory containing samples
SAMPLE_DIR="" # Update with full path

# Get curr sample name
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}


echo "Running bowtie2 against control gene sequences for ${SAMPLE_NAME}...."


# Set output dir
OUT_DIR="/vsv_seq/ctrl_genes2/${SAMPLE_NAME}/$GENE" 
INDEX_DIR="/pos_control/${GENE}_WGS"

OUT_DIR2="/vsv_seq/ctrl_genes2/${SAMPLE_NAME}/$GENE2" 
INDEX_DIR2="/pos_control/${GENE2}_WGS"

# ensure dir exists
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR2


# Set path to input files
READ1="${SAMPLE_DIR}/${SAMPLE_NAME}_GEX_FL_S1_R1_001.fastq.gz"
READ2="${SAMPLE_DIR}/${SAMPLE_NAME}_GEX_FL_S1_R2_001.fastq.gz"


# Run bowtie2
echo "Running bowtie2 for $GENE"
bowtie2 -p 8 -x $INDEX_DIR -1 $READ1 -2 $READ2 -S ${OUT_DIR}/${SAMPLE_NAME}.sam
echo "Running bowtie2 for $GENE2"
bowtie2 -p 8 -x $INDEX_DIR2 -1 $READ1 -2 $READ2 -S ${OUT_DIR2}/${SAMPLE_NAME}.sam


