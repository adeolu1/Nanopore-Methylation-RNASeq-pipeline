#!/bin/bash
##
# STAR genome index
#GENOME=/path/to/ref
GENOME=$1
# Directory containing FASTQ files
#FASTQDIR=/media/adeolu/DATA/EpiTrans/20dpf
FASTQDIR=$2

# Output directory
#OUTDIR=/media/adeolu/DATA/EpiTrans/20_dpf_analysis
OUTDIR=$3
THREADS="${4:-10}"   # Default to 10 threads if not provided

mkdir -p "$OUTDIR"

run_star () {
    R1="$1"
    R2="${R1%_1.fq.gz}_2.fq.gz"
    SAMPLE=$(basename "$R1" _1.fq.gz)
    mkdir -p "$OUTDIR/$SAMPLE"
    echo "Starting $SAMPLE"
    STAR \
        --runMode alignReads \
        --genomeDir "$GENOME" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$OUTDIR/$SAMPLE/${SAMPLE}_" \
        --runThreadN "$THREADS" \
        --quantMode GeneCounts
}

export -f run_star
export GENOME OUTDIR THREADS

find "$FASTQDIR" -name "*_1.fq.gz" | parallel -j 4 run_star {}

echo "All alignmnet Job finished on $(date)"
