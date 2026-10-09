#!/bin/bash

GTF="$1"
BAMDIR="$2"
OUTDIR="$3"
JOBS="${4:-5}"

echo "GTF: $GTF"
echo "BAMDIR: $BAMDIR"
echo "OUTDIR: $OUTDIR"
echo "JOBS: $JOBS"

if [ ! -f "$GTF" ]; then
    echo "ERROR: GTF file not found!"
    exit 1
fi

if [ ! -d "$BAMDIR" ]; then
    echo "ERROR: BAM directory not found!"
    exit 1
fi

mkdir -p "$OUTDIR"

export GTF OUTDIR

find "$BAMDIR" -name "*_Aligned.sortedByCoord.out.bam" | \
parallel -j "$JOBS" '
BAM={}
SAMPLE=$(basename "$BAM" _Aligned.sortedByCoord.out.bam)

echo "Counting $SAMPLE"

htseq-count \
    -f bam \
    -r pos \
    -s reverse \
    -t exon \
    -i gene_id \
    "$BAM" \
    "$GTF" \
    > "$OUTDIR/${SAMPLE}.counts.txt"
'

echo "HTSeq counting completed."
