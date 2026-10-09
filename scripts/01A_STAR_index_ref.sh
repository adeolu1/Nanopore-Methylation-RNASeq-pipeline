#!/bin/bash

## Reference .fa
FASTA=$1
## Annotation .gtf
GTF=$2
## output directory for the generated index
INDEX=$3
## Number of CPU to be used 
THREADS=${4:-10}

echo "FASTA: $FASTA"
echo "GTF: $GTF"
echo "INDEX: $INDEX"

if [ ! -f "$FASTA" ]; then
    echo "ERROR: FASTA file not found!"
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "ERROR: GTF file not found!"
    exit 1
fi

##### make output directory for reference index  file
mkdir -p "$INDEX"

STAR \
    --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX" \
    --genomeFastaFiles "$FASTA" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang 99 \
    --genomeSAindexNbases 13

echo "All Job finished on $(date)"
