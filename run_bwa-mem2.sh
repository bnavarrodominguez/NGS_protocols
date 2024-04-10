#!/bin/bash
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 reference fastq1 fastq2 output nthreads"
	exit 2
fi


# Define variables
REFERENCE=$1
FASTQ1=$2
FASTQ2=$3
OUT=$4
THREADS=$5

# Check if reference is indexed, and if not, index it
if [ ! -e "$REFERENCE.bwt.2bit.64" ]; then
    echo "Indexing reference genome..."
    bwa-mem2 index $REFERENCE
fi


# Map reads to reference using BWA-MEM and output SAM file
RG="@RG\tID:NL\tSM:$OUT\tPL:illumina\tLB:lib1\tPU:unit1"
#bwa mem -t $THREADS -M -R $RG $REFERENCE $FASTQ1 $FASTQ2 | samtools sort -o ${OUT}.sorted.bam -
bwa-mem2 mem -t $THREADS -M -R $RG $REFERENCE $FASTQ1 $FASTQ2 | samtools sort - -o ${OUT}.sorted.bam

# Index sorted BAM file${OUT}.sorted.bam

#bamtools index $(basename $FASTQ1 .fastq.gz).sorted.bam
samtools index ${OUT}.sorted.bam

## reduce bam size
#reduce_bam.py ${OUT}.sorted.bam

