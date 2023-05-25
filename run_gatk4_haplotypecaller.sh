#!/bin/bash
if [ "$#" -ne 4 ]; then
	    echo "Usage: $0 sample.bam reference.fasta region[chr or chr:start-end] ploidy"
	            exit 2
fi



bam=$1
ref=$2
chr=$3
p=$4

##### check reference indexes

if [ -f ${ref}.fai ]; then
	echo "${ref}.fai exists, moving on..."
else	
	echo "Generating ${ref}.fai index ..."
	        samtools faidx $ref
fi

# Check if reference is indexed, and if not, index it
if [ -f "$(basename $ref .fasta).dict" ]; then
	echo "$(basename $ref .fasta).dict exists, moving on..."
else
	echo "Creating $ref dictionary..."
	    picard CreateSequenceDictionary R=$ref O=$(basename $ref .fasta).dict
fi



#### Dedup
#echo "


if [ -f "$(basename $bam .bam).dedup.bam" ]; then
	echo "$(basename $bam .bam).dedup.bam exists, moving on ..."
else
	echo "Marking duplicates..."
picard MarkDuplicates \
	I=${bam} \
	O=$(basename $bam .bam).dedup.bam \
	METRICS_FILE= $(basename $bam .bam).dup_metrics.txt\
	REMOVE_DUPLICATES=true \
	VALIDATION_STRINGENCY=LENIENT AS=true 
samtools index $(basename $bam .bam).dedup.bam
fi

## haplotype caller

if [ -f "$(basename $bam .bam).${chr}.vcf.gz" ]; then
	echo "$(basename $bam .bam).${chr}.vcf.gz exists, moving on ..."
else 
	echo "Running GATK Haplotype Caller ..."
	gatk --java-options "-Xmx10g" HaplotypeCaller \
	-R $ref \
	-I $(basename $bam .bam).dedup.bam \
	-O $(basename $bam .bam).vcf.gz \
	-ploidy $p \
	-L $chr
	-ERC GVCF 

fi


