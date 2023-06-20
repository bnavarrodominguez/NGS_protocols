#!/bin/bash
if [ "$#" -ne 7 ]; then
		    echo "Usage: $0 map_file.txt reference.fasta out_prefix interval[chr or chr:start-stop] ploidy nthreads filterHet[yes/no]"
		    	            exit 2
fi

#### The sample map (first argument) is a tab-delimited text file with sample_name\tpath_to_sample_vcf per line.
files=$1
###3 the reference fasta file
ref=$2
#### a name for the database
out=$3
####3 region to analyze in format chr (for the full chromosome) or chr:start-stop (for region within a chromosome)
region=$4
### number of threads
nthr=$6
#### ploidy
p=$5

#### filter het
fhet=$7


## Make database

if [ -d ${out}_database ]; then
	echo "${out}_database already exists; moving on ..." 
else
echo "Importing interval $region to ${out}_database"
	gatk --java-options "-Xmx40g -Xms40g" GenomicsDBImport \
	            --sample-name-map $files \
	                --genomicsdb-workspace-path ${out}_database \
	                --reader-threads $nthr \
			-L $region
fi

#### Genotype population

if [ -f ${out}.all_variants.vcf.gz ]; then
	echo "${out}.vcf.gz already exists, moving on..."
else
	echo "Genotyping cohort... ${out} "
	gatk --java-options "-Xmx4g" GenotypeGVCFs \
	   -R $ref \
	   -V gendb://${out}_database \
	   -O ${out}.all_variants.vcf.gz \
	   --ploidy $p \
	   -L $region \
	   --all-sites
fi

########### Filtering
	echo "Excluding indels and low quality SNPs.."
	gatk SelectVariants \
		-R $ref \
		-V ${out}.all_variants.vcf.gz \
		--select-type-to-include SNP \
		-O ${out}.all_snps.vcf.gz
	echo "Variant filtering..."
	########### variant filtering based on https://gatk.broadinstitute.org/hc/en-us/articles/360035532412-Can-t-use-VQSR-on-non-model-organism-or-small-dataset
	##### do not remove filtered snps, just flag them
	gatk VariantFiltration \
		-R $ref \
		-V ${out}.all_snps.vcf.gz \
		-O ${out}.snps_qc.vcf.gz 	\
		--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 ||  MappingQualityRankSum < -12.5" \
		--filter-name "gatk_filter"
	
if [ $fhet == "yes" ]; then
	echo "Filtering heterozygous calls"

	gatk VariantFiltration \
		-R $ref \
		-V ${out}.snps_qc.vcf.gz \
		-O ${out}.snps_qc_fhetVF.vcf.gz \
		--genotype-filter-expression "isHet == 1" \
		--genotype-filter-name "isHetFilter"
	
	#### Convert heterozygous calls to N
	gatk SelectVariants \
		-V ${out}.snps_qc_fhetVF.vcf.gz \
		--set-filtered-gt-to-nocall \
		-O ${out}.snps_qc_fhetSV.vcf.gz


	soft=${out}.snps_qc_fhetSV.vcf.gz

elif [ $fhet == "no" ]; then
       echo "Not filtering heterozygous calls" 
	soft=${out}.snps_qc.vcf.gz
else 
	echo "Filter heterozygous calls value not valid"	
	echo "Not filtering heterozygous calls" 
	soft=${out}.snps_qc.vcf.gz

fi
	
echo "Hard filtering: retain only filter=PASS, biallelic and min ac =1"
bcftools view -v snps -m2 -M2 -i 'F_MISSING<0.1' -f 'PASS,.' --min-ac 1:minor -I $soft -O z -o $(basename $soft .vcf.gz).bial-mac1.hard.vcf.gz
	tabix $(basename $soft .vcf.gz).bial-mac1.hard.vcf.gz




