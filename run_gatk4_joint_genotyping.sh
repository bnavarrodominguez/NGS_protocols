#!/bin/bash
if [ "$#" -ne 6 ]; then
		    echo "Usage: $0 map_file.txt reference.fasta out_prefix interval[chr or chr:start-stop] ploidy nthreads"
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


## Make database

if [ -d ${out}_database ]; then
	echo "${out}_database already exists; moving on ..." 
else
echo "Importing interval $region to ${out}_database"
	gatk --java-options "-Xmx16g -Xms16g" GenomicsDBImport \
	            --sample-name-map $files \
	                --genomicsdb-workspace-path ${out}_database \
	                --reader-threads $nthr \
			--intervals $region
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
	   -all-sites
fi

########### Filtering
if [ -f ${out}.filtered.vcf.gz ]; then
	echo "${out}.filtered_snps.vcf.gz already exists, moving on ..."
else
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
		-O ${out}.filtered_snps.vcf.gz 	\
		--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 ||  MappingQualityRankSum < -12.5" \
		--filter-name "gatk_filter"
	###### retain only biallelic snps with a min depth of 4 and max depth of 200 and where the minor allele is present in at least one individual (discard "SNPS" where all individuals are different from the reference, yet equal to each other"
	bcftools view -v snps -m2 -M2 --min-ac 1:minor -I ${out}.filtered_snps.vcf.gz -O z -o ${out}.bialelic_mac1_filtered_snps.vcf.gz
	tabix ${out}.bialelic_mac1_filtered_snps.vcf.gz



fi
