#!/bin/bash 
#####Dependencies: sratoolkit, gnu-parallel
 
display_usage() { 
	echo "Parallel download of fastq.gz from a list of SRA accession numbers" 
	echo -e "\nUsage:$0 SRA_accesions.txt \n" 
	}

if [  $# -ne 1 ] 
	then 
		display_usage
		exit 1
	fi 
 

accno=$1


while IFS='' read -r line || [[ -n "$line" ]]; do 
	if [ -s fastq/$line"_pass_1.fastq.gz" ]
	then
		echo "# $line already downladed" 
		#exit 1
	else
	echo "fastq-dump --gzip --readids --origfmt --skip-technical --clip --read-filter pass --split-3 $line";  
	#i=$((i+1)); 
	fi
done < $accno > commands_"$accno".txt; 

parallel < commands_"$accno".txt
