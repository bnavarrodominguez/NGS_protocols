#!/bin/bash

display_usage() {
       echo -e "\nUsage:$0 LibraryName fastq.gz/fq.gz adapters.fa threads \n
	Check adapters files in /usr/local/lib/Trimmomatic-0.33/adapters/ \n"
        }


if [  $# -ne 4 ]
        then
                display_usage
                exit 1
        else

                lib=$1
                thr=$4
		suf=$2
                reads1=${lib}_1.${suf}
                #reads1=${lib}_1.fq.gz
                reads2=${lib}_2.${suf}
                #reads2=${lib}_2.fq.gz
		adapters=$3
        fi


#cp $adapters adapters.fa
#echo "
trimmomatic PE -threads 12 $reads1 $reads2 ${lib}_paired_1.fastq.gz ${lib}_unpaired_1.fastq.gz ${lib}_paired_2.fastq.gz ${lib}_unpaired_2.fastq.gz ILLUMINACLIP:${adapters}:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:10:30 MINLEN:100

fastqc ${lib}_paired_1.fastq.gz ${lib}_paired_2.fastq.gz


#"


