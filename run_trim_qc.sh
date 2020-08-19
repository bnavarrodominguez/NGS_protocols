#!/bin/bash

display_usage() {
        echo -e "\nUsage:$0 LibraryName threads \n" 
        }

if [  $# -ne 2 ]
        then
                display_usage
                exit 1
        else

                lib=$1
                thr=$2
                reads1=${lib}_1.fastq.gz
                reads2=${lib}_2.fastq.gz
        fi

if [ ! -f TruSeq3-PE.fa ]; then
        wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
fi



trimmomatic PE -threads 12 $reads1 $reads2 ${lib}_paired_1.fastq.gz ${lib}_unpaired_1.fastq.gz ${lib}_paired_2.fastq.gz ${lib}_unpaired_2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100


fastqc ${lib}_paired_1.fastq.gz ${lib}_paired_2.fastq.gz
