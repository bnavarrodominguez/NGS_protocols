#!/bin/bash

### parameters from https://github.com/rancilhac/Museoscript/blob/master/museoscript.sh

display_usage() {
        echo -e "\nUsage:$0 LibraryName\n" 
        }

if [  $# -ne 1 ]
        then
                display_usage
                exit 1
        else

                lib=$1
                #thr=$2
                reads1=${lib}_1.fastq.gz
                #reads1=${lib}_1.fq.gz
                reads2=${lib}_2.fastq.gz
                #reads2=${lib}_2.fq.gz
        fi

#if [ ! -f TruSeq3-PE.fa ]; then
#        wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
#fi

adapters=/usr/local/lib/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa

java -jar /usr/local/lib/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 $reads1 $reads2 ${lib}_paired_1.fastq.gz ${lib}_unpaired_1.fastq.gz ${lib}_paired_2.fastq.gz ${lib}_unpaired_2.fastq.gz ILLUMINACLIP:$adapters:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:2:25 MINLEN:36

fastqc ${lib}_paired_1.fastq.gz ${lib}_paired_2.fastq.gz
