#!/bin/bash

display_usage() {
<<<<<<< HEAD
       echo -e "\nUsage:$0 LibraryName adapters.fa threads \n
	Check adapters files in /usr/local/lib/Trimmomatic-0.33/adapters/ \n"
        }

if [  $# -ne 3 ]
=======
        echo -e "\nUsage:$0 LibraryName\n" 
        }

if [  $# -ne 1 ]
>>>>>>> 9d14234d6600afb717a4bca61ffaf5a15be1c988
        then
                display_usage
                exit 1
        else

                lib=$1
<<<<<<< HEAD
                thr=$3
=======
                #thr=$2
>>>>>>> 9d14234d6600afb717a4bca61ffaf5a15be1c988
                reads1=${lib}_1.fastq.gz
                #reads1=${lib}_1.fq.gz
                reads2=${lib}_2.fastq.gz
                #reads2=${lib}_2.fq.gz
<<<<<<< HEAD
		adapters=$2
        fi

if [ ! -f TruSeq3-PE.fa ]; then
cp $adapters adapters.fa
fi

trimmomatic PE -threads 12 $reads1 $reads2 ${lib}_paired_1.fastq.gz ${lib}_unpaired_1.fastq.gz ${lib}_paired_2.fastq.gz ${lib}_unpaired_2.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:10:30 MINLEN:100
#
#


=======
        fi

#if [ ! -f TruSeq3-PE.fa ]; then
#        wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
#fi

adapters=/usr/local/lib/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa

java -jar /usr/local/lib/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 $reads1 $reads2 ${lib}_paired_1.fastq.gz ${lib}_unpaired_1.fastq.gz ${lib}_paired_2.fastq.gz ${lib}_unpaired_2.fastq.gz ILLUMINACLIP:$adapters:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:10:30 MINLEN:100
>>>>>>> 9d14234d6600afb717a4bca61ffaf5a15be1c988

fastqc ${lib}_paired_1.fastq.gz ${lib}_paired_2.fastq.gz
