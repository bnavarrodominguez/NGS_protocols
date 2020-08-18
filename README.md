<h2>NGS protocols</h2>

General scripts for processing and analyzing NGS data
* **parallel_fastqdump.sh**: download fastq.gz files in parallel from a list of NCBI SRA accessions. Dependencies: sratoolkit, gnu-parallel
* **reduce_bam.py**: keeps only mapped reads and their mates (paired or unpaired) in a BAM file. Modified from https://github.com/fjruizruano/ngs-protocols/blob/master/reduce_bam.py 
