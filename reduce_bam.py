#!/usr/bin/python

from subprocess import call
import sys

print "Usage: reduce_bam.py BamFile"

try:
    bamfile = sys.argv[1]
except:
    bamfile = raw_input("Introduce name of your bam file: ")

name = bamfile.split(".")
name = name[0]

print "\nReducing file %s..." % (bamfile)

#to get all the reads where both mapped.
call("samtools view -b -F 12 %s > %s_temp1.bam" % (bamfile,name), shell=True)

#to get all the reads that did not map, but whose mate mapped
call("samtools view -b -f 4 -F 8 %s > %s_temp2.bam" % (bamfile,name), shell=True)

#to get all the reads that mapped, but whose mates did not.
call("samtools view -b -f 8 -F 4 %s > %s_temp3.bam" % (bamfile,name), shell=True)

#to merge bam files and sort by name
#call("samtools merge -u - temp[123].bam | samtools sort - %s_mapped" % (name), shell=True)
call("samtools merge -u - %s_temp[123].bam | samtools sort -T %s_aln.sorted - -o %s_mapped.bam" % (name,name,name), shell=True)


#to index sorted bam file
call("samtools index %s_mapped.bam" % (name), shell=True)

call("rm %s_temp1.bam %s_temp2.bam %s_temp3.bam %s_aln.sorted" % (name,name,name), shell=True)
