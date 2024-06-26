#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $last_seq_id = 0;
my $fetching_seq = 0;
my $seq = undef;
my $seq_id = undef;
my %seq_lengths = ();

my $print_usage = 0;

my $usage = <<USAGE;

  This script converts GFF formatted annotation into a NCBI "Feature Table" 
  for submission of sequences to Genbank.

  It expects GFF3 formatted input (http://gmod.org/wiki/GFF3#GFF3_Annotation_Section)
  This GFF3 file should not contain sequences. 

  It's been tested with GFF3 files exported from Geneious v2023.2 (without sequences).

  Reads from stdin and writes to stdout.

  Mark Stenglein,  3/6/2018
  Modified Beatriz Navarro, 4/25/2024

  Usage: $0 [-h] 

   [-h]          print this message

USAGE

if ((scalar @ARGV == 0) and -t STDIN) { print $usage and exit; }

GetOptions ("h" => \$print_usage);

# read stdin
while (<>)
{
  chomp;

  # try to parse out sequence length info, if present
  if (/^##/)
  {
     # read in sequence length if present (this information is present in Geneious-exported GFF3)
     if (/^##sequence-region/)
     {
	my @fields = split "\t";
        my $seq_id = $fields[1];
        my $seq_length = $fields[3];
        $seq_lengths{$seq_id} = $seq_length;
     }
   }
   else
   {
      chomp;
      # GFF is tab-delimited, split into fields
      my @fields = split "\t";
      my $seq_id = $fields[0];
      my $type = $fields[2];
      my $start = $fields[3];
      my $end = $fields[4];
      my $strand = $fields[6];
      my $name = $fields[8];
      my $partial=0;

      if ($strand eq '-')
      {
         die "script doesn't support -strand features yet\nline: $_\n";
      }

      # first line for a entry for a new sequence
      if ($seq_id ne $last_seq_id)
      {
         $last_seq_id = $seq_id;
         print ">Feature $seq_id\n";
      }

      # match any text up to a ; character
      if ($name =~ /Name=([^;]+)/)
      {
         $name = $1;
	 # strip "CDS" or "_CDS" or "gene" from name
         $name =~ s/_CDS//;
         $name =~ s/ CDS//;
         $name =~ s/ gene//;
	 # if name contains 
         if ($name =~ /partial/)
         {
            $partial = 1;
         }
      }
      else
      {
         die "error: couldn't parse feature name from $name\nline: $_\n";
      }

      # handle CDS or mat_peptide feature types
      if (($type eq "CDS") or ($type eq "tRNA") or ($type eq "gene") or ($type eq "rRNA"))
      {
         my $seq_length = length($seq_lengths{$seq_id});

	 # handle CDS that extend outside of sequence (
	 # CDS that start at beginning of sequence (position 1)
         if ($start == 1) 
	 { 
	    print "<".$start."\t"; 
	 }
         else             
	 { 
	    # CDS starts within sequence
	    print  "$start\t"; 
	 }
	 # CDS that extend past end of sequence 
	 # (only possible to know if sequence length defined)
         if (defined $seq_length and ($end == $seq_length))
	 { 
	    print $end.">\t"; 
	 }
         else             
	 { 
	    # CDS ends within sequence
	    print  "$end\t"; 
	 }
	 

	 # type of feature
         print "$type\n";

	 # genes
	 if ($type eq "gene")
	 {
	 print "\t\t\tgene\t$name\n"
	 }
	 else
	 {
	 # name of product
         print "\t\t\tproduct\t$name\n";
         }
	 # if the CDS name includes "partial" make a note of that
         if ($partial) 
         { 
            print "\t\t\tnote\tpartial\n";
         }
      }
      else
      {
         warn "unsupported feature type: $type\nline: $_\n";
      }
   }
}

