#!/usr/bin/perl -w

#USAGE: ppDesigner.pl targetInfo.txt > primers.txt

#Copyright 2009, 2010 The Regents of the University of California.  All Rights Reserved.

use File::Temp qw/ tempfile tempdir /;
use Sort::Array qw(
      Sort_Table
      Discard_Duplicates
  );


my $input_file = $ARGV[0];

eval `cat $input_file` or die "couldn't parse input file";

$| = 1;

#Pick input fasta file (optional)
my $input_chr = $ARGV[1];

#Load targets
my %targetTable = load_inputs($exon_info_file);

our $maxCG = 0 if(!$maxCG);
	
#Cycle through each FASTA file 
foreach my $target_chrom (sort keys(%targetTable))
{
	#Have we been told to only run a single FASTA file? If so, are we on this one?
	next if($input_chr && $input_chr ne $target_chrom);
	
	#Cycle through each site
	foreach my $target_start (sort keys(%{$targetTable{$target_chrom}}))
	{
		my $target_end = $targetTable{$target_chrom}->{$target_start}->{'end'};
		my $target_id = $targetTable{$target_chrom}->{$target_start}->{'exid'};
		my $target_strand = $targetTable{$target_chrom}->{$target_start}->{'use_strand'};
		#Obtain reference sequence for the target
		my $target_sequence = get_sequence($target_chrom, $target_start, $target_end);
		$last_chr = $target_chrom;

		#Get probes to capture the target
		my $time_start = time();
		#print "$target_id\t$target_chrom:$target_start-$target_end\n";
		my @probeInfo = get_probes($target_sequence,$target_start,$target_end,$target_strand);
		my $time_used = time() - $time_start;
		$last_chr = $target_chrom;
		#Output probe data
		if(scalar(@probeInfo) == 0)
		{
				output_text('',$time_used, $target_id, $target_chrom,$target_start, $target_end);
		}
		else
		{
			for(my $index=0;$index<scalar(@probeInfo);$index++)
			{
				output_text($probeInfo[$index],$time_used, $target_id, $target_chrom,$target_start, $target_end);
			}
		}
	}
}
#Output probe statistics
output_statistics();
