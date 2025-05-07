#!/usr/bin/perl -w

#Copyright 2009, 2010 The Regents of the University of California.  All Rights Reserved.

#loadInputs.pl
#Loads a standard-format input file into a hash table.
#Input file format:
#exonID	chr	start	end	strand

sub load_inputs(){
	#Open the target file
	my $exon_info_file = shift;
	open(INFILE, "$exon_info_file") || die("Error in opening gene file $exon_info_file!");
	while(my $line = <INFILE>){
	        chop($line);
		my @fields = split(/[ \t]+/, $line);
		
		next if($fields[1] =~ /random/);

		#Skip if not the correct fasta file
		next if($target_chr && $target_chr ne $fields[1]);
		my $start;	

		$start = sprintf("%09d", $fields[2]);#-$splicingRegions+1);
		$exonTable{$fields[1]}->{$start}->{'end'}=$fields[3];#+$splicingRegions;
		$exonTable{$fields[1]}->{$start}->{'exid'}=$fields[0];
		if (($fields[4]) && ($fields[4] eq '+' || $fields[4] eq '-'))
		{
			$exonTable{$fields[1]}->{$start}->{'use_strand'}=$fields[4];
		}
		else
		{
			$exonTable{$fields[1]}->{$start}->{'use_strand'}=' ';
		}
	}
	close(INFILE);
	return %exonTable;
}

1;
