#!/usr/bin/perl -w

#Copyright 2009, 2010 The Regents of the University of California.  All Rights Reserved.

#Parameters for probe design.  Replace with desired parameter values.

our $gSeqLen = 0;

our $last_chr;

our %cmpTable = ("A", "T", "T", "A", "G", "C", "C", "G", "N", "N", "M", "K", "R", "Y", "W", "W", "S", "S", "Y", "R", "K", "M", "V", "B", "H", "D", "D", "H", "B", "V");

our @candidateProbeSets;

#Some constants
our $R = 1.987;
our %dH_buffer;
our %dS_buffer;
our %deltaH;
our %deltaS;
our @kMerDist;
our %oligoFreqBuffer;
our ($exonStart, $exonEnd);

&readNNParam($NNParamFile);

sub readNNParam(){
	my $NNParamFile = shift;
	open(NNFILE, $NNParamFile) || die("Error in opening NN parameter file!");
	my $line = <NNFILE>;
	while($line = <NNFILE>){
		$line =~ s/\r//;
		chop($line);
		my ($seqF, $seqR, $dH, $dS, $mismatch) = split(/[:\t]/, $line);
		if(!$mismatch){
			$deltaH{'pm'}->{$seqF} = $dH;
			$deltaS{'pm'}->{$seqF} = $dS;
		}else{
			$deltaH{'mm'}->{$seqF}->{$seqR} = $dH;
			$deltaS{'mm'}->{$seqF}->{$seqR} = $dS;			
		}
	}
	close(NNFILE);
}

1;
