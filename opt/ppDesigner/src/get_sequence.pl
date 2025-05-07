#!/usr/bin/perl -w

#Copyright 2009, 2010 The Regents of the University of California.  All Rights Reserved.

#get_sequence.pl
#Extracts a nucleotide sequence.  Also keeps track of the last chromosome used in the global variable "$last_chr"
#Uses global variable $flankinglen.  Returns the basic sequence only.

$last_chr = 'chr';
our ($chrSeqName, $chrSeq);

#Obtain a genomic sequence for probe creation
sub get_sequence()
{
	my $chrom = shift;
	my $start = shift;
	my $end = shift;

	if ($chrom ne $last_chr)
	{
		($chrSeqName, $chrSeq) = &read_fa($HsDir . $chrom.".fa", 1);
		#if(!$chrSeq || length($chrSeq) < 10)
		#{
		#	($chrSeqName, $chrSeq) = &read_fa($HsDir . $chrom.".fa", 1);
		#}
		#if(!$chrSeq || length($chrSeq) < 10)
		#{
		#	($chrSeqName, $chrSeq) = &read_fa($HsDir . $chrom.".fa.masked", 1);
		#}
		die("Very small chromosome detected.  FASTA file $chrom sequence most likely not correctly loaded!\n") if(!$chrSeq || length($chrSeq) < 100);
		$last_chr = $chrom;
	}

	my $substrstart;
	if ($start - $flankingLen < 0)
	{
		$substrstart = 0;
	}
	else
	{
		$substrstart = $start-$flankingLen;
	} 
	
	my $gSeq = substr($chrSeq,$substrstart, $end-$start+1+2*$flankingLen);
	
	#lower case letters are repeats, replaced to Ns
	#$gSeq =~ s/[atgc]/N/g;

	$gSeqLen = length($gSeq);
	return $gSeq;
}

#Read the nth sequence from a Fasta file.
sub read_fa(){
        my $fileName = shift;
        my $seqNumber = shift;
        open(FA_IN, "$fileName") || return("Error in opening the Fasta file $fileName!");
        my ($seqName, $seq, $nSeq);
        $nSeq = 0;
        while(my $line = <FA_IN>){
                chop($line);
                if($line =~ /^>/){
                        if($nSeq == $seqNumber){
                                close(FA_IN);
                                return ($seqName, $seq);
                        }
                        $seqName = substr($line,1,length($line)-1);
                        $seq = "";
                        $nSeq++;
                }elsif($line =~ /[ATGCN]/i){
                        $line =~ s/[\t ]+//;
                        $seq = $seq . $line;
                }
        }
        close(FA_IN);
        return ($seqName, $seq);
}

1;
