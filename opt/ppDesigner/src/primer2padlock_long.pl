#!/usr/bin/perl -w
# usage:
# primer2padlock_long [maxH1H2len] < [primers.txt] > [probes.seq]


use strict;
my %cmpTable = ("A", "T", "T", "A", "G", "C", "C", "G","R","Y","Y", "R", "M", "K", "K", "M", "S", "S", "W", "W",  "N", "N");
my $AP1 = "GTCATATCGGTCACTGTT"; #V6s
my $AP2 = "GATCAGGATACACACTACCC"; #V6s
my $linkerL = "GTTGGAGGCTCATCGTTCCTATT";
my $linkerR = "CAGGCAGATGTTATCGAGGTCCGAC"; #V8
my $maxH1H2Len = $ARGV[0];
my $longProbes=0;
my @bases=("A", "T", "C", "G");

while(my $line = <STDIN>){
	next if($line !~ /chr/ || $line =~ /good/);
	chop($line);
	my @fields = split(/[\t ]+/, $line);
	my $LP = $fields[5];
	$LP =~ s/C/Y/ig;
	my $RP = $fields[7];
	$RP =~ s/C/Y/ig;
	my $DmrID= $fields[0];
	my ($chr,$chr_start,$chr_end) = split(/[:\-]/, $fields[1]);
	my $offset = $fields[2];
	my $target_start = $fields[3]+$offset+$chr_start;
	my $target_end = $fields[4]+$offset+$chr_start;
	my $len = length($LP.$RP);

	if($len > $maxH1H2Len){
		$longProbes++;
		next;
	}
	my $fill_in="";
	while($len < $maxH1H2Len){
		$fill_in = $fill_in . "C";
		$len ++; 
	}
	my $padlock = $AP1. &RevComp($LP) . $linkerL. $fill_in . $linkerR . &RevComp($RP) . $AP2;
	my @fragments = split(/R/, $padlock);
	my @variants;
	if(scalar(@fragments)<4){
		push(@variants, $fragments[0]);
		for(my $i=0; $i<scalar(@fragments)-1; $i++){
			my $size = scalar(@variants);
			for(my $j=0; $j < $size; $j++){
				push(@variants, $variants[$j] . 'G' . $fragments[$i+1]);
			}
			for(my $j=0; $j < $size; $j++){
				$variants[$j] = $variants[$j] . 'A' . $fragments[$i+1];
			}
		}
	}else{
		$variants[0] = $padlock;
		$variants[0] =~ s/R/G/g;
		$variants[1] = $padlock;
		$variants[1] =~ s/R/A/g;
	}
#	print $line, "\n";
    for(my $i=1; $i<=scalar(@variants); $i++){
	    print "$DmrID:$chr:$target_start-$target_end:$i\t", $variants[$i-1], "\n";
    }
}
warn "$longProbes probes are too long and not included.\n";

#-------------------
sub RevComp(){
#-------------------
    my $seq = shift;
	my $seqLen = length($seq);
	my $revcom = $cmpTable{substr($seq,0,1)};
	for(my $i = 1; $i < $seqLen; $i++){
		$revcom = $cmpTable{substr($seq, $i, 1)} . $revcom;
	}
    return $revcom;
}
