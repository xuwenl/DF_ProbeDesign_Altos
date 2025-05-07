#!/usr/bin/perl -w

#Copyright 2009, 2010 The Regents of the University of California.  All Rights Reserved.

#Output data to screen.  This format can be easily converted to oligos.

#Statistical variables
our $totalTargetSize = 0;
our $totalUniqueTargetSize = 0;

sub output_text()
{
	my @probeInfo = shift;
	my $time_used = shift;
	my $exonExid = shift;
	my $chr = shift;
	my $exonStart = shift;
	my $exonEnd = shift;

	#Were any probes found?  If so, output their data
	if($probeInfo[0])
	{
		#Calculate the % covereage for this particular target (useful to know)
		my $searchStart = $flankingLen;
		my $searchEnd = $gSeqLen - $flankingLen;
		if($searchEnd <= $searchStart) {
			$searchEnd = $searchStart + 1;
		}
		my $setCoverage = &calcSetCoverage(\@probeInfo, $searchStart, $searchEnd);
		$totalCapturedTargets += $setCoverage;
		my $pctCoverred = sprintf("%3.1f",$setCoverage/($searchEnd-$searchStart)*100);
		$pctCoverred = 100 if($pctCoverred>100);
		foreach my $probe (@probeInfo)
		{		
			print "$exonExid\t$chr:$exonStart-$exonEnd\t", -$flankingLen,"\t$probe\t$pctCoverred\t",$time_used, "s\n";			
		}
	}
	else
	{
		#No probes
		print "$exonExid\t$chr:$exonStart-$exonEnd\t#no good probes found.\n";
	}
}

sub output_statistics()
{
	print "total_exon_size=$totalTargetSize\ttotal_captured_size=$totalCapturedTargets\n";
}

1;
