#!/usr/bin/perl -w

#Copyright 2009, 2010 The Regents of the University of California.  All Rights Reserved.

our %proxLigDincTable = ("TA", 9, "GA", 2, "AA", 7, "TG", 1, "GG", 6, "AG", 4, "AT", 8, "GT", 3, "TT", 5);
our %proxExtDincTable = ("GA", 1, "TA", 5, "GT", 8, "AT", 7, "AA", 2, "GG", 3, "TG", 4, "AG", 9, "TT", 6);

sub get_probes()
{
	my ($target_sequence,$target_start,$target_end,$use_strand) = @_;

	my @probeInfo;
	
	$gSeqLen = length($target_sequence);
	
	#Calculate some target statistics
	$totalTargetSize += $target_end-$target_start+1;
	
	#If the target is more than 90% N's, reject it.
	if ((uc($target_sequence) =~ tr/N//)/length($target_sequence) >= 0.9)
	{
		#Return an empty probeset
		return @probeInfo;
	}
	
	# print "In get_probes().\n";
	#Small target?  If so, we just need one probe.
	if($target_end-$target_start+1 <= $targetMaxLen) #160
	{
		# print "In small targets.\n";
		my $searchStart = $flankingLen < $gSeqLen-$targetMinLen ? $flankingLen : $gSeqLen-$targetMinLen;
		my $searchEnd = $gSeqLen-$flankingLen > $targetMinLen ? $gSeqLen-$flankingLen : $targetMinLen;
		if($arewebisulfite == 1)
		{
			my $minCG = 0;
			while($minCG <= $maxCG and scalar(@probeInfo) == 0)
			{
				@probeInfo = &getProbes_smallExon($target_sequence, $searchStart, $searchEnd,$use_strand, $minCG);
				$minCG++;
			}
		}else{
			@probeInfo = &getProbes_smallExon($target_sequence, $searchStart, $searchEnd,$use_strand, 0);
		}
	}
	elsif($target_end-$target_start+1 > $largeExonSize)
	{
		#Large target?  If so, we need to split it, or the search will take an extremely long time.
		my $seed_probe="";		
		my $frag_start = 0;
		my $old_frag_start = 0;
		do
		{
			#Break the target into a chunk and find probes.
			my $subSeq = substr($target_sequence,$frag_start,$largeExonSize);
			my $subSeqLen =length($subSeq);
			my $searchStart = $seed_probe =~ /[a-z0-9]/i ? 0: $flankingLen ;
			my $searchEnd = $subSeqLen-$flankingLen; # ($subSeqLen == $largeExonSize) ? $subSeqLen-80 :  $subSeqLen-$flankingLen; <-221007: Kian
			my @fragProbeInfo;
			if($arewebisulfite == 1)
			{
				my $minCG = 0;
				while($minCG <= $maxCG and scalar(@fragProbeInfo) == 0)
				{
					@fragProbeInfo = &getProbes($subSeq, $seed_probe, $searchStart, $searchEnd,$use_strand, $minCG);
					$minCG++;
				}
			}else{
				@fragProbeInfo = &getProbes($subSeq, $seed_probe, $searchStart, $searchEnd,$use_strand, 0);
			}
			$seed_probe='';

			#Obtain the last probe in the set and keep it for the next chunk (so we don't get overlap)
			for(my $k=0; $k<scalar(@fragProbeInfo); $k++)
			{
				my @probe_fields = split(/\t/, $fragProbeInfo[$k]);
				my ($targetStart,$b_start,$LP,$TmLP,$RP,$TmRP,$score,$strand,$target)= @probe_fields;
				#adjust the offset of probe position.
				$probe_fields[0] += $frag_start;
				$probe_fields[1] += $frag_start;
				push(@probeInfo, join("\t", @probe_fields));
				#Use the last probe to seed the search in the next fragment
				if($k ==scalar(@fragProbeInfo)-1){
					my ($a2_start, $a2_len, $a2_Tm, $a2_oligo, $b2_len, $b2_start, $b2_Tm, $b2_oligo);
					#Need to consider strand for seeding procedure
					if($strand =~ /w/i)
					{
						$a2_len = length($LP);
						$b2_len = length($RP);
						$a2_start = -$a2_len;					
						$a2_Tm = $TmLP;
						$a2_oligo=$LP;
						$b2_Tm = $TmRP;
						$b2_oligo=$RP;
					}
					else
					{
						$a2_len = length($RP);
						$b2_len = length($LP);
						$a2_start = -$a2_len;
						$a2_Tm = $TmRP;
						$a2_oligo=$RP;
						$b2_Tm = $TmLP;
						$b2_oligo=$LP;
					}
					$b2_start = $b_start - $targetStart;
					#Generate a seed probe output from the last probe in the probeset
					$seed_probe = "$a2_start,$a2_len,$a2_Tm,$a2_oligo,$b2_start,$b2_len,$b2_Tm,$b2_oligo,$score,$strand";
					$frag_start += $targetStart+$jumpFromPrevious; # Kian added 221007. Previously +80
				}							
			}
		
			if($frag_start == $old_frag_start)
			{
				$frag_start += $jumpIf0Found; # Kian added 221007. Previously +80
			}
			$old_frag_start = $frag_start;
		
		} while(($frag_start+$targetMinLen < $target_end-$target_start+2*$flankingLen));			
		#Keep finding new probes and storing them in probeInfo until we have enough
	}
	else
	{
		#Medium sized probe; get probes that are highest scoring and maximize coverage
		my $searchStart = $flankingLen;
		my $searchEnd = $gSeqLen - $flankingLen;
		if($arewebisulfite == 1){
			my $minCG = 0;
			while($minCG <= $maxCG and scalar(@probeInfo) == 0)
			{
				@probeInfo = &getProbes($target_sequence, '', $searchStart, $searchEnd,$use_strand, $minCG);
				$minCG++;
			}
		}else
		{
			@probeInfo = &getProbes($target_sequence, '', $searchStart, $searchEnd,$use_strand, 0);		
		}
	}

	return @probeInfo;
}

#Recursively generates probe sets.
sub compSet{
	my $h_input_set = shift;
	my $h_allProbeList = shift;
	my $startvalue = shift;
	my $max_setsize = shift;
	my $regionStart = shift;
	my $regionEnd =shift;
	my $regionLen = shift;	
	
	my $foundone = 0;	
	if (scalar(@{$h_input_set}) < $max_setsize)
	{
		my $i=-1;
		my @input_set = @{$h_input_set};
		foreach my $this_probe (@{$h_allProbeList})
		{
			$i++;
			next if ($i<$startvalue);
			$foundone =isCompatibleProbes(\@input_set,$this_probe);
			next if($foundone != 1);
			push(@{$h_input_set}, $this_probe );
			compSet($h_input_set, $h_allProbeList, $startvalue+1, $max_setsize, $regionStart, $regionEnd, $regionLen);
			@{$h_input_set}=@input_set;
		}
	}
	my $setScore = sprintf("%4.1f", &calcSetScore($h_input_set, $regionStart, $regionEnd, $regionLen));
	push(@main::candidateProbeSets, join(";", @{$h_input_set}) . "\t$setScore");
}

sub getProbes_smallExon(){
	my $fwdSeq = shift;
	my $regionStart = shift; #start position of the region where probes will be designed.
	my $regionEnd = shift; #end position of the region where probes will be designed.
	my $regionLen = length($fwdSeq);

	my $revSeq = &RevComp($fwdSeq);

	my $use_strand = shift;
	my $minCG = shift;	
	#If we are creating probes for methylation analysis, we need to convert C->T (except in the case of CG)
	if($arewebisulfite == 1)
	{
		$fwdSeq = bisulfiteCvt($fwdSeq);
		$revSeq = bisulfiteCvt($revSeq);
	}

	my $bestProbeRev;
	my $bestProbeFwd;
	
	if($use_strand eq '+')
	{
		$bestProbeFwd=&getProbeList_smallExon($fwdSeq, 'W',$regionStart,$regionEnd, $minCG);	
	}
	elsif($use_strand eq '-')
	{
		$bestProbeRev=&getProbeList_smallExon($revSeq, 'C',$regionLen-$regionEnd,$regionLen-$regionStart, $minCG);	
	}
	else
	{	
		$bestProbeFwd=&getProbeList_smallExon($fwdSeq, 'W',$regionStart,$regionEnd, $minCG);	
	
		$bestProbeRev=&getProbeList_smallExon($revSeq, 'C',$regionLen-$regionEnd,$regionLen-$regionStart, $minCG);
	}
	my @probeListAll;
	push(@probeListAll,$bestProbeFwd);
	push(@probeListAll,$bestProbeRev);

	my @sorted_probeList = Sort_Table(
                            	cols      => '10',
                                field     => '9',
                                sorting   => 'descending',
                                structure => 'csv',
                                separator => ',',
                                data      => \@probeListAll);

   	my $seqLen = length($revSeq);
   	my @probeInfo;	
	my $probe = $sorted_probeList[0];
	return if(!$probe || $probe !~/[0-9ATGCW]/);
	   	my ($a_start, $a_len, $a_Tm, $a_oligo, $b_start, $b_len, $b_Tm, $b_oligo, $score, $strand) = split(/,/, $probe);
	   	my ($LP, $TmLP, $RP, $TmRP, $target);
		if($strand =~ /w/i){
			$LP = $a_oligo;
			$RP = $b_oligo;
			$target = substr($fwdSeq, $a_start, $b_start-$a_start+$b_len);
			$TmLP = $a_Tm;
			$TmRP = $b_Tm;
		}else{
			$LP = $b_oligo;
			$RP = $a_oligo;
			$target = substr($revSeq, $seqLen-$b_start-$b_len, $b_start-$a_start+$b_len);
			$TmLP = $b_Tm;
			$TmRP = $a_Tm;
		}
		my $targetStart = $a_start+$a_len;
	push(@probeInfo,  "$targetStart\t$b_start\t$LP\t$TmLP\t$RP\t$TmRP\t$score\t$strand\t$target");
	return @probeInfo;
}

# Find a set of probes that maximizes the coverage of the input sequence. 
sub getProbes(){
	my $fwdSeq = shift;
	my $seed_probe = shift;
	my $regionStart = shift; #start position of the region where probes will be design.
	my $regionEnd = shift; #end position of the region where probes will be design.
	return if($regionEnd-$regionStart<20);
	
	my $regionLen = length($fwdSeq);	
	my $revSeq = &RevComp($fwdSeq);

	my $use_strand = shift;
	my $minCG = shift;	
	#If we are creating probes for methylation analysis, we need to convert C->T (except in the case of CG)
	if($arewebisulfite == 1)
	{
		$fwdSeq = bisulfiteCvt($fwdSeq);
		$revSeq = bisulfiteCvt($revSeq);
	}
	
	my @probeListRev;
	my @probeListFwd;
	
	if($use_strand eq '+')
	{
		@probeListFwd=&getProbeList($fwdSeq, 'w', $minCG);
	}
	elsif($use_strand eq '-')
	{
		@probeListRev=&getProbeList($revSeq, 'c', $minCG);
	}
	else
	{	
		@probeListFwd=&getProbeList($fwdSeq, 'w', $minCG);
	
		@probeListRev=&getProbeList($revSeq, 'c', $minCG);	
	}

	return if(scalar(@probeListFwd)<1 && scalar(@probeListRev) < 1);

	undef @candidateProbeSets;
	
	my @probeListAll;	
	push(@probeListAll,@probeListFwd);
	push(@probeListAll,@probeListRev);
	undef(@probeListFwd);
	undef(@probeListRev);	
	return if (scalar(@probeListAll) == 0);
	my $target_sized = length($fwdSeq) - 2*$flankingLen;
	my $max_probeSetSize = int($target_sized*1.6 / $targetMinLen) + 1;
	$max_probeSetSize ++ if($seed_probe);

	for(my $i=0; $i < scalar(@probeListAll); $i++){
		my @list;		
		next if(!$probeListAll[$i]);
		my $isCompatible; 
		if($seed_probe) {
			push (@list, $seed_probe);
			$isCompatible = isCompatibleProbes(\@list,$probeListAll[$i]);
		}		
		if(scalar(@list)<1 || $isCompatible==1){
			push (@list, $probeListAll[$i]);	
			&compSet(\@list, \@probeListAll, $i+1,$max_probeSetSize, $regionStart, $regionEnd, $regionLen);
		}
		undef(@list);
	}

	#Create sets of probes
	return if(scalar(@candidateProbeSets)<1);
	
	my @sorted_probeSets = Sort_Table(
                            	cols      => '2',
                                field     => '2',
                                sorting   => 'descending',
                                structure => 'csv',
                                separator => '\t',
                                data      => \@candidateProbeSets);	
    #extract primers
   	my ($bestSet, $bestScore) = split(/\t/, $sorted_probeSets[0]);
	undef(@sorted_probeSets);
   	my @bestProbes = split(/;/, $bestSet);
   	my $seqLen = length($revSeq);
   	my @probeInfo;
	foreach my $probe (@bestProbes){
		next if($seed_probe && $probe eq $seed_probe);
	   	my ($a_start, $a_len, $a_Tm, $a_oligo, $b_start, $b_len, $b_Tm, $b_oligo, $score, $strand) = split(/,/, $probe);
	   	my ($LP, $TmLP, $RP, $TmRP, $target);
		if($strand =~ /w/i){
			$LP = $a_oligo;
			$RP = $b_oligo;
			$target = substr($fwdSeq, $a_start, $b_start-$a_start+$b_len);
			$TmLP = $a_Tm;
			$TmRP = $b_Tm;
		}else{
			$LP = $b_oligo;
			$RP = $a_oligo;
			$target = substr($revSeq, $seqLen-$b_start-$b_len, $b_start-$a_start+$b_len);
			$TmLP = $b_Tm;
			$TmRP = $a_Tm;
		}
		my $targetStart = $a_start+$a_len;
		my $targetEnd = $b_start;
		push(@probeInfo,  "$targetStart\t$targetEnd\t$LP\t$TmLP\t$RP\t$TmRP\t$score\t$strand\t$target");
   	}
   	@probeInfo = Sort_Table(
                            	cols      => '9',
                                field     => '1',
                                sorting   => 'ascending',
                                structure => 'csv',
                                separator => '\t',
                                data      => \@probeInfo);	
	return @probeInfo;
}

# Calculate the coverage for a set of probes.
sub calcSetCoverage(){
	my $h_probeInfo = shift;
	my $regionStart = shift;
	my $regionEnd =shift;
	my @baseCoverred;
	my $totalCoverage = 0;
	foreach my $probe (@{$h_probeInfo}){
	   	my ($a_start, $b_start, $a_seq, $a_Tm, $b_seq, $b_Tm, $score, $strand, $target) = split(/\t/, $probe);
		for(my $i=$a_start; $i<=$b_start; $i++){
				if($i>=$regionStart && $i<=$regionEnd){
					$baseCoverred[$i]=1; 
				}
		}	
	}
	foreach my $baseStat (@baseCoverred){
		$totalCoverage ++ if($baseStat);
	}
	return $totalCoverage;	
}

# Calculate the combined score for a set of probes.
sub calcSetScore(){
	my $h_probeSets = shift;
	my $regionStart = shift;
	my $regionEnd =shift;
	my $regionLen = shift;	
	my @baseCoverred;
	my $thescores;
	my $totalScore = 0;
	my $totalCaptureSeqLen =0;
	
	my $found = 0;
	
	foreach my $probe (@{$h_probeSets}){
		next if(!$probe);
		
	   	my ($a_start, $a_len, $a_Tm, $a_oligo, $b_start, $b_len, $b_Tm, $b_oligo, $score, $strand) = split(/,/, $probe);
	   	next if(!$a_len);
		$totalCaptureSeqLen  += $b_start-$a_start-$a_len;
		for(my $i=$a_start+$a_len; $i<$b_start; $i++){
			next if($i < $regionStart || $i > $regionEnd);
			if($strand =~ /w/i){
				$baseCoverred[$i]++;
			}else{
				$baseCoverred[$regionLen - $i]++;
			}
		}	
		$totalScore += $score;
	}	
	$totalScore/=scalar(@{$h_probeSets});
	my $uniqueBpCoverred=0;
	my $redundantBases=0;
	for(my $i=0; $i<scalar(@baseCoverred); $i++){
		$uniqueBpCoverred++ if($baseCoverred[$i]);
		$redundantBases++ if($baseCoverred[$i] && $baseCoverred[$i]>1);
	}
	my $pctCoverred = $uniqueBpCoverred/($regionEnd - $regionStart);
	
	$pctCoverred = 1 if($pctCoverred>1);
	my $redundancy = $redundantBases/($regionEnd - $regionStart);
	# A good probe set should cover more region with little redundancy.

	$totalScore += $pctCoverred*$pctCoverred*2000;
	$totalScore -= $redundancy*500;
	#Give extra points if a probe set covers >95% of the target region
	$totalScore += 500 if($pctCoverred>0.90);
	return $totalScore;	
}

#rules:
#  1. No overlap between two probes on the same strand
###  2. A probe should cover a gap of at least $minBaseCoverred bases on the opposite strand
sub isCompatibleProbes(){
	my $h_probeSets = shift;
	my $new_probe = shift;	
	my $minBaseCoverred = 0; # 221007: changed from 50 to 0

	return 0 if(!$new_probe);
	my (@fwdBaseCoverred, @revBaseCoverred, @primerBaseCoverred);
	my $found = 0;
	foreach my $probe (@{$h_probeSets}){
		next if(!$probe);
		return 0 if($probe eq $new_probe);
		my ($a_start, $a_len, $a_Tm, $a_oligo, $b_start, $b_len, $b_Tm, $b_oligo, $score, $strand) = split(/,/, $probe);
		next if($strand !~ /[wc]/i || $a_start<0 || $b_start<0);
		for(my $i=0; $i<$a_len; $i++){
				$primerBaseCoverred[$i+$a_start]=1;
		}	
		for(my $i=0; $i<$b_len; $i++){
				$primerBaseCoverred[$i+$b_start]=1;
		}		
	
		if($strand =~ /w/i){
			for(my $i=$a_start+$a_len; $i<$b_start; $i++){
				$fwdBaseCoverred[$i]=1;
			}	
		}else{
			for(my $i=$a_start+$a_len; $i<$b_start; $i++){
				$revBaseCoverred[$i]=1;
			}	
		}
	} 
	
	my ($a2_start, $a2_len, $a2_Tm, $a2_oligo, $b2_start, $b2_len, $b2_Tm, $b2_oligo, $score2, $strand2) = split(/,/, $new_probe);
	return 0 if($strand2 !~ /[wc]/i);
	my $base_coverred=0;
	my $base_overlap=0; # 221007 Kian: counter for the overlapping bases
	for(my $i=0; $i<$a2_len; $i++){
		$base_overlap++ if($primerBaseCoverred[$i+$a2_start]); # 221007 Kian
		# return 0 if($primerBaseCoverred[$i+$a2_start]);
	}	
	for(my $i=0; $i<$b2_len; $i++){
		$base_overlap++ if($primerBaseCoverred[$i+$b2_start]); # 221007 Kian
		# return 0 if($primerBaseCoverred[$i+$b2_start]);
	}		

	if($strand2 =~ /w/i){
		for(my $i=$a2_start+$a2_len; $i<$b2_start; $i++){ 
			$base_overlap++ if($fwdBaseCoverred[$i]); # 221007 Kian
			# return 0 if($fwdBaseCoverred[$i]);
			$base_coverred++ if(!$revBaseCoverred[$i]);
		}		
	}else{
		for(my $i=$a2_start+$a2_len; $i<$b2_start; $i++){
			$base_overlap++ if($revBaseCoverred[$i]); # 221007 Kian
			# return 0 if($revBaseCoverred[$i]);
			$base_coverred++ if(!$fwdBaseCoverred[$i] );
		}  
	}
	return 0 if($base_overlap > $maxOverlap); # 221007 Kian
	return 0 if($base_coverred < $minBaseCoverred);
	return 1;		
}

#For a small exon, only one probe is enough. Therefore the search of candidate primers is simpler and faster.
sub getProbeList_smallExon(){
	my ($template_seq, $strand, $targetStart, $targetEnd, $minCG) = @_;	
	my $seqLen = length($template_seq);	
	my $gapMinLen;
	my $gapMaxLen;
	$gapMinLen = $targetMinLen;
	$gapMaxLen = $targetMaxLen;

	my @oligoListA;
	my @oligoListB;
	my %bestProbes;
	my @probeList;

	for(my $i=0; $i<$targetStart-$primerMinLen; $i++){
		for(my $j=$primerMinLen; $j<=$primerMaxLen; $j++){			
			my $oligo = uc(substr($template_seq, $i, $j));
			last if($oligo =~ /[NnMRWSYKVHDBmrwsykvhdb]/          #no Ns (repeat masked sequences)
					|| $oligo =~ /GGGGGG/
                                        || $oligo =~ /TTTTTT/
					|| $oligo =~ /CCCCCC/
					|| $oligo =~ /AAAAAA/);       #no five A/T/G/Cs in a run
			next if($arewebisulfite == 1 and ($oligo =~ s/CG/CG/g) > $minCG);
			my $Tm = shortOligoTm($oligo, 200, 1.5, 50, 0.2, 0.0, 50);
			next if($Tm < $primerMinTm);
			last if($Tm > $primerMaxTm);
			my $score =  0; #-8*abs($H1OptLen-$j); 
			#nucleotide biases in A and C have stronger effects than G and T.
			my %nuclDistTable = &nuclDist($oligo);
			$score -= abs($nuclDistTable{'A'}-0.25)*200 if($nuclDistTable{'A'});
			$score -= abs($nuclDistTable{'C'}-0.25)*200 if($nuclDistTable{'C'});
			push(@oligoListA, "$i,$j,$score,$Tm,$oligo");
		}
	}
	
	for(my $i=$targetEnd; $i<$seqLen-$primerMinLen; $i++){
		for(my $j=$primerMinLen; $j<=$primerMaxLen; $j++){			
			my $oligo = uc(substr($template_seq, $i, $j));
			last if($oligo =~ /[NnMRWSYKVHDBmrwsykvhdb]/          #no Ns (repeat masked sequences)
					|| $oligo =~ /GGGGGG/
                                        || $oligo =~ /TTTTTT/
                                        || $oligo =~ /CCCCCC/
                                        || $oligo =~ /AAAAAA/);       #no five A/T/G/Cs in a run
			next if($arewebisulfite == 1 and ($oligo =~ s/CG/CG/g) > $minCG);
			my $Tm = shortOligoTm($oligo, 200, 1.5, 50, 0.2, 0.0, 50);
			next if($Tm < $primerMinTm);
			last if($Tm > $primerMaxTm);
			my $score =  0; #-8*abs($H2OptLen-$j); 
			#nucleotide biases in A and C have stronger effects than G and T.
			my %nuclDistTable = &nuclDist($oligo);
			$score -= abs($nuclDistTable{'A'}-0.25)*200 if($nuclDistTable{'A'});
			$score -= abs($nuclDistTable{'C'}-0.25)*200 if($nuclDistTable{'C'});			

			push(@oligoListB, "$i,$j,$score,$Tm,$oligo");
		}
	}

	#Optional folding energy calculation
	my $target_folden;
	if ($using_unafold)
	{
		$target_folden = calc_folding_energy_unafold($template_seq);
	}
	else
	{
		$target_folden = -5;
	}

	foreach my $candidateA (@oligoListA){
		my ($a_start, $a_len, $a_score, $a_Tm, $a_oligo) = split(/,/, $candidateA);		
		#next if ($a_Tm < 55 || ($a_start+$gapMinLen+$a_len) >= $seqLen); # to speed up the search a little bit.
		next if (($a_start+$gapMinLen+$a_len) >= $seqLen); # to speed up the search a little bit.
		foreach my $candidateB (@oligoListB){
			my ($b_start, $b_len, $b_score, $b_Tm, $b_oligo) = split(/,/, $candidateB);
			my $gap = $b_start-$a_start-$a_len;
			
			#Constraints like this can make experiments easier, but reduce number of probes found.
			#next if(! ($a_len + $b_len == 40) );
			#next if(! (length($a_oligo)+length($b_oligo) == $H1_plus_H2_Len) );
			next if( length($a_oligo)+length($b_oligo) > $H1_plus_H2_Len );
			#next if(($a_len + $gap + $b_len) != 235);

			
			next if($gap < $gapMinLen 
					|| $gap>$gapMaxLen
					|| $a_oligo =~ /AAAAAA$/
					|| $a_oligo =~ /[NMRWSYKVHDB]/
					|| $a_oligo =~ /TTTTTT$/
					|| $b_oligo =~ /^AAAAAA/
					|| $b_oligo =~ /^TTTTTT/
					|| $b_oligo =~ /[NMRWSYKVHDB]/);
			
			#Calculate score
			my $p_score = 0;
			$p_score = calc_probe_score($target_folden, $template_seq, $a_start, $a_start+$a_len, $b_start, $b_start+$b_len,0); #last argument is zero: no bisulfite conversion was done.

			my $candidate = ($strand =~ /w/i) ? 
			      sprintf("%d,%d,%d,%s,%d,%d,%d,%s,%f,W", $a_start, $a_len, $a_Tm, $a_oligo, $b_start, $b_len, $b_Tm, $b_oligo, $p_score) 
			     :sprintf("%d,%d,%d,%s,%d,%d,%d,%s,%f,C", $seqLen-$b_start-$b_len, $b_len, $b_Tm, $b_oligo, $seqLen-$a_start-$a_len, $a_len, $a_Tm, $a_oligo, $p_score);
			my $bin = int(($a_start+$b_start+$b_len)/10);
			if(!$bestProbes{$bin} || $bestProbes{$bin}->{'s'} < $p_score){
				$bestProbes{$bin}->{'p'}=$candidate;
				$bestProbes{$bin}->{'s'}=$p_score;
			}
		}
	}
	my @array_of_keys = keys(%bestProbes);
	foreach my $bin (keys(%bestProbes)){
		push(@probeList, $bestProbes{$bin}->{'p'});
	}
	my @sorted_probeList = Sort_Table(
                            	cols      => '10',
                                field     => '9',
                                sorting   => 'descending',
                                structure => 'csv',
                                separator => ',',
                                data      => \@probeList);
	return $sorted_probeList[0];
}

# Generate a list of candidate probes which contain two capturing arms of 
#   appropriate Tm within a appropriate distance
sub getProbeList(){
	my $template_seq = shift;
	my $strand = shift;
	my $minCG = shift;
	my $seqLen = length($template_seq);	
	my $gapMinLen = $targetMinLen;
	my $gapMaxLen = $targetMaxLen;

	my @oligoList;
	my %bestProbes;

	my @oligoStarts;
	for(my $i=0; $i<$seqLen; $i++){
		for(my $j=$primerMinLen; $j<=$primerMaxLen; $j++){
			next if($i+$j>=$seqLen);
			my $oligo = uc(substr($template_seq, $i, $j));
			last if($oligo =~ /[NnMRWSYKVHDBmrwsykvhdb]/          #no Ns (repeat masked sequences)
				|| $oligo =~ /GGGGGG/
                                        || $oligo =~ /TTTTTT/
                                        || $oligo =~ /CCCCCC/
                                        || $oligo =~ /AAAAAA/);       #no five A/T/G/Cs in a run
			next if($arewebisulfite == 1 and ($oligo =~ s/CG/CG/g) > $minCG);
			my $Tm = shortOligoTm($oligo, 200, 1.5, 50, 0.2, 0.0, 50);
			next if($Tm < $primerMinTm);
			last if($Tm > $primerMaxTm);
			my $score = 0;
			
			#nucleotide biases in A and C have stronger effects than G and T.
			my %nuclDistTable = &nuclDist($oligo);
			$score -= abs($nuclDistTable{'A'}-0.25)*200 if($nuclDistTable{'A'});
			$score -= abs($nuclDistTable{'C'}-0.25)*200 if($nuclDistTable{'C'});			
			push(@oligoList, "$i,$j,$score,$Tm,$oligo");
			push(@oligoStarts,sprintf("%05d",$i));
		}
	}
	
	@oligoList = Sort_Table(
                            	cols      => '5',
                                field     => '1',
                                sorting   => 'ascending',
                                structure => 'csv',
                                separator => ',',
                                data      => \@oligoList);
	@oligoStarts = sort(@oligoStarts);
	my $totalOligos = scalar(@oligoList);
	my @probeList;
	
	#Optional folding energy calculation
	my $target_folden;
	if ($using_unafold)
	{
		$target_folden = calc_folding_energy_unafold($template_seq);
	}
	else
	{
		$target_folden = -5;
	}
		
	for(my $i=0; $i<$totalOligos-1; $i++){
		last if($oligoStarts[$i]+$gapMinLen >$seqLen);
		my ($a_start, $a_len, $a_score, $a_Tm, $a_oligo) = split(/,/, $oligoList[$i]);		
		#next if ($a_Tm <55); # speeds up the search a little bit, but less probes found.
		for(my $j=$i+1; $j<scalar(@oligoList); $j++){
			next if($oligoStarts[$j]<$a_start+$gapMinLen+$a_len);
			last if($oligoStarts[$j]>$a_start+$gapMaxLen+$a_len);
			my ($b_start, $b_len, $b_score, $b_Tm, $b_oligo) = split(/,/, $oligoList[$j]);
			my $gap = $b_start-$a_start-$a_len;			
			
			#next if(($a_len + $b_len) > 51);
			#next if(($a_len + $gap + $b_len) != 235);
			#next if(! (length($a_oligo)+length($b_oligo) == $H1_plus_H2_Len) );
			next if( length($a_oligo)+length($b_oligo) > $H1_plus_H2_Len );
			
			next if($gap < $gapMinLen 
					|| $gap>$gapMaxLen
					|| $a_oligo =~ /AAAAAA$/
					|| $a_oligo =~ /TTTTTT$/
					|| $b_oligo =~ /^AAAAAA/
					|| $b_oligo =~ /^TTTTTT/);
		#Calculate score		
		my $p_score = 0;
		$p_score = calc_probe_score($target_folden, $template_seq, $a_start, $a_start+$a_len, $b_start, $b_start+$b_len,0); #last argument is zero: no bisulfite conversion was done.

			my $candidate = ($strand =~ /w/i) ? 
			      sprintf("%d,%d,%d,%s,%d,%d,%d,%s,%f,W", $a_start, $a_len, $a_Tm, $a_oligo, $b_start, $b_len, $b_Tm, $b_oligo, $p_score) 
			     :sprintf("%d,%d,%d,%s,%d,%d,%d,%s,%f,C", $seqLen-$b_start-$b_len, $b_len, $b_Tm, $b_oligo, $seqLen-$a_start-$a_len, $a_len, $a_Tm, $a_oligo, $p_score);
			
			my $bin = int(($a_start+$b_start+$b_len)/50)+1;
			if(!$bestProbes{$bin} || $bestProbes{$bin}->{'s'} < $p_score){
				$bestProbes{$bin}->{'p'}=$candidate;
				$bestProbes{$bin}->{'s'}=$p_score;
			}
		}
	}
	foreach my $bin (keys(%bestProbes)){
		push(@probeList, $bestProbes{$bin}->{'p'});
	}
	my @sorted_probeList = Sort_Table(
                            	cols      => '10',
                                field     => '9',
                                sorting   => 'descending',
                                structure => 'csv',
                                separator => ',',
                                data      => \@probeList);
	return @sorted_probeList;
}

sub RevComp(){
    my $seq = shift;
	my $seqLen = length($seq);
	my $revcom = '';
	for(my $i = 0; $i < $seqLen; $i++){
		#For any base that is not A/T/G/C, such as a degenerate base, use "N".
		my $rcBase = $cmpTable{substr($seq, $i, 1)} ? $cmpTable{substr($seq, $i, 1)} :'N';
		$revcom = $rcBase . $revcom;
	}
    return $revcom;
}

sub bisulfiteCvt()
{
	my $input_sequence = shift;
	
	#Protect CG dinucleotides...however, LATER in probe->oligo conversion, remember to make degenerate probes
	$input_sequence =~ s/CG/QX/g; #Convert to a random impossible string so we don't have any errors
	
	#Bisulfite Conversion
	$input_sequence =~ tr/C/T/;
	
	#Restore protected dinucleotides
	$input_sequence =~ s/QX/CG/g;

	return $input_sequence;
}


sub calc_probe_score()
{

	my $target_folden = shift;
	my $total_seq = shift;
	my $a_start = shift;
	my $a_end = shift;
	my $b_start = shift;
	my $b_end = shift;
	my $bisulfite = shift;
	
	my $H1_oligo = uc(substr($total_seq,$a_start,$a_end-$a_start));
	my $H2_oligo = uc(substr($total_seq,$b_start,$b_end-$b_start));
	my $target_oligo = uc(substr($total_seq,$a_end,$b_start-$a_end));
	
	# Length of Target
	my $target_length = length($target_oligo);
	
	# GC Content of Target
	my $target_gc = ((uc($target_oligo) =~ tr/[GC]//) + (uc($target_oligo) =~ tr/N//)/2)/ $target_length; #Count half the N's.
	
	# Melting temperature of H1
	my $H1_Tm = shortOligoTm($H1_oligo, 200, 1.5, 50, 0.2, 0.0, 50);
	# Melting temperature of H2
	my $H2_Tm = shortOligoTm($H2_oligo, 200, 1.5, 50, 0.2, 0.0, 50);
	
	# Length of H1
	my $H1_length = length($H1_oligo);
	# Length of H2
	my $H2_length = length($H2_oligo);
	
	# Proximal Ligation Base
	my $proxligdinc = substr($H1_oligo,-1,1); # Xuwen changed this from -2,2 to -1,1 for splintR PLP design, 3' end for splintR
	# $proxligdinc =~ s/C/T/g;
	
	# Proximal Extension Base
	my $proxextdinc = substr($H2_oligo,0,1); # Xuwen changed this from 2 to 1 for splintR PLP design, , 5' end for splintR
	# $proxextdinc =~ s/C/T/g;
	
	my $proxlig53 = $proxextdinc.$proxligdinc;
	my @liglisttoavoid = ('GC','GG','GT','CA','CC','CG','CT');
	# H1 GC Content
	my $H1_gc = ((uc($H1_oligo) =~ tr/[GC]//) + (uc($H1_oligo) =~ tr/N//)/2)/ $H1_length;
	# H2 GC Content
	my $H2_gc = ((uc($H2_oligo) =~ tr/[GC]//) + (uc($H2_oligo) =~ tr/N//)/2)/ $H2_length;

	my ($x1,$x2,$x3,$x4,$x5,$x6,$x7,$x8,$x9,$x10,$x11,$probe_score);
	
	#Now, calculate score with neural network
	if($arewebisulfite == 0)
	{
		$x1 = $target_folden;
		$x2 = $target_length;
		$x3 = $target_gc;
		$x4 = $H1_Tm;
		$x5 = $H2_Tm;
		$x6 = $H1_length;
		$x7 = $H2_length;
		if (grep { $_ eq $proxlig53 } @liglisttoavoid) { # Xuwen added this rule for splintR PLP design
			# print 'The string is in the list.';
    		$probe_score = neural_net_capture($x1,$x2,$x3,$x4,$x5,$x6,$x7)/2;	
			} else {
    			# print 'The string is not in the list.';
				$probe_score = neural_net_capture($x1,$x2,$x3,$x4,$x5,$x6,$x7);	
				}
		# $probe_score = neural_net_capture($x1,$x2,$x3,$x4,$x5,$x6,$x7);	
	}
	else
	{
		$x1 = $H1_length;
		$x2 = $H2_length;
		$x3 = $H1_Tm;
		$x4 = $H2_Tm;
		$x5 = $target_length;
		$x6 = $target_gc;
		$x7 = $proxLigDincTable{$proxligdinc};
		$x8 = $proxExtDincTable{$proxextdinc};
		$x9 = $target_folden;

		$probe_score = neural_net_capture_LC($x1,$x2,$x3,$x4,$x5,$x6,$x7,$x8,$x9);
		$probe_score = $probe_score - 0.01*(uc($H1_oligo.'NN'.$H2_oligo) =~ tr/C//);
	}

	return $probe_score;
}

sub calc_folding_energy_unafold()
{
	$target_oligo = shift();
	
	#Output oligo to UNAFOLD and run hybrid-ss-min to find folding energy
	my $random_number = int(rand(100000));
	open(UNAFILE,">temp_unafold.$random_number.seq");
	print UNAFILE "$target_oligo;\n";
	close(UNAFILE);
	`/usr/local/bin/hybrid-ss-min --NA=DNA --energyOnly --tmin=37 --tmax=37 temp_unafold.$random_number.seq`;
	
	#Read folding energy from file
	open(UNAOUTPUT,"temp_unafold.$random_number.dG");
	#Skip over first line (column headers)
	my $line = <UNAOUTPUT>;
	$line = <UNAOUTPUT>;
	my @fields = split(/\t/,$line);
	my $foldingenergy = $fields[1];
	close(UNAOUTPUT);

	unlink("temp_unafold.$random_number.seq");
	unlink("temp_unafold.$random_number.run");
	unlink("temp_unafold.$random_number.dG");
	
	#Return folding energy
	return $foldingenergy;
}


sub shortOligoTm(){
        my $seq = shift;
        my $C_primer = shift;           # nM
        my $C_Mg = shift;               # mM
        my $C_MonovalentIon = shift;    #mM
        my $C_dNTP = shift;             #mM
        my $percentage_DMSO = shift;
        my $percentage_annealed = shift; #percentage of templates that anneal to primers

        #$seq =~ s/[ \t\n]+//g;
	$seq = uc($seq);
        $percentage_annealed = 50.0 if (!$percentage_annealed);
        $percentage_annealed /= 100.0;

        my $C_SodiumEquivalent = $C_MonovalentIon + 120 * sqrt($C_Mg-$C_dNTP);
        my $seqLength = length($seq);
        my $dH = $deltaH{'pm'}->{substr($seq, 0, 1)} + $deltaH{'pm'}->{substr($seq, $seqLength-1, 1)};
        my $dS = $deltaS{'pm'}->{substr($seq, 0, 1)} + $deltaS{'pm'}->{substr($seq, $seqLength-1, 1)};
		
	my ($internal_dH, $internal_dS) = &quickInternaldHdS($seq);
	$dH+= $internal_dH;
	$dS+= $internal_dS;	
	$dS += 0.368 * $seqLength * log($C_SodiumEquivalent/1000.0);
        my $Tm = sprintf("%5.2f", ($dH * 1000) / ($dS + $R * (log($C_primer*(1-$percentage_annealed)/$percentage_annealed)-21.4164)) - 273.15 - 0.75*$percentage_DMSO);
        return $Tm;
}

# use two hash table for dH and dS to speed up the calculation.
sub quickInternaldHdS(){
	my $seq = shift;
	my ($dH, $dS);
	my $seqLength = length($seq);
	my $kMerLen = 8;
	for(my $i = 0; $i < $seqLength-1; $i +=$kMerLen){
		my $subSeq = substr($seq,$i,$kMerLen);
		my $subSeqLen = length($subSeq);
		if(!$dH_buffer{$subSeq}){
		for(my $i = 0; $i < $subSeqLen-1; $i ++){
                	$dH_buffer{$subSeq} += $deltaH{'pm'}->{substr($subSeq, $i, 2)};
	                $dS_buffer{$subSeq} += $deltaS{'pm'}->{substr($subSeq, $i, 2)};
		}
			$dH_buffer{&RevComp($subSeq)}=$dH_buffer{$subSeq};
			$dS_buffer{&RevComp($subSeq)}=$dS_buffer{$subSeq};
		}
		$dH += $dH_buffer{$subSeq};
		$dS += $dS_buffer{$subSeq};
		if($i+$kMerLen < $seqLength){
			$dH += $deltaH{'pm'}->{substr($seq, $i+$kMerLen-1, 2)};
			$dS += $deltaS{'pm'}->{substr($seq, $i+$kMerLen-1, 2)};
		}
	}
	return($dH, $dS);
}

sub nuclDist(){
	my $seq = shift;
	my %distTable;
	my $seqLen = length($seq);
	while(length($seq)>0){
		my $base = chop($seq);
		$distTable{$base}+=1/$seqLen;
	}
	return %distTable;
}

sub flipDNA()
{
	my $oligo = shift;
	
	my $UColigo = uc($oligo);
	$UColigo =~ tr/G/c/;
	$UColigo =~ tr/C/g/;
	$UColigo =~ tr/A/t/;
	$UColigo =~ tr/T/a/;
	my $newoligo = uc($UColigo);
	
	return $newoligo;
}

sub convertBaseToNumATCG()
{
	my $base = shift;
	my $outnum;
	
	if ($base eq 'A')
	{
		$outnum = 0;
	}
	if ($base eq 'T')
	{
		$outnum = 1;
	}
	if ($base eq 'C')
	{
		$outnum = 2;
	}
	if ($base eq 'G')
	{
		$outnum = 3;
	}
	
	return $outnum;
}

1;



1;
