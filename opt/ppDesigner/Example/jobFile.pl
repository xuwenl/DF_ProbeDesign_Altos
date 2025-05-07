#!/usr/bin/perl -w
#Use probe script with this file as ARGV[0].

our $softwareDir='/opt/ppDesigner/src';
our $HsDir='/GenomeDB/HsGenome/hg18/';

our $exon_info_file='/opt/ppDesigner/Example/targetFile_LC.txt';

our $primerMaxLen=25;
our $primerMinLen=15;
our $H1_plus_H2_Len = 40;

our $targetMinLen=178;
our $targetMaxLen=182;

our $arewebisulfite=1;
our $using_unafold=1;


#DO NOT CHANGE THE BELOW LINES.
eval `cat $softwareDir/probe_parameters.pl` or die 'couldnt parse file';
eval `cat $softwareDir/load_inputs.pl` or die 'couldnt parse file';
eval `cat $softwareDir/get_sequence.pl` or die 'couldnt parse file';
eval `cat $softwareDir/get_probes.pl` or die 'couldnt parse file';
eval `cat $softwareDir/neural_net.pl` or die 'couldnt parse file';
eval `cat $softwareDir/output_text.pl` or die 'couldnt parse file';
