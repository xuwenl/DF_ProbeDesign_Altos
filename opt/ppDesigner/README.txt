Padlock Probe Designer (ppDesigner)

Installation instructions
-------------------------
ppDesigner runs on any UNIX-like OS.  It has been tested on Ubuntu Linux 9.04 - 10.04 and MacOS X 10.5 to 10.6.  It should work on any modern Mac or Linux-based system.

Prerequisites:
	Perl (included with most UNIX-like systems)
	Ensure that your Perl install has the File::Temp and Sort::Array modules.  Most should.
	BioPerl (download and install from here: http://www.bioperl.org/)

Recommended Optional Components:
	UNAFold (for better probe efficiency prediction; download and install from here: http://dinamelt.bioinfo.rpi.edu/download.php)
	(Using UNAFold will result in a more accurate prediction of probe efficiency.  It is not required, and will not change the probe sequence.)

To install, simply copy the ppDesigner directory to a directory on your machine.  We recommend using /opt/ppDesigner/ as the target.

Make all scripts executable by typing the command:
	chmod +x /opt/ppDesigner/src/*.pl

If you chose to use a different location than /opt/ppDesigner, you will need to open the probe_parameters.pl file and change line 4 to point to the new location of the NN_param.txt file.

If you are using UNAFold, you will need to make sure that the UNAFold executables are in your $PATH variable.


Usage instructions
------------------
The general use case for ppDesigner is as follows:
1) Create a job file and a target file
2) Obtain target genome in FASTA format
3) Run the job file through the probe designer
4) Convert the probe designer output file to a list of synthesizable oligos

* Job and Target File Creation *
To create a job file, we suggest you use the template provided in the Example directory.  Modify the padlock probe parameters as needed.
$softwareDir - location of the software
$HsDir - location of genome
$exon_info_file - This variable points to your target file.
$primerMaxLen - Maximum length of a capturing arm.
$primerMinLen - Minimum length of a capturing arm.
$H1_plus_H2_len - Total maximum length of the two capturing arms together.  If zero, total arm maximum length will be 2*$primerMaxLen.
$targetMinLen - Minimum length of the target.
$targetMaxLen - Maximum length of the target.
$arewebisulfite - Set to 1 if bisulfite conversion is to be performed.  If designing genomic probes, set to 0.
$using_unafold - Set to 1 if UNAFold is being used.  While UNAFold will not change the choice of probes, it will increase the prediction of probe efficiency.  Generally 1.


To create a target file, we suggest you use the template provided in the Example directory.  The target file is a tab-separated text file with the following columns:
TargetID	FASTA_filename	Start	End	Strand(optional)

TargetID should be alphanumeric, with underscores allowed.
FASTA_filename should be the name of the genome file, minus the trailing '.fa' extension.
Start should be the numerical target start position.
End should be the numerical target end position.
Strand is optional and allows strand-specific probe generation; if set to either + or -, only probes from the + or - strand respectively will be returned.  This allows targeting of single-stranded molecules.

* Obtain Target Genome *
Place the genome in FASTA format in the $HsDir directory specified by the job file.  Each FASTA file should contain one DNA molecule.
The file should be named XXXX.fa; XXXX should be listed in the target file's column 2.

* Run the job file through the probe designer *
This step can take a very long time (hours to days, depending on the number of probes and type of targets).

From a terminal, execute the probe designer with the following syntax:
/opt/ppDesigner/src/ppDesigner.pl /opt/ppDesigner/Example/jobFile.pl > /opt/ppDesigner/Example/Results/outputFile.txt
If you would like to only process one genome out of many listed in the target file, run the program as follows:
/opt/ppDesigner/src/ppDesigner.pl /opt/ppDesigner/Example/jobFile.pl NC_011750 > /opt/ppDesigner/Example/Results/NC_011750_outputFile.txt
This will only generate probes for the NC_011750 FASTA file.

* Convert the probe designer output file to a list of synthesizable oligos *
To generate a list of oligos, simply run either long or short padlock generation scripts:

/opt/ppDesigner/src/primer2padlock_short.pl [max_H1_plus_H2_len] ["IDT" or "notIDT"] < /opt/ppDesigner/Example/outputFile.txt > /opt/ppDesigner/Example/outputFile.probes.txt
/opt/ppDesigner/src/primer2padlock_long.pl [max_H1_plus_H2_len] ["IDT" or "notIDT"] < /opt/ppDesigner/Example/outputFile.txt > /opt/ppDesigner/Example/outputFile.probes.txt

For bisulfite probes: note that the H1 and H2 output from the previous step will have many G's and only few C's in CG positions. The actual probes will contain the reverse complement of these sequences, thus, provide correct annealing.

This will generate the file outputFile.probes.txt.  For ordering these oligos (for example from IDT), the suggested oligo name is in column 15 and the oligo sequence is in column 16 (the last two columns).  These oligos can be synthesized and used for capture. 
For "notIDT" synthesis, the suggested probe sequences wil contain adapters on both ends to allow for parallel oligonucleotides synthesis on microarrays, PCR amplification, and probes preparation.

If you have any questions, please contact Kun Zhang (kzhang@bioeng.ucsd.edu) or Athurva Gore (ajgore@eng.ucsd.edu).
