# Spatial ppDesigner
sppDesigner is a software for designing padlock probes for in situ profiling. It is based on ppDesginer from Diep et al. *Nature Methods* (2012). 

## Requirements
- Perl interpreter which is preinstalled on many linux systems.
- [SortArray](https://metacpan.org/pod/Sort::Array) package on perl
```
sudo perl -MCPAN -e shell
install Sort::Array
```
- a reference genome with one fasta file per chromosome
## How to run
To call ppDesigner you need two files 1) a file specifying regions on the reference genome for which you want to target with padlock probes, 2) a "job file" with all the settings. 
We have included an example run with t1.target as the target file, "j1.jfile" as the job file. "run_ppd.sh" runs ppDesigner and creates the output "testrun.ppd". 

The notebook "run_sppd.ipynb" runs sppDesigner for all genes in the genome and creates a resource from which probes can be pulled for every new probe set. To run that, you need a 
genome annotation .gtf file (Gencode or RefSeq), [BWA](https://bio-bwa.sourceforge.net/) aligner, with a BWA genome index. 

# License
sppDesigner is a modified version of ppDesigner (PMID: 22306810). The original software is covered by a copyright owned by UC San Diego (invent@ucsd.edu)
