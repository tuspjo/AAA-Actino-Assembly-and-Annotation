# AAA-Actino-Assembly-and-Annotation
the ass_ann script goes to github

This script uses some of the most well known tools to assemble and annotate the genomes of actinobacteria from nanopore fastqs and illumina data. Normally, an filtering step before this pipeline would weed out datasets which does not assemble fully/sufficiently.  

So the current script requires a very particular folder structure, which needs to look like this 
-> STAINNAME
    -> illumina
         somethingsomething.1.fq.gz
         somethingsomething.2.fq.gz
    -> nanopore
         fastqfile1
         fastqfile2
         fastqfile3
         ...
(huh, that does not look good) 

Further, it requires a set of files containing the version numbers for each tool in use, which on the machine NBC-shared can be found in this location    
  /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1

The following tools are being used: 
antiSMASH
assembly-stats
autoMLST
bandage
blastn
busco
filtlong
flye2.9
porechop
prokka
trimgalore
unicycler

TBC

