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

Further, it requires a set of files containing the version numbers for each tool in use, which on the machine NBC-shared can be found in this location    
  /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1

The following tools are being used (again, locations on NBCshared): 
	
	Tool		Version		Location
	antiSMASH	6.1.0-4de2be2	module
	autoMLST	NA		module
	bandage		0.8.1		/home/WIN.DTU.DK/tuspjo/install/Bandage
	blastn		2.10.0+		/opt/blast/2.10/bin/blastn (loaded with a module)
	busco		BUSCO 5.1.2	module
	filtlong	v0.2.0		/home/WIN.DTU.DK/tuspjo/install/Filtlong/bin/filtlong
	flye2.9		2.9-b1768	module
	porechop	0.2.4		module
	prokka		1.14.6		module
	trimgalore	0.6.4_dev	module
	unicycler	v0.4.8-beta	module

Certain versions of some of the tools needs to be loaded : 
	module load antismash/dev porechop flye unicycler autoMLST prokka/1.14.6 busco trimgalore
(ok, apparently only prokka, previously also BUSCO)
