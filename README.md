# AAA-Actino-Assembly-and-Annotation
the ass_ann script goes to github!

This script uses some of the most well known tools to assemble and annotate the genomes of actinobacteria from nanopore fastqs and illumina data. Normally, a filtering step before this pipeline would weed out datasets which does not assemble fully/sufficiently.  

Usage:
./aaa_v1.6.4.sh -c threads -s strain_name

So the current script requires a very particular folder structure, which needs to look like this 

	-> STAINNAME
    	-> illumina
        		somethingsomething.1.fq.gz
        		somethingsomething.2.fq.gz
			someotherthing.1.fq.gz
			someotherthing.2.fq.gz
    	-> nanopore
        		file1fastq.gz
        		file2fastq.gz
        		file3fastq.gz
         		...

Further, it requires a set of files containing the version numbers for each tool in use, which is included in the folder "versions"

The following tools are being used (locations on NBCshared): 
	
	Tool		Version		Location
	antiSMASH	6.1.0-4de2be2	module
	autoMLST	NA		module
	bandage		0.8.1		/home/WIN.DTU.DK/tuspjo/install/Bandage
	blastn		2.10.0+		/opt/blast/2.10/bin/blastn (loaded with a module)
	busco		BUSCO 5.1.2	module
	flye2.9		2.9-b1768	module
	prokka		1.14.6		module
	trimgalore	0.6.4_dev	module
	masurca/polca	v4.0.5		module only to be loaded when used bc loads a bunch of stuff which might interfere
	polypolish	v0.5.0		module
Certain versions of some of the tools needs to be loaded : 
	module load antismash/dev flye autoMLST prokka/1.14.6 busco trimgalore polypolish #porechop filtlong unicycler
(ok, apparently only prokka, previously also BUSCO)

Prokka setup

Prokka annotates a lot more genes if given extra databases to the default ones.

Download PFAM HMM and PGAP HMM and put them in your prokka location in db/hmm/ and press them.

Currently, Pfam-A v.32.0 and PGAP HMMs released 2021-11-10 are used (http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/ https://ftp.ncbi.nlm.nih.gov/hmm/current/ 
