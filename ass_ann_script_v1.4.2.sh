#!/bin/bash

# exit on the first error encountered
set -e

usage() {
    echo -e "Usage:\n -c threads -s strain_name"
}

# if no args, abort with help
[ $# -eq 0 ] && usage && exit 1

while getopts ":s:c:" arg; do
  case $arg in
    c) # number of threads
      THREADS=${OPTARG}
      ;;
    s) # strain name
      STRAINNAME=${OPTARG}
      ;;
    h) # help
      usage
      exit 0
      ;;
    *) # help then error
      usage
      exit 1
      ;;
  esac
done
if [ -z $STRAINNAME ] ; then
 usage
 exit 1
fi


#check dependencies (makes failing faster)
autoMLST 2>&1 | grep 'ref REFDB'
ALE  > tmp
Bandage --version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/bandage
blastn -version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/blastn
busco --version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/busco
filtlong --version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/filtlong
flye --version | cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/flye2.9
porechop --version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/porechop
prokka --version 2>&1|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/prokka
trim_galore --version | cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/trimgalore
unicycler --version | cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/unicycler
antismash --version |cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/antiSMASH

#clean nanopore data
porechop -i nanopore/ -o nanopore/porechop.fq
filtlong --min_length 1000 nanopore/porechop.fq > nanopore/min1kb.fq

#assembly nanopore data
flye -t $THREADS -i 5 -o flye --plasmids --nano-hq nanopore/min1kb.fq
cat flye/flye.log|grep 'Reads N50.*' -o|cut -f 3 -d ' '|printf 'nanopore N50: %s\n' "$(cat)" > n50
rm nanopore/porechop.fq nanopore/min1kb.fq


#the illumina read naming needs to exactly match the below for trimgalore to recognize it
trim_galore -j 8 --length 100  -o illumina --paired --quality 20 --fastqc --gzip illumina/*1.fq.gz illumina/*2.fq.gz
zcat illumina/*val_1.fq.gz |grep ^@|wc -l|printf 'illumina read pairs: %s\n' "$(cat)" > ill_pairs

#failsafe check that the illumina and nanopore data match
bowtie2-build flye/assembly.fasta flye/assemby.index
bowtie2 -p $THREADS -x flye/assemby.index -1 illumina/*_val_1.fq.gz -2 illumina/*_val_2.fq.gz -S tmp 2>&1|grep overall > overall_ill_on_flye_mapping_percent
rm tmp
cat overall_ill_on_flye_mapping_percent|sed 's/\..*//' > overall_ill_on_flye_mapping_short
line=$(cat overall_ill_on_flye_mapping_short)
target=60
if [ $line -ge $target ]
then
   echo "the percent of illumina reads mapping on the assembly is $line which is nice"
else
   echo "ERROR: mapped reads indicate that there is an issue with the illumina and/or nanopore dataset used" >> $STRAINNAME.AA.log
   echo "ERROR: percent of mapped reads is $line . At least 60, but usually >95 is expected for datasets from the same source" >> $STRAINNAME.AA.log
   echo "ERROR: mapped reads indicate that there is an issue with the illumina and/or nanopore dataset used"
   echo "ERROR: percent of mapped reads is $line . At least 60, but usually >95 is expected for datasets from the same source" 

   exit 1
fi


#polish assembly
unicycler_polish --threads $THREADS -a flye/assembly.fasta -1 illumina/*_val_1.fq.gz -2 illumina/*_val_2.fq.gz
rm illumina/*val_1.fq.gz illumina/*val_2.fq.gz

#phylogenetic placement
mkdir automlst
cp *final_polish.fasta automlst/.
autoMLST --cpu $THREADS automlst/ automlst/
cat automlst/mash_distances.txt | grep -v \# | sort -k5 | tail -n 1 | cut -f 3,5,7,8 | sed 's/_/\n/1' | sed 's/\t/\n/g' | split -l 1 -a 1 -d - automlst/result
cat automlst/result0 > automlst/genus
#paste automlst/result1 automlst/result2 | sed 's/\t/_/' > automlst/species
paste automlst/result0 automlst/result1 automlst/result2 | sed 's/\t/_/g' > automlst/genus_species_ANI

#prepare for annotation
cat *final_polish.fasta | sed "s/contig/$STRAINNAME/" | sed "s/scaffold/${STRAINNAME}_scaf/" > prokkainput.fna

#annotate: note that the 6 actinobacrterial strains as well as PFA should be included
prokka --outdir ${STRAINNAME}_prokka_actinoannotPFAM --prefix $STRAINNAME --genus `cat automlst/genus` --species sp.$STRAINNAME --cdsrnaolap --cpus $THREADS --rnammer --increment 10 --evalue 1e-05 prokkainput.fna
antismash --output-dir ${STRAINNAME}_antiSMASH --cb-general --cb-subclusters --cb-knownclusters -c $THREADS ${STRAINNAME}_prokka_actinoannotPFAM/$STRAINNAME.gbk  --clusterhmmer --cc-mibig --asf --tigr --pfam2go --html-description `cat automlst/genus_species_ANI`


#make short results log of the assembly
busco -i ${STRAINNAME}_prokka_actinoannotPFAM/$STRAINNAME.faa -l actinobacteria_class_odb10 -o busco -m prot -c $THREADS -q 
Bandage image flye/assembly_graph.gfa $STRAINNAME.jpg
cp prokkainput.fna $STRAINNAME.contigs.fasta
cp flye/assembly_graph.gfa $STRAINNAME.graph.gfa
cp flye/flye.log $STRAINNAME.flye.log
cp ${STRAINNAME}_antiSMASH/$STRAINNAME.zip ${STRAINNAME}.antiSMASH.zip
echo $STRAINNAME|sed 's/^/Strain name: /' >> $STRAINNAME.AA.log
cat automlst/result0|printf 'genus: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/result1|printf 'species: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/result2|printf 'estimated ANI: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/genus_species_ANI|printf 'genus_species_ANI: %s\n' "$(cat)" >> $STRAINNAME.AA.log
assembly-stats -s $STRAINNAME.contigs.fasta|grep number|cut -f 3 |printf 'number of contigs: %s\n' "$(cat)" >> $STRAINNAME.AA.log
assembly-stats -s $STRAINNAME.contigs.fasta|grep total_length|cut -f 3|printf 'total assembly length: %s\n' "$(cat)" >> $STRAINNAME.AA.log
assembly-stats -s $STRAINNAME.contigs.fasta|grep longest|cut -f 3|printf 'longest contig: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat flye/assembly_info.txt >> $STRAINNAME.AA.log
#echo n50 >> $STRAINNAME.AA.log
#echo ill_pairs >> $STRAINNAME.AA.log
cat n50 >> $STRAINNAME.AA.log
cat ill_pairs >> $STRAINNAME.AA.log
cat *changes|wc -l|printf 'illumina polished bases: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat overall_ill_on_flye_mapping_percent|printf 'percent illumina reads mapping on flye assembly: %s\n' "$(cat)" >> $STRAINNAME.AA.log
ls *antiSMASH/*region*gbk|wc -l |sed 's/^/BGC Regions: /' >> $STRAINNAME.AA.log
cat ${STRAINNAME}_prokka_actinoannotPFAM/${STRAINNAME}.log |grep unann|tail -n 1|sed 's/.* There are //1' >> $STRAINNAME.AA.log
cat busco/short*|tail -n 7 >> $STRAINNAME.AA.log

#document tool versions
echo $0 >> $STRAINNAME.AA.log
which autoMLST|sed 's/.*autoMLST\//autoMLST v/1'|sed 's/\/.*//' >> $STRAINNAME.AA.log
Bandage --version|printf 'Bandage v%s\n' "$(cat)"|sed 's/Version: //' >> $STRAINNAME.AA.log
bowtie2 --version |grep version >> $STRAINNAME.AA.log
blastn -version|grep -v Pack|sed 's/: / v/' >> $STRAINNAME.AA.log
busco --version|sed 's/ / v/' >> $STRAINNAME.AA.log
filtlong --version >> $STRAINNAME.AA.log
flye --version 2>&1|printf 'flye v%s\n' "$(cat)" >> $STRAINNAME.AA.log
porechop --version|printf 'porechop v%s\n' "$(cat)" >> $STRAINNAME.AA.log
prokka --version 2>&1|grep prokka_|sed 's/ / v/' >> $STRAINNAME.AA.log
trim_galore --version |grep version|sed 's/.*version /v/'|printf 'trimgalore %s\n' "$(cat)" >> $STRAINNAME.AA.log
unicycler --version 2>&1 >> $STRAINNAME.AA.log
antismash --version >> $STRAINNAME.AA.log 

#clean up temporary files
rm -r  automlst busco* flye ${STRAINNAME}_antiSMASH
rm   prokkainput.fna 0* overall_ill_on_flye_mapping_percent overall_ill_on_flye_mapping_short n50 ill_pairs

#make little celebratory statement marking the finishing of the pipeline
echo 'Your assembly and annotation of'
echo $STRAINNAME  
echo 'has successfully finished' 
 
