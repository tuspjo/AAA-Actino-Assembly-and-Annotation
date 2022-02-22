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
Bandage --version|cmp - `dirname "$0"`/versions/bandage
blastn -version|cmp - `dirname "$0"`/versions/blastn
busco --version|cmp - `dirname "$0"`/versions/busco
#filtlong --version|cmp - `dirname "$0"`/versions/filtlong
flye --version | cmp - `dirname "$0"`/versions/flye2.9
#porechop --version|cmp - `dirname "$0"`/versions/porechop
prokka --version 2>&1|cmp - `dirname "$0"`/versions/prokka
trim_galore --version | cmp - `dirname "$0"`/versions/trimgalore
#unicycler --version | cmp - `dirname "$0"`/versions/unicycler
antismash --version |cmp - `dirname "$0"`/versions/antiSMASH
polypolish --version |cmp - `dirname "$0"`/versions/polypolish
module load masurca
masurca --version|cmp - `dirname "$0"`/versions/masurca
module unload masurca

#Bandage --version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/bandage
#blastn -version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/blastn
#busco --version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/busco
##filtlong --version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/filtlong
#flye --version | cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/flye2.9
##porechop --version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/porechop
#prokka --version 2>&1|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/prokka
#trim_galore --version | cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/trimgalore
##unicycler --version | cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/unicycler
#antismash --version |cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/antiSMASH
#polypolish --version |cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/polypolish
#module load masurca
#masurca --version|cmp - /home/WIN.DTU.DK/tuspjo/install/ass_ann_script_dependencies/v.1.4.1/masurca
#module unload masurca

##clean nanopore data #no longer necessary after upgrades to guppy. also porechop hardly knows the RBK110-96 barcodes as development stopped years before the kit was released.
#porechop -i nanopore/ -o nanopore/porechop.fq
#filtlong --min_length 1000 nanopore/porechop.fq > nanopore/min1kb.fq

#assembly nanopore data
flye -t $THREADS -i 5 -o flye --nano-raw nanopore/*gz 
cat flye/flye.log|grep 'Reads N50.*' -o|cut -f 3 -d ' '|printf 'nanopore N50: %s\n' "$(cat)" > n50
#rm nanopore/porechop.fq nanopore/min1kb.fq


#the illumina read naming needs to exactly match the below for trimgalore to recognize it
trim_galore -j 8 --length 100  -o illumina --paired --quality 20 --gzip illumina/*1.fq.gz illumina/*2.fq.gz
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


#polypolish
bwa index flye/assembly.fasta
bwa mem -t $THREADS -a flye/assembly.fasta illumina/*val_1.fq.gz > alignments_1.sam
bwa mem -t $THREADS -a flye/assembly.fasta illumina/*val_2.fq.gz > alignments_2.sam
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish flye/assembly.fasta filtered_1.sam filtered_2.sam > polypolished.fasta 2> polypolish.log
rm alignments_1.sam alignments_2.sam filtered_1.sam filtered_2.sam 

#POLCA
echo blaa
module load masurca
polca.sh -a polypolished.fasta -r 'illumina/*_val_1.fq.gz illumina/*_val_2.fq.gz' -t $THREADS
module unload masurca

#phylogenetic placement
mkdir automlst
cp polypolished.fasta.PolcaCorrected.fa automlst/.
autoMLST --cpu $THREADS automlst/ automlst/
cat automlst/mash_distances.txt | grep -v \# | sort -k5 | tail -n 1 | cut -f 2,3,5,7,8 | sed 's/_/\n/2' | sed 's/\t/\n/g' | split -l 1 -a 1 -d - automlst/result
cat automlst/result1 > automlst/genus
paste automlst/result1 automlst/result2 | sed 's/\t/_/' > automlst/species
paste automlst/result0 automlst/result1 automlst/result2 | sed 's/\t/_/g' > automlst/genus_species_ANI
paste automlst/result1 automlst/result2 automlst/result0 automlst/result3| sed 's/\t/_/g' > automlst/genus_species_ref_ANI

#prepare for annotation
cat polypolished.fasta.PolcaCorrected.fa | sed "s/contig/$STRAINNAME/" | sed "s/scaffold/${STRAINNAME}_scaf/" |sed 's/_polypolish//' > $STRAINNAME.contigs.fasta

#annotate: note that the 6 actinobacrterial strains as well as PFA should be included
cat `dirname "$0"`/trusted_annotations/*gbff >> trusted.gbff
prokka --outdir ${STRAINNAME}_prokka_actinoannotPFAM --prefix $STRAINNAME --genus `cat automlst/genus` --species sp.$STRAINNAME --cdsrnaolap --cpus $THREADS --rnammer --increment 10 --proteins trusted.gbff --evalue 1e-05 $STRAINNAME.contigs.fasta
rm trusted.gbff
antismash --output-dir ${STRAINNAME}_antiSMASH --cb-general --cb-subclusters --cb-knownclusters -c $THREADS ${STRAINNAME}_prokka_actinoannotPFAM/$STRAINNAME.gbk  --clusterhmmer --cc-mibig --asf --tigr --pfam2go --html-description `cat automlst/genus_species_ref_ANI`

#make short results log of the assembly
busco -i ${STRAINNAME}_prokka_actinoannotPFAM/$STRAINNAME.faa -l actinobacteria_class_odb10 -o busco -m prot -c $THREADS -q 
Bandage image flye/assembly_graph.gfa $STRAINNAME.jpg
cp flye/assembly_graph.gfa $STRAINNAME.graph.gfa
cp flye/flye.log $STRAINNAME.flye.log
cp ${STRAINNAME}_antiSMASH/$STRAINNAME.zip ${STRAINNAME}.antiSMASH.zip
echo $STRAINNAME|sed 's/^/Strain name: /' >> $STRAINNAME.AA.log
cat automlst/result1|printf 'genus: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/result2|printf 'species: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/result0|printf 'reference acc: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/result3|printf 'estimated ANI: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/genus_species_ref_ANI|printf 'genus_species_ANI: %s\n' "$(cat)" >> $STRAINNAME.AA.log
assembly-stats -s $STRAINNAME.contigs.fasta|grep number|cut -f 3 |printf 'number of contigs: %s\n' "$(cat)" >> $STRAINNAME.AA.log
assembly-stats -s $STRAINNAME.contigs.fasta|grep total_length|cut -f 3|printf 'total assembly length: %s\n' "$(cat)" >> $STRAINNAME.AA.log
assembly-stats -s $STRAINNAME.contigs.fasta|grep longest|cut -f 3|printf 'longest contig: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat flye/assembly_info.txt >> $STRAINNAME.AA.log
cat n50 >> $STRAINNAME.AA.log
cat ill_pairs >> $STRAINNAME.AA.log
cat ${STRAINNAME}.graph.gfa |grep ^S|cut -f 3|awk '{ print length }'|sed 's/$/ nt/' > ${STRAINNAME}.edgelength
cat ${STRAINNAME}.graph.gfa |grep ^S|cut -f 2,4|sed 's/\t/\tcov /'|paste - ${STRAINNAME}.edgelength >> $STRAINNAME.AA.log
rm ${STRAINNAME}.edgelength
cat ${STRAINNAME}.graph.gfa |grep -P '^L|^P' >> $STRAINNAME.AA.log
cat polypolish.log |grep -P 'Polishing contig' -A 4 --no-group-separator >> $STRAINNAME.AA.log
cat *report >>$STRAINNAME.AA.log
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
#filtlong --version >> $STRAINNAME.AA.log
flye --version 2>&1|printf 'flye v%s\n' "$(cat)" >> $STRAINNAME.AA.log

module load masurca
masurca --version |printf 'masurca/polca v%s\n' "$(cat)" >> $STRAINNAME.AA.log
module unload masurca

polypolish --version 2>&1|printf '%s\n' "$(cat)" >> $STRAINNAME.AA.log
#porechop --version|printf 'porechop v%s\n' "$(cat)" >> $STRAINNAME.AA.log
prokka --version 2>&1|grep prokka_|sed 's/ / v/' >> $STRAINNAME.AA.log
trim_galore --version |grep version|sed 's/.*version /v/'|printf 'trimgalore %s\n' "$(cat)" >> $STRAINNAME.AA.log
#unicycler --version 2>&1 >> $STRAINNAME.AA.log
antismash --version >> $STRAINNAME.AA.log

#clean up temporary files
rm -r  automlst busco* flye ${STRAINNAME}_antiSMASH
rm overall_ill_on_flye_mapping_percent overall_ill_on_flye_mapping_short n50 ill_pairs
rm illumina/*val_1.fq.gz illumina/*val_2.fq.gz
rm polypolished.fasta bwa.err polypolished.fasta.alignSorted.bam polypolished.fasta.alignSorted.bam.bai polypolished.fasta.batches polypolished.fasta.bwa.amb polypolished.fasta.bwa.ann polypolished.fasta.bwa.bwt polypolished.fasta.bwa.pac polypolished.fasta.bwa.sa polypolished.fasta.fai polypolished.fasta.fix.success polypolished.fasta.index.success polypolished.fasta.map.success polypolished.fasta.names polypolished.fasta.PolcaCorrected.fa polypolished.fasta.report polypolished.fasta.report.success polypolished.fasta.sort.success polypolished.fasta.unSorted.sam polypolished.fasta.vcf polypolished.fasta.vc.success polypolish.log samtools.err

#make little celebratory statement marking the finishing of the pipeline
echo 'Your assembly and annotation of'
echo $STRAINNAME
echo 'has successfully finished'

