#!/bin/bash

# exit on the first error encountered
set -e

usage() {
    echo -e "Usage:\n -c threads -s strain_name"
}

# if no args, abort with help
[ $# -eq 0 ] && usage && exit 1

while getopts ":s:c:f:" arg; do
  case $arg in
    c) # number of threads
      THREADS=${OPTARG}
      ;;
    s) # strain name
      STRAINNAME=${OPTARG}
      ;;
    f) # filtlong percentage
      FILT=${OPTARG}
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

# Check for strain_name
if [ -z $STRAINNAME ] ; then
 usage
 exit 1
fi

# check folder structure
if [ ! -d "nanopore" ]
then
	echo "Error: No Nanopore folder found/check directory"
	exit
fi

# check dependencies (makes failing faster)
(autoMLST 2>&1 | grep 'ref REFDB') || (echo "correct version of autoMLST not found/loaded"; exit 1)
(Bandage --version|cmp - `dirname "$0"`/versions/bandage) || (echo "correct version of Bandage not found/loaded"; exit 1)
(blastn -version|cmp - `dirname "$0"`/versions/blastn) || (echo "correct version blastn not found/loaded"; exit 1)
(busco --version|cmp - `dirname "$0"`/versions/busco) || (echo "correct version of Busco not found/loaded"; exit 1)
#filtlong --version|cmp - `dirname "$0"`/versions/filtlong 
(flye --version | cmp - `dirname "$0"`/versions/flye2.9) || (echo "correct version of flye not found/loaded"; exit 1)
(prokka --version 2>&1|cmp - `dirname "$0"`/versions/prokka) || (echo "correct version of Prokka not found/loaded"; exit 1)
(antismash --version |cmp - `dirname "$0"`/versions/antiSMASH) ||(echo "correct version of antiSMASH not found/loaded"; exit 1)
(medaka --version | cmp - `dirname "$0"`/versions/medaka) ||(echo "correct version of medaka not found/loaded"; exit 1)

# assemble nanopore data
zcat nanopore/*gz|gzip > allnp.fq.gz

if [ ! -z $FILT ]
then
	filtlong -p $FILT allnp.fq.gz|gzip > filt.fq.gz
	echo "filtlong -p $FILT was used to filter the raw nanopore data" >>$STRAINNAME.AA.log
	flye -t $THREADS -i 5 -o flye --nano-hq  filt.fq.gz
else
	flye -t $THREADS -i 5 -o flye --nano-hq allnp.fq.gz

fi

##npgm-contigger
#python3 /bigdata/dalofa/installs/npgm-contigger-main/contigger/contigger.py --infile flye/assembly_graph.gfa  --output npgm-contigged.fa 2>npgm-contigger.err
#flye -t $THREADS --polish-target npgm-contigged.fa -o polish -i 5 --nano-raw allnp.fq.gz
#cat polish/polished_1.fasta|awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > polish.singleline.fa #from user ljq on https://www.biostars.org/p/9262/
#echo "the npgm-contigger was used to complete the repeat graph edges into replicons" >> $STRAINNAME.AA.log
#sed -i 's/^contig_/unused_flye_contig_/' flye/assembly_info.txt
#cat polish.singleline.fa |awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > $STRAINNAME.lensort1.fa
#cat $STRAINNAME.lensort1.fa|tr '\n' ' '|sed 's/$/\n/'|sed 's/>/\n>/g'|sed '/^$/d' > $STRAINNAME.lensort2.fa
#cat $STRAINNAME.lensort2.fa|awk '{ print length }'|paste - $STRAINNAME.lensort2.fa |sort -n -r |cut -f 2|sed 's/ /\n/'|sed 's/ //g'|sed '/^$/d' > $STRAINNAME.lensort3.fa
#rm $STRAINNAME.lensort1.fa $STRAINNAME.lensort2.fa

#cat $STRAINNAME.lensort3.fa |sed '/^$/d'|sed 'N;s/\n/ /'|cat -n - |sed 's/^     //' |sed 's/\t>/ /'|sed 's/^/>contig_/'|sed 's/ /\n/2' > oneline
#cat oneline|grep -v \> |awk '{ print length }' |sed 's/^/length /'|sed 's/$/ nt/'> seqlen
#cat oneline| grep \>|paste - seqlen > npgm-stats.txt
#cat oneline|sed 's/ /\n/2' > npgm-contigger.fa
#cat polish/flye.log >> flye/flye.log
#rm allnp.fq.gz 
#rm -r oneline seqlen npgm-contigged.fa polish polish.singleline.fa $STRAINNAME.lensort3.fa

#cat npgm-contigger.fa|sed 's/ .*//' > flye/assembly.fasta

cat flye/flye.log|grep 'Reads N50.*' -o|cut -f 3 -d ' '|printf 'nanopore N50: %s\n' "$(cat)" > n50
cat flye/assembly.fasta |awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > $STRAINNAME.lensort1.fa
cat $STRAINNAME.lensort1.fa|tr '\n' ' '|sed 's/$/\n/'|sed 's/>/\n>/g'|sed '/^$/d' > $STRAINNAME.lensort2.fa
cat $STRAINNAME.lensort2.fa|awk '{ print length }'|paste - $STRAINNAME.lensort2.fa |sort -n -r |cut -f 2|sed 's/ /\n/'|sed 's/ //g' > flye/for_polishing.fa
rm $STRAINNAME.lensort1.fa $STRAINNAME.lensort2.fa

###Polish with medaka
medaka_consensus  -d flye/for_polishing.fa -i allnp.fq.gz -o medaka -m r1041_e82_260bps_sup_v4.0.0 -f -t $THREADS


#phylogenetic placement
mkdir automlst
cp medaka/consensus.fasta automlst/.
autoMLST --cpu $THREADS automlst/ automlst/
cat automlst/mash_distances.txt | grep -v \# | sort -k5 | tail -n 1 | cut -f 2,3,5,7,8 | sed 's/_/\n/2' | sed 's/\t/\n/g' | split -l 1 -a 1 -d - automlst/result
cat automlst/result1 > automlst/genus
paste automlst/result1 automlst/result2 | sed 's/\t/_/' > automlst/species
paste automlst/result0 automlst/result1 automlst/result2 | sed 's/\t/_/g' > automlst/genus_species_ANI
paste automlst/result1 automlst/result0 automlst/result3| sed 's/\t/_/g' > automlst/genus_ref_ANI
paste automlst/result1 automlst/result2 automlst/result0 automlst/result3| sed 's/\t/_/g' > automlst/genus_species_ref_ANI

#prepare for annotation
#cat $STRAINNAME.lensort3.fa | sed "s/contig/$STRAINNAME/" | sed "s/scaffold/${STRAINNAME}_scaf/" |sed 's/_polypolish//' > $STRAINNAME.contigs.fasta
#rm $STRAINNAME.lensort3.fa
cat medaka/consensus.fasta | sed "s/contig/$STRAINNAME/" | sed "s/scaffold/${STRAINNAME}_scaf/" |sed 's/_polypolish//' > $STRAINNAME.contigs.fasta
#rm $STRAINNAME.lensort3.fa 


#annotate: note that the 6 actinobacrterial strains as well as PFA should be included
cat `dirname "$0"`/trusted_annotations/*gbff >> trusted.gbff
prokka --outdir ${STRAINNAME}_prokka_actinoannotPFAM --prefix $STRAINNAME --genus `cat automlst/genus` --species sp. --strain $STRAINNAME --cdsrnaolap --cpus $THREADS --rnammer --increment 10 --proteins trusted.gbff --evalue 1e-05 $STRAINNAME.contigs.fasta
rm trusted.gbff
basename $0|sed 's/^/tool_version:/' > aaaversion
echo $STRAINNAME > strainname
cat automlst/genus|sed 's/^/autoMLST_genus:/' > genus
paste  strainname genus aaaversion |sed 's/\t/_/g'> antiSMASH_html_comment
antismash --output-dir ${STRAINNAME}_antiSMASH --cb-general --cb-subclusters --cb-knownclusters -c $THREADS ${STRAINNAME}_prokka_actinoannotPFAM/$STRAINNAME.gbk  --genefinding-tool none --clusterhmmer --cc-mibig --asf --tigr --pfam2go --html-description `cat antiSMASH_html_comment` 2>antiSMASH.errorlog
rm aaaversion strainname genus antiSMASH_html_comment 

#make short results log of the assembly
busco -i ${STRAINNAME}_prokka_actinoannotPFAM/$STRAINNAME.faa -l actinobacteria_class_odb10 -o busco -m prot -c $THREADS -q 
Bandage image flye/assembly_graph.gfa $STRAINNAME.jpg --names
cp flye/assembly_graph.gfa $STRAINNAME.graph.gfa
cp flye/flye.log $STRAINNAME.flye.log
cp ${STRAINNAME}_antiSMASH/$STRAINNAME.zip ${STRAINNAME}.zip
echo $STRAINNAME|sed 's/^/Strain name: /' >> $STRAINNAME.AA.log
cat automlst/result1|printf 'autoMLST genus: %s\n' "$(cat)" >> $STRAINNAME.AA.log
#cat automlst/result2|printf 'autoMLST species: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/result0|printf 'autoMLST reference acc: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/result3|printf 'autoMLST estimated ANI: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat automlst/genus_species_ref_ANI|printf 'autoMLST genus_species_ANI: %s\n' "$(cat)" >> $STRAINNAME.AA.log
assembly-stats -s $STRAINNAME.contigs.fasta|grep number|cut -f 3 |printf 'number of contigs: %s\n' "$(cat)" >> $STRAINNAME.AA.log
assembly-stats -s $STRAINNAME.contigs.fasta|grep total_length|cut -f 3|printf 'total assembly length: %s\n' "$(cat)" >> $STRAINNAME.AA.log
assembly-stats -s $STRAINNAME.contigs.fasta|grep longest|cut -f 3|printf 'longest contig: %s\n' "$(cat)" >> $STRAINNAME.AA.log
cat flye/assembly_info.txt >> $STRAINNAME.AA.log

#cat npgm-stats.txt >> $STRAINNAME.AA.log
#rm npgm-stats.txt 
cat n50 >> $STRAINNAME.AA.log
cat ${STRAINNAME}.graph.gfa |grep ^S|cut -f 3|awk '{ print length }'|sed 's/$/ nt/' > ${STRAINNAME}.edgelength
cat ${STRAINNAME}.graph.gfa |grep ^S|cut -f 2,4|sed 's/\t/\tcov /'|paste - ${STRAINNAME}.edgelength >> $STRAINNAME.AA.log
rm ${STRAINNAME}.edgelength
cat ${STRAINNAME}.graph.gfa |grep -P '^L|^P' >> $STRAINNAME.AA.log
ls *antiSMASH/*region*gbk|wc -l |sed 's/^/BGC Regions: /' >> $STRAINNAME.AA.log
cat ${STRAINNAME}_prokka_actinoannotPFAM/${STRAINNAME}.log |grep unann|tail -n 1|sed 's/.* There are //1' >> $STRAINNAME.AA.log
cat busco/short*|tail -n 7 >> $STRAINNAME.AA.log

#document tool versions
basename $0 >> $STRAINNAME.AA.log
which autoMLST|sed 's/.*autoMLST\//autoMLST v/1'|sed 's/\/.*//' >> $STRAINNAME.AA.log
Bandage --version|printf 'Bandage v%s\n' "$(cat)"|sed 's/Version: //' >> $STRAINNAME.AA.log
blastn -version|grep -v Pack|sed 's/: / v/' >> $STRAINNAME.AA.log
busco --version|sed 's/ / v/' >> $STRAINNAME.AA.log
flye --version 2>&1|printf 'flye v%s\n' "$(cat)" >> $STRAINNAME.AA.log

#POLCA-version for non-module
#echo "POLCA start"
#OLD_PATH=$PATH
#export PATH=`get_masurca_path.sh`:$PATH
#OLD_LD_PATH=LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=`get_masurca_lib_path.sh`:$LD_LIBRARY_PATH
#polca.sh --version |printf 'masurca/polca v%s\n' "$(cat)" >> $STRAINNAME.AA.log 
#export PATH=$OLD_PATH
#export LD_LIBRARY_PATH=$OLD_LD_PATH

prokka --version 2>&1|grep prokka_|sed 's/ / v/' >> $STRAINNAME.AA.log
antismash --version >> $STRAINNAME.AA.log

#clean up temporary files
rm -r  automlst busco* flye ${STRAINNAME}_antiSMASH
#rm npgm-contigged.fa.fai npgm-contigger.fa
rm n50

#make little celebratory statement marking the finishing of the pipeline
echo 'Your assembly and annotation of'
echo $STRAINNAME
echo 'has successfully finished'

