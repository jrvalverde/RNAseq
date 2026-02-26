#!/bin/bash
#
# Compute Feature counts using HTSeq.
#

alignment_dir=/data/jr/work/adrian/rnaseq/alignments
gff=/data/jr/work/adrian/reference_genome/Gg_GRGc7b.gtf


# Check if reference is in GFF or GTF format

if [[ "$gff" =~ .gff$ ]] ; then gene_id='gene' ; fi
if [[ "$gff" =~ .gtf$ ]] ; then gene_id='gene_id' ; fi

for i in $alignment_dir/*bam ; do
   
   bam=`realpath $i`
   gff=`realpath $gff`
   sample=`basename $i .bam | sed -e 's/.sorted//g'`  # Remove ".bam" or ".sorted.bam"

   if [ -s "$sample.htseq.counts" ] ; then
   	continue
   fi

   echo -e "\nProcessing BAM file $sample ..."
	
   (	
  	htseq-count -f 'bam' -r 'pos' -t 'exon' -i $gene_id -m 'union' \
 	$bam $gff > $sample.htseq.counts
   	
	echo "Finished Feature Count: $sample"
	
   ) 2> log.err &
   
done
