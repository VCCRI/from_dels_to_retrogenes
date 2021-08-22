#!/bin/bash

indir=$1

outdir=$indir

outfile="${outdir}"/All__insertion_points_for_deletions_that_are_retrocopied_genes_inspt_refseq.sort_by_sample.txt
outfile2="${outdir}"/All__insertion_points_for_deletions_that_are_retrocopied_genes_inspt_refseq.sort_by_gene.txt

:>$outfile
i=1
for infile in "${indir}"/*.insertion_points_for_deletions_that_are_retrocopied_genes_with_insertion_point_reference_sequence.tsv; do

  if [[ $i -eq 1 ]]; then
    head -n 1 $infile > $outfile
  fi
  ((i=i+1))

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"
  echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $infile "| grep -v '^$' >>" $outfile
  awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $infile | grep -v '^$' >> $outfile

done

# sort by gene, then by sample, then by cohort
echo ''
head -n 1 $outfile > $outfile2
echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $outfile '| sort -k3,3 -k2,2 -k1,1 >>' $outfile2
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $outfile | sort -k3,3 -k2,2 -k1,1 >> $outfile2

echo ''
echo 'outfile:' $outfile
echo 'outfile2:' $outfile2
echo ''


