#!/bin/bash

indir=$1

outdir=$indir

outfile="${outdir}"/All__insertion_points_for_deletions_that_are_retrocopied_genes.sort_by_sample.txt
outfile2="${outdir}"/All__insertion_points_for_deletions_that_are_retrocopied_genes.sort_by_gene.txt

file1=$(ls -1 $indir | head -n 1)
read -r hdr1 < "${indir}"/"${file1}"
echo -e "$hdr1" > $outfile

:>$outfile
for infile in "${indir}"/*.insertion_points_for_deletions_that_are_retrocopied_genes*.tsv; do

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


