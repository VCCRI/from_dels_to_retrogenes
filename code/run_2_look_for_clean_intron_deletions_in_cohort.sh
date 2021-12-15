#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
dummy_sample=$3 # this parameter is not used, fill it with a dot, this script will extract each sample id from each file name and pass it on to next script.
outdir=$4 # output files will be written here and input files will be found here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline

for infile in "${outdir}"/*.structural_variants_and_BND.tsv; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  outfile="${outdir}"/"${sample}".clean_intron_deletions_with_strand_and_vaf.tsv

  ./run_2_look_for_clean_intron_deletions.sh $genome $cohort $sample $outdir $tmpdir $infile $outfile

done


