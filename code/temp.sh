#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
sample=$3 # will appear in output data in sample field and will appear in output file names
outdir=$4 # output files will be written here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
in_bam_dir=$6 # the bam file for this sample, for records not having depth info, will use samtools to fetch the depth 10bp either side of each breakend of the clean-intron-deletion

for infile in "${outdir}"/A020049117311.clean_intron_deletions_with_strand_and_vaf_and_depth.tsv; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  if [[ $in_bam_dir == "." ]]; then
    in_bam="."
  else
    in_bam="${in_bam_dir}"/"${sample}".bam
  fi

  outfile="${outdir}"/"${sample}".clean_intron_deletions_with_strand_and_vaf_and_depth_and_numExons.tsv

  echo './temp_run_2b.sh' $genome $cohort $sample $outdir $tmpdir $in_bam $infile $outfile
  ./temp_run_2b.sh $genome $cohort $sample $outdir $tmpdir $in_bam $infile $outfile
  echo ''

done


