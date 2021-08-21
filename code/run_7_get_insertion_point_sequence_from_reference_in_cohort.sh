#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
dummy_sample=$3 # this parameter is not used, fill it with a dot, this script will extract each sample id from each file name and pass it on to next script.
outdir=$4 # output files will be written here and input files will be found here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
ref_fasta=$6
python_path=$7

for infile in "${outdir}"/*.insertion_points_for_deletions_that_are_retrocopied_genes_using_all_samples_insertion_points.tsv; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  outfile="${outdir}"/"${sample}".insertion_points_for_deletions_that_are_retrocopied_genes_with_insertion_point_reference_sequence.tsv

  ./run_7_get_insertion_point_sequence_from_reference.sh $genome $cohort $sample $outdir $tmpdir $ref_fasta $python_path $infile $outfile

done
