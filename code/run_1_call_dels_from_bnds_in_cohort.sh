#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
dummy_sample=$3 # this parameter is not used, fill it with a dot, this script will extract each sample id from each file name and pass it on to next script.
outdir=$4 # output files will be written here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
infile_from_gridss=$6 # the input gridss vcf files, with an asterisk as the first character of the basename, which is where the sample id will be and will be extracted from.
# infile_from_gridss is expected to be /my/directory/of/gridss/output/*.gridss.vcf or *.gridss.vcf.gz, where * is the sample id.

for infile in $infile_from_gridss; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  outfile="${outdir}"/"${sample}".structural_variants_and_BND.tsv

  run_1_call_dels_from_bnds.sh $genome $cohort $sample $outdir $tmpdir $infile_from_gridss $outfile

done


