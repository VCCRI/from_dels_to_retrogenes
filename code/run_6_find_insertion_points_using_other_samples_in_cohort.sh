#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
dummy_sample=$3 # this parameter is not used, fill it with a dot, this script will extract each sample id from each file name and pass it on to next script.
outdir=$4 # output files will be written here and input files will be found here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
indir_infile_from_gridss=$6 # the input gridss vcf file, will be used to create file names of output data in this pipeline
python_path=$7
indir_in_bam=$8
other_samples_retrocopies=$9

# infile_from_gridss is expected to be /my/directory/of/gridss/output/*.gridss.vcf or *.gridss.vcf.gz, where * is the sample id.
# indir_in_bam is expected to be /my/directory/of/bam/file/*.bam, where * is the sample id, or is a representative sample bam file to be used for all samples, or is "." to denote no bam available.

other_samples_retrocopies=/g/data/jb96/emmrat/clean_intron_deletions/work_2021_august/results/All/All__insertion_points_for_deletions_that_are_retrocopied_genes_N.sort_by_sample.txt

for infile in "${outdir}"/*.insertion_points_for_deletions_that_are_retrocopied_genes_N.tsv; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  outfile="${outdir}"/"${sample}".insertion_points_for_deletions_that_are_retrocopied_genes_using_all_samples_insertion_points.tsv

  infile_from_gridss="${indir_infile_from_gridss/\*/$sample}"
  in_bam="${indir_in_bam/\*/$sample}"

  ./run_6_find_insertion_points_using_other_samples.sh $genome $cohort $sample $outdir $tmpdir $infile_from_gridss $python_path $in_bam $infile $other_samples_retrocopies $outfile

done

