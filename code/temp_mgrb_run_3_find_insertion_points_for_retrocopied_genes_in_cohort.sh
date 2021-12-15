#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
dummy_sample=$3 # this parameter is not used, fill it with a dot, this script will extract each sample id from each file name and pass it on to next script.
outdir=$4 # output files will be written here and input files will be found here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
indir_for_gridss=$6 # the input gridss vcf file. It needs to be bgzip compressed as <name>.vcf.gz and tabix -p vcf indexed <name>.vcf.gz.tbi, and sample is the only sample in this vcf file.
python_path=$7
indir_for_in_bam=$8

# infile_from_gridss is expected to be /my/directory/of/gridss/output/*.gridss.vcf or *.gridss.vcf.gz, where * is the sample id.
# indir_in_bam is expected to be /my/directory/of/bam/file/*.bam, where * is the sample id, or is a representative sample bam file to be used for all samples, or is "." to denote no bam available.

for infile in "${outdir}"/BAAOC*.clean_intron_deletions_with_strand_and_vaf_and_depth_and_numExons.tsv; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  outfile="${outdir}"/"${sample}".insertion_points_for_deletions_that_are_retrocopied_genes_1_new.tsv

  infile_from_gridss="${indir_for_gridss}"/"${sample}".gridss.vcf.gz
  if [[ $indir_for_in_bam == "." ]]; then
    in_bam=../input_files/MGRB/AAABL.bam
  else
    in_bam="${indir_for_in_bam}"/"${sample}".bam
    if [ -f "$in_bam" ]; then
      do_nothing=1
    else
      in_bam=../input_files/MGRB/AAABL.bam
    fi
  fi

  ./run_3_find_insertion_points_for_retrocopied_genes.sh $genome $cohort $sample $outdir $tmpdir $infile_from_gridss $python_path $in_bam $infile $outfile

done


