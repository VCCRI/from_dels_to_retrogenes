#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
sample=$3 # will appear in output data in sample field and will appear in output file names
outdir=$4 # output files will be written here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
infile_from_gridss=$6 # the input gridss vcf file, will be used to create file names of output data in this pipeline
python_path=$7 # needed for the python program to use python packages, will be used in command: export PYTHONPATH=$python_path
in_bam=$8 # not used by this script
infile=$9
other_samples_retrocopies="${10}"
outfile="${11}"

in_sv=$infile_from_gridss

if [[ $genome == "hg19" ]]; then
  gene_regions=../reference_data/hg19_UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq_and_Gencode_20211116.bed
  blacklist_regions=../reference_data/hg19_insertion_points_blacklist_hg19.txt
else
  gene_regions=../reference_data/hg38_UCSC_GRCh38_GenesAndGenePredictions_genes_RefSeq_and_Gencode_20211116.bed
  blacklist_regions=../reference_data/hg38_insertion_points_blacklist_hg38.txt
fi

export PYTHONPATH=$python_path

sw="."

# Use insertion points of other samples to find insertion points not yet identified in a given sample,
# also using incomplete evidence in this given sample for the insertion point.

echo 'python3' $sw'/find_insertion_points_using_insertion_points_of_other_samples.py --in_retrocopies' $infile '--in_sv' $in_sv '--in_retrocopies_other_samples' $other_samples_retrocopies '-o' $outfile
python3 $sw/find_insertion_points_using_insertion_points_of_other_samples.py --in_retrocopies $infile --in_sv $in_sv --in_retrocopies_other_samples $other_samples_retrocopies -o $outfile
echo ''

echo ''
echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''

