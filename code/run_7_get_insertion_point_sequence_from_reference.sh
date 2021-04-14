#!/bin/bash
set -euo pipefail

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
sample=$3 # will appear in output data in sample field and will appear in output file names
outdir=$4 # output files will be written here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
ref_fasta=$6 # the reference fasta that the bam file was aligned to, the bam file used to call structural variants by gridss
python_path=$7 # needed for the python program to use python packages, will be used in command: export PYTHONPATH=$python_path

export PYTHONPATH=$python_path

sw="."

infile="${outdir}"/"${sample}".insertion_points_for_retrocopied_genes_including_using_other_samples.tsv
outfile="${outdir}"/"${sample}".insertion_points_for_retrocopied_genes_including_using_other_samples.insertion_pt_ref_seq.tsv

echo ''
echo 'python3 get_insertion_point_sequence_from_reference.py --input' $infile '--output' $outfile '--ref_fasta' $ref_fasta
python3 get_insertion_point_sequence_from_reference.py --input $infile --output $outfile --ref_fasta $ref_fasta
echo ''

echo 'outfile:' $outfile
echo ''



