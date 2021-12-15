#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
sample=$3 # will appear in output data in sample field and will appear in output file names
outdir=$4 # output files will be written here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
in_bam=$6 # the bam file for this sample, for records not having depth info, will use samtools to fetch the depth 10bp either side of each breakend of the clean-intron-deletion
infile=$7
outfile=$8

echo 'genome' $genome
echo 'cohort' $cohort
echo 'sample' $sample
echo 'outdir' $outdir
echo 'tmpdir' $tmpdir
echo 'in_bam' $in_bam
echo 'infile' $infile
echo 'outfile' $outfile
echo ''

gene_max_num_exons=../reference_data/hg38_UCSC_GRCh38_GenesAndGenePredictions_numExons_RefSeq_and_Gencode_20211116.txt
if [[ $genome == "hg19" ]]; then
  gene_max_num_exons=../reference_data/hg19_UCSC_GRCh37_GenesAndGenePredictions_numExons_RefSeq_and_Gencode_20211116.txt
fi

echo 'Rscript add_gene_max_num_exons_to_clean_intron_deletions.R' $infile $outfile $gene_max_num_exons
Rscript add_gene_max_num_exons_to_clean_intron_deletions.R $infile $outfile $gene_max_num_exons
echo ''

echo ''
echo 'outfile:' $outfile
echo ''

