#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
sample=$3 # will appear in output data in sample field and will appear in output file names
outdir=$4 # output files will be written here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
infile_from_gridss=$6 # the input gridss vcf file, will be used to create file names of output data in this pipeline
python_path=$7 # needed for the python program to use python packages, will be used in command: export PYTHONPATH=$python_path
in_bam=$8 # preferably the bam file for this sample, sequenced to 30x or 40x depth, or else a representative bam file of similar depth and same genome version, 
          # to see regions of unusually high depth that must be avoided when looking for insertion sites
infile=$9
other_samples_retrocopies="${10}"
outfile="${11}"

in_sv=$infile_from_gridss

if [[ $genome == "hg19" ]]; then
  gene_regions=../reference_data/hg19_UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq_20200324.txt
  blacklist_regions=../reference_data/hg19_insertion_points_blacklist_hg19.txt
else
  gene_regions=../reference_data/hg38_UCSC_GRCh38_GenesAndGenePredictions_genes_RefSeq_20200901.txt
  blacklist_regions=../reference_data/hg38_insertion_points_blacklist_hg38.txt
fi

export PYTHONPATH=$python_path

sw="."

# In subsequent processing, extra fields need to be added. Add them now so that input and output files have same format.

infile_basename=$(basename $infile)
tmp_infile="${tmpdir}"/"${infile_basename}"

echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {left="is_left_BND_present";right="is_right_BND_present"} else { left="";right=""; if ($10!=""){left="left_BND_is_present"}; if ($11!=""){right="right_BND_is_present"} } printf $1; for (i=2; i<=NF; ++i) {printf OFS $i}; printf OFS left OFS right RS }' $infile '| tr -d $''\r' '>' $tmp_infile
#
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {left="is_left_BND_present";right="is_right_BND_present"} else { left="";right=""; if ($10!=""){left="left_BND_is_present"}; if ($11!=""){right="right_BND_is_present"} } printf $1; for (i=2; i<=NF; ++i) {printf OFS $i}; printf OFS left OFS right RS }' $infile | tr -d $'\r' > $tmp_infile
#
echo ''

other_basename=$(basename $other_samples_retrocopies)
tmp_other="${tmpdir}"/"${other_basename}"

echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {left="is_left_BND_present";right="is_right_BND_present"} else { left="";right=""; if ($10!=""){left="left_BND_is_present"}; if ($11!=""){right="right_BND_is_present"} } printf $1; for (i=2; i<=NF; ++i) {printf OFS $i}; printf OFS left OFS right RS }' $other_samples_retrocopies '| tr -d $''\r' '>' $tmp_other
#
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {left="is_left_BND_present";right="is_right_BND_present"} else { left="";right=""; if ($10!=""){left="left_BND_is_present"}; if ($11!=""){right="right_BND_is_present"} } printf $1; for (i=2; i<=NF; ++i) {printf OFS $i}; printf OFS left OFS right RS }' $other_samples_retrocopies | tr -d $'\r' > $tmp_other
#
echo ''

echo 'python3' $sw'/find_insertion_points_using_insertion_points_of_other_samples.py --in_retrocopies' $tmp_infile '--in_sv' $in_sv '--in_retrocopies_other_samples' $tmp_other '-o' $outfile
python3 $sw/find_insertion_points_using_insertion_points_of_other_samples.py --in_retrocopies $tmp_infile --in_sv $in_sv --in_retrocopies_other_samples $tmp_other -o $outfile
echo ''

echo ''
echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''

