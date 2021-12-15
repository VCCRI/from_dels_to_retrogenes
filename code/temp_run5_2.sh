#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
sample=$3 # will appear in output data in sample field and will appear in output file names
outdir=$4 # output files will be written here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
infile_from_gridss=$6 # the input gridss vcf file, will be used to create file names of output data in this pipeline. It needs to be bgzip compressed as <name>.vcf.gz and tabix -p vcf indexed <name>.vcf.gz.tbipython_path=$7 # needed for the python program to use python packages, will be used in command: export PYTHONPATH=$python_path
python_path=$7 # needed for the python program to use python packages, will be used in command: export PYTHONPATH=$python_path
in_bam=$8 # preferably the bam file for this sample, sequenced to 30x or 40x depth, or else a representative bam file of similar depth and same genome version, 
          # to see regions of unusually high depth that must be avoided when looking for insertion sites
in_dels=$9
in_prev_outfile="${10}"
outfile="${11}"

echo 'genome' $genome
echo 'cohort' $cohort
echo 'sample' $sample
echo 'outdir' $outdir
echo 'tmpdir' $tmpdir
echo 'infile_from_gridss' $infile_from_Gridss
echo 'python_path' $python_path
echo 'in_bam' $in_bam
echo 'in_dels' $in_dels
echo 'in_prev_outfile' $in_prev_outfile
echo 'outfile' $outfile
echo ''

currdir=$(pwd)

in_sv=$infile_from_gridss

if [[ $genome == "hg19" ]]; then
  gene_regions="${currdir}"/../reference_data/hg19_UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq_and_Gencode_20211116.bed
  overly_mapped_regions="${currdir}"/../reference_data/hg19_insertion_points_overly_mapped_hg19.txt
else
  gene_regions="${currdir}"/../reference_data/hg38_UCSC_GRCh38_GenesAndGenePredictions_genes_RefSeq_and_Gencode_20211116.bed
  overly_mapped_regions="${currdir}"/../reference_data/hg38_insertion_points_overly_mapped_hg38.txt
fi

export PYTHONPATH=$python_path

sw="."

# Sort by gene name, then by clean_intron_deletion coords, because the one clean_intron_deletion can belong to 2 genes
# and thus leaving the introns sorted by coords means that the introns of one gene can be interrupted by the introns of another gene.
# eg. 12:124457849-124458431 is an intron of ZNF664 and ZNF664-RFLNA
in_dels_basename=$(basename $in_dels)
tmp_in_dels="${tmpdir}"/"${in_dels_basename}"
tmp_hdr="${tmpdir}"/tmp_hdr.txt
head -n 1 $in_dels > $tmp_hdr
echo 'grep -v' '^cohort' $in_dels '| sort -k7,7 -k3,3V -k4,4V -k5,5V | cat' $tmp_hdr '- >' $tmp_in_dels
grep -v '^cohort' $in_dels | sort -k7,7 -k3,3V -k4,4V -k5,5V | cat $tmp_hdr - > $tmp_in_dels
echo ''

keep_looking_for_more_insertion_points=1
max_num_rounds_of_looking_for_insertion_points=4
i=1

basename_in_prev_outfile=$(basename $in_prev_outfile)
prev_outfile="${tmpdir}"/prev_outfile_"${basename_in_prev_outfile}"
echo 'cp' $in_prev_outfile $prev_outfile
cp $in_prev_outfile $prev_outfile
echo ''

outfile_basename=$(basename $outfile)

while [[ $keep_looking_for_more_insertion_points == "1" ]]; do

  blacklist_regions="${tmpdir}/blacklist_regions_for_${outfile_basename}"
  tmp_outfile="${tmpdir}/tmp_${outfile_basename}"

  # add existing insertion_points to blacklist_regions

  echo 'cut -d$''\t' '-f17-19' $prev_outfile '| grep -v' '^insertion' '| sort -k1,1V -k2,2V -k3,3V | grep -P -v' '^\t' '>' $blacklist_regions
  cut -d$'\t' -f17-19 $prev_outfile | grep -v '^insertion' | sort -k1,1V -k2,2V -k3,3V | grep -P -v '^\t' > $blacklist_regions
  echo ''

  # look for more insertion points

  if [[ $in_bam == "." ]]; then

    echo 'python3' $sw'/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels' $tmp_in_dels '--in_sv' $in_sv '--gene_regions' $gene_regions '-o' $tmp_outfile '--gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --overly_mapped_regions' $overly_mapped_regions '--overly_mapped_depth 500 --max_dist_btwn_ins_pts 500 --blacklist_regions' $blacklist_regions
    python3 $sw/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels $tmp_in_dels --in_sv $in_sv --gene_regions $gene_regions -o $tmp_outfile --gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --overly_mapped_regions $overly_mapped_regions --overly_mapped_depth 500 --max_dist_btwn_ins_pts 500 --blacklist_regions $blacklist_regions
    echo ''

  else

    echo 'python3' $sw'/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels' $tmp_in_dels '--in_sv' $in_sv '--gene_regions' $gene_regions '-o' $tmp_outfile '--gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --overly_mapped_regions' $overly_mapped_regions '--bam' $in_bam '--overly_mapped_depth 500 --max_dist_btwn_ins_pts 500 --blacklist_regions' $blacklist_regions
    python3 $sw/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels $tmp_in_dels --in_sv $in_sv --gene_regions $gene_regions -o $tmp_outfile --gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --overly_mapped_regions $overly_mapped_regions --bam $in_bam --overly_mapped_depth 500 --max_dist_btwn_ins_pts 500 --blacklist_regions $blacklist_regions
    echo ''

  fi

  # Don't output genes for which an insertion point was not found.
  # Add the insertion points we have just found to the previously found insertion points.

  echo 'cat' $tmp_outfile '| awk' 'BEGIN {FS="\t";OFS="\t"} {if (($17 != "") && ($17 != ".") && (NR>1)) {print $0}}' '| cat' $prev_outfile '- >' $outfile
  cat $tmp_outfile | awk 'BEGIN {FS="\t";OFS="\t"} {if (($17 != "") && ($17 != ".") && (NR>1)) {print $0}}' | cat $prev_outfile - > $outfile
  echo ''

  echo 'cp' $outfile $prev_outfile
  cp $outfile $prev_outfile
  echo ''

  # If we didn't find any new insertion points, then stop looking.

  echo 'num_found=$(cat' $tmp_outfile '| awk' 'BEGIN {FS="\t";OFS="\t"} {if (($17 != "") && ($17 != ".") && (NR>1)) {print $0}}' '| wc -l | cut -d'':' '-f1)'
  num_found=$(cat $tmp_outfile | awk 'BEGIN {FS="\t";OFS="\t"} {if (($17 != "") && ($17 != ".") && (NR>1)) {print $0}}' | wc -l | cut -d':' -f1)
  if [[ $num_found -eq 0 ]]; then
    keep_looking_for_more_insertion_points=0
  fi
  echo 'num_found:' $num_found 'keep_looking_for_more_insertion_points:' $keep_looking_for_more_insertion_points

  # To ensure that we don't go into an infinite loop, there is a maximum number of loops we will do.

  ((i=i+1))
  if [[ $i -gt $max_num_rounds_of_looking_for_insertion_points ]]; then 
    keep_looking_for_more_insertion_points=0
  fi
  echo 'i:' $i 'keep_looking_for_more_insertion_points:' $keep_looking_for_more_insertion_points

done

echo ''
echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''


