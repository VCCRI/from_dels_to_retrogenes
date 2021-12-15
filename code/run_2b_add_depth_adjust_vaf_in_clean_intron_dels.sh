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

infile_basename=$(basename $infile)
tmp_outfile="${tmpdir}"/"${infile_basename%.tsv}"_and_depth.tsv

awk -v in_bam="$in_bam" '
function max(a, b) {
  if (a>b) {
    return a
  } else {
    return b
  }
}
function fix_depth_pair(depth_pair) {
  if ((depth_pair == ".") || (depth_pair == "")) {
    depth_pair = "."
  } else {
    num_bits = split(depth_pair, array, ",")
    if (num_bits >= 2) {
      # The input data contains two depth numbers (or more, which would be unexpected), separated by a comma.
      do_nothing=1 # depth_pair doesnt need fixing
    } else { # (num_bits == 1)
      # There is a bug in the input data. Sometimes the comma between the two numbers is missing.
      # If the length of the depth_pair is an even number, then the comma was in the middle.
      # Otherwise, figure out where the comma would have been, such that the depths are as low as possible.
      if ((length(depth_pair) % 2) == 0) { # length is an even number
        halfway = length(depth_pair) / 2
        num1 = substr(depth_pair, 1, halfway)
        num2 = substr(depth_pair, (halfway+1))
        depth_pair = num1","num2
      } else { # length is an odd number
        halfway = int(length(depth_pair) / 2)
        num1 = substr(depth_pair, 1, halfway)
        num2 = substr(depth_pair, (halfway+1))
        max_depth_12 = max(num1, num2)
        num3 = substr(depth_pair, 1, (halfway+1))
        num4 = substr(depth_pair, (halfway+2))
        max_depth_34 = max(num3, num4)
        if (max_depth_12 > max_depth_34) {
          depth_pair = num3","num4
        } else {
          depth_pair = num1","num2
        }
      }
    }
  }
  return depth_pair
}
BEGIN {FS="\t";OFS="\t"}
{
  batch=$1
  sample=$2
  clean_intron_del_chrom=$3
  clean_intron_del_start=$4
  clean_intron_del_end=$5
  clean_intron_del_svtype=$6
  clean_intron_del_gene=$7
  gene_intron_chrom=$8
  gene_intron_start=$9
  gene_intron_end=$10
  gene_intron_strand=$11
  vaf=$12
  depth1=$13
  depth2=$14
  depth3=$15
  depth4=$16
  depth5=$17

  # output to be determined by this script
  max_bnd1_depth="."
  max_bnd2_depth="."
  bam_bnd1_depths="."
  bam_bnd2_depths="."
  new_vaf=vaf

  # depth1 = gridssASRP
  # depth2 = gridssREF
  # depth3 = gridssREFPAIR
  # depth4 = gridssRP
  # depth5 = gridssSR

  depth1 = fix_depth_pair(depth1)
  depth2 = fix_depth_pair(depth2)
  depth3 = fix_depth_pair(depth3)
  depth4 = fix_depth_pair(depth4)
  depth5 = fix_depth_pair(depth5)

  if (depth1 != ".") {
    split(depth1, array, ",")
    max_bnd1_depth = array[1]
    max_bnd2_depth = array[2]
  }
  if (depth2 != ".") {
    split(depth2, array, ",")
    max_bnd1_depth = max(array[1], max_bnd1_depth)
    max_bnd2_depth = max(array[2], max_bnd2_depth)
  }
  if (depth3 != ".") {
    split(depth3, array, ",")
    max_bnd1_depth = max(array[1], max_bnd1_depth)
    max_bnd2_depth = max(array[2], max_bnd2_depth)
  }
  if (depth4 != ".") {
    split(depth4, array, ",")
    max_bnd1_depth = max(array[1], max_bnd1_depth)
    max_bnd2_depth = max(array[2], max_bnd2_depth)
  }
  if (depth5 != ".") {
    split(depth5, array, ",")
    max_bnd1_depth = max(array[1], max_bnd1_depth)
    max_bnd2_depth = max(array[2], max_bnd2_depth)
  }

  if (in_bam != ".") {

    pos = clean_intron_del_start - 10
    coords = clean_intron_del_chrom":"pos"-"pos
    cmd = "samtools view "in_bam" "coords" | wc -l"
    cmd | getline bam_bnd1_depths1
    close(cmd)

    pos = clean_intron_del_start + 10
    coords = clean_intron_del_chrom":"pos"-"pos
    cmd = "samtools view "in_bam" "coords" | wc -l"
    cmd | getline bam_bnd1_depths2
    close(cmd)

    bam_bnd1_depths=bam_bnd1_depths1","bam_bnd1_depths2

    pos = clean_intron_del_end - 10
    coords = clean_intron_del_chrom":"pos"-"pos
    cmd = "samtools view "in_bam" "coords" | wc -l"
    cmd | getline bam_bnd2_depths1
    close(cmd)

    pos = clean_intron_del_end + 10
    coords = clean_intron_del_chrom":"pos"-"pos
    cmd = "samtools view "in_bam" "coords" | wc -l"
    cmd | getline bam_bnd2_depths2
    close(cmd)

    bam_bnd2_depths=bam_bnd2_depths1","bam_bnd2_depths2

    if (vaf == 0) {

      bam_bnd1_vaf1 = 0
      if (bam_bnd1_depths1 > bam_bnd1_depths2) {
        bam_bnd1_vaf = 1.0 - (bam_bnd1_depths2/bam_bnd1_depths1)
      } else {
        bam_bnd1_vaf = 1.0 - (bam_bnd1_depths1/bam_bnd1_depths2)
      }

      bam_bnd2_vaf2 = 0
      if (bam_bnd2_depths1 > bam_bnd2_depths2) {
        bam_bnd2_vaf = 1.0 - (bam_bnd2_depths2/bam_bnd2_depths1)
      } else {
        bam_bnd2_vaf = 1.0 - (bam_bnd2_depths1/bam_bnd2_depths2)
      }

      if (bam_bnd1_vaf > bam_bnd2_vaf) {
        new_vaf = bam_bnd1_vaf
      } else {
        new_vaf = bam_bnd2_vaf
      }
    }
  }

  print $0, max_bnd1_depth, max_bnd2_depth, bam_bnd1_depths, bam_bnd2_depths, new_vaf

}' $infile > $tmp_outfile

gene_max_num_exons=../reference_data/hg38_UCSC_GRCh38_GenesAndGenePredictions_numExons_RefSeq_and_Gencode_20211116.txt
if [[ $genome == "hg19" ]]; then
  gene_max_num_exons=../reference_data/hg19_UCSC_GRCh37_GenesAndGenePredictions_numExons_RefSeq_and_Gencode_20211116.txt
fi

echo 'Rscript add_gene_max_num_exons_to_clean_intron_deletions.R' $tmp_outfile $outfile $gene_max_num_exons
Rscript add_gene_max_num_exons_to_clean_intron_deletions.R $tmp_outfile $outfile $gene_max_num_exons
echo ''

echo ''
echo 'outfile:' $outfile
echo ''

