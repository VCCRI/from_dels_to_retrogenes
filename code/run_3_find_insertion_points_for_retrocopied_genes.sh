#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
sample=$3 # will appear in output data in sample field and will appear in output file names
outdir=$4 # output files will be written here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
infile_from_gridss=$6 # the input gridss vcf file, will be used to create file names of output data in this pipeline. It needs to be bgzip compressed as <name>.vcf.gz and tabix -p vcf indexed <name>.vcf.gz.tbi
python_path=$7 # needed for the python program to use python packages, will be used in command: export PYTHONPATH=$python_path
in_bam=$8 # preferably the bam file for this sample, sequenced to 30x or 40x depth, or else a representative bam file of similar depth and same genome version, 
          # to see regions of unusually high depth that must be avoided when looking for insertion sites
infile=$9
outfile="${10}"

currdir=$(pwd)

in_dels=$infile
in_sv=$infile_from_gridss

if [[ $genome == "hg19" ]]; then
  gene_regions="${currdir}"/../reference_data/hg19_UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq_20200324.txt
  overly_mapped_regions="${currdir}"/../reference_data/hg19_insertion_points_overly_mapped_hg19.txt
else
  gene_regions="${currdir}"/../reference_data/hg38_UCSC_GRCh38_GenesAndGenePredictions_genes_RefSeq_20200901.txt
  overly_mapped_regions="${currdir}"/../reference_data/hg38_insertion_points_overly_mapped_hg38.txt
fi

export PYTHONPATH=$python_path

sw="."

# Sort by gene name, then by clean_intron_deletion coords, because the one clean_intron_deletion can belong to 2 genes
# and thus leaving the introns sorted by coords means that the introns of one gene can be interrupted by the introns of another gene.
# eg. 12:124457849-124458431 is an intron of ZNF664 and ZNF664-RFLNA
in_dels_basename=$(basename $in_dels)
tmp_in_dels="${tmpdir}"/"${in_dels_basename}"
echo 'sort -k7,7 -k3,3V -k4,4V -k5,5V' $in_dels '>' $tmp_in_dels
sort -k7,7 -k3,3V -k4,4V -k5,5V $in_dels > $tmp_in_dels
echo ''

if [[ $in_bam == "." ]]; then

  echo 'python3' $sw'/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels' $tmp_in_dels '--in_sv' $in_sv '--gene_regions' $gene_regions '-o' $outfile '--gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --overly_mapped_regions' $overly_mapped_regions '--overly_mapped_depth 140'
  python3 $sw/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels $tmp_in_dels --in_sv $in_sv --gene_regions $gene_regions -o $outfile --gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --overly_mapped_regions $overly_mapped_regions --overly_mapped_depth 140
  echo ''

else

  echo 'python3' $sw'/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels' $tmp_in_dels '--in_sv' $in_sv '--gene_regions' $gene_regions '-o' $outfile '--gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --overly_mapped_regions' $overly_mapped_regions '--bam' $in_bam '--overly_mapped_depth 140'
  python3 $sw/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels $tmp_in_dels --in_sv $in_sv --gene_regions $gene_regions -o $outfile --gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --overly_mapped_regions $overly_mapped_regions --bam $in_bam --overly_mapped_depth 140
  echo ''

fi

echo ''
echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''


