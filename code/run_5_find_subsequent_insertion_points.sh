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
round_of_insertion_point_output=$9

in_insertpts="${outdir}"/All__insertion_points_for_deletions_that_are_retrocopied_genes.sort_by_sample.txt
in_dels="${outdir}"/"${sample}".clean_intron_deletions_with_strand_and_vaf.tsv
in_sv=$infile_from_gridss

outfile="${outdir}"/"${sample}".insertion_points_for_deletions_that_are_retrocopied_genes_"${round_of_insertion_point_output}".tsv

if [[ $genome == "hg19" ]]; then
  gene_regions=../reference_data/hg19_UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq_20200324.txt
else
  gene_regions=../reference_data/hg38_UCSC_GRCh38_GenesAndGenePredictions_genes_RefSeq_20200901.txt
fi

if [[ $genome == "hg19" ]]; then
  blacklist_regions=../reference_data/hg19_insertion_points_blacklist_hg19.txt
else
  blacklist_regions=../reference_data/hg38_insertion_points_blacklist_hg38.txt
fi

export PYTHONPATH=$python_path

sw="."

if [ -f "$in_sv" ]; then
  do_nothing=1
else
  in_sv_vcf="${in_sv%.gz}"
  echo 'bgzip -f' $in_sv_vcf
  bgzip -f $in_sv_vcf
  echo ''
fi

if [ -f "$in_sv".tbi ]; then
  do_nothing=1
else
  echo 'tabix -p vcf' $in_sv
  tabix -p vcf $in_sv
  echo ''
fi

in_dels_basename=$(basename $in_dels)
blacklist_regions_basename=$(basename $blacklist_regions)
tmp_blacklist_regions="${tmpdir}/${blacklist_regions_basename}"
outfile_basename=$(basename $outfile)
tmp_outfile="${tmpdir}/${outfile_basename}"

# add existing insertion_points to blacklist_regions

echo 'cut -d$''\t' '-f9-11' $in_insertpts '| grep -v' '^insertion' '| cat' $blacklist_regions '- | sort -k1,1V -k2,2V -k3,3V | grep -P -v' '^\t' '>' $tmp_blacklist_regions
cut -d$'\t' -f9-11 $in_insertpts | grep -v '^insertion' | cat $blacklist_regions - | sort -k1,1V -k2,2V -k3,3V | grep -P -v '^\t' > $tmp_blacklist_regions
echo ''

# Sort by gene name, then by clean_intron_deletion coords, because the one clean_intron_deletion can belong to 2 genes
# and thus leaving the introns sorted by coords means that the introns of one gene can be interrupted by the introns of another gene.
# eg. 12:124457849-124458431 is an intron of ZNF664 and ZNF664-RFLNA
in_dels_basename=$(basename $in_dels)
tmp_in_dels="${tmpdir}"/"${in_dels_basename}"
echo 'sort -k7,7 -k3,3V -k4,4V -k5,5V' $in_dels '>' $tmp_in_dels
sort -k7,7 -k3,3V -k4,4V -k5,5V $in_dels > $tmp_in_dels
echo ''

# look for more insertion points

echo 'python3' $sw'/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels' $tmp_in_dels '--in_sv' $in_sv '--gene_regions' $gene_regions '-o' $tmp_outfile '--gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --bam' $in_bam '--blacklist_regions' $tmp_blacklist_regions '--blacklist_depth 140'
python3 $sw/find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels $tmp_in_dels --in_sv $in_sv --gene_regions $gene_regions -o $tmp_outfile --gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --bam $in_bam --blacklist_regions $tmp_blacklist_regions --blacklist_depth 140
echo ''

# Don't output genes for which an insertion point was not found.

echo 'cat' $tmp_outfile '| awk' 'BEGIN {FS="\t";OFS="\t"} {if ($9 != "") {print $0}}' '>' $outfile
cat $tmp_outfile | awk 'BEGIN {FS="\t";OFS="\t"} {if ($9 != "") {print $0}}' > $outfile
echo ''

echo ''
echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''


