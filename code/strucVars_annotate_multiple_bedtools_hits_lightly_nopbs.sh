#!/bin/bash

set -euo pipefail

sample=$1
infile_tsv=$2
outdir=$3
outfile=$4
outprefix=$5
sw_and_refs=$6
cohort=$7
tmpdir=$8

# This script uses bedtools to find the intersection between each row of the input data and the reference tables.
# Bedtools produces multiple output rows for a given input row when it overlaps with more than one entry in the reference table.
# Often the reference table value is the same value for all those overlaps. Eg. the reference table contains multiple transcript for the same region for the same gene.
# This script then collapses those multiple rows to one row for each different value it overlaps with.
# It does this by sorting, then running a script to fill in the unique value(s) of the last column(s) added by the bedtools intersect.
# Please note that this script uses linux sort on the first 5 columns instead of bedtools.
# Rows having the same beginning key are then collapsed to one row.
# For multiple rows of the following records where columns 1 to 3 are the same:
# chr17   77522483        81970118        T       <DUP:TANDEM>
# chr17   77522483        77522536        T       TATACACATATATATATATATATATATATACACACATATATATATATATATATAC
# bedtools sort does not sort on the fifth column.
# The consequences of this is that multiple rows remains, which can exponentially amplify the number of rows for each subsequent reference table comparison.

# set environment variables for softwares and references.
# they will be set by export command in the file sw_and_refs.
. "${sw_and_refs}"

outprefix_sorted="${outprefix}".sorted.tsv
outprefix_CDSexons="${outprefix}".CDSexons.tsv
outprefix_CDSexons_collapsed="${outprefix}".CDSexons.collapsed.tsv
outprefix_genes="${outprefix}".CDSexons.entireGene.tsv
outprefix_genes_collapsed="${outprefix}".CDSexons.entireGene.collapsed.tsv
outprefix_bndInsideExon="${outprefix}".CDSexons.entireGene.bndInsideExon.tsv
outprefix_gencode_CDSexons="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.tsv
outprefix_gencode_CDSexons_collapsed="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.collapsed.tsv
outprefix_gencode_genes="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.gencodeEntireGene.tsv
outprefix_gencode_genes_collapsed="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.gencodeEntireGene.collapsed.tsv
outprefix_gencode_bndInsideExon="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.gencodeEntireGene.gencodeBndInsideExon.tsv
outprefix_EHRFr99="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.gencodeEntireGene.gencodeBndInsideExon.EHRFr99.tsv
outprefix_EHRFr99_collapsed="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.gencodeEntireGene.gencodeBndInsideExon.EHRFr99.collapsed.tsv
outprefix_segdup="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.gencodeEntireGene.gencodeBndInsideExon.EHRFr99.gnomadSV.DGV.SegDup.tsv
outprefix_cleanIntronDels="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.gencodeEntireGene.gencodeBndInsideExon.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.HGMD.FANTOM5_TSS.cleanIntronDels.tsv
outprefix_cleanIntronDels_collapsed="${outprefix}".CDSexons.entireGene.bndInsideExon.gencodeCDSexons.gencodeEntireGene.gencodeBndInsideExon.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.HGMD.FANTOM5_TSS.cleanIntronDels.collapsed.tsv

random_number=$RANDOM
script_name=$(basename $0)
outfile_basename=$(basename $outfile)



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: sort the input'
echo ''

this_input="${infile_tsv}"
echo 'sed s/^#Chr/0chrom/' "${this_input}" '| sed s/^#CHROM/0chrom/ | sed s/^CHROM/0chrom/ | sed s/^chrom/0chrom/ | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_sorted}"
sed 's/^#Chr/0chrom/' "${this_input}" | sed 's/^#CHROM/0chrom/' | sed 's/^CHROM/0chrom/' | sed 's/^chrom/0chrom/' | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_sorted}"
echo ''
this_input="${outprefix_sorted}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with ucsc_refseq_cdsexons'
echo ''

# head $ucsc_refseq_cdsexons
# chr1	67096251	67096321	C1orf141	-
# chr1	67103237	67103382	C1orf141	-
# chr1	67111576	67111644	C1orf141	-

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".CDSexons."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
echo 'head -n 1' "${this_input}" '| sed -e s/$/\t.\t.\t.\tgene_CDSexons\tgene_CDSexons_strand/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed -e 's/$/\t.\t.\t.\tgene_CDSexons\tgene_CDSexons_strand/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${ucsc_refseq_cdsexons}" '| \'
echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- | cut -d$\t -f"1-'$cut1','$cut1_plus_4 '>' "${outprefix_CDSexons}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" | bedtools intersect -wao -a - -b "${ucsc_refseq_cdsexons}" | \
  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat "${this_hdr}" - | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" > "${outprefix_CDSexons}"
rm "${this_hdr}"
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It might be needed for any duplicate exons, otherwise same gene might appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${outprefix_CDSexons}" '>' "${outprefix_CDSexons_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${outprefix_CDSexons}" > "${outprefix_CDSexons_collapsed}"
echo ''
this_input="${outprefix_CDSexons_collapsed}"
rm $outprefix_sorted
rm $outprefix_CDSexons



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with ucsc_refseq_genes'
echo ''

# head $ucsc_refseq_genes
# chr1	11873	14409	DDX11L1	+
# chr1	14361	29370	WASH7P	-
# chr1	17368	17436	MIR6859-1	-
# chr1	29925	31295	MIR1302-2HG	+
# chr1	30365	30503	MIR1302-2	+

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".genes."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l) # 23
cut1_plus_4=$((cut1+4)) # 27
echo 'head -n 1' "${this_input}" '| sed -e s/$/\t.\t.\t.\tgene_entireGene\tgene_entireGene_strand/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed -e 's/$/\t.\t.\t.\tgene_entireGene\tgene_entireGene_strand/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${ucsc_refseq_genes}" '| \'
echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- | cut -d$\t -f"1-'$cut1','$cut1_plus_4 '>' "${outprefix_genes}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" | bedtools intersect -wao -a - -b "${ucsc_refseq_genes}" | \
  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat "${this_hdr}" - | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" > "${outprefix_genes}"
rm "${this_hdr}"
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${outprefix_genes}" '>' "${outprefix_genes_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${outprefix_genes}" > "${outprefix_genes_collapsed}"
echo ''
this_input="${outprefix_genes_collapsed}"
rm $outprefix_CDSexons_collapsed
rm $outprefix_genes



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do left BND fall inside an exon?'
echo ''

# head $ucsc_refseq_cdsexons
# chr1	67096251	67096321	C1orf141	-
# chr1	67103237	67103382	C1orf141	-
# chr1	67111576	67111644	C1orf141	-

bndInsideExon_left_bnd="${tmpdir}"/temp."${sample}".bndInsideExon_left_bnd."${outfile_basename}".bed
bndInsideExon_left_bnd_intersect="${tmpdir}"/temp."${sample}".bndInsideExon_left_bnd_intersect."${outfile_basename}".bed
bndInsideExon_left_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".bndInsideExon_left_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".bndInsideExon_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tleft_BND_inside_exon\texon_containing_left_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tleft_BND_inside_exon\texon_containing_left_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${bndInsideExon_left_bnd}"
awk 'BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" | cat "${this_hdr}" - > "${bndInsideExon_left_bnd}" || true
echo 'bedtools intersect -wao -a' "${bndInsideExon_left_bnd}" '-b' "${ucsc_refseq_cdsexons}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8} | cat' "${this_hdr}" '- >' "${bndInsideExon_left_bnd_intersect}"
bedtools intersect -wao -a "${bndInsideExon_left_bnd}" -b "${ucsc_refseq_cdsexons}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8}' | cat "${this_hdr}" - > "${bndInsideExon_left_bnd_intersect}"
echo ''
rm $bndInsideExon_left_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It is needed for the duplicate exons, otherwise same gene will appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${bndInsideExon_left_bnd_intersect}" '>' "${bndInsideExon_left_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${bndInsideExon_left_bnd_intersect}" > "${bndInsideExon_left_bnd_intersect_collapsed}"
echo ''
rm $bndInsideExon_left_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do right BND fall inside an exon?'
echo ''

bndInsideExon_right_bnd="${tmpdir}"/temp."${sample}".bndInsideExon_right_bnd."${outfile_basename}".bed
bndInsideExon_right_bnd_intersect="${tmpdir}"/temp."${sample}".bndInsideExon_right_bnd_intersect."${outfile_basename}".bed
bndInsideExon_right_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".bndInsideExon_right_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".bndInsideExon_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tright_BND_inside_exon\texon_containing_right_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tright_BND_inside_exon\texon_containing_right_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${bndInsideExon_right_bnd}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" | cat "${this_hdr}" - > "${bndInsideExon_right_bnd}" || true
echo 'bedtools intersect -wao -a' "${bndInsideExon_right_bnd}" '-b' "${ucsc_refseq_cdsexons}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8} | cat' "${this_hdr}" '- | sed s/:-1--1//g >' "${bndInsideExon_right_bnd_intersect}"
bedtools intersect -wao -a "${bndInsideExon_right_bnd}" -b "${ucsc_refseq_cdsexons}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8}' | cat "${this_hdr}" - | sed 's/:-1--1//g' > "${bndInsideExon_right_bnd_intersect}"
echo ''
rm $bndInsideExon_right_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It is needed for the duplicate exons, otherwise same gene will appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${bndInsideExon_right_bnd_intersect}" '>' "${bndInsideExon_right_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${bndInsideExon_right_bnd_intersect}" > "${bndInsideExon_right_bnd_intersect_collapsed}"
echo ''
rm $bndInsideExon_right_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: merge BND_inside_exon left BND and right BND into one file'
echo ''

echo 'awk -v left_bnd_file='"${bndInsideExon_left_bnd_intersect_collapsed}" '\'
echo '    -v right_bnd_file='"${bndInsideExon_right_bnd_intersect_collapsed}"' -f' $sw'/victorchang_scripts/combine_bnd_files_key5cols.awk \'
echo '    '"${this_input}" '>' "${outprefix_bndInsideExon}"
awk -v left_bnd_file="${bndInsideExon_left_bnd_intersect_collapsed}" \
    -v right_bnd_file="${bndInsideExon_right_bnd_intersect_collapsed}" -f $sw/victorchang_scripts/combine_bnd_files_key5cols.awk \
    "${this_input}" > "${outprefix_bndInsideExon}"
echo ''
rm $bndInsideExon_left_bnd_intersect_collapsed
rm $bndInsideExon_right_bnd_intersect_collapsed
this_input="${outprefix_bndInsideExon}"
rm $outprefix_genes_collapsed



if [[ "$ucsc_gencode_cdsexons" != "" ]]; then

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with ucsc_gencode_cdsexons'
echo ''

# head $ucsc_gencode_cdsexons
# chr1	67096251	67096321	C1orf141	-
# chr1	67103237	67103382	C1orf141	-
# chr1	67111576	67111644	C1orf141	-

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".gencode_CDSexons."${random_number}".txt
tmpfile11="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile11."${random_number}".txt
tmpfile12="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile12."${random_number}".txt
tmpfile13="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile13."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
echo 'head -n 1' "${this_input}" '| sed -e s/$/\t.\t.\t.\tgencode_gene_CDSexons\tgencode_gene_CDSexons_strand/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed -e 's/$/\t.\t.\t.\tgencode_gene_CDSexons\tgencode_gene_CDSexons_strand/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${ucsc_gencode_cdsexons}" '| \'
echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- | cut -d$\t -f"1-'$cut1','$cut1_plus_4 '>' "${outprefix_gencode_CDSexons}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" > $tmpfile11
bedtools intersect -wao -a $tmpfile11 -b "${ucsc_gencode_cdsexons}" > $tmpfile12
rm $tmpfile11
sort $tmpfile12 | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 > $tmpfile13
rm $tmpfile12
cat "${this_hdr}" $tmpfile13 | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" > "${outprefix_gencode_CDSexons}"
rm "${this_hdr}"
rm $tmpfile13
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It might be needed for any duplicate exons, otherwise same gene might appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${outprefix_gencode_CDSexons}" '>' "${outprefix_gencode_CDSexons_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${outprefix_gencode_CDSexons}" > "${outprefix_gencode_CDSexons_collapsed}"
echo ''
this_input="${outprefix_gencode_CDSexons_collapsed}"

fi



if [[ "$ucsc_gencode_genes" != "" ]]; then

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with ucsc_gencode_genes'
echo ''

# head $ucsc_refseq_genes
# chr1	11873	14409	DDX11L1	+
# chr1	14361	29370	WASH7P	-
# chr1	17368	17436	MIR6859-1	-
# chr1	29925	31295	MIR1302-2HG	+
# chr1	30365	30503	MIR1302-2	+

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".genes."${random_number}".txt
tmpfile21="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile21."${random_number}".txt
tmpfile22="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile22."${random_number}".txt
tmpfile23="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile23."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l) # 23
cut1_plus_4=$((cut1+4)) # 27
echo 'head -n 1' "${this_input}" '| sed -e s/$/\t.\t.\t.\tgencode_gene_entireGene\tgencode_gene_entireGene_strand/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed -e 's/$/\t.\t.\t.\tgencode_gene_entireGene\tgencode_gene_entireGene_strand/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${ucsc_gencode_genes}" '| \'
echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- | cut -d$\t -f"1-'$cut1','$cut1_plus_4 '>' "${outprefix_gencode_genes}"
echo 'awk' $this_input 'to' $tmpfile21
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" > $tmpfile21
bedtools intersect -wao -a $tmpfile21 -b "${ucsc_gencode_genes}" > $tmpfile22
rm $tmpfile21
sort $tmpfile22 | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | uniq > $tmpfile23
rm $tmpfile22
cat "${this_hdr}" $tmpfile23 | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" > "${outprefix_gencode_genes}"
rm "${this_hdr}"
rm $tmpfile23
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${outprefix_gencode_genes}" '>' "${outprefix_gencode_genes_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${outprefix_gencode_genes}" > "${outprefix_gencode_genes_collapsed}"
echo ''
this_input="${outprefix_gencode_genes_collapsed}"

fi



if [[ "$ucsc_gencode_cdsexons" != "" ]]; then

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do left BND fall inside a gencode exon?'
echo ''

# head $ucsc_refseq_cdsexons
# chr1	67096251	67096321	C1orf141	-
# chr1	67103237	67103382	C1orf141	-
# chr1	67111576	67111644	C1orf141	-

gencode_bndInsideExon_left_bnd="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_left_bnd."${outfile_basename}".bed
gencode_bndInsideExon_left_bnd_intersect="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_left_bnd_intersect."${outfile_basename}".bed
gencode_bndInsideExon_left_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_left_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".gencode_bndInsideExon_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tgencode_left_BND_inside_exon\tgencode_exon_containing_left_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tgencode_left_BND_inside_exon\tgencode_exon_containing_left_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${gencode_bndInsideExon_left_bnd}"
awk 'BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" | cat "${this_hdr}" - > "${gencode_bndInsideExon_left_bnd}" || true
echo 'bedtools intersect -wao -a' "${gencode_bndInsideExon_left_bnd}" '-b' "${ucsc_gencode_cdsexons}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8} | cat' "${this_hdr}" '- >' "${gencode_bndInsideExon_left_bnd_intersect}"
bedtools intersect -wao -a "${gencode_bndInsideExon_left_bnd}" -b "${ucsc_gencode_cdsexons}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8}' | cat "${this_hdr}" - > "${gencode_bndInsideExon_left_bnd_intersect}"
echo ''
rm $gencode_bndInsideExon_left_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It is needed for the duplicate exons, otherwise same gene will appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${gencode_bndInsideExon_left_bnd_intersect}" '>' "${gencode_bndInsideExon_left_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${gencode_bndInsideExon_left_bnd_intersect}" > "${gencode_bndInsideExon_left_bnd_intersect_collapsed}"
echo ''
rm $gencode_bndInsideExon_left_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do right BND fall inside a gencode exon?'
echo ''

gencode_bndInsideExon_right_bnd="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_right_bnd."${outfile_basename}".bed
gencode_bndInsideExon_right_bnd_intersect="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_right_bnd_intersect."${outfile_basename}".bed
gencode_bndInsideExon_right_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_right_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".gencode_bndInsideExon_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tgencode_right_BND_inside_exon\tgencode_exon_containing_right_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tgencode_right_BND_inside_exon\tgencode_exon_containing_right_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${gencode_bndInsideExon_right_bnd}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" | cat "${this_hdr}" - > "${gencode_bndInsideExon_right_bnd}" || true
echo 'bedtools intersect -wao -a' "${gencode_bndInsideExon_right_bnd}" '-b' "${ucsc_gencode_cdsexons}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8} | cat' "${this_hdr}" '- | sed s/:-1--1//g >' "${gencode_bndInsideExon_right_bnd_intersect}"
bedtools intersect -wao -a "${gencode_bndInsideExon_right_bnd}" -b "${ucsc_gencode_cdsexons}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8}' | cat "${this_hdr}" - | sed 's/:-1--1//g' > "${gencode_bndInsideExon_right_bnd_intersect}"
echo ''
rm $gencode_bndInsideExon_right_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It is needed for the duplicate exons, otherwise same gene will appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${gencode_bndInsideExon_right_bnd_intersect}" '>' "${gencode_bndInsideExon_right_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${gencode_bndInsideExon_right_bnd_intersect}" > "${gencode_bndInsideExon_right_bnd_intersect_collapsed}"
echo ''
rm $gencode_bndInsideExon_right_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: merge BND_inside_exon left BND and right BND into one file'
echo ''

echo 'awk -v left_bnd_file='"${gencode_bndInsideExon_left_bnd_intersect_collapsed}" '\'
echo '    -v right_bnd_file='"${gencode_bndInsideExon_right_bnd_intersect_collapsed}"' -f' $sw'/victorchang_scripts/combine_bnd_files_key5cols.awk \'
echo '    '"${this_input}" '>' "${outprefix_gencode_bndInsideExon}"
awk -v left_bnd_file="${gencode_bndInsideExon_left_bnd_intersect_collapsed}" \
    -v right_bnd_file="${gencode_bndInsideExon_right_bnd_intersect_collapsed}" -f $sw/victorchang_scripts/combine_bnd_files_key5cols.awk \
    "${this_input}" > "${outprefix_gencode_bndInsideExon}"
echo ''
rm $gencode_bndInsideExon_left_bnd_intersect_collapsed
rm $gencode_bndInsideExon_right_bnd_intersect_collapsed
this_input="${outprefix_gencode_bndInsideExon}"

fi



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with EHRFr99 (Ensemble Human Regulatory Features)'
echo ''

# head $refdata_EHRFr99
# chrom	start	end	feature_type
# chr2	113379801	113380200	Promoter Flanking Region
# chr18	32661402	32662400	Promoter Flanking Region
# chr3	41288801	41289200	CTCF Binding Site

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".EHRFr99_header."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l) # 17
cut1_plus_4=$((cut1+4)) # 21
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\tensemble_human_regulatory_features/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\tensemble_human_regulatory_features/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${refdata_EHRFr99}" '| cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_EHRFr99}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" | bedtools intersect -wao -a - -b "${refdata_EHRFr99}" | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_EHRFr99}"
rm "${this_hdr}"
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${outprefix_EHRFr99}" '>' "${outprefix_EHRFr99_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${outprefix_EHRFr99}" > "${outprefix_EHRFr99_collapsed}"
echo ''
this_input="${outprefix_EHRFr99_collapsed}"
if [[ "$ucsc_gencode_cdsexons" != "" ]]; then
  rm $outprefix_gencode_bndInsideExon
fi



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do left BND fall in a segmental duplication region?'
echo ''

# head $refdata_segmental_duplications
# chr1	10169	37148	chr1:180723	0	+	chr1	180723	207666	26943	1	1000	N/A	N/A	N/A	N/A	align_both/0014/both0071547	27025	30	128	26897	26628	269	164	105	0.989998884633974	0.988895903739741	0.0100683956884972	0.0100742688854825
# chr1	180723	207666	chr1:10169	0	+	chr1	10169	37148	26979	1	1000	N/A	N/A	N/A	N/A	align_both/0014/both0071547	27025	30	128	26897	26628	269	164	105	0.989998884633974	0.988895903739741	0.0100683956884972	0.0100742688854825
# chr1	88000	121417	chr1:265774	0	+	chr1	265774	297956	32182	2	1000	N/A	N/A	N/A	N/A	align_both/0014/both0071548	33449	25	1299	32150	31941	209	133	76	0.993499222395023	0.992727272727273	0.00652911487615718	0.00653207256400941

left_bnd="${tmpdir}"/temp."${sample}".left_bnd."${outfile_basename}".bed
left_bnd_intersect="${tmpdir}"/temp."${sample}".left_bnd_intersect."${outfile_basename}".bed
left_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".left_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".segdup_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tleft_BND_in_segmental_duplication_region\tsegmental_duplication_region_hit_by_left_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tleft_BND_in_segmental_duplication_region\tsegmental_duplication_region_hit_by_left_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS=""} {if (NR>1) print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${left_bnd}"
awk 'BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" | cat "${this_hdr}" - > "${left_bnd}" || true
echo 'bedtools intersect -wao -a' "${left_bnd}" '-b' "${refdata_segmental_duplications}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $5 ":" $6 "-" $7} | cat' "${this_hdr}" '- | sed s/:-1--1//g >' "${left_bnd_intersect}"
bedtools intersect -wao -a "${left_bnd}" -b "${refdata_segmental_duplications}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $5 ":" $6 "-" $7}' | cat "${this_hdr}" - | sed 's/:-1--1//g' > "${left_bnd_intersect}"
echo ''
rm $left_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${left_bnd_intersect}" '>' "${left_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${left_bnd_intersect}" > "${left_bnd_intersect_collapsed}"
echo ''
rm $left_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do right BND fall in a segmental duplication region?'
echo ''

right_bnd="${tmpdir}"/temp."${sample}".right_bnd."${outfile_basename}".bed
right_bnd_intersect="${tmpdir}"/temp."${sample}".right_bnd_intersect."${outfile_basename}".bed
right_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".right_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".segdup_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tright_BND_in_segmental_duplication_region\tsegmental_duplication_region_hit_by_right_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tright_BND_in_segmental_duplication_region\tsegmental_duplication_region_hit_by_right_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${right_bnd}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" | cat "${this_hdr}" - > "${right_bnd}" || true
echo 'bedtools intersect -wao -a' "${right_bnd}" '-b' "${refdata_segmental_duplications}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $5 ":" $6 "-" $7} | cat' "${this_hdr}" '- | sed s/:-1--1//g >' "${right_bnd_intersect}"
bedtools intersect -wao -a "${right_bnd}" -b "${refdata_segmental_duplications}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $5 ":" $6 "-" $7}' | cat "${this_hdr}" - | sed 's/:-1--1//g' > "${right_bnd_intersect}"
echo ''
rm $right_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${right_bnd_intersect}" '>' "${right_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${right_bnd_intersect}" > "${right_bnd_intersect_collapsed}"
echo ''
rm $right_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: merge segmental duplication left BND and right BND into one file'
echo ''

echo 'awk -v left_bnd_file='"${left_bnd_intersect_collapsed}" '\'
echo '    -v right_bnd_file='"${right_bnd_intersect_collapsed}"' -f' $sw'/victorchang_scripts/combine_bnd_files_key5cols.awk \'
echo '    '"${this_input}" '>' "${outprefix_segdup}"
awk -v left_bnd_file="${left_bnd_intersect_collapsed}" \
    -v right_bnd_file="${right_bnd_intersect_collapsed}" -f $sw/victorchang_scripts/combine_bnd_files_key5cols.awk \
    "${this_input}" > "${outprefix_segdup}"
echo ''
rm $left_bnd_intersect_collapsed
rm $right_bnd_intersect_collapsed
this_input="${outprefix_segdup}"
rm $outprefix_EHRFr99_collapsed



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with introns and identify clean_intron_deletions'
echo ''

# head 
# chr1    67093605        67096250        C1orf141        -
# chr1    67096322        67103236        C1orf141        -
# chr1    33540699        33541128        CSMD2   -
# chr1    33541310        33542718        CSMD2   -

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".cleanIntronDels_header."${random_number}".txt
hdr_bit1=`head -n 1 "${this_input}" | sed 's/^chrom/0chrom/'`
hdr_bit2=$(echo -e "chrom\tstart\tend\tis_clean_intron_deletion")
echo -e "${hdr_bit1}\t${hdr_bit2}\toverlap" > "${this_hdr}"

chrom1=1
pos1=2
end1=3
alt1=5
chrom1=$(( $chrom1 + 3 ))
pos1=$(( $pos1 + 3 ))
end1=$(( $end1 + 3 ))
alt1=$(( $alt1 + 3 ))
num_cols1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
num_cols2=$(head -n 1 "${ucsc_refseq_introns}" | sed -e 's/\t/\n/g' | wc -l)
overlap=$(( $num_cols1 + 3 + $num_cols2 + 1 ))
chrom2=$(( overlap - 5 ))
pos2=$(( overlap - 4 ))
end2=$(( overlap - 3 ))
gene=$(( overlap - 2 ))

this_input_basename=$(basename $this_input)
tmpfile_cleanIntronDels="${tmpdir}"/"${this_input_basename}"
echo 'Create a temporary file of only the clean_intron_deletions for this sample.'
echo 'grep -v' 'SVTYPE=BND' "${this_input}" '| awk' 'BEGIN {FS="\t";OFS="\t"} {print $1, $2, $end1, $0}' '| awk' 'BEGIN {FS="\t";OFS="\t"} {if ($2 != $3) {print $0}}' '| bedtools intersect -a - -b' $ucsc_refseq_introns '-wao | \'
echo '  awk -v overlap='"$overlap" '-v pos1='"$pos1" '-v end1='"$end1" '-v alt1='"$alt1" '-v gene='"$gene" '-v chrom2='"$chrom2" '-v pos2='"$pos2" '-v end2='"$end2" 'function abs(v) {return v < 0 ? -v : v} BEGIN {FS="\t";OFS="\t"} {if ( ($overlap>0) && (abs($overlap-($end2-$pos2))<=5) && (abs(($end1-$pos1)-($end2-$pos2))<=5) && (abs($pos1-$pos2)<=5) && (abs($end1-$end2)<=5) ) {print $1, $pos1, $end1, "is_clean_intron"}}' '>' $tmpfile_cleanIntronDels
grep -v 'SVTYPE=BND' "${this_input}" | awk 'BEGIN {FS="\t";OFS="\t"} {print $1, $2, $3, $0}' | awk 'BEGIN {FS="\t";OFS="\t"} {if ($2 != $3) {print $0}}' | bedtools intersect -a - -b $ucsc_refseq_introns -wao | \
  awk -v overlap="$overlap" -v pos1="$pos1" -v end1="$end1" -v alt1="$alt1" -v gene="$gene" -v chrom2="$chrom2" -v pos2="$pos2" -v end2="$end2" 'function abs(v) {return v < 0 ? -v : v} BEGIN {FS="\t";OFS="\t"} {if ( ($overlap>0) && (abs($overlap-($end2-$pos2))<=5) && (abs(($end1-$pos1)-($end2-$pos2))<=5) && (abs($pos1-$pos2)<=5) && (abs($end1-$end2)<=5) ) {print $1, $pos1, $end1, "is_clean_intron"}}' > $tmpfile_cleanIntronDels || true
echo ''

num_clean_intron_deletions=$(wc -l "${tmpfile_cleanIntronDels}" | cut -d" " -f1)

if [[ $num_clean_intron_deletions -eq 0 ]]; then
  echo 'There are no clean_intron_deletions for this sample.'
  echo -e "${hdr_bit1}\tis_clean_intron_deletion" > "${this_hdr}"
  echo 'grep -v' '^chrom' "${this_input}" '| awk' 'BEGIN {FS="\t";OFS="\t"} {print $0, "."}' '| cat' "${this_hdr}" '- | sed s/^0chrom/chrom/ >' "${outprefix_cleanIntronDels}"
  grep -v '^chrom' "${this_input}" | awk 'BEGIN {FS="\t";OFS="\t"} {print $0, "."}' | cat "${this_hdr}" - | sed 's/^0chrom/chrom/' > "${outprefix_cleanIntronDels}" || true
else
  echo 'Use the temporary file of only clean_intron_deletions as a reference for this sample.'
  echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${tmpfile_cleanIntronDels}" '| grep -v ^chrom | cat' "${this_hdr}" '- | cut -d$\t -f"1-$cut1,$cut1_plus_4" | \'
  echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_cleanIntronDels}"
  bedtools intersect -wao -a "${this_input}" -b "${tmpfile_cleanIntronDels}" | grep -v '^chrom' | cat "${this_hdr}" - | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" | \
    sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_cleanIntronDels}" || true
  echo ''
fi
rm $outprefix_segdup

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${outprefix_cleanIntronDels}" '>' "${outprefix_cleanIntronDels_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${outprefix_cleanIntronDels}" > "${outprefix_cleanIntronDels_collapsed}"
echo ''
echo 'Output:' $outprefix_cleanIntronDels_collapsed
echo ''
rm $outprefix_cleanIntronDels

rm -rf $outprefix_sorted
rm -rf $outprefix_CDSexons
rm -rf $outprefix_CDSexons_collapsed
rm -rf $outprefix_genes
rm -rf $outprefix_genes_collapsed
rm -rf $outprefix_bndInsideExon
rm -rf $outprefix_gencode_CDSexons
rm -rf $outprefix_gencode_CDSexons_collapsed
rm -rf $outprefix_gencode_genes
rm -rf $outprefix_gencode_genes_collapsed
rm -rf $outprefix_gencode_bndInsideExon
rm -rf $outprefix_EHRFr99
rm -rf $outprefix_EHRFr99_collapsed
rm -rf $outprefix_segdup
rm -rf $outprefix_cleanIntronDels

echo ''
echo 'Finished!'
echo ''

