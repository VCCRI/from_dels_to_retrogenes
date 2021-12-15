#!/bin/bash

genome=$1 # hg19 or hg38
cohort=$2 # will appear in output data in cohort field
sample=$3 # will appear in output data in sample field and will appear in output file names
outdir=$4 # output files will be written here
tmpdir=$5 # temporary files will be written here, some will be deleted, user can delete this directory after running the pipeline
infile=$6
outfile=$7

currdir=$(pwd)

if [[ $genome == "hg19" ]]; then
  introns="${currdir}"/../reference_data/hg19_UCSC_GRCh37_GenesAndGenePredictions_introns_RefSeq_and_Gencode_20211116.bed
else
  introns="${currdir}"/../reference_data/hg38_UCSC_GRCh38_GenesAndGenePredictions_introns_RefSeq_and_Gencode_20211116.bed
fi

echo 'infile:' $infile
echo ''

tmpfile="${tmpdir}"/"${cohort}"."${sample}".tmpfile.txt

echo 'awk print $1, $2, $3, vaf, depth1, depth2, depth3, depth4, depth5, $0, "."' $infile '\'
echo 'sort -k1,1 -k2,2n -k3,3n | \'
echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if ($2 != $3) {print $0}}' '| \'
echo 'bedtools intersect -a - -b' $introns '-wao >' $tmpfile
echo ''

awk '
BEGIN {FS="\t";OFS="\t"} 
{
  if (NR==1) {
    vaf_i=-1
    svtype_i=-1
    depth1_i=-1
    depth2_i=-1
    depth3_i=-1
    depth4_i=-1
    depth5_i=-1
    for (i=1; i<=NF; ++i) {
      if ($i=="VAF") { vaf_i=i }
      if (length($i)>=4) { if (substr($i,1,4)=="VAF_") { vaf_i=i } }
      if ($i=="SVTYPE") { svtype_i=i }
      if ($i=="ASRP") { depth1_i=i }
      if ($i=="gridssASRP") { depth1_i=i }
      if ($i=="REF") { if (i>5) { depth2_i=i } }
      if ($i=="gridssREF") { depth2_i=i }
      if ($i=="REFPAIR") { depth3_i=i }
      if ($i=="gridssREFPAIR") { depth3_i=i }
      if ($i=="RP") { depth4_i=i }
      if ($i=="gridssRP") { depth4_i=i }
      if ($i=="SR") { depth5_i=i }
      if ($i=="gridssSR") { depth5_i=i }
    }
  } else {
    if (svtype_i != -1) {
      if (($svtype_i=="DEL") || ($svtype_i=="INDEL")) {
        vaf="."
        depth1="."
        depth2="."
        depth3="."
        depth4="."
        depth5="."
        if (vaf_i != -1) { vaf=$vaf_i } 
        if (depth1_i != -1) { depth1=$depth1_i } 
        if (depth2_i != -1) { depth2=$depth2_i } 
        if (depth3_i != -1) { depth3=$depth3_i } 
        if (depth4_i != -1) { depth4=$depth4_i } 
        if (depth5_i != -1) { depth5=$depth5_i } 
        print $1, $2, $3, vaf, depth1, depth2, depth3, depth4, depth5, $0, "."
      }
    }
  }
}' $infile | \
sort -k1,1 -k2,2n -k3,3n | \
awk 'BEGIN {FS="\t";OFS="\t"} {if ($2 != $3) {print $0}}' | \
bedtools intersect -a - -b $introns -wao > $tmpfile
echo ''

overlap=$(head -n 1 "${tmpfile}" | sed -e 's/\t/\n/g' | wc -l) # 85

echo 'awk -v cohort='"$cohort" '-v sample='"$sample" '-v overlap='"$overlap"
echo 'print cohort, sample, $chrom1, $pos1, $end1, $alt1, $gene, $chrom2, $pos2, $end2, $strand, $vaf, $depth1, $depth2, $depth3, $depth4, $depth5'
echo '}' $tmpfile '| uniq >' $outfile
echo ''

awk -v cohort="$cohort" -v sample="$sample" -v overlap="$overlap" '
function abs(v) {return v < 0 ? -v : v}
BEGIN {
  FS="\t";OFS="\t"
  chrom2=overlap-5
  pos2=overlap-4
  end2=overlap-3
  gene=overlap-2
  strand=overlap-1
  chrom1=1
  pos1=2
  end1=3
  vaf=4
  depth1=5
  depth2=6
  depth3=7
  depth4=8
  depth5=9
  alt1=14  
}
{
  if ( ($overlap>0) && (abs($overlap-($end2-$pos2))<=5) && (abs(($end1-$pos1)-($end2-$pos2))<=5) && (abs($pos1-$pos2)<=5) && (abs($end1-$end2)<=5) ) {
    print cohort, sample, $chrom1, $pos1, $end1, $alt1, $gene, $chrom2, $pos2, $end2, $strand, $vaf, $depth1, $depth2, $depth3, $depth4, $depth5
  }
}' $tmpfile | uniq > $outfile
echo ''

#rm $tmpfile

echo 'Finished!'
echo ''


