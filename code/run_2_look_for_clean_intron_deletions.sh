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
  introns="${currdir}"/../reference_data/hg19_UCSC_GRCh37_GenesAndGenePredictions_introns_RefSeq_20200324.txt
else
  introns="${currdir}"/../reference_data/hg38_UCSC_GRCh38_GenesAndGenePredictions_introns_RefSeq_20200901.txt
fi

echo 'infile:' $infile
echo ''

tmpfile="${tmpdir}"/"${cohort}"."${sample}".tmpfile.txt

vaf=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n -i '^VAF$' | cut -d':' -f1`
svlen=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n -i '^SVLEN$' | cut -d':' -f1`
svtype=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n -i '^SVTYPE$' | cut -d':' -f1`
if [[ $vaf == "" ]]; then
  vaf=0
fi
if [[ $svlen == "" ]]; then
  svlen=0
fi
echo 'vaf:' $vaf 'svlen:' $svlen 'svtype:' $svtype
echo ''

echo 'awk -v vaf='"$vaf" '-v svlen='"$svlen" '-v svtype='"$svtype" 'BEGIN {FS="\t";OFS="\t"} {if ((NR>1) && (($svtype=="DEL")||($svtype=="INDEL"))) {print $1, $2, $3, $vaf, $0, "."}}' $infile '| sort -k1,1 -k2,2n -k3,3n | awk' 'BEGIN {FS="\t";OFS="\t"} {if ($2 != $3) {print $0}}' '| bedtools intersect -a - -b' $introns '-wao >' $tmpfile

awk -v vaf="$vaf" -v svlen="$svlen" -v svtype="$svtype" 'BEGIN {FS="\t";OFS="\t"} {if ((NR>1) && (($svtype=="DEL")||($svtype=="INDEL"))) {print $1, $2, $3, $vaf, $0, "."}}' $infile | sort -k1,1 -k2,2n -k3,3n | awk 'BEGIN {FS="\t";OFS="\t"} {if ($2 != $3) {print $0}}' | bedtools intersect -a - -b $introns -wao > $tmpfile
echo ''

overlap=$(head -n 1 "${tmpfile}" | sed -e 's/\t/\n/g' | wc -l) # 85
chrom2=$(( overlap - 5 ))
pos2=$(( overlap - 4 ))
end2=$(( overlap - 3 ))
gene=$(( overlap - 2 ))
strand=$(( overlap - 1 ))
chrom1=1
pos1=2
end1=3
vaf=4
alt1=9

echo 'awk -v cohort='"$cohort" '-v sample='"$sample" '-v overlap='"$overlap" '-v chrom1='"$chrom1" '-v pos1='"$pos1" '-v end1='"$end1" '-v alt1='"$alt1" '-v gene='"$gene" '-v chrom2='"$chrom2" '-v pos2='"$pos2" '-v end2='"$end2" '-v strand='"$strand" '-v vaf='"$vaf" 'function abs(v) {return v < 0 ? -v : v} BEGIN {FS="\t";OFS="\t"} {if ( ($overlap>0) && (abs($overlap-($end2-$pos2))<=5) && (abs(($end1-$pos1)-($end2-$pos2))<=5) && (abs($pos1-$pos2)<=5) && (abs($end1-$end2)<=5) ) {print cohort, sample, '$chrom1', '$pos1', '$end1', '$alt1', '$gene', '$chrom2', '$pos2', '$end2', '$strand', '$vaf'}}' $tmpfile '| uniq >' $outfile

awk -v cohort="$cohort" -v sample="$sample" -v overlap="$overlap" -v chrom1="$chrom1" -v pos1="$pos1" -v end1="$end1" -v alt1="$alt1" -v gene="$gene" -v chrom2="$chrom2" -v pos2="$pos2" -v end2="$end2" -v strand="$strand" -v vaf="$vaf" 'function abs(v) {return v < 0 ? -v : v} BEGIN {FS="\t";OFS="\t"} {if ( ($overlap>0) && (abs($overlap-($end2-$pos2))<=5) && (abs(($end1-$pos1)-($end2-$pos2))<=5) && (abs($pos1-$pos2)<=5) && (abs($end1-$end2)<=5) ) {print cohort, sample, $chrom1, $pos1, $end1, $alt1, $gene, $chrom2, $pos2, $end2, $strand, $vaf}}' $tmpfile | uniq > $outfile
echo ''

rm $tmpfile

echo 'Finished!'
echo ''


