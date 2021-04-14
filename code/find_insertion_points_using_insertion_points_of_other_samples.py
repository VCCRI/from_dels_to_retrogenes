# python3 find_insertion_points_using_insertion_points_of_other_samples.py --in_retrocopies temp4.txt --in_sv A0813K1N62883440T1D1.gridss_2_7_3.vcf.gz --in_retrocopies_other_samples TGA__insertion_points_for_deletions_that_are_retrocopied_genes.insertion_pt_flags.sort_by_sample.unique.txt -o temp4_out.txt

# in_retrocopies format.
# Contains the left-most and right-most clean-intron-deletion coordinates, the point before the left-most and the point after the right-most coordinates that connect to the insertion-point (called retrocopy_start and retrocopy_end), and the insertion-points left and right positions (usually approx 11 bp apart).
# Insertion point can be on a different chromosome. The clean-intron-deletions and retrocopy start/end are on the same chromosome as each other.
#
# cohort	sample	gene	chrom	clean_intron_del_start	clean_intron_del_end	retrocopy_start	retrocopy_end	insertion_chrom	left_insertion_point	right_insertion_point	insertion_sequence	retrocopy_insertion_direction	is_left_BND_present	is_right_BND_present
# mycohort	mysample	RBP5	12	7281658	7281688						
# mycohort	mysample	TYRO3	15	41866014	41870084	41851396	41871536	13	44069827	44069838	AAAAAAAAAAAAAAAAAAAAAAAAA	FORWARD	left_BND_is_present	right_BND_is_present
# mycohort	mysample	SMAD4	18	48556993	48604626	48556624	48606219	9	127732636	127732713	FORWARD	left_BND_is_present	right_BND_is_present	
# mycohort	mysample	PRKRA	2	179296982	179315693						

# in_retrocopies_other_samples format, it is the same as the in_retrocopies format:
#
# cohort	sample	gene	chrom	clean_intron_del_start	clean_intron_del_end	retrocopy_start	retrocopy_end	insertion_chrom	left_insertion_point	right_insertion_point	insertion_sequence	retrocopy_insertion_direction	is_left_BND_present	is_right_BND_present
# mycohort	sample111	PER3	1	7890012	7890066								
# mycohort	sample111	RBP5	12	7281658	7281688								
# mycohort	sample311	RBP5	12	7281658	7281688								
# mycohort	sample311	INTS10	8	19675175	19679949	19674931	19680587	3	129510413	129510428	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	FORWARD	left_BND_is_present	right_BND_is_present
# mycohort	sample211	SKA3	13	21729289	21746478	21727743	21750661	11	108585748	108585765	TGTTTTTTTTTTTTTTTTTTTTTTAGTTTTTA	FORWARD	left_BND_is_present	right_BND_is_present
# mycohort	sample111	SKA3	13	21729289	21746478	21727734	21750661	11	108585748	108585765	TGTTTTTTTTTTTTTTTTTTTTT	FORWARD	left_BND_is_present	right_BND_is_present
# mycohort	sampleT1D1	CBX3	7	26241446	26242590	26241365	26252971	15	40854180	40854194	AAAAAAAA	FORWARD	left_BND_is_present	right_BND_is_present
# mycohort	sampleT4D1	RBP5	12	7281658	7281688								


# The following can convert find_insertion_points_using_insertion_points_of_other_samples.py output to this find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py input:
# awk 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {direction="retrocopy_insertion_direction";left="is_left_BND_present";right="is_right_BND_present"} else { left="";right="";direction=""; if ($9!=""){left="left_BND_is_present";direction="FORWARD"}; if ($10!=""){right="right_BND_is_present";direction="FORWARD"} } printf $1; for (i=2; i<=NF; ++i) {printf OFS $i}; printf OFS direction OFS left OFS right RS }' TGA__insertion_points_for_deletions_that_are_retrocopied_genes.sort_by_sample.txt | tr -d $'\r' > TGA__insertion_points_for_deletions_that_are_retrocopied_genes.insertion_pt_flags.sort_by_sample.txt

# Have one copy of each insertion point, even if there are multiple samples that have the same insertion point.
# cohort	sample	gene	chrom	clean_intron_del_start	clean_intron_del_end	retrocopy_start	retrocopy_end	insertion_chrom	left_insertion_point	right_insertion_point	insertion_sequence	retrocopy_insertion_direction	is_left_BND_present	is_right_BND_present
# head -n 1 TGA__insertion_points_for_deletions_that_are_retrocopied_genes.insertion_pt_flags.sort_by_sample.txt > temp_hdr.txt
# awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print "cohort", "sample", $3, $4, $5, $6, $7, $8, $9, $10, $11, "", $13, $14, $15}}' TGA__insertion_points_for_deletions_that_are_retrocopied_genes.insertion_pt_flags.sort_by_sample.txt | sort | uniq | cat temp_hdr.txt - > TGA__insertion_points_for_deletions_that_are_retrocopied_genes.insertion_pt_flags.sort_by_sample.unique.txt

# Program specification:
#
# For each clean_intron_deletion gene having blank insertion_point
#   To identify an insertion point, usually need two different BNDs, inserted from points are on either side of the clean_intron_deletion, and inserted to points are approx 11 bp apart.
#   Get the insertion_point ranges of other samples having same clean_intron_deletion gene. both where inserted to and inserted from.
#   If there is at least one BND connecting the inserted_from and inserted_to ranges then
#     Call this insertion point, say whether the left BND or right BND were seen (ie, inserted_from is before or after the clean_intron_deletion gene).

import sys
import csv
import subprocess
import vcf # pip3 install --install-option="--prefix=/g/data/jb96/software/python_packages" pyvcf
import argparse


# cd /tmp
# virtualenv test_env
# source test_env/bin/activate
#
# python3
# help('modules')

def parse_arguments(args):

    parser = argparse.ArgumentParser(description='Read in retrocopied genes and any insertion points already found for them. An insertion point consists of two points - the left and right sides of the insertion region, usually separated by approx 11 bp. For those that do not have an insertion point found yet, look at the retrocopied genes and insertion points of other samples. If another sample has the same retrocopied gene and the 2 sides of the insertion point, and if this sample that is missing an insertion point does have one of the sides (left or right) insertion point connections as a BND record, then we assume that it has the same insertion point, even though we have not seen the second side of the insertion point.')
    parser.add_argument('-i', '--in_retrocopies', dest='in_retrocopies_file', 
                        help='Tab-delimited list of retrocopied genes for which this program will try to fill in any missing insertion points, in format: gene   chrom   clean_intron_del_start   clean_intron_del_end   retrocopy_start   retrocopy_end   insertion_chrom   left_insertion_point   right_insertion_point   insertion_sequence   is_left_BND_present   is_right_BND_present')
    parser.add_argument('--in_sv', dest='in_sv_file', 
                        help='Input VCF file of structural variants containing BND records that represent potential retrocopied gene insertion points. Must be bgzip compressed and indexed with tabix so this program can retrieve by co-ordinates.')
    parser.add_argument('--in_retrocopies_other_samples', dest='in_retrocopies_other_samples_file', 
                        help='Tab-delimited list of retrocopied genes of other samples that this program will consult to look for already-existing already-known insertion points')
    parser.add_argument('-o', '--output', dest='output_file', 
                        help='Output file, one line per input clean_intron_deletion, specifying the insertion points. The multiple records of the same gene will have the same insertion points.')

    args = parser.parse_args()

    return args


def is_integer(s):
  try:
    int(s)
    return True
  except ValueError:
    return False


def do_regions_overlap( chrom1, left1, right1, chrom2, left2, right2 ):

  result = False

  chrom1 = str(chrom1)
  chrom2 = str(chrom2)
  left1 = int(left1)
  right1 = int(right1)
  left2 = int(left2)
  right2 = int(right2)
  if (left1 > right1):
    temp = left1
    left1 = right1
    right1 = temp
  if (left2 > right2):
    temp = left2
    left2 = right2
    right2 = temp

  if (chrom1 == chrom2):
    if ((left1 <= left2) and (right1 >= left2)):
      result = True
    if ((left1 <= right2) and (right1 >= right2)):
      result = True
    if ((left1 >= left2) and (right1 >= right2)):
      result = True

  return result


def is_it_same_point( chrom1, pos1, chrom2, pos2 ):
  it_is_same_point = False
  if (chrom1 == chrom2):
    pos1 = int(pos1)
    pos2 = int(pos2)
    diff = abs(pos1 - pos2)
    if (diff <= 10):
      it_is_same_point = True
  return it_is_same_point


def is_point2_before_point1( chrom1, pos1, chrom2, pos2 ):
  yes_it_is = False
  if (chrom1 == chrom2):
    pos1 = int(pos1)
    pos2 = int(pos2)
    diff = abs(pos1 - pos2)
    if (pos2 <= pos1):
      yes_it_s = True
  return yes_it_is


def is_point2_after_point1( chrom1, pos1, chrom2, pos2 ):
  yes_it_is = False
  if (chrom1 == chrom2):
    pos1 = int(pos1)
    pos2 = int(pos2)
    diff = abs(pos1 - pos2)
    if (pos2 >= pos1):
      yes_it_is = True
  return yes_it_is


class InsertionPoint:
    """
    This object represents the insertion point of a retrocopied gene.
    """
    def __init__(self, initial_retrocopy_chrom, initial_retrocopy_start, initial_retrocopy_end, initial_insertion_chrom, initial_left_insertion_point, initial_right_insertion_point, initial_insertion_sequence, initial_retrocopy_insertion_direction, initial_is_left_BND_present, initial_is_right_BND_present):
        self.retrocopy_chrom = initial_retrocopy_chrom
        self.retrocopy_start = 0
        if is_integer(initial_retrocopy_start):
          self.retrocopy_start = int(initial_retrocopy_start)
        self.retrocopy_end = 0
        if is_integer(initial_retrocopy_end):
          self.retrocopy_end = int(initial_retrocopy_end)
        self.insertion_chrom = initial_insertion_chrom
        self.left_insertion_point = int(initial_left_insertion_point)
        self.right_insertion_point = int(initial_right_insertion_point)
        self.insertion_sequence = initial_insertion_sequence
        self.retrocopy_insertion_direction = initial_retrocopy_insertion_direction
        self.is_left_BND_present = initial_is_left_BND_present
        self.is_right_BND_present = initial_is_right_BND_present
    def print(self, **kwargs):
        print( "retrocopy: " + self.retrocopy_chrom + ":" + str(self.retrocopy_start) + "-" + str(self.retrocopy_end) + " insertion: " + self.insertion_chrom + ":" + str(self.left_insertion_point) + "-" + str(self.right_insertion_point) + " seq: " + str(self.insertion_sequence) + " direction: " + str(self.retrocopy_insertion_direction) + " is_left_BND_present: " + str(self.is_left_BND_present) + " is_right_BND_present: " + str(self.is_right_BND_present))
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)


def read_retrocopies_with_insertion_points( input_file_path ):

  # Read through all the retrocopied genes and their insertion points for many samples.
  # Keep only those that have insertion points. 
  # We will be using those insertion points to see whether the sample that we are processing has an incomplete insertion point
  # at the same positions as these insertion points of other samples, for the same gene.
  # There are many old retrocopied genes and thus many samples have the exact same retrocopied gene and insertion points.

  retrocopies = {}

  with open(input_file_path, newline='') as csvfile:
    retrocopies_data = list(csv.reader(csvfile, delimiter='\t'))

  is_header = True
  for data in retrocopies_data:
    cohort = str(data[0])
    sample = str(data[1])
    gene = str(data[2])
    chrom = str(data[3])
    clean_intron_del_start = str(data[4])
    clean_intron_del_end = str(data[5])
    retrocopy_start = str(data[6])
    retrocopy_end = str(data[7])
    insertion_chrom = str(data[8])
    left_insertion_point = str(data[9])
    right_insertion_point = str(data[10])
    insertion_sequence = str(data[11])
    retrocopy_insertion_direction = str(data[12])
    is_left_BND_present = str(data[13])
    is_right_BND_present = str(data[14])

    keep_this_record = False
    if (is_header):
      keep_this_record = False
      is_header = False
    else:
      if ((is_left_BND_present != ".") and (is_left_BND_present != "") and (is_right_BND_present != ".") and (is_right_BND_present != "")):
        keep_this_record = True

    # This retrocopied gene has a complete insertion point.
    # Let's save it for comparisons with the main sample we are processing.
    # Ignore the header.

    if (keep_this_record):

      left_insertion_point = int(left_insertion_point)
      right_insertion_point = int(right_insertion_point)

      this_insertion_point_is_already_in_list = False

      if gene in retrocopies:
        retrocopies_for_this_gene = retrocopies[gene]
        for this_retrocopy_insertion_point in retrocopies_for_this_gene:
          this_chrom = this_retrocopy_insertion_point.insertion_chrom
          if (chrom == this_chrom):
            this_left_insertion_point = this_retrocopy_insertion_point.left_insertion_point
            this_right_insertion_point = this_retrocopy_insertion_point.right_insertion_point
            overlap_result = do_regions_overlap( chrom, left_insertion_point, right_insertion_point, this_chrom, this_left_insertion_point, this_right_insertion_point )
            if (overlap_result == True):
              this_insertion_point_is_already_in_list = True # This insertion point is already saved in the comparison list. Must have previously seen it for a different sample.

      if (this_insertion_point_is_already_in_list == False):

        new_retrocopy_insertion_point = InsertionPoint( chrom, retrocopy_start, retrocopy_end, insertion_chrom, left_insertion_point, right_insertion_point, insertion_sequence, retrocopy_insertion_direction, is_left_BND_present, is_right_BND_present )

        if gene in retrocopies:
          retrocopies_for_this_gene = retrocopies[gene]
          retrocopies_for_this_gene.append( new_retrocopy_insertion_point )
        else:
          retrocopies_for_this_gene = [ new_retrocopy_insertion_point ]
        retrocopies[gene] = retrocopies_for_this_gene

  #for gene in retrocopies:
  #  print(gene)
  #  for gene_insertion_point in retrocopies[gene]:
  #    print(gene_insertion_point)

  return retrocopies


def does_sample_have_this_insertion_point( vcf_reader, clean_intron_del_chrom, clean_intron_del_start, clean_intron_del_end, one_insertion_point ):

  sample_does_have_this_insertion_point = False
  new_insertion_point = ""

  insertion_chrom = one_insertion_point.insertion_chrom
  left_insertion_point = one_insertion_point.left_insertion_point - 1 - 10
  right_insertion_point = one_insertion_point.right_insertion_point + 10

  # hg19 contig hs37d5 is making this program crash with:
  #
  #  File "/g/data/jb96/emmrat/clean_intron_deletions/code/find_insertion_points_using_insertion_points_of_other_samples.py", line 259, in does_sample_have_this_insertion_point
  #    for record in vcf_reader.fetch( insertion_chrom, left_insertion_point, right_insertion_point ):
  #  File "/g/data/jb96/software/python_packages/lib/python3.7/site-packages/vcf/parser.py", line 631, in fetch
  #    self.reader = self._tabix.fetch(chrom, start, end)
  #  File "pysam/libctabix.pyx", line 509, in pysam.libctabix.TabixFile.fetch
  #ValueError: could not create iterator for region 'hs37d5:7314085-7314154'
  #
  # It occurs when there are no entries at all in the vcf file for hs37d5 even though this contig is defined in the header of the vcf file.

  interrogate_this_position = True
  try:
    vcf_reader.fetch( insertion_chrom, left_insertion_point, right_insertion_point )
  except:
    interrogate_this_position = False

  if (interrogate_this_position):
    #print('tabix ' + str(in_sv_path) + ' record.var_subtype=' + str(insertion_chrom) + ':' + str(left_insertion_point) + '-' + str(right_insertion_point))
    for record in vcf_reader.fetch( insertion_chrom, left_insertion_point, right_insertion_point ):

     #print('record.var_type=' + str(record.var_type) + ' record.var_subtype=' + str(record.var_subtype) + ' record.ALT=' + str(record.ALT))
     if ((record.var_type=='sv') and (record.var_subtype=='complex')):
       for this_alt in record.ALT:
         if (this_alt.type == 'BND'):
           bnd_chrom = record.CHROM
           bnd_pos = record.POS
           bnd_alt_chrom = this_alt.chr
           bnd_alt_pos = this_alt.pos
           bnd_connectingSequence = this_alt.connectingSequence
           #print('found: ' + str(bnd_chrom) + ':' + str(bnd_pos) + ' ' + str(bnd_alt_chrom) + ':' + str(bnd_alt_pos) )

           it_is_a_candidate_point_1 = is_point2_before_point1( clean_intron_del_chrom, clean_intron_del_start, bnd_alt_chrom, bnd_alt_pos )
           it_is_a_candidate_point_2 = is_point2_after_point1( clean_intron_del_chrom, clean_intron_del_end, bnd_alt_chrom, bnd_alt_pos )

           if ((it_is_a_candidate_point_1) or (it_is_a_candidate_point_2)):

             sample_does_have_this_insertion_point = True

             new_retrocopy_chrom = one_insertion_point.retrocopy_chrom # which is the same as clean_intron_del_chrom
             new_retrocopy_start = one_insertion_point.retrocopy_start
             new_retrocopy_end = one_insertion_point.retrocopy_end
             new_insertion_chrom = one_insertion_point.insertion_chrom
             new_left_insertion_point = one_insertion_point.left_insertion_point
             new_right_insertion_point = one_insertion_point.left_insertion_point
             new_insertion_sequence = bnd_connectingSequence
             new_retrocopy_insertion_direction = one_insertion_point.retrocopy_insertion_direction
             new_is_left_BND_present = ""
             new_is_right_BND_present = ""

             if (it_is_a_candidate_point_1):
               new_is_left_BND_present = "left_BND_is_present"
               new_retrocopy_start = bnd_alt_pos
             else:
               new_is_right_BND_present = "right_BND_is_present"
               new_retrocopy_end = bnd_alt_pos

             diff1 = abs(int(one_insertion_point.left_insertion_point) - int(bnd_pos))
             diff2 = abs(int(one_insertion_point.right_insertion_point) - int(bnd_pos))
             if (diff1 <= diff2):
               new_left_insertion_point = bnd_pos
               if (new_is_left_BND_present == "left_BND_is_present"):
                 new_retrocopy_insertion_direction = "FORWARD"
               else:
                 new_retrocopy_insertion_direction = "REVERSE"
             else:
               new_right_insertion_point = bnd_pos
               if (new_is_right_BND_present == "right_BND_is_present"):
                 new_retrocopy_insertion_direction = "FORWARD"
               else:
                 new_retrocopy_insertion_direction = "REVERSE"

             new_insertion_point = InsertionPoint(new_retrocopy_chrom, new_retrocopy_start, new_retrocopy_end, new_insertion_chrom, new_left_insertion_point, new_right_insertion_point, new_insertion_sequence, new_retrocopy_insertion_direction, new_is_left_BND_present, new_is_right_BND_present)

  return sample_does_have_this_insertion_point, new_insertion_point


def find_insertion_points_using_insertion_points_of_other_samples(in_retrocopies_path, in_sv_path, in_retrocopies_other_samples_path, output_file_path):

  writer = csv.writer(open(output_file_path, 'w'), delimiter='\t')

  # Read in all the retrocopies and their insertion points of other samples,
  # so that insertion points are ready to be used by the sample being processed by this program.
  other_samples_retrocopies_all_genes = read_retrocopies_with_insertion_points( in_retrocopies_other_samples_path )

  # Open the tabix indexed VCF file of BND records, ready to look for the insertion point BND records for each clean_intron_deletion.
  sv_vcf_reader = vcf.Reader(filename=in_sv_path)

  # Read in all retrocopied genes of the sample being processed.
  with open(in_retrocopies_path, newline='') as csvfile:
    retrocopies_data = list(csv.reader(csvfile, delimiter='\t'))
    #print(dir(retrocopied_genes_data))

  # Process each retrocopied gene. For those that don't have insertion points yet, try to find the 2 insertion points.
  #for this_retrocopied_gene in list_of_retrocopied_genes:

  is_header = True
  for data in retrocopies_data:
    cohort = str(data[0])
    sample = str(data[1])
    gene = str(data[2])
    chrom = str(data[3])
    clean_intron_del_start = str(data[4])
    clean_intron_del_end = str(data[5])
    retrocopy_start = str(data[6])
    retrocopy_end = str(data[7])
    insertion_chrom = str(data[8])
    left_insertion_point = str(data[9])
    right_insertion_point = str(data[10])
    insertion_sequence = str(data[11])
    retrocopy_insertion_direction = str(data[12])
    is_left_BND_present = str(data[13])
    is_right_BND_present = str(data[14])

    this_record_needs_processing = False
    if (is_header):
      this_record_needs_processing = False
      is_header = False
    else:
      if ((is_left_BND_present == ".") or (is_left_BND_present == "")):
        this_record_needs_processing = True
      if ((is_right_BND_present == ".") or (is_right_BND_present == "")):
        this_record_needs_processing = True

    # If this record is the header or if it already has its two insertion points, then simply copy this record to output.
    if (this_record_needs_processing == False):
      writer.writerow(data)

    # For the records that don't have two insertion points, process it.
    # Process it to see if it has one insertion point of another sample.
    # If it does, then take that other sample's two insertion points.
    else:

      found_insertion_point = False

      if gene in other_samples_retrocopies_all_genes:

        other_samples_retrocopies = other_samples_retrocopies_all_genes[gene]

        keep_looking_for_insertion_points = True
        i = 0
        if (i >= len(other_samples_retrocopies)):
          keep_looking_for_insertion_points = False

        while (keep_looking_for_insertion_points):

          data_of_another_sample = other_samples_retrocopies[i]

          sample_does_have_this_insertion_point, insertion_point = does_sample_have_this_insertion_point( sv_vcf_reader, chrom, clean_intron_del_start, clean_intron_del_end, data_of_another_sample )

          # This retrocopied gene doesn't have 2 insertion points, but it does have one insertion point found in another sample.
          # So let's assign that other sample's insertion points to this sample.
          if (sample_does_have_this_insertion_point):

            retrocopy_chrom = insertion_point.retrocopy_chrom
            retrocopy_start = insertion_point.retrocopy_start
            retrocopy_end = insertion_point.retrocopy_end
            insertion_chrom = insertion_point.insertion_chrom
            left_insertion_point = insertion_point.left_insertion_point
            right_insertion_point = insertion_point.right_insertion_point
            insertion_sequence = insertion_point.insertion_sequence
            retrocopy_insertion_direction = insertion_point.retrocopy_insertion_direction
            is_left_BND_present = insertion_point.is_left_BND_present
            is_right_BND_present = insertion_point.is_right_BND_present

            # We found a new insertion point, so write out this retrocopied gene giving it this new insertion point.
            outline = [cohort, sample, gene, chrom, clean_intron_del_start, clean_intron_del_end, retrocopy_start, retrocopy_end, insertion_chrom, left_insertion_point, right_insertion_point, insertion_sequence, retrocopy_insertion_direction, is_left_BND_present, is_right_BND_present]
            writer.writerow(outline)

            found_insertion_point = True
            keep_looking_for_insertion_points = False

          i = i + 1
          if (i >= len(other_samples_retrocopies)):
            keep_looking_for_insertion_points = False

      # If we didn't find a new insertion point, then write out this retrocopied gene as is.
      if (found_insertion_point == False):
        writer.writerow(data)

  return


def main(args):

  args = parse_arguments(args)
  in_retrocopies_path = args.in_retrocopies_file
  in_sv_path = args.in_sv_file
  in_retrocopies_other_samples_path = args.in_retrocopies_other_samples_file
  output_file_path = args.output_file
  print(' ')

  find_insertion_points_using_insertion_points_of_other_samples(in_retrocopies_path, in_sv_path, in_retrocopies_other_samples_path, output_file_path)

  print(' ')
  print('find_insertion_points_using_insertion_points_of_other_samples.py')


if __name__ == "__main__":
    main(sys.argv)



