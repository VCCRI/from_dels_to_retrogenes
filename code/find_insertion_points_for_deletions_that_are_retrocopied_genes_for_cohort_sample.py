# python3 find_insertion_points_for_deletions_that_are_retrocopied_genes_for_cohort_sample.py --in_dels temp1_SCHN_1891.txt --in_sv SCHN_1891.gridss_2_7_3.vcf.gz --gene_regions /home/emma/emma/clean_intron_deletions/UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq_20200324.txt --gene_region_extension_for_start_of_gene 20 --gene_region_extension 100000 --blacklist_regions insertion_points_blacklist_hg19.txt --bam SCHN_1891.bam --blacklist_depth 200 --args.max_dist_btwn_ins_pts 35 -o temp1_out_SCHN_1891.txt

# This program reads in clean_intron_deletions (there can be multiple ones for a given sample and gene, one for each intron of the gene).
# For each sample gene, this program looks for the insertion point.
# An insertion point will manifest as a BND from the left side of the gene, and another BND from the right side of the gene,
# that both connect to the same place somewhere else in the genome (usually on another chromosome).
# Actually, those two BNDs don't connect to exactly the same place. Their two connection points will be approx. 11 bp appart.
# One of the connection may have an inserted sequence between the gene and the insertion point, usually poly-A (or poly-T).
#
# If a blacklist_regions file is provided, then this program will avoid those regions when looking for insertion_points.
# If a bam file is provided, then this program will avoid regions that have a high depth of blacklist_depth or above, default is 200 if blacklist_depth not provided.
# The program avoids such regions by ignoring them when it finds an insertion point in the region, continuing to look further for insertion_points.
#
# This program finds only the first insertion point (or none if there aren't any) for a given sample and gene.
# If a sample has more than one insertion point (for example, on different chromosomes), then this program finds only one of them.
# To find another insertion point (if there are any more for the same sample and gene), then you could add the found insertion point to the blacklist regions and run this program again.
# Thus, the program will not find the first insertion point (because it is now in a blacklisted region) and instead find the next one (if there are any)
# or find no insertion point (when this sample does not have any more insertion points for this gene).

# in_dels format,
# contains clean_intron_deletion (chrom, start_pos, end_pos), gene, gene_intron (chrom, start_pos, end_pos), gene_strand, vaf of the clean_intron_deletion
# mycohort     mysample        12      7281658 7281688 <DEL>   RBP5    12      7281661 7281689 -       0.63
# mycohort     mysample        15      41866014        41870084        <DEL>   TYRO3   15      41866014        41870082        +       0.3
# mycohort     mysample        18      48556993        48573289        <DEL>   SMAD4   18      48556994        48573288        +       0.39
# mycohort     mysample        18      48573665        48575055        <DEL>   SMAD4   18      48573666        48575054        +       0.29
# mycohort     mysample        18      48581364        48584495        <DEL>   SMAD4   18      48581364        48584493        +       0.36
# mycohort     mysample        18      48584826        48586235        <DEL>   SMAD4   18      48584827        48586234        +       0.56
# mycohort     mysample        18      48586285        48591791        <DEL>   SMAD4   18      48586287        48591791        +       0.54
# mycohort     mysample        18      48591975        48593387        <DEL>   SMAD4   18      48591977        48593387        +       0.47
# mycohort     mysample        18      48593556        48603006        <DEL>   SMAD4   18      48593558        48603006        +       0.39
# mycohort     mysample        18      48603147        48604626        <DEL>   SMAD4   18      48603147        48604624        +       0.43
# mycohort     mysample        2       179296982       179300872       <DEL>   PRKRA   2       179296982       179300870       -       0.24
# mycohort     mysample        2       179301046       179306336       <DEL>   PRKRA   2       179301047       179306335       -       0.29
# mycohort     mysample        2       179306431       179307993       <DEL>   PRKRA   2       179306432       179307992       -       0.32
# mycohort     mysample        2       179308111       179309148       <DEL>   PRKRA   2       179308112       179309147       -       0.42
# mycohort     mysample        2       179309227       179312231       <DEL>   PRKRA   2       179309228       179312230       -       0.45
# mycohort     mysample        2       179312313       179314968       <DEL>   PRKRA   2       179312314       179314967       -       0.42
# mycohort     mysample        2       179315139       179315693       <DEL>   PRKRA   2       179315139       179315691       -       0.4

# output format, insertion_point fields are filled only when insertion_point was found.
# Contains the left-most and right-most clean-intron-deletion coordinates, the point before the left-most and the point after the right-most coordinates that connect to the insertion-point (called retrocopy_start and retrocopy_end), and the insertion-points left and right positions (usually approx 11 bp apart).
# Insertion point can be on a different chromosome. The clean-intron-deletions and retrocopy start/end are on the same chromosome as each other.
# cohort       sample          gene	chrom	clean_intron_del_start	clean_intron_del_end	retrocopy_start	retrocopy_end	insertion_chrom	left_insertion_point	right_insertion_point	insertion_sequence	retrocopy_insertion_direction
# mycohort     mysample        RBP5	12	7281658	7281688						
# mycohort     mysample        TYRO3	15	41866014	41870084	41851396	41871536	13	44069827	44069838	AAAAAAAAAAAAAAAAAAAAAAAAA	FORWARD
# mycohort     mysample        SMAD4	18	48556993	48604626	48556624	48606219	9	127732636	127732713		REVERSE
# mycohort     mysample        PRKRA	2	179296982	179315693				

# gene_regions format:
# 1	11873	14409	DDX11L1	+
# 1	14361	29370	WASH7P	-
# 1	17368	17436	MIR6859-1	-
# 1	30365	30503	MIR1302-2	+

# blacklist_regions format:
# 2	33141300	33141700
# Y	17994840	17994950
# 13	44069825	44069826


# There are genes (eg. TYRO3) whose hg19 definition is too short because it doesn't include some longer transcripts defined in hg38.
# In that case, when the insertion point's discordant pair includes the longer transcript region that falls outside the official gene region, we fail to find the insertion point.
# In that case, when the insertion point's discordant pair does not extend all the way to the region that falls outside the official gene region
# and thus falls entirely within the official gene region, we do find the insertion point.
# Let's expand the region that we consider to belong to a gene in order to account for genes whose official transcripts are missing some extended region.
# How much to add? Let's see what the difference in gene length is between hg19 and hg38, for each gene.
# There are 129 genes whose hg19 and hg38 lengths are different by over 100,000 bp. eg.
# gene   hg19_length hg38_length   diff
# RAB27B       67041      177676 110635
# RAB38        62194      366016 303822
# RBFOX1     1694319     2473593 779274
# RBFOX3      426805      576257 149452

# There are some simple genomic regions such as GGGGGGGGGGGGGGGGGGGGGG where 100,000 reads might be mapped, 
# for which most of the read pair mates map elsewhere in the genome.
# Do not allow these regions to be identified as insertion points because they would be false positive insertion points.
# Format of the optional blacklist_regions file:
# chrom		start_pos	end_pos

# This program does not yet add the following fields to its output. They will be added by a script afterwards:
# cohort, sample, retrocopy_insertion_direction, is_left_BND_present, is_right_BND_present.
# If retrocopy_start is connected by a BND to left_insertion_point and retrocopy_end is connected by a BND to right_insertion_point, then retrocopy_insertion_direction="FORWARD".
# If retrocopy_start is connected by a BND to right_insertion_point and retrocopy_end is connected by a BND to left_insertion_point, then retrocopy_insertion_direction="REVERSE".
# The subsequent script will set retrocopy_insertion_direction="FORWARD" because the subsequent script doesn't know what it should be.
# The current program needs to be modified to correctly set retrocopy_insertion_direction.


import sys
import csv
import subprocess
import vcf # pip3 install --install-option="--prefix=/g/data/jb96/software/python_packages" pyvcf
import argparse
import pysam

# cd /tmp
# virtualenv test_env
# source test_env/bin/activate
#
# python3
# help('modules')

def parse_arguments(args):

    parser = argparse.ArgumentParser(description='Read in clean_intron_deletions that represent retrocopied genes. For each gene in a set of one of more of those deletions, look for the two insertion points in the input vcf of structural variants near the ends of the gene. Input reference gene_regions defines the gene ends. Output the insertion points for each gene.')
    parser.add_argument('-i', '--in_dels', dest='in_dels_file', 
                        help='Tab-delimited list of input clean_intron_deletions, sorted by gene, in format: CHROM POS END SVTYPE GENE')
    parser.add_argument('--in_sv', dest='in_sv_file', 
                        help='Input VCF file of structural variants containing BND records that represent potential retrocopied gene insertion points. Must be bgzip compressed and indexed with tabix so this program can retrieve by co-ordinates.')
    parser.add_argument('--gene_regions', dest='gene_regions_file', 
                        help='Tab-delimited list of gene regions, one line per gene, in format: CHROM START_POS END_POS GENE STRAND')
    parser.add_argument('--gene_region_extension_for_start_of_gene', dest='gene_region_extension_for_start_of_gene', required=False, 
                        help='Add this many basepairs before the beginning of the gene region to look for its insertion point')
    parser.add_argument('--gene_region_extension', dest='gene_region_extension', required=False, 
                        help='Add this many basepairs on after the end of the gene region to look for its insertion point')
    parser.add_argument('--blacklist_regions', dest='blacklist_regions', required=False, 
                        help='Do not find any insertion points in these blacklisted regions. If any found then keep looking for insertion points elsewhere.')
    parser.add_argument('--bam', dest='bam', required=False, 
                        help='Do not find any insertion points in regions that have a high read depth in this bam file. Defaults to 200x.')
    parser.add_argument('--blacklist_depth', dest='blacklist_depth', required=False, 
                        help='Do not find any insertion points in regions whose bam depth is this value or higher.')
    parser.add_argument('--max_dist_btwn_ins_pts', dest='max_dist_btwn_ins_pts', required=False, 
                        help='The maximum distances allowed between the 2 ends of an insertion point. Defaults to 35 bp.')
    parser.add_argument('-o', '--output', dest='output_file', 
                        help='Output file, one line per retrocopied gene, specifying the insertion points if found. The multiple input records of the same gene produce one output line.')

    args = parser.parse_args()

    if args.gene_region_extension is None:
      args.gene_region_extension = 0
    else:
      args.gene_region_extension = int(args.gene_region_extension)
    if args.blacklist_depth is None:
      args.blacklist_depth = 140
    else:
      args.blacklist_depth = int(args.blacklist_depth)

    #print('args.gene_region_extension')
    #print(args.gene_region_extension)

    return args


def is_integer(s):
  try:
    int(s)
    return True
  except ValueError:
    return False


def min_start(val1, val2):
  min_val = None
  if val1 is not None:
    if val2 is not None:
      min_val = min(val1, val2)
    else:
      min_val = val1
  else:
    min_val = val2
  return min_val


def max_end(val1, val2):
  max_val = None
  if val1 is not None:
    if val2 is not None:
      max_val = max(val1, val2)
    else:
      max_val = val1
  else:
    max_val = val2
  return max_val


def reverse_complement(inseq):
  outseq = ''
  complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
  for i in range(0, len(inseq)):
    this_char = inseq[i:(i+1)]
    this_char_compl = this_char
    if this_char in complement:
      this_char_compl = complement[this_char]
    outseq = this_char_compl + outseq
  return outseq


def longest_AAA_or_TTT(inseq):

  longest = 0
  A_length = 0
  T_length = 0
  prev_char = ''
  inseq = inseq.upper()
  for this_char in inseq:

    if (this_char == 'A'):
      if (A_length == 0):
        A_length = 1
      else:
        A_length = A_length + 1
      longest = max(longest, A_length)
    else:
      A_length = 0

    if (this_char == 'T'):
      if (T_length == 0):
        T_length = 1
      else:
        T_length = T_length + 1
      longest = max(longest, T_length)
    else:
      T_length = 0

    prev_char = this_char

  return longest


class CleanIntronDeletion:
    """
    This object represents a clean_intron_deletion.
    """
    def __init__(self, initial_chrom, initial_start, initial_end, initial_gene):
        self.chrom = initial_chrom
        self.start = initial_start
        self.end = initial_end
        self.gene = initial_gene

    def print(self, **kwargs):
        print( self.chrom + ":" + str(self.start) + "-" + str(self.end) + " " + self.gene)

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)


class InsertionPointPair:
    """
    This object represents the left and right insertion points of a retrocopied gene.
    """
    def __init__(self, initial_left_chrom, initial_left_pos, initial_left_alt_chrom, initial_left_alt_pos, initial_left_seq, initial_right_chrom, initial_right_pos, initial_right_alt_chrom, initial_right_alt_pos, initial_right_seq, initial_distance_between_insertion_points, initial_sequence, initial_sequence_longest_AAA_or_TTT, initial_retrocopy_insertion_direction):
        self.left_chrom = initial_left_chrom
        self.left_pos = initial_left_pos
        self.left_alt_chrom = initial_left_alt_chrom
        self.left_alt_pos = initial_left_alt_pos
        self.left_seq = initial_left_seq
        self.right_chrom = initial_right_chrom
        self.right_pos = initial_right_pos
        self.right_alt_chrom = initial_right_alt_chrom
        self.right_alt_pos = initial_right_alt_pos
        self.right_seq = initial_right_seq
        if initial_distance_between_insertion_points is not None:
          self.distance_between_insertion_points = initial_distance_between_insertion_points
        else:
          self.distance_between_insertion_points = -1
        if initial_sequence is not None:
          self.sequence = initial_sequence
        else:
          self.sequence = ''
        if initial_sequence_longest_AAA_or_TTT is not None:
          self.sequence_longest_AAA_or_TTT = initial_sequence_longest_AAA_or_TTT
        else:
          self.sequence_longest_AAA_or_TTT = 0
        self.retrocopy_insertion_direction = initial_retrocopy_insertion_direction

    def print(self, **kwargs):
        if self.left_chrom is not None:
          print( "left_chrom = " + str(self.left_chrom) )
        if self.left_pos is not None:
          print( "left_pos = " + str(self.left_pos) )
        if self.left_alt_chrom is not None:
          print( "left_alt_chrom = " + str(self.left_alt_chrom) )
        if self.left_alt_pos is not None:
          print( "left_alt_pos = " + str(self.left_alt_pos) )
        if self.left_seq is not None:
          print( "left_seq = " + str(self.left_seq) )
        if self.right_chrom is not None:
          print( "right_chrom = " + str(self.right_chrom) )
        if self.right_pos is not None:
          print( "right_pos = " + str(self.right_pos) )
        if self.right_alt_chrom is not None:
          print( "right_alt_chrom = " + str(self.right_alt_chrom) )
        if self.right_alt_pos is not None:
          print( "right_alt_pos = " + str(self.right_alt_pos) )
        if self.right_seq is not None:
          print( "right_seq = " + str(self.right_seq) )
        if self.distance_between_insertion_points is not None:
          print( "distance_between_insertion_points = " + str(self.distance_between_insertion_points) )
        if self.sequence is not None:
          print( "sequence =" + str(self.sequence) )
        if self.sequence_longest_AAA_or_TTT is not None:
          print( "sequence_longest_AAA_or_TTT = " + str(self.sequence_longest_AAA_or_TTT) )
        if self.retrocopy_insertion_direction is not None:
          print( "retrocopy_insertion_direction = " + str(self.retrocopy_insertion_direction) )

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)


class GeneAndItsCleanIntronDeletions:
    """
    This object represents one gene and its multiple (0, 1, or more) clean_intron_deletions.
    """
    def __init__(self, initial_gene):
        self.gene = initial_gene
        self.chrom = None
        self.min_start = None
        self.max_end = None
        self.gene_start = None
        self.gene_end = None
        self.list_of_clean_intron_deletions = []
        self.insertion_point_pair = None

    def print(self, **kwargs):
        print( "gene=" + self.gene )
        if self.chrom is not None:
          print( "chrom=" + str(self.chrom) )
        if self.min_start is not None:
          print( "min_start=" + str(self.min_start) )
        if self.max_end is not None:
          print( "max_end  =" + str(self.max_end) )
        if self.gene_start is not None:
          print( "gene_start=" + str(self.gene_start) )
        if self.gene_end is not None:
          print( "gene_end  =" + str(self.gene_end) )
        print("clean_intron_deletions:")
        for this_clean_intron_deletion in self.list_of_clean_intron_deletions:
          this_clean_intron_deletion.print()
        print("insertion_point_pairs:")
        for this_insertion_point_pair in self.list_of_insertion_point_pairs:
          this_insertion_point_pair.print()

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)


def read_clean_intron_deletions_for_genes( input_file_path ):

  # mycohort	sample222	12	7281658	7281688	<DEL>	RBP5	12	7281661	7281689	-	0.49
  # mycohort	sample222	15	41853506	41853734	<DEL>	TYRO3	15	41853509	41853735	+	0.38
  # mycohort	sample222	15	41866014	41870084	<DEL>	TYRO3	15	41866014	41870082	+	0.48
  # mycohort	sample222	22	30163536	30165665	<DEL>	UQCR10	22	30163538	30165665	+	0.47

  list_of_genes_and_their_clean_intron_deletions = []

  with open(input_file_path, newline='') as csvfile:
    clean_intron_deletion_data = list(csv.reader(csvfile, delimiter='\t'))
    #print(dir(clean_intron_deletion_data))

  # Record each gene in the loop
  this_gene = ''
  for this_clean_intron_deletion in clean_intron_deletion_data:
    next_gene = str(this_clean_intron_deletion[6])
    # Start recording a new gene
    if (next_gene != this_gene):
      if (this_gene != ''):
        list_of_genes_and_their_clean_intron_deletions.append( this_gene_and_its_dels_obj )
        #this_gene_and_its_dels_obj.print()
      this_gene = next_gene
      this_gene_and_its_dels_obj = GeneAndItsCleanIntronDeletions(this_gene)
    # Record this clean_intron_deletion for this gene
    this_cohort = str(this_clean_intron_deletion[0])
    this_sample = str(this_clean_intron_deletion[1])
    this_chrom = str(this_clean_intron_deletion[2])
    this_start = int(this_clean_intron_deletion[3])
    this_end = int(this_clean_intron_deletion[4])
    this_clean_intron_deletion = CleanIntronDeletion( this_chrom, this_start, this_end, this_gene )
    this_gene_and_its_dels_obj.list_of_clean_intron_deletions.append( this_clean_intron_deletion )
    this_gene_and_its_dels_obj.cohort = this_cohort
    this_gene_and_its_dels_obj.sample = this_sample
    this_gene_and_its_dels_obj.chrom = this_chrom
    this_gene_and_its_dels_obj.min_start = min_start( this_start, this_gene_and_its_dels_obj.min_start )
    this_gene_and_its_dels_obj.max_end = max_end( this_end, this_gene_and_its_dels_obj.max_end )

  # Make sure to record the last gene in the loop
  if (this_gene != ''):
    list_of_genes_and_their_clean_intron_deletions.append( this_gene_and_its_dels_obj )
    #this_gene_and_its_dels_obj.print()

  return list_of_genes_and_their_clean_intron_deletions


class GeneRegion:
    """
    This object represents information about a gene from a reference file
    """
    def __init__(self, initial_gene, initial_start, initial_end):
        self.gene = initial_gene
        self.start = initial_start
        self.end = initial_end

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)


class ChromRegion:
    """
    This object represents a chromosomal region of chrom:start-end
    """
    def __init__(self, initial_chrom, initial_start, initial_end):
        self.chrom = initial_chrom
        self.start = initial_start
        self.end = initial_end

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)


def read_all_genes_start_and_end( input_file_path, gene_region_extension, gene_region_extension_for_start_of_gene ):

  # 1	11873	14409	DDX11L1	+
  # 1	14361	29370	WASH7P	-
  # 1	17368	17436	MIR6859-1	-
  # 1	30365	30503	MIR1302-2	+

  gene_regions = {}

  with open(input_file_path, newline='') as csvfile:
    gene_data = list(csv.reader(csvfile, delimiter='\t'))

  for data in gene_data:
    chrom = str(data[0])
    start = int(data[1])
    end = int(data[2])
    gene = str(data[3])
    strand = str(data[4])
    if gene in gene_regions:
      start = min( start, gene_regions[gene].start )
      end = max( start, gene_regions[gene].end )
    gene_region = GeneRegion( gene, start, end )
    gene_regions[gene] = gene_region
  #print("FBRSL1 "+ str(gene_regions["FBRSL1"].start) + " " + str(gene_regions["FBRSL1"].end))

  for gene in gene_regions:
    gene_regions[gene].start = gene_regions[gene].start - gene_region_extension
    if (strand == "+"):
      gene_regions[gene].start = gene_regions[gene].start - gene_region_extension_for_start_of_gene
    if (gene_regions[gene].start < 1):
      gene_regions[gene].start = 1
    gene_regions[gene].end = gene_regions[gene].end + gene_region_extension
    if (strand == "-"):
      gene_regions[gene].end = gene_regions[gene].end + gene_region_extension_for_start_of_gene

  return gene_regions


def find_gene_start_and_end_for_gene( gene, gene_regions ):

  gene_info = gene_regions[ gene ]
  gene_start = gene_info.start
  gene_end = gene_info.end
  return (gene_start, gene_end)


def find_gene_start_and_end_for_list_of_genes_and_their_clean_intron_deletions( list_of_genes_and_their_clean_intron_deletions, gene_regions ):

  new_list_of_genes_and_their_clean_intron_deletions  = []

  for one_gene_and_its_clean_intron_deletions in list_of_genes_and_their_clean_intron_deletions:
    this_gene = one_gene_and_its_clean_intron_deletions.gene
    (gene_start, gene_end) = find_gene_start_and_end_for_gene( this_gene, gene_regions )
    one_gene_and_its_clean_intron_deletions.gene_start = gene_start
    one_gene_and_its_clean_intron_deletions.gene_end = gene_end
    new_list_of_genes_and_their_clean_intron_deletions.append( one_gene_and_its_clean_intron_deletions )

  return new_list_of_genes_and_their_clean_intron_deletions


def list_BNDs_by_alt_chrom( list_bnds ):

  list_bnds_by_chrom = {}

  for this_bnd in list_bnds:
    if (this_bnd.alt_chrom in list_bnds_by_chrom):
      chrom_list = list_bnds_by_chrom[ this_bnd.alt_chrom ]
    else:
      chrom_list = []
    chrom_list.append( this_bnd )
    list_bnds_by_chrom[ this_bnd.alt_chrom ] = chrom_list

  return list_bnds_by_chrom


class InsertionPointBndRecord:
    """
    This object represents a BND record that is a candidate for being one of the insertion points
    """
    def __init__(self, initial_chrom, initial_pos, initial_alt_chrom, initial_alt_pos, initial_sequence, initial_sequence_length, initial_sequence_longest_AAA_or_TTT):
        self.chrom = initial_chrom
        self.pos = initial_pos
        self.alt_chrom = initial_alt_chrom
        self.alt_pos = initial_alt_pos
        self.sequence = initial_sequence
        self.sequence_length = initial_sequence_length
        self.sequence_longest_AAA_or_TTT = initial_sequence_longest_AAA_or_TTT

    def print(self, **kwargs):
        print( str(self.pos) + " " + str(self.alt_chrom) + ":" + str(self.alt_pos) + " " + str(self.sequence_longest_AAA_or_TTT) + " " + str(self.sequence) + " " + str(self.sequence_length) )

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)


def does_region_overlap_regions( chrom, pos1, pos2, regions ):

  region_overlaps_regions = False

  chrom = str(chrom)
  pos1 = int(pos1)
  pos2 = int(pos2)
  start_pos = min(pos1, pos2)
  end_pos = max(pos1, pos2)

  keep_looking = True
  if (len(regions) == 0):
    keep_looking = False
  i = 0
  while (keep_looking):
    this_region = regions[i]
    if (chrom == this_region.chrom):
      if ((start_pos <= this_region.start) and (end_pos >= this_region.start)):
        region_overlaps_regions = True
      if ((start_pos <= this_region.end) and (end_pos >= this_region.end)):
        region_overlaps_regions = True
      if ((start_pos >= this_region.start) and (end_pos <= this_region.end)):
        region_overlaps_regions = True
    if (region_overlaps_regions):
      keep_looking = False
    i = i + 1
    if (i >= len(regions)):
      keep_looking = False

  return region_overlaps_regions


def is_bam_region_too_deep( args, chrom, pos1, pos2, retrogene_chrom, retrogene_pos1, retrogene_pos2, bam_pysam_handle ):

  bam_region_is_too_deep = False

  if args.bam is not None:

    chrom = str(chrom)
    pos1 = int(pos1)
    pos2 = int(pos2)
    start_pos = min(pos1, pos2)
    end_pos = max(pos1, pos2)
    if (start_pos == end_pos):
      end_pos = end_pos + 1

    max_depth_of_all_readings = 0

    get_depth_every_N_nucleotides = 10
    if ((end_pos - start_pos) > 100):
      get_depth_every_N_nucleotides = int((end_pos - start_pos) / 10) # no more than 10 interrogations
    num_interrogations = int((end_pos - start_pos) / get_depth_every_N_nucleotides)
    pos_upto = start_pos
    while (pos_upto < end_pos):

      pos_upto_plus_1 = pos_upto + 1

      # dir(pysam)
      # ['AlignedRead', 'AlignedSegment', 'AlignmentFile', 'AlignmentHeader', 'BGZFile', 'CBACK', 'CDEL', 'CDIFF', 'CEQUAL', 'CHARD_CLIP', 'CINS', 'CMATCH', 'CPAD', 'CREF_SKIP', 'CSOFT_CLIP', 'FDUP', 'FMREVERSE', 'FMUNMAP', 'FPAIRED', 'FPROPER_PAIR', 'FQCFAIL', 'FREAD1', 'FREAD2', 'FREVERSE', 'FSECONDARY', 'FSUPPLEMENTARY', 'FUNMAP', 'FastaFile', 'Fastafile', 'FastqFile', 'FastqProxy', 'FastxFile', 'FastxRecord', 'GZIterator', 'GZIteratorHead', 'HFile', 'HTSFile', 'IndexedReads', 'IteratorColumn', 'IteratorRow', 'KEY_NAMES', 'Pileup', 'PileupColumn', 'PileupRead', 'Samfile', 'SamtoolsError', 'TabixFile', 'Tabixfile', 'VCF', 'VCFRecord', 'VariantFile', 'VariantHeader', 'VariantHeaderRecord', 'VariantRecord', '__all__', '__builtins__', '__cached__', '__doc__', '__file__', '__loader__', '__name__', '__package__', '__path__', '__samtools_version__', '__spec__', '__version__', 'addreplacerg', 'array_to_qualitystring', 'asBed', 'asGFF3', 'asGTF', 'asTuple', 'asVCF', 'bam2fq', 'bamshuf', 'bedcov', 'calmd', 'cat', 'collate', 'config', 'depad', 'depth', 'dict', 'faidx', 'fasta', 'fastq', 'fixmate', 'flagstat', 'get_defines', 'get_include', 'get_libraries', 'get_verbosity', 'idxstats', 'index', 'libcalignedsegment', 'libcalignmentfile', 'libcbcf', 'libcbcftools', 'libcbgzf', 'libcfaidx', 'libchtslib', 'libcsamfile', 'libcsamtools', 'libctabix', 'libctabixproxies', 'libcutils', 'libcvcf', 'merge', 'mpileup', 'os', 'pad2unpad', 'phase', 'py_bcftools', 'py_samtools', 'pysam', 'qualities_to_qualitystring', 'qualitystring_to_array', 'quickcheck', 'reheader', 'rmdup', 'samimport', 'samtools', 'set_verbosity', 'sort', 'split', 'stats', 'sys', 'sysconfig', 'tabix_compress', 'tabix_file_iterator', 'tabix_generic_iterator', 'tabix_index', 'tabix_iterator', 'targetcut', 'tview', 'utils', 'version', 'view']

      # If this position has a lot of bad quality reads, we want to know, so that we can reject this region as having too many reads.
      #coverage_array = bam_pysam_handle.count_coverage( chrom, start=pos_upto, stop=pos_upto_plus_1, quality_threshold=0 )
      #print("coverage_array")
      #this_depth = coverage_array[0][0] + coverage_array[1][0] + coverage_array[2][0] + coverage_array[3][0]

      iter = bam_pysam_handle.fetch( chrom, int(pos_upto), int(pos_upto_plus_1) )
      this_depth = 0
      for this_read in iter:
        #print(this_read)
        #print(this_read.to_string())
        #print(this_read.next_reference_name)
        #print(this_read.next_reference_start)
        # If this read's mate pointing to the retrocopied_gene location, then don't count it in the depth,
        # because it is a read supporting this position being the insertion_point of the retrocopied_gene,
        # and thus we allow many reads.
        # We will not disqualify this insertion_point position for having lots of read depth that supports it being the insertion_point.
        # The depth check will only disqualify this insertion_point position when it has lots of reads whose mates are all over the place,
        # because that is obviously a region that has had all sorts of orphans reads erroneously mapped to it.
        # If there is a lot of read depth of reads whose mate is not the retrocopied_gene position,
        # then we assume that it is probably got read mates mapped all over the place, without actually checking for that.
        this_read_supports_this_insertion_point = False
        max_possible_dist_between_ins_pt_read_mate_position_and_retrogene_position = 3000
        if (this_read.next_reference_name == retrogene_chrom):
          if (abs(this_read.next_reference_start - retrogene_pos1) <= max_possible_dist_between_ins_pt_read_mate_position_and_retrogene_position):
            this_read_supports_this_insertion_point = True
          if (abs(this_read.next_reference_start - retrogene_pos2) <= max_possible_dist_between_ins_pt_read_mate_position_and_retrogene_position):
            this_read_supports_this_insertion_point = True
        if (this_read_supports_this_insertion_point == False):
          this_depth = this_depth + 1

      max_depth_of_all_readings = max(this_depth, max_depth_of_all_readings)

      pos_upto = pos_upto + get_depth_every_N_nucleotides

    if (max_depth_of_all_readings >= args.blacklist_depth):
      bam_region_is_too_deep = True

  return bam_region_is_too_deep


def find_candidate_insertion_point_BND_records( args, look_chrom, look_start, look_end, avoid_chrom, avoid_start, avoid_end, left_look_chrom, left_look_start, left_look_end, right_look_chrom, right_look_start, right_look_end, vcf_reader ):

  # This function looks at BND records between the gene_start and the first clean_intron_deletion,
  # and between the last clean_intron_deletion and the gene_end.
  # If this clean_intron_deletion is a retrocopied gene, then insertion points are expected just before the first and just after the last clean_intron_deletions.
  # So we look for BND records in the before-first region and after-last region, connecting to some other region outside of the gene.
  # An insertion points will not be from before-first to after-last positions, so do not pick up any of those BNDs.

  # For each insertion point, there really should be 2 BND records - one on this chromosome and one on the insertion chromosome.
  # However, it is possible that the quality of one of the BND records was not enough to be actually called.
  # We are looking for BND records on the clean_intron_deletion chromosome.
  # Really we should also look for BND records elsewhere that connect back to near the clean_intron_deletion.
  # We are not doing that though because it will be too time-consuming. 
  # Thus, if the clean_intron_deletion BND for the insertion is not present on the clean_intron_deletion chromosome,
  # then we will fail to find the insertion_point, even if there is the other BND present elsewhere that connects to this clean_intron_deletion chromosome at the insertion point.

  # Do not allow insertion points candidates in any blacklisted regions.

  list_of_candidate_insertion_point_BND_records = []

  # Make sure start < end, swap them around if necessary.
  if (look_start > look_end):
    temp_start = look_end
    look_end = look_start
    look_start = temp_start
  if (left_look_start > left_look_end):
    temp_start = left_look_end
    left_look_end = left_look_start
    left_look_start = temp_start
  if (right_look_start > right_look_end):
    temp_start = right_look_end
    right_look_end = right_look_start
    right_look_start = temp_start

  # This program might crash if the contig has no entries, with error:
  #
  #    for record in vcf_reader.fetch( look_chrom, look_start, look_end ):
  #  File "/g/data/jb96/software/python_packages/lib/python3.7/site-packages/vcf/parser.py", line 631, in fetch
  #    self.reader = self._tabix.fetch(chrom, start, end)
  #  File "pysam/libctabix.pyx", line 509, in pysam.libctabix.TabixFile.fetch
  #ValueError: could not create iterator for region 'hs37d5:7314085-7314154'

  interrogate_this_position = True
  try:
    vcf_reader.fetch( look_chrom, look_start, look_end )
  except:
    interrogate_this_position = False

  if (interrogate_this_position):
    for record in vcf_reader.fetch( look_chrom, look_start, look_end ):
      if ((record.var_type=='sv') and (record.var_subtype=='complex')):
        for this_alt in record.ALT:
          if (this_alt.type == 'BND'):
            bnd_chrom = record.CHROM
            bnd_pos = record.POS
            bnd_alt_chrom = this_alt.chr
            bnd_alt_pos = this_alt.pos
            #print('try ' + str(bnd_chrom) + ':' + str(bnd_pos) + ' ' + str(bnd_alt_chrom) + ':' + str(bnd_alt_pos) )
            if ((bnd_alt_chrom is None) or (bnd_alt_pos is None)):
              do_nothing = 1
            elif ((bnd_chrom == avoid_chrom) and (bnd_pos >= avoid_start) and (bnd_pos <= avoid_end)): # the BND is inside the gene
              do_nothing = 1
            elif ((bnd_alt_chrom == avoid_chrom) and (bnd_alt_pos >= avoid_start) and (bnd_alt_pos <= avoid_end)): # the BND partner is inside the gene
              do_nothing = 1
            else:

              # This BND record has one end inside a region of interest (either before-first or after-last), but make sure both its ends are not inside those regions.
              bnd_pos_is_inside_gene = False
              bnd_alt_is_inside_gene = False
              if ((bnd_chrom == left_look_chrom) and (bnd_pos >= left_look_start) and (bnd_pos <= left_look_end)):
                bnd_pos_is_inside_gene = True
              if ((bnd_alt_chrom == left_look_chrom) and (bnd_alt_pos >= left_look_start) and (bnd_alt_pos <= left_look_end)):
                bnd_alt_is_inside_gene = True
              if ((bnd_chrom == right_look_chrom) and (bnd_pos >= right_look_start) and (bnd_pos <= right_look_end)):
                bnd_pos_is_inside_gene = True
              if ((bnd_alt_chrom == right_look_chrom) and (bnd_alt_pos >= right_look_start) and (bnd_alt_pos <= right_look_end)):
                bnd_alt_is_inside_gene = True
              if (bnd_pos_is_inside_gene and bnd_alt_is_inside_gene):
              # BND and its partner position are both inside the region between the gene_start and first clean_intron_deletion, or last clean_intron_deletion and gene_end.
                do_nothing = 1
              else:

                #print( 'left ' + str(left_look_chrom) + ':' + str(left_look_start) + '-' + str(left_look_end) )
                #print( 'right ' + str(right_look_chrom) + ':' + str(right_look_start) + '-' + str(right_look_end) )
                #print( 'bnd ' + str(bnd_chrom) + ':' + str(bnd_pos) + ' ' + str(bnd_alt_chrom) + ':' + str(bnd_alt_pos) )
                connecting_sequence = this_alt.connectingSequence
                connecting_sequence_length = len(connecting_sequence)
                connecting_sequence_longest_AAA_or_TTT = longest_AAA_or_TTT(connecting_sequence)
                #print(record)
                #print(dir(record))
                #print(this_alt)
                #print(dir(this_alt))
                #print('chr=' + str(this_alt.chr))
                #print('connectingSequence=' + str(this_alt.connectingSequence))
                #print('connecting_sequence_longest_AAA_or_TTT=' + str(connecting_sequence_longest_AAA_or_TTT))
                #print('pos=' + str(this_alt.pos))
                #print('remoteOrientation=' + str(this_alt.remoteOrientation))
                #print('type=' + str(this_alt.type))
                #print('withinMainAssembly=' + str(this_alt.withinMainAssembly))
                #print('bnd_alt_chrom=' + str(bnd_alt_chrom) + ' bnd_alt_pos=' + str(bnd_alt_pos))
                this_candidate_insertion_point_BND_record = InsertionPointBndRecord(bnd_chrom, bnd_pos, bnd_alt_chrom, bnd_alt_pos, connecting_sequence, connecting_sequence_length, connecting_sequence_longest_AAA_or_TTT)
                list_of_candidate_insertion_point_BND_records.append( this_candidate_insertion_point_BND_record )
                #this_candidate_insertion_point_BND_record.print()

  return list_of_candidate_insertion_point_BND_records


def identify_insertion_points_from_candidates( args, left_insertion_point_candidates, right_insertion_point_candidates, blacklisted_regions, bam_pysam_handle, max_distance_between_insertion_points ):

  # Find the left BND closest to the first clean_intron_deletion that connects to the same chromosome that a right BND connects to, not far (within 35 bp) of the right BND connection position.
  # Find the right BND closest to the last clean_intron_deletion that connects to the same chromosome that a left BND connects to, not far (within 35 bp) of the left BND connection position.

  list_of_insertion_point_pairs = []

  left_insertion_point_candidates_by_chrom = list_BNDs_by_alt_chrom( left_insertion_point_candidates )
  right_insertion_point_candidates_by_chrom = list_BNDs_by_alt_chrom( right_insertion_point_candidates )

  max_ij = max( len(left_insertion_point_candidates), len(right_insertion_point_candidates) )
  found_insertion_points = False
  min_AAA_or_TTT_length_to_take_precedence_over_distance_between_insertion_points = 6
  best_insertion_point_pair = InsertionPointPair(None, None, None, None, None, None, None, None, None, None, None, None, None, None)

  i = len(left_insertion_point_candidates)
  for j in range( 0, max_ij ):
    i = i - 1

    # consider the next candidate BND on the left of the clean_intron_deletions
    if ( (i >= 0) and (i < len(left_insertion_point_candidates)) ):
      this_left_BND = left_insertion_point_candidates[i]

      if this_left_BND.alt_chrom in right_insertion_point_candidates_by_chrom:
        right_BNDs_having_alt_on_this_same_chrom_as_left_BND_alt_chrom = right_insertion_point_candidates_by_chrom[ this_left_BND.alt_chrom ]

        for this_right_BND in right_BNDs_having_alt_on_this_same_chrom_as_left_BND_alt_chrom:
          #this_right_BND.print()

          this_distance_between_insertion_points = abs( this_right_BND.alt_pos - this_left_BND.alt_pos )
          this_retrocopy_insertion_direction = "FORWARD" # if (this_left_BND.alt_pos < this_right_BND.alt_pos)
          if (this_right_BND.alt_pos < this_left_BND.alt_pos):
            this_retrocopy_insertion_direction = "REVERSE"

          # consider only BND pairs that are close to each other
          if (this_distance_between_insertion_points <= max_distance_between_insertion_points):

            is_this_putative_insertion_point_in_blacklisted_region = does_region_overlap_regions( this_left_BND.alt_chrom, this_left_BND.alt_pos, this_right_BND.alt_pos, blacklisted_regions )
            # is the putative insertion point in a blacklisted region?
            if (is_this_putative_insertion_point_in_blacklisted_region == False):

              is_this_putative_insertion_point_in_bam_region_having_too_much_depth = is_bam_region_too_deep( args, this_left_BND.alt_chrom, this_left_BND.alt_pos, this_right_BND.alt_pos, this_left_BND.chrom, this_left_BND.pos, this_right_BND.pos, bam_pysam_handle )
              # is the putative insertion point in a region having lots of depth and lots of reads erroneously mapped here and thus read mate BNDs are false positives?
              if (is_this_putative_insertion_point_in_bam_region_having_too_much_depth == False):

                if (this_left_BND.sequence_longest_AAA_or_TTT >= min_AAA_or_TTT_length_to_take_precedence_over_distance_between_insertion_points):

                  if (this_left_BND.sequence_longest_AAA_or_TTT > best_insertion_point_pair.sequence_longest_AAA_or_TTT):
                    best_insertion_point_pair = InsertionPointPair( this_left_BND.chrom, this_left_BND.pos, this_left_BND.alt_chrom, this_left_BND.alt_pos, this_left_BND.sequence, this_right_BND.chrom, this_right_BND.pos, this_right_BND.alt_chrom, this_right_BND.alt_pos, this_right_BND.sequence, this_distance_between_insertion_points, this_left_BND.sequence, this_left_BND.sequence_longest_AAA_or_TTT, this_retrocopy_insertion_direction)

                  elif (this_left_BND.sequence_longest_AAA_or_TTT == best_insertion_point_pair):
                    if (this_distance_between_insertion_points < best_insertion_point_pair.distance_between_insertion_points):
                      best_insertion_point_pair = InsertionPointPair( this_left_BND.chrom, this_left_BND.pos, this_left_BND.alt_chrom, this_left_BND.alt_pos, this_left_BND.sequence, this_right_BND.chrom, this_right_BND.pos, this_right_BND.alt_chrom, this_right_BND.alt_pos, this_right_BND.sequence, this_distance_between_insertion_points, this_left_BND.sequence, this_left_BND.sequence_longest_AAA_or_TTT, this_retrocopy_insertion_direction)

                  else: # (this_left_BND.sequence_longest_AAA_or_TTT < best_insertion_point_pair):
                    do_nothing = 1

                else: # (this_left_BND.sequence_longest_AAA_or_TTT < min_AAA_or_TTT_length_to_take_precedence_over_distance_between_insertion_points):

                  if (this_distance_between_insertion_points < best_insertion_point_pair.distance_between_insertion_points):

                    this_sequence = ''
                    this_sequence_longest_AAA_or_TTT = 0
                    if (this_right_BND.sequence_longest_AAA_or_TTT >= min_AAA_or_TTT_length_to_take_precedence_over_distance_between_insertion_points):
                      this_sequence = this_right_BND.sequence
                      this_sequence_longest_AAA_or_TTT = this_right_BND.sequence_longest_AAA_or_TTT

                    best_insertion_point_pair = InsertionPointPair( this_left_BND.chrom, this_left_BND.pos, this_left_BND.alt_chrom, this_left_BND.alt_pos, this_left_BND.sequence, this_right_BND.chrom, this_right_BND.pos, this_right_BND.alt_chrom, this_right_BND.alt_pos, this_right_BND.sequence, this_distance_between_insertion_points, this_sequence, this_sequence_longest_AAA_or_TTT, this_retrocopy_insertion_direction)

    # consider the next candidate BND on the right of the clean_intron_deletions
    if ( (j >= 0) and (j < len(right_insertion_point_candidates)) ):
      this_right_BND = right_insertion_point_candidates[j]

      if this_right_BND.alt_chrom in left_insertion_point_candidates_by_chrom:
        left_BNDs_having_alt_on_this_same_chrom_as_right_BND_alt_chrom = left_insertion_point_candidates_by_chrom[ this_right_BND.alt_chrom ]

        for this_left_BND in left_BNDs_having_alt_on_this_same_chrom_as_right_BND_alt_chrom:
          #this_left_BND.print()

          this_distance_between_insertion_points = abs( this_left_BND.alt_pos - this_right_BND.alt_pos )
          this_retrocopy_insertion_direction = "FORWARD" # if (this_left_BND.alt_pos < this_right_BND.alt_pos)
          if (this_right_BND.alt_pos < this_left_BND.alt_pos):
            this_retrocopy_insertion_direction = "REVERSE"

          # consider only BND pairs that are close to each other
          if (this_distance_between_insertion_points <= max_distance_between_insertion_points):

            is_this_putative_insertion_point_in_blacklisted_region = does_region_overlap_regions( this_left_BND.alt_chrom, this_left_BND.alt_pos, this_right_BND.alt_pos, blacklisted_regions )
            # is the putative insertion point in a blacklisted region?
            if (is_this_putative_insertion_point_in_blacklisted_region == False):

              is_this_putative_insertion_point_in_bam_region_having_too_much_depth = is_bam_region_too_deep( args, this_left_BND.alt_chrom, this_left_BND.alt_pos, this_right_BND.alt_pos, this_left_BND.chrom, this_left_BND.pos, this_right_BND.pos, bam_pysam_handle )
              # is the putative insertion point in a region having lots of depth and lots of reads erroneously mapped here and thus read mate BNDs are false positives?
              if (is_this_putative_insertion_point_in_bam_region_having_too_much_depth == False):

                if (this_right_BND.sequence_longest_AAA_or_TTT >= min_AAA_or_TTT_length_to_take_precedence_over_distance_between_insertion_points):

                  if (this_right_BND.sequence_longest_AAA_or_TTT > best_insertion_point_pair.sequence_longest_AAA_or_TTT):
                    best_insertion_point_pair = InsertionPointPair( this_right_BND.chrom, this_right_BND.pos, this_right_BND.alt_chrom, this_right_BND.alt_pos, this_right_BND.sequence, this_left_BND.chrom, this_left_BND.pos, this_left_BND.alt_chrom, this_left_BND.alt_pos, this_left_BND.sequence, this_distance_between_insertion_points, this_right_BND.sequence, this_right_BND.sequence_longest_AAA_or_TTT, this_retrocopy_insertion_direction)

                  elif (this_right_BND.sequence_longest_AAA_or_TTT == best_insertion_point_pair):
                    if (this_distance_between_insertion_points < best_insertion_point_pair.distance_between_insertion_points):
                      best_insertion_point_pair = InsertionPointPair( this_right_BND.chrom, this_right_BND.pos, this_right_BND.alt_chrom, this_right_BND.alt_pos, this_right_BND.sequence, this_left_BND.chrom, this_left_BND.pos, this_left_BND.alt_chrom, this_left_BND.alt_pos, this_left_BND.sequence, this_distance_between_insertion_points, this_right_BND.sequence, this_right_BND.sequence_longest_AAA_or_TTT, this_retrocopy_insertion_direction)

                  else: # (this_right_BND.sequence_longest_AAA_or_TTT < best_insertion_point_pair):
                    do_nothing = 1

                else: # (this_right_BND.sequence_longest_AAA_or_TTT < min_AAA_or_TTT_length_to_take_precedence_over_distance_between_insertion_points):

                  if (this_distance_between_insertion_points < best_insertion_point_pair.distance_between_insertion_points):

                    this_sequence = ''
                    this_sequence_longest_AAA_or_TTT = 0
                    if (this_left_BND.sequence_longest_AAA_or_TTT >= min_AAA_or_TTT_length_to_take_precedence_over_distance_between_insertion_points):
                      this_sequence = this_left_BND.sequence
                      this_sequence_longest_AAA_or_TTT = this_left_BND.sequence_longest_AAA_or_TTT

                    best_insertion_point_pair = InsertionPointPair( this_right_BND.chrom, this_right_BND.pos, this_right_BND.alt_chrom, this_right_BND.alt_pos, this_right_BND.sequence, this_left_BND.chrom, this_left_BND.pos, this_left_BND.alt_chrom, this_left_BND.alt_pos, this_left_BND.sequence, this_distance_between_insertion_points, this_sequence, this_sequence_longest_AAA_or_TTT, this_retrocopy_insertion_direction)

  list_of_insertion_point_pairs.append( best_insertion_point_pair )

  return list_of_insertion_point_pairs


def write_out_header(writer):

  outline = ["cohort", "sample", "gene", "chrom", "clean_intron_del_start", "clean_intron_del_end", "retrocopy_start", "retrocopy_end", "insertion_chrom", "left_insertion_point", "right_insertion_point", "insertion_sequence", "retrocopy_insertion_direction" ]
  writer.writerow(outline)
  return


def write_out_gene_and_its_insertion_points(writer, obj):

  obj2 = obj.list_of_insertion_point_pairs[0]
  retrocopy_start = obj2.left_pos
  retrocopy_end = obj2.right_pos
  if (obj2.left_pos is not None) and (obj2.right_pos is not None):
    if (obj2.left_pos > obj2.right_pos):
      retrocopy_start = obj2.right_pos
      retrocopy_end = obj2.left_pos
  left_insertion_point = obj2.left_alt_pos
  right_insertion_point = obj2.right_alt_pos
  if (obj2.left_alt_pos is not None) and (obj2.right_alt_pos is not None):
    if (obj2.left_alt_pos > obj2.right_alt_pos):
      left_insertion_point = obj2.right_alt_pos
      right_insertion_point = obj2.left_alt_pos

  outline = [obj.cohort, obj.sample, obj.gene, obj.chrom, obj.min_start, obj.max_end, retrocopy_start, retrocopy_end, obj2.left_alt_chrom, left_insertion_point, right_insertion_point, obj2.sequence, obj2.retrocopy_insertion_direction ]
  writer.writerow(outline)
  return


def find_insertion_points_for_deletions_that_are_retrocopied_genes(args, in_dels_path, in_sv_path, gene_regions_path, output_file_path, gene_region_extension, gene_region_extension_for_start_of_gene, blacklisted_regions, bam_pysam_handle, max_distance_between_insertion_points):

  writer = csv.writer(open(output_file_path, 'w'), delimiter='\t')
  write_out_header(writer)

  # Read in all the reference genes, find the start and end positions of each gene, ready to be used for any clean_intron_deletions in any gene.
  gene_regions = read_all_genes_start_and_end( gene_regions_path, gene_region_extension, gene_region_extension_for_start_of_gene )

  # Read in all the clean_intron_deletions, ready to process soon one by one.
  list_of_genes_and_their_clean_intron_deletions = read_clean_intron_deletions_for_genes(in_dels_path)
  list_of_genes_and_their_clean_intron_deletions = find_gene_start_and_end_for_list_of_genes_and_their_clean_intron_deletions( list_of_genes_and_their_clean_intron_deletions, gene_regions )

  # Open the tabix indexed VCF file of BND records, ready to look for the insertion point BND records for each clean_intron_deletion.
  sv_vcf_reader = vcf.Reader(filename=in_sv_path)

  # Process each gene having clean_intron_deletions to find the 2 insertion points.
  for this_gene_and_its_dels in list_of_genes_and_their_clean_intron_deletions:

    left_look_chrom = this_gene_and_its_dels.chrom
    left_look_start = this_gene_and_its_dels.gene_start
    left_look_end = this_gene_and_its_dels.min_start
    right_look_chrom = this_gene_and_its_dels.chrom
    right_look_start = this_gene_and_its_dels.max_end
    right_look_end = this_gene_and_its_dels.gene_end

    avoid_chrom = this_gene_and_its_dels.chrom
    avoid_start = this_gene_and_its_dels.min_start
    avoid_end = this_gene_and_its_dels.max_end
    left_insertion_point_candidates = find_candidate_insertion_point_BND_records( args, left_look_chrom, left_look_start, left_look_end, avoid_chrom, avoid_start, avoid_end, left_look_chrom, left_look_start, left_look_end, right_look_chrom, right_look_start, right_look_end, sv_vcf_reader )
    #print('left_insertion_point_candidates')
    #for this_one in left_insertion_point_candidates:
    #    print(this_one)

    right_insertion_point_candidates = find_candidate_insertion_point_BND_records( args, right_look_chrom, right_look_start, right_look_end, avoid_chrom, avoid_start, avoid_end, left_look_chrom, left_look_start, left_look_end, right_look_chrom, right_look_start, right_look_end, sv_vcf_reader )
    #print('right_insertion_point_candidates')
    #for this_one in right_insertion_point_candidates:
    #    print(this_one)

    list_of_insertion_point_pairs = identify_insertion_points_from_candidates( args, left_insertion_point_candidates, right_insertion_point_candidates, blacklisted_regions, bam_pysam_handle, max_distance_between_insertion_points )
    this_gene_and_its_dels.list_of_insertion_point_pairs = list_of_insertion_point_pairs

    write_out_gene_and_its_insertion_points(writer, this_gene_and_its_dels)

    print(" ")
    print("a gene and its clean_intron_deletions:")
    this_gene_and_its_dels.print()
    print(" ")

  return


def get_blacklist_regions(args):

  blacklisted_regions = []
  if args.blacklist_regions is not None:

    with open(args.blacklist_regions, newline='') as csvfile:
      region_data = list(csv.reader(csvfile, delimiter='\t'))

    for data in region_data:
      chrom = str(data[0])
      start = int(data[1])
      end = int(data[2])
      this_gene_region = ChromRegion( chrom, start, end )
      blacklisted_regions.append( this_gene_region )

  return blacklisted_regions


def open_bam_file(args):

  if args.bam is None:
    bam_pysam_handle = ''
  else:
    bam_pysam_handle = pysam.AlignmentFile( args.bam, "rb" )

  return bam_pysam_handle


def main(args):

  args = parse_arguments(args)
  in_dels_path = args.in_dels_file
  in_sv_path = args.in_sv_file
  gene_regions_path = args.gene_regions_file
  output_file_path = args.output_file
  max_distance_between_insertion_points = 35
  if args.max_dist_btwn_ins_pts is not None:
    max_distance_between_insertion_points = int(args.max_dist_btwn_ins_pts)
  gene_region_extension = 100000
  if args.gene_region_extension is not None:
    gene_region_extension = int(args.gene_region_extension)
  gene_region_extension_for_start_of_gene = 20
  if args.gene_region_extension_for_start_of_gene is not None:
    gene_region_extension_for_start_of_gene = int(args.gene_region_extension_for_start_of_gene)
  print(' ')

  blacklisted_regions = get_blacklist_regions(args)
  bam_pysam_handle = open_bam_file(args)

  find_insertion_points_for_deletions_that_are_retrocopied_genes(args, in_dels_path, in_sv_path, gene_regions_path, output_file_path, gene_region_extension, gene_region_extension_for_start_of_gene, blacklisted_regions, bam_pysam_handle, max_distance_between_insertion_points)

  print(' ')
  print('find_insertion_points_for_deletions_that_are_retrocopied_genes.py')


if __name__ == "__main__":
    main(sys.argv)


