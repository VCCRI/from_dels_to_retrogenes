# python3 get_insertion_point_sequence_from_reference.py --input infile.txt --output outfile.txt --ref_fasta /home/emma/emma/reference_genomes/hs38dhx/hs38dhx.fa

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2021, Victor Chang Medical Research Institute'

import sys
import csv
import pysam
import subprocess
import argparse

######################################################

def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def parse_arguments(args):
    parser = argparse.ArgumentParser(description='The input file contains the two insertion points of a retrocopied gene. This program gets the reference sequence of that interval and put it in the output file.')
    parser.add_argument('-i', '--input', dest='input_file', 
                        help='Input file')
    parser.add_argument('-o', '--output', dest='output_file', 
                        help='Output file')
    parser.add_argument('--ref_fasta', dest='ref_fasta', 
                        help='The path to the reference genome fasta from which the reference sequence will be obtained.')
    return parser.parse_args()


def read_infile(input_file_path):
  with open(input_file_path, newline='') as csvfile:
    data = list(csv.reader(csvfile, delimiter='\t'))
  return data


def write_out_line(writer, line, new_col):

  outline = []
  for this_col in line:
    outline.append(str(this_col))
  outline.append(str(new_col))
  writer.writerow(outline)

  return


def get_ref_seq( chrom, start_pos, end_pos, ref_fasta ):

  fasta_seq = ''

  idx = chrom.find("*") # * for unmapped, HLA-DRB1*12:01:01. Don't deal with this, too complicated, just avoid it for now.
  if (idx > -1):
    do_nothing = 1
  else:

    if (int(start_pos) > int(end_pos)):
      temp = start_pos
      start_pos = end_pos
      end_pos = temp

    locus = chrom + ':' + str(start_pos) + '-' + str(end_pos)
    fasta_result = pysam.faidx(ref_fasta, locus)
    fasta_seq = ''
    fasta_result_lines = fasta_result.split("\n")
    if (len(fasta_result_lines) > 1):
      for i in range( 1, len(fasta_result_lines) ):
        fasta_seq = fasta_seq + fasta_result_lines[i]

  return fasta_seq


def get_insertion_point_sequence_from_reference(input_file_path, output_file_path, ref_fasta):

  infile = read_infile(input_file_path)

  writer = csv.writer(open(output_file_path, 'w'), delimiter='\t')

  in_header = True
  for inline in infile:
    # print(inline)

    if (in_header):

      write_out_line(writer, inline, "insertion_point_reference_sequence")
      in_header = False

    else:

      insertion_pt_ref_seq = ""
      ins_pt_chrom = str(inline[16])
      ins_pt_start_pos = str(inline[17])
      ins_pt_end_pos = str(inline[18])

      if ((is_integer(ins_pt_start_pos)) and (is_integer(ins_pt_end_pos))):
        insertion_pt_ref_seq = get_ref_seq(ins_pt_chrom, ins_pt_start_pos, ins_pt_end_pos, ref_fasta)

      write_out_line(writer, inline, insertion_pt_ref_seq)

  return


def main(args):

  args = parse_arguments(args)
  input_file_path = args.input_file
  output_file_path = args.output_file
  ref_fasta = args.ref_fasta

  get_insertion_point_sequence_from_reference(input_file_path, output_file_path, ref_fasta)

  print(' ')
  print('End of get_insertion_point_sequence_from_reference')


if __name__ == "__main__":
    main(sys.argv)


