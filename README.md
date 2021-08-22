# from_dels_to_retrogenes

# Identify retrocopied genes from whole genome sequencing structural variant clean intron deletions

This small pipeline of scripts identifies retrocopied genes in whole genome sequencing data (WGS) from structural variant (SV) deletion variants (DEL or INDEL) that correspond to an entire intron, referred to here as clean intron deletions.

Prior to running this small pipeline, Gridss ([Cameron et al. 2017](https://genome.cshlp.org/content/27/12/2050), [PMID: 29097403](https://pubmed.ncbi.nlm.nih.gov/29097403/), software available at [https://github.com/PapenfussLab/gridss](https://github.com/PapenfussLab/gridss)) or equivalent software that generates VCF BND records for deletions and translocations needs to be run, to generate the SV BND records that will be used by this pipeline to identify retrocopied genes and their insertion points.

## Usage

```bash
# Make sure the reference files provided with this pipeline are decompressed.
cd from_dels_to_retrogenes/reference_data
gunzip *.gz
cd ../code

# Identify retrocopied genes and their insertion points in every sample in the input directory by running the following scripts.
# If your data is hg19 and not hg38, then genome parameters should be hg19 instead of hg38.
# It is mandatory for each sample to have a gridss output vcf file or similar file with BND records representing deletions and translocation.
# These scripts assume that the gridss output files are bgzip compressed (.vcf.gz) and tabix indexed (.tbi) with file name <SAMPLE_ID>.*.vcg.gz.
# The BND records are used to identify the clean-intron-deletions of retrocopied genes and to identify insertion-points.
# It is optional to provide a bam file for each sample, which can be a representative bam file of the same sequencing depth and reference genome but doesn't have to be for the same sample.
# If provided, the bam file will be consulted to determine whether insertion points fall in a region having too much sequencing depth.
# These scripts assume that bam files as named as <SAMPLE_ID>.bam, and indexed .bai file.
# These scripts assume that the same output directory is used for all scripts and thus, with the exception of the first script, that the input directory is the same as the output directory.
#
# Identify retrocopied genes and their insertion points in every sample in the input directory.
./run_1_call_dels_from_bnds_in_cohort.sh hg38 MY_COHORT_ID . /my/output/directory /my/temp/directory /my/directory/of/gridss/files
./run_2_look_for_clean_intron_deletions_in_cohort.sh hg38 MY_COHORT_ID . /my/output/directory /my/temp/directory
./run_3_find_insertion_points_for_retrocopied_genes_in_cohort.sh hg38 MY_COHORT_ID . /my/output/directory /my/temp/directory /my/directory/of/gridss/files /my/PYTHON_PATH /my/directory/of/bam/files
#
# Concatenate the results from all samples.
./run_4_concat_insertion_points_for_all_samples.sh /my/output/directory
#
# Look for subsequent insertion points for the retrocopied genes that already have insertion points,
# because it is possible for a parental gene to be retrocopied more than once in the same sample.
./run_5_find_subsequent_insertion_points_in_cohort.sh hg38 MY_COHORT_ID . /my/output/directory /my/temp/directory /my/directory/of/gridss/files /my/PYTHON_PATH /my/directory/of/bam/files
./run_5b_concat_insertion_points_for_all_samples.sh /my/output/directory
#
# Use the previously identified insertion points to identify insertion points in retrocopied genes that do not yet have insertion points identified.
./run_6_find_insertion_points_using_other_samples_in_cohort.sh hg38 MY_COHORT_ID . /my/output/directory /my/temp/directory /my/directory/of/gridss/files /my/PYTHON_PATH /my/directory/of/bam/files /output/file/from/run_5b_concat_insertion_points_for_all_samples.sh/script
#
# Annotate insertion points with the reference sequence of the insertion point.
./run_7_get_insertion_point_sequence_from_reference_in_cohort.sh hg38 MY_COHORT_ID . /my/output/directory /my/temp/directory /my/reference/fasta/file /my/PYTHON_PATH
#
# Concatenate the retrogene results from the multiple samples merged output file, regardless of whether insertion_point was found or not.
# The files are:
# /my/output/directory/All__insertion_points_for_deletions_that_are_retrocopied_genes_inspt_refseq.sort_by_gene.txt
# /my/output/directory/All__insertion_points_for_deletions_that_are_retrocopied_genes_inspt_refseq.sort_by_sample.txt
./run_7b_concat_insertion_points_for_all_samples.sh /my/output/directory
#
```

## Dependencies

```bash
bash
awk
python3
bedtools
bgzip (part of htslib)
tabix (part of htslib)

python3 dependencies, to be present in PYTHONPATH:
import argparse
import commands
import csv
import datetime
import math
import os
import pysam
import random
import re
import subprocess
import sys
import vcf # pip3 install --install-option="--prefix=/my/python_packages" pyvcf
```

## Citation

Rath EM et al. 2021. Manuscript in preparation.


