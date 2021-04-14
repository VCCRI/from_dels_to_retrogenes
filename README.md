# from_dels_to_retrogenes

# Identify retrocopied genes from whole genome sequencing structural variant clean intron deletions

This small suite of scripts identifies retrocopied genes in whole genome sequencing data (WGS) from structural variant (SV) deletion variants (DEL or INDEL) that correspond to an entire intron.
Prior to running this small pipeline, Gridss ([Cameron et al. 2017](https://genome.cshlp.org/content/27/12/2050), [PMID: 29097403](https://pubmed.ncbi.nlm.nih.gov/29097403/), software available at [https://github.com/PapenfussLab/gridss](https://github.com/PapenfussLab/gridss)) or equivalent software needs to be run, to generate the SV BND records that will be used by this pipeline to identify retrocopied genes.

## Usage

```bash
# Identify retrocopied genes and their insertion points in each sample,
# by running the following for each sample.
# If your data is hg19 and not hg38, then genome parameters should be hg19.
#
./run_1_call_dels_from_bnds.sh hg38 cohort_id sample_id /my/output/directory /my/temp/directory /my/gridss/output.vcf.gz
./run_2_look_for_clean_intron_deletions.sh hg38 cohort_id sample_id /my/output/directory /my/temp/directory /my/gridss/output.vcf.gz
./run_3_find_insertion_points_for_retrocopied_genes.sh hg38 cohort_id sample_id /my/output/directory /my/temp/directory /my/gridss/output.vcf.gz /my/PYTHONPATH/python_packages/lib/python3.N/site-packages /a/representative/bam/file
#
# After running all samples individually, concatenate the results from all samples,
#
# by running the following once for all samples together.
#
./run_4_concat_insertion_points_for_all_samples.sh /my/output/directory
#
# Look for subsequent insertion points for the retrocopied genes that already have insertion points,
# by running the following one script for each sample and then running the concat script for all samples.
# Repeat these two scripts, incrementing the last parameter each round so that the output file names can be different to previous rounds,
# until no new insertion points are identified, which will manifest as an empty output file.
#
./run_5_find_subsequent_insertion_points.sh /my/output/directory hg38 cohort_id sample_id /my/output/directory /my/temp/directory /my/gridss/output.vcf.gz /my/PYTHONPATH /a/representative/bam/file 2
./run_4_concat_insertion_points_for_all_samples.sh /my/output/directory
./run_5_find_subsequent_insertion_points.sh /my/output/directory hg38 cohort_id sample_id /my/output/directory /my/temp/directory /my/gridss/output.vcf.gz /my/PYTHONPATH /a/representative/bam/file 3
./run_4_concat_insertion_points_for_all_samples.sh /my/output/directory
./run_5_find_subsequent_insertion_points.sh /my/output/directory hg38 cohort_id sample_id /my/output/directory /my/temp/directory /my/gridss/output.vcf.gz /my/PYTHONPATH /a/representative/bam/file 4
./run_4_concat_insertion_points_for_all_samples.sh /my/output/directory
#
# Use the previously identified insertion points to identify insertion points in retrocopied genes that do not yet have insertion points identified,
# by running the following for each sample.
#
./run_6_find_insertion_points_using_other_samples.sh hg38 cohort_id sample_id /my/output/directory /my/temp/directory /my/gridss/output.vcf.gz /my/PYTHONPATH /a/representative/bam/file
#
# For each sample having insertion points including using other samples,
# annotate with the reference sequence of the insertion point.
#
./run_7_get_insertion_point_sequence_from_reference.sh hg38 cohort_id sample_id /my/output/directory /my/temp/directory reference_fasta /my/PYTHONPATH
#
# The output produced by the above, *.insertion_points_for_retrocopied_genes_including_using_other_samples.insertion_pt_ref_seq.tsv, does include the initially found insertion points and those found using the insertion points of other samples, but does not include the subsequently found insertion points for genes for which insertion points were initially found.
# The subsequently found insertion points for genes for which insertion points were intially found are output in *.insertion_points_for_deletions_that_are_retrocopied_genes_N.tsv where N is the round that it was found in.
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

Rath EM et al. Manuscript in preparation.


