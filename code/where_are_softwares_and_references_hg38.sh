#!/bin/bash

export genome_version='hg38'
export ref_fasta=/g/data/jb96/References_and_Databases/hg38.noalt.decoy.bwa/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
export ref_fasta_fa=/g/data/jb96/References_and_Databases/hg38.noalt.decoy.bwa/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa
export ref_fasta_with_alts=/g/data/jb96/References_and_Databases/Homo_sapiens_assembly38.fasta/Homo_sapiens_assembly38.fasta
export ref_fasta_chroms=/g/data/jb96/References_and_Databases/hg38.noalt.decoy.bwa.chroms
export sw=/g/data/jb96/software
export humandb=$sw/annovar/humandb
export gatk_path=$sw/GATK/gatk-4.2.0.0/gatk
export gatk3_path=$sw/GATK/GenomeAnalysisTK-3.8-1-0/GenomeAnalysisTK.jar
export bwa=$sw/bwa-0.7.17/bwa
export picard_jar=$sw/picard/picard-2.18.26/picard.jar
# For GATK HaplotypeCaller:
export gatk_dbsnp=/g/data3/a32/References_and_Databases/GATK_bundle/hg38/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz
# For GATK VariantRecalibrator:
export gatk_db_dir=/g/data/jb96/References_and_Databases/GATK_bundle/hg38
export mills_resource_vcf=$gatk_db_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
export axiomPoly_resource_vcf=$gatk_db_dir/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
export dbsnp_resource_vcf=$gatk_db_dir/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz
export hapmap_resource_vcf="$gatk_db_dir/hapmap_3.3.hg38.vcf.gz"
export omni_resource_vcf="$gatk_db_dir/1000G_omni2.5.hg38.vcf.gz"
export one_thousand_genomes_resource_vcf="$gatk_db_dir/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
# For annovar annotation
export wgEncodeGencodeBasic=wgEncodeGencodeBasicV33
export EHRF_regions=EHRF_r104_regions
export exac03_ed=exac03_ed1_high
# For vep annotation
export vep_assembly=GRCh38
export vep_loftee=/g/data/jb96/software/vep_plugins/LOFTEE_GRCh38/loftee
export vep_loftee_human_ancestor_fa=/g/data/jb96/software/vep_plugins/LOFTEE_GRCh38/human_ancestor/human_ancestor.fa.gz
export vep_loftee_conservation_file=/g/data/jb96/software/vep_plugins/LOFTEE_GRCh38/conservation_file/phylocsf_gerp.sql
export vep_utrannotator=/g/data/jb96/software/vep_plugins/UTRannotator/UTRannotator/uORF_starts_ends_GRCh38_PUBLIC.txt
# For gridss structural variant calling
export gridss_blacklist="${sw}"/gridss_resources/hg38-blacklist.v2.bed
export regions_of_LowComplexity_SimpleRepeats=/g/data/jb96/software/public_domain_MGRB_structural_variant_annotation/structural_variant_annotation/reference_data/UCSC_GRCh38_Repeats_20170921_LowComplexity_SimpleRepeats.bed
# For structural variant annotation
export ucsc_refseq_cdsexons=/g/data/jb96/software/public_domain_MGRB_structural_variant_annotation/structural_variant_annotation/reference_data/UCSC_GRCh38_GenesAndGenePredictions_CDSexons_RefSeq_20211116.txt
export ucsc_refseq_genes=/g/data/jb96/software/public_domain_MGRB_structural_variant_annotation/structural_variant_annotation/reference_data/UCSC_GRCh38_GenesAndGenePredictions_genes_RefSeq_20211116_plus_TAPVR1.txt
export ucsc_refseq_introns=/g/data/jb96/software/public_domain_MGRB_structural_variant_annotation/structural_variant_annotation/reference_data/UCSC_GRCh38_GenesAndGenePredictions_introns_RefSeq_20211116.txt
export ucsc_gencode_cdsexons=/g/data/jb96/software/public_domain_MGRB_structural_variant_annotation/structural_variant_annotation/reference_data/UCSC_GRCh38_GenesAndGenePredictions_CDSexons_wgEncodeGencodeCompV38_20211116.txt
export ucsc_gencode_genes=/g/data/jb96/software/public_domain_MGRB_structural_variant_annotation/structural_variant_annotation/reference_data/UCSC_GRCh38_GenesAndGenePredictions_genes_wgEncodeGencodeCompV38_20211116.txt
export ucsc_gencode_introns=/g/data/jb96/software/public_domain_MGRB_structural_variant_annotation/structural_variant_annotation/reference_data/UCSC_GRCh38_GenesAndGenePredictions_introns_wgEncodeGencodeCompV38_20211116.txt
#export refdata_EHRFr90=/g/data/jb96/software/annovar/humandb/hg38_EHRF_chr.txt
export refdata_EHRFr99=/g/data/jb96/software/annovar/humandb/hg38_EHRF_r104_chr.txt
export refdata_gnomadsv=/g/data/jb96/software/annovar/humandb/hg38_gnomAD-SV_short_ordered_chr.txt
export refdata_DGVgoldStdVar_r20160516=/g/data/jb96/software/annovar/humandb/hg38_DGVgoldStdVars_r20160516.txt
export refdata_segmental_duplications=/g/data/jb96/software/annovar/humandb/GRCh38GenomicSuperDup.tab
export refdata_dbscSNV=/g/data/jb96/software/annovar/humandb/hg38_splice_dbscSNV_1_0p6.txt
export refdata_FANTOM5_TSS=/g/data/jb96/software/annovar/humandb/hg38_FANTOM5_CAGEr_humanHeart_transcriptStartSites_1column.txt
export refdata_FANTOM5_TSS_5cols=/g/data/jb96/software/annovar/humandb/hg38_FANTOM5_CAGEr_humanHeart_transcriptStartSites_5columns.txt
export pext_artery_aorta=/g/data/jb96/software/annovar/humandb/hg38_pext_sv_for_exons_of_Artery_Aorta.bed
export pext_artery_coronary=/g/data/jb96/software/annovar/humandb/hg38_pext_sv_for_exons_of_Artery_Coronary.bed
export pext_heart_atrial_appendage=/g/data/jb96/software/annovar/humandb/hg38_pext_sv_for_exons_of_Heart_AtrialAppendage.bed
export pext_heart_left_ventricle=/g/data/jb96/software/annovar/humandb/hg38_pext_sv_for_exons_of_Heart_LeftVentricle.bed
export pext_genes=/g/data/jb96/software/annovar/humandb/hg38_pext_sv_for_genes.bed
export hotspot_exons=/g/data/jb96/software/annovar/humandb/hg38_hotspot_exons_sv.bed
# For splicing/splice sites
export spliceogen_sw="${sw}"/spliceogen_2020march_updated_2021march # This crashes on non-standard nucleotide in reference fasta: Running MaxEntScan... Can't take log of 0 at score5.pl line 96, <FILE> line 144119. java.lang.ArrayIndexOutOfBoundsException: Index 1 out of bounds for length 1. at processScoresMES.main(processScoresMES.java:27). MaxEntScan returned non-zero exit status. It is likely not all variants were processed. Exiting...
#export spliceogen_sw="${sw}"/spliceogen_2020july # This crashes with unsolved error: gzip: stdout: Broken pipe. java.lang.ArrayIndexOutOfBoundsException: Index 1 out of bounds for length 1 at mergeOutput.main(mergeOutput.java:49). sort: write failed: 'standard output': Broken pipe. sort: write error
export fasta_for_spliceogen=/g/data/a32/References_and_Databases/hg38.fa_faidx/hg38.fa
export gtf_for_spliceogen=/g/data/jb96/software/spliceogen_2020july/resources/gencode.v33.basic.annotation.gtf.gz
export ese_for_spliceogen=/g/data/jb96/software/spliceogen_2020july/resources/ESE.txt
export ess_for_spliceogen=/g/data/jb96/software/spliceogen_2020july/resources/ESS.txt
export refseq_gene_protein_lengths=/g/data/jb96/software/public_domain_MGRB_structural_variant_annotation/structural_variant_annotation/reference_data/UCSC_GRCh38_GenesAndGenePredictions_geneProteinLengths_RefSeq_20200321.txt
# For spliceai tensorflow program
export spliceai_reference=/g/data/jb96/References_and_Databases/grch38.p12/GRCh38.p12.genome.fa
# For ConanVarvar
#export R_LIBS_USER=/g/data/jb96/software/conanvarvar/R_libraries_for_conanvarvar # conanvarvar_1_run_conanvarvar_on_directory_of_bams.pbs needs this
# For many other scripts, so can't have these 2 different R_LIBS_USER in this one where_are_softwares_and_references.sh script
export R_LIBS_USER=/g/data/jb96/software/victorchang_scripts/R_libraries_for_victorchang_scripts # strucVars_extract_and_view_annotations.pbs and many scripts need this
# For liftover
export liftover_chain=/g/data/jb96/software/liftover/hg38ToHg19.over.chain.gz
export liftover_to_genome=hg19
# For atlantool
export atlantool_path=/g/data/jb96/software/atlantool/atlantool/releases/atlantool-linux
# For annotation of genes by bedtools intersect of gene regions, has multiple annotation columns
export gnomad_constraints=/g/data/jb96/software/annovar/humandb/hg38_gnomadv211_lof_metrics_by_gene.txt
