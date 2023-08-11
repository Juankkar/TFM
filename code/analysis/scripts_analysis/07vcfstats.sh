#!/usr/bin/env bash

## Path with the files to use
path_logs=../../../metadata/logs/vcfstats/
path_basic_stats=../results_analysis/tables/vcf_basic_stats.tsv
path_add_stats=../results_analysis/tables/additional_vcf_stats.tsv

function vcf_basicstats_filter(){

  operation=$(sed -E 's/\s:\s/\t/g' ${path_logs}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E $1 \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')

      echo $operation
}

function vcf_stats_filter() {

  operation=$(sed -E 's/\s:\s/\t/g' ${path_logs}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E $1 \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

    echo $operation
}


## Headers of the first TSV file (basic VCF statistics)
echo -e "pass_fail_filter\t\
snps\t\
mnps\t\
insertions\t\
deletions\t\
indels\t\
same_reference\t\
sample" > $path_basic_stats

for sample in ERR696683 ERR7533{68..78}
do

  passed_filters=$(vcf_basicstats_filter "^Pass")
  failed_filters=$(vcf_basicstats_filter "^Fail")
  snps=$(vcf_basicstats_filter "^SNPs")
  mnps=$(vcf_basicstats_filter "^MNPs")
  inser=$(vcf_basicstats_filter "^Insertions")
  del=$(vcf_basicstats_filter "^Deletions")
  indel=$(vcf_basicstats_filter "^Indels")
  same_ref=$(vcf_basicstats_filter "^Same")

  # Values of the TSV file (basic VCF statistics)
  echo -e "$passed_filters|$failed_filters\t\
  $snps\t\
  $mnps\t\
  $inser\t\
  $del\t\
  $indel\t\
  $same_ref\t\
  $sample"

done >> $path_basic_stats

## Headers of the second TSV file (additional VCF statistics)
echo -e "snp_transition_transversion\t\
total_het_hom\t\
snp_het_hom\t\
mnp_het_hom\t\
Insert_het_hom\t\
del_het_hom\t\
Indel_het_hom\t\
insert_del\t\
indel_snp_plus_mnp\t\
sample" > $path_add_stats

for sample in ERR696683 ERR7533{68..78}
do

  snp_trans=$(vcf_stats_filter "^SNP\sT")
  tot_het_hom=$(vcf_stats_filter "^Total\sH")
  snp_het_hom=$(vcf_stats_filter "^SNP\sH")
  mnp_het_hom=$(vcf_stats_filter "^MNP\sH")
  insert_het_hom=$(vcf_stats_filter "^Insertion\sH")
  del_het_hom=$(vcf_stats_filter "^Deletion\sH")
  indel_het_hom=$(vcf_stats_filter "^Indel\sH")
  insert_del=$(vcf_stats_filter "^Insertion\/Deletion")
  indel_snp_plus_mnp=$(vcf_stats_filter "^Indel\/SNP")

  # Values of the second VCF file (additional VCF statistics)
  echo -e "$snp_trans\t\
  $tot_het_hom\t\
  $snp_het_hom\t\
  $mnp_het_hom\t\
  $insert_het_hom\t\
  $del_het_hom\t\
  $indel_het_hom\t\
  $insert_del\t\
  $indel_snp_plus_mnp\t\
  $sample"

done >> $path_add_stats


