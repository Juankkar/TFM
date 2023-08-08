#!/usr/bin/env bash

function vcf_basicstats_filter(){

  operation=$()$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E $1 \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')

      echo $operation
}

function vcf_stats_filter() {

  operation=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E $1 \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

    echo $operation
}


path=../../../metadata/logs/vcfstats/ 

echo -e "pass_fail_filter\tsnps\tmnps\tinsertions\tdeletions\tindels\tsame_reference\tsample" \
  > ../results_analysis/tables/vcf_basic_stats.tsv

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

  echo -e "$passed_filters|$failed_filters\t$snps\t$mnps\t$inser\t$del\t$indel\t$same_ref\t$sample"

done >> ../results_analysis/tables/vcf_basic_stats.tsv 

echo -e "snp_transition_transversion\ttotal_het_hom\tsnp_het_hom\tmnp_het_hom\tInsert_het_hom\tdel_het_hom\tIndel_het_hom\tinsert_del\tindel_snp_plus_mnp\tsample" > ../results_analysis/tables/additional_vcf_stats.tsv 

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

  echo -e "$snp_trans\t$tot_het_hom\t$snp_het_hom\t$mnp_het_hom\t$insert_het_hom\t$del_het_hom\t$indel_het_hom\t$insert_del\t$indel_snp_plus_mnp\t$sample"

done >> ../results_analysis/tables/additional_vcf_stats.tsv 


