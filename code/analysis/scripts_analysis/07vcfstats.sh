#!/usr/bin/env bash

path=../../../metadata/logs/vcfstats/ 

echo -e "pass_fail_filter\tsnps\tmnps\tinsertions\tdeletions\tindels\tsame_reference\tsample" \
  > ../results_analysis/tables/vcf_basic_stats.tsv

for sample in ERR696683 ERR7533{68..78}
do

  passed_filters=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Pass" \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')

  failed_filters=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Fail" \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')

  snps=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^SNPs" \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')


  mnps=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^MNPs" \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')

  inser=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Insertions" \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')


  inser=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Insertions" \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')

  del=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Deletions" \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')

  indel=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Indels" \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')

  same_ref=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Same" \
    | cut -f 2 \
    | awk '{sum += $1} END {print sum}')

  echo -e "$passed_filters|$failed_filters\t$snps\t$mnps\t$inser\t$del\t$indel\t$same_ref\t$sample"

done >> ../results_analysis/tables/vcf_basic_stats.tsv 

echo -e "snp_transition_transversion\ttotal_het_hom\tsnp_het_hom\tmnp_het_hom\tInsert_het_hom\tdel_het_hom\tIndel_het_hom\tinsert_del\tindel_snp_plus_mnp\tsample" > ../results_analysis/tables/additional_vcf_stats.tsv 

for sample in ERR696683 ERR7533{68..78} 
do

  snp_trans=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^SNP\sT" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

  tot_het_hom=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Total\sH" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

  tot_het_hom=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Total\sH" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

  snp_het_hom=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^SNP\sH" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

  mnp_het_hom=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^MNP\sH" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

  insert_het_hom=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Insertion\sH" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

  del_het_hom=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Deletion\sH" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

  indel_het_hom=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Indel\sH" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}')

  insert_del=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Insertion\/Deletion" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}') 

  indel_snp_plus_mnp=$(sed -E 's/\s:\s/\t/g' ${path}${sample}* \
    | sed 's/:\s/\t/g' \
    | grep -E "^Indel\/SNP" \
    | cut -f 2 \
    | sed 's/)//' \
    | sed 's/\s(/\t/' \
    | sed 's/\//\t/' \
    | awk '{sum_col2 += $2; sum_col3 += $3} END {print sum_col2"|"sum_col3}') 


  echo -e "$snp_trans\t$tot_het_hom\t$snp_het_hom\t$mnp_het_hom\t$insert_het_hom\t$del_het_hom\t$indel_het_hom\t$insert_del\t$indel_snp_plus_mnp\t$sample"

done >> ../results_analysis/tables/additional_vcf_stats.tsv 


