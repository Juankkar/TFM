#!/usr/bin/env bash

echo -e "num_seqs_afterQC\tmapped_reads\tproperly_mapped\tsingletons\tsample" \
  > ../results_analysis/tables/flagstats.tsv

for sample in ERR696683 ERR7533{68..78}
do
  num_seqs=$(head -n 1 ../../../metadata/logs/flagstats/${sample}* \
    | grep -E [[:digit:]] \
    | tr " " "\t" \
    | sed 's/\+//' \
    | cut -f 1 \
    | awk '{sum += $1} END {print sum}')

   mapped_reads=$(cat ../../../metadata/logs/flagstats/${sample}* \
    | grep "mapped (" \
    | sed 's/ /\t/g' \
    | cut -f 1\
    | awk '{sum += $1} END {print sum}')

  properly_paired=$(cat ../../../metadata/logs/flagstats/${sample}* \
   | grep "properly" \
   | sed 's/ /\t/g' \
   | cut -f 1 \
   | awk '{sum += $1} END {print sum}')

  singletons=$(cat ../../../metadata/logs/flagstats/${sample}* \
   | grep "singletons" \
   | sed 's/ /\t/g' \
   | cut -f 1 \
   | awk '{sum += $1} END {print sum}')

   echo -e "$num_seqs\t$mapped_reads\t$properly_paired\t$singletons\t$sample"

done >> ../results_analysis/tables/flagstats.tsv