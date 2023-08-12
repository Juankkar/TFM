#!/usr/bin/env bash

## Path with the files to use
path_logs=../../../metadata/logs/flagstats/
path_results=../results_analysis/tables/flagstats.tsv

#######################################
# Filter values from flagstat files
# Arguments:
#   Pattern to filter  
# Return
#   Numerical value of the stats    
#######################################
function flagstats_filter() {

  operation=$(cat ${path_logs}${sample}* \
    | grep -E $1 \
    | sed 's/ /\t/g' \
    | cut -f 1\
    | awk '{sum += $1} END {print sum}')

    echo $operation
}

echo -e "num_seqs_afterQC\t\
mapped_reads\t\
properly_mapped\t\
singletons\t\
sample" > $path_results

for sample in ERR696683 ERR7533{68..78}
do
  num_seqs=$(head -n 1 ${path_logs}${sample}* \
    | grep -E [[:digit:]] \
    | tr " " "\t" \
    | sed 's/\+//' \
    | cut -f 1 \
    | awk '{sum += $1} END {print sum}')

  mapped_reads=$(flagstats_filter "mapped\s\(")
  properly_paired=$(flagstats_filter "properly")
  singletons=$(flagstats_filter "singletons")

   echo -e "$num_seqs\t\
   $mapped_reads\t\
   $properly_paired\t\
   $singletons\t\
   $sample"

done >> $path_results

