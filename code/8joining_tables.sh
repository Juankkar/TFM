#!/usr/bin/env bash

#################
## DIRECTORIES ##
#################

if [[ ! -d results/biostatistics/joined_tables/ ]]
then
    mkdir results/biostatistics/joined_tables/
fi

#################
##  EXECUTION  ##
#################

for pattern in biotype clin_sig clinvar_clnsig consequence location polyphen protein_position sift variant_class num_variants
do
    cat results/biostatistics/tables/*$pattern* \
        | awk "!/^${pattern}/ || NR == 1" \
        > results/biostatistics/joined_tables/${pattern}.tsv
done  