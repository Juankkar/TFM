#!/usr/bin/env bash

if [[ ! -d results/biostatistics/joined_tables/ ]]
then
    mkdir results/biostatistics/joined_tables/
fi

for pattern in clin_sig biotype consequence polyphen pubmed clinvar_clnsig gene
do
    cat results/biostatistics/tables/*$pattern* \
        | awk "!/^${pattern}/ || NR == 1" \
        > results/biostatistics/joined_tables/${pattern}.tsv
done  