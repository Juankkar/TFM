#!usr/bin/env bash

## Path to files
path_draw=../../../data/raw/
path_dpro=../../../data/processed/ 
path_result=../results_analysis/tables/fastp_analysis.tsv 

## Headers of the TSV file
echo -e "nseqs_before_fastp_1\t\
nseqs_after_fastp_1\t\
nseqs_before_fastp_2\t\
nseqs_after_fastp_2\t\
sample" > $path_result

for sample in ERR696683 ERR7533{68..78}
do
    before_fastp_1=$(zcat ${path_draw}${sample}*_1* | grep @ -c)
    after_fastp_1=$(zcat ${path_dpro}${sample}*_1* | grep @ -c)

    before_fastp_2=$(zcat ${path_draw}${sample}*_2* | grep @ -c)
    after_fastp_2=$(zcat ${path_dpro}${sample}*_2* | grep @ -c)

    # Values of the TSV files 
    echo -e "$before_fastp_1\t\
    $after_fastp_1\t\
    $before_fastp_2\t\
    $after_fastp_2\t\
    $sample"

done >> $path_result
