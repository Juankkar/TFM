#!/usr/bin/env bash

##################
## DIREECTORIES ##
##################

if [[ ! -d results/variants/ ]]
then
    mkdir results/variants/
fi

##################
##  VARIABLES   ##
##################

read -p "Write the reference genome (including file termination, example: .fa) -> " reference
read -p "Write the sample name -> " sample
read -p "Min reads that has to be mapping to considerate a variant -> " n_C

##################
##  EXECUTION   ##
##################

## Calling variants
freebayes \
    -C $n_C \
    -f data/reference/$reference \
    results/mapped_reads/bam_files/${sample}_dedup.bam \
    > results/variants/${sample}.vcf
