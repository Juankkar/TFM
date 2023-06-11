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

sample_list=$(echo $(cat your_sample_list.txt))

##################
##  EXECUTION   ##
##################
for sample in ${sample_list[*]}
do
    ## Calling variants
    freebayes \
        -C $2 \
        -f data/reference/$1 \
        results/mapped_reads/bam_files/${sample}_dedup.bam \
        > results/variants/${sample}.vcf
done
