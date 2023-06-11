#!/usr/bin/env bash

###############
## VARIABLES ##
###############

sample_list=$(echo $(cat your_sample_list.txt))  

###############
## EXECUTION ##
###############
for sample in ${sample_list[*]}
do

    ## first we will sort the files
    samtools sort -n \
        results/mapped_reads/${sample}${1}.sam \
        -o results/mapped_reads/${sample}${1}_sorted.sam

    samtools sort -n \
        results/mapped_reads/${sample}${2}.sam \
        -o results/mapped_reads/${sample}${2}_sorted.sam

    ## merging the files
    samtools merge results/mapped_reads/${sample}.sam \
        results/mapped_reads/${sample}${1}_sorted.sam \
        results/mapped_reads/${sample}${2}_sorted.sam

    ## sorting the SAM
    samtools sort results/mapped_reads/${sample}.sam \
        -o results/mapped_reads/${sample}_sorted.sam
  
done

## Delete the samples not joined
rm results/mapped_reads/*${1}.sam
rm results/mapped_reads/*${2}.sam

rm results/mapped_reads/*${1}_sorted.sam
rm results/mapped_reads/*${2}_sorted.sam