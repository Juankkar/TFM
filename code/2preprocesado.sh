#!/usr/bin/env bash

## Comenzaremos filtrando los cromosomas que estemos
## intersados en estudiar

read -p "> Path/Name filter BAM -> " filter_bam
# read -p "> Path/Name sorted BAM? -> " sorted_bam
# read -p "> Path/Name FASTQ file? -> " fastq_name
 
samtools view $1 \
    | awk '{if($3 == chr1 && $3 == ch5)}{print $0}' \
    | samtools view -Sb > $filter_bam

 