#!/usr/bin/env bash

#################
## DIRECTORIES ##
#################

if [[ ! -d results/mapped_reads/bam_files/markduplicates ]]
then
    mkdir results/mapped_reads/bam_files/markduplicates/
fi

#################
##  VARIABLES  ##
#################

sample_list=$(echo $(cat your_sample_list.txt))  

#################
##  EXECUTION  ##
#################

for sample in ${sample_list[*]}
do
    ## Mark duplicates
    picard MarkDuplicates --INPUT results/mapped_reads/bam_files/${sample}_sorted.bam \
        --OUTPUT results/mapped_reads/bam_files/${sample}_dedup.bam \
        --METRICS_FILE results/mapped_reads/bam_files/markduplicates/${sample}_markDuplicatesMetrics.txt \
        --ASSUME_SORTED True

    ## Indexamos de nuevo el bam de duplicados en este caso
    samtools index results/mapped_reads/bam_files/${sample}_dedup.bam
done

