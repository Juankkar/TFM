#!/usr/bin/env bash

#################
## DIRECTORIES ##
#################

if [[ ! -d metadata/logs/markduplicates ]]
then
    mkdir metadata/logs/markduplicates/
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
        --METRICS_FILE metadata/logs/markduplicates/${sample}_markDuplicatesMetrics.txt \
        --ASSUME_SORTED True

    ## Indexing the new BAM file generated
    samtools index results/mapped_reads/bam_files/${sample}_dedup.bam
done

