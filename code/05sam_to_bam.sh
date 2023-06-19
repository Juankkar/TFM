#!/usr/bin/env bash

#################
## DIRECTORIES ##
#################

for dir in results/mapped_reads/bam_files/ metadata/logs/flagstats/
do 
    if [[ ! -d $dir ]]
    then 
        mkdir $dir 
    fi
done

#################
##  VARIABLES  ##
#################

sample_list=$(echo $(cat your_sample_list.txt))  

#################
##  EXECUTION  ##
#################

for sample in ${sample_list[*]}
do
    ## From SAM to BAM
    samtools view -bS \
        results/mapped_reads/${sample}_sorted.sam \
        > results/mapped_reads/bam_files/${sample}.bam

    ## Sorting
    samtools sort \
        results/mapped_reads/bam_files/${sample}.bam \
        > results/mapped_reads/bam_files/${sample}_sorted.bam

    ## Index
    samtools index results/mapped_reads/bam_files/${sample}_sorted.bam

    ## Summary, basic statistics
    samtools flagstats \
        results/mapped_reads/bam_files/${sample}_sorted.bam \
        > metadata/logs/flagstats/${sample}.flagstats
done