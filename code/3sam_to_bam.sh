#!/usr/bin/env bash

#################
## DIRECTORIES ##
#################

if [[ ! -d results/mapped_reads/bam_files/ || 
        ! -d results/mapped_reads/bam_files/info/ ]]
then 
    mkdir results/mapped_reads/bam_files/
    mkdir results/mapped_reads/bam_files/info/
fi

#################
##  VARIABLES  ##
#################

read -p "Write the name of your sample -> " sample

#################
##  EXECUTION  ##
#################

## From SAM to BAM
samtools view -bS \
    results/mapped_reads/${sample}.sam \
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
    > results/mapped_reads/bam_files/info/${sample}.flagstats