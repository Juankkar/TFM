#!/usr/bin/env bash

if [[ ! -d results/mapped_reads/bam_files/ ]]
then 
    mkdir results/mapped_reads/bam_files/
fi

read -p "Write rute of your sam file -> " sam
read -p "Write rute for your bam file -> " bam
read -p "Write rute for your sorted bam file -> " sorted_bam
read -p "Write rute for your flagstats file -> " flagstats_file

## From SAM to BAM
samtools view -bS $sam > $bam
## Sorting
samtools sort $bam > $sorted_bam
## Index
samtools index $sorted_bam
## Summary, basic statistics
samtools flagstats $sorted_bam > $flagstats_file