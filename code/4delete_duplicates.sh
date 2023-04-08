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

read -p "Write the name fo your sample -> " sample

#################
##  EXECUTION  ##
#################

## Mark duplicates
picard MarkDuplicates --INPUT results/mapped_reads/bam_files/${sample}_sorted.bam \
    --OUTPUT results/mapped_reads/bam_files/${sample}_dedup.bam \
    --METRICS_FILE results/mapped_reads/bam_files/markduplicates/markDuplicatesMetrics.txt \
    --ASSUME_SORTED True

# picard AddOrReplaceReadGroups -I results/mapped_reads/bam_files/${sample}_dedup.bam \
#     -O results/mapped_reads/bam_files/${sample}_dedup_RG.bam \
#     -RGLB lib1 \
#     -RGPU unit1 \
#     -RGID M02899 \
#     -RGPL ILLUMINA \

# ## Indexamos de nuevo el bam de duplicados en este caso
samtools index results/mapped_reads/bam_files/${sample}_dedup.bam

