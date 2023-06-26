#!/usr/bin/env bash

###################
##  DIRECTORIES  ##
###################

if [[ ! -d data/original_bam/filtering ]]
then
    mkdir data/original_bam/filtering
fi

#################
##  VARIABLES  ##
#################


original_bam_list=$(ls data/original_bam/*.bam -1 | \
                    cut -d "/" -f 3 | \
                    cut -d "." -f 1)

#################
##  EXECUTION  ##
#################


for bam in ${original_bam_list[*]}
do

    ## Indexing the BAM files
    samtools index data/original_bam/${bam}.bam

    ## Chromosome filter
    samtools view -b data/original_bam/${bam}.bam $1 \
        > data/original_bam/filtering/${bam}_${1}.bam

    ### Sorted filtered BAM file
    samtools sort -n data/original_bam/filtering/${bam}_${1}.bam \
        -o data/original_bam/filtering/${bam}_${1}_sorted.bam

    ### BAM to FASTQ
    samtools fastq -@ 8 data/original_bam/filtering/${bam}_${1}_sorted.bam \
        -1 data/raw/${bam}_${1}_1.fastq.gz \
        -2 data/raw/${bam}_${1}_2.fastq.gz \
        -0 /dev/null -s /dev/null -n

done
