#!/usr/bin/env bash

## Path files
path_origbam=data/original_bam/ 
path_ob_filt=data/original_bam/filtering/ 

###################
##  DIRECTORIES  ##
###################

if [[ ! -d $path_ob_filt ]]
then
    mkdir $path_ob_filt
fi

#################
##  VARIABLES  ##
#################


original_bam_list=$(ls ${path_origbam}*.bam -1 | \
                    cut -d "/" -f 3 | \
                    cut -d "." -f 1)

#################
##  EXECUTION  ##
#################


for bam in $original_bam_list
do

    ## Indexing the BAM files
    samtools index ${path_origbam}${bam}.bam

    ## Chromosome filter ($1 = the chr)
    samtools view -b ${path_origbam}${bam}.bam $1 \
        > ${path_ob_filt}${bam}_${1}.bam

    ### Sorted filtered BAM file
    samtools sort -n ${path_ob_filt}${bam}_${1}.bam \
        -o ${path_ob_filt}${bam}_${1}_sorted.bam

    ### BAM to FASTQ
    samtools fastq -@ 8 ${path_ob_filt}${bam}_${1}_sorted.bam \
        -1 data/raw/${bam}_${1}_1.fastq.gz \
        -2 data/raw/${bam}_${1}_2.fastq.gz \
        -0 /dev/null -s /dev/null -n

done
