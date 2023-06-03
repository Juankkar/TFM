#!/usr/bin/env bash

## ftp from the files bams -> add the path of the report.tsv of the project 
## when execute the script

#################
## DIRECTORIES ##
#################

if [[ ! -d data/original_bam/ ]]
then
    mkdir data/original_bam/
fi

#################
##  VARIABLES  ##
#################

## Extracting to a variable the bam files from report.tsv
bam_files=$(cat metadata/report.tsv \
    | cut -f 8 \
    | grep -v "submitted_ftp")

## We transform the bam files list into an array
bam_array=$(echo $bam_files)

#################
##  EXECUTION  ##
#################

## This loop will dowload the bam files from the project to the
## "original_bam" directory
for bam in ${bam_array[*]}
do
    if [ ! -f data/original_bam/$(echo $bam \
        | cut -d "/" -f6) ]
    then
        wget -P data/original_bam/ $bam
    elif [ -f data/original_bam/$(echo $bam \
        | cut -d "/" -f6) ]
    then
        echo "===>>> THE FILE data/original_bam/$(echo $bam \
        | cut -d "/" -f6) ALREADY EXIST!!! <<<==="
    else
        echo "==>>> ERROR <<<=="
    fi
done
