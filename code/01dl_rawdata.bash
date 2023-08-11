#!/usr/bin/env bash

## ftp from the files bams -> add the path of the report.tsv of the project 
## when execute the script

path_origbam=data/original_bam/ 
path_report=metadata/report.tsv 

#################
## DIRECTORIES ##
#################

if [[ ! -d $path_origbam ]]
then
    mkdir $path_origbam
fi

#################
##  VARIABLES  ##
#################

## Extracting to a variable the bam files from report.tsv
bam_files=$(cat $path_report \
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
    if [ ! -f ${path_origbam}$(echo $bam \
        | cut -d "/" -f6) ]
    then
        wget -P ${path_origbam}/ $bam
    elif [ -f ${path_origbam}$(echo $bam \
        | cut -d "/" -f6) ]
    then
        echo "===>>> THE FILE \
              ${path_origbam}$(echo $bam | cut -d "/" -f6) \
              ALREADY EXIST!!! <<<==="
    else
        echo "==>>> ERROR <<<=="
    fi
done
