#!/usr/bin/env bash

## ftp from the files bams -> add the path of the report.tsv of the project 
## when execute the script

## Extracting to a variable the bam files from report.tsv
bam_files=$(cat $1 \
    | cut -f8 \
    | grep -v "submitted_ftp")

## We transform the bam files list into an array
bam_array=$(echo $bam_files)

read -p "Ruta donde quiere las lecturas crudas: " rawdata


## This loop will dowload the bam files from the project to the
## "original_bam" directory
for bam in ${bam_array[*]}
do
    if [ ! -f ${rawdata}$(echo $bam \
        | cut -d "/" -f6) ]
    then
        wget -P $rawdata $bam
    elif [ -f ${rawdata}$(echo $bam \
        | cut -d "/" -f6) ]
    then
        echo "===>>> THE FILE ${rawdata}$(echo $bam \
        | cut -d "/" -f6) ALREADY EXIST!!! <<<==="
    else
        echo "==>>> ERROR <<<=="
    fi
done
