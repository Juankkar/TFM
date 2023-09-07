#!/usr/bin/env bash

## ftp from the files bams -> add the path of the report.tsv 
## of the project when execute the script

## Path files
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

## Extracting to a variable the bam ftps from report.tsv
bam_files=$(cat $path_report \
    | cut -f 8 \
    | grep -v "submitted_ftp")

#################
##  EXECUTION  ##
#################

## This loop will dowload the bam files from the project to the
## "original_bam" directory
for bam in $bam_files      
do            # path + (echo = original BAM name from report) 
    if [ ! -f ${path_origbam}$(echo $bam | cut -d "/" -f6) ]
    then
        wget -P ${path_origbam}/ $bam

    elif [ -f ${path_origbam}$(echo $bam | cut -d "/" -f6) ]
    then
        echo -n "==> THE FILE " 
        echo -n "${path_origbam}$(echo $bam | cut -d "/" -f6) "
        echo -e "ALREADY EXIST!!! <==\n"

    else
        echo "==>>> ERROR <<<=="
    fi
done
