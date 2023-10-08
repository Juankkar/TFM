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

### This was the original code used to download the data for the publication -->

## This loop will dowload the bam files from the project to the
## "original_bam" directory
# for bam in $bam_files      
# do            # path + (echo = original BAM name from report) 
#     if [ ! -f ${path_origbam}$(echo $bam | cut -d "/" -f6) ]
#     then
#         wget -P ${path_origbam}/ $bam
# 
#     elif [ -f ${path_origbam}$(echo $bam | cut -d "/" -f6) ]
#     then
#         echo -n "==> THE FILE " 
#         echo -n "${path_origbam}$(echo $bam | cut -d "/" -f6) "
#         echo -e "ALREADY EXIST!!! <==\n"
# 
#     else
#         echo "==>>> ERROR <<<=="
#     fi
# done

#----------------------------------------------------------#
## GNU Parallel: "Im about to end this man's whole career"
## This will download all data using parallelization :/ 
## The most heavy file is 14 Gb (~93 Gb in total all files)
parallel wget -P ${path_origbam} ::: $bam_files 
#----------------------------------------------------------#

## Unfortunately there are problems with both methods, some of the BAM
## files can get corrupted, lets try to solve this problem

## ERROR check

raw_bams=$(cat $path_report \
        | cut -f 8 \
        | grep -v "submitted_ftp" \
        | cut -d "/" -f 6)

while true;do

    corrupted=0

    ## Indexing  the bam files, in case of corruption, the respective files would
    ## be downloaded again
    parallel samtools index ::: $(for bam in $raw_bams;do echo ${path_origbam}$bam;done)
    
    ## This checks if some of the BAM files are not able to be indexed -> corrupted
    ## in that case the corrupted will be downloaded again
    for run in $raw_bams;do
        if [[ ! -f $(echo ${path_origbam}$(echo ${run}.bai)) ]];then
            echo "==> The file $run is corrupted, downloading again <=="
            rm ${path_origbam}${run}
            wget -P $path_origbam $(echo $(cat $path_report \
                                            | cut -f 8 \
                                            | grep "$run"))
            corrupted=1
        fi
    done 

    rm $path_origbam*.bai

    if [[ $corrupted -eq 0 ]];then
        break
    fi

done
