#!/usr/bin/env bash

## Path to files
path_annot=../../../results/annotations/ 
 
for sample in ERR696683 ERR7533{68..78}
do
    # Joining DNA sequences from CHRs to their ERR sample
    cat ${path_annot}${sample}* \
        > ${path_annot}${sample}.txt
done  

rm ${path_annot}*fastq.gz.txt

for sample in ERR696683 ERR7533{68..78}
do
    ## Calculates the number of pb in each sequences
    while read -r line
    do 
        echo ${#line}
    done < ${path_annot}${sample}.txt \
        > ${path_annot}${sample}_len.txt 
done

rm ${path_annot}ERR696683.txt \
    ${path_annot}ERR7533{68..78}.txt