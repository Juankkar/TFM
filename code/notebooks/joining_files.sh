#!/usr/bin/env bash

for sample in ERR696683 ERR7533{68..78}
do
    cat ../../results/annotations/${sample}* \
        > ../../results/annotations/${sample}.txt
done  

rm ../../results/annotations/*fastq.gz.txt

for sample in ERR696683 ERR7533{68..78}
do
    while read -r line
    do 
        echo ${#line}
    done < ../../results/annotations/${sample}.txt > ../../results/annotations/${sample}_len.txt 
done

rm ../../results/annotations/ERR696683.txt \
    ../../results/annotations/ERR7533{68..78}.txt