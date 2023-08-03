#!/usr/bin/env bash

for sample in ERR696683 ERR7533{68..78}
do
    cat ../../results/annotations/${sample}* \
        > ../../results/annotations/${sample}.txt
done  

rm ../../results/annotations/*fastq.gz.txt