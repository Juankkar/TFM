#!/usr/bin/env bash

###############
## VARIABLES ##
###############

echo "Select the reads that you want to merge:"
read -p "Sample name for forward reads -> " forward
read -p "Sample name for reverse reads -> " reverse
read -p "Sample name for the merging file -> " merging

vector=($forward $reverse)

###############
## EXECUTION ##
###############

## first we will sort the files
for ends in ${vector[*]}
do
    samtools sort -n \
        results/mapped_reads/${ends}.sam \
        -o results/mapped_reads/${ends}_sorted.sam
done

## merging the files
samtools merge results/mapped_reads/${merging}.sam \
    results/mapped_reads/${forward}.sam \
    results/mapped_reads/${reverse}.sam

## sorting the SAM
samtools sort results/mapped_reads/${merging}.sam \
    -o results/mapped_reads/${merging}_sorted.sam

## Delete other sam files not joined to prevent excessive space usage
for ends in ${vector[*]}
do
    rm results/mapped_reads/${ends}.sam \
       results/mapped_reads/${ends}_sorted.sam
done