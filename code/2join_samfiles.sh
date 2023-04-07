#!/usr/bin/env bash

echo "Select the reads that you want to merge:"
read -p "Write the rute for fowrward reads -> " forward
read -p "Write the rute for fowrward sorted reads ->" forward_sorted
read -p "Write the rute for reverse reads -> " reverse
read -p "Write the rute for reverse reads_sorted -> " reverse_sorted
read -p "Write the rute for the merging file -> " merging
read -p "Write the rute for the merging_sorted file -> " merging_sorted

## first we will sort the files
samtools sort -n $forward -o $forward_sorted
samtools sort -n $reverse -o $reverse_sorted

## merging the files
samtools merge $merging $forward_sorted $reverse_sorted

## sorting the SAM
samtools sort $merging -o $merging_sorted

## Delete other sam files not joined to prevent excessive space usage
rm $forward $forward_sorted $reverse $reverse_sorted