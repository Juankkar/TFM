#!/usr/bin/env bash

## Patth to files
path_annot=../../../results/annotations/ 
path_draw=../../../data/raw/ 

if [[ ! -d $path_annot ]]
then
    mkdir $path_annot 
fi

for reads in _1 _2
do
    # $1 = chromosome that you want to use
    for file in $(echo $(ls -1 $path_draw \
                            | grep -v README \
                            | grep $reads \
                            | grep $1)) 
    do
    # Extracting all DNA sequences from the reads and passing them to a TXT file
        zcat ${path_draw}${file} \
            | sed -n '2~4p' \
            > ${path_annot}${file}.txt 
    done
done
