#!/usr/bin/env bash

if [[ ! -d ../../../results/annotations/ ]]
then
    mkdir ../../../results/annotations/ 
fi

for reads in _1 _2
do
    ## $1 = chromosome that you want to use
    for file in $(echo $(ls -1 ../../../data/raw/ | grep -v README | grep $reads | grep $1)) 
        do
            zcat ../../../data/raw/$file | sed -n '2~4p' \
                > ../../../results/annotations/${file}.txt 
        done
done
