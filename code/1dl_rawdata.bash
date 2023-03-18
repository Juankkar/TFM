#!/usr/bin/env bash

## Obtenemos los runs de reports como un array

## ftp de las Lecturas 1
reads1=$(cat $1 \
| cut -f7 \
| cut -d ";" -f 1 \
| tail -n 1)

## ftp de las lecturas 2
reads2=$(cat $1 \
| cut -f7 \
| cut -d ";" -f 2 \
| tail -n 1)

## Los convertimos a un array
reads1_array=$(echo $reads1)
reads2_array=$(echo $reads2)

## Desacargamos las lecturas
for r1 in ${reads1_array[*]}
do
    for r2 in ${reads2_array[*]}
    do
        wget -P ../raw/ $r1
        wget -P ../raw/ $r2
    done
done

