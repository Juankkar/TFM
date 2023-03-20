#!/usr/bin/env bash

## Obtenemos los runs de reports como un array

## ftp de las Lecturas 1
reads1=$(cat $1 \       ## $1 es el archivo report almacenado en ../data/report
| cut -f7 \
| cut -d ";" -f 1 \
| grep -v "fastq_ftp")

## ftp de las lecturas 2
reads2=$(cat $1 \
| cut -f7 \
| cut -d ";" -f 2 \
| grep -v "fastq_ftp")

## Los convertimos a un array
reads1_array=$(echo $reads1)
reads2_array=$(echo $reads2)

read -p "Ruta donde quiere las lecturas crudas: " rawdata

## Desacargamos las lecturas
for r1 in ${reads1_array[*]}
do
    wget -P $rawdata $r1
done

for r2 in ${reads2_array[*]}
do
    wget -P $rawdata $r2
done
