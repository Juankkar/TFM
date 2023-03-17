#!/usr/bin/env bash

## Obtenemos los runs de reports como un array
runs=$(cat $1 \
| cut -f4 | grep -v "^run")

echo ${runs[*]}

## Usamos fastq-dump para obtener los reads
for run in ${runs[*]}
do
    fastq-dump --outdir ../data/raw $run --no-verify-ssl
done

