#!/usr/bin/env bash

#################
## DIRECTORIES ##
#################

if [[ ! -d results/variants/vep/ ]]
then
    mkdir results/variants/vep/
fi

#################
##  VARIABLES  ##
#################

read -p "Select your sample VCF file for VEP -> " sample
read -p "Select assembly (example: GRCh38) -> " assembly
read -p "Select your species (example: homo_sapiens) -> " specie 

#################
##  EXECUTION  ##
#################

# vep -i results/variants/${sample}.vcf \
#     --cache \
#     --assembly ${assembly} \
#     --offline \
#     --species ${specie} \
#     --force_overwrite \
#     -o results/variants/vep/${sample}.txt \
#     --vcf \
#     --fork 4

## No se si me va a dar tiempo de subir esto antes de que lo veas, pero creo que asi esta mas completo :)
## Supongo que sobraran cosas 
vep -i results/variants/${sample}.vcf \
    --offline \
    --force_overwrite \
    --assembly ${assembly} \
    --af \
    --af_1kg \
    --appris \
    --biotype \
    --buffer_size 500 \
    --check_existing \
    --distance 5000 \
    --filter_common \
    --mane \
    --polyphen b \
    --protein \
    --pubmed \
    --regulatory \
    --sift b \
    --species $specie \
    --symbol \
    --transcript_version \
    --tsl \
    --uniprot \
    --cache \
    -o results/variants/vep/${sample}.txt