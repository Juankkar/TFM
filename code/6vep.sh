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

vep -i results/variants/${sample}.vcf \
    --cache \
    --assembly ${assembly} \
    --offline \
    --species ${specie} \
    --force_overwrite \
    -o results/variants/vep/${sample}.txt \
    --vcf \
    --fork 4
    
