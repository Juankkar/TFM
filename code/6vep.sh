#!/usr/bin/env bash

#################
## DIRECTORIES ##
#################

if [[ ! -d results/variants/vep/ ]]
then
    mkdir results/variants/vep/
fi

if [[ ! -f data/ClinVar/clinvar.vcf.gz ]]
then
    # Compressed VCF file
    wget -P data/ClinVar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
fi

if [[ ! -f data/ClinVar/clinvar.vcf.gz.tbi ]]
then
    # Index file
    wget -P data/ClinVar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
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
    --offline \
    --force_overwrite \
    --assembly ${assembly} \
    --appris \
    --biotype \
    --variant_class \
    --check_existing \
    --filter_common \
    --mane \
    --polyphen b \
    --pubmed \
    --regulatory \
    --sift b \
    --species $specie \
    --symbol \
    --transcript_version \
    --tsl \
    --cache \
    --tab \
    -o results/variants/vep/${sample}.txt \
    --custom data/ClinVar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN