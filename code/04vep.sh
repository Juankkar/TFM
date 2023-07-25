#!/usr/bin/env bash

#################
## DIRECTORIES ##
#################

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
