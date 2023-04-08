#!/usr/bin/env bash

##################
## DIREECTORIES ##
##################

if [[ ! -d results/variants/ ||
      ! -d results/variants/stats ]]
then
    mkdir results/variants/
    mkdir results/variants/stats
fi

##################
##  VARIABLES   ##
##################

read -p "Write the reference genome name -> " reference
read -p "Write the sample name -> " sample
read -p "Min reads that has to be mapping to considerate a variant -> " n_C

##################
##  EXECUTION   ##
##################

## Calling variants
freebayes \
    -C $n_C \
    -f data/reference/$reference \
    results/mapped_reads/bam_files/${sample}_dedup.bam \
    > results/variants/${sample}.vcf

## Vemos las estadísticas sobre las variantes encontradas
rtg /results/variants/${sample}.vcf \
    > results/variants/stats/${sample}.vcfstats

# Filtering 1; only indels
vcftools --vcf results/variants/${sample}.vcf \
    --keep-only-indels \
    --recode \
    --recode-INFO-all \
    --out results/variants/${sample}_indels.vcf

# Filtering 2; only SNPs
vcftools --vcf results/variants/${sample}.vcf \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --out results/variants/${sample}_snvs.vcf