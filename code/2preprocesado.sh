#!/usr/bin/env bash

## Comenzaremos filtrando los cromosomas que estemos
## intersados en estudiar

# samtools index C2466_GCCAAT.bwa.sorted.rmdup.recal.realigned.fixed.bam
# samtools view -b C2466_GCCAAT.bwa.sorted.rmdup.recal.realigned.fixed.bam chr{1,5,17} > filtering/IV_method1.bam
# samtools sort -n IV_method1.bam -o IV_method1_sorted.bam
# samtools fastq -@ 8 IV_method1.bam \
#     -1 ../../raw/IV_method1_1.fastq.gz \
#     -2 ../../raw/IV_method1_2.fastq.gz \
#     -0 /dev/null -s /dev/null -n

