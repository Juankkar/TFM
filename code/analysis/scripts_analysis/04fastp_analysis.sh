#!usr/bin/env bash

echo -e "nseqs_before_fastp_1\tnseqs_after_fastp_1\tnseqs_before_fastp_2\tnseqs_after_fastp_2\tsample" \
    > ../results_analysis/tables/fastp_analysis.tsv

for sample in ERR696683 ERR7533{68..78}
do
    before_fastp_1=$(zcat ../../../data/raw/${sample}*_1* | grep @ -c)
    after_fastp_1=$(zcat ../../../data/processed/${sample}*_1* | grep @ -c)

    before_fastp_2=$(zcat ../../../data/raw/${sample}*_2* | grep @ -c)
    after_fastp_2=$(zcat ../../../data/processed/${sample}*_2* | grep @ -c)

    echo -e "$before_fastp_1\t$after_fastp_1\t$before_fastp_2\t$after_fastp_2\t$sample"
done >> ../results_analysis/tables/fastp_analysis.tsv
