configfile: "config.yaml"

rule all:
    input:
        expand("results/fastqc_result/{sample}_fastqc.html", 
               sample=config["samples"]),
        expand("data/processed/fastp_processed/{sample}_fastp.fastq.gz",
               sample=config["samples"]),
        expand("results/mapped_reads/{sample}.sam", 
               sample=config["samples"])

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

## View the quality of the samples
rule fastqc:
    input: 
        "data/raw/{sample}.fastq.gz"
    output:
        "results/fastqc_result/{sample}_fastqc.html"
    shell:
        """
        fastqc {input} -o results/fastqc_result/
        """

## Pre-processed of the data
rule fastp:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "data/processed/fastp_processed/{sample}_fastp.fastq.gz"
    params:
        cut_tail=config["fastp_cuttail"],
        cut_front=config["fastp_cutfront"],
        cut_meanq=config["fastp_cutmeanq"]
    shell:
        """
        fastp -i {input} -o {output} \
        --cut_tail '{params.cut_tail}' \
        --cut_front '{params.cut_front}' \
        --cut_mean_quality '{params.cut_meanq}'
        mv *.json data/processed/fastp_processed
        mv fastp.html data/processed/fastp_processed
        """

## Creating sam files for forward and reverse reads
rule bwa_mapping:
    input:
        reference = "data/reference/genome.fa",
        files = get_bwa_map_input_fastqs
    output:
        "results/mapped_reads/{sample}.sam"
    shell:
        """
        bwa index {input.reference}
        bwa mem -a {input.reference} {input.files} -o {output} 2> info.out
        mv info.out results/mapped_reads/
        """

## script for joining the SAM files
rule merge_sam_files:
    shell:
        "bash code/2join_samfiles.sh"

## Transforming SAM to BAM and sorting
rule sam_to_bam:
    shell:
        "bash code/3sam_to_bam.sh"

## Delete duplicates
rule delete_duplicates:
    shell:
        "bash code/4delete_duplicates.sh"

## Extracting variants
rule extracting_variants:
    shell:
        "bash code/5extracting_variants.sh"