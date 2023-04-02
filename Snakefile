configfile: "config.yaml"

rule all:
    input:
        expand("results/fastqc_result/{sample}_fastqc.html", sample=config["samples"]),
        expand("data/processed/fastp_processed/{sample}_fastp.fastq.gz", sample=config["samples"])

# def get_bwa_map_input_fastqs(wildcards):
#     return config["samples"][wildcards.sample]

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
        """

## Mapping the samples with the reference genome
rule bwa_mapping:
    input:
        reference = "data/reference/genome.fa",
        samp = "data/processed/fastp_processed/{sample}_fastp.fastq.gz"
    output:
        "results/mapped_reads/{sample}.bam"
    shell:
        """
        bwa mem {input.reference} {input.samp} | samtools view -Sb - > {output}
        """


