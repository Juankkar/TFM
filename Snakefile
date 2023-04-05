configfile: "config.yaml"

rule all:
    input:
        expand("results/fastqc_result/{sample}_fastqc.html", sample=config["samples"]),
        expand("data/processed/fastp_processed/{sample}_fastp.fastq.gz", sample=config["samples"])

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
        "data/processed/fastp_processed/"
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

## La idea era crear una regla que creara un archivo SAM
## pero al usar lo que sé lo que pasa es que me crean dos
## archivos un sam con las lecturas forward y otra con las
## lecturas reverse, lo cual es un problema porque no sé
## como juntarlos en un único SAM, pinta bien...
rule bwa_mapping:
    input:
        reference = "data/reference/genome.fa",
        files = get_bwa_map_input_fastqs
    output:
        "results/mapped_reads/{sample}.sam"
    shell:
        """
        bwa index {input.reference}
        bwa mem -a {input.reference} {input.files} -o {output} 2> chr16.out
        mv chr16.out results/mapped_reads/
        """

## Transform SAM file into BAM


