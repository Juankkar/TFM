configfile: "config.yaml"

rule all:
    input:
        expand("results/fastqc_result/{sample}_fastqc.html", 
               sample=config["samples"]),
        expand("data/processed/fastp_processed/{sample}_fastp.fastq.gz",
               sample=config["samples"]),
        expand("results/fastqc_result/trimmed/{sample}_fastp_fastqc.html", 
               sample=config["samples"]),
        expand("results/mapped_reads/{sample}.sam", 
               sample=config["samples"])

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

## Downloading the data
rule download_data:
    shell:
        "code/1dl_rawdata.bash"

## View the quality of the samples
rule fastqc:
    input: 
        "data/raw/{sample}.fastq.gz"
    output:
        "results/fastqc_result/{sample}_fastqc.html"
    conda:
        "code/enviroments/TFM.yml"
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
        cut_meanq=config["fastp_cutmeanq"],
        length=config["fastp_length"]
    conda:
        "code/enviroments/TFM.yml"
    shell:
        """
        fastp -i {input} -o {output} \
        --cut_tail '{params.cut_tail}' \
        --cut_front '{params.cut_front}' \
        --cut_mean_quality '{params.cut_meanq}' \
        -l {params.length}
        mv *.json data/processed/fastp_processed
        mv fastp.html data/processed/fastp_processed
        """

## View the quality of the trimmed samples
rule fastqc_trimmed:
    input: 
        "data/processed/fastp_processed/{sample}_fastp.fastq.gz"
    output:
        "results/fastqc_result/trimmed/{sample}_fastp_fastqc.html"
    conda:
        "code/enviroments/TFM.yml"
    shell:
        """
        fastqc {input} -o results/fastqc_result/trimmed/
        """

## Creating sam files for forward and reverse reads
rule bwa_mapping:
    input:
        reference = "data/reference/genome.fa",
        files = get_bwa_map_input_fastqs
    output:
        "results/mapped_reads/{sample}.sam"
    conda:
        "code/enviroments/TFM.yml"
    shell:
        """
        bwa index {input.reference}
        bwa mem -a {input.reference} {input.files} -o {output} 2> info.out
        mv info.out results/mapped_reads/
        """

## script for joining the SAM files
rule merge_sam_files:
    conda:
        "code/enviroments/TFM.yml"
    shell:
        "bash code/2join_samfiles.sh"

## Transforming SAM to BAM and sorting
rule sam_to_bam:
    conda:
        "code/enviroments/TFM.yml"
    shell:
        "bash code/3sam_to_bam.sh"

## Delete duplicates
rule delete_duplicates:
    conda:
        "code/enviroments/TFM.yml"
    shell:
        "bash code/4delete_duplicates.sh"

## Extracting variants
rule extracting_variants:
    conda:
        "code/enviroments/TFM.yml"
    shell:
        "bash code/5extracting_variants.sh"

## Variant Effect Prediction DB
rule vep:
    params:
        species = config["vep_species"],
        assembly = config["vep_assembly"]
    conda: 
        "code/enviroments/vep.yml"
    shell:
        """
        vep_install -a cf \
            -s {params.species} \
            ls-y {params.assembly}
        """

## Running VEP
rule vep_cli:
    conda: 
        "code/enviroments/vep.yml"
    shell:
        "bash code/6vep.sh"

## Doing some biostatistics in R
rule biostatisticsR_tables:
    params:
        dir1 = "results/biostatistics/",
        dir2 = "results/biostatistics/tables",
        dir3 = "results/biostatistics/plots"
    conda:
        "code/enviroments/biostatistics.yml"
    shell:
        """
        for dir in {params.dir1} {params.dir2} {params.dir3}
        do
            if [[ ! -d $dir ]]
            then
                mkdir $dir
            fi
        done

        Rscript code/7biostatistics_tables.R
        """

## Joining tables
rule joining_tables:
    conda:
        ## It can be any of them for this one really
        "code/enviroments/TFM.yml"
    shell:
        """
        bash code/8joining_tables.sh     
        """

## Finally we will plot the data
rule R_ploting:
    conda:
        "code/enviroments/biostatistics.yml"
    shell:
        """
        Rscript code/9ploting.R
        """


        
