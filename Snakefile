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
        "code/01dl_rawdata.bash"

## View the quality of the samples
rule fastqc:
    input: 
        "data/raw/{sample}.fastq.gz"
    output:
        "results/fastqc_result/{sample}_fastqc.html"
    conda:
        "code/enviroments/Greference_tools.yml"
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
        "code/enviroments/Greference_tools.yml"
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
        "code/enviroments/Greference_tools.yml"
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
        "code/enviroments/Greference_tools.yml"
    shell:
        """
        bwa index {input.reference}
        bwa mem -a {input.reference} {input.files} -o {output} 2> info.out
        mv info.out results/mapped_reads/
        """

## script for joining the SAM files
rule merge_sam_files:
    input:
        script = "code/02join_samfiles.sh" 
    params:
        ## For example mine are _1 and _2, but could be _R1 _R2
        ends_1 = config["reads_fordward_termination"],
        ends_2 = config["reads_reversed_termination"]
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        "bash {input.script} {params.ends_1} {params.ends_2}"

## Transforming SAM to BAM and sorting
rule sam_to_bam:
    input:
        script = "code/03sam_to_bam.sh"
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        "bash {input.script}"

## Delete duplicates
rule delete_duplicates:
    input:
        script = "code/04delete_duplicates.sh"
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        "bash {input.script}"

## Extracting variants
rule extracting_variants:
    input:
        script = "code/05extracting_variants.sh" 
    params:
        ref_genome = config["ref_genome_name_file"],
        min_reads= config["min_reads_variant"]
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        "bash {input.script} {params.ref_genome} {params.min_reads}"

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

## Running VEP in the command line
rule vep_cli:
    input:
        script = "code/06vep.sh"
    params:
        species = config["vep_species"],
        assembly=config["vep_assembly"]
    conda: 
        "code/enviroments/vep.yml"
    shell:
        "bash {input.script} {params.species} {params.assembly}"

## Doing some biostatistics in R
rule biostatisticsR_tables:
    params:
        dir1 = "results/biostatistics/",
        dir2 = "results/biostatistics/tables",
        dir3 = "results/biostatistics/plots",
    conda:
        "code/enviroments/biostatisticsR.yml"
    shell:
        """
        for dir in {params.dir1} {params.dir2} {params.dir3}
        do
            if [[ ! -d $dir ]]
            then
                mkdir $dir
            fi
        done

        Rscript code/07biostatistics_tables.R
        """

## Joining tables
rule joining_tables:
    input:
        script = "code/08joining_tables.sh" 
    conda:
        ## It can be any of them for this one really
        "code/enviroments/Greference_tools.yml"
    shell:
        """
        bash {input.script}
        """

## We will plot the data
rule R_ploting:
    conda:
        "code/enviroments/biostatisticsR.yml"
    shell:
        """
        Rscript code/09ploting.R
        """
