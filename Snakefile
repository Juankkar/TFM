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

##-----------------------------------##
##     SETUP FOR THE WORKFLOW        ##
##-----------------------------------##

## These first 3 rules are necessary to set everything ready

## 1 Downloading the data == "no_pain_no_gain"
rule download_data:
    input:
        script_download = "code/01dl_rawdata.bash",
        script_rename="code/02rename.py"
    output:
        touch("tasks/01download_data.done")
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        """
        bash {input.script_download} 

        python {input.script_rename} 
        """


## 2 Preprocessing the data (Executing this twice will generate an error)
rule pre_processing:
    input:
        script = "code/03extracting_fastq.sh"
    output:
        touch("tasks/02pre_processing.done")
    params:
        chr_choosed = config["chromosome"]
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        """
        bash {input.script} {params.chr_choosed}

        ## This is to prevent excessive space usage!!!
        rm data/original_bam/filtering/*
        """


## 3 Downloading the reference genome:
rule reference_genome:
    output:
        "data/reference/genome.fa"
    params:
        url = config["url_reference_genome"]
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        """
        wget -O {output}.gz {params.url}
        gzip -d {output}.gz 
        """

##-----------------------------------------##
##     REAL STARTING POINT WORKFLOW        ##
##-----------------------------------------##

## Now you can check if everything is ready to go using the command:
## snakemake -n 
## If this doesn't work, something bad happened

## 4 View the quality of the samples
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


## 5 Pre-processed of the data
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


## 6 View the quality of the trimmed samples
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


## 7 Creating sam files for forward and reverse reads
rule bwa_mapping:
    input:
        reference = "data/reference/genome.fa",
        files = get_bwa_map_input_fastqs
    output:
        "results/mapped_reads/{sample}.sam"
    log:
        "metadata/logs/sam/{sample}_infosam.out"
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        """
        bwa index {input.reference}
        bwa mem -a {input.reference} {input.files} \
            -o {output} \
            2> {log}
        """


## 8 script for joining the SAM files
rule merge_sam_files:
    input:
        script = "code/04join_samfiles.sh"
    output:
        touch("tasks/08merged_sam.done")
    params:
        # For example mine are _1 and _2, but could be _R1 _R2
        ends_1 = config["reads_forward_termination"],
        ends_2 = config["reads_reverse_termination"]
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        "bash {input.script} {params.ends_1} {params.ends_2}"


## 9 Transforming SAM to BAM and sorting
rule sam_to_bam:
    input:
        script = "code/05sam_to_bam.sh",
    output:
        touch("tasks/09sam_to_bam.done")
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        "bash {input.script}"


## 10 Deleting duplicates
rule delete_duplicates:
    input:
        script = "code/06delete_duplicates.sh",
    output:
        touch("tasks/10deleted_duplicates.done")
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        "bash {input.script}"


## 11 Extracting variants
rule extracting_variants:
    input:
        script = "code/07extracting_variants.sh"
    output:
        touch("tasks/11extracting_variants.done")
    params:
        ref_genome = config["ref_genome_name_file"],
        min_reads= config["min_reads_variant"]
    conda:
        "code/enviroments/Greference_tools.yml"
    shell:
        "bash {input.script} {params.ref_genome} {params.min_reads}"


## 12 Variant Effect Prediction DB
rule vep_install_db:
    output:
        touch("tasks/12vep_dependencies.done")
    params:
        species = config["vep_species"],
        assembly = config["vep_assembly"]
    conda: 
        "code/enviroments/vep.yml"
    shell:
        """
        vep_install -a cf \
            -s {params.species} \
            ls-y {params.assembly} \
            --ASSEMBLY {params.assembly}
        """


## 13 Running VEP in the command line
rule vep_cli:
    input:
        script = "code/08vep.sh"
    output:
        touch("tasks/13vep_cli.done")
    params:
        species = config["vep_species"],
        assembly=config["vep_assembly"]
    conda: 
        "code/enviroments/vep.yml"
    shell:
        "bash {input.script} {params.species} {params.assembly}"


## 14 Parsing data from the VCF files with R
rule parsing_dataR:
    input:
        script = "code/09parsing_vep_data.R"
    output:
        touch("tasks/14parsing_dataR.done")
    params:
        dir1 = "results/biostatistics/",
        dir2 = "results/biostatistics/tables",
        dir3 = "results/biostatistics/plots",
        dir4 = "results/biostatistics/joined_tables",
        gene_filter=config["gene_to_filterR"],
        chr_choosed=config["chromosome"],
        gene=config["gene_to_filterR"]
    conda:
        "code/enviroments/biostatisticsR.yml"
    shell:
        """
        for dir in {params.dir1} {params.dir2} {params.dir3} {params.dir4}
        do
            if [[ ! -d $dir ]]
            then
                mkdir $dir
            fi
        done

        Rscript {input.script} {params.gene_filter}

        ## Joining parsed tables for each sample
        cat results/biostatistics/tables/*{params.chr_choosed}* \
            | awk "!/^$(cat results/biostatistics/tables/*{params.chr_choosed}* \
            | cut -f 1 | head -n 1)/ || NR == 1" \
            > results/biostatistics/joined_tables/{params.gene}.tsv
        """


##---------------------------------------------------------##
##     FINAL BOSS, YOU NEED TO HAVE THE 5 GENE TABLES      ##
##---------------------------------------------------------##

## 15
rule R_plotting:
    input:
        script = "code/10final_plot.R"
    output:
        png1="results/biostatistics/plots/final_plot.png",
        png2="results/biostatistics/plots/other_plot.png"
    params:
        gene=config["gene_to_filterR"]
    conda:
        "code/enviroments/biostatisticsR.yml"
    shell:
        """
        Rscript {input.script} && echo "THE SCRIPT RAN WELL congrats :)"        

            # If this is red you didn't made it yet, maybe this will encourage you :)
            # https://www.youtube.com/watch?v=tYzMYcUty6s&ab_channel=TeamPsycosmos
        """
