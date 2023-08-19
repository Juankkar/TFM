configfile: "config.yaml"

rule all:
    input:
        expand("data/original_bam/filtering/{sample}_chr7_sorted.bam",
               sample=config["samples"]),
        expand("results/fastqc_result/{sample}_1_fastqc.html", 
               sample=config["samples"]),
        expand("results/fastqc_result/{sample}_2_fastqc.html", 
               sample=config["samples"]),
        expand("data/processed/{sample}_1_fastp.fastq.gz",
               sample=config["samples"]),
        expand("data/processed/{sample}_2_fastp.fastq.gz",
               sample=config["samples"]),
        expand("results/fastqc_result/trimmed/{sample}_1_fastp_fastqc.html", 
               sample=config["samples"]),
        expand("results/fastqc_result/trimmed/{sample}_2_fastp_fastqc.html", 
               sample=config["samples"]),            
        expand("results/mapped_reads/{sample}.sam", 
               sample=config["samples"]),
        expand("results/mapped_reads/{sample}_sorted.sam",
               sample=config["samples"]),
        expand("results/mapped_reads/bam_files/{sample}.bam",
               sample=config["samples"]),
        expand("results/mapped_reads/bam_files/{sample}_sorted.bam",
               sample=config["samples"]),
        expand("results/mapped_reads/bam_files/{sample}_dedup.bam",
               sample=config["samples"]),
        expand("results/variants/{sample}.vcf",
               sample=config["samples"]),
        expand("results/variants/vep/{sample}.txt",
               sample=config["samples"])


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
        "code/environments/Greference_tools.yml"
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
        ## choose one sample, if you execute either sample, it will generate all results
        "data/original_bam/filtering/{sample}_sorted.bam"
    params:
        chr_choosed = config["chromosome"]
    conda:
        "code/environments/Greference_tools.yml"
    shell:
        """
        bash {input.script} {params.chr_choosed}
        """


## 3 Downloading the reference genome:
# if this is your 2º run, you might need to use the option --force
rule reference_genome:
    output:
        "data/reference/genome.fa"
    params:
        url = config["url_reference_genome"]
    conda:
        "code/environments/Greference_tools.yml"
    shell:
        """
        rm {output}* || \
            echo "==> If you hadn't a previous reference in the directory and there is an ERROR, is NORMAL <==\n"

        wget -O {output}.gz {params.url}
        gzip -d {output}.gz 
        """

##-----------------------------------------##
##     REAL STARTING POINT WORKFLOW        ##
##-----------------------------------------##

## 4 View the quality of the samples
rule fastqc:
    input: 
        read1 = "data/raw/{sample}_1.fastq.gz",
        read2 = "data/raw/{sample}_2.fastq.gz"
    output:
        read1 = "results/fastqc_result/{sample}_1_fastqc.html",
        read2 = "results/fastqc_result/{sample}_2_fastqc.html"
    conda:
        "code/environments/Greference_tools.yml"
    shell:
        """
        for read in {input.read1} {input.read2}
        do 
            fastqc $read -o results/fastqc_result/
        done
        """


## 5 Pre-processed of the data
rule fastp:
    input:
        read1 = "data/raw/{sample}_1.fastq.gz",
        read2 = "data/raw/{sample}_2.fastq.gz"
    output:
        read1 = "data/processed/{sample}_1_fastp.fastq.gz",
        read2 = "data/processed/{sample}_2_fastp.fastq.gz"
    params:
        cut_tail=config["fastp_cuttail"],
        cut_front=config["fastp_cutfront"],
        cut_meanq=config["fastp_cutmeanq"],
        length=config["fastp_length"]
    conda:
        "code/environments/Greference_tools.yml"
    shell:
        """
        fastp -i {input.read1} -I {input.read2} \
            -o {output.read1} -O {output.read2} \
            --cut_tail '{params.cut_tail}' \
            --cut_front '{params.cut_front}' \
            --cut_mean_quality '{params.cut_meanq}' \
            -l {params.length}
            mv *.json data/processed/ 
            mv fastp.html data/processed
        """


## 6 View the quality of the trimmed samples
rule fastqc_trimmed:
    input: 
        read1 = "data/processed/{sample}_1_fastp.fastq.gz",
        read2 = "data/processed/{sample}_2_fastp.fastq.gz"
    output:
        read1 = "results/fastqc_result/trimmed/{sample}_1_fastp_fastqc.html",
        read2 = "results/fastqc_result/trimmed/{sample}_2_fastp_fastqc.html"
    conda:
        "code/environments/Greference_tools.yml"
    shell:
        """
        for read in {input.read1} {input.read2}
        do
        fastqc $read -o results/fastqc_result/trimmed/
        done
        """


## 7 Creating sam files for forward and reverse reads
rule bwa_mapping:
    input:
        reference = "data/reference/genome.fa",
        read1 = "data/processed/{sample}_1_fastp.fastq.gz",
        read2 = "data/processed/{sample}_2_fastp.fastq.gz"
    output:
        sam = "results/mapped_reads/{sample}.sam",
        sam_sorted = "results/mapped_reads/{sample}_sorted.sam"
    log:
        "metadata/logs/sam/{sample}_infosam.out"
    conda:
        "code/environments/Greference_tools.yml"
    shell:
        """
        ## Mapping
        bwa index {input.reference}
        bwa mem -a {input.reference} \
        {input.read1} {input.read2} \
            -o {output.sam} \
            2> {log}
        
        ## sorting the SAM files
        samtools sort {output.sam} \
            -o {output.sam_sorted} 
        """


# 8 Transforming SAM to BAM and sorting
rule sam_to_bam:
    input:
        "results/mapped_reads/{sample}_sorted.sam",
    output:
       bam = "results/mapped_reads/bam_files/{sample}.bam",
       bam_sorted = "results/mapped_reads/bam_files/{sample}_sorted.bam"
    log:
        "metadata/logs/flagstats/{sample}.flagstats"
    conda:
        "code/environments/Greference_tools.yml"
    shell:
        """
        ## From SAM to BAM
        samtools view -bS \
            {input} \
            > {output.bam}

        ## Sorting
        samtools sort \
            {output.bam} \
            > {output.bam_sorted}

        ## Index
        samtools index {output.bam_sorted}

        ## Summary, basic statistics
        samtools flagstat \
            {output.bam_sorted} \
            > {log}
        """


## 9 Deleting duplicates
rule delete_duplicates:
    input:
        "results/mapped_reads/bam_files/{sample}_sorted.bam"
    output: 
        "results/mapped_reads/bam_files/{sample}_dedup.bam"
    conda:
        "code/environments/Greference_tools.yml"
    log:
        "metadata/logs/markduplicates/{sample}_markDuplicatesMetrics.txt"
    shell:
        """
        ## Mark duplicates
        picard MarkDuplicates --INPUT {input} \
            --OUTPUT {output} \
            --METRICS_FILE {log} \
            --ASSUME_SORTED True
    
        ## Indexing the new BAM file generated
        samtools index {output}
        """


## 10 Extracting variants
rule extracting_variants:
    input:
        reference = "data/reference/genome.fa",
        bam = "results/mapped_reads/bam_files/{sample}_dedup.bam" 
    output:
        "results/variants/{sample}.vcf"
    params:
        ref_genome = config["ref_genome_name_file"],
        min_reads= config["min_reads_variant"]
    log:
        "metadata/logs/vcfstats/{sample}.vcfstats"
    conda:
        "code/environments/Greference_tools.yml"
    shell:
        """
        freebayes \
            -C {params.min_reads} \
            -f {input.reference} \
            {input.bam} \
            > {output}
        
        rtg vcfstats {output} > {log}
        """


## 11 Variant Effect Prediction DB (version of GRCh38 109)
rule vep_install_db:
    output:
        touch("tasks/11vep_dependencies.done")
    params:
        species = config["vep_species"],
        assembly = config["vep_assembly"]
    conda: 
        "code/environments/vep.yml"
    shell:
        """
        vep_install -a cf \
            -s {params.species} \
            --ASSEMBLY {params.assembly} \
            --NO_UPDATE
        """


## 12 Running VEP in the command line 
rule vep_cli:
    input:
        script_dl_clivar = "code/04vep.sh",
        vcf = "results/variants/{sample}.vcf"
    output:
        "results/variants/vep/{sample}.txt"
    params:
        species = config["vep_species"],
        assembly=config["vep_assembly"]
    conda: 
        "code/environments/vep.yml"
    shell:
        """
        bash {input.script_dl_clivar}
        
        vep -i {input.vcf} \
                --offline \
                --force_overwrite \
                --assembly {params.assembly} \
                --appris \
                --biotype \
                --variant_class \
                --check_existing \
                --filter_common \
                --mane \
                --polyphen b \
                --pubmed \
                --regulatory \
                --sift b \
                --species {params.species} \
                --symbol \
                --transcript_version \
                --tsl \
                --cache \
                --tab \
                -o {output} \
                --custom data/ClinVar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN
        """


## 13 Parsing data from the VCF files with R
rule parsing_dataR:
    input:
        script = "code/05parsing_vep_data.R"
    output:
        touch("tasks/13parsing_dataR.done")
    params:
        sample_list=config["samples"],
        dir1 = "results/biostatistics/",
        dir2 = "results/biostatistics/tables",
        dir3 = "results/biostatistics/plots",
        dir4 = "results/biostatistics/joined_tables",
        gene_filter=config["gene_to_filterR"],
        chr_choosed=config["chromosome"],
    conda:
        "code/environments/biostatisticsR.yml"
    shell:
        """
        for dir in {params.dir1} {params.dir2} {params.dir3} {params.dir4}
        do
            if [[ ! -d $dir ]]
            then
                mkdir $dir
            fi
        done

        Rscript {input.script} {params.gene_filter} {params.sample_list} 

        ## Joining parsed tables for each sample
        cat results/biostatistics/tables/*{params.chr_choosed}* \
            | awk "!/^$(cat results/biostatistics/tables/*{params.chr_choosed}* \
            | cut -f 1 | head -n 1)/ || NR == 1" \
            > results/biostatistics/joined_tables/{params.gene_filter}.tsv
        """


##---------------------------------------------------------##
##     FINAL BOSS, YOU NEED TO HAVE THE 5 GENE TABLES      ##
##---------------------------------------------------------##

## 14 plotting with R
rule R_plotting:
    input:
        script = "code/06final_plot.R"
    output:
        png1="results/biostatistics/plots/final_plot.png",
        png2="results/biostatistics/plots/other_plot.png"
    conda:
        "code/environments/biostatisticsR.yml"
    shell:
        """
        Rscript {input.script} && echo "THE SCRIPT RAN WELL congrats :)"        

            # If this is red you didn't made it yet, maybe this will encourage you :)
            # https://www.youtube.com/watch?v=tYzMYcUty6s&ab_channel=TeamPsycosmos
        """
