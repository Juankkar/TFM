#!/usr/bin/env python

import pygraphviz as pgv

# Definir las reglas del grafo
graph_content = """
digraph snakemake_rules {{
    rankdir=LR;
    size="15,6";
    label="Rule 1 to Rule 15"
    
    // Node styling
    node [style=filled];
    "rule download_data" [fillcolor=orange];
    "rule pre_processing" [fillcolor=green];
    "rule reference_genome" [fillcolor=green];
    "rule fastqc" [fillcolor=green];
    "rule fastp" [fillcolor=green];
    "rule fastqc_trimmed" [fillcolor=green];
    "rule bwa_mapping" [fillcolor=green];
    "rule sam_to_bam" [fillcolor=green];
    "rule delete_duplicates" [fillcolor=green];
    "rule extracting_variants" [fillcolor=green];
    "rule vep_install_db" [fillcolor=orange];
    "rule vep_cli" [fillcolor=green];
    "rule parsing_dataR" [fillcolor=green];
    "rule R_plotting" [fillcolor=orange];
    
    // Edge connections
    
    "1. Preparación del flujo" -> "rule download_data" -> "rule pre_processing"  -> "rule reference_genome"; 
    "2. Calidad de las secuencias" -> "rule fastqc" -> "rule fastp" -> "rule fastqc_trimmed" -> "3. Mapeado, SAM/BAM" -> "rule bwa_mapping" -> "rule sam_to_bam" -> "rule delete_duplicates" -> "4. Estudio de variantes" -> "rule extracting_variants" -> "rule vep_install_db" -> "rule vep_cli" -> "rule parsing_dataR" -> "rule R_plotting";
}}
"""

# Crear un objeto de gráfico
graph = pgv.AGraph(string=graph_content)

# Establecer las opciones de diseño
graph.layout(prog='dot')

# Guardar el gráfico como una imagen
graph.draw('../results/snakemake_rules_linear.png')
