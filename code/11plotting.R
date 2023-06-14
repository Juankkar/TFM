#!/usr/bin/env Rscript

#################
##  LIBRARIES  ##
#################

suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(gt)))
suppressPackageStartupMessages(suppressWarnings(library(glue)))

###############
##  SYSTEM   ##
###############

nrows_location_plot <- system('read -p "For the location plot you have to choose the number of rows: " rows ; echo $rows',
                              intern=TRUE)

## The table metadata/severe_consequences.csv was got from ENSEMBLE using --> code/annotations.ipynb 

high_consequence <- 'echo $(cat metadata/severe_consequences.csv | grep "HIGH" | cut -d "," -f 1) | tr " " "|"'
consequence_level_high <- c(system(high_consequence, intern=TRUE))

moderate_consequence <- 'echo $(cat metadata/severe_consequences.csv | grep "MODERATE" | cut -d "," -f 1) | tr " " "|"'
consequence_level_moderate <- c(system(moderate_consequence,intern=TRUE))

low_consequence <- 'echo $(cat metadata/severe_consequences.csv | grep "LOW" | cut -d "," -f 1) | tr " " "|"'
consequence_level_low <- c(system(low_consequence,intern=TRUE))

modifier_consequece <- 'echo $(cat metadata/severe_consequences.csv | grep "MODIFIER" | cut -d "," -f 1) | tr " " "|"'
consequence_level_modifier <- c(system(modifier_consequece,intern=TRUE))

## Algorithm to add a new character value column based on a list given separate with "|" data from a list 
## throwback to my TFG :)
selec <-function(ord,lista_tokens,var) {
  paste(lista_tokens[-ord],collapse="|")
  if(!is.na(ord)) return(grepl(lista_tokens[ord],tolower(var)) & !grepl(paste(lista_tokens[-ord],collapse="|"),tolower(var)))
  else return(grepl(paste(lista_tokens,collapse="|"),tolower(var)))
}

list_level_consequence <- tolower(c(consequence_level_high,consequence_level_moderate,consequence_level_low,consequence_level_modifier))

list_coding_consequence <- c("splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant|stop_lost|start_lost",
                             "inframe_insertion|inframe_deletion|missense_variant|protein_altering_variant|coding_sequence_variant|incomplete_terminal_codon_variant",
                             "start_retained_variant|stop_retained_variant|synonymous_variant")
#####################
##  READING DATA   ## 
##       AND       ##
##  PRE-PROCESSED  ##
#####################

biotype <- read_tsv("results/biostatistics/joined_tables/biotype.tsv")

## This code it is necessary in order to add a "0" when a value of the first column (biotype,
## consequence...) does not appear in one of the samples
df_biotype <- expand.grid(biotype=as.character(unique(as.character(biotype$biotype))),
                          sample=as.character(unique(as.character(biotype$sample)))) %>% 
    left_join(
        biotype, by=c("biotype","sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n)) 

clin_sig <- read_tsv("results/biostatistics/joined_tables/clin_sig.tsv") 
df_clin_sig <- expand.grid(clin_sig=as.character(unique(as.character(clin_sig$clin_sig))),
                           sample=as.character(unique(as.character(clin_sig$sample)))) %>% 
    left_join(
        clin_sig, by=c("clin_sig","sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n)) %>%
    separate_longer_delim(clin_sig, delim=",") %>%
    group_by(clin_sig, sample) %>%
    summarise(n=sum(n)) %>%
    ungroup()  

clinvar <- read_tsv("results/biostatistics/joined_tables/clinvar_clnsig.tsv") 
df_clinvar <- expand.grid(clinvar_clnsig=as.character(unique(as.character(clinvar$clinvar_clnsig))),
                          sample=as.character(unique(as.character(clinvar$sample)))) %>% 
    left_join(
        clinvar, by=c("clinvar_clnsig","sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n)) 

consequence <- read_tsv("results/biostatistics/joined_tables/consequence.tsv") 
df_consequence <- expand.grid(consequence=as.character(unique(as.character(consequence$consequence))),
                              sample=as.character(unique(as.character(consequence$sample)))) %>% 
    left_join(
        consequence, by=c("consequence","sample")
        ) %>% 
    mutate(
        n=ifelse(is.na(n),0,n),
        consequence_lower=tolower(consequence),
        level= case_when(selec(1,list_level_consequence,consequence_lower)~"high" ,
                         selec(2,list_level_consequence,consequence_lower)~"moderate",
                         selec(3,list_level_consequence,consequence_lower)~"low",
                         selec(4,list_level_consequence,consequence_lower)~"modifier"),
        coding= as.character(case_when(selec(1,list_coding_consequence,consequence_lower)~"coding",
                             selec(2,list_coding_consequence,consequence_lower)~"coding",
                             selec(3,list_coding_consequence,consequence_lower)~"coding")),
        coding=ifelse(is.na(coding),"not_coding",coding)
    ) %>% 
    select(-"consequence_lower") 

location <- read_tsv("results/biostatistics/joined_tables/location.tsv") 
df_location <- expand.grid(location=as.character(unique(as.character(location$location))),
                           sample=as.character(unique(as.character(location$sample)))) %>% 
    left_join(
        location %>% mutate(location=as.character(location)), by=c("location","sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n),
           location=as.numeric(location)) 

num_variants <- read_tsv("results/biostatistics/joined_tables/num_variants.tsv")

# prot_position <- read_tsv("results/biostatistics/joined_tables/protein_position.tsv") %>%
#     group_by(sample) %>%
#     mutate(percentage=(protein_position/max(protein_position))*100,
#            percentiles=case_when(percentage >= 0 & percentage <= 10 ~ "0-10%",
#                                  percentage > 10 & percentage <= 20 ~ "10-20%",
#                                  percentage > 20 & percentage <= 30 ~ "20-30%",
#                                  percentage > 30 & percentage <= 40 ~ "30-40%",
#                                  percentage > 40 & percentage <= 50 ~ "40-50%",
#                                  percentage > 50 & percentage <= 60 ~ "50-60%",
#                                  percentage > 60 & percentage <= 70 ~ "60-70%",
#                                  percentage > 70 & percentage <= 80 ~ "70-80%",
#                                  percentage > 80 & percentage <= 90 ~ "80-90%",
#                                  percentage > 90 & percentage <= 100 ~ "90-100%")) %>%
#     select(-protein_position) %>%
#     ungroup() %>%
#     group_by(sample, percentiles) %>%
#     summarise(n=sum(n))

# prot_position %>% print(n=Inf)


polyphen <- read_tsv("results/biostatistics/joined_tables/polyphen.tsv") 
df_polyphen <- expand.grid(polyphen=as.character(unique(as.character(polyphen$polyphen))),
                           sample=as.character(unique(as.character(polyphen$sample)))) %>% 
    left_join(
        polyphen, by=c("polyphen","sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n)) 

sift <- read_tsv("results/biostatistics/joined_tables/sift.tsv") 
df_sift <- expand.grid(sift=as.character(unique(as.character(sift$sift))),
                           sample=as.character(unique(as.character(sift$sample)))) %>% 
    left_join(
        sift %>% mutate(sift=as.character(sift)), by=c("sift","sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n)) 

variant_class <- read_tsv("results/biostatistics/joined_tables/variant_class.tsv") 
df_variant_class <- expand.grid(variant_class=as.character(unique(as.character(variant_class$variant_class))),
                       sample=as.character(unique(as.character(variant_class$sample)))) %>% 
    left_join(
        variant_class, by=c("variant_class","sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n)) 

#################
##  EXECUTION  ##
#################

#----------------------#
# Studying the biotype #
#----------------------#
biotype_wider <- df_biotype %>%
    pivot_wider(biotype, names_from=sample,values_from=n)

biotype_wider_fixed <- biotype_wider %>%
    mutate(row_sum = rowSums(biotype_wider[,-1])
        #    biotype = case_when(row_sum < 100 ~ "Others",
        #                        row_sum >= 100 ~ as.character(biotype))
                               ) %>% 
    select(-"row_sum") %>%
    pivot_longer(-biotype, names_to="sample", values_to="n") %>%
    group_by(biotype,sample) %>%
    summarize(n = sum(n)) %>%
    arrange(desc(n)) %>%
    pivot_wider(biotype, names_from="sample", values_from="n") %>%
    ungroup()

ggbiotype <- biotype_wider_fixed %>% 
    pivot_longer(-biotype, names_to="sample", values_to="n") 

heatmap_biotyoe <- ggbiotype %>%
    ggplot(aes(reorder(sample,n,decreasing=TRUE), 
               reorder(biotype,n, decreasing=FALSE), 
               fill = n)) +
    geom_tile(color="black") +
    geom_text(aes(label=n)) +
    scale_fill_gradient(name="BIOTYPE\nvs\nSAMPLE",
                        low="white", high="red") +
    labs(
        title="Biotype of the variants for each sample",
        x=NULL,
        y=NULL
    ) +
    theme(
        panel.background=element_blank(),
        axis.line.x=element_line(),
        plot.title=element_text(size=14, face="bold", hjust=.5),
        axis.title.x = element_text(margin = margin(t = 10),size = 12),
        axis.title.y = element_text(margin = margin(r = 10), size = 12),
        axis.text.x = element_text(angle = 90,vjust = .05),
        axis.text = element_text(size=10, color="black")
    )

# max_biotype <- max(ggbiotype$n) 

# biotype_plot <- ggbiotype %>% 
#     ggplot(aes(n, reorder(sample,n), fill=biotype)) +
#     geom_bar(stat="identity", position="dodge") +
#     scale_x_continuous(expand=expansion(0),
#                        limits=c(0,max_biotype+100)) +
#     labs(
#         title = "Biotype of the variants for each sample",
#         x = "Number of variants",
#         y = "Samples",
#         fill="BIOTYPE"
#     ) +
#     theme_classic() +
#     theme(
#         plot.title=element_text(hjust=.5, size=14, face="bold"),
#         axis.title=element_text(size=12, face="bold"),
#         axis.text=element_text(size=11, color="black")
#     )

## Saving Processed data
gt_biotype <- biotype_wider_fixed %>%
    gt() %>%
    tab_header(title=md("Biotype product from the variation"))

# ggsave(filename = "results/biostatistics/plots/biotype.png",
#        plot = biotype_plot,
#        height = 5,
#        width = 10)

ggsave(filename = "results/biostatistics/plots/biotype_heatmap.png",
       plot = heatmap_biotyoe,
       height = 5,
       width = 10)

gtsave(data=gt_biotype,
       filename=glue("results/biostatistics/plots/biotype_plot.html"))

#---------------------------------------------#
#  Studying the Clinical Significance status  #
#---------------------------------------------#
clin_sig_wider <- df_clin_sig %>%
    pivot_wider(clin_sig, names_from=sample,values_from=n) 

clinsig_wider_fixed <- clin_sig_wider %>%
    pivot_longer(-clin_sig, names_to="sample", values_to="n") %>%
    group_by(clin_sig,sample) %>%
    summarize(n = sum(n)) %>%
    arrange(desc(n)) %>%
    pivot_wider(clin_sig, names_from="sample", values_from="n") %>%
    ungroup()

ggclinsig <- clinsig_wider_fixed %>% 
    pivot_longer(-clin_sig, names_to="sample", values_to="n") 

max_clinsig <- max(ggclinsig$n) 

clin_sig_plot <- ggclinsig %>% 
    ggplot(aes(n, reorder(sample,n), fill=clin_sig)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_clinsig+100)) +
    labs(
        title = "Clinical Significance of the variants",
        x = "Number of variants",
        y = "Samples",
        fill="Clinical Significance"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

## Saving Processed data
gt_clin_sig <- clinsig_wider_fixed %>%
    gt() %>%
    tab_header(title=md("Clinical Significance for the variants"))

gtsave(data=gt_clin_sig,
       filename=glue("results/biostatistics/plots/clin_sig_plot.html"))

ggsave(filename = "results/biostatistics/plots/clin_sig.png",
       plot = clin_sig_plot,
       height = 5,
       width = 10)

#-------------------------------#
#  Studying the ClinVar status  #
#-------------------------------#

max_clinvar_clnsig <- max(df_clinvar$n) 

clinvar_plot <- df_clinvar %>%
    ggplot(aes(n, reorder(sample,n), fill=clinvar_clnsig)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_clinvar_clnsig+100)) +
    labs(
        title = "ClinVar Status of the variants for each sample",
        x = "Number of variants",
        y = "Samples"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

gt_clinvar <- df_clinvar %>%
    pivot_wider(clinvar_clnsig, names_from="sample", values_from="n") %>%
    gt() %>%
    tab_header(title=md("CliVar for the variants"))

gtsave(data=gt_clinvar,
       filename=glue("results/biostatistics/plots/clinvar.html"))

ggsave(filename = "results/biostatistics/plots/clinvar.png",
       plot = clinvar_plot,
       height = 5,
       width = 10)

#----------------------------------------------#
# Studying the consequences for the variation  #
#----------------------------------------------#

consequence_wider <- df_consequence %>%
    select("consequence", "sample", "n") %>%
    pivot_wider(consequence, names_from=sample,values_from=n)

consequence_level <- df_consequence %>% 
    group_by(sample, level) %>%
    summarise(n=sum(n))

max_consequence_level <- max(consequence_level$n) 

consequence_level_plot <- consequence_level %>%
    mutate(level=factor(level,
        levels=c("high", "low", "moderate", "modifier")
    )) %>%
    ggplot(aes(n, reorder(sample,n), fill=level)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_consequence_level+100)) +
    labs(
        title = "Consequences of the variation (severe levels)",
        x = "Number of variants",
        y = "Samples",
        fill="SEVERE"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

consequence_coding <- df_consequence %>% 
    group_by(sample, coding) %>%
    summarise(n=sum(n))

max_consequence_coding <- max(consequence_coding$n) 

consequence_coding_plot <- consequence_coding %>%
    ggplot(aes(n, reorder(sample,n), fill=coding)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_consequence_coding+100)) +
    labs(
        title = "Consequences of the variation (Coding/not Coding)",
        x = "Number of variants",
        y = "Samples",
        fill="VARIANT"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

consequence_coding_only <- df_consequence %>% 
    filter(coding == "coding")

max_con_coding_only <- max(consequence_coding_only$n) 

con_coding_only_plot <- consequence_coding_only %>%
    # mutate(consequence=case_when(n < 50 ~ "other",
    #                              n > 50 ~ as.character(consequence))) %>%
    ggplot(aes(n, reorder(sample,n), fill=consequence)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_con_coding_only+100)) +
    labs(
        title = "Consequences of the variation (coding only)",
        x = "Number of variants",
        y = "Samples",
        fill="CODING VARIANT"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

gt_consequence <- consequence_wider %>%
    gt() %>%
    tab_header(title=md("Consequence product from the variation"))

## Saving Processed data
ggsave(filename = "results/biostatistics/plots/consequence_level.png",
       plot = consequence_level_plot,
       height = 5,
       width = 10)

ggsave(filename = "results/biostatistics/plots/consequence_coding.png",
       plot = consequence_coding_plot,
       height = 5,
       width = 10)

ggsave(filename = "results/biostatistics/plots/consequence_coding_only.png",
       plot = con_coding_only_plot,
       height = 5,
       width = 10)

gtsave(data=gt_consequence,
           filename=glue("results/biostatistics/plots/consequence.html"))

#----------------------------------------#
# Studying the location of the variants  #
#----------------------------------------#

location_plot <- df_location %>%
    mutate(location=location/10^6) %>%
    ggplot(aes(location, n, color=sample)) +
    geom_line(show.legend=FALSE) +
    facet_wrap(~sample, nrow=as.numeric(glue("{nrows_location_plot}"))) +
    labs(
         title = "Distribution of variants on the chromosome",
         x="Position (mb)",
         y="Number of variants"
    ) +
    theme_test() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

ggsave(filename = "results/biostatistics/plots/location.png",
       plot = location_plot,
       height = 5,
       width = 10)

#--------------------------------------#
# Studying the number of the variants  #
#--------------------------------------#

max_number <- max(num_variants$num_variants) 

nvariants_plot <- num_variants %>%
    ggplot(aes(num_variants, reorder(sample,num_variants))) +
    geom_bar(stat="identity", fill="#1e81b0", color="black") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_number+1000)) +
    labs(
        title = "Number of variants of echa samples",
        x = "Number of variants",
        y = "Samples"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

gt_number <- num_variants %>%
    gt() %>%
    tab_header(title=md("Number of variants in each sample"))

## Saving Processed data
ggsave(filename = "results/biostatistics/plots/num_variants.png",
       plot = nvariants_plot,
       height = 5,
       width = 10)

gtsave(data=gt_number,
           filename=glue("results/biostatistics/plots/num_variants.html"))

#--------------------------------#
#  Studying the PolyPhen status  #
#--------------------------------#

max_polyphen <- max(df_polyphen$n) 

polyphen_plot <- df_polyphen %>%
    ggplot(aes(n, reorder(sample,n), fill=polyphen)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_polyphen+100)) +
    labs(
        title = "PolyPhen Status of the variants for each sample",
        x = "Number of variants",
        y = "Samples"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

gt_polyphen <- df_polyphen %>%
    pivot_wider(polyphen, names_from="sample", values_from="n") %>%
    gt() %>%
    tab_header(title=md("PolyPhen for the variants"))


gtsave(data=gt_polyphen,
       filename=glue("results/biostatistics/plots/polyphen.html"))

ggsave(filename = "results/biostatistics/plots/polyphen.png",
       plot = polyphen_plot,
       height = 5,
       width = 10)

#----------------------------#
#  Studying the Sift status  #
#----------------------------#

max_sift <- max(df_sift$n) 

sift_plot <- df_sift %>%
    ggplot(aes(n, reorder(sample,n), fill=sift)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_sift+100)) +
    labs(
        title = "Sift Status of the variants",
        x = "Number of variants",
        y = "Samples",
        fill="SIFT"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

gt_sift <- df_sift %>%
    pivot_wider(sift, names_from="sample", values_from="n") %>%
    gt() %>%
    tab_header(title=md("Sift for the variants"))

gtsave(data=gt_sift,
       filename=glue("results/biostatistics/plots/sift.html"))

ggsave(filename = "results/biostatistics/plots/sift.png",
       plot = sift_plot,
       height = 5,
       width = 10)

#--------------------------------#
#  Studying the Variant Classes  #
#--------------------------------#

max_variant_class <- max(df_variant_class$n) 

variant_class_plot <- df_variant_class %>%
    ggplot(aes(n, reorder(sample,n), fill=variant_class)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_variant_class+100)) +
    labs(
        title = "Variants Classes",
        x = "Number of variants",
        y = "Samples",
        fill="CLASSES"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

gt_variant_class <- df_variant_class %>%
    pivot_wider(variant_class, names_from="sample", values_from="n") %>%
    gt() %>%
    tab_header(title=md("variant_class for the variants"))

gtsave(data=gt_variant_class,
       filename=glue("results/biostatistics/plots/variant_class.html"))

ggsave(filename = "results/biostatistics/plots/variant_class.png",
       plot = variant_class_plot,
       height = 5,
       width = 10)