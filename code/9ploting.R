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

#####################
##  READING DATA   ## 
##       AND       ##
##  PRE-PROCESSED  ##
#####################

biotype <- read_tsv("results/biostatistics/joined_tables/biotype.tsv")

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
    mutate(n=ifelse(is.na(n),0,n)) 

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
    mutate(n=ifelse(is.na(n),0,n)) 

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
    mutate(row_sum = rowSums(biotype_wider[,-1]),
           biotype = case_when(row_sum < 100 ~ "Others",
                               row_sum >= 100 ~ as.character(biotype))) %>% 
    select(-"row_sum") %>%
    pivot_longer(-biotype, names_to="sample", values_to="n") %>%
    group_by(biotype,sample) %>%
    summarize(n = sum(n)) %>%
    arrange(desc(n)) %>%
    pivot_wider(biotype, names_from="sample", values_from="n") %>%
    ungroup()

ggbiotype <- biotype_wider_fixed %>% 
    pivot_longer(-biotype, names_to="sample", values_to="n") 

max_biotype <- max(ggbiotype$n) 

biotype_plot <- ggbiotype %>% 
    ggplot(aes(n, reorder(sample,n), fill=biotype)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_biotype+100)) +
    labs(
        title = "Biotype of the variants for each sample",
        x = "Number of variants",
        y = "Samples",
        fill="BIOTYPE"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

## Saving Processed data
gt_biotype <- biotype_wider_fixed %>%
    gt() %>%
    tab_header(title=md("Biotype product from the variation"))

ggsave(filename = "results/biostatistics/plots/biotype.png",
       plot = biotype_plot,
       height = 5,
       width = 10)

for(type_file in c("html", "docx")){
    gtsave(data=gt_biotype,
           filename=glue("results/biostatistics/plots/biotype_plot.{type_file}"))
}

#---------------------------------------------#
#  Studying the Clinical Significance status  #
#---------------------------------------------#
clin_sig_wider <- df_clin_sig %>%
    pivot_wider(clin_sig, names_from=sample,values_from=n) 

clinsig_wider_fixed <- clin_sig_wider %>%
    mutate(row_sum = rowSums(clin_sig_wider[,-1]),
           clin_sig = case_when(row_sum < 20 ~ "other",
                                row_sum >= 20 ~ as.character(clin_sig))) %>%
    select(-"row_sum") %>%
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

for(type_file in c("html", "docx")){
    gtsave(data=gt_clin_sig,
           filename=glue("results/biostatistics/plots/clin_sig_plot.{type_file}"))
}

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

for(type_file in c("html", "docx")){
    gtsave(data=gt_clinvar,
           filename=glue("results/biostatistics/plots/clinvar.{type_file}"))
}

ggsave(filename = "results/biostatistics/plots/clinvar.png",
       plot = clinvar_plot,
       height = 5,
       width = 10)

#----------------------------------------------#
# Studying the consequences for the variation  #
#----------------------------------------------#

consequence_wider <- df_consequence %>%
    pivot_wider(consequence, names_from=sample,values_from=n)

consequence_wider_fixed <- consequence_wider %>%
    mutate(row_sum = rowSums(consequence_wider[,-1]),
           consequence = case_when(row_sum < 200 ~ "other",
                                   row_sum >= 200 ~ as.character(consequence))) %>%
    select(-"row_sum") %>%
    pivot_longer(-consequence, names_to="sample", values_to="n") %>%
    group_by(consequence,sample) %>%
    summarize(n = sum(n)) %>%
    arrange(desc(n)) %>%
    pivot_wider(consequence, names_from="sample", values_from="n") %>%
    ungroup()

ggconsequence <- consequence_wider_fixed %>% 
    pivot_longer(-consequence, names_to="sample", values_to="n") 

max_consequence <- max(ggconsequence$n) 

consequence_plot <- ggconsequence %>% 
    ggplot(aes(n, reorder(sample,n), fill=consequence)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(expand=expansion(0),
                       limits=c(0,max_consequence+100)) +
    labs(
        title = "Consequence of the variants",
        x = "Number of variants",
        y = "Samples",
        fill="CONSEQUENCE"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, size=14, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        axis.text=element_text(size=11, color="black")
    )

gt_consequence <- consequence_wider_fixed %>%
    gt() %>%
    tab_header(title=md("Consequence product from the variation"))

## Saving Processed data
ggsave(filename = "results/biostatistics/plots/consequence.png",
       plot = consequence_plot,
       height = 5,
       width = 10)

for(type_file in c("html", "docx")){
    gtsave(data=gt_consequence,
           filename=glue("results/biostatistics/plots/consequence.{type_file}"))
}

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

for(type_file in c("html", "docx")){
    gtsave(data=gt_number,
           filename=glue("results/biostatistics/plots/num_variants.{type_file}"))
}

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

for(type_file in c("html", "docx")){
    gtsave(data=gt_polyphen,
           filename=glue("results/biostatistics/plots/polyphen.{type_file}"))
}

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

for(type_file in c("html", "docx")){
    gtsave(data=gt_sift,
           filename=glue("results/biostatistics/plots/sift.{type_file}"))
}

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

for(type_file in c("html", "docx")){
    gtsave(data=gt_variant_class,
           filename=glue("results/biostatistics/plots/variant_class.{type_file}"))
}

ggsave(filename = "results/biostatistics/plots/variant_class.png",
       plot = variant_class_plot,
       height = 5,
       width = 10)