#!/usr/bin/env Rscript

suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(gt)))
suppressPackageStartupMessages(suppressWarnings(library(glue)))

## Tables
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

polyphen <- read_tsv("results/biostatistics/joined_tables/polyphen.tsv") 
df_polyphen <- expand.grid(polyphen=as.character(unique(as.character(polyphen$polyphen))),
                           sample=as.character(unique(as.character(polyphen$sample)))) %>% 
    left_join(
        polyphen, by=c("polyphen","sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n)) 




## Studying the biotype
biotype_wider <- df_biotype %>%
    pivot_wider(biotype, names_from=sample,values_from=n) 

gt_biotype <- biotype_wider %>%
    gt() %>%
    tab_header(title=md("Biotype product from the variation"))

for(type_file in c("html", "docx")){
    gtsave(data=gt_biotype,
           filename=glue("results/biostatistics/plots/biotype_plot.{type_file}"))
}

## Studying the Clinical Significance status

clin_sig_wider <- df_clin_sig %>%
    pivot_wider(clin_sig, names_from=sample,values_from=n) 

gt_clin_sig <- clin_sig_wider %>%
    gt() %>%
    tab_header(title=md("Clinical Significance for the variants"))

for(type_file in c("html", "docx")){
    gtsave(data=gt_clin_sig,
           filename=glue("results/biostatistics/plots/clin_sig_plot.{type_file}"))
}

## Studying the ClinVar status

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


ggsave(filename = "results/biostatistics/plots/clinvar.png",
       plot = clinvar_plot,
       height = 5,
       width = 10)

## Studying the PolyPhen status

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

ggsave(filename = "results/biostatistics/plots/polyphen.png",
       plot = polyphen_plot,
       height = 5,
       width = 10)