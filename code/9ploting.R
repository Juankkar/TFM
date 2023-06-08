#!/usr/bin/env Rscript

###############
##  SYSTEM   ##
###############

nrows_location_plot <- system('read -p "For the location plot you have to choose the number of rows: " rows ; echo $rows',
                              intern=TRUE)

#################
##  LIBRARIES  ##
#################

suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(gt)))
suppressPackageStartupMessages(suppressWarnings(library(glue)))

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


polyphen <- read_tsv("results/biostatistics/joined_tables/polyphen.tsv") 
df_polyphen <- expand.grid(polyphen=as.character(unique(as.character(polyphen$polyphen))),
                           sample=as.character(unique(as.character(polyphen$sample)))) %>% 
    left_join(
        polyphen, by=c("polyphen","sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n)) 

#################
##  EXECUTION  ##
#################

## Studying the biotype
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


#biotype_wider %>% mutate(sumatorio = rowSums(biotype_wider[,-1])) %>% select(sumatorio) %>% arrange(desc(sumatorio)) %>% 

for(type_file in c("html", "docx")){
    gtsave(data=gt_biotype,
           filename=glue("results/biostatistics/plots/biotype_plot.{type_file}"))
}

## Studying the Clinical Significance status

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

# Studying the ClinVar status

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

# Studying the consequences for the variation
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

## Studying the location of the variants
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