#!/usr/bin/env Rscript

suppressMessages(suppressWarnings({
        library(tidyverse)
}))

## Reading the data
suppressMessages(suppressWarnings({
    apc <- read_tsv("results/biostatistics/joined_tables/APC.tsv")
    braf <- read_tsv("results/biostatistics/joined_tables/BRAF.tsv")
    kras <- read_tsv("results/biostatistics/joined_tables/KRAS.tsv")
    pik3ca <- read_tsv("results/biostatistics/joined_tables/PIK3CA.tsv")
    tp53 <- read_tsv("results/biostatistics/joined_tables/TP53.tsv")
}))

## joining tables and processing
df_genes <- rbind(apc,braf,kras,pik3ca,tp53) %>% 
    separate_wider_delim(sample, delim="_", names=c("sample","chr")) %>% 
    select(-"chr")

suppressMessages(suppressWarnings({
df_nvariants <- df_genes %>%
    select("sample","symbol","variant_class") %>%
    group_by(sample,symbol,variant_class) %>%
        summarise(n=n())

df_clinvar <- df_genes %>%
    select("clinvar_clnsig", "clinvar_clndn", "sample") %>%
    separate_longer_delim(clinvar_clndn, delim="|") %>%
    separate_longer_delim(clinvar_clnsig, delim="/") %>%
    group_by(clinvar_clnsig, clinvar_clndn,sample) %>%
        summarise(n=n()) %>%
    mutate(clinvar_clnsig = str_replace_all(clinvar_clnsig, pattern="_", replacement=" "),
           clinvar_clndn = str_replace_all(clinvar_clndn, pattern="_", replacement=" ")) %>%
    filter(clinvar_clnsig != "-" & 
           clinvar_clndn %in% c("APC-Associated Polyposis Disorders", "Carcinoma of colon",
                                "Colorectal cancer", "Familial multiple polyposis syndrome",
                                "Familial colorectal cancer"))
}))

## Varaintes totales
df_genes_processed <- expand.grid(sample=as.character(unique(as.character(df_nvariants$sample))),
                                  symbol=as.character(unique(as.character(df_nvariants$symbol))),
                                  variant_class=as.character(unique(as.character(df_nvariants$variant_class)))) %>% 
left_join(
    df_nvariants, by=c("sample","symbol","variant_class")
    ) %>% 
    mutate(n=ifelse(is.na(n),0,n))
    
    ## Variantes de ClinVar
df_clinvar <- expand.grid(clinvar_clnsig=as.character(unique(as.character(df_clinvar$clinvar_clnsig))),
                          clinvar_clndn=as.character(unique(as.character(df_clinvar$clinvar_clndn))),
                          sample=as.character(unique(as.character(df_clinvar$sample)))) %>% 
    left_join(
        df_clinvar, by=c("clinvar_clnsig", "clinvar_clndn", "sample")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n))

## Plotting the data

variant_class_pretty <- c("indel" = "Copy number variation (CNV)",
                          "SNV" = "Single nucleotide variation (SNV)")

final_plot <- df_genes_processed %>%
    mutate(
        sample=factor(sample,
                      levels=sort(unique(df_genes$sample))),
        symbol=factor(symbol,
                      levels=c("PIK3CA","APC","BRAF","KRAS","TP53")),
        variant_class=factor(variant_class,
                             levels=c("indel","SNV"),
                             labels=c("Copy number variation (CNV)",
                                      "Single nucleotide variant (SNV)"))
    ) %>%
    ggplot(aes(sample, n, fill=symbol)) +
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    scale_fill_manual(values=c("white","lightgray",
                               "#99DBF5","#FFEEBB",
                               "#EA906C")) +
    scale_y_continuous(expand=expansion(0),
                       limits=c(0,260),
                       breaks=seq(0,250,50)) +
    facet_wrap(~variant_class, ncol=1, 
               strip.position="left", scales="free_y") +
    labs(
        title = "Distribution of the variants per genes and samples",
        x = "Samples",
        y = "Number of variants",
        fill="Genes"
    ) +
    theme_test() +
    theme(
        plot.title=element_text(hjust=.5,size=16,face="bold"),
        axis.title=element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=10),
        axis.text.x=element_text(angle = 20, hjust = .9, vjust = .9),
        axis.ticks.x = element_line(linewidth=.5),
        axis.ticks.length.x = unit(.2,"cm"),
        legend.position="top",
        legend.title=element_text(size=12, face="bold"),
        strip.background=element_blank(),
        strip.text=element_text(size=10.5, face="bold"),
        strip.placement="outside"
    )
    
ggsave(file="results/biostatistics/plots/final_plot.png",
       plot=final_plot,
       heigh=7,
       width=9)


o <- df_clinvar %>% 
    mutate(clinvar_clndn = factor(clinvar_clndn,
                                  levels = c("APC-Associated Polyposis Disorders", "Carcinoma of colon",
                                             "Colorectal cancer", "Familial multiple polyposis syndrome",
                                             "Familial colorectal cancer"),
                                  labels = c("Associated\nPolyposis\nDisorders\n(APC)", "Carcinoma\ncolon",
                                             "Colorectal\ncancer", "Familial\nmultiple\npolyposis\nsyndrome",
                                             "Familial\ncolorectal\ncancer"))) %>%
    ggplot(aes(clinvar_clndn, n, fill = clinvar_clnsig)) +
    geom_bar(stat="identity", position="dodge") + 
    facet_wrap(~sample, ncol=2)


ggsave(file="results/biostatistics/plots/other_plot.png",
       plot=o,
       heigh=7,
       width=9)

print("##==========================##")
print("##===> WORK FINISHED!!! <===##")
print("##==========================##")
