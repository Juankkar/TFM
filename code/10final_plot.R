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

## Processing to get the values that we want

df_apc <- apc %>%
    select("sample","symbol","variant_class") %>% 
    separate_wider_delim(sample, delim="_", names=c("sample","chr")) %>% 
    select(-"chr") 

df_braf <- braf %>%
    select("sample","symbol","variant_class") %>% 
    separate_wider_delim(sample, delim="_", names=c("sample","chr")) %>% 
    select(-"chr") 

df_kras <- kras %>%
    select("sample","symbol","variant_class") %>% 
    separate_wider_delim(sample, delim="_", names=c("sample","chr")) %>% 
    select(-"chr") 

df_pik3ca <- pik3ca %>%
    select("sample","symbol","variant_class") %>% 
    separate_wider_delim(sample, delim="_", names=c("sample","chr")) %>% 
    select(-"chr") 

df_tp53 <- tp53 %>%
    select("sample","symbol","variant_class") %>% 
    separate_wider_delim(sample, delim="_", names=c("sample","chr")) %>% 
    select(-"chr")

## Joining the data together
suppressMessages(suppressWarnings({
    df_genes <- rbind(df_apc,df_braf,df_kras,df_pik3ca,df_tp53) %>%
        group_by(sample,symbol,variant_class) %>%
        summarise(n=n())
    
    df_genes_processed <- expand.grid(sample=as.character(unique(as.character(df_genes$sample))),
                                      symbol=as.character(unique(as.character(df_genes$symbol))),
                                      variant_class=as.character(unique(as.character(df_genes$variant_class)))) %>% 
    left_join(
        df_genes, by=c("sample","symbol","variant_class")
        ) %>% 
    mutate(n=ifelse(is.na(n),0,n))
}))

## Plotting the data

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
                                      "Single nucleotide polimorfism (SNP)"))
    ) %>%
    ggplot(aes(sample, n, fill=symbol)) +
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    scale_fill_manual(values=c("white","lightgray",
                               "#99DBF5","#FFEEBB",
                               "#EA906C")) +
    scale_y_continuous(expand=expansion(0),
                       limits=c(0,250),
                       breaks=seq(0,250,50)) +
    facet_wrap(~variant_class, ncol=1, scales="free_y") +
    labs(
        title = "Distribution of the variants per genes and samples",
        x = "Samples",
        y = "Number of variants",
        fill="Genes"
    ) +
    theme_test() +
    theme(
        plot.title=element_text(hjust=.5,size=16,face="bold"),
        plot.background=element_rect(linewidth=1, color="black"),
        axis.title=element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=10),
        axis.text.x=element_text(angle = 20, hjust = .9, vjust = .9),
        axis.ticks.x = element_line(linewidth=.5),
        axis.ticks.length.x = unit(.2,"cm"),
        legend.position="top",
        legend.title=element_text(size=12, face="bold"),
        strip.background=element_blank(),
        strip.text=element_text(size=11)
    )
    
ggsave(file="results/biostatistics/plots/final_plot.png",
       plot=final_plot,
       heigh=7,
       width=9)

print("##==========================##")
print("##===> WORK FINISHED!!! <===##")
print("##==========================##")
