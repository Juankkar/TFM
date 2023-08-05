#!/usr/bin/env Rscript

library(tidyverse)

flagstats <- read_tsv("../results_analysis/tables/flagstats.tsv")

max_num_reads <- max(flagstats$num_seqs_afterQC )
y_axis_labels <- c("FALSE" = "Probabilidad de precipitación", 
                   "TRUE" = "Promedio de precipitación\npor evento (mm)")
flagstats %>%
    select(num_seqs_afterQC, mapped_reads_per, properly_paired_per, singletons_per, sample) %>% 
    pivot_longer(-sample, names_to="names", values_to="values") %>% 
    mutate(facet = case_when(names == "num_seqs_afterQC" ~ "Sequences after QC",
                             names != "num_seqs_afterQC" ~ "Percentage"),
           names = factor(names,
                          levels=c("num_seqs_afterQC", "mapped_reads_per",
                                   "properly_paired_per", "singletons_per")),
           facet = factor(facet,
                          levels=c("Sequences after QC", "Percentage"))) %>% 
    ggplot(aes(sample, values, fill = names)) +
    geom_bar(stat="identity", color = "black",
             width = .5, position="dodge") +
    facet_wrap(~facet, scale="free", ncol=1, 
               strip.position="left") +
    scale_y_continuous(expand=expansion(.05)) +
    scale_fill_manual(values=c("lightgray","white", "#99DBF5","#FFEEBB", "#EA906C"),
                      labels = c("Seqs after QC", "% Mapped reads", 
                                 "% Properly paired", "% Singletons")) +
    labs(
        title = "Statistics from flagstats",
        x = "Samples",
        y = NULL,
        fill=NULL
    ) + 
    theme_minimal() +
    theme(
        plot.background=element_rect(fill="white", color="white"),
        panel.background=element_rect(fill="white", color="white"),
        plot.title=element_text(hjust=.5, face="bold", size = 14),
        axis.title=element_text(face="bold", size = 12),
        axis.text= element_text(color="black",size=11),
        axis.text.x= element_text(angle=35, vjust=1, hjust=1),
        legend.position="top",
        strip.placement="outside",
        strip.text=element_text(face="bold", size=12)
    )
    
ggsave("../results_analysis/plots/flagstats.png",
       width=6.5, height=6)
