#!/usr/bin/env Rscript

suppressMessages(suppressWarnings({
        library(tidyverse)
}))

vcf_stats <- read_tsv("../results_analysis/tables/vcf_basic_stats.tsv") 
ratios_vcf <- read_tsv("../results_analysis/tables/ratios_vcfstats.tsv")

basic_stats <- vcf_stats %>%
    select(-"pass_fail_filter") %>%
    pivot_longer(-sample, 
                 names_to="names", 
                 values_to="values")  %>% 
    mutate(names = factor(names,
                          levels=c("snps",
                                   "mnps",
                                   "insertions",
                                   "deletions",
                                   "indels",
                                   "same_reference"),
                          labels=c("SNPs",
                                   "MNPs",
                                   "Insertions",
                                   "Deletions",
                                   "Indels",
                                   "Same as\nreference")),
            values=values/10^3) %>% 
    ggplot(aes(sample, values, fill = names)) +
    geom_bar(stat="identity", 
             position="dodge", 
             color="black", 
             size=.25) +
    scale_y_continuous(expand=expansion(0),
                      limits=c(0,32),
                      breaks=seq(0,30,5)) +
    scale_fill_manual(values=c("white","black",
                               "lightgray","#99DBF5",
                               "#FFEEBB","#EA906C")) +
    labs(
        title = "Variants from all the samples (VCF)",
        x = "Samples",
        y = "Nº of variants x 10³",
        fill= NULL
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

ggsave("../results_analysis/plots/vcf_basic_stats.png",
       plot=basic_stats,
       width=7,
       height=5)


ratios1 <- ratios_vcf  %>% 
    select(-"...1", 
           "r_snp_transition_transversion", 
           "r_insert_del", 
           "r_indel_snp_plus_mnp")  

plot_ratios1 <- ratios1 %>%
    pivot_longer(-sample) %>% 
    filter(name %in% c("r_snp_transition_transversion", 
                       "r_insert_del", 
                       "r_indel_snp_plus_mnp")) %>% 
    mutate(name=factor(name,
                       levels=c("r_snp_transition_transversion", 
                                "r_insert_del", 
                                "r_indel_snp_plus_mnp"),
                       labels=c("SNP Transition/Transversion", 
                                "Insertion/Deletion", 
                                "Indel/SNP+MNP")))  %>% 
    ggplot(aes(sample, value, fill=name)) +
    geom_bar(stat="identity", 
             position="dodge", 
             color="black") +
    scale_y_continuous(expand=expansion(0),
                      limits=c(0,3),
                      breaks=seq(0,3,.5)) +
    scale_fill_manual(values=c("white", 
                               "lightgray", 
                               "#99DBF5")) +
    labs(
        title = "Ratios from the VCF",
        x = "Samples",
        y = "Ratio",
        fill= NULL
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

ggsave("../results_analysis/plots/ratios_vcf.png",
       plot=plot_ratios1,
       width=7,
       height=5) 


ratios2 <- ratios_vcf  %>% 
    select(-"...1", 
           -"r_snp_transition_transversion", 
           -"r_insert_del", 
           -"r_indel_snp_plus_mnp")

plot_ratios2_total <- ratios2 %>%
    pivot_longer(-sample) %>% 
    filter(name == "r_total_het_hom") %>%
    ggplot(aes(sample, value)) +
    geom_bar(stat="identity", 
             color="black", 
             fill="black", 
             width=.5) +
    scale_y_continuous(expand=expansion(0),
                      limits=c(0,2),
                      breaks=seq(0,2,.5)) +
    labs(
        title = "Total ratio from the Heterozygotes/Homozygotes VCF",
        x = "Samples",
        y = "Ratio"
    ) +
    theme_classic() +
    theme(
        plot.background=element_rect(fill="white", color="white"),
        panel.background=element_rect(fill="white", color="white"),
        plot.title=element_text(hjust=.5, face="bold", size = 14),
        axis.title=element_text(face="bold", size = 12),
        axis.text= element_text(color="black",size=11),
        axis.text.x= element_text(angle=35, vjust=1, hjust=1),
        legend.position="top",
        strip.placement="outside",
        strip.text=element_text(face="bold", size=12),
        axis.ticks.x=element_blank()
    )


ggsave("../results_analysis/plots/ratios_het_hom_total_vcf.png",
       plot=plot_ratios2_total,
       width=6,
       height=5) 

plot_ratios2 <- ratios2 %>%
    pivot_longer(-sample) %>% 
    filter(name != "r_total_het_hom") %>%
    mutate(name=factor(name,
                       levels=c("r_snp_het_hom", 
                                "r_mnp_het_hom", 
                                "r_Insert_het_hom",
                                "r_del_het_hom", 
                                "r_Indel_het_hom"),
                       labels=c("SNP", 
                                "MNP", 
                                "Insert", 
                                "Deletion", 
                                "Indel"))) %>%
    ggplot(aes(sample, value, fill=name)) +
    geom_bar(stat="identity", 
             position="dodge", 
             color="black") +
    scale_y_continuous(expand=expansion(0),
                      limits=c(0,4),
                      breaks=seq(0,4,.5)) +
    scale_fill_manual(values=c("white",
                               "lightgray",
                               "#99DBF5",
                               "#FFEEBB",
                               "#EA906C")) +
    labs(
        title = "Ratios from the VCF, Heterozygotes/Homozygotes",
        x = "Samples",
        y = "Ratio",
        fill= NULL
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


ggsave("../results_analysis/plots/ratios_het_hom_vcf.png",
       plot=plot_ratios2,
       width=7,
       height=6) 
