#!/usr/bin/env Rscript
suppressMessages(suppressWarnings({
        library(tidyverse)
}))

suppressMessages(suppressWarnings({
    vcf_stats <- read_tsv("../results_analysis/tables/ratios_vcfstats.tsv")  %>% 
        select(-"...1", -"r_snp_transition_transversion", 
               -"r_insert_del", 
               -"r_indel_snp_plus_mnp") %>%
        pivot_longer(-sample, names_to="het_hom", values_to="num") %>% 
        filter(het_hom != "r_total_het_hom")   
}))    

print("===> Total values")
vcf_stats %>% nrow()
print("===> Total values pre group")
vcf_stats %>% 
    group_by(het_hom) %>% 
    summarise(n=n())

print("=------------------------------------=")
print("===>           Normality          <===")
print("=------------------------------------=")

shapiro.test(vcf_stats$num)

vcf_stats %>% 
    group_by(het_hom) %>% 
    rstatix::shapiro_test(num) 

print("=------------------------------------=")
print("===>      Kruskal_Wallis/Dunnet   <===")
print("=------------------------------------=")

vcf_stats %>%
    rstatix::kruskal_test(num ~ het_hom) 

vcf_stats %>%
    rstatix::dunn_test(num ~ het_hom, p.adj = "bonf")  %>%  
    write_tsv("../results_analysis/tables/invference_vcf.tsv") %>%
    print(n=Inf)


vcf_stats %>% 
    mutate(sig=case_when(het_hom == "r_mnp_het_hom" ~ "Sig",
                         het_hom != "r_mnp_het_hom" ~ "no.sig"),
           het_hom = factor(het_hom,
                            levels = c("r_Indel_het_hom","r_del_het_hom",
                                       "r_Insert_het_hom", "r_mnp_het_hom",
                                       "r_snp_het_hom"),
                            labels = c("Indel", "Deletion", "Insertion",
                                       "MNP","SNP"))) %>%
    ggplot(aes(num, het_hom, fill=sig)) +
    geom_boxplot(show.legend=FALSE, width=.65) +
    scale_x_continuous(expand = expansion(0),
                      limits = c(0,4)) +
    scale_fill_manual(values = c("#73d206", "#d00808")) +
    labs(
        title = "Differences between ratios",
        subtitle = "X² Kruskal-Wallis (*p* < 0.05) | Post-Hoc: Dunnet Test (Bonferroni)",
        x = "Ratio",
        y = "Variation type"
    ) +
    theme_classic() +
     theme(
        plot.background=element_rect(fill="white", color="white"),
        panel.background=element_rect(fill="white", color="white"),
        plot.title=element_text(hjust=.5, face="bold", size = 14),
        plot.subtitle=ggtext::element_markdown(hjust=.5, face="bold", size = 11.5),
        axis.title=element_text(face="bold", size = 12),
        axis.text= element_text(color="black",size=11),
    )

ggsave("../results_analysis/plots/statistics_vcf_inference.png",
       width=6,
       height=4)