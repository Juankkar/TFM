#!/usr/bin/env Rscript

library(tidyverse)

vcf_stats <- read_tsv("../results_analysis/tables/ratios_vcfstats.tsv")  %>% 
    select(-"...1", -"r_snp_transition_transversion", 
           -"r_insert_del", 
           -"r_indel_snp_plus_mnp") %>%
    pivot_longer(-sample, names_to="het_hom", values_to="num") %>% 
    filter(het_hom != "r_total_het_hom")   

print("===> Total values")
vcf_stats %>% nrow()
print("===> Total values pre group")
vcf_stats %>% 
    group_by(het_hom) %>% 
    summarise(n=n())

vcf_stats %>% 
   ggplot(aes(het_hom, num, fill=het_hom)) +
   geom_boxplot()

ggsave("../results_analysis/plots/statistics_vcf_inference.png",
       width=7,
       height=5)

print("=------------------------------------=")
print("===>           Normality          <===")
print("=------------------------------------=")

shapiro.test(vcf_stats$num)

vcf_stats %>% 
    group_by(het_hom) %>% 
    rstatix::shapiro_test(num) 

vcf_stats %>%
    rstatix::kruskal_test(num ~ het_hom) 

vcf_stats %>%
    rstatix::dunn_test(num ~ het_hom, p.adj = "bonf")
