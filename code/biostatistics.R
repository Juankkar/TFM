#!/usr/bin/env Rscript


## Packages 
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(library(glue))

## Selecting our sample
sample <- system(
    'read -p "Select your sample -> " sample ; echo $sample',
    intern=TRUE
)

## calculating the rows from the variants TXT to avoid them in order
## to read our data 
rows_skip <- system(
    glue("grep ^## results/variants/vep/{sample}.txt | wc -l"),
    intern = TRUE
    )

## Reading the data
variants <- read_table(glue("results/variants/vep/{sample}.txt"), 
                       skip=as.numeric(rows_skip))

## We will select the variables that we are interesting on
variants_filt <- variants %>% 
    select(Uploaded_variation="#Uploaded_variation", "Gene", 
           "Location","Allele","Feature","BIOTYPE",
           "CLIN_SIG", "PUBMED", "PHENO", "Feature_type",
           "Consequence", "PolyPhen")

print("===> Genes with variations <===")
genes <- variants %>%
  mutate(Gene = case_when(Gene == "-" ~ "UC",
                          !(Gene == "-") ~ as.character(Gene))) %>%
  group_by(Gene) %>%
  count() %>%
  arrange(desc(n))

genes

genes_table <- genes %>%
    mutate(sample=sample)
  
write.table(x = genes_table,
            file = glue("results/biostatistics/tables/{sample}_genes.txt"),
            sep="\t",
            row.names=F)

print("#########################################")
print("#########################################")

print("===> Biotype with variations <===")
biotype <- variants %>%
  mutate(BIOTYPE = case_when(BIOTYPE == "-" ~ "UC",
                             !(BIOTYPE == "-") ~ as.character(BIOTYPE))) %>%
  group_by(BIOTYPE) %>%
  count() %>%
  arrange(desc(n))

biotype

biotype_table <- biotype %>%
    mutate(sample=sample)
  
write.table(x = biotype_table,
            file = glue("results/biostatistics/tables/{sample}_biotype.txt"),
            sep="\t",
            row.names=F)

print("#########################################")
print("#########################################")

print("===> Clinical Significance <===")
clin_sig <- variants %>%
  mutate(CLIN_SIG = case_when(CLIN_SIG == "-" ~ "UC",
                          !(CLIN_SIG == "-") ~ as.character(CLIN_SIG))) %>%
  group_by(CLIN_SIG) %>%
  count() %>%
  arrange(desc(n))

clin_sig

clin_sig_table <- clin_sig %>% 
  mutate(sample=sample)

write.table(x = clin_sig_table,
            file = glue("results/biostatistics/tables/{sample}_clin_sig.txt"),
            sep="\t",
            row.names=F)

print("#########################################")
print("#########################################")

print("===> Consequence of the variations <===")
consequence <- variants %>%
  mutate(Consequence = case_when(Consequence == "-" ~ "UC",
                          !(Consequence == "-") ~ as.character(Consequence))) %>%
  group_by(Consequence) %>%
  count() %>%
  arrange(desc(n))

consequence

consequence_table <- biotype %>%
    mutate(sample=sample)
  
write.table(x = biotype_table,
            file = glue("results/biostatistics/tables/{sample}_consequence.txt"),
            sep="\t",
            row.names=F)


print("#########################################")
print("#########################################")

print("===> PUBMED ID Literature <===")
pubmed <- variants %>%
  mutate(PUBMED = case_when(PUBMED == "-" ~ "UC",
                          !(PUBMED == "-") ~ as.character(PUBMED))) %>%
  group_by(PUBMED) %>%
  count() %>%
  arrange(desc(n))

pubmed

pubmed_table <- pubmed %>%
    mutate(sample=sample)
  
write.table(x = pubmed_table,
            file = glue("results/biostatistics/tables/{sample}_pubmed.txt"),
            sep="\t",
            row.names=F)