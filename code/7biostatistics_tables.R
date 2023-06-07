#!/usr/bin/env Rscript


## Packages 
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(library(glue))

## Selecting our sample
sample <- system(
    'read -p "Select your sample -> " sample ; echo $sample',
    intern=TRUE
)

## calculating the rows from the variants tsv to avoid them in order
## to read our data 
rows_skip <- system(
    glue("grep ^## results/variants/vep/{sample}.txt | wc -l"),
    intern = TRUE
    )

## Reading the data
variants <- read_table(glue("results/variants/vep/{sample}.txt"), 
                       skip=as.numeric(rows_skip))

print("===> Genes with variations <===")
genes <- variants %>%
  filter(Gene != "-") %>%
  group_by(Gene) %>%
  count() %>%
  arrange(desc(n))

colnames(genes) <- tolower(colnames(genes))

genes

genes_table <- genes %>%
    mutate(sample=sample)

col1_genes_name <- colnames(genes_table)[1]
  
write_tsv(x = genes_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_genes_name}.tsv"))

print("#########################################")
print("#########################################")

print("===> Biotype with variations <===")
biotype <- variants %>%
  filter(BIOTYPE != "-") %>%
  group_by(BIOTYPE) %>%
  count() %>%
  arrange(desc(n))

colnames(biotype) <- tolower(colnames(biotype))

biotype

biotype_table <- biotype %>%
    mutate(sample=sample)

col1_biotype_name <- colnames(genes_table)[1]

write_tsv(x = biotype_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_biotype_name}.tsv"))

print("#########################################")
print("#########################################")

print("===> Clinical Significance <===")
clin_sig <- variants %>%
  filter(CLIN_SIG != "-") %>%
  group_by(CLIN_SIG) %>%
  count() %>%
  arrange(desc(n))

colnames(clin_sig) <- tolower(colnames(clin_sig))

clin_sig

clin_sig_table <- clin_sig %>% 
  mutate(sample=sample)

col1_clin_sig_name <- colnames(clin_sig_table)[1]

write_tsv(x = clin_sig_table,
            file = glue("results/biostatistics/tables/{sample}_clin_sig.tsv"))

print("#########################################")
print("#########################################")

print("===> Consequence of the variations <===")
consequence <- variants %>%
  filter(Consequence != "-") %>%
  group_by(Consequence) %>%
  count() %>%
  arrange(desc(n))

colnames(consequence) <- tolower(colnames(consequence))

consequence

consequence_table <- consequence %>%
    mutate(sample=sample)
  
col1_consequence_name <- colnames(consequence_table)[1]
  
write_tsv(x = consequence_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_consequence_name}.tsv"))


print("#########################################")
print("#########################################")

print("===> PUBMED ID Literature <===")
pubmed <- variants %>%
  filter(PUBMED != "-") %>%
  group_by(PUBMED) %>%
  count() %>%
  arrange(desc(n))

colnames(pubmed) <- tolower(colnames(pubmed))

pubmed

pubmed_table <- pubmed %>%
    mutate(sample=sample)

col1_pubmed_name <- colnames(pubmed_table)[1]
  
write_tsv(x = pubmed_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_pubmed_name}.tsv"))

print("#########################################")
print("#########################################")

print("===> PolyPhen status <===")
polyphen <- variants %>%
  filter(PolyPhen != "-") %>%
  mutate(PolyPhen = str_remove_all(PolyPhen, pattern = "[:punct:]"),
         PolyPhen = str_remove_all(PolyPhen, pattern = "\\d")) %>%
  group_by(PolyPhen) %>%
  count()

colnames(polyphen) <- tolower(colnames(polyphen))

polyphen

polyphen_table <- polyphen %>%
    mutate(sample=sample)

col1_polyphen_name <- colnames(polyphen_table)[1]
  
write_tsv(x = polyphen_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_polyphen_name}.tsv"))

print("#########################################")
print("#########################################")

print("===> ClinVar status <===")
clinvar <- variants %>%
  select("ClinVar", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN") %>%
  filter(ClinVar_CLNSIG != "-") %>%
  group_by(ClinVar_CLNSIG) %>%
  count()

colnames(clinvar) <- tolower(colnames(clinvar))

clinvar

clinvar_table <- clinvar %>%
    mutate(sample=sample)

col1_clinvar_name <- colnames(clinvar_table)[1]
  
write_tsv(x = clinvar_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_clinvar_name}.tsv"))