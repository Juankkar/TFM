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

print("===> Genes with variations <===")
genes <- variants %>%
  filter(Gene != "-") %>%
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
  filter(BIOTYPE != "-") %>%
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
  filter(CLIN_SIG != "-") %>%
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
  filter(Consequence != "-") %>%
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
  filter(PUBMED != "-") %>%
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

print("#########################################")
print("#########################################")

print("===> PolyPhen status Literature <===")
polyphen <- variants %>%
  filter(PolyPhen != "-") %>%
  mutate(PolyPhen = str_remove_all(PolyPhen, pattern = "[:punct:]"),
         PolyPhen = str_remove_all(PolyPhen, pattern = "\\d")) %>%
  group_by(PolyPhen) %>%
  count()

polyphen

polyphen_table <- polyphen %>%
    mutate(sample=sample)
  
write.table(x = polyphen_table,
            file = glue("results/biostatistics/tables/{sample}_polyphen.txt"),
            sep="\t",
            row.names=F)

print("#########################################")
print("#########################################")

print("===> ClinVar status Literature <===")
clinvar <- variants %>%
  select("ClinVar", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN") %>%
  filter(ClinVar_CLNSIG != "-") %>%
  group_by(ClinVar_CLNSIG) %>%
  count()

clinvar

clinvar_table <- clinvar %>%
    mutate(sample=sample)
  
write.table(x = clinvar_table,
            file = glue("results/biostatistics/tables/{sample}_clinvar.txt"),
            sep="\t",
            row.names=F)