#!/usr/bin/env Rscript


#################
##  LIBRARIES  ##
#################

suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(library(glue))

################
##  SAMPLES   ##
################

sample_list <- scan(
  file="your_sample_list.txt",
  what = character(),
  quiet = TRUE,
  sep = "\n"
)

sample_list

param <- commandArgs(trailingOnly = TRUE)

gene_filter <- param[1]

gene_filter

for(sample in sample_list){

#####################
##  READING DATA   ## 
##       AND       ##
##  PRE-PROCESSED  ##
#####################

sample_file <- glue("results/variants/vep/{sample}.txt")

rows_skip <- tibble(
  our_vector = scan(
    file=sample_file,
    what = character(),
    quiet = TRUE,
    sep = "\n")) %>%
  filter(str_detect(our_vector,"^##")) %>%
  nrow()

variants <- read_table(sample_file, 
                       skip=as.numeric(rows_skip)) %>%
  filter(SYMBOL == gene_filter)

#################
##  EXECUTION  ##
#################

variant_class <- variants %>%
  filter(VARIANT_CLASS != "-") %>%
  group_by(VARIANT_CLASS) %>%
  count() %>%
  arrange(desc(n))

colnames(variant_class) <- tolower(colnames(variant_class))

variant_class

variant_class_table <- variant_class %>%
    mutate(sample=sample)

col1_variant_class_name <- colnames(variant_class_table)[1]
  
write_tsv(x = variant_class_table,
          file = glue("results/biostatistics/tables/{sample}_{col1_variant_class_name}.tsv"))

print(glue("===> VARIATION CLASSES {sample} DONE <==="))

print("#########################################")
print("#########################################")

biotype <- variants %>%
  filter(BIOTYPE != "-") %>%
  group_by(BIOTYPE) %>%
  count() %>%
  arrange(desc(n))

colnames(biotype) <- tolower(colnames(biotype))

biotype

biotype_table <- biotype %>%
    mutate(sample=sample)

col1_biotype_name <- colnames(biotype_table)[1]

write_tsv(x = biotype_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_biotype_name}.tsv"))

print(glue("===> BIOTYPE OF THE VARIATIONS {sample} DONE <==="))

print("#########################################")
print("#########################################")

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

print(glue("===> CLINICAL SIGNIFICANCE {sample} DONE <==="))

print("#########################################")
print("#########################################")

consequence <- variants %>%
  filter(Consequence != "-") %>%
  separate_longer_delim(Consequence, delim=",") %>%
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

print(glue("===> CONSEQUENCE OF THE VARIATIONS {sample} DONE <==="))

print("#########################################")
print("#########################################")

sift <- variants %>%
  filter(SIFT != "-") %>%
  mutate(SIFT = str_remove_all(SIFT, pattern = "[:punct:]"),
         SIFT = str_remove_all(SIFT, pattern = "\\d")) %>%
  group_by(SIFT) %>%
  count() %>%
  arrange(desc(n))

colnames(sift) <- tolower(colnames(sift))

sift

sift_table <- sift %>%
    mutate(sample=sample)

col1_sift_name <- colnames(sift)[1]
  
write_tsv(x = sift_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_sift_name}.tsv"))

print(glue("===> SIFT OF THE VARIATION {sample} DONE <==="))

print("#########################################")
print("#########################################")

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

print(glue("===> POLYPHEN STATUS {sample} DONE <==="))

print("#########################################")
print("#########################################")

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

print(glue("===> CLINVAR STATUS {sample} DONE <==="))

print("#########################################")
print("#########################################")

clinvar_clndn <- variants %>%
  select("ClinVar", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN") %>%
  filter(ClinVar_CLNDN != "-") %>%
  group_by(ClinVar_CLNDN) %>%
  count()

colnames(clinvar_clndn) <- tolower(colnames(clinvar_clndn))

clinvar_clndn_table <- clinvar_clndn %>%
    mutate(sample=sample)

col1_clinvar_clndn_name <- colnames(clinvar_clndn_table)[1]
  
write_tsv(x = clinvar_clndn_table,
          file = glue("results/biostatistics/tables/{sample}_{col1_clinvar_clndn_name}.tsv"))

print(glue("===> CLINVAR DIAGNOSIS {sample} DONE <==="))

print("#########################################")
print("#########################################")

num_variants <- tibble(num_variants=nrow(variants),
                       sample=sample)

num_variants

col1_rows_name <- colnames(num_variants)[1]

write_tsv(x = num_variants,
          file = glue("results/biostatistics/tables/{sample}_{col1_rows_name}.tsv"))

print(glue("===> NUMBER OF VARIANTS {sample} DONE <==="))

print("#########################################")
print("#########################################")

location <- variants %>%
  select(Location) %>%
  mutate(Location = str_replace(Location, 
                                pattern = "^\\d+:",
                                replacement="")) %>%
  separate_longer_delim(Location, "-") %>%
  mutate(Location=as.numeric(Location)) %>%
  group_by(Location) %>%
  summarise(n=n()) %>%
  arrange(Location) %>%
  drop_na()

colnames(location) <- tolower(colnames(location))

location 

location_table <- location %>%
    mutate(sample=sample) 

col1_location_name <- colnames(location_table)[1]
  
write_tsv(x = location_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_location_name}.tsv"))

variants %>% select("Location")

print(glue("===> DISTIBUTION CHROMOSOME {sample} DONE <==="))

print("#########################################")
print("#########################################")

protein_pos <- variants %>%
  mutate(Protein_position = as.numeric(Protein_position)) %>%
  group_by(Protein_position) %>%
  summarise(n=n()) %>%
  arrange(Protein_position) %>%
  drop_na()

colnames(protein_pos) <- tolower(colnames(protein_pos))

protein_pos 

protein_pos_table <- protein_pos %>%
    mutate(sample=sample) 

col1_protein_pos_name <- colnames(protein_pos_table)[1]
  
write_tsv(x = protein_pos_table,
            file = glue("results/biostatistics/tables/{sample}_{col1_protein_pos_name}.tsv"))

print(glue("===> PROTEIN POSITION {sample} DONE <==="))

}

print("##==========================##")
print("##===> WORK FINISHED!!! <===##")
print("##==========================##")