#!/usr/bin/env Rscript


#################
##  LIBRARIES  ##
#################

suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(library(glue))

########################
##  SAMPLES AND GENE  ##
########################

sample_list <- scan(
  file="your_sample_list.txt",
  what = character(),
  quiet = TRUE,
  sep = "\n"
)

print("this is your sample list: ")
sample_list

param <- commandArgs(trailingOnly = TRUE)

gene_filter <- param[1]

print("This is your gene: ")
gene_filter

####################
## EXECUTION LOOP ##
####################

for(sample in sample_list){

sample_file <- glue("results/variants/vep/{sample}.txt")

#####################
##  READING DATA   ## 
##       AND       ##
##   PROCESSESING  ##
#####################

rows_skip <- tibble(
  our_vector = scan(
    file=sample_file,
    what = character(),
    quiet = TRUE,
    sep = "\n")) %>%
  filter(str_detect(our_vector,"^##")) %>%
  nrow()

variants_selec <- read_table(sample_file, 
                       skip=as.numeric(rows_skip)) %>%
  filter(SYMBOL == glue("{gene_filter}")) %>%
  mutate(sample=glue("{sample}"))

colnames(variants_selec) <- tolower(colnames(variants_selec))  

########################
## SAVING PARSED DATA ##
########################

write_tsv(x = variants_selec,
          file = glue("results/biostatistics/tables/{sample}.tsv"))
variants_selec %>% print(n=6)

}

print("##==========================##")
print("##===> WORK FINISHED!!! <===##")
print("##==========================##")