#!/usr/bin/env Rscript


#################
##  LIBRARIES  ##
#################

suppressPackageStartupMessages(suppressWarnings({
  library(dplyr)
  library(glue)
  library(stringr)
  library(readr)
}))


########################
##  SAMPLES AND GENE  ##
########################

params <- commandArgs(trailingOnly = TRUE)

sample_list <- params[2:length(params)]
print(glue(">>> This is your sample list: {sample_list}"))


gene_filter <- params[1]
print(glue(">>> This is your gene: {gene_filter}"))

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