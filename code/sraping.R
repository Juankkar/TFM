library(tidyverse)
library(glue)
library(rentrez)


entrez_dbs()


entrez_db_summary(db="pubmed")
entrez_db_searchable(db="sra")


busqueda <- entrez_search(db ="pubmed", term = "human genome")
busqueda %>% glimpse()

x <- entrez_fetch(db="pubmed", id=busqueda, rettype = "abstract")

anio <- 1950:2022
mus_musculus <- glue("mice genome AND {anio}[PDAT]")
human_genome <- glue("human genome AND {anio}[PDAT]")

human_genome_search <- tibble(
  anio=anio,
  human_genome = human_genome,
  mus_musculus = mus_musculus
) %>%
  mutate(
    mice_pubmed = map_dbl(mus_musculus, 
                          ~entrez_search(db="pubmed", term=.x)$count),
    human_pubmed = map_dbl(human_genome, 
                           ~entrez_search(db="pubmed", term=.x)$count),
    # mice_genome = map_dbl(mus_musculus, 
    #                       ~entrez_search(db="sra", term=.x)$count),
    # human_genome2 = map_dbl(human_genome, 
    #                         ~entrez_search(db="sra", term=.x)$count)
    )


human_genome_search %>% 
  select(-c("mus_musculus", "human_genome")) %>%
  pivot_longer(-anio, 
               values_to = "valores", 
               names_to = "names") %>% 
  ggplot(aes(anio, valores, fill=names, color=names)) +
  geom_bar(stat="identity",size=1, position = position_dodge()) +
  scale_fill_manual(values = c("pink", "skyblue","green", "yellow")) +
  scale_color_manual(values = c("red", "blue","forestgreen", "orange")) +
  theme_test() +
  theme(legend.position = "top")
