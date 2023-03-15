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
dna <- glue("dna sequencing AND {anio}[PDAT]")
rna <- glue("RNA-seq AND {anio}[PDAT]")

human_genome_search <- tibble(
  anio=anio,
  human_genome = human_genome,
  mus_musculus = mus_musculus
) %>%
  mutate(
    dna_seq = map_dbl(dna, 
                      ~entrez_search(db="pubmed", term=.x)$count),
    rna_seq = map_dbl(rna, 
                      ~entrez_search(db="pubmed", term=.x)$count)
    )


human_genome_search %>% 
  select(-c("mus_musculus", "human_genome")) %>%
  pivot_longer(-anio, 
               values_to = "valores", 
               names_to = "names") %>% 
  mutate(names=factor(names,
                      levels = c("dna_seq", "rna_seq"),
                      labels = c("DNA-seq", "RNA-seq"))) %>% 
  ggplot(aes(anio, valores, fill=names, color=names)) +
  geom_bar(stat="identity",size=1, position = position_dodge()) +
  scale_fill_manual(values = c("pink", "skyblue")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(0,30000),
                     breaks = seq(0,30000,5000)) +
  scale_x_continuous(expand = expansion(0),
                     limits = c(1950, 2022),
                     breaks = seq(1949,2022, 10)) + 
  labs(
    title = "Publicaciones DNA-seq vs RNA-seq, PUBMED",
    x="Año",
    y="Número de publicaciones",
    fill = NULL,
    color=NULL
  ) +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = .5, size = 18))
