#!usr/bin/env Rscript

suppressMessages(suppressWarnings({
        library(tidyverse)
}))

read_csv("len_mean.csv") %>%
    ggplot(aes(sample, mean_len)) +
    geom_bar(stat="identity", width=.75, color="black", fill="lightgray") +
    geom_errorbar(aes(ymin= mean_len - sd_len,
                      ymax= mean_len + sd_len), width=.3) +
    scale_y_continuous(expand=expansion(0),
                       limits=c(0, 110),
                       breaks=seq(0,105, 15)) +

    theme_classic()

ggsave("mean_len.png")
