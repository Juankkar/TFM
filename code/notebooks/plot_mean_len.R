#!usr/bin/env Rscript

suppressMessages(suppressWarnings({
        library(tidyverse)
}))

read_csv("len_mean.csv") %>%
    ggplot(aes(sample, mean_len)) +
    geom_bar(stat="identity", width=.5, color="black", fill="lightgray") +
    geom_errorbar(aes(ymin= mean_len - sd_len,
                      ymax= mean_len + sd_len), width=.3) +
    scale_y_continuous(expand=expansion(0),
                       limits=c(0, 135),
                       breaks=seq(0,130, 20)) +
    labs(
        title = "Avarage length of the sequences for each sample",
        x = "Sample",
        y = "Sequence length (pb)"
    ) +
    theme_classic() +
    theme(
        plot.title=element_text(hjust=.5, face="bold", size = 14),
        axis.title=element_text(face="bold", size = 12),
        axis.text= element_text(color="black",size=11),
        axis.text.x= element_text(angle=35, vjust=1, hjust=1)
    )

ggsave("../../results/objective/mean_len.png",
       width=6,
       height=5)
