#!/usr/bin/env Rscript

library(tidyverse)

metrFile <- "${metricsFile}"

metrTab <- read_csv(metrFile)

pdf("${plotFile}", width = 11, height = 8.5)
metrTab %>%
    ggplot(aes(x = Barcode, y = `Number of Reads`)) +
    geom_col(fill = "#ec008c", colour = "#666666") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

metrTab %>%
    ggplot(aes(x = Barcode, y = `Estimated Number of Cells`)) +
    geom_col(fill = "#00b6ed", colour = "#666666") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()