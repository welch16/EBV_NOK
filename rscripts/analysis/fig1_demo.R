
rm(list = ls())

library(patchwork)
library(scales)
library(tidyverse)
library(DESeq2)

## From Mark's description

## Compare NOKS vs NOKS+EBV (undifferentiated) with Scott's data
##
## Volcano plot
## Pathway analysis with GSEA
## Possibly other pathway analysis

datadr = "data/RSEM/hg19"
figsdr = "figs/fig2"

files = datadr %>%
    list.files(full.names = TRUE) %>%
    str_subset("genes.results") %>%
    str_subset("Nok")

rsem = files %>%
    set_names(
        basename(.) %>%
        str_replace(".genes.results","")) %>%
    map(read_tsv) %>%
    tibble(
        file = names(.),
        rsem = .)  %>%
    mutate(
        cell = if_else(file %>% str_detect("EBV"),"NOKS+EBV","NOKS"),
        treat = if_else(file %>% str_detect("mono"),"None","MC"))

## quick check of count histograms
library(ggjoy)

theme_set(theme_bw())

log2_count_histogram = rsem %>%
    unnest() %>%
    ggplot(aes(expected_count))+
    geom_histogram(aes(fill = interaction(cell,treat)),
                   bins = 45,colour = "black")+
    scale_fill_brewer(palette = "Pastel1",name = "")+
    facet_grid( file ~ .)+
    theme(
        legend.position = "top",
        strip.text.y = element_text(angle = 0))+
    scale_x_log10()+
    geom_vline(xintercept = 100,colour = "red",
               linetype = 2)

ggsave(file = file.path(figsdr,"Histograms_log2_expected_counts_scott.pdf"),
       plot = log2_count_histogram,
       height = unit(8,"in"),
       width = unit(5,"in"))

## RLE analysis

medians = rsem %>%
    unnest() %>%
    group_by(gene_id) %>%
    summarize(
        median = median(expected_count)) %>%
    ungroup()

rsem = rsem %>%
    unnest() %>%
    inner_join(medians,by = "gene_id") %>%
    mutate(
        rle = log2( expected_count / median)) %>%
    mutate(median = NULL)

rle_histogram = rsem %>% 
    ggplot(aes(rle))+
    geom_histogram(aes(fill = interaction(cell,treat)),
                   bins = 45,colour = "black")+
    scale_fill_brewer(palette = "Pastel1",name = "")+
    facet_grid( file ~ .)+
    theme(
        legend.position = "top",
        strip.text.y = element_text(angle = 0))+
    geom_vline(xintercept = 0,colour = "red",
               linetype = 2)+
    scale_x_continuous(limits = c(-5,5))

ggsave(file = file.path(figsdr,"Histogram_RLE_scott.pdf"),
       plot = rle_histogram,
       height = unit(8,"in"),
       width = unit(5,"in"))
       
