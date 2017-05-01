rm(list = ls())

library(tidyverse)

dr = "data/summary/hg19/RNAseq"

files = list.files(dr,full.names = TRUE)

ff = files %>% map(read_tsv) %>%
    map(.f = function(x) x %>% filter(seqnames != "chrM") %>%
                         mutate(prop = nreads / sum(nreads))) %>%
    map2(.f = function(x,y) x %>% mutate(file = y),basename(files)) %>%
    bind_rows




ff %>% filter(seqnames == "chr17") %>%
    ggplot(
        aes(file,prop)) + geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))+
    coord_flip()
