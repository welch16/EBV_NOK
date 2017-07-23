
rm(list =ls())

library(tidyverse)

rawdr = "./data/RSEM/hg19"
files = rawdr %>%
    list.files(pattern = "result",full.names = TRUE)

files = files %>%
    split(
        files %>%
        basename() %>%
        strsplit("\\.") %>%
        map_chr( ~ .[1]))


files = files[ grep("clone1",names(files))][[1]]

example_data = files %>%
    map(read_tsv)

names(example_data) = c("gene",
                        "iso")

ex1 = "ENSG00000000003_TSPAN6"

my_example = example_data %>%
    map(filter,gene_id == ex1)

## $gene
##                  gene_id
## 1 ENSG00000000003_TSPAN6
##                                                                                              transcript_id(s)
## 1 ENST00000373020_TSPAN6-001,ENST00000431386_TSPAN6-201,ENST00000494424_TSPAN6-002,ENST00000496771_TSPAN6-003
##   length effective_length expected_count   TPM FPKM
## 1 1434.1           1384.1         367.34 19.53 15.6

## $iso
##                transcript_id                gene_id length effective_length
## 1 ENST00000373020_TSPAN6-001 ENSG00000000003_TSPAN6   2206             2156
## 2 ENST00000431386_TSPAN6-201 ENSG00000000003_TSPAN6   1099             1049
## 3 ENST00000494424_TSPAN6-002 ENSG00000000003_TSPAN6    820              770
## 4 ENST00000496771_TSPAN6-003 ENSG00000000003_TSPAN6   1025              975
##   expected_count  TPM FPKM IsoPct
## 1         210.84 7.20 5.75  36.85
## 2          72.81 5.11 4.08  26.15
## 3          45.29 4.33 3.46  22.16
## 4          38.39 2.90 2.31  14.84

## RSEM evaluates at the isoform levels and then computes the results for a gene_id,
## for each file we processed at the gene level, we have a similar file at the
## isoform level.
