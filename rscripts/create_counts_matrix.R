#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--metadatafile",action = "store_true",type = "character",
                help = "Input tsv file consisting of the metadata used to generate the PCA plot:
                          file | cell | treatment | lab | rep "),
    make_option("--outfile",action = "store_true",type = "character",default = tempfile(),
                help = "File where the count matrix is saved")
    )

opt = parse_args(OptionParser(option_list = optList))

opt$metadatafile = "data/metadata/PCA_definition.tsv"

stopifnot(file.exists(opt$metadatafile))

library(base,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(tximport,quietly = TRUE)
library(readr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(tidyr,quietly = TRUE)

metadata = read_tsv(opt$metadatafile)

txdata = tximport(metadata$file,
                  type = "rsem" , importer = read_tsv)

countdata = txdata[["counts"]]
genes = rownames(countdata)

files = metadata$file %>% basename
files = gsub(".genes.results","",files)

out = countdata %>% as.data.frame %>% as.tbl %>%
    mutate(gene = genes) %>%
    select(gene,everything())
names(out) = c("gene",files)

write_tsv(out,opt$outfile)
