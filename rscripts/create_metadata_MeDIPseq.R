
rm(list = ls())

library(magrittr)
library(dplyr)
library(readr)

indir = "data/BAM/hg19/bowtie_dip"
files = list.files(indir,full.names = TRUE,pattern = "sort.bam")
files = files[grep("bai",files,invert = TRUE)]

cell = ifelse(grepl("akata",basename(files)),"EBV","NOK")

treat = ifelse(grepl("mono",basename(files)),"NoTrt","CaFBS")

repl = basename(files) %>% strsplit("rep") %>% unlist
repl = repl[grepl("sort",repl)]
repl = paste0("Rep", gsub(".sort.bam","",repl) )

type = ifelse(grepl("Input",basename(files)),"Input",
              ifelse(grepl("hmC",basename(files)),"hmC","mC"))

df  = tibble( file = files,
              cell , treatment = treat, type,  rep = repl)

write_delim(df , path = "data/metadata/MeDIPseq_definition.tsv" ,
            delim = "\t")
