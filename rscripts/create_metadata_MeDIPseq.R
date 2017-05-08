
rm(list = ls())

library(magrittr)
library(dplyr)
library(readr)

indir = "data/BAM/hg19/bowtie_dip"
files = list.files(indir,full.names = TRUE,pattern = "sort.bam")
files = files[grep("bai",files,invert = TRUE)]

filesInput_pool = files[grep("rep0",files)]
files = files[grep("rep0",files,invert = TRUE)]

filesInput_pool = c(filesInput_pool,files[grep("Input",files,invert = TRUE)])

create_structure <- function(files)
{
    cell = ifelse(grepl("akata",basename(files)),"EBV","NOK")
    treat = ifelse(grepl("mono",basename(files)),"NoTrt","CaFBS")

    repl = basename(files) %>% strsplit("rep") %>% unlist
    repl = repl[grepl("sort",repl)]
    repl = paste0("Rep", gsub(".sort.bam","",repl) )

    type = ifelse(grepl("Input",basename(files)),"Input",
           ifelse(grepl("hmC",basename(files)),"hmC","mC"))

    tibble( file = files,cell , treatment = treat, type,  rep = repl)
}

df = create_structure(files) %>% filter(!grepl("Input",file))

dfInput_pool = create_structure(filesInput_pool) %>%
    filter(grepl("Input",file))

write_delim(df , path = "data/metadata/MeDIPseq_definition.tsv" ,
            delim = "\t")

write_delim(dfInput_pool , path = "data/metadata/MeDIPseq_definition_InputPooled.tsv" ,
            delim = "\t")
