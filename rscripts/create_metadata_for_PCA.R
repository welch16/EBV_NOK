
rm(list = ls())

library(magrittr)
library(dplyr)
library(readr)

indir = "data/RSEM/hg19"
files = list.files(indir,full.names = TRUE,pattern = "genes.results")

lab = ifelse(grepl("Nok",basename(files)),"Scott","Rona")

cell = ifelse(lab == "Rona",
       ifelse(grepl("akata",basename(files)),"EBV","NOK"),
       ifelse(grepl("EBV",basename(files)),"EBV","NOK"))

treat = ifelse(lab == "Rona",
        ifelse(grepl("no_treatment",basename(files)),"None",
                     ifelse(grepl("methyl",basename(files)),"MC","CaFBS")),
        ifelse(grepl("mono",basename(files)),"None","MC"))

repl = basename(files) %>% strsplit("rep") %>% unlist
repl = repl[grepl("genes.results",repl)]
repl = paste0("Rep", gsub(".genes.results","",repl) )

df  = tibble( file = files,
             cell , treatment = treat,lab , rep = repl)

write_delim(df , path = "data/metadata/PCA_definition.tsv" ,
            delim = "\t")
