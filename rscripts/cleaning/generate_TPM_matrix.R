
rm(list = ls())


library(readr)
library(matrixStats)
library(tidyverse)
library(purrr)


## first step is to load all the results obtained by RSEM

indr = "data/RSEM/hg19"
files = list.files(indr,full.names = TRUE,pattern = "genes.results")

RSEM_results = map(files,read_delim,"\t") 


## define variable in each files and names of files

vars = RSEM_results[[1]] %>% colnames
vars = vars[grep("id",vars,invert = TRUE)]
samples = files %>% basename %>% {gsub(".genes.results","",.)}

genes = RSEM_results[[1]][,1:2]

metadata = tibble(sample = samples) %>%
    mutate(
        EBV = ifelse(grepl("akata",sample) | grepl("EBV",sample),"EBV","NOKS"),
        Treatment = ifelse(grepl("CaFBS",sample),"CaFBS",
                    ifelse(grepl("mono",sample) | grepl("no_treat",sample),"NoTr",
                           "MC")),
        Replicate = strsplit(sample,"-") %>% map_chr(.f = function(x)x[length(x)]) %>%
            {gsub("rep","Rep",.)},
        lab = ifelse(grepl("Noks",sample),"Scott","Ben"),
        name = paste(lab,EBV,Treatment,Replicate,sep = ".")
        
    ) %>% select(name,sample,everything())


output = lapply(vars,
                function(x){
                    RSEM_results %>%
                        map(.f = function(z){
                            z %>% select_(x)}) %>% bind_cols}) %>%
    map(.f = function(z){
        colnames(z) = metadata$name
        z %>% bind_cols(genes) %>% select(contains("id"),everything())
        }  )
names(output) = vars


outdr = "data/metadata"
outfiles = file.path(outdr,paste0("RSEM_gene_",vars,"_matrix.tsv"))


map2(output,outfiles,write_tsv)





               

