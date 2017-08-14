#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--all_files",action = "store_true",type = "character",
                help = "All RSEM outcomes to add into the matrix"),
    make_option("--out_file",action = "store_true",type = "character",
                default = tempfile(),
                help = "Name of the file where the TPM matrix is saved")
    
  )

opt = parse_args(OptionParser(option_list = optList))

## opt$all_files = "data/RSEM/hg19/*.genes.results"

separate_files <- function(ff)
{
  if(grepl(",",ff)){
      ff = ff %>% strsplit(',') %>% unlist     
  }else{
      ff = Sys.glob(ff)      
  }
  ff
}

opt$all_files = separate_files(opt$all_files)

library(tidyverse)

clean_names <- function(files)
{
    files %>%
        basename() %>%
        strsplit("\\.") %>%
        map_chr( ~ .[1])
}    



## load all files
rsem_data = opt$all_files %>%
    map( read_tsv)
names(rsem_data) = clean_names(opt$all_files)

## parse them into one matrix
rsem_data = rsem_data %>%
    map(select, contains("id"),TPM) %>%
    map2( names(.) , ~ mutate(.x , set = .y)) %>%
    bind_rows() %>%
    spread(set,TPM)


write_tsv(rsem_data,
          path = opt$out_file)

