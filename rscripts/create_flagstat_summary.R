#!/usr/bin/env Rscript

## this script is used to compare a genomic region

library(optparse,quietly = TRUE)

optList = list(
    make_option("--logsdir",action = "store_true",type = "character",
                help = "Log directory with the files to summarized with samtools flagstat output"),
    make_option("--outfile",action = "store_true",type = "character",default = tempfile(),
                help = "Output file with a table with the following columns:
                          Total reads | Mapped reads | Perc. mapped reads")
  )

opt = parse_args(OptionParser(option_list = optList))

library(readr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(tidyr,quietly = TRUE)

files = list.files(opt$logsdir,full.names = TRUE)

parse_file <- function(file)
{
    string = (file %>% read_file %>% strsplit("\\n"))[[1]]   

    tbl =tibble(string) %>% separate(string,into = c("q","rest"),sep = "\\+") %>%
        mutate(q = as.numeric(q)) %>% filter(q > 0)

    out = tibble(file = basename(file),total = tbl$q[1],mapped = tbl$q[2]) %>%
        mutate( mapped_perc = mapped / total,
               file = gsub(".logs","",file))      
    out 
}

parses = lapply(files,parse_file) %>% bind_rows

write_delim(parses,path = opt$outfile,delim = ",")
