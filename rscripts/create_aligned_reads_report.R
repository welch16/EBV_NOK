#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-v","--verbose"),action = "store_true",default = TRUE,type = "logical",
              help = "Print extra output"),
  make_option(c("-d","--directory"),action = "store_true",default = "",type = "character",
              help = "Directory where the logs are stored"),
  make_option(c("-o","--outputfile"),action = "store_true",default = tempfile(),type = "character",
              help = "Output file with the alignment results table"),
  make_option(c("-a","--aligner"),action = "store_true",default = "bowtie",type = "character",
              help = "Aligner used. The only one supported is 'bowtie'")
)


opt <- parse_args(OptionParser(option_list = option_list))

if(opt$verbose){
  message("Using file in directory: ",opt$directory)
  message("Saving results in ",opt$outputfile)
  message("RSEM processed reads aligned with ",opt$aligner)
  message("Loading packages")

}
  
stopifnot(file.exists(opt$directory))
stopifnot(opt$aligner %in% c("bowtie"))

library(readr,quietly = opt$verbose)
library(dplyr,quietly = opt$verbose)
library(magrittr,quietly = opt$verbose)
library(tidyr,quietly = opt$verbose)
source("rfuns/parse_RSEM_bowtie_logs.R")

files = list.files(opt$directory,full.names = TRUE)

if(opt$aligner == "bowtie"){

  dt_list = lapply(files,parse_RSEM_bowtie_log)  

}

dt <- do.call(rbind,dt_list)

files <- basename(files)
files <- sapply(strsplit(files,".",fixed = TRUE),function(x)x[1])

dt$file = files

dt = dt %>% select(file,everything()) %>%
  mutate(aligned_prop = round(100 * aligned / processed,2))


write_delim(dt , opt$outputfile,delim = "\t")
