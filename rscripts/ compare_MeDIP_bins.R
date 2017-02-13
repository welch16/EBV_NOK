#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
  make_option("--samples", action = "store_true",type = "character",
              help = "Files used to make the plots"),
  make_option("--sample_names", action = "store_true", type = "character",default = NULL,
              help = "Names of the files in '--samples'"),
  make_option("--bin_size",action = "store_true",type = "numeric",default = 200,
              help = "Bin size used to compare the files in samples"),
  make_option("--frag_len",action = "store_true",type = "numeric",default = 200,
              help = "Fragment length used to extend the reads"),
  make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
              help = "Directory where the figures are being saved")
)

opt = parse_args(OptionParser(option_list = optList))

library(base,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(GenomicRanges,quietly = TRUE)
library(GenomicRanges,quietly = TRUE)
library(parallel,quietly = TRUE)

separateFiles <- function(ff)
{
  if(grepl(",",ff)){
    ff = ff %>% strsplit(',') %>% unlist     
  }else{
    ff = Sys.glob(ff)      
  }
  ff
}


files = separateFiles(opt$samples)

if(is.null(opt$sample_names)){
  opt$sample_names = paste("Sample",seq_along(files) , sep = "-")
}

reads = files %>% lapply(readGAlignments,param = NULL) %>% lapply(as,"GRanges")

sizes = reads %>% lapply(as.data.frame) %>% lapply(as.tbl) %>% bind_rows %>%
  group_by(seqnames) %>% summarise(min = min(start),
                                   max = max(end))








