#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--design_file",action = "store_true",type = "character",
                help = "tsv file with the following columns: File | AKATA | Treatment | Type | Rep "),
  make_option("--bin_size",action = "store_true",type = "numeric",default = 200,
              help = "Bin size used to compare the files in samples"),
  make_option("--frag_len",action = "store_true",type = "numeric",default = 200,
              help = "Fragment length used to extend the reads"),
  make_option("--size_file",action = "store_true",type = "character",
              help = "Name of chromosome sizes file"),  
  make_option("--outfile",action = "store_true",type = "character",default = tempfile(),
              help = "Name of the file where the bin counts are stored."),
  make_option("--cores",action = "store_true",type = "numeric",default = 8,
              help = "Number of cores used for parallel")
)

opt = parse_args(OptionParser(option_list = optList))

library(base,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(GenomicAlignments,quietly = TRUE)
library(GenomicRanges,quietly = TRUE)
library(parallel,quietly = TRUE)
library(readr,quietly = TRUE)
library(tidyr,quietly = TRUE)
library(purrr,quietly = TRUE)

source("rfuns/MeDIP_analysis.R")

options(cores = opt$cores)

sizes = read_delim(opt$size_file,delim ="\t",col_names = c("seqnames","size"))

bins = sizes %>% split(.$seqnames) %>% map(create_bins,opt$bin_size) %>%
  as.list %>% GRangesList %>% unlist
names(bins) = NULL



design = read_delim(opt$design_file,delim = "\t")

split_design = design %>% split(interaction(.$rep,.$cell))

design = design %>% mutate(name = "MeDIPseq") %>%
    unite("name",name,cell,treatment,type,rep,sep = "_")

split_design = split_design %>% lapply(merge,design,by = "file")

##
##split_design = split_design %>% lapply(function(x)x[1:2,])
##split_design = split_design[1:2]


calculate_bin_counts <- function(desi,bins)
{
    print(desi)
    reads = desi$file %>% mclapply(readGAlignments,param = NULL) %>%
        mclapply(as,"GRanges") %>% mclapply(resize,opt$frag_len)
    names(reads) = desi$name

    ## get counts
    reads %>% mclapply(function(x)countOverlaps(bins,x)) %>%
        as.data.frame %>% as.tbl    
}    

split_bins = split_design %>% lapply(calculate_bin_counts,bins)

bin_counts = bind_cols(split_bins)


bindata = bind_cols(
    bins %>% as.data.frame %>% as.tbl,
    bin_counts
   ) %>%
    mutate(strand = NULL,
           width = NULL,
           end = NULL)

write_delim(bindata ,path = opt$outfile ,delim = "\t")
