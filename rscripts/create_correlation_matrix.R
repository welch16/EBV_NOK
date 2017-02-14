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
  make_option("--size_file",action = "store_true",type = "character",
              help = "Name of chromosome sizes file"),
  make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
              help = "Directory where the figures are being saved"),
  make_option("--use_log",action = "store_true",type = "logical",default = FALSE,
              help = "Logical value indicating if the data should be plotted in log-scale"),
  make_option("--cores",action = "store_true",type = "numeric",default = 8,
              help = "Number of cores used for parallel")
)

opt = parse_args(OptionParser(option_list = optList))

library(base,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(purrr,quietly = TRUE)
library(GenomicAlignments,quietly = TRUE)
library(GenomicRanges,quietly = TRUE)
library(parallel,quietly = TRUE)
library(readr,quietly = TRUE)

options(cores = opt$cores)

source("rfuns/MeDIP_analysis.R")

opt$samples = "data/BAM/hg19/bowtie_dip/MeDIPseq-NOKS-akata-mono-Input-rep?.sort.bam"
opt$size_file = "/p/keles/SOFTWARE/hg19.chrom.sizes"

opt$samples = separateFiles(opt$samples)

## opt$sample_names = paste0(
##     paste0("Input",seq_len(length(files))),collapse = ",")

if(is.null(opt$sample_names)){
  opt$sample_names = paste0("Replicate",seq_along(opt$samples) )
}else{
  opt$sample_names = strsplit(opt$sample_names,",") %>% unlist
  
}

reads = opt$samples %>% mclapply(readGAlignments,param = NULL) %>% mclapply(as,"GRanges")
names(reads) = opt$sample_names

sizes = read_delim(opt$size_file,delim ="\t",col_names = c("seqnames","size"))

bins = sizes %>% split(.$seqnames) %>% map(create_bins,opt$bin_size) %>%
  as.list %>% GRangesList %>% unlist
names(bins) = NULL

reads = reads %>% mclapply(resize,opt$frag_len)

bin_counts = reads %>% mclapply(function(x)countOverlaps(bins,x)) %>%
  as.data.frame  %>% as.tbl

if(opt$use_log){
  bin_counts = bin_counts %>%
    mutate_all(funs(log10(1 + . )))
}

N = ncol(bin_counts)

library(ggplot2,quietly = TRUE)
library(gridExtra,quietly = TRUE)
library(grid,quietly = TRUE)
library(gtable,quietly = TRUE)
library(hexbin,quietly = TRUE)


nms = opt$sample_names
