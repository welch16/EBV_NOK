#!/usr/bin/env Rscript

library(optparse)

## defines option parse and help caption for script
optList = list(
  make_option("--A" , action = "store_true",type = "character",
              help = "Files corresponding to A cell lines and undifferentiated genes"),
  make_option("--A_diff", action = "store_true",type = "character",
              help = "Files corresponding to A cell lines and differentiated genes"),
  make_option("--B", action = "store_true",type = "character",
              help = "Files corresponding to B cell lines and undifferentiated genes"),
  make_option("--B_diff",action = "store_true",type = "character",
              help = "Files corresponding to B cell lines and differentiated genes"),
  make_option(c("-t","--type"),action = "store_true",default = "rsem",type = "character",
              help = "Tool used to quantify transcripts. The default value is 'rsem'"),
  make_option(c("--xlab","-x"),action = "store_true",default = "",type = "character",
              help = "x-lab used in the plot"),
  make_option(c("--ylab","-y"),action = "store_true",default = "",type = "character",
              help = "y-lab used in the plot"),  
  make_option(c("--outputfile","-o"),action = "store_true",default = tempfile(),type = "character",
              help = "Directory and prefix used to store the plots of the analysis")
)

opt = parse_args(OptionParser(option_list = optList))

library(base)
opt$A = Sys.glob(opt$A)
opt$B = Sys.glob(opt$B)
opt$A_diff = Sys.glob(opt$A_diff)
opt$B_diff = Sys.glob(opt$B_diff)

# stops if files dont exists
stopifnot(all(file.exists(opt$A)),
          all(file.exists(opt$A_diff)),
          all(file.exists(opt$B)),
          all(file.exists(opt$B_diff)))

stopifnot(opt$type %in% c("none", "kallisto", "salmon", "sailfish","rsem"))

source("rfuns/geneExpression_visualization.R")

names(opt$A) = getRep(opt$A)
names(opt$A_diff) = getRep(opt$A_diff)
names(opt$B) = getRep(opt$B)
names(opt$B_diff) = getRep(opt$B_diff)

## load packages
library(tximport)
library(magrittr)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(readr)
library(scales)
library(hexbin)


## load files
A = tximport(opt$A,type = opt$type,importer = read_tsv)
A_diff = tximport(opt$A_diff,type = opt$type,importer = read_tsv)
B = tximport(opt$B,type = opt$type,importer = read_tsv)
B_diff = tximport(opt$B_diff,type = opt$type,importer = read_tsv)

countDT = createDataTable(A,A_diff,B,B_diff,"counts")
abundanceDT = createDataTable(A,A_diff,B,B_diff,"abundance")

theme_set(theme_bw())

sc = 5

pdf(file = paste0(opt$outputfile, "_abundance_hexbin_plot.pdf"))
geneExpression_count_plot(abundanceDT,"log2FC_A","log2FC_B",opt,sc)
dev.off()

pdf(file = paste0(opt$outputfile, "_counts_hexbin_plot.pdf"))
geneExpression_count_plot(countDT,"log2FC_A","log2FC_B",opt,sc)
dev.off()

