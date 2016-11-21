#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--A_noTr", action = "store_true",type = "character",
                help = "Files corresponding to A cell line and no treatment was applied"),
    make_option("--B_noTr", action = "store_true", type = "character",
                help = "Files corresponding to B cell line and no treatment was applied"),
    make_option("--A_Tr", action = "store_true", type = "character",
                help = "Files corresponding to A cell line and the treatment was applied"),
    make_option("--B_Tr", action = "store_true", type = "character",
                help = "Files corresponding to B cell line and the treatment was applied"),
    make_option(c("-o","--outfile"), action = "store_true",default = tempfile(),type = "character",
                help = "Name of the outfile where the output list of genes is saved"),
  make_option(c("-t","--type"),action = "store_true",default = "rsem",type = "character",
              help = "Tool used to quantify transcripts. The default value is 'rsem'"),
  make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
              help = "Directory where the figures are being saved"),
  make_option("--fdr",action = "store_true",type = "numeric",default = 0.05,
              help = "FDR used to call differentially expressed genes")
  )

opt = parse_args(OptionParser(option_list = optList))

library(base,quietly = TRUE)

opt$A_noTr = "data/RSEM/hg19/RNAseq-noks-no_treatment-rep?.genes.results"
opt$A_Tr = "data/RSEM/hg19/RNAseq-noks-methyl_cell-rep?.genes.results"
opt$B_noTr = "data/RSEM/hg19/RNAseq-akata-noks-no_treatment-rep?.genes.results"
opt$B_Tr = "data/RSEM/hg19/RNAseq-akata-noks-methyl_cell-rep?.genes.results"


opt$A_noTr = Sys.glob(opt$A_noTr)
opt$B_noTr = Sys.glob(opt$B_noTr)

opt$A_Tr = Sys.glob(opt$A_Tr)
opt$B_Tr = Sys.glob(opt$B_Tr)

library(tximport,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(readr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(tidyr,quietly = TRUE)
library(ggplot2,quietly = TRUE)

source("rfuns/geneExpression_analysis.R")

names(opt$A_noTr) = getRep(opt$A_noTr,"A_noTr")
names(opt$B_noTr) = getRep(opt$B_noTr,"B_noTr")
names(opt$A_Tr) = getRep(opt$A_noTr,"A_Tr")
names(opt$B_Tr) = getRep(opt$B_noTr,"B_Tr")

stopifnot(all(file.exists(opt$A_noTr)),
          all(file.exists(opt$B_noTr)),
          all(file.exists(opt$A_Tr)),
          all(file.exists(opt$B_Tr)))

stopifnot(opt$type %in% c("none", "kallisto", "salmon", "sailfish","rsem"))

## load files with tximport
A_noTr = tximport(opt$A_noTr,type = opt$type,importer = read_tsv)
B_noTr = tximport(opt$B_noTr,type = opt$type,importer = read_tsv)
A_Tr = tximport(opt$A_Tr,type = opt$type,importer = read_tsv)
B_Tr = tximport(opt$B_Tr,type = opt$type,importer = read_tsv)

theme_set(theme_bw())

library(DESeq2,quietly = TRUE)
