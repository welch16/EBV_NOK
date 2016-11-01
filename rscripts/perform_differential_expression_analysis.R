#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

## defines option parse and help caption for script
optList = list(
  make_option("--base" , action = "store_true",type = "character",
              help = "Files corresponding to A cell lines and undifferentiated genes"),
  make_option("--cond", action = "store_true",type = "character",
              help = "Files corresponding to A cell lines and differentiated genes"),
  make_option(c("-o","--outfile"),action = "store_true",type = "character",
              default = tempfile(),help = "File where the results are saved"),
  make_option(c("-t","--type"),action = "store_true",default = "rsem",type = "character",
              help = "Tool used to quantify transcripts. The default value is 'rsem'"),  
  make_option("--package",action = "store_true",type = "character",default = "ebseq",              
              help = "Package used to perform the analysis, the default value is 'ebseq'"),
  make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
              help = "Directory where the figures are being saved"),
  make_option("--fdr",action = "store_true",type = "numeric",default = 0.05,
              help = "FDR used to call differentially expressed genes")
 )


opt = parse_args(OptionParser(option_list = optList))

message("**************************************",opt$outfile)

library(base,quietly = TRUE)
opt$base = Sys.glob(opt$base)
opt$cond = Sys.glob(opt$cond)

library(tximport,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(readr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(tidyr,quietly = TRUE)
library(ggplot2,quietly = TRUE)

source("rfuns/geneExpression_analysis.R")

names(opt$base) = getRep(opt$base,"base")
names(opt$cond) = getRep(opt$cond,"treatment")

# stops if files dont exists
stopifnot(all(file.exists(opt$base)),
          all(file.exists(opt$cond)))
stopifnot(opt$type %in% c("none", "kallisto", "salmon", "sailfish","rsem"))
stopifnot(opt$package %in% c("deseq","ebseq","limma","edger"))


## load files with tximport
base = tximport(opt$base,type = opt$type,importer = read_tsv)
cond = tximport(opt$cond,type = opt$type,importer = read_tsv)

theme_set(theme_bw())

if(opt$package == "deseq2"){
  library(DESeq2)

  

}else if(opt$package == "deseq"){
  library(DESeq,quietly = TRUE)

  ## remove genes with low read count in treated and
  ## untreated condition
  toRemove = union( which(base[["counts"]] %>% rowGeoMeans == 0),
                    which(cond[["counts"]] %>% rowGeoMeans == 0)) %>% unique
  
  countData = cbind(base[["counts"]],cond[["counts"]]) %>% round
  condition = factor(c(rep("untreated",length(opt$base)),rep("treated",length(opt$cond))))

  ## remove genes with low reads
  countData = countData[-toRemove,]  
  
  ## specific to A samples
  countDataSet = newCountDataSet(countData,condition)
  countDataSet = estimateSizeFactors(countDataSet)
  countDataSet = estimateDispersions(countDataSet)

  ## tests to get differentially expressed genes
  results = nbinomTest(countDataSet,"untreated","treated")
  results = results %>% as.tbl

  plots = DESeq_plots(countDataSet , results,opt$fdr)
  dir.create(dirname(opt$figs),recursive = TRUE,showWarnings = FALSE)
  pdf(paste0(opt$figs,"_DESeq_baseMean_vs_dispersion.pdf"))
  u = print(plots[[1]])
  dev.off()
  pdf(paste0(opt$figs,"_DESeq_MA_plot.pdf"))
  u = print(plots[[2]])
  dev.off()
  pdf(paste0(opt$figs,"_DESeq_pvalue_histogram.pdf"),height = 5)
  u = print(plots[[3]])
  dev.off()
  
  results = cbind(results,countData)  
  results = results %>% as.tbl %>%
    mutate(dispersion = fData(countDataSet)[[1]])

  dir.create(dirname(opt$outfile),recursive = TRUE,showWarnings = FALSE)  
  write_delim(results,path = opt$outfile,delim = "\t")
  
}else if(opt$package == "ebseq"){
  library(EBSeq)
  
}else if(opt$package == "edger"){
  library(edgeR)

}else if(opt$package == "limma"){
  library(limma)


}else{
  stop("Package provided by user is not available")
}

