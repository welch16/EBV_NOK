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

getRep <- function(files)
{
  paste0("rep",seq_along(files))
}

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
A = tximport(opt$A,type = opt$type,reader = read_tsv)
A_diff = tximport(opt$A_diff,type = opt$type,reader = read_tsv)
B = tximport(opt$B,type = opt$type,reader = read_tsv)
B_diff = tximport(opt$B_diff,type = opt$type,reader = read_tsv)

## using geometric mean, we average per replicate and create a tibble
gemean <- function(x,...)exp(mean(log(x),...))
rowGeoMeans <- function(mat)apply(mat,1,gemean)

createDataTable <- function(A,A_diff,B,B_diff,what = c("abundance","counts"))
{

  genes = A[[what]] %>% rownames %>%
    strsplit("_",fixed = TRUE) %>%
    sapply(function(x)x[2])
  
  DT = tibble(
    genes = genes,
    A = rowGeoMeans(A[[what]]),
    A_diff = rowGeoMeans(A_diff[[what]]),
    B = rowGeoMeans(B[[what]]),
    B_diff = rowGeoMeans(B_diff[[what]]))

  DT = DT %>%
    mutate(log2FC_A = log2(1 + A_diff) - log2(1 + A),
           log2FC_B = log2(1 + B_diff) - log2(1 + B))

  DT
}

countDT = createDataTable(A,A_diff,B,B_diff,"counts")
abundanceDT = createDataTable(A,A_diff,B,B_diff,"abundance")

pal = viridis(1e3, option = "D")

theme_set(theme_bw())

sc = 5

pdf(file = paste0(opt$outputfile, "_abundance_hexbin_plot.pdf"))
abundanceDT %>% ggplot(aes(log2FC_A,log2FC_B))+stat_binhex(bins = 140)+
  scale_fill_gradientn(colours = pal,trans = "log10",
                labels = trans_format('log10',math_format(10^.x)),
                guide = guide_colourbar(title = NULL,
                  barheight = unit(.92,"npc"),
                  barwidth = unit(0.01,"npc")))+
  xlim(-sc,sc)+ylim(-sc,sc)+xlab(opt$xlab)+ylab(opt$ylab)
dev.off()

pdf(file = paste0(opt$outputfile, "_counts_hexbin_plot.pdf"))
countDT %>% ggplot(aes(log2FC_A,log2FC_B))+stat_binhex(bins = 140)+
  scale_fill_gradientn(colours = pal,trans = "log10",
                labels = trans_format('log10',math_format(10^.x)),
                guide = guide_colourbar(title = NULL,
                  barheight = unit(.92,"npc"),
                  barwidth = unit(0.01,"npc")))+
  xlim(-sc,sc)+ylim(-sc,sc)+xlab(opt$xlab)+ylab(opt$ylab)
dev.off()

