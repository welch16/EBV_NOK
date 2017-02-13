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

separateFiles <- function(ff)
{
  if(grepl(",",ff)){
    ff = ff %>% strsplit(',') %>% unlist     
  }else{
    ff = Sys.glob(ff)      
  }
  ff
}

options(cores = opt$cores)


files = separateFiles(opt$samples)

if(is.null(opt$sample_names)){
  opt$sample_names = paste("Sample",seq_along(files) , sep = "-")
}

reads = files %>% mclapply(readGAlignments,param = NULL) %>% mclapply(as,"GRanges")

sizes = read_delim(opt$size_files,delim ="\t",col_names = c("seqnames","coord"))

create_bins <- function(.)
{
  GRanges(seqnames = .$seqnames,
          ranges = IRanges(
            start = seq(1,.$size,by = opt$bin_size),
            width = opt$bin_size
          ))
}

bins = sizes %>% split(.$seqnames) %>% map(create_bins) %>%
  as.list %>% GRangesList %>% unlist


# 
# library(ggplot2)
# library(gridExtra)
# library(grid)
# library(gtable)
# 
# corr_1 = rnorm(100)
# corr_2 = rnorm(100)
# corr_12 = rnorm(100)
# corr_list = list(corr_1, corr_2, corr_12)
# ttls = c(
#   'variance within variable 1',
#   'correlation within variable 1 & 2',
#   'variance within variable 2'
# )
# plots = list()
# for (i in 1:3) {
#   temp_df = data.frame(x = corr_list[[i]])
#   temp = ggplot(data = temp_df, aes(x = x)) +
#     geom_density() +
#     
#     theme(axis.title = element_blank()) #+
#   
#   #  ggtitle(ttls[i])
#   plots[[i]] = temp
# }
# 
# ng <- nullGrob()
# gp <- arrangeGrob(plots[[1]], plots[[2]],
#                   ng, plots[[3]])
# 
# # The gp object is a gtable;
# # thus gtable functions can be applied to add the the necessary labels
# 
# # A list of text grobs - the labels
# vars <- list(textGrob("Variable 1"), textGrob("Variable 2"))
# 
# # So that there is space for the labels,
# # add a row to the top of the gtable,
# # and a column to the left of the gtable.
# gp <- gtable_add_cols(gp, unit(1.5, "lines"), 0)
# gp <- gtable_add_rows(gp, unit(1.5, "lines"), 0)
# 
# # Add the label grobs.
# # The labels on the left should be rotated; hence the edit.
# # t and l refer to cells in the gtable layout.
# # gtable_show_layout(gp) shows the layout.
# gp <-
#   gtable_add_grob(gp,
#                   lapply(vars, editGrob, rot = 90),
#                   t = 2:3,
#                   l = 1)
# gp <- gtable_add_grob(gp, vars, t = 1, l = 2:3)
# 
# # Draw it
# grid.newpage()
# grid.draw(gp)
