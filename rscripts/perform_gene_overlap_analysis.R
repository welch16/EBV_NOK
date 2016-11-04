#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
  make_option("--data_1A", action = "store_true", type = "character",
              help = "File with the differentially expressed genes of dataset 1 in cell A"),
  make_option("--data_2A", action = "store_true", type = "character",
              help = "File with the differentially expressed genes of dataset 2 in cell A"),
  make_option("--data_1B", action = "store_true", type = "character",
              help = "File with the differentially expressed genes of dataset 1 in cell B"),
  make_option("--data_2B", action = "store_true", type = "character",
              help = "File with the differentially expressed genes of dataset 2 in cell B"),
  make_option("--outfile", action = "store_true", type = "character", default = tempfile(),
              help = "File where the results are going to be saved"),
  make_option("--cells", action = "store_true", type = "character",default = "A,B",
              help = "Names of the cell lines used separated by a ','. The default value is
                 'A,B'"),
  make_option("--datasets", action = "store_true", type = "character", default = "1,2",
              help = "Names of the datasets used separated by a ','. The default value is '1,2'"),
  make_option("--figs", action = "store_true", type = "character", default = "./Rplots",
              help = "Dir and prefix used to save the figures generated"))

opt = parse_args(OptionParser(option_list = optList))

library(base,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(readr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(tidyr,quietly = TRUE)
library(ggplot2,quietly = TRUE)

files = opt[grep("data_",names(opt))] %>% unlist

stopifnot(all(file.exists(files)))

results = lapply(files,read_tsv)

cells = strsplit(opt$cells,",") %>% unlist %>% rep(each = 2)
datasets = strsplit(opt$datasets,",") %>% unlist %>% rep(2)

fdr_values = c(0.01,0.05,0.1)
mean_count_values = c(10,50,100,250,500,750,1000)

clean_results <- function(result,cell_label,data_label)
{
  result %>% separate(id,into = c("ensembl","gene"),sep = "_",
                     extra = "merge") %>%
    mutate(cell = cell_label,
           dataset = data_label ) %>%
    select(gene,ensembl,baseMean,log2FoldChange,pval,padj,cell,
         dataset)
}

results = mapply(clean_results,results,cells,datasets,SIMPLIFY = FALSE)
results = do.call(rbind,results)

my_grid = expand.grid(fdr = fdr_values,min_count = mean_count_values) %>% as.tbl

apply_gene_op <- function(set1,set2,opfun)
{
  opfun(set1,set2)
}
  
overlap_by_cell <- function(result)
{
  sets = result %>% split(result$dataset)
  out = list()
  out[["both"]] = apply_gene_op(sets[[1]]$gene,
       sets[[2]]$gene,intersect) %>% length
  out[[paste0("only-",names(sets)[1])]] =
    apply_gene_op(sets[[1]]$gene,sets[[2]]$gene,setdiff) %>%
    length
  un = apply_gene_op(sets[[1]]$gene,sets[[2]]$gene,union) %>%
    length
  out[[paste0("only-",names(sets)[2])]] =
    un - Reduce(sum,out)
  unlist(out)
}


gene_overlap <- function(row,results,grid)
{
  
  fdr = my_grid[row,]$fdr
  count = my_grid[row,]$min_count

  results = results %>%
    filter(padj <= fdr & baseMean >= count)

  cell_results = split(results,results$cell)
  lapply(cell_results,overlap_by_cell) %>%
    unlist

}

process_overlap <- function(overlap,my_grid)
{
  overlap = do.call(rbind,overlap) %>% as.data.frame
  cbind(my_grid,overlap) %>% as.tbl %>%
    gather(venn,ngenes,c(-fdr,-min_count)) %>%
    mutate(venn = gsub(".","=",venn,fixed = TRUE)) %>%
    separate(venn,c("cell","venn"),"=")
}

overlap_all_genes = lapply(my_grid %>% nrow %>% seq_len,
  gene_overlap,results,my_grid) %>% process_overlap(my_grid)

overlap_up_genes = lapply(my_grid %>% nrow %>% seq_len,
  gene_overlap,results %>% filter(log2FoldChange >=  0), my_grid) %>%
  process_overlap(my_grid)
overlap_down_genes = lapply(my_grid %>% nrow %>% seq_len,
  gene_overlap,results %>% filter(log2FoldChange <  0), my_grid) %>%
  process_overlap(my_grid)


theme_set(theme_bw())

gene_plot <- function(overlap_genes)
{
  ggplot(overlap_genes,aes(min_count,ngenes,colour = cell)) + geom_line()+
    facet_grid(venn ~ fdr,scales = "free_y")+
      theme(legend.position = "top")+
        scale_color_brewer(palette = "Set1",name = "Cell")
}

overlap = do.call(rbind,
  list(
    overlap_all_genes %>% mutate(genes = "all"),
    overlap_up_genes %>% mutate(genes = "up"),
    overlap_down_genes %>% mutate(genes = "down")))

write_tsv(overlap,opt$outfile)

pdf(paste0(opt$figs,"_all_genes.pdf"))
u = print(gene_plot(overlap_all_genes))
dev.off()

pdf(paste0(opt$figs,"_up_genes.pdf"))
u = print(gene_plot(overlap_up_genes))
dev.off()

pdf(paste0(opt$figs,"_down_genes.pdf"))
u = print(gene_plot(overlap_down_genes))
dev.off()





