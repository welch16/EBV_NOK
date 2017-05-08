
rm(list = ls())

library(tidyverse)
library(ChIPseeker)
library(Vennerable)

indr =  "data/Diff.Genes/hg19/DESeq2_contrasts"
files = list.files(indr,
                   full.names = TRUE,
                   pattern = "cell")

genes = files %>% map(read_tsv)
names(genes) = c("Ben_CaFBS","Ben_MC","Scott_MC")


# Summary stats -----------------------------------------------------------


incidence_matrix = function(genes,padj_thresh)
{
  diff_genes = genes %>%
    map(filter,padj <= padj_thresh)
  
  all_genes = diff_genes %>% bind_rows() %>%
    select(gene_id) %>% unique()
  
  adj = diff_genes %>% map(.f = function(x,all_genes){
    ifelse(all_genes[[1]] %in% x$gene_id,
           1, 0)
  },all_genes) %>% bind_cols() %>% as.tbl()
  
  bind_cols(all_genes,adj) %>%
    mutate( total = rowSums(.[,-1]))
}


# Plots -------------------------------------------------------------------

vennplot_th = function(genes,padj_thresh)
{
  diff_genes = genes %>%
    map(filter,padj <= padj_thresh)

  diff_genes %>% map(select,gene_id) %>% 
    map(.f = function(x)x$gene_id) %>%
    vennplot(by = "Vennerable")

}


# run ---------------------------------------------------------------------

figsdr = "figs/Venn_genes"
dir.create(figsdr,showWarnings = TRUE)

create_venn_plot = function(genes,padj_thresh)
{
  pdf(file.path(figsdr,paste0("Venn_diagram_padj",padj_thresh ,".pdf")))
  vennplot_th(genes,padj_thresh) %>% print()
  dev.off()
}

pvals = c(1e-4,1e-3,1e-2,5e-2)

pvals %>% map(.f = function(x)create_venn_plot(genes,x))



# google_sheets -----------------------------------------------------------

library(googlesheets)

diff_genes = genes %>% map(filter,padj <= 1e-1)
nms = names(diff_genes)

dr = "data/local"

ff = basename(files) %>% {gsub(".tsv",".csv",.)} %>% {file.path(dr,.)}

diff_genes %>% map2(ff,.f = function(x,y)write_csv(x,path =y ))


gs_upload(ff[1],"Ben_CaFBS_diff_EBV")
gs_upload(ff[2],"Ben_MC_diff_EBV")
gs_upload(ff[3],"Scott_MC_diff_EBV")



# incidence matrices ------------------------------------------------------

incs  = pvals %>% map(.f = function(p)incidence_matrix(diff_genes,p))
pfiles = file.path(dr,paste0("Incidence_matrix_",pvals,".csv"))

incs %>% map2(pfiles,.f = function(x,y)write_csv(x,path = y))

gs_upload(pfiles[1],paste0("Incidence_matrix_",pvals[1]))
gs_upload(pfiles[2],paste0("Incidence_matrix_",pvals[2]))
gs_upload(pfiles[3],paste0("Incidence_matrix_",pvals[3]))
gs_upload(pfiles[4],paste0("Incidence_matrix_",pvals[4]))

