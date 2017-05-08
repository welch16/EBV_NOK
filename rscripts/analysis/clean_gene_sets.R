
rm(list = ls())

library(tidyverse)

genedr = "./data/Diff.Genes/hg19/DESeq2_contrasts/"
model = "cell:"

files = list.files(genedr,full.names = TRUE,recursive = TRUE)

files = files %>%
    {.[grep("DESeq2_c",.)]}

files = files[-3]  ## remove the dataset that considers EBV-mono-rep1 sample

library(parallel)

options(mc.cores = 12)

genes = files %>% map(read_tsv)

genes = genes %>%
  map(select,gene,contains(model))

nms = names(genes[[1]])
nms_new = gsub(model,"",nms)

genes = genes %>% 
  map(.f = function(x){
    names(x) = nms_new
    x
  }) %>% 
  map(separate,gene,into = c("ensembl_id","gene_id"),sep = "\\_")

which.genes = genes %>% 
  map(select,gene_id) %>% 
  bind_rows() %>% 
  group_by(gene_id) %>% 
  summarize(n = n()) %>% filter(n == length(genes)) %>% 
  select(gene_id)

genes = genes %>%
  map(right_join,which.genes,by = "gene_id") %>%
  map(mutate,padj = p.adjust(pvalue,method = "BH"))

list_genes = left_join(which.genes,genes[[1]][,1:2],by = "gene_id") %>%
  select(ensembl_id,gene_id)

files = gsub(".tsv","only_cell.tsv",files)

genes = genes %>% map2(files,write_tsv)


write_tsv(list_genes,
          path = "data/Diff.Genes/hg19/DESeq2_contrasts_genes_list.tsv")







