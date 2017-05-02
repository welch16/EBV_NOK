
rm(list = ls())

library(tidyverse)

genedr = "/store01/Collab_EJohannsen_RWelch/EBV_NOK/data/Diff.Genes/hg19"

files = list.files(genedr,full.names = TRUE,recursive = TRUE)

files = files %>%
    {.[grep("DESeq2_c",.)]}

library(parallel)

options(mc.cores = 12)

genes = mclapply(files[-3],read_tsv)

genes = genes %>%
    map(.f = function(x) x %>% select(gene)) %>%
    {do.call(intersect,.)} %>%
    separate(gene,into = c("ensembl_id","gene_id"),sep = "\\_")
    
write_tsv(genes,
          path = "data/Diff.Genes/hg19/DESeq2_contrasts_genes_list.tsv")
