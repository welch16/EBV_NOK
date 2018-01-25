
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(allez)

load(here::here("apps/Fig2_dashboad/fig2data.RData"))

## tidy outcome of DESeq2 hypothesis tests
diff_genes = diff_genes %>% 
  separate(gene_id , into = c("ENSEMBL_ID","SYMBOL"), sep = "\\_") 

## remove genes with duplicated symbols
dup_symbols = diff_genes %>% 
  group_by(SYMBOL) %>% 
  summarize(n = n()) %>% 
  filter(n > 1) %>% 
  pluck("SYMBOL")

diff_genes = diff_genes %>% 
  filter(!SYMBOL %in% dup_symbols)

## sanity check
diff_genes %>%
  ggplot(aes(pvalue))+
  geom_histogram(binwidth = 0.01,boundary = 0,
                 fill = "white",colour = "navyblue")

## ---------------------------------------------- allez

gene_scores1 = diff_genes %>% 
  pluck("stat") %>% 
  set_names(
    pluck(diff_genes,"SYMBOL")
  )

allez_go = allez(gene_scores1, "org.Hs.eg", idtype = "SYMBOL", library.loc = NULL,
      sets = "GO", locallist = NULL,
      universe = c("global", "local"),
      transform = c("none"), # c("none", "binary", "rank","nscore"),
      annotate = TRUE)

allez_kegg = allez(gene_scores1, "org.Hs.eg", idtype = "SYMBOL", library.loc = NULL,
                 sets = "KEGG", locallist = NULL,
                 universe = c("global", "local"),
                 transform = c("none"), # c("none", "binary", "rank","nscore"),
                 annotate = TRUE)





## ---------------------------------------------- clusterProfiler

