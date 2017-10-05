
rm(list = ls())

library(tidyverse)
library(scales)

indr = "data/Diff.Genes/hg19/Sept17"
files = list.files(indr,full.names = TRUE)

meta = tibble(
  file = files %>% basename()
) %>% 
  mutate(
    cell = if_else(grepl("dRdZ",file),"EBV_dRdZ",
                   if_else(grepl("NOKS",file),"NOKS","EBV")),
    unit = if_else(grepl("genes",file),"gene","isoform"),
    file = map_chr(file , ~ gsub(paste0(".",tools::file_ext(.)),"",.))
  ) 
  

diff_genes = files %>% 
  set_names(meta$file) %>% 
  map(read_tsv) %>% 
  map2( names(.), ~ mutate(.x,file = .y))


theme_set(theme_bw())
pal = viridis::viridis(1e3)

## p.value histograms
diff_genes %>% 
  map(select,file,pvalue) %>% 
  bind_rows() %>% 
  inner_join(meta,by = "file") %>% 
  ggplot(aes(pvalue))+
  geom_histogram(boundary = 0,colour = "blue",fill = "white",bins = 31)+
  facet_grid( unit ~ cell , scales = "free_y")


fdr = 1e-10
## volcano plots
diff_genes %>% 
  map(select,file,log2FC,padj,pvalue) %>% 
  bind_rows() %>% 
  inner_join(meta, by = "file") %>%
  mutate(log10pval = -log10(pvalue),
         DE = if_else(padj <= fdr,"yes","no")) %>% 
  filter(!is.na(DE)) %>% 
  ggplot(aes(log2FC,log10pval,colour = DE))+
  geom_point(alpha = 1/2)+
  geom_vline(xintercept = 0,linetype = 2)+
  facet_grid(unit  ~ cell)+ylim(0,30)+
  theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = paste0("DE: padj <=",fdr))+
  ylab(expression(-log[10](p.value)))
  
  
