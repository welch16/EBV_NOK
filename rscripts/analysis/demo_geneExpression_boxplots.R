
rm(list = ls())

library(tidyverse)

tpmmat = read_tsv("data/metadata/RSEM_gene_TPM_matrix.tsv")


tpmmat = tpmmat %>% 
  separate(gene_id ,into = c("ensembl","gene"),sep = "\\_") %>%
  dplyr::rename(transcript = `transcript_id(s)`)

tpmmat %>% filter(gene == "BATF2") %>% select(contains("Scott")) %>% as.data.frame()


gene_exp <- function(gg,tpmmat,title = FALSE)
{
  mat = tpmmat %>% filter(gene == gg) %>% select(contains("Scott"))
  mat = mat %>% t 
  nms = mat %>% rownames
  mat = tibble(sample = gsub("Scott.","",nms), tpm = mat[,1] ) %>% 
    mutate(cell = ifelse(grepl("EBV",sample),"EBV","NOKS"))

  out = ggplot(mat,aes(cell,tpm))+geom_boxplot() +
    geom_point()+geom_text_repel(aes(label = sample),show.legend = FALSE)
  
  if(title) out = out + ggtitle(gg)

  out
}  

library(ggrepel)

# gene_exp("BATF2",tpmmat)
# gene_exp("AEN",tpmmat)
# gene_exp("ZNF790",tpmmat)
# gene_exp("PLA2G16",tpmmat)
# gene_exp("BBC3",tpmmat)
# gene_exp("FHAD1",tpmmat)
# gene_exp("PRKCH",tpmmat)
# 
# gene_exp("BRMS1",tpmmat)
# gene_exp("ZNF263",tpmmat)
# gene_exp("DPP3",tpmmat)
# gene_exp("C11orf84",tpmmat)

dipdr = "data/MeDIPseq_results"
dipmat = read_tsv(file.path(dipdr,"MeDIPseq_PromotersCounts_upstr500_downstr1000fraglen300.tsv"))
colnames(dipmat)[1] = "ensembl"

dipdepths = read_tsv("data/metadata/MeDIPseq_sequencingDepth.tsv")


## clean MeDIP-seq, removing input1,...,input3 because we
## pooled them into input0

dipdepths = dipdepths %>%
  filter(!grepl("Input-rep1",file),
         !grepl("Input-rep2",file),
         !grepl("Input-rep3",file))

dipmat = dipmat %>% select(-contains("Input-rep1"),
                           -contains("Input-rep2"),
                           -contains("Input-rep3"))

## remove samples with very low sequencing depths

minDepth = 10e6

w = which(dipdepths$depth > minDepth)
dipdepths = dipdepths %>% filter(depth > minDepth) %>%
  mutate(file = gsub(".sort.bam","",file))

dipmat = dipmat[, c(seq_len(4),4 + w)]

dipmat = dipmat %>% right_join(tpmmat[,1:2],by = "ensembl")


gene_methyl <- function(gg,dipmat,title = FALSE)
{
  mat = dipmat %>% filter(gene == gg) %>% select(-contains("Input")) %>%
    select(contains("hmC"))
  mat = mat %>% t 
  nms = mat %>% rownames
  mat = tibble(sample = gsub("MeDIPseq-","",nms), tpm = mat[,1] ) %>% 
    mutate(cell = ifelse(grepl("akata",sample),"EBV","NOKS"))
  
  out = ggplot(mat,aes(cell,tpm))+geom_boxplot() +
    geom_point()+geom_text_repel(aes(label = sample),show.legend = FALSE)
  
  if(title) out = out + ggtitle(gg)
  
  out

}  


# gene_methyl("BATF2",dipmat)
# gene_methyl("AEN",dipmat)
# gene_methyl("ZNF790",dipmat)
# gene_methyl("PLA2G16",dipmat)
# gene_methyl("BBC3",dipmat)
# gene_methyl("FHAD1",dipmat)
# gene_methyl("PRKCH",dipmat)
# 
# gene_methyl("BRMS1",dipmat)
# gene_methyl("ZNF263",dipmat)
# gene_methyl("DPP3",dipmat)
# gene_methyl("C11orf84",dipmat)

library(grid)
library(gridExtra)

gene_plot = function(gene,tpmmat,dipmat)
{
  p1 = gene_exp(gene,tpmmat) + ylab("Gene expression")
  p2 = gene_methyl(gene,dipmat) + ylab("Methylation")
  
  grid.arrange(p1,p2,top = gene,nrow = 1)
  
  
}

# gene_plot("FHAD1",tpmmat,dipmat)
# gene_plot("PRKCH",tpmmat,dipmat)
# 
# 
# gene_plot("C11orf24",tpmmat,dipmat)
# gene_plot("YIF1A",tpmmat,dipmat)
# gene_plot("EIF3M",tpmmat,dipmat)
# gene_plot("C11orf84",tpmmat,dipmat)
# gene_plot("C11orf84",tpmmat,dipmat)
# 
# 
# 
