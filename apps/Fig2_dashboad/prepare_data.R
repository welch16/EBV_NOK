


library(tidyverse)


diff_genes = "data/Diff.Genes/hg19/Gene_Isoform/Scott_Genes_mono:EBV_vs_NOKS.tsv" %>% 
  read_tsv() %>% 
  separate(
    gene_id, into = c("ensembl_id","gene_id"), sep = "\\_"
  ) %>% 
  select(-contains("transcript"))

tpm_mat = "data/TPM_matrices/Genes_TPM_matrix.tsv" %>% 
  read_tsv() %>% 
  select(gene_id,contains("Noks",ignore.case = FALSE)) %>% 
  separate(
    gene_id, into = c("ensembl_id","gene_id"), sep = "\\_"
  ) 

genes = intersect(
  diff_genes %>% pluck("ensembl_id"),
  tpm_mat %>% pluck("ensembl_id")
)

diff_genes = diff_genes %>% 
  filter(ensembl_id %in% genes)

tpm_mat = tpm_mat %>% 
  filter(ensembl_id %in% genes)

save(diff_genes,tpm_mat,file = "./apps/Fig2_dashboad/fig2data.RData")