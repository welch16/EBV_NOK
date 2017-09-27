
rm(list = ls())

library(tidyverse)

isoform_tpm = "data/TPM_matrices/Isoforms_TPM_matrix.tsv" %>%
    read_tsv()

isoform_tpm %>%
    filter(grepl("ACTG1",gene_id)) %>%
    select(contains("id"),contains("Noks",ignore.case = FALSE)) %>%
    as.data.frame()

##               transcript_id               gene_id RNAseq-Noks_EBV-MC-rep1
## 1 ENST00000331925_ACTG1-201 ENSG00000184009_ACTG1                 3709.72
## 2 ENST00000447294_ACTG1-202 ENSG00000184009_ACTG1                    0.00
##   RNAseq-Noks_EBV-MC-rep2 RNAseq-Noks_EBV-MC-rep3 RNAseq-Noks_EBV-MC-rep4
## 1                  3430.7                 2357.74                 2705.01
## 2                     0.0                    0.00                    0.00
##   RNAseq-Noks_EBV-mono-rep1 RNAseq-Noks_EBV-mono-rep2 RNAseq-Noks_EBV-mono-rep3
## 1                   4405.00                   4147.43                   4564.56
## 2                      2.33                      0.00                      0.40
##   RNAseq-Noks_EBV-mono-rep4 RNAseq-Noks-MC-rep1 RNAseq-Noks-MC-rep2
## 1                   4356.51             2432.43             2651.05
## 2                      0.00              500.68                0.19
##   RNAseq-Noks-MC-rep3 RNAseq-Noks-MC-rep4 RNAseq-Noks-mono-rep1
## 1             2170.06             2495.06               4616.10
## 2                0.00                0.00                589.71
##   RNAseq-Noks-mono-rep2 RNAseq-Noks-mono-rep3 RNAseq-Noks-mono-rep4
## 1               4869.39               3368.44               3590.53
## 2                  2.28                  0.00                  0.00

results_iso = "data/Diff.Genes/hg19/Gene_Isoform/Scott_Isoforms_MC:EBV_vs_NOKS.tsv" %>%
    read_tsv()

results_gene = "data/Diff.Genes/hg19/Gene_Isoform/Scott_Genes_MC:EBV_vs_NOKS.tsv" %>%
    read_tsv()

results_iso %>%
    filter(grepl("ACTG1",gene_id)) %>%
    as.data.frame()

##               transcript_id               gene_id      log2FC log2FC_se
## 1 ENST00000331925_ACTG1-201 ENSG00000184009_ACTG1   0.1704526 0.2363626
## 2 ENST00000447294_ACTG1-202 ENSG00000184009_ACTG1 -28.8305987 4.0805438
##         stat       pvalue         padj rank_mean_tpm:EBV_MC
## 1  0.7211485 4.708181e-01 8.774822e-01                   23
## 2 -7.0653815 1.601753e-12 1.659529e-09                14594
##   rank_mean_tpm:NOKS_MC
## 1                    31
## 2                   888


results_gene %>%
    filter(grepl("ACTG1",gene_id)) %>%
    as.data.frame()


## > results_iso %>% nrow
## [1] 103727



## > results_gene %>%
## +     filter(grepl("ACTG1",gene_id)) %>%
## +     as.data.frame()
##                 gene_id                                    transcript_id(s)
## 1 ENSG00000184009_ACTG1 ENST00000331925_ACTG1-201,ENST00000447294_ACTG1-202
##      log2FC log2FC_se     stat    pvalue    padj rank_mean_tpm:EBV_MC
## 1 0.1138943 0.2204438 0.516659 0.6053942 0.77858                   44
##   rank_mean_tpm:NOKS_MC
## 1                    49
