
rm( list = ls())

library(Vennerable)
library(tidyverse)

indr = "data/Diff.Genes/hg19/DESeq2_marginal"
files = indr %>% list.files(full.names = TRUE) %>%
    {.[grep("MC",.)]}

gene_list = files %>%
    map(read_tsv)
names(gene_list) = c("EBV","EBV-","EBV-dRdZ","EBV (Scott)","EBV- (Scott)")


create_venn_diagram <- function(thresh,gene_list)
{
    v1 = Venn(
        gene_list %>%
        map(filter,padj <= thresh) %>%
        map(select,gene_id) %>%
        map(.f = function(x)x[[1]])
    )


    plot(v1,show = list(Faces = FALSE))
              
}


thresh = c(1e-10,1e-8,1e-5,1e-3,1e-2,5e-2)

figsdr = "figs/diff_expression/DESeq2_marginal/Venn"

## 2 sets

thresh %>%
    map(.f = function(x,genes){
        pdf(file.path(figsdr,paste0("Venn_Diff_EBV_vs_NOKS_padj",x,".pdf")))
        create_venn_diagram(x,genes)
        dev.off()},
        gene_list[1:2])
        

thresh %>%
    map(.f = function(x,genes){
        pdf(file.path(figsdr,paste0("Venn_Diff_Up_EBV_vs_NOKS_padj",x,".pdf")))
        create_venn_diagram(x,genes)
        dev.off()},
        gene_list[1:2] %>%
        map(filter,log2FoldChange > 0))

thresh %>%
    map(.f = function(x,genes){
        pdf(file.path(figsdr,paste0("Venn_Diff_Down_EBV_vs_NOKS_padj",x,".pdf")))
        create_venn_diagram(x,genes)
        dev.off()},
        gene_list[1:2] %>%
        map(filter,log2FoldChange < 0))

## 3 sets

thresh %>%
    map(.f = function(x,genes){
        pdf(file.path(figsdr,paste0("Venn_Diff_EBV_vs_NOKS_vs_dRdZ_padj",x,".pdf")))
        create_venn_diagram(x,genes)
        dev.off()},
        gene_list[1:3])
        

thresh %>%
    map(.f = function(x,genes){
        pdf(file.path(figsdr,paste0("Venn_Diff_Up_EBV_vs_NOKS_vs_dRdZ_padj",x,".pdf")))
        create_venn_diagram(x,genes)
        dev.off()},
        gene_list[1:3] %>%
        map(filter,log2FoldChange > 0))

thresh %>%
    map(.f = function(x,genes){
        pdf(file.path(figsdr,paste0("Venn_Diff_Down_EBV_vs_NOKS_vs_dRdZ_padj",x,".pdf")))
        create_venn_diagram(x,genes)
        dev.off()},
        gene_list[1:3] %>%
        map(filter,log2FoldChange < 0))
