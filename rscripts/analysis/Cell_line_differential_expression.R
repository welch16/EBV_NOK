
rm(list = ls())

library(tidyverse)
library(DESeq2)
library(BiocParallel)
library(tximport)
library(pheatmap)

options(mc.cores = 12)

indr = "./data/RSEM/hg19"
files = indr %>%
    list.files(full.names = TRUE,pattern = "genes.results")

## use Scott data and new mutant type


files = c(
    files %>%
        {.[grep("Noks",files)]},
    files  %>%
        {.[grep("clone",files)]})


txdata = tximport(files,type = "rsem" , importer = read_tsv)

bfiles = files %>% basename()

coldata = data.frame(
    strain = ifelse(grepl("Noks_EBV",bfiles),"EBV",
             ifelse(grepl("Noks",bfiles),"NOKS","dRdZ")),
    treat = ifelse(grepl("MC",bfiles) | grepl("methy",bfiles),
                   "methyl","mono")) %>%
    mutate(
       interac = factor(paste(strain,treat,sep = "_"))) %>%
    arrange(strain,treat)

nms = coldata %>%
    group_by(interac) %>%
    summarize( n = n()) %>%
    ungroup()

rownames(coldata) = map2( nms$interac,
                         nms$n,
                         .f = function(x,i) paste0(x,seq_len(i))) %>%
    unlist()


DD = DESeqDataSetFromMatrix(
    round(txdata[["counts"]]),colData = coldata,
                               design= ~ strain + treat + strain:treat )

DD = DD[rowSums(counts(DD)) > 1,]
DDmodel1 = DESeq(DD,minReplicatesForReplace = Inf)
design(DD) = ~ interac
DDmodel2 = DESeq(DD,minReplicatesForReplace = Inf)


strain_results = list(
    ## dRdZ vs NOKS
    "strain:dRdZ_vs_NOKS" = results(DDmodel1,cooksCutoff = FALSE,
                                    contrast = c("strain","dRdZ","NOKS"),
                                    tidy = TRUE,
                                    parallel = TRUE),
    ## EBV vs NOKS
    "strain:EBV_vs_NOKS" = results(DDmodel1,cooksCutoff = FALSE,
                                   contrast = c("strain","EBV","NOKS"),
                                   tidy = TRUE,
                                   parallel = TRUE),
    ## dRdZ vs NOKS when treat = methyl
    "methyl:dRdZ_vs_NOKS" = results(DDmodel2,cooksCutoff = FALSE,
                                    contrast = c("interac","dRdZ_methyl","NOKS_methyl"),
                                    tidy = TRUE,
                                    parallel = TRUE),
    ## EBV vs NOKS when treat = methyl
    "methyl:EBV_vs_NOKS" = results(DDmodel2,cooksCutoff = FALSE,
                                   contrast = c("interac","EBV_methyl","NOKS_methyl"),
                                   tidy = TRUE,
                                   parallel = TRUE),
    ## dRdZ vs NOKS without treatment
    "mono:dRdZ_vs_NOKS" = results(DDmodel2,cooksCutoff = FALSE,
                                  contrast = c("interac","dRdZ_mono","NOKS_mono"),
                                  tidy = TRUE,
                                  parallel = TRUE),
    ## EBV vs NOKS without treatment
    "mono:EBV_vs_NOKS" = results(DDmodel2,cooksCutoff = FALSE,
                                 contrast = c("interac","EBV_mono","NOKS_mono"),
                                 tidy = TRUE,
                                 parallel = TRUE))

## clean format
strain_results = strain_results %>%
    map(as_tibble) %>%
    map(separate,row,into = c("ENSEMBL_ID","GENE_ID"),sep = "\\_")




## plot

pvalue_histogram <- function(res,mod)
{
  
    res %>% ggplot(aes(pvalue)) +
        geom_histogram(bins = 50,colour = "blue",fill = "white",boundary = 0) +
        ggtitle(mod)

}


figsdr = "figs/diff_expression/DESeq2_strain"
ff = file.path(figsdr,"DESeq2_strain")



pdf(paste0(ff,"_pvalue_histograms.pdf"),height = 5)
u = strain_results %>%
    map2(names(.),pvalue_histogram) %>%
    map(print)
dev.off()


volcano_plot <- function(res,mod,fdr = 0.05)
{

    res = res %>%
        mutate(col = ifelse(padj <= fdr,"yes","no")) %>%
        filter(!is.na(padj))

    mm = 1.3 * min(min(res$log2FoldChange),max(res$log2FoldChange)) %>%
        abs()

    res %>% ggplot(aes(log2FoldChange,-log10(pvalue),colour = col))+geom_point()+
        scale_colour_brewer(name = paste("adjusted pval <=",fdr),palette = "Set1")+
        theme(legend.position = "top")+ggtitle(mod)+xlim(-mm,mm)+
        geom_vline(xintercept = 0,linetype = 2)


}

pdf(paste0(ff,"_volcano_plots.pdf"))
u = strain_results %>%
    map2(names(.),
         volcano_plot,fdr = 1e-5) %>%
    map(print)
dev.off()

signal2noise_barplot <- function(res,mod)
{    

    res = res %>%
        arrange(stat) %>%
        mutate(rank = seq_len(nrow(.)))

    res %>% ggplot(aes(rank,stat))+geom_line()+
        ggtitle(mod)+geom_abline(slope = 0,intercept = 0,linetype = 2)
}

pdf(paste0(ff,"_signal_to_noise_plot.pdf"))
u = strain_results %>%
    map2(names(.),signal2noise_barplot) %>%
    map(print)
dev.off()


rld = rlog(DDmodel1,blind = FALSE)

pdf(paste0(ff,"_PCA_analysis.pdf"))
u = plotPCA(rld,intgroup = c("strain","treat")) %>%
    print()
dev.off()


pheat_pval <- function(rld,my_results,K,suffix)
{
    mat = assay(rld)[,]
    mat = mat - rowMeans(mat)
    rownames(mat) = my_results$GENE_ID
    r = viridis::viridis(1e3,option = "D")

    coldata = colData(rld)
    coldata$sizeFactor = NULL

    coldata = coldata %>% as.data.frame
    
    topgenes = head(order(my_results$pvalue),K)
    pheatmap(mat[topgenes,],color = r,
             annotation_col = coldata,
             show_rownames = K <= 75,
             filename = file.path(figsdr,paste0("Heatmap_bottom",K,"genes_pval_",suffix,".pdf")),
             height = 9,width = 6)
    
}    


nms = strain_results %>%
    names()
Ks = c(20,50,100,5e2,1e3,5e3,1e4)

for( i in seq_along(nms)){
    a = map(Ks,function(x)pheat_pval(rld,strain_results[[i]],x,nms[i]))
}    

