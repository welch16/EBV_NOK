
rm(list = ls())

library(tidyverse)
library(BiocParallel)

bp = MulticoreParam(workers = 10,bpprogressbar = TRUE)


upstream = 500
downstream = 1000
fraglen = 300


dr = "data/MeDIPseq_results"
promoters = read_tsv(file.path(dr,
                           paste0("MeDIPseq_PromotersCounts_upstr",upstream,
                                  "_downstr",downstream,"fraglen",fraglen,".tsv")))
    

genes = read_tsv("data/Diff.Genes/hg19/DESeq2_contrasts_genes_list.tsv")

genes = genes %>%
    rename(ensembl_gene_id = ensembl_id)

promoters = promoters %>%
    right_join(genes,by = "ensembl_gene_id")

dipdepths = read_tsv("data/metadata/MeDIPseq_sequencingDepth.tsv")

mindepth = 1e7

dipdepths = dipdepths %>%
    filter(depth > mindepth) %>%
    filter(!grepl("Input",file)) %>%
    filter(!grepl("rep0",file)) %>%
    mutate(
        file = gsub(".sort.bam","",file)
    )

base = promoters %>%
    select(ensembl_gene_id,
           gene_id,
           seqnames,
           start,
           end)

promoters = bind_cols(
    base,
    promoters[,dipdepths$file]
    )


clean_cols <- function(geneList)
{

    X = geneList[,-seq_len(5)]
    
    toremove = X %>% as.list %>%
        map(.f = function(x)is.nan(x) | is.na(x)) %>%
        map(which) %>% do.call(c,.) %>% unique

    geneList[-toremove,]


}

promoters = promoters %>%
    clean_cols()


base = promoters[,seq_len(5)]
count_mats = list()
count_mats[["mC"]] = promoters[,-seq_len(5)] %>%
    select(-contains("hmC")) %>% as.matrix()
count_mats[["hmC"]] = promoters[,-seq_len(5)] %>%
    select(contains("hmC")) %>% as.matrix()

create_coldata <- function(count_mat)
{
    nms = count_mat %>% colnames

    out = data.frame(
        cell = ifelse(grepl("akata",nms),"EBV","NOKS"))
    rownames(out) = nms

    out


}

coldata_list = count_mats %>%
    map(create_coldata)

## DESeq2 analysis

library(DESeq2)

create_deseqData <- function(count_mat,coldata)
{   
    deseqData = DESeqDataSetFromMatrix(countData = count_mat,
                                 colData = coldata,
                                 design = ~ cell)

    deseqData
    

}


deseq_data_list = count_mats %>%
    map2(coldata_list,create_deseqData)



deseq_analysis <- function(dds)
{
    res = results(dds,cooksCutoff = FALSE,contrast = c("cell","EBV","NOKS"),
                  tidy = FALSE)

    res %>%  as.data.frame %>%
        as.tbl %>%
        dplyr::rename(p.value = pvalue)  

    

}

deseq_model_list = deseq_data_list %>%
    map(DESeq)

deseq_results_list = deseq_model_list %>%
    map(deseq_analysis) %>%
    map(.f = function(x)bind_cols(base,x))

save(deseq_model_list,file = file.path(dr,"DESeq2_diff_methylation_model.RData"))
save(deseq_results_list,file = file.path(dr,"DESeq2_diff_methylation_results.RData"))


## edgeR analysis
library(edgeR)


create_edgeRData <- function(count_mat,coldata)
{


    DGEList(count_mat,group = coldata$cell)
    

}


edgeR_data_list = count_mats %>%
    map2(coldata_list,create_edgeRData)

edgeR_model_list = edgeR_data_list %>%
    map(estimateCommonDisp) %>%
    map(estimateTagwiseDisp,prior.df = 10)


edgeR_results_list = edgeR_model_list %>%
    map(exactTest,pair = c("EBV","NOKS")) %>%
    map(topTags,n = nrow(base)) %>%
    map(.f = function(x)x[[1]] %>% as.tbl)
        

save(edgeR_model_list,file = file.path(dr,"edgeR_diff_methylation_model.RData"))
save(edgeR_results_list,file = file.path(dr,"edgeR_diff_methylation_results.RData"))




## mC_deseqResults_list %>% map2(names(mC_deseqResults_list),
##                               .f = function(x,y){
##                                   write_tsv(x,
##                                             path = file.path(outdr,
##                                                              paste0(y,"_diffMethyl_mC.tsv")))})
 
## hmC_deseqResults_list %>% map2(names(hmC_deseqResults_list),
##                               .f = function(x,y){
##                                   write_tsv(x,
##                                             path = file.path(outdr,
##                                                              paste0(y,"_diffMethyl_hmC.tsv")))})



## rld_mC = mC_deseqModel_list %>% map(rlog,blind = FALSE)
## rld_hmC = hmC_deseqModel_list %>% map(rlog,blind = FALSE)

## theme_set(theme_bw())


## rowVars <- function (x,na.rm = TRUE)
## {
##     sqr = function(x) x * x
##     n = rowSums(!is.na(x))
##     n[n <= 1] = NA
##     return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
## }


## PCA_plot <- function(rl,tit,ntop = 500,x = 1,y = 2,tab)
## {

##     rv = rowVars(assay(rl))

##     select = order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]

##     pca = prcomp(t(assay(rl)[select, ]))
##     percentVar = pca$sdev^2/sum(pca$sdev^2)

##     intgroup.df =  as.data.frame(colData(rl)[, ,drop = FALSE]) %>% as.tbl

##     dt = intgroup.df %>% mutate(file = tab[[1]],
##                                 x = pca$x[,x], y = pca$x[,y]) %>%
##         mutate(cond = ifelse(grepl("mono",file),"No-Tr","CaFBS")) %>% split(.$cell)  %>%
##         map(.f = function(z)
##             z %>% mutate(Rep = paste("Rep",seq_len(nrow(z)),sep = "-"))) %>%
##         bind_rows  

##     dt %>%
##         ggplot(aes(x,y,colour = cell))+
##         geom_point(size = 2,shape = 1)+coord_fixed()+
##         xlab(paste0("PC",x,": ",round(percentVar[x] * 100), "% variance"))+
##         ylab(paste0("PC",y,": ",round(percentVar[y] * 100), "% variance"))+
##         scale_color_brewer(palette = "Set1",name = "Cell")+        
##         geom_label_repel(aes(x,y,label = file),size = 2,
##                          box.padding = unit(.3,"lines"),
##                          point.padding = unit(.7,"lines"),
##                          show.legend = FALSE)+ ggtitle(tit)+
##         theme(legend.position = "top")

## }


## library(ggrepel)

## mC_PCA_plots = rld_mC %>%
##     map2(names(rld_mC),PCA_plot,
##          tab = dipdepths %>% filter(!grepl("Input",file)) %>%
##              filter(!grepl("hmC",file))
##          )

         
## hmC_PCA_plots = rld_hmC %>%
##     map2(names(rld_hmC),PCA_plot,
##          tab = dipdepths %>% filter(!grepl("Input",file)) %>%
##              filter(grepl("hmC",file))
##          )



## figsdr = "figs/diff_methylation/DESeq2_results"


## pdf(file.path(figsdr,"DESeq2_PCA_mC.pdf"))
## u = mC_PCA_plots %>% map(print)
## dev.off()



## pdf(file.path(figsdr,"DESeq2_PCA_hmC.pdf"))
## u = hmC_PCA_plots %>% map(print)
## dev.off()


## library(pheatmap)

## rowVars <- function (x,na.rm = TRUE)
## {
##     sqr = function(x) x * x
##     n = rowSums(!is.na(x))
##     n[n <= 1] = NA
##     return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
## }



## pheat_pval <- function(rld,my_results,K,suffix)
## {
##     mat = assay(rld)[,]
##     mat = mat - rowMeans(mat)
##     rownames(mat) = my_results$gene
##     r = viridis::viridis(1e3,option = "D")

##     coldata = colData(rld)
##     coldata$sizeFactor = NULL

##     coldata = coldata %>% as.data.frame
    
##     topgenes = head(order(my_results$p.value),K)
##     pheatmap(mat[topgenes,],color = r,
##              annotation_col = coldata,
##              show_rownames = K <= 75,
##              filename = file.path(figsdr,paste0("Heatmap_bottom",K,"genes_pval_",suffix,".pdf")),
##              height = 9,width = 6)
    
## }    

## pheat_var <- function(rld,my_results,K,suffix)
## {
##     mat = assay(rld)[,]
##     mat = mat - rowMeans(mat)
##     rownames(mat) = my_results$gene
##     r = viridis::viridis(1e3,option = "D")

##     coldata = colData(rld)
##     coldata$sizeFactor = NULL

##     coldata = coldata %>% as.data.frame

    

##     topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)

##     pheatmap(mat[topgenes,],annotation_col = coldata,color = r,
##              show_rownames = K <= 75,
##              filename = file.path(figsdr,paste0("Heatmap_top",K,"genes_var_",suffix,".pdf")),
##              height = 9,width = 6)    
## }    


## Ks = c(20,50,100,5e2,1e3,5e3,1e4)

## a = lapply(Ks,function(x)pheat_pval(rld_mC[[3]],mC_deseqResults_list[[3]],x,"mC"))
## b = lapply(Ks,function(x)pheat_var(rld_mC[[3]],mC_deseqResults_list[[3]],x,"mC"))

## a = lapply(Ks,function(x)pheat_pval(rld_hmC[[3]],hmC_deseqResults_list[[3]],x,"hmC"))
## b = lapply(Ks,function(x)pheat_var(rld_hmC[[3]],hmC_deseqResults_list[[3]],x,"hmC"))


