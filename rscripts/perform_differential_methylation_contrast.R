
rm(list = ls())

library(ggplot2)
library(superheat)
library(tidyverse)
library(rtracklayer)
library(BiocParallel)

bp = MulticoreParam(workers = 10,bpprogressbar = TRUE)

## load DeSeq2 results
deseqdr = "data/Diff.Genes/hg19/DESeq2_contrasts"
deseqfiles = list.files(deseqdr,full.names = TRUE)
deseqfiles = c(deseqfiles[grepl("ben",deseqfiles)],
               deseqfiles[grepl("wo",deseqfiles)])


deseqresults = lapply(deseqfiles,read_tsv)


## D.E. genes by cell
diffgenes = deseqresults %>%
    map( .f = function(x){
        out = x %>% dplyr::select(gene,dplyr::contains("cell:"))
        colnames(out) = gsub("cell:","",colnames(out))
        out %>% separate(gene , into = c("ensembl","gene"),sep = "\\_") %>%
            dplyr::rename(p.value = pvalue)
        })

names(diffgenes) = c("ben_CaFBS","ben_MC","scott_MC")

## load DIP-seq count matrix
dipdr = "data/MeDIPseq_results"
dipmat = read_tsv(file.path(dipdr,"MeDIPseq_PromotersCounts_upstr500_downstr1000fraglen300.tsv"))
colnames(dipmat)[1] = "ensembl"


## load RSEM TPM matrix
dr = "data/metadata"
rsemfile = list.files(dr,full.names = TRUE,pattern = "TPM")
tpmmat = read_tsv(rsemfile)
tpmmat = tpmmat %>% separate(gene_id, into = c("ensembl","gene"),sep = "\\_")

## load CpG islands
cpgmat = read_tsv("data/metadata/UCSC_genomeBrowser_CpGislands_hg19.tsv")




## We pooled the Input samples, therefore we are removing the base replicates that we used
dipdepths = read_tsv("data/metadata/MeDIPseq_sequencingDepth.tsv")

dipdepths = dipdepths %>%
    filter(!grepl("Input-rep1",file),
           !grepl("Input-rep2",file),
           !grepl("Input-rep3",file))

dipmat = dipmat %>% select(-contains("Input-rep1"),
                           -contains("Input-rep2"),
                           -contains("Input-rep3"))

## filter samples with very low sequencing depth. We considered 10M
## as in the QC analysis, the samples with less than 10M aligned reads
## seemed to be lower quality than the rest

minDepth = 10e6

w = which(dipdepths$depth > minDepth) %>% {dipdepths$file[.]} %>% {gsub(".sort.bam","",.)}

dipdepths = dipdepths %>% filter(depth > minDepth) %>%
    mutate(file = gsub(".sort.bam","",file))


dipmat = bind_cols(dipmat %>% select(-contains("MeDIPseq")),
                   dipmat[,w])

diplist = diffgenes %>% map(left_join,dipmat,by = "ensembl")

clean_cols <- function(geneList)
{

    X = geneList[,-seq_len(10)]
    
    toremove = X %>% as.list %>%
        map(.f = function(x)is.nan(x) | is.na(x)) %>%
        map(which) %>% do.call(c,.) %>% unique

    geneList[-toremove,]


}

diplist = diplist %>% map(clean_cols)


## this test is underpowered, hence we are going to try to
## test differential methylation with Deseq2


count_mats = diffgenes %>% map(.f = function(x){
    dipmat %>% dplyr::filter(ensembl %in% x$ensembl) %>%
        select(-contains("rep0"))})

transform2matrix <- function(mat_tbl,which )
{

    mat = mat_tbl %>% select(ensembl,contains("MeDIP"))

    ##

    toremove = mat %>% as.list %>%
        map(.f = function(x)is.nan(x) | is.na(x)) %>%
        map(which) %>% do.call(c,.) %>% unique

    if(length(toremove) > 0){
        mat = mat[-toremove,]
    }


    
    empty = mat[,-1] %>% rowSums
    empty = which(empty == 0)

    if(length(empty) > 0){
        mat = mat[-empty,]
    }


    if(which == "mC"){
        mat = mat %>% select(-contains("hmC"),ensembl)             
    }else{
        mat = mat %>% select(contains("hmC"),ensembl)        
    }

    mat = mat %>% select(ensembl,everything())

    nms = mat$ensembl
    mat = mat[,-1] %>% as.matrix
    rownames(mat) = nms

    out = list()

    out[["count_mat"]] = mat

    cols = colnames(mat)
    colData = data.frame(cell = ifelse(grepl("akata",cols),"EBV","NOKS"))
    rownames(colData) = cols

    out[["col_data"]] = colData
    out
    
}   

mC_countData_list = count_mats %>% map(transform2matrix,"mC")
hmC_countData_list = count_mats %>% map(transform2matrix,"hmC")


library(DESeq2)

deseqData <- function(countData)
{

    deseqData = DESeqDataSetFromMatrix(countData = countData[[1]],
                                 colData = countData[[2]],
                                 design = ~ cell)

    deseqData
    


}

mC_deseqData_list = mC_countData_list %>% map(deseqData)
hmC_deseqData_list = hmC_countData_list %>% map(deseqData)

mC_deseqModel_list = mC_deseqData_list %>% map(DESeq)
hmC_deseqModel_list = hmC_deseqData_list %>% map(DESeq)




deseq_analysis <- function(dds)
{
    res = results(dds,cooksCutoff = FALSE,contrast = c("cell","EBV","NOKS"),
                  tidy = TRUE,parallel = TRUE,BPPARAM = bp)

    res %>%  as.tbl %>%
        dplyr::rename(p.value = pvalue,ensembl = row)
     

}


mC_deseqResults_list = mC_deseqModel_list %>% map(deseq_analysis)
hmC_deseqResults_list = hmC_deseqModel_list %>% map(deseq_analysis)


mC_deseqResults_list = mC_deseqResults_list %>% map(.f = function(x)
    x %>% left_join( tpmmat[,1:2],by = "ensembl") %>%
    select(ensembl,gene,everything()))

hmC_deseqResults_list = hmC_deseqResults_list %>% map(.f = function(x)
    x %>% left_join( tpmmat[,1:2],by = "ensembl") %>%
    select(ensembl,gene,everything()))


save(mC_deseqResults_list,file = "data/MeDIPseq_results/DESeq2_diffMethyl_MC.RData")
save(hmC_deseqResults_list,file = "data/MeDIPseq_results/DESeq2_diffMethyl_HMC.RData")

outdr = "data/MeDIPseq_results"

mC_deseqResults_list %>% map2(names(mC_deseqResults_list),
                              .f = function(x,y){
                                  write_tsv(x,
                                            path = file.path(outdr,
                                                             paste0(y,"_diffMethyl_mC.tsv")))})
 
hmC_deseqResults_list %>% map2(names(hmC_deseqResults_list),
                              .f = function(x,y){
                                  write_tsv(x,
                                            path = file.path(outdr,
                                                             paste0(y,"_diffMethyl_hmC.tsv")))})



rld_mC = mC_deseqModel_list %>% map(rlog,blind = FALSE)
rld_hmC = hmC_deseqModel_list %>% map(rlog,blind = FALSE)

theme_set(theme_bw())


rowVars <- function (x,na.rm = TRUE)
{
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}


PCA_plot <- function(rl,tit,ntop = 500,x = 1,y = 2,tab)
{

    rv = rowVars(assay(rl))

    select = order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]

    pca = prcomp(t(assay(rl)[select, ]))
    percentVar = pca$sdev^2/sum(pca$sdev^2)

    intgroup.df =  as.data.frame(colData(rl)[, ,drop = FALSE]) %>% as.tbl

    dt = intgroup.df %>% mutate(file = tab[[1]],
                                x = pca$x[,x], y = pca$x[,y]) %>%
        mutate(cond = ifelse(grepl("mono",file),"No-Tr","CaFBS")) %>% split(.$cell)  %>%
        map(.f = function(z)
            z %>% mutate(Rep = paste("Rep",seq_len(nrow(z)),sep = "-"))) %>%
        bind_rows  

    dt %>%
        ggplot(aes(x,y,colour = cell))+
        geom_point(size = 2,shape = 1)+coord_fixed()+
        xlab(paste0("PC",x,": ",round(percentVar[x] * 100), "% variance"))+
        ylab(paste0("PC",y,": ",round(percentVar[y] * 100), "% variance"))+
        scale_color_brewer(palette = "Set1",name = "Cell")+        
        geom_label_repel(aes(x,y,label = file),size = 2,
                         box.padding = unit(.3,"lines"),
                         point.padding = unit(.7,"lines"),
                         show.legend = FALSE)+ ggtitle(tit)+
        theme(legend.position = "top")

}


library(ggrepel)

mC_PCA_plots = rld_mC %>%
    map2(names(rld_mC),PCA_plot,
         tab = dipdepths %>% filter(!grepl("Input",file)) %>%
             filter(!grepl("hmC",file))
         )

         
hmC_PCA_plots = rld_hmC %>%
    map2(names(rld_hmC),PCA_plot,
         tab = dipdepths %>% filter(!grepl("Input",file)) %>%
             filter(grepl("hmC",file))
         )



figsdr = "figs/diff_methylation/DESeq2_results"


pdf(file.path(figsdr,"DESeq2_PCA_mC.pdf"))
mC_PCA_plots %>% map(print)
dev.off()



pdf(file.path(figsdr,"DESeq2_PCA_hmC.pdf"))
hmC_PCA_plots %>% map(print)
dev.off()


library(pheatmap)

rowVars <- function (x,na.rm = TRUE)
{
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}



pheat_pval <- function(rld,my_results,K,suffix)
{
    mat = assay(rld)[,]
    mat = mat - rowMeans(mat)
    rownames(mat) = my_results$gene
    r = viridis::viridis(1e3,option = "D")

    coldata = colData(rld)
    coldata$sizeFactor = NULL

    coldata = coldata %>% as.data.frame
    
    topgenes = head(order(my_results$p.value),K)
    pheatmap(mat[topgenes,],color = r,
             annotation_col = coldata,
             show_rownames = K <= 75,
             filename = file.path(figsdr,paste0("Heatmap_bottom",K,"genes_pval_",suffix,".pdf")),
             height = 9,width = 6)
    
}    

pheat_var <- function(rld,my_results,K,suffix)
{
    mat = assay(rld)[,]
    mat = mat - rowMeans(mat)
    rownames(mat) = my_results$gene
    r = viridis::viridis(1e3,option = "D")

    coldata = colData(rld)
    coldata$sizeFactor = NULL

    coldata = coldata %>% as.data.frame

    

    topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)

    pheatmap(mat[topgenes,],annotation_col = coldata,color = r,
             show_rownames = K <= 75,
             filename = file.path(figsdr,paste0("Heatmap_top",K,"genes_var_",suffix,".pdf")),
             height = 9,width = 6)    
}    


Ks = c(20,50,100,5e2,1e3,5e3,1e4)

a = lapply(Ks,function(x)pheat_pval(rld_mC[[3]],mC_deseqResults_list[[3]],x,"mC"))
b = lapply(Ks,function(x)pheat_var(rld_mC[[3]],mC_deseqResults_list[[3]],x,"mC"))

a = lapply(Ks,function(x)pheat_pval(rld_hmC[[3]],hmC_deseqResults_list[[3]],x,"hmC"))
b = lapply(Ks,function(x)pheat_var(rld_hmC[[3]],hmC_deseqResults_list[[3]],x,"hmC"))


