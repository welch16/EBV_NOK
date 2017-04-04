
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

deseqresults = lapply(deseqfiles,read_tsv)

## D.E. genes by cell
diffgenes = deseqresults %>%
    map( .f = function(x){
        out = x %>% dplyr::select(gene,dplyr::contains("cell:"))
        colnames(out) = gsub("cell:","",colnames(out))
        out %>% separate(gene , into = c("ensembl","gene"),sep = "\\_") %>%
            dplyr::rename(p.value = pvalue)
        })

names(diffgenes) = c("ben_CaFBS","ben_MC","scott_MC","scott_MC_woRep1")

## load DIP-seq count matrix
dipdr = "data/MeDIPseq_results"
dipmat = read_tsv(file.path(dipdr,"MeDIPseq_PromotersCounts_upstr500_downstr1000fraglen300.tsv"))
colnames(dipmat)[1] = "ensembl"


## load RSEM TPM matrix
dr = "data/metadata"
rsemfile = list.files(dr,full.names = TRUE,pattern = "TPM")
tpmmat = read_tsv(rsemfile)
tpmmat = tpmmat %>% separate(gene_id, into = c("ensembl","gene"),sep = "\\_")


## load CpG islands from UCSC genome browser
## session =  browserSession()
## genome(session) = "hg19"
## cpg = session[["CpG Islands"]]

## cpgmat = cpg %>% as.data.frame %>% as.tbl %>%
##     mutate(strand = NULL) %>% select(name,everything()) %>%
##     filter(!grepl("random",seqnames))
## write_tsv(cpgmat,"data/metadata/UCSC_genomeBrowser_CpGislands_hg19.tsv")

cpgmat = read_tsv("data/metadata/UCSC_genomeBrowser_CpGislands_hg19.tsv")

## obtain nr of aligned reads in each of the MeDIP-seq files
## bamfiles = list.files("data/BAM/hg19/bowtie_dip/",
##                       pattern = "sort",full.names = TRUE) %>%
##     .[grep("bai",.,invert= TRUE)]

## ## dipdepths = bamfiles %>% map(Rsamtools::countBam) %>% map_int(.f = function(.).$records)

## dipdepths = tibble(file = basename(bamfiles),depth = dipdepths)


## write_tsv(dipdepths,
##           "data/metadata/MeDIPseq_sequencingDepth.tsv")

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

## filter dipmat into genes,

diplist = diffgenes %>% map(left_join,dipmat,by = "ensembl")

sF = 1e9

normalize_dip <- function(dip,dipdepth,scaleFactor )
{
    base = dip[,seq_len(10)]

    counts = dip[,-seq_len(10)] %>% as.list
    nms = counts %>% names

    rpkm = counts %>% map2(dipdepths$depth,
                           .f = function(count,nreads) scaleFactor * count / nreads) %>%
        bind_cols

    bind_cols(base,rpkm)

}

## this contains RPKM because of the factor
diplist = diplist %>% map(normalize_dip,dipdepths,sF)

clean_cols <- function(geneList)
{

    X = geneList[,-seq_len(10)]
    
    toremove = X %>% as.list %>%
        map(.f = function(x)is.nan(x) | is.na(x)) %>%
        map(which) %>% do.call(c,.) %>% unique

    geneList[-toremove,]


}



divide_by_Input <- function(dip)
{

    base = dip[,seq_len(10)]

    input = dip %>% select(contains("Input"))
    dip = dip[,-seq_len(10)] %>% select(-contains("Input"))

    f = 1e-8
    akata_treat = dip %>%
        select(contains("akata")) %>% select(contains("CaFBS")) %>%
        mutate_each( funs( (f + .) /( f +  input[[1]])))
    akata_mono =  dip %>%
        select(contains("akata")) %>% select(contains("mono")) %>%
        mutate_each( funs( (f + .) / (f + input[[2]])))
    noks_treat =  dip %>%
        select(-contains("akata")) %>% select(contains("CaFBS")) %>%
        mutate_each( funs( (f + .) / (f + input[[3]])))
    noks_mono = dip %>%
        select(-contains("akata")) %>% select(contains("mono")) %>%
        mutate_each( funs( (f + .) / (f + input[[4]])))

    out = bind_cols(akata_treat,
              akata_mono,
              noks_treat,
              noks_mono) %>% mutate_all(funs(log10(.)))
    bind_cols(base,out)
}


diplist = diplist %>% map(clean_cols)
diplist = diplist %>% map(divide_by_Input)

heatmap_wrap <- function(geneList,methTr,scale)
{
    stopifnot(methTr %in% c("mC","hmC"))

    base = geneList[,seq_len(10)]
    X = geneList[,-seq_len(10)]
    
    if(methTr == "mC"){
        X = X %>% select(-contains("hmC"))
    }else{
        X = X  %>% select(contains("hmC"))
    }

    colnames(X) = gsub("MeDIPseq-","",colnames(X))
    
    superheat(X,
              pretty.order.rows = TRUE,
              pretty.order.cols = TRUE,
              left.label = "none",
              grid.vline.col = "white",
              bottom.label.text.angle = 80,
              bottom.label.text.size = 5,
              grid.hline = FALSE,
              scale = scale,
              heat.na.col = "white")
    
}

## plot MeDIP-seq heatmaps

figsdr = "figs/joinGeneExpression_Methylation"
    
png(file.path(figsdr,"methylHeatmaps","Heatmap_benCaFBS_mC.png"),
    width = 1200,height = 1400)
heatmap_wrap( diplist[[1]],"mC",TRUE)
dev.off()

png(file.path(figsdr,"methylHeatmaps","Heatmap_benMC_mC.png"),
    width = 1200,height = 1400)
heatmap_wrap( diplist[[2]],"mC",TRUE)
dev.off()

png(file.path(figsdr,"methylHeatmaps","Heatmap_scottMC_mC.png"),
    width = 1200,height = 1400)
heatmap_wrap( diplist[[3]],"mC",TRUE)
dev.off()

png(file.path(figsdr,"methylHeatmaps","Heatmap_scottMC_woRep1_mC.png"),
    width = 1200,height = 1400)
heatmap_wrap( diplist[[4]],"mC",TRUE)
dev.off()

png(file.path(figsdr,"methylHeatmaps","Heatmap_benCaFBS_hmC.png"),
    width = 1200,height = 1400)
heatmap_wrap( diplist[[1]],"hmC",TRUE)
dev.off()

png(file.path(figsdr,"methylHeatmaps","Heatmap_benMC_hmC.png"),
    width = 1200,height = 1400)
heatmap_wrap( diplist[[2]],"hmC",TRUE)
dev.off()

png(file.path(figsdr,"methylHeatmaps","Heatmap_scottMC_hmC.png"),
    width = 1200,height = 1400)
heatmap_wrap( diplist[[3]],"hmC",TRUE)
dev.off()

png(file.path(figsdr,"methylHeatmaps","Heatmap_scottMC_woRep1_hmC.png"),
    width = 1200,height = 1400)
heatmap_wrap( diplist[[4]],"hmC",TRUE)
dev.off()

t.test_wrap <- function(...)
{
    args = list(...)

    x = args[[1]] %>% t %>% .[,1]
    y = args[[2]] %>% t %>% .[,1]
    
    t.test(x,y,var.equal = TRUE) %>% broom::glance() %>% as.tbl
    
}   

## wilcox.test_wrap <- function(...)
## {
##     args = list(...)

##     x = args[[1]] %>% t %>% .[,1]
##     y = args[[2]] %>% t %>% .[,1]

##     wilcox.test(x,y) %>% broom::glance() %>% as.tbl %>%
##         mutate(estimate1 = median(x),
##                estimate2 = median(y))
     
## }


test_geneSet <- function(geneList,methTr,test = "t")
{
    stopifnot(methTr %in% c("mC","hmC"))

    dt = geneList[,-seq_len(10)]

    if(methTr == "mC"){
        dt = dt %>% select(-contains("hmC"))
    }else{
        dt = dt  %>% select(contains("hmC"))
    }
    
    X = dt %>% select(contains("akata")) 
    Y = dt %>% select(-contains("akata"))

    if(test == "t"){
        out = bpmapply(t.test_wrap,
                       X %>% split(geneList$ensembl),
                       Y %>% split(geneList$ensembl),SIMPLIFY = FALSE,BPPARAM = bp)

    }else{
        out = bpmapply(wilcox.test_wrap,
                       X %>% split(geneList$ensembl),
                     Y %>% split(geneList$ensembl),SIMPLIFY = FALSE,BPPARAM = bp)                          
    }

    out = out %>% bind_rows %>% mutate(ensembl = geneList$ensembl,
                                 gene = geneList$gene) %>%
        select(ensembl,gene,everything())


    out %>%
        dplyr::rename(EBV = estimate1,NOKS = estimate2)
      
}

diffmethMC = diplist %>% map(test_geneSet,"mC","t")
diffmethHMC = diplist %>% map(test_geneSet,"hmC","t")

## diffmethMC_wilcox = diplist %>% map(test_geneSet,"mC","wilcox")
## diffmethHMC_wilcox = diplist %>% map(test_geneSet,"hmC","wilcox")


## plot p.value histograms and adjust p.values

library(gridExtra)
library(grid)


pval_hist <- function(dts,titles)
{
    plots = map2(dts,titles,
                 .f = function(x,y){
                     x %>%
                         ggplot(aes(p.value))+
                         geom_histogram(bins = 20,
                                        boundary = 0,
                                        fill = "white",
                                        colour = "black")+
                         ggtitle(y)})
    arrangeGrob(plots[[1]],plots[[2]],plots[[3]],plots[[4]],
                 nrow = 2)                     
}

pp = pval_hist(diffmethMC,names(diffmethMC))
ggsave( filename = file.path(figsdr,"diffmethylPvalues","HistPval_mC.pdf"),
       plot = pp)

pp = pval_hist(diffmethHMC,names(diffmethHMC))

ggsave( filename = file.path(figsdr,"diffmethylPvalues","HistPval_hmC.pdf"),
       plot = pp)


diffmethMC = diffmethMC %>%
    map(.f = function(x)x %>%
                        mutate(padj = p.adjust(p.value,method = "BH")))

diffmethHMC = diffmethHMC %>%
    map(.f = function(x)x %>%
                        mutate(padj = p.adjust(p.value,method = "BH")))

save(diffmethMC,file = "data/MeDIPseq_results/T_diffMethyl_MC.RData")
save(diffmethHMC, file = "data/MeDIPseq_results/T_diffMethyl_HMC.RData")

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

deseq_input <- function(countData)
{

    deseqData = DESeqDataSetFromMatrix(countData = countData[[1]],
                                 colData = countData[[2]],
                                 design = ~ cell)

    deseqData
    


}

mC_deseqData_list = mC_countData_list %>% map(deseq_input)
hmC_deseqData_list = hmC_countData_list %>% map(deseq_input)



deseq_analysis <- function(deseqData)
{
    dds = DESeq(deseqData)

    res = results(dds,cooksCutoff = FALSE,contrast = c("cell","EBV","NOKS"),
                  tidy = TRUE,parallel = TRUE,BPPARAM = bp)

    res %>%  as.tbl %>%
        dplyr::rename(p.value = pvalue,ensembl = row)
     

}


mC_deseqResults_list = mC_deseqData_list %>% map(deseq_analysis)
hmC_deseqResults_list = hmC_deseqData_list %>% map(deseq_analysis)


pp = pval_hist(mC_deseqResults_list,names(mC_deseqResults_list))

ggsave( filename = file.path(figsdr,"diffmethylPvalues","HistPval_mC_DESeq2.pdf"),
       plot = pp)

pp = pval_hist(hmC_deseqResults_list,names(hmC_deseqResults_list))

ggsave( filename = file.path(figsdr,"diffmethylPvalues","HistPval_hmC_DESeq2.pdf"),
       plot = pp)


mC_deseqResults_list = mC_deseqResults_list %>% map(.f = function(x)
    x %>% left_join( tpmmat[,1:2],by = "ensembl") %>%
    select(ensembl,gene,everything()))

hmC_deseqResults_list = hmC_deseqResults_list %>% map(.f = function(x)
    x %>% left_join( tpmmat[,1:2],by = "ensembl") %>%
    select(ensembl,gene,everything()))

save(mC_deseqResults_list,file = "data/MeDIPseq_results/DESeq2_diffMethyl_MC.RData")
save(hmC_deseqResults_list,file = "data/MeDIPseq_results/DESeq2_diffMethyl_HMC.RData")


## overlap of diff. expressed and diff. methyl genes,

mC_jointAnalysis_list = diffgenes %>% map2(mC_deseqResults_list,
                                     left_join,
                                     by = c("ensembl","gene"),
                                     suffix = c(".exp",".methyl"))

hmC_jointAnalysis_list = diffgenes %>% map2(hmC_deseqResults_list,
                                     left_join,
                                     by = c("ensembl","gene"),
                                     suffix = c(".exp",".methyl"))




## save(mC_deseqResults_list,file = "data/MeDIPseq_results/DESeq2_diffMethyl_MC.RData")
## save(hmC_deseqResults_list,file = "data/MeDIPseq_results/DESeq2_diffMethyl_HMC.RData")




## ## overlap with CpG islands

## rd$cpg = ifelse(countOverlaps(rd,cpg) > 0, TRUE, FALSE)
## auxrd = rd %>% mcols %>% as.data.frame %>% as.tbl %>% dplyr::rename(ensembl = ensembl_gene_id)

## diffgenes = diffgenes %>% separate(gene,into = c("ensembl","name"),remove = FALSE,sep = "\\_")

## diffgenes = left_join(diffgenes,auxrd,by ="ensembl")

## ## only take no-treatment

## genefiles = genefiles[grep("noks",genefiles)]

## ## read count data only
## txdata = tximport(genefiles,type = "rsem" , importer = read_tsv)
## txdata = txdata[["counts"]]

## nms = txdata %>% rownames
## colnames(txdata) = gsub(".genes.results","",basename(genefiles))


## ## indexes
## akatamono = which(grepl("akata",colnames(txdata)) & grepl("no_treat",colnames(txdata)))
## noksmono = which(!grepl("akata",colnames(txdata)) & grepl("no_treat",colnames(txdata)))
## akataCaFBS = which(grepl("akata",colnames(txdata)) & grepl("CaFBS",colnames(txdata)))
## noksCaFBS = which(!grepl("akata",colnames(txdata)) & grepl("CaFBS",colnames(txdata)))
## akataMC = which(grepl("akata",colnames(txdata)) & grepl("methyl",colnames(txdata)))
## noksMC = which(!grepl("akata",colnames(txdata)) & grepl("methyl",colnames(txdata)))

## akataCaFBS = rowMeans(txdata[,c(akatamono,akataCaFBS)])
## akataMC = rowMeans(txdata[,c(akatamono,akataMC)])
## noksCaFBS = rowMeans(txdata[,c(noksmono,noksCaFBS)])
## noksMC = rowMeans(txdata[,c(noksmono,noksMC)])


## counts = tibble(genes = nms,
##                 akataCaFBS,
##                 akataMC,
##                 noksCaFBS,
##                 noksMC) %>%
##     separate(genes,into = c("ensembl","name"),remove = FALSE,sep = "\\_")

## cumulativeDist = counts %>%
##     mutate_if(is.numeric,cume_dist)



## diffgenes = diffgenes %>% left_join(cumulativeDist,by = "ensembl") %>%
##     mutate(EBV_dist = ifelse(model == "noneVsCaFBS",akataCaFBS,akataMC),
##            NOKS_dist = ifelse(model == "noneVsCaFBS",noksCaFBS,noksMC),
##            genes = NULL,name.y = NULL,
##            akataCaFBS = NULL,akataMC = NULL,noksCaFBS = NULL,noksMC = NULL) %>%
##     dplyr::rename(name = name.x)


## ## dip-seq reads in genes

## fraglen = 150

## dipdr = "data/BAM/hg19/bowtie_dip"
## dipfiles = list.files(dipdr,full.names = TRUE,pattern = "sort")
## dipfiles = dipfiles[grep("bai",dipfiles,invert = TRUE)]

## ## dipfiles = dipfiles[grep("mono",dipfiles,invert = TRUE)] 

## inputfiles = dipfiles[grepl("Input",dipfiles) & grepl("rep0",dipfiles)]
## dipfiles = dipfiles[grep("Input",dipfiles,invert = TRUE)]
## dipfiles = dipfiles[grep("hmC",dipfiles,invert = TRUE)]


## process_aligned <- function(files,fraglen){
##     reads = files %>% mclapply(readGAlignments,param = NULL)
##     reads = reads %>% mclapply(as,"GRanges")
##     names(reads) = gsub(".sort.bam","",basename(files))

##     reads = reads %>% mclapply(resize,fraglen)
##     reads}

## dip = process_aligned(dipfiles,fraglen)
## input = process_aligned(inputfiles,fraglen)

## dipcounts = dip %>% mclapply(function(x)countOverlaps(rd,x))
## inputcounts = input %>% mclapply(function(x)countOverlaps(rd,x))

## dipdepth = dip %>% sapply(length)
## inputdepth = input %>% sapply(length) 


## dipcounts = mapply(function(x,y)1e9 *  x /y ,dipcounts , dipdepth,SIMPLIFY = FALSE)
## inputcounts = mapply(function(x,y)1e9 *  x /y ,inputcounts , inputdepth,SIMPLIFY = FALSE)

## dipmat = dipcounts %>% as.data.frame %>% as.matrix
## inputmat = inputcounts %>% as.data.frame %>% as.matrix

## w1 = ifelse(dipdepth > 10e6,1,0)
## w2 = ifelse(inputdepth > 10e6,1,0)

## expt_counts = list()
## expt_counts[["EBV"]] = matrixStats::rowWeightedMeans(dipmat,w = ifelse(grepl("akata",names(w1)),w1,0))
## expt_counts[["NOKS"]] = matrixStats::rowWeightedMeans(dipmat,w = ifelse(!grepl("akata",names(w1)),w1,0))
## expt_counts[["EBV_Input"]] = matrixStats::rowWeightedMeans(inputmat,w = ifelse(grepl("akata",names(w2)),w2,0))
## expt_counts[["NOKS_Input"]] = matrixStats::rowWeightedMeans(inputmat,w = ifelse(!grepl("akata",names(w2)),w2,0))

## mcols(rd) = cbind(mcols(rd),DataFrame(width = width(rd)),DataFrame(expt_counts))

## auxrd = rd %>% mcols %>% as.data.frame %>% as.tbl %>% dplyr::rename(ensembl = ensembl_gene_id) %>%
##     mutate(cpg = NULL)



## diffgenes = diffgenes %>% left_join(auxrd,by = "ensembl")


## diffgenes = diffgenes %>% mutate(DIP_log2FC =log2( NOKS / EBV),
##                                  DIPInput_log2FC = log2(NOKS_Input / EBV_Input))






## geneExpression_count_plot <- function(DT,xvar,yvar,sc)
## {
##   require(viridis)
##   require(scales)
##   require(ggplot2)
##   pal = viridis(1e3, option = "D")
  
##   ggplot(DT,aes_string(xvar,yvar))+stat_binhex(bins = 100) +
##     scale_fill_gradientn(colours = pal,trans = "log10",
##                          labels = trans_format('log10',math_format(10^.x)),
##                          guide = guide_colourbar(title = NULL,
##                            barheight = unit(0.92,"npc"),
##                            barwidth = unit(0.01,"npc")))+
##     xlim(-sc,sc)+ylim(-sc,sc)
## }


## library(ggplot2)

## theme_set(theme_bw()+theme(legend.position = "top"))

## pdf(width = 10,height = 10)
## diffgenes %>% filter(!is.na(cpg) & !is.na(padj)) %>%
##     ggplot(aes(log2FC,DIPInput_log2FC,colour = padj <= 1e-3))+
##     geom_point(alpha = I(1/3))+ylim(-2,2)+xlim(-10,10)+
##     geom_smooth(method = "lm",se =FALSE,colour = "orange")+
##     facet_grid(model ~ cpg)+
##     scale_color_manual(values = c("navyblue","red"))
## diffgenes %>% filter(!is.na(cpg) & !is.na(padj)) %>%
##     ggplot(aes(log2FC,DIP_log2FC,colour = padj <= 1e-3))+
##     geom_point(alpha = I(1/3))+ylim(-2,2)+xlim(-10,10)+
##     geom_smooth(method = "lm",se =FALSE,colour = "orange")+
##     facet_grid(model ~ cpg)+
##     scale_color_manual(values = c("navyblue","red"))
## dev.off()


## write_tsv(diffgenes,path = "data/integration_RNAseq_MeDIPseq.tsv")

##  ## diffgenes %>% filter(padj <= 1e-7 & EBV_dist > .75 & NOKS_dist > .75 & cpg) %>% filter(log2FC < -1)  %>% filter(model == "noneVsMC") %>% arrange(log2FC) %>% dplyr::select(dplyr::contains("MeDIP")) %>% dplyr::select(-dplyr::contains("hmC")) %>% as.data.frame %>% head(5)


