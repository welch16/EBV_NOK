
rm(list = ls())

library(tidyverse)
library(rtracklayer)


library(dplyr)
library(purrr)
library(tidyr)
library(magrittr)
library(readr)

library(tximport)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(IRanges)
library(biomaRt)

library(parallel)
options(mc.cores = 12)

genedr = "data/RSEM/hg19"
genefiles = list.files(genedr,full.names = TRUE,pattern = "genes.results")
    
## CpG islands
session =  browserSession()
genome(session) = "hg19"
cpg = session[["CpG Islands"]]

## Gene name - positions
mart <- useMart('ensembl', 'hsapiens_gene_ensembl')
attr <- c("chromosome_name", "start_position", "end_position",
          "ensembl_gene_id")
bm <- getBM(attr, "chromosome_name", c(as.character(1:22),"X","Y"), mart)
rd <- with(bm,
           RangedData(IRanges(start_position, end_position),
                      space=chromosome_name, ensembl_gene_id))

rd = as(rd,"GRanges")
seqlevels(rd) = paste0("chr",seqlevels(rd))

rd = promoters(rd,upstream = 1e3,downstream = 500)

## DESeq2 results

deseqdr = "data/Diff.Genes/hg19/DESeq2_contrasts"
deseqfiles = list.files(deseqdr,full.names = TRUE,pattern = "ben")

results = lapply(deseqfiles,read_tsv)


## D.E. genes by cell
diffgenes = results %>%
    lapply(
        function(x){
            out = x %>% dplyr::select(gene,dplyr::contains("cell:"))
            colnames(out) = gsub("cell:","",colnames(out))
            out
        })


diffgenes = mapply(function(x,y)x %>% mutate(model = y),diffgenes,c("noneVsCaFBS","noneVsMC"),SIMPLIFY = FALSE) %>%
    bind_rows

## overlap with CpG islands

rd$cpg = ifelse(countOverlaps(rd,cpg) > 0, TRUE, FALSE)
auxrd = rd %>% mcols %>% as.data.frame %>% as.tbl %>% dplyr::rename(ensembl = ensembl_gene_id)

diffgenes = diffgenes %>% separate(gene,into = c("ensembl","name"),remove = FALSE,sep = "\\_")

diffgenes = left_join(diffgenes,auxrd,by ="ensembl")

## only take no-treatment

genefiles = genefiles[grep("noks",genefiles)]

## read count data only
txdata = tximport(genefiles,type = "rsem" , importer = read_tsv)
txdata = txdata[["counts"]]

nms = txdata %>% rownames
colnames(txdata) = gsub(".genes.results","",basename(genefiles))


## indexes
akatamono = which(grepl("akata",colnames(txdata)) & grepl("no_treat",colnames(txdata)))
noksmono = which(!grepl("akata",colnames(txdata)) & grepl("no_treat",colnames(txdata)))
akataCaFBS = which(grepl("akata",colnames(txdata)) & grepl("CaFBS",colnames(txdata)))
noksCaFBS = which(!grepl("akata",colnames(txdata)) & grepl("CaFBS",colnames(txdata)))
akataMC = which(grepl("akata",colnames(txdata)) & grepl("methyl",colnames(txdata)))
noksMC = which(!grepl("akata",colnames(txdata)) & grepl("methyl",colnames(txdata)))

akataCaFBS = rowMeans(txdata[,c(akatamono,akataCaFBS)])
akataMC = rowMeans(txdata[,c(akatamono,akataMC)])
noksCaFBS = rowMeans(txdata[,c(noksmono,noksCaFBS)])
noksMC = rowMeans(txdata[,c(noksmono,noksMC)])


counts = tibble(genes = nms,
                akataCaFBS,
                akataMC,
                noksCaFBS,
                noksMC) %>%
    separate(genes,into = c("ensembl","name"),remove = FALSE,sep = "\\_")

cumulativeDist = counts %>%
    mutate_if(is.numeric,cume_dist)



diffgenes = diffgenes %>% left_join(cumulativeDist,by = "ensembl") %>%
    mutate(EBV_dist = ifelse(model == "noneVsCaFBS",akataCaFBS,akataMC),
           NOKS_dist = ifelse(model == "noneVsCaFBS",noksCaFBS,noksMC),
           genes = NULL,name.y = NULL,
           akataCaFBS = NULL,akataMC = NULL,noksCaFBS = NULL,noksMC = NULL) %>%
    dplyr::rename(name = name.x)


## dip-seq reads in genes

fraglen = 150

dipdr = "data/BAM/hg19/bowtie_dip"
dipfiles = list.files(dipdr,full.names = TRUE,pattern = "sort")
dipfiles = dipfiles[grep("bai",dipfiles,invert = TRUE)]

## dipfiles = dipfiles[grep("mono",dipfiles,invert = TRUE)] 

inputfiles = dipfiles[grepl("Input",dipfiles) & grepl("rep0",dipfiles)]
dipfiles = dipfiles[grep("Input",dipfiles,invert = TRUE)]
dipfiles = dipfiles[grep("hmC",dipfiles,invert = TRUE)]


process_aligned <- function(files,fraglen){
    reads = files %>% mclapply(readGAlignments,param = NULL)
    reads = reads %>% mclapply(as,"GRanges")
    names(reads) = gsub(".sort.bam","",basename(files))

    reads = reads %>% mclapply(resize,fraglen)
    reads}

dip = process_aligned(dipfiles,fraglen)
input = process_aligned(inputfiles,fraglen)

dipcounts = dip %>% mclapply(function(x)countOverlaps(rd,x))
inputcounts = input %>% mclapply(function(x)countOverlaps(rd,x))

dipdepth = dip %>% sapply(length)
inputdepth = input %>% sapply(length) 


dipcounts = mapply(function(x,y)1e9 *  x /y ,dipcounts , dipdepth,SIMPLIFY = FALSE)
inputcounts = mapply(function(x,y)1e9 *  x /y ,inputcounts , inputdepth,SIMPLIFY = FALSE)

dipmat = dipcounts %>% as.data.frame %>% as.matrix
inputmat = inputcounts %>% as.data.frame %>% as.matrix

w1 = ifelse(dipdepth > 10e6,1,0)
w2 = ifelse(inputdepth > 10e6,1,0)

expt_counts = list()
expt_counts[["EBV"]] = matrixStats::rowWeightedMeans(dipmat,w = ifelse(grepl("akata",names(w1)),w1,0))
expt_counts[["NOKS"]] = matrixStats::rowWeightedMeans(dipmat,w = ifelse(!grepl("akata",names(w1)),w1,0))
expt_counts[["EBV_Input"]] = matrixStats::rowWeightedMeans(inputmat,w = ifelse(grepl("akata",names(w2)),w2,0))
expt_counts[["NOKS_Input"]] = matrixStats::rowWeightedMeans(inputmat,w = ifelse(!grepl("akata",names(w2)),w2,0))

mcols(rd) = cbind(mcols(rd),DataFrame(width = width(rd)),DataFrame(expt_counts))

auxrd = rd %>% mcols %>% as.data.frame %>% as.tbl %>% dplyr::rename(ensembl = ensembl_gene_id) %>%
    mutate(cpg = NULL)



diffgenes = diffgenes %>% left_join(auxrd,by = "ensembl")


diffgenes = diffgenes %>% mutate(DIP_log2FC =log2( NOKS / EBV),
                                 DIPInput_log2FC = log2(NOKS_Input / EBV_Input))






geneExpression_count_plot <- function(DT,xvar,yvar,sc)
{
  require(viridis)
  require(scales)
  require(ggplot2)
  pal = viridis(1e3, option = "D")
  
  ggplot(DT,aes_string(xvar,yvar))+stat_binhex(bins = 100) +
    scale_fill_gradientn(colours = pal,trans = "log10",
                         labels = trans_format('log10',math_format(10^.x)),
                         guide = guide_colourbar(title = NULL,
                           barheight = unit(0.92,"npc"),
                           barwidth = unit(0.01,"npc")))+
    xlim(-sc,sc)+ylim(-sc,sc)
}


library(ggplot2)

theme_set(theme_bw()+theme(legend.position = "top"))

pdf(width = 10,height = 10)
diffgenes %>% filter(!is.na(cpg) & !is.na(padj)) %>%
    ggplot(aes(log2FC,DIPInput_log2FC,colour = padj <= 1e-3))+
    geom_point(alpha = I(1/3))+ylim(-2,2)+xlim(-10,10)+
    geom_smooth(method = "lm",se =FALSE,colour = "orange")+
    facet_grid(model ~ cpg)+
    scale_color_manual(values = c("navyblue","red"))
diffgenes %>% filter(!is.na(cpg) & !is.na(padj)) %>%
    ggplot(aes(log2FC,DIP_log2FC,colour = padj <= 1e-3))+
    geom_point(alpha = I(1/3))+ylim(-2,2)+xlim(-10,10)+
    geom_smooth(method = "lm",se =FALSE,colour = "orange")+
    facet_grid(model ~ cpg)+
    scale_color_manual(values = c("navyblue","red"))
dev.off()


write_tsv(diffgenes,path = "data/integration_RNAseq_MeDIPseq.tsv")

 ## diffgenes %>% filter(padj <= 1e-7 & EBV_dist > .75 & NOKS_dist > .75 & cpg) %>% filter(log2FC < -1)  %>% filter(model == "noneVsMC") %>% arrange(log2FC) %>% dplyr::select(dplyr::contains("MeDIP")) %>% dplyr::select(-dplyr::contains("hmC")) %>% as.data.frame %>% head(5)


