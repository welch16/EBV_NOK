rm(list = ls())

library(GenomicFeatures)
library(GenomicRanges)
library(MEDIPS)
library(tidyverse)
library(biomaRt)
library(viridis)
library(scales)
library(ChIPseeker)
library(clusterProfiler)
library(Vennerable)

## qi's packages
library(org.Hs.eg.db)
library(parallel)
options(mc.cores = 20)

## general parameters, to match with the used in promoter matrix
ups = 500
dns = 1e3

load(file = "./data/MEDIPS/hmC_edgeR.RData") # hmC_edgeR_results
load(file = "./data/MEDIPS/mC_edgeR.RData")  # mC_edgeR_results

write_tsv(mC_edgeR_results,path = "data/MEDIPS/mC_bins.tsv")
write_tsv(hmC_edgeR_results,path = "data/MEDIPS/hmC_bins.tsv")

pval = 0.1
mC_sig = MEDIPS.selectSig( mC_edgeR_results,p.value = pval)
hmC_sig = MEDIPS.selectSig( hmC_edgeR_results,p.value = pval)

mC_sig_merge = MEDIPS.mergeFrames(mC_sig,distance = 250)
hmC_sig_merge = MEDIPS.mergeFrames(hmC_sig,distance = 250)

convert2GRanges <- function(x)
{
    names(x)[1:3] = c("seqnames","start","end")

    x$start = x$start %>% as.numeric()
    x$end = x$end %>% as.numeric()
        
    gr = makeGRangesFromDataFrame(x)
    names(gr) = NULL

    mcols(gr) = x[,-(1:3)]
    
    gr
    
}   

diff_methyl_bins  =
    list(
        mC = mC_sig,
        hmC = hmC_sig
    ) %>%
    map(convert2GRanges)

diff_methyl_regions =
    list(
        mC = mC_sig_merge,
        hmC = hmC_sig_merge
    ) %>%
    map(convert2GRanges)

    
mergeSummaryStats <- function(regions,bins)
{

    ov = findOverlaps(regions,bins)
    ov_list = ov %>%
        as.data.frame() %>%
        as.tbl() %>%
        split(.$queryHits) %>%
        map(.f = function(x) x$subjectHits)

    bins_tbl = bins %>% as.data.frame() %>% as.tbl() %>%
        dplyr::select(-contains("counts")) %>%
        dplyr::select(-contains("sort.bam")) %>%
        mutate(
            regionID = queryHits(ov),
            width = NULL,
            strand = NULL,
            edgeR.adj.p.value = NULL
        ) %>%
        dplyr::select(regionID,
                      everything())

    bin_grp = bins_tbl %>%
        group_by(regionID)
    
    regions_tbl = bin_grp %>%
        summarize(
            seqnames = seqnames %>%  unique(),
            start = min(start),
            end = max(end),
            CF = sum(CF)
        )

    rpkm  = bin_grp %>%
        summarize_at( vars(matches("mean")), mean)

    logpart = bin_grp %>%
        summarize_at(vars(matches("log")),mean)


    pval = bin_grp %>%
        summarize(
            min.pval = min(edgeR.p.value),
            med.pval = median(edgeR.p.value)
        )
    

    regions_tbl %>%
        left_join(rpkm,by = "regionID") %>%
        left_join(logpart,by ="regionID") %>%
        left_join(pval, by = "regionID") %>%
        mutate(
            minlog10pval = -log10(min.pval),
            medlog10pval = -log10(med.pval)
            )




}    

diff_methyl_regions =
    diff_methyl_regions %>%
    map2(
        diff_methyl_bins,
        mergeSummaryStats)

write_tsv(diff_methyl_regions[["mC"]],"data/MEDIPS/Diff_mC_regions_merge.tsv")
write_tsv(diff_methyl_regions[["hmC"]],"data/MEDIPS/Diff_hmC_regions_merge.tsv")

diff_methyl_gr = diff_methyl_regions %>%
    map(.f = function(x)x %>% mutate(regionID = NULL)) %>%
    map(convert2GRanges)


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene



diff_methyl_annot = diff_methyl_gr %>%
    map(
        annotatePeak,
        tssRegion=c(-1000, 500),
        TxDb=txdb, annoDb="org.Hs.eg.db")



figsdr = "figs/diff_methylation/annotation"

methyl_split = diff_methyl_gr %>%
    map(.f = function(x){
        S4Vectors::split(x , ifelse(x$edgeR.logFC > 0 , "up","down")) %>%
            as.list()
    }) %>%
    unlist()

methyl_split_annot = methyl_split %>%
    map(
        annotatePeak,
        tssRegion = c(-dns,ups),
        TxDb = txdb, annoDb = "org.Hs.eg.db")


a = plotDistToTSS(methyl_split_annot,title = "Distance from Diff. Meth region to TSS")

pdf(file.path(figsdr,"mC_hmC_annotation_comparison.pdf"))
plotAnnoBar(diff_methyl_annot)
plotAnnoBar(methyl_split_annot)
plotDistToTSS(diff_methyl_annot,title = "Distance from Diff. Meth region to TSS")
plotDistToTSS(methyl_split_annot,title = "Distance from Diff. Meth region to TSS")
dev.off()

pdf(file.path(figsdr,"mC_annotation_summaries.pdf"))
plotAnnoPie(diff_methyl_annot[[1]])
upsetplot(diff_methyl_annot[[1]])
dev.off()

pdf(file.path(figsdr,"hmC_annotation_summaries.pdf"))
plotAnnoPie(diff_methyl_annot[[2]])
upsetplot(diff_methyl_annot[[2]])
dev.off()

pdf(file.path(figsdr,"mC_up_annotation_summaries.pdf"))
plotAnnoPie(methyl_split_annot[["mC.up"]])
upsetplot(methyl_split_annot[["mC.up"]])
dev.off()

pdf(file.path(figsdr,"mC_down_annotation_summaries.pdf"))
plotAnnoPie(methyl_split_annot[["mC.down"]])
upsetplot(methyl_split_annot[["mC.down"]])
dev.off()

pdf(file.path(figsdr,"hmC_up_annotation_summaries.pdf"))
plotAnnoPie(methyl_split_annot[["hmC.up"]])
upsetplot(methyl_split_annot[["hmC.up"]])
dev.off()

pdf(file.path(figsdr,"hmC_down_annotation_summaries.pdf"))
plotAnnoPie(methyl_split_annot[["hmC.down"]])
upsetplot(methyl_split_annot[["hmC.down"]])
dev.off()


gene_list = read_tsv("data/Diff.Genes/hg19/DESeq2_contrasts_genes_list.tsv")

diff_genes = list.files("data/Diff.Genes/hg19",full.names = TRUE,recursive = TRUE) %>%
    {.[grep("DESeq2_con",.)]} %>%
    {.[grep("list",.,invert = TRUE)]} %>% {.[-3]}

diff_genes = diff_genes %>% map(read_tsv) %>%
    map(.f = function(x) x %>% tidyr::separate(gene,into = c("ensembl_id","gene_id"), sep = "\\_"))

diff_genesP = diff_genes[[3]][!diff_genes[[3]]$ensembl_id %in% gene_list$ensembl_id,]

diff_genes = diff_genes %>%
    map(.f = function(x) right_join(x,gene_list,by = c("ensembl_id","gene_id")) %>%
                         dplyr::select(ensembl_id,gene_id,contains("cell:"))) %>%
    map(.f = function(x){
        names(x) = gsub("cell:","",names(x))
        x 
    }) %>%
    map(.f = function(x) x %>% mutate(reg = ifelse(log2FC > 0 , "Up","Down")))
   

names(diff_genes) = c("ben_CaFBS","ben_mC","scott_mC")


gene_overlap_tab <- function(meth_annot,gene_set,pval = 5e-2,distn = 10e3)
{

    my_genes = gene_set %>% dplyr::filter(padj <= pval) %>% {.$ensembl_id}
    w = meth_annot %>%
        map(.f = function(x)
            which(
                x@anno$ENSEMBL %in% my_genes &
                abs(x@anno$distanceToTSS) <= distn)
            )

    dets = meth_annot %>% map(.f = function(x) x@detailGenomicAnnotation)
    dets = dets %>% map2(w, .f = function(x,y)x[y,]) %>% map(as.tbl)

    nms = names(dets[[1]])
    dets %>% map(colSums) %>%
        map2(names(dets),
             .f = function(x,y)
                 tibble(region = nms,
                        counts = x,
                        meth = ifelse(grepl("down",y),"Hypo-methy","Hyper-methyl"))) %>%
        bind_rows() %>%
        mutate(pval_th = pval,
               distn)
          
}

theme_set(theme_bw())

dist = 1e3
pval=0.01
pdf(file.path(figsdr,paste0("Barplot_diff_exp_geneSets_dist",dist,"_pval",pval,".pdf")),height =5)
gene_overlap_tab(methyl_split_annot[1:2],
                 diff_genes[[1]],
                 pval,dist) %>%
    ggplot(aes(region,counts,fill = meth))+
        geom_bar(stat = "identity",position = "dodge",colour = "black")+
        scale_fill_brewer(palette = "Pastel1",name = "") +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 60,hjust = 1),
              axis.title.x = element_blank())+
    ggtitle("Ben_CaFBS")+ylab("Number of genes")
gene_overlap_tab(methyl_split_annot[1:2],
                 diff_genes[[2]],
                 pval,dist) %>%
    ggplot(aes(region,counts,fill = meth))+
        geom_bar(stat = "identity",position = "dodge",colour = "black")+
        scale_fill_brewer(palette = "Pastel1",name = "") +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 60,hjust = 1),
              axis.title.x = element_blank())+
    ggtitle("Ben_MC")+ylab("Number of genes")
gene_overlap_tab(methyl_split_annot[1:2],
                 diff_genes[[3]],
                 pval,dist) %>%
    ggplot(aes(region,counts,fill = meth))+
        geom_bar(stat = "identity",position = "dodge",colour = "black")+
        scale_fill_brewer(palette = "Pastel1",name = "") +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 60,hjust = 1),
              axis.title.x = element_blank())+
    ggtitle("Scott_MC")+ylab("Number of genes")
dev.off()



pdf(file.path(figsdr,paste0("Barplot_UpReg_exp_geneSets_dist",dist,"_pval",pval,".pdf")),height =5)
gene_overlap_tab(methyl_split_annot[1:2],
                 diff_genes[[1]] %>% dplyr::filter(log2FC > 0),
                 pval,dist) %>%
    ggplot(aes(region,counts,fill = meth))+
        geom_bar(stat = "identity",position = "dodge",colour = "black")+
        scale_fill_brewer(palette = "Pastel1",name = "") +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 60,hjust = 1),
              axis.title.x = element_blank())+
    ggtitle("Ben_CaFBS")+ylab("Number of genes")
gene_overlap_tab(methyl_split_annot[1:2],
                 diff_genes[[2]] %>% dplyr::filter(log2FC > 0),
                 pval,dist) %>%
    ggplot(aes(region,counts,fill = meth))+
        geom_bar(stat = "identity",position = "dodge",colour = "black")+
        scale_fill_brewer(palette = "Pastel1",name = "") +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 60,hjust = 1),
              axis.title.x = element_blank())+
    ggtitle("Ben_MC")+ylab("Number of genes")
gene_overlap_tab(methyl_split_annot[1:2],
                 diff_genes[[3]] %>% dplyr::filter(log2FC > 0),
                 pval,dist) %>%
    ggplot(aes(region,counts,fill = meth))+
        geom_bar(stat = "identity",position = "dodge",colour = "black")+
        scale_fill_brewer(palette = "Pastel1",name = "") +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 60,hjust = 1),
              axis.title.x = element_blank())+
    ggtitle("Scott_MC")+ylab("Number of genes")
dev.off()


pdf(file.path(figsdr,paste0("Barplot_DownReg_exp_geneSets_dist",dist,"_pval",pval,".pdf")),height =5)
gene_overlap_tab(methyl_split_annot[1:2],
                 diff_genes[[1]] %>% dplyr::filter(log2FC < 0),
                 pval,dist) %>%
    ggplot(aes(region,counts,fill = meth))+
        geom_bar(stat = "identity",position = "dodge",colour = "black")+
        scale_fill_brewer(palette = "Pastel1",name = "") +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 60,hjust = 1),
              axis.title.x = element_blank())+
    ggtitle("Ben_CaFBS")+ylab("Number of genes")
gene_overlap_tab(methyl_split_annot[1:2],
                 diff_genes[[2]] %>% dplyr::filter(log2FC < 0),
                 pval,dist) %>%
    ggplot(aes(region,counts,fill = meth))+
        geom_bar(stat = "identity",position = "dodge",colour = "black")+
        scale_fill_brewer(palette = "Pastel1",name = "") +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 60,hjust = 1),
              axis.title.x = element_blank())+
    ggtitle("Ben_MC")+ylab("Number of genes")
gene_overlap_tab(methyl_split_annot[1:2],
                 diff_genes[[3]] %>% dplyr::filter(log2FC < 0),
                 pval,dist) %>%
    ggplot(aes(region,counts,fill = meth))+
        geom_bar(stat = "identity",position = "dodge",colour = "black")+
        scale_fill_brewer(palette = "Pastel1",name = "") +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 60,hjust = 1),
              axis.title.x = element_blank())+
    ggtitle("Scott_MC")+ylab("Number of genes")
dev.off()

pdf(file.path(figsdr,paste0("Venn_diagram_diff_genes_pval",pval,".pdf")))
vennplot(
    diff_genes %>% map(.f = function(x) x %>%
                                        dplyr::filter(padj <= pval) %>%
                                        {.$ensembl_id}),by = "Vennerable")
dev.off()

pdf(file.path(figsdr,paste0("Venn_diagram_UpReg_diff_genes_pval",pval,".pdf")))
vennplot(
    diff_genes %>% map(.f = function(x) x %>%
                                        dplyr::filter(padj <= pval & log2FC > 0) %>%
                                        {.$ensembl_id}),by = "Vennerable")
dev.off()

pdf(file.path(figsdr,paste0("Venn_diagram_DownReg_diff_genes_pval",pval,".pdf")))
vennplot(
    diff_genes %>% map(.f = function(x) x %>%
                                        dplyr::filter(padj <= pval & log2FC < 0) %>%
                                        {.$ensembl_id}),by = "Vennerable")
dev.off()

## Venn diagram's Diff Expressed genes, Diff Methyl Genes

pdf(file.path(figsdr,paste0("VennDiagram_DiffExpBenCaFBS_DiffMethMC_genes_dist",dist,"_pval",pval,".pdf")))
genes = list(
    DiffMeth = diff_methyl_annot[["mC"]]@anno %>%
        as.data.frame() %>%
        as.tbl() %>%
        dplyr::filter(distanceToTSS <= dist) %>%
        dplyr::select(ENSEMBL) %>%
        {.[[1]]} %>% {.[!is.na(.)]},
    DiffExp = diff_genes[["ben_CaFBS"]] %>%
        dplyr::filter(padj <= pval) %>% {.$ensembl_id}
    )
vennplot( genes,by = "Vennerable")
dev.off()

pdf(file.path(figsdr,paste0("VennDiagram_DiffExpBenMC_DiffMethMC_genes_dist",dist,"_pval",pval,".pdf")))
genes = list(
    DiffMeth = diff_methyl_annot[["mC"]]@anno %>%
        as.data.frame() %>%
        as.tbl() %>%
        dplyr::filter(distanceToTSS <= dist) %>%
        dplyr::select(ENSEMBL) %>%
        {.[[1]]} %>% {.[!is.na(.)]},
    DiffExp = diff_genes[["ben_mC"]] %>%
        dplyr::filter(padj <= pval) %>% {.$ensembl_id}
    )
vennplot( genes,by = "Vennerable")
dev.off()

pdf(file.path(figsdr,paste0("VennDiagram_DiffExpScottMC_DiffMethMC_genes_dist",dist,"_pval",pval,".pdf")))
genes = list(
    DiffMeth = diff_methyl_annot[["mC"]]@anno %>%
        as.data.frame() %>%
        as.tbl() %>%
        dplyr::filter(distanceToTSS <= dist) %>%
        dplyr::select(ENSEMBL) %>%
        {.[[1]]} %>% {.[!is.na(.)]},
    DiffExp = diff_genes[["scott_mC"]] %>%
        dplyr::filter(padj <= pval) %>% {.$ensembl_id}
    )
vennplot( genes,by = "Vennerable")
dev.off()


