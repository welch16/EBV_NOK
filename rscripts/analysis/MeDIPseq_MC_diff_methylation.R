
rm(list = ls())

library(GenomicAlignments)
library(MEDIPS)
library(MEDIPSData)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(BiocParallel)

## general parameters, to match with the used in promoter matrix


sizes = readr::read_tsv("/p/keles/SOFTWARE/hg19.chrom.sizes",col_names = FALSE)
chr = sizes$X1

mindepth = 10e6


indr = "data/BAM/hg19/bowtie_dip"
dipdepths = readr::read_tsv("data/metadata/MeDIPseq_sequencingDepth.tsv") %>%
    filter(depth > mindepth)

input_files = dipdepths %>% filter(grepl("rep0",file)) %>% {file.path(indr,.$file)}

mC_files = dipdepths %>% filter(!grepl("Input",file) & !grepl("hmC",file)) %>%
    {file.path(indr,.$file)}

hmC_files = dipdepths %>% filter(!grepl("Input",file) & grepl("hmC",file)) %>%
    {file.path(indr,.$file)}


uniq = 1e-3
extend = 300
shift = 0
ws = 200
BSgenome = "BSgenome.Hsapiens.UCSC.hg19"

bp = MulticoreParam(workers = 16)
    
message("Loading input files ...................")
input_sets = bplapply(input_files,
                      MEDIPS.createSet,BSgenome = BSgenome,
                      extend = extend, shift = shift, uniq = uniq,
                      window_size = ws, chr.select = chr,
                      BPPARAM = bp)

message("Loading treatment files ...................")
mC_sets = bplapply(mC_files,
                   MEDIPS.createSet,BSgenome = BSgenome,
                   extend = extend, shift = shift,uniq = uniq,
                   window_size = ws,chr.select =chr,
                   BPPARAM = bp)

## message("Loading treatment files ...................")
## hmC_sets = bplapply(hmC_files,
##                    MEDIPS.createSet,BSgenome = BSgenome,
##                    extend = extend, shift = shift,uniq = uniq,
##                    window_size = ws,chr.select =chr,
##                    BPPARAM = bp)

names(input_sets) = basename(input_files)
names(mC_sets) = basename(mC_files)
## names(hmC_sets) = basename(hmC_files)


message("Calculating coupling vector..................")
couplingVec = MEDIPS.couplingVector(pattern = "CG",refObj = input_sets[[1]])



## MEDIPS.plotCalibrationPlot(CSet = couplingVec, main = "Calibration Plot", MSet = DE_sets[[1]],
##                                  plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)

p.adj =1e-1
minRowSum = 10 
  
message("Fitting edgeR model........................")
mC_edgeR_results = MEDIPS.meth(MSet1 = mC_sets[grepl("akata",names(mC_sets))],
                                 MSet2 = mC_sets[!grepl("akata",names(mC_sets))],
                                 CSet = couplingVec,
                                 ISet1 = input_sets[grepl("akata",names(input_sets))],
                                 ISet2 = input_sets[!grepl("akata",names(input_sets))],
                                 p.adj = "BH", diff.method = "edgeR",
                                MeDIP = TRUE, CNV = FALSE, minRowSum = minRowSum)

save(mC_edgeR_results,file = "./mC_edgeR.RData")

## hmC_edgeR_results = MEDIPS.meth(MSet1 = hmC_sets[grepl("akata",names(hmC_sets))],
##                                  MSet2 = hmC_sets[!grepl("akata",names(hmC_sets))],
##                                  CSet = couplingVec,
##                                  ISet1 = input_sets[grepl("akata",names(input_sets))],
##                                  ISet2 = input_sets[!grepl("akata",names(input_sets))],
##                                  p.adj = "BH", diff.method = "edgeR",
##                                 MeDIP = TRUE, CNV = FALSE, minRowSum = minRowSum)

## save(hmC_edgeR_results,file = "hmC_edgeR_t.RData")


    
##     message("Finding differentially methylated regions ....................")
##     edgeR_Sig = MEDIPS.selectSig(results = edgeR_results, p.value = pval, adj = TRUE,
##                               ratio = NULL, bg.counts = NULL, CNV = FALSE)



## chr.select = chr[23]

## MEDIPS_diff_methylation <- function(chr.select,                                    
##                                     DE_files,
##                                     input_files,
##                                     BSgenome,
##                                     extend,
##                                     shift,
##                                     uniq,
##                                     ws,
##                                     pval = 0.1,
##                                     p.adj = "BH",
##                                     diff.meth = "edgeR",
##                                     minRowSum = 10)
## {


##     browser()

    
##     message("Loading input files ...................")
##     input_sets = bplapply(input_files,
##                         MEDIPS.createSet,BSgenome = BSgenome,
##                                extend = extend, shift = shift, uniq = uniq,
##                         window_size = ws, chr.select = chr.select,
##                         BPPARAM = bp)

##     message("Loading treatment files ...................")
##     DE_sets = bplapply(DE_files,
##                      MEDIPS.createSet,BSgenome = BSgenome,
##                      extend = extend, shift = shift,uniq = uniq,
##                      window_size = ws,chr.select =chr.select,BPPARAM = bp)

##     names(input_sets) = basename(input_files)
##     names(DE_sets) = basename(DE_files)

##     message("Calculating coupling vector..................")
##     couplingVec = MEDIPS.couplingVector(pattern = "CG",refObj = DE_sets[[3]])


##     MEDIPS.plotCalibrationPlot(CSet = couplingVec, main = "Calibration Plot", MSet = DE_sets[[1]],
##                                  plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
    
    
##     message("Fitting edgeR model........................")
##     edgeR_results = MEDIPS.meth(MSet1 = DE_sets[grepl("akata",names(DE_sets))],
##                                  MSet2 = DE_sets[!grepl("akata",names(DE_sets))],
##                                  CSet = couplingVec,
##                                  ISet1 = input_sets[grepl("akata",names(input_sets))],
##                                  ISet2 = input_sets[!grepl("akata",names(input_sets))],
##                                  p.adj = p.adj, diff.method = diff.meth,
##                                 MeDIP = TRUE, CNV = FALSE, minRowSum = minRowSum)


    
    
##     message("Finding differentially methylated regions ....................")
##     edgeR_Sig = MEDIPS.selectSig(results = edgeR_results, p.value = pval, adj = TRUE,
##                               ratio = NULL, bg.counts = NULL, CNV = FALSE)

## }


## aa = MEDIPS_diff_methylation(chr,
##                         mC_files,input_files,BSgenome,
##                         extend,shift,uniq,ws = 1e2,pval = 0.1,minRowSum = 10)




## ## mC_edgeR_Sig_gain = mC_edgeR_Sig[which(mC_edgeR_Sig[, grep("logFC", colnames(mC_edgeR_Sig))] > 0), ]
## ## mr.edgeR.s.gain.m = MEDIPS.mergeFrames(frames = mr.edgeR.s.gain, distance = 1)


## ## columns = names(mr.edgeR)[grep("counts", names(mr.edgeR))]
## ## rois = MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.gain.m,columns = columns, summarize = NULL)

## ## rois.s = MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.gain.m,columns = columns, summarize = "avg")

## ## anno.mart.gene = MEDIPS.getAnnotation(dataset = c("hsapiens_gene_ensembl"),annotation = c("GENE"), chr = "chr22")





## ## mC_sets = lapply(mC_files,
## ##                    MEDIPS.createSet,BSgenome = BSgenome,
## ##                                extend = extend, shift = shift, uniq = uniq,
## ##                                window_size = ws, chr.select = chr.select)

## ## hmC_sets = lapply(hmC_files,
## ##                    MEDIPS.createSet,BSgenome = BSgenome,
## ##                                extend = extend, shift = shift, uniq = uniq,
## ##                                window_size = ws, chr.select = chr.select)

## ## names(input_sets) = basename(input_files)
## ## names(mC_sets) = basename(mC_files)
## ## names(hmC_sets) = basename(mC_files)
