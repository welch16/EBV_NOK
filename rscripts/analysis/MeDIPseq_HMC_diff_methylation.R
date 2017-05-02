
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
hmC_sets = bplapply(hmC_files,
                   MEDIPS.createSet,BSgenome = BSgenome,
                   extend = extend, shift = shift,uniq = uniq,
                   window_size = ws,chr.select =chr,
                   BPPARAM = bp)

names(input_sets) = basename(input_files)
names(hmC_sets) = basename(hmC_files)


message("Calculating coupling vector..................")
couplingVec = MEDIPS.couplingVector(pattern = "CG",refObj = input_sets[[1]])



## MEDIPS.plotCalibrationPlot(CSet = couplingVec, main = "Calibration Plot", MSet = DE_sets[[1]],
##                                  plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)

p.adj =1e-1
minRowSum = 10 

EBV_curves = lapply(hmC_sets[grepl("akata",names(hmC_sets))][-5],
           function(x)MEDIPS.calibrationCurve(MSet = x, CSet = couplingVec))


NOKS_curves = lapply(hmC_sets[!grepl("akata",names(hmC_sets))][-5],
           function(x)MEDIPS.calibrationCurve(MSet = x, CSet = couplingVec))


## plot diagnostics
## Removed the hmC_mono rep2 for NOKS and EBV as the calibration
## curve was not possible to be calculated




message("Fitting edgeR model........................")
hmC_edgeR_results = MEDIPS.meth(MSet1 = hmC_sets[grepl("akata",names(hmC_sets))][-5],
                                 MSet2 = hmC_sets[!grepl("akata",names(hmC_sets))][-5],
                                 CSet = couplingVec,
                                ISet1 = input_sets[grepl("akata",names(input_sets))],
                                ISet2 = input_sets[!grepl("akata",names(input_sets))],
                                 p.adj = "BH", diff.method = "edgeR",
                                MeDIP = TRUE, CNV = FALSE, minRowSum = minRowSum)

save(hmC_edgeR_results,file = "./data/MEDIPS/hmC_edgeR.RData")

