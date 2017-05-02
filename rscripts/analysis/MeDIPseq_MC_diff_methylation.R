
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

names(input_sets) = basename(input_files)
names(mC_sets) = basename(mC_files)


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

save(mC_edgeR_results,file = ".data/MEDIPS/mC_edgeR.RData")

