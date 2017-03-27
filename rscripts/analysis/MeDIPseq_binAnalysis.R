
rm(list = ls())

library(GenomicAlignments)
library(GenomicRanges)
library(tidyverse)

library(parallel)
options(mc.cores = 20)

devtools::load_all("~/Desktop/Docs/Code/ChIPUtils")

indr = "data/BAM/hg19/bowtie_dip"
files = list.files(indr,full.names = TRUE,pattern = "sort")

files = files[grep("bai",files,invert = TRUE)]
sizes = read_tsv("/p/keles/SOFTWARE/hg19.chrom.sizes",col_names = FALSE) %>%
    filter(X1 != "chrM")

binsize = 500
fraglen = 300
shift = 0


chr = "chr10"

reads_chr <- function(chr,files)
{
    sizes = sizes %>% filter(X1 == chr)
    gr = GRanges(sizes$X1,ranges = IRanges(start = 1, width = sizes$X2))
    reads = mclapply(files,readGAlignments,param = ScanBamParam(which = gr))
    reads = map(reads,as,"GRanges")
    reads = reads %>% map(resize,fraglen)
    reads
}    


reads = reads_chr(chr,files)


bins = create_bins(binsize,sizes %>% filter(X1 == chr))

depth  = reads %>% map_dbl(length)

counts = reads %>% map(.f = function(.)countOverlaps(bins,.))

rkm = map2(counts,depth,.f = function(x,y) 1e9 * x / y)

