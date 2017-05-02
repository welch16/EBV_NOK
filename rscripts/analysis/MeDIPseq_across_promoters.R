
## this script builds a count matrix with the MeDIP-seq data, across a set of genes
##
## 1 - gets the ensembl hg19 genes
## 2 - converts them to promoters
## 3 - reads the MeDIP-seq files, and count the number of extended reads in
##     the promoters regions

rm(list = ls())

library(biomaRt)
library(tidyverse)
library(GenomicAlignments)
library(ChIPpeakAnno)

## Functions to search for ensembl genes
#    listMarts()
#    grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice")
#    listDatasets(grch37)


## Get gene positions for all ensembl genes in hg19

data("TSS.human.GRCh37")
tss = TSS.human.GRCh37
seqlevels(tss) = paste0("chr",seqlevels(tss))
seqlevels(tss,force = TRUE) = paste0("chr",c(seq_len(22),"X","Y"))

## define promoters according to Grimm et al, 2012
upstream = 500
downstream = 1000
fraglen = 300

promoters = promoters(tss, upstream,downstream)
promoters$description = NULL

sizes = read_tsv("/p/keles/SOFTWARE/hg19.chrom.sizes",col_names= FALSE)
genome = GRanges(seqnames = sizes$X1,
                 ranges = IRanges(
                     start = 1,
                     width = sizes$X2))


promoters = subsetByOverlaps(promoters,genome)
promoters$ensembl_gene_id = names(promoters)
names(promoters) = NULL

promoterfile = "data/BED/hg19/MeDIPseq_promoters.bed"

bedfile  = promoters %>% as.data.frame %>% as.tbl %>% .[,c(seq_len(3),5)]

fileConn = file(promoterfile)
writeLines(c("browser position chr22:20100000-20140000",
             'track name=Promoters description="Promoter regions" color=0,0,255'),
           fileConn)
close(fileConn)

write_tsv( bedfile,promoterfile,append = TRUE,
          col_names = FALSE)


## get MeDIP-seq data, we are going to extend the reads by 300 bps
## and count the number of reads in the promoter regions

dipdr = "data/BAM/hg19/bowtie_dip"
dipfiles = list.files(dipdr,pattern = "sort",full.names = TRUE) %>%
    .[grep("bai",.,invert = TRUE)]

## count overlapping  reads into promoters regions
process_reads <- function(file,promoters,fraglen)
{
    w = width(promoters) %>% max
    which = resize(promoters, w + 2 * fraglen,fix = "center")
    reads = readGAlignments(file,param = ScanBamParam(which = which))

    reads = reads %>% as("GRanges") %>% resize(fraglen)
    countOverlaps(promoters,reads) 
    
}    


library(BiocParallel)

bp = MulticoreParam(workers = 16,bpprogressbar = TRUE)

countvectors = bplapply(dipfiles,process_reads,promoters,fraglen,BPPARAM = bp)
names(countvectors) = gsub(".sort.bam","",basename(dipfiles))

outdr = "data/MeDIPseq_results"
save(countvectors,file = file.path(outdr,"MeDIPseq_counts.RData"))

## tidying them up into wide data
count_matrix = bind_cols(
    promoters %>% as.data.frame %>% as.tbl %>%
      select(ensembl_gene_id,everything()) %>%
      mutate(width = NULL,strand = NULL),
    countvectors %>% do.call(cbind,.) %>% as.data.frame %>% as.tbl)

write_tsv(count_matrix ,
          path = file.path(outdr,
                           paste0("MeDIPseq_PromotersCounts_upstr",upstream,
                                  "_downstr",downstream,"fraglen",fraglen,".tsv")))
    
