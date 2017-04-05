#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--bamfile",action = "store_true",type = "character",
                help = "Name of the BAM file to convert to wig"),
    make_option("--fraglen",action = "store_true",type = "numeric",default = 300,
                help = "Fragment length used to extend the aligned reads in the BAM file"),
    make_option("--bigwigfile",action = "store_true",
                type = "character",
                help = "Name of the BW file where results are saved")
    )

opt = parse_args(OptionParser(option_list = optList))

## opt$bamfile = "data/BAM/hg19/bowtie_dip/MeDIPseq-NOKS-akata-CaFBS-hmC-rep1.sort.bam"
## opt$bigwigfile = "./out.bw"

suppressMessages(library(GenomicAlignments))
suppressMessages(library(rtracklayer))

reads = readGAlignments(opt$bamfile,param = NULL)
reads = as(reads,"GRanges")

## removing chrM as we never use it 
lvs = seqlevelsInUse(reads)
seqlevels(reads,force = TRUE) = lvs[lvs != "chrM"]

## resizing reads by fraglen and normalizing the coverage
reads = resize(reads ,opt$fraglen)
nreads = length(reads)
cover = coverage(reads)
normFactor = 1e9 / nreads
out = as(cover,"RangedData")
out$score = round(normFactor * out$score,2)


## formatting and saving as bigwig
out = as(out,"GRanges")
seqlengths(out) = seqlengths(reads)
bw = BigWigFile(opt$bigwigfile)
export.bw(out,bw)
