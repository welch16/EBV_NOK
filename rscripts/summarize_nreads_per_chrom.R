#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
  make_option("--genomeBam", action = "store_true", type = "character",
              help = "Name of the bam file with the reads aligned to the genome"),
  make_option("--outfile",action = "store_true",type = "character",
              help = "File where the nr. of reads per chr are saved"))

opt = parse_args(OptionParser(option_list = optList))

## opt$genomeBam = "data/BAM/hg19/bowtie_genome/RNAseq-Noks-mono-rep3_rsem_default.bam"

library(GenomicAlignments)
library(tidyverse)


reads = readGAlignments(opt$genomeBam,param = NULL)
reads = as(reads,"GRanges")

summaryTable = reads %>% as.data.frame %>% as.tbl %>%
    group_by(seqnames) %>% summarize(nreads = n(),
                                     freads = sum(strand == "+"),
                                     readlength = median(width))

write_tsv(summaryTable,opt$outfile)
