
rm(list = ls())

library(tidyverse)
library(MEDIPS)

load("data/MEDIPS/sig_diff_methyl_regions_mC.RData") ## sig_regions

sig_regions_merged = sig_regions %>% map(MEDIPS.mergeFrames,distance = 250) %>%
  map(as.tbl) %>% map2(c("down","up"),.f = function(x,y)x %>% mutate(methyl = y)) %>%
  bind_rows() %>% mutate(start = as.numeric(start),
                         stop = as.numeric(stop),
                         width = stop - start + 1)


library(ChIPpeakAnno)
library(biomaRt)


data("TSS.human.GRCh37")
tss = TSS.human.GRCh37
seqlevels(tss) = paste0("chr",seqlevels(tss))
seqlevels(tss,force = TRUE) = paste0("chr",c(1:22,"X","Y"))

gr = GRanges(seqnames = sig_regions_merged$chr,
             ranges = IRanges(
               start = sig_regions_merged$start,
               end = sig_regions_merged$stop
             ),methyl = sig_regions_merged$methyl)

gr = annotatePeakInBatch(gr,AnnotationData = tss)

# data("TSS.human.GRCh38")
# tss2 = TSS.human.GRCh38
# 
# 
# attr <- c("chromosome_name", "start_position", "end_position",
#           "ensembl_gene_id","external_gene_name","strand","transcription_start_site"
#           )
# 
# bm = getBM(attr, "chromosome_name", c(as.character(1:22),"X","Y"), ensembl)
# 
# 
# rd = with(bm,
#             RangedData(IRanges(start_position, end_position),
#                        space=chromosome_name, ensembl_gene_id,gene = external_gene_name,strand,
#                        tss = transcription_start_site))
# rd = rd %>% as("GRanges")
# 
# 
# rd = with(bm,
#             tibble(seqnames = chromosome_name,start  = start_position,
#                    end =  end_position,strand,
#                    ensembl_gene_id,gene = external_gene_name,
#                    tss = transcription_start_site))
# 
# 
# CLTCL1 = rd  %>% filter(gene == "CLTCL1")
# annot[CLTCL1$ensembl_gene_id]
# tss[CLTCL1$ensembl_gene_id]
# 

