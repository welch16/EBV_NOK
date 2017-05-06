
rm(list = ls())

library(tidyverse)

load("data/MEDIPS/DiffMethyl_annotation.RData")

ucsc_url = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=welch16&hgS_otherUserSessionName=MC_tracks&position="
context = 1e3

diff_meth_anno = diff_meth_anno %>%
  filter(methSet %in% c("mC.up","mC.down")) %>%
  mutate(
    edgeR.logCPM = NULL,
    strand = NULL,
    width = NULL,min.pval = NULL,med.pval = NULL,
    geneId = NULL
  ) %>% 
  mutate_if(is.logical,funs(ifelse(.,1,0))) %>%
  mutate(url = paste0(ucsc_url, seqnames,":",start - context ,"-",end + context)) %>%
  select(url , everything()) %>% 
  select(-contains("rpkm")) %>%
  select(-contains("rms")) %>%
  rename(chrom = seqnames) %>% 
  mutate(geneStrand = ifelse(geneStrand == "+","F","R")) %>%
  split(.$methSet) %>%
  map(mutate,methSet = NULL) %>%
  map( .f = function(x)x %>% arrange(desc(medlog10pval)))

dr = "data/MEDIPS/UCSC_link"

dir.create(dr,showWarnings = FALSE)

diff_meth_anno %>% 
  map2(
    names(diff_meth_anno) %>% 
      {gsub(".up","_EBV_dominant",.)} %>%
      {gsub(".down","_NOKS_dominant",.)},
    .f = function(tab,pref){
      write_tsv(tab,path = file.path(dr,paste0(pref,".tsv")))
    })
  
  
  
load("data/MEDIPS/DiffMethyl_annotation.RData")

ucsc_url = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=welch16&hgS_otherUserSessionName=hMC_tracks&position="
context = 1e3

diff_meth_anno = diff_meth_anno %>%
  filter(methSet %in% c("hmC.up","hmC.down")) %>%
  mutate(
    edgeR.logCPM = NULL,
    strand = NULL,
    width = NULL,min.pval = NULL,med.pval = NULL,
    geneId = NULL
  ) %>% 
  mutate_if(is.logical,funs(ifelse(.,1,0))) %>%
  mutate(url = paste0(ucsc_url, seqnames,":",start - context ,"-",end + context)) %>%
  select(url , everything()) %>% 
  select(-contains("rpkm")) %>%
  select(-contains("rms")) %>%
  rename(chrom = seqnames) %>% 
  mutate(geneStrand = ifelse(geneStrand == "+","F","R")) %>%
  split(.$methSet) %>%
  map(mutate,methSet = NULL) %>%
  map( .f = function(x)x %>% arrange(desc(medlog10pval)))

dr = "data/MEDIPS/UCSC_link"

dir.create(dr,showWarnings = FALSE)

diff_meth_anno %>% 
  map2(
    names(diff_meth_anno) %>% 
    {gsub(".up","_EBV_dominant",.)} %>%
    {gsub(".down","_NOKS_dominant",.)},
    .f = function(tab,pref){
      write_tsv(tab,path = file.path(dr,paste0(pref,".tsv")))
    })







