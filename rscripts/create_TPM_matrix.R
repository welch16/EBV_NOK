#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--all_files",action = "store_true",type = "character",
                help = "All RSEM outcomes to add into the matrix"),
    make_option("--outfile",action = "store_true",type = "character",
                default = tempfile(),
                help = "Name of the file where the TPM matrix is saved")
    
  )

opt = parse_args(OptionParser(option_list = optList))


## opt$all_files = "data/RSEM/hg19/RNAseq-akata-noks-*clone*.genes.results"

separateFiles <- function(ff)
{
  if(grepl(",",ff)){
      ff = ff %>% strsplit(',') %>% unlist     
  }else{
      ff = Sys.glob(ff)      
  }
  ff
}

opt$all_files = separateFiles(opt$all_files)

library(tidyverse)
library(tximport)


## load files with tximport
txdata = tximport(opt$all_files,
                  type = 'rsem' , importer = read_tsv)


rpkmat = txdata$counts / ( txdata$length / 1e3)


scaling_factor = rpkmat %>% colSums(na.rm = TRUE) %>%
    {./1e6}

tpmmat = rpkmat %>% split(rownames(.)) %>%
    map(.f = function(x) x / scaling_factor) %>%
    {do.call(rbind,.)}

tpmmat = tpmmat %>% as.data.frame() %>%
    as.tbl() %>% mutate(genes = rownames(.)) %>%
    select(genes,everything())

colnames(tpmmat) = c("genes",gsub(".genes.results","",basename(opt$all_files)))

tpmmat = tpmmat %>% separate(genes,into = c("ensembl_id","gene_id"), sep = "\\_")

write_tsv(tpmmat,path = opt$outfile)
