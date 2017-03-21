
rm(list = ls())

devtools::load_all("~/Desktop/Docs/Code/ChIPUtils")


library(GenomicAlignments)

dr = "/store01/Collab_EJohannsen_RWelch/data_plos_genetics_paper/BAM"
files = list.files(dr,full.names= TRUE)

files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]

library(parallel)
library(magrittr)

options(mc.cores = 22)

reads = mclapply(files,readGAlignments,param = NULL)
reads = reads %>% mclapply(as,"GRanges")

sizes = readr::read_tsv("/p/keles/SOFTWARE/mm9.chrom.sizes",col_names=FALSE)

library(dplyr)

pbc = reads %>% sapply(PBC)

library(BiocParallel)

scc = lapply(reads,strand_cross_corr,seq_len(500),sizes)

names(scc) = paste0("Rep",seq_along(scc))


scc = mapply(function(x,y)x %>% mutate(repl = y),scc,names(scc),SIMPLIFY = FALSE) %>% bind_rows


library(ggthemes)
library(ggplot2)

ggplot(scc,aes(shift,cross.corr,colour = repl))+
    geom_line()+
    scale_color_ptol("cyl")+theme_minimal()

ggsave("scc_mouse_data.png"   )
