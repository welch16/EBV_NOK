
rm(list = ls())

library(ChIPUtils)
library(parallel)
library(magrittr)

dipdr = "data/BAM/hg19/bowtie_dip"
dipfiles = list.files(dipdr,full.names = TRUE,pattern = "sort")
dipfiles = dipfiles[grep("bai",dipfiles,invert = TRUE)]
dipfiles = dipfiles[grep("Input",dipfiles,invert = TRUE)]

options(mc.cores = 20)

reads = dipfiles %>% mclapply(create_reads)


pbc = reads %>% mclapply(PBC) %>% unlist

depth = reads %>% sapply(nreads)


library(dplyr)
tt = tibble(pbc,depth ,sample = gsub(".sort.bam","",basename(dipfiles)))

##          pbc    depth                             sample
## 1  0.9897847 28371167 MeDIPseq-NOKS-akata-CaFBS-hmC-rep1
## 2  0.9857791 19151870 MeDIPseq-NOKS-akata-CaFBS-hmC-rep2
## 3  0.9820665 21491618 MeDIPseq-NOKS-akata-CaFBS-hmC-rep3
## 4  0.9870850 17646546  MeDIPseq-NOKS-akata-CaFBS-mC-rep1
## 5  0.9896884 21918960  MeDIPseq-NOKS-akata-CaFBS-mC-rep2
## 6  0.9836471 19720048  MeDIPseq-NOKS-akata-CaFBS-mC-rep3
## 7  0.9900538 26701337  MeDIPseq-NOKS-akata-mono-hmC-rep1
## 8  0.9853771 21766903  MeDIPseq-NOKS-akata-mono-hmC-rep2
## 9  0.9861336 12447951  MeDIPseq-NOKS-akata-mono-hmC-rep3
## 10 0.9884849 17538040   MeDIPseq-NOKS-akata-mono-mC-rep1
## 11 0.9933624  4548734   MeDIPseq-NOKS-akata-mono-mC-rep2
## 12 0.9858711 16257730   MeDIPseq-NOKS-akata-mono-mC-rep3
## 13 0.9882710 42746559       MeDIPseq-NOKS-CaFBS-hmC-rep1
## 14 0.9820927 20537563       MeDIPseq-NOKS-CaFBS-hmC-rep2
## 15 0.9795356 15998230       MeDIPseq-NOKS-CaFBS-hmC-rep3
## 16 0.9871846 17761213        MeDIPseq-NOKS-CaFBS-mC-rep1
## 17 0.9893404 12066021        MeDIPseq-NOKS-CaFBS-mC-rep2
## 18 0.9838283 16470868        MeDIPseq-NOKS-CaFBS-mC-rep3
## 19 0.9897919 25958512        MeDIPseq-NOKS-mono-hmC-rep1
## 20 0.9869623 28298338        MeDIPseq-NOKS-mono-hmC-rep2
## 21 0.9832744 16605381        MeDIPseq-NOKS-mono-hmC-rep3
## 22 0.9893766 15868553         MeDIPseq-NOKS-mono-mC-rep1
## 23 0.9966170    91351         MeDIPseq-NOKS-mono-mC-rep2
## 24 0.9942169  1937624         MeDIPseq-NOKS-mono-mC-rep3

library(readr)
library(data.table)

sizes = read_tsv("/p/keles/SOFTWARE/hg19.chrom.sizes",col_names= FALSE) %>%
    as.data.frame %>% as.data.table %>% rename(V1 = X1,V2 = X2)

shift = seq_len(300)
scc = reads %>% lapply(strand_cross_corr,shift,chrom.sizes = sizes,parallel = TRUE)


scc = scc %>% lapply(as.tbl)

pdf("figs/MeDIPseq_SCC.pdf")
plots = mapply(function(x,y)ggplot(x,aes(shift,cross.corr))+geom_line()+ggtitle(y),scc,gsub(".sort.bam","",basename(dipfiles)),SIMPLIFY = FALSE)
u = plots %>% lapply(print)
dev.off()


scc = mapply(function(x,y)x %>% mutate(file = y),scc,gsub(".sort.bam","",basename(dipfiles)),
             SIMPLIFY = FALSE) %>% bind_rows

write_tsv(scc,"data/quality_control/MeDIPseq_scc.tsv")
