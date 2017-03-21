
rm(list = ls())

#library(ChIPUtils)

library(GenomicAlignments)
library(magrittr)
library(BiocParallel)
library(matrixStats)

devtools::load_all("~/Desktop/Docs/Code/ChIPUtils")

dipdr = "data/BAM/hg19/bowtie_dip"
dipfiles = list.files(dipdr,full.names = TRUE,pattern = "sort")
dipfiles = dipfiles[grep("bai",dipfiles,invert = TRUE)]
dipfiles = dipfiles[grep("Input",dipfiles,invert = TRUE)]

options(mc.cores = 20)


rr = dipfiles[4:6] %>% lapply(readGAlignments,param = NULL)

library(readr)
library(dplyr)
library(data.table)

sizes = read_tsv("/p/keles/SOFTWARE/hg19.chrom.sizes",col_names= FALSE)

rr = rr %>% lapply(as,"GRanges")

aa = do.call(getMethod(c,"GenomicRanges"),rr)


scc = strand_cross_corr(aa,shift = seq_len(500),sizes)




rr %>% lapply(PBC)


## > PBC(aa)
## [1] 0.9777749


system.time(
aa <- scc(rr[[1]],seq_len(300),sizes)
)


system.time(
    bb <- PBC(rr[[1]])
    )




reads = dipfiles %>% mclapply(create_reads)


pbc = reads %>% mclapply(PBC) %>% unlist

depth = reads %>% sapply(nreads)


##          min        max     nsc
##        <dbl>      <dbl>   <dbl>
## 1 0.02246122 0.03398293 1.51296

library(dplyr)
tt = tibble(pbc,depth ,sample = gsub(".sort.bam","",basename(dipfiles)))

##          pbc      depth                             sample
## 1  0.9897847 28,371,167 MeDIPseq-NOKS-akata-CaFBS-hmC-rep1
## 2  0.9857791 19,151,870 MeDIPseq-NOKS-akata-CaFBS-hmC-rep2
## 3  0.9820665 21,491,618 MeDIPseq-NOKS-akata-CaFBS-hmC-rep3

## 4  0.9870850 17,646,546  MeDIPseq-NOKS-akata-CaFBS-mC-rep1
## 5  0.9896884 21,918,960  MeDIPseq-NOKS-akata-CaFBS-mC-rep2
## 6  0.9836471 19,720,048  MeDIPseq-NOKS-akata-CaFBS-mC-rep3

## 7  0.9900538 26,701,337  MeDIPseq-NOKS-akata-mono-hmC-rep1
## 8  0.9853771 21,766,903  MeDIPseq-NOKS-akata-mono-hmC-rep2
## 9  0.9861336 12,447,951  MeDIPseq-NOKS-akata-mono-hmC-rep3

## 10 0.9884849 17,538,040   MeDIPseq-NOKS-akata-mono-mC-rep1
## 11 0.9933624  4,548,734   MeDIPseq-NOKS-akata-mono-mC-rep2
## 12 0.9858711 16,257,730   MeDIPseq-NOKS-akata-mono-mC-rep3

## 13 0.9882710 42,746,559       MeDIPseq-NOKS-CaFBS-hmC-rep1
## 14 0.9820927 20,537,563       MeDIPseq-NOKS-CaFBS-hmC-rep2
## 15 0.9795356 15,998,230       MeDIPseq-NOKS-CaFBS-hmC-rep3
## 16 0.9871846 17,761,213        MeDIPseq-NOKS-CaFBS-mC-rep1
## 17 0.9893404 12,066,021        MeDIPseq-NOKS-CaFBS-mC-rep2
## 18 0.9838283 16,470,868        MeDIPseq-NOKS-CaFBS-mC-rep3
## 19 0.9897919 25,958,512        MeDIPseq-NOKS-mono-hmC-rep1
## 20 0.9869623 28,298,338        MeDIPseq-NOKS-mono-hmC-rep2
## 21 0.9832744 16,605,381        MeDIPseq-NOKS-mono-hmC-rep3
## 22 0.9893766 15,868,553         MeDIPseq-NOKS-mono-mC-rep1
## 23 0.9966170     91,351         MeDIPseq-NOKS-mono-mC-rep2
## 24 0.9942169  1,937,624         MeDIPseq-NOKS-mono-mC-rep3

library(readr)
library(data.table)

sizes = read_tsv("/p/keles/SOFTWARE/hg19.chrom.sizes",col_names= FALSE) %>%
    as.data.frame %>% as.data.table %>% rename(V1 = X1,V2 = X2)

shift = seq_len(500)
scc = reads %>% lapply(strand_cross_corr,shift,chrom.sizes = sizes,parallel = TRUE)


scc = scc %>% lapply(as.tbl)

pdf("figs/MeDIPseq_SCC.pdf")
plots = mapply(function(x,y)ggplot(x,aes(shift,cross.corr))+geom_line()+ggtitle(y),scc,gsub(".sort.bam","",basename(dipfiles)),SIMPLIFY = FALSE)
u = plots %>% lapply(print)
dev.off()


scc = mapply(function(x,y)x %>% mutate(file = y),scc,gsub(".sort.bam","",basename(dipfiles)),
             SIMPLIFY = FALSE) %>% bind_rows

write_tsv(scc,"data/quality_control/MeDIPseq_scc.tsv")

##                                  file           min         max        nsc
## 1  MeDIPseq-NOKS-akata-CaFBS-hmC-rep1  5.053178e-03 0.012706688   2.514593
## 2  MeDIPseq-NOKS-akata-CaFBS-hmC-rep2  3.251908e-03 0.007327670   2.253345
## 3  MeDIPseq-NOKS-akata-CaFBS-hmC-rep3  2.615224e-03 0.006611483   2.528075
## 4   MeDIPseq-NOKS-akata-CaFBS-mC-rep1  8.660414e-03 0.014555929   1.680743
## 5   MeDIPseq-NOKS-akata-CaFBS-mC-rep2  7.794165e-03 0.014067769   1.804910
## 6   MeDIPseq-NOKS-akata-CaFBS-mC-rep3  1.116062e-02 0.015895212   1.424222
## 7   MeDIPseq-NOKS-akata-mono-hmC-rep1  4.372276e-03 0.010538927   2.410398
## 8   MeDIPseq-NOKS-akata-mono-hmC-rep2  4.440353e-03 0.008012284   1.804425
## 9   MeDIPseq-NOKS-akata-mono-hmC-rep3  1.520429e-03 0.003818908   2.511732
## 10   MeDIPseq-NOKS-akata-mono-mC-rep1  8.360911e-03 0.014133506   1.690427
## 11   MeDIPseq-NOKS-akata-mono-mC-rep2  1.710977e-03 0.002879452   1.682929
## 12   MeDIPseq-NOKS-akata-mono-mC-rep3  9.789695e-03 0.013689649   1.398373
## 13       MeDIPseq-NOKS-CaFBS-hmC-rep1  5.872459e-03 0.011156080   1.899729
## 14       MeDIPseq-NOKS-CaFBS-hmC-rep2  2.563588e-03 0.005019506   1.958001
## 15       MeDIPseq-NOKS-CaFBS-hmC-rep3  1.469457e-03 0.002968420   2.020079
## 16        MeDIPseq-NOKS-CaFBS-mC-rep1  7.777039e-03 0.012869587   1.654818
## 17        MeDIPseq-NOKS-CaFBS-mC-rep2  4.254199e-03 0.006711667   1.577657
## 18        MeDIPseq-NOKS-CaFBS-mC-rep3  9.533423e-03 0.014640391   1.535691
## 19        MeDIPseq-NOKS-mono-hmC-rep1  3.319376e-03 0.006438094   1.939550
## 20        MeDIPseq-NOKS-mono-hmC-rep2  6.518762e-03 0.009345933   1.433698
## 21        MeDIPseq-NOKS-mono-hmC-rep3  1.179685e-03 0.002719148   2.304978
## 22         MeDIPseq-NOKS-mono-mC-rep1  6.703965e-03 0.011003823   1.641390
## 23         MeDIPseq-NOKS-mono-mC-rep2 -1.677234e-05 0.000177101 -10.559110
## 24         MeDIPseq-NOKS-mono-mC-rep3  1.136981e-03 0.001913271   1.682764
