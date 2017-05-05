
rm(list = ls())

library(tidyverse)
library(viridis)
library(hexbin)
library(scales)
library(gridExtra)

clean_results <- function(x)
{
  x = x %>% as.tbl()
  base = x %>% 
    select(chr,start,stop,CF) 
  
  counts = x %>% select(contains("count"))
  
  edgeR = x %>% select(contains("edge"))
  
  bind_cols(base,counts,edgeR) 
  
}

load("data/MEDIPS/mC_edgeR.RData") ## the file is too big, we need to remove tons of columns 
## before being able to work with it
mC_results = clean_results(mC_edgeR_results)
rm(mC_edgeR_results)

load("data/MEDIPS/hmC_edgeR.RData")

hmC_results = clean_results(hmC_edgeR_results)
rm(hmC_edgeR_results)

MEDIPS_plot <- function(x,y,bb = 50)
{
  pal = viridis(1e3, option = "D")
  
  DT  = tibble(x,y)
  
  ggplot(DT,aes_string("x","y"))+stat_binhex(bins = bb) +
    scale_fill_gradientn(colours = pal,trans = "log10",
                         labels = trans_format('log10',math_format(10^.x)),
                         guide = guide_colourbar(title = NULL,
                                                 barheight = unit(0.92,"npc"),
                                                 barwidth = unit(0.01,"npc")))
  
}

redline = geom_abline(slope = 1 , intercept = 0,linetype = 2,colour = "red")

figsdr = "figs/diff_methylation/MEDIPS_results"

pdf(file.path(figsdr,"mC_EBV_NOKS_meanCounts.pdf"))
MEDIPS_plot(mC_results$MSets1.counts.mean,
            mC_results$MSets2.counts.mean)+
  ggtitle("mC")+xlab("EBV")+ylab("NOKS")+
  ylim(0,125)+xlim(0,125)+redline
dev.off()

pdf(file.path(figsdr,"input_EBV_NOKS_meanCounts.pdf"))
MEDIPS_plot(mC_results$ISets1.counts.mean,
            mC_results$ISets2.counts.mean)+
  ggtitle("Input")+xlab("EBV")+ylab("NOKS")+
  ylim(0,125)+xlim(0,125)+redline
dev.off()

pdf(file.path(figsdr,"hmC_EBV_NOKS_meanCounts.pdf"))
MEDIPS_plot(hmC_results$MSets1.counts.mean,
            hmC_results$MSets2.counts.mean)+
  ggtitle("hmC")+xlab("EBV")+ylab("NOKS")+
  ylim(0,125)+xlim(0,125)+redline
dev.off()


pdf(file.path(figsdr,"all_input_mc.pdf"))
MEDIPS_plot(mC_results$MeDIPseq.NOKS.akata.CaFBS.Input.rep0.sort.bam.counts,
            mC_results$MeDIPseq.NOKS.CaFBS.Input.rep0.sort.bam.counts)+
  ggtitle("CaFBS")+xlab("EBV")+ylab("NOKS")+xlim(0,100)+ylim(0,100)
MEDIPS_plot(mC_results$MeDIPseq.NOKS.akata.mono.Input.rep0.sort.bam.counts,
            mC_results$MeDIPseq.NOKS.mono.Input.rep0.sort.bam.counts)+
  ggtitle("mono")+xlab("EBV")+ylab("NOKS")+xlim(0,100)+ylim(0,100)
MEDIPS_plot(mC_results$MeDIPseq.NOKS.akata.mono.Input.rep0.sort.bam.counts,
            mC_results$MeDIPseq.NOKS.akata.CaFBS.Input.rep0.sort.bam.counts)+
  ggtitle("EBV")+xlab("mono")+ylab("CaFBS")+xlim(0,100)+ylim(0,100)
MEDIPS_plot(mC_results$MeDIPseq.NOKS.mono.Input.rep0.sort.bam.counts,
            mC_results$MeDIPseq.NOKS.CaFBS.Input.rep0.sort.bam.counts)+
ggtitle("NOKS")+xlab("mono")+ylab("CAFBS")+xlim(0,100)+ylim(0,100)
dev.off()

pdf(file.path(figsdr,"all_input_hmc.pdf"))
MEDIPS_plot(hmC_results$MeDIPseq.NOKS.akata.CaFBS.Input.rep0.sort.bam.counts,
            hmC_results$MeDIPseq.NOKS.CaFBS.Input.rep0.sort.bam.counts)+
  ggtitle("CaFBS")+xlab("EBV")+ylab("NOKS")+xlim(0,100)+ylim(0,100)
MEDIPS_plot(hmC_results$MeDIPseq.NOKS.akata.mono.Input.rep0.sort.bam.counts,
            hmC_results$MeDIPseq.NOKS.mono.Input.rep0.sort.bam.counts)+
ggtitle("mono")+xlab("EBV")+ylab("NOKS")+xlim(0,100)+ylim(0,100)
MEDIPS_plot(hmC_results$MeDIPseq.NOKS.akata.mono.Input.rep0.sort.bam.counts,
            hmC_results$MeDIPseq.NOKS.akata.CaFBS.Input.rep0.sort.bam.counts)+
  ggtitle("EBV")+xlab("mono")+ylab("CaFBS")+xlim(0,100)+ylim(0,100)
MEDIPS_plot(hmC_results$MeDIPseq.NOKS.mono.Input.rep0.sort.bam.counts,
            hmC_results$MeDIPseq.NOKS.CaFBS.Input.rep0.sort.bam.counts)+
  ggtitle("NOKS")+xlab("mono")+ylab("CAFBS")+xlim(0,100)+ylim(0,100)
dev.off()




