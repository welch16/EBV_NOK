
rm(list = ls())

library(tidyverse)
library(MEDIPS)
library(GenomicRanges)

load(file = "./data/MEDIPS/mC_edgeR.RData") ## mC_edgeR_results

names(mC_edgeR_results)[1:3] = c("seqnames","start","end")

mC_edgeR_results = DataFrame(mC_edgeR_results)

gr = makeGRangesFromDataFrame(mC_edgeR_results )

results = mC_edgeR_results %>% as.data.frame %>% as.tbl


results %>% ggplot(aes(edgeR.logFC)) + geom_histogram(bins = 50) + xlim(-5,5)
## there are tails, slightly upmethylated, centered around zero


results %>% filter(CF > 0) %>% ggplot(aes(edgeR.p.value))+geom_histogram(bins = 50)

library(viridis)
library(hexbin)
library(scales)

pal = viridis(1e3,option = "D")


MEDIPS_plot <- function(DT,xvar,yvar)
{
  require(viridis)
  require(scales)
  require(ggplot2)
  pal = viridis(1e3, option = "D")
  
  ggplot(DT,aes_string(xvar,yvar))+stat_binhex(bins = 80) +
    scale_fill_gradientn(colours = pal,trans = "log10",
                         labels = trans_format('log10',math_format(10^.x)),
                         guide = guide_colourbar(title = NULL,
                           barheight = unit(0.92,"npc"),
                           barwidth = unit(0.01,"npc")))
                          
}

results = results %>% mutate(MSets.log.mean =  log(.5 * (MSets1.rpkm.mean + MSets2.rpkm.mean)))
nms = names(results)


ident = geom_abline(slope = 1,intercept = 0,linetype = 2,colour = "red")

theme_set(theme_bw())


## mC_edgeR_results = MEDIPS.meth(MSet1 = mC_sets[grepl("akata",names(mC_sets))],
##                                  MSet2 = mC_sets[!grepl("akata",names(mC_sets))],
##                                  CSet = couplingVec,
##                                  ISet1 = input_sets[grepl("akata",names(input_sets))],
##                                  ISet2 = input_sets[!grepl("akata",names(input_sets))],
##                                  p.adj = "BH", diff.method = "edgeR",
##                                 MeDIP = TRUE, CNV = FALSE, minRowSum = minRowSum)

figsdr = "figs/diff_methylation/MEDIPS_results"

pdf(file.path(figsdr,"aveRPKM_EBV_vs_NOKS_bin200_fragLen300.pdf"))
MEDIPS_plot(results,nms[41],nms[44])+ident+xlim(0,20)+ylim(0,20)+
    xlab("mean RPKM EBV")+ylab("mean RPKM NOKS")
dev.off()



pdf(file.path(figsdr,"aveRPKM_EBV_vs_NOKS_bin200_fragLen300_Input.pdf"))
MEDIPS_plot(results,nms[47],nms[48])+ident+xlim(0,20)+ylim(0,70)+
    xlab("mean RPKM EBV Input")+ylab("mean RPKM NOKS Input")
dev.off()


pdf(file.path(figsdr,"aveRPKM_Input_vs_MC_bin200_fragLen300_EBV.pdf"))
MEDIPS_plot(results,nms[47],nms[41])+ident+xlim(0,20)+ylim(0,20)+
    xlab("mean RPKM EBV Input")+ylab("mean RPKM EBV")
dev.off()


pdf(file.path(figsdr,"aveRPKM_Input_vs_MC_bin200_fragLen300_NOKS.pdf"))
MEDIPS_plot(results,nms[48],nms[44])+ident+xlim(0,40)+ylim(0,20)+
    xlab("mean RPKM NOKS Input")+ylab("mean RPKM NOKS")
dev.off()


pdf(file.path(figsdr,"mC_EBV_NOKS_MAplot.pdf"))
MEDIPS_plot(results,nms[54],nms[50])+ylim(-5,5)+xlim(0,5)+geom_smooth(se = FALSE)
dev.off()

annot = MEDIPS.getAnnotation(dataset = "hsapiens_gene_ensembl",
                             annotation = c("TSS","EXON","GENE"))

sig_regions = MEDIPS.selectSig( mC_edgeR_results,p.value = .1) %>%
    mutate( upmeth = edgeR.logFC > 0) %>% split(.$upmeth) 

load("data/MEDIPS/sig_diff_methyl_regions_mC.RData")


## down_methyl

down_bed = MEDIPS.mergeFrames(sig_regions[[1]],distance = 250)

down_meth = annot %>% names %>% map( .f = function(x){
    MEDIPS.setAnnotation(down_bed,annot[[x]]) %>% as.tbl})


down_meth = down_meth %>% map(.f = function(x){
    x = x %>% mutate(ID = as.numeric(gsub("ID_","",ID)),
                     start = as.numeric(start),
                     stop = as.numeric(stop))
    ids = x %>% split(.$ID) %>% map_chr(.f = function(y){
        y %>% select(contains("id",ignore.case = FALSE)) %>% do.call(c,.) %>%
            {.[!is.na(.)]} %>% paste(collapse = ",")
    })
    nids = ids %>% map_dbl(.f = function(z){
        z %>% strsplit(",") %>% .[[1]] %>% length
    })
   x %>% select(-contains("id",ignore.case = FALSE)) %>%
       mutate(annot = ifelse(ids == "",NA,ids),
              nannot = nids,
              width = stop -start + 1)  
})



up_bed = MEDIPS.mergeFrames(sig_regions[[2]],distance = 250)


up_meth = annot %>% names %>% map( .f = function(x){
    MEDIPS.setAnnotation(up_bed,annot[[x]]) %>% as.tbl})


up_meth = up_meth %>% map(.f = function(x){
    x = x %>% mutate(ID = as.numeric(gsub("ID_","",ID)),
                     start = as.numeric(start),
                     stop = as.numeric(stop))
    ids = x %>% split(.$ID) %>% map_chr(.f = function(y){        
        y %>% select(contains("id",ignore.case = FALSE)) %>% do.call(c,.) %>%
            {.[!is.na(.)]} %>% paste(collapse = ",")
    })
    nids = ids %>% map_dbl(.f = function(z){
        z %>% strsplit(",") %>% .[[1]] %>% length
    })
   x %>% select(-contains("id",ignore.case = FALSE)) %>%
       mutate(annot = ifelse(ids == "",NA,ids),
              nannot = nids,
              width = stop - start + 1)  
})


outdr = "data/BED/hg19"

upfile = file.path(outdr,"MEDIPS_upmethyl_mC.bed")
dofile = file.path(outdr,"MEDIPS_downmethyl_mC.bed")

bedfiles  = list()
bedfiles[[1]] = up_bed %>% as.tbl %>% select(chr,start,stop)
bedfiles[[2]] = down_bed %>% as.tbl %>% select(chr,start,stop)

fileConn = file(upfile)
writeLines(c("browser position chr22:20100000-20140000",
             'track name=Upmethy description="MEDIPS Up-methylated regions" color=102,0,102'),
           fileConn)
close(fileConn)

write_tsv(bedfiles[[1]],upfile,append = TRUE,
          col_names = FALSE)

fileConn = file(dofile)
writeLines(c("browser position chr22:20100000-20140000",
             'track name=Downmethy description="MEDIPS Down-methylated regions" color=255,153,153'),
           fileConn)
close(fileConn)

write_tsv(bedfiles[[2]],dofile,append = TRUE,
          col_names = FALSE)

save(sig_regions,file = "data/MEDIPS/sig_diff_methyl_regions_mC.RData")
