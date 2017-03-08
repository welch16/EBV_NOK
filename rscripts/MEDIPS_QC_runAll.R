
rm(list = ls())

library(base)
library(methods)
library(stats4)
library(stats)
library(parallel)
library(compiler)

library(BSgenome)
library(Rsamtools)
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicAlignments)
library(SummarizedExperiment)
library(Biobase)


library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
    
theme_set(theme_bw())


indir = "data/BAM/hg19/bowtie_dip/"
files = list.files(indir,full.names = TRUE)

files = files[grep("sort.bam",files)]
files = files[grep("bai",files,invert = TRUE)]


## parameters
extend = 200
shift = 0
window_size = 200
uniq = 0
genome = "BSgenome.Hsapiens.UCSC.hg19"

theme_set(theme_bw())


plot_saturation <- function(medips,name)
{
    ds = medips$distinctSets
    est = medips$estimation
    
    df = bind_rows(
        tibble(nreads = ds[,1],corr = ds[,2],ps = "Observed"),
        tibble(nreads = est[,1],corr = est[,2],ps = "Estimated"))

    ggplot(df,aes(nreads,corr,colour = ps)) + geom_line() +
        scale_color_brewer(palette = "Set1",name = "")+ylim(0,1)+
        theme(legend.position = "top")+xlab("Number of reads")+
        scale_x_continuous(labels = comma)+
        ylab("Pearson Correlation")+ggtitle(name)
}


plot_coverage <- function(medips,name)
{

    df = tibble(cover = medips$cov.res)

    df %>% ggplot(aes(cover)) +
        geom_histogram(aes(y = ..count../sum(..count..)),binwidth = 1,fill= "white",colour = "black")+
        scale_x_continuous(limits = c(0,30))+
        scale_y_continuous(labels = comma)+
        xlab("Number of CpG's covered")+ylab("Frequency")+ggtitle(name)       
}



medips_analysis <- function(readsfile,genome,extend,shift,window_size,uniq)
{
    name = gsub(".sort.bam","",readsfile) %>% basename
    message("Calculating saturation ",name)
    sr = MEDIPS.saturation(file = readsfile, BSgenome = genome,
                           uniq = uniq, extend = extend,
                           shift = shift, window_size = window_size,
                           nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)

    ## CpG coverage
    message("Calculating coverage ",name)
    cr = MEDIPS.seqCoverage(file = readsfile , pattern = "CG",
                            BSgenome = genome,
                            extend = extend, shift = shift,
                            uniq = uniq)

    ## enrichment analysis
    message("Calculating enrichment ",name)
    er = MEDIPS.CpGenrich(file = readsfile, BSgenome = genome,
                          extend = extend, shift = shift, uniq = 0)

    out = list()
    
    out[["satu"]] = plot_saturation(sr,name)
    out[["cover"]] = plot_coverage(cr,name)

    summary = tibble(name,
                 depth = sr$numberReads,
                 maxEstCor = sr$maxEstCor[2],
                 maxTruCor = sr$maxTruCor[2],
                 pattern = cr$pattern,
                 depthPattern = cr$numberReadsWO)

    er = er %>% as.data.frame %>% as.tbl
    out[["summary"]] = bind_cols(summary,er)

    out   
}

analysis = files %>% lapply(medips_analysis,genome,extend,shift,window_size,uniq)


saturat = analysis %>% lapply(function(x)x[[1]])
cover = analysis %>% lapply(function(x)x[[2]])
summar = analysis %>% lapply(function(x)x[[3]])


figsdr = "figs/methylation/medips"

pdf(file.path(figsdr,"saturation.pdf"),height = 5)
u = saturat %>% lapply(print)
dev.off()

pdf(file.path(figsdr,"coverage.pdf"),height = 5)
u = cover %>% lapply(print)
dev.off()

readr::write_tsv(summar %>% bind_rows,
           "data/quality_control/MEDIPS/MEDIPseq_summary.tsv")







