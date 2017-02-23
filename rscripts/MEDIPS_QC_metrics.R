#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--readsfile",action = "store_true",type = "character",
                help = "File with MeDIP-seq aligned reads"),
    make_option("--outfile",action = "store_true",type = "character",default = tempfile(),
                help = "Output file with summary QC metrics"),
    make_option("--shift",action = "store_true",default = 0 , type = "numeric",
                help = "Nr. of bps that the reads are going to be shifted toward the 3' end"),
    make_option("--extend",action = "store_true",default = 200,type = "numeric",
                help = "Nr. of bps that the reads are going to be extended toward the 3' end"),
    make_option("--window_size",action = "store_true",default = 100,type = "numeric",
                "Window size, used to partition the genome into bins of size 'window_size'"),
    make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
                help = "Directory and prefix for all the figures to be generated"),
    make_option("--cores",action = "store_true",type = "numeric",default = 8,
              help = "Number of cores used for parallel")
)

opt = parse_args(OptionParser(option_list = optList))

## opt$readsfile = "data/BAM/hg19/bowtie_dip/MeDIPseq-NOKS-mono-hmC-rep1.sort.bam"
## opt$window_size = 200

library(base)

suppressMessages( library(BSgenome))
suppressMessages(library(MEDIPS))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(parallel))

options(mc.cores = opt$cores)

## saturation analysis


sr = MEDIPS.saturation(file = opt$readsfile, BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
                       uniq = 0, extend = opt$extend, shift = opt$shift, window_size = opt$window_size,
                       nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)

## CpG coverage
cr = MEDIPS.seqCoverage(file = opt$readsfile , pattern = "CG",
                        BSgenome = "BSgenome.Hsapiens.UCSC.hg19", extend = opt$extend, shift = opt$shift,
                        uniq = 0)


## enrichment analysis
er = MEDIPS.CpGenrich(file = opt$readsfile, BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
                      extend = opt$extend, shift = opt$shift, uniq = 0)

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(scales))
suppressMessages(library(tidyr))
    
name = gsub(".sort.bam","",opt$readsfile) %>% basename
theme_set(theme_bw())


plot_saturation <- function(medips)
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


plot_coverage <- function(medips)
{

    df = tibble(cover = medips$cov.res)

    df %>% ggplot(aes(cover)) +
        geom_histogram(aes(y = ..count../sum(..count..)),binwidth = 1,fill= "white",colour = "black")+
        scale_x_continuous(limits = c(0,30))+
        scale_y_continuous(labels = comma)+
        xlab("Number of CpG's covered")+ylab("Frequency")+ggtitle(name)
        

}


pdf(paste0(opt$figs,"_saturation.pdf"),height = 4)
plot_saturation(sr) %>% print
dev.off()


pdf(paste0(opt$figs,"_coverage.pdf"))
plot_coverage(cr) %>% print
dev.off()

summary = tibble(name,
                 depth = sr$numberReads,
                 maxEstCor = sr$maxEstCor[2],
                 maxTruCor = sr$maxTruCor[2],
                 pattern = cr$pattern,
                 depthPattern = cr$numberReadsWO)

er = er %>% as.data.frame %>% as.tbl

summary = bind_cols(summary,er)



readr::write_tsv(summary, opt$outfile)                 

