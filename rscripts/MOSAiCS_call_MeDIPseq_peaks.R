#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
  make_option("--bindata_file", action = "store_true",type = "character",
              help = "Files used to make the plots"),
  make_option("--peakfile",action = "store_true",type = "character",default = tempfile(),
              help = "Name of the output file where the peaks are saved."),
  make_option("--peak_sample",action = "store_true",type = "character",
              help = "Name of the sample used to call peaks, the script will determine the Input
                      sample by itself."),
  make_option("--read_dir",action = "store_true",type = "character",
              help = "Directory where the aligned read files are saved"),
  make_option("--fdr",action = "store_true",type = "numeric",default = 0.05,
              help = "Empirical FDR used to call peaks"),
  make_option("--maxgap",action = "store_true",type = "numeric",default = 200,
              help = "Maximum gap considered to merge adjacent enriched bins"),
  make_option("--thresh",action = "store_true",type = "numeric",default = 10,
              help = "Minimum number of extended reads per bin"),
  make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
              help = "Directory and prefix for all the figures to be generated"),              
  make_option("--cores",action = "store_true",type = "numeric",default = 8,
              help = "Number of cores used for parallel")
)

opt = parse_args(OptionParser(option_list = optList))

library(base,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(mosaics,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(purrr,quietly = TRUE)
library(GenomicAlignments,quietly = TRUE)
library(GenomicRanges,quietly = TRUE)
library(parallel,quietly = TRUE)
library(readr,quietly = TRUE)
library(Rsamtools,quietly = TRUE)

options(mc.cores = opt$cores)

## opt$bindata_file = "data/bins/MeDIPseq/BinMatrix_binsize200_fragLen200_InputPooled.tsv"
## opt$peak_sample = "EBV_CaFBS_hmC_Rep1"
## opt$read_dir = "data/BAM/hg19/bowtie_dip"

bin_matrix = read_tsv(opt$bindata_file)
characteristics = strsplit(opt$peak_sample,"_")[[1]]
samples = bin_matrix %>% colnames


get_Input <- function(samples,characteristics)
{
    samples = samples[grepl("Input",samples)]
    if(characteristics[1] == "EBV"){
        samples = samples[grepl("EBV",samples)]
    }else{
        samples = samples[grepl("NOK",samples)]
    }

    if(characteristics[2] == "CaFBS"){
        samples[grepl("CaFBS",samples)]
    }else{
        samples[grepl("NoTrt",samples)]
    }
        
        
}

get_seqDepth <- function(dr,charac)
{
    files = list.files(dr,full.names = TRUE)
    files = files[grep("sort",files)]
    files = files[grep("bai",files,invert = TRUE)]

    if(charac[1] == "EBV"){
        files = files[grep("akata",files)]
    }else{
        files = files[grep("akata",files,invert = TRUE)]
    }

    files = files[grep(charac[2],files)]
    files = files[grep(charac[3],files)]
    if(charac[3] != "hmC"){
        files = files[grep("hmC",files,invert = TRUE)]
    }
    files = files[grep(charac[4],files,ignore.case = TRUE)]
    
    bf = BamFile(files)

    sl = seqlengths(bf)
    chr = sl %>% names

    counts = chr %>% mclapply(function(x){
        countBam(bf,param = ScanBamParam(which = GRanges(x,IRanges(1,sl[x]))))
        })

    counts = counts %>% sapply(function(x)x$records)
    names(counts) = chr

    counts
                               
}
                        

## binSize = get_binSize(bin_matrix)
input_name = get_Input(samples,characteristics)
seqdepth = get_seqDepth(opt$read_dir,characteristics)

## exclude chrM
bin_matrix = bin_matrix %>% filter(seqnames != "chrM")
seqdepth = seqdepth[names(seqdepth) != "chrM"]

## seqdepth = seqdepth[names(seqdepth) %in% c("chr10","chr11")]
## bin_matrix = bin_matrix %>% filter(seqnames %in% c("chr10","chr11"))


coord = bin_matrix %>% select(seqnames,start) %>%
    dplyr::rename(chrID =seqnames,coord = start) %>%
    mutate(coord = coord - 1)

input = bin_matrix %>% select(dplyr::contains(input_name))

chip = bin_matrix %>%
    select(dplyr::contains(ifelse(characteristics[1] == "EBV","EBV","NOK"))) %>%
    select(dplyr::contains(ifelse(characteristics[2] == "mono","NoTrt","CaFBS"))) %>%
    select(dplyr::contains(characteristics[4]))

if(characteristics[3] == "hmC"){
    chip = chip %>% select(dplyr::contains("hmC"))
}else{
    chip = chip %>% select(-dplyr::contains("hmC"))
}


bindata = new("BinData",
              chrID = coord$chrID,
              coord = coord$coord %>% as.integer,
              tagCount = chip[[1]],
              input = input[[1]],
              dataType = c("chip","input"),
              seqDepth = seqdepth %>% sum)

df = bind_cols(chip,input)
names(df) = c("chip","input")

fit = mosaicsFit(bindata,analysisType = "IO",parallel = TRUE,
                 nCore = opt$cores)

pdf(paste0(opt$figs,"_mosaicsGOF.pdf"))
plot(fit)
dev.off()

source("rfuns/MeDIP_analysis.R")

library(ggplot2,quietly = TRUE)

theme_set(theme_bw())

df = df %>% mutate_all(funs(log10(1  +. )))

pdf(paste0(opt$figs,"_hexbin_Input_vs_ChIP.pdf"))
hexbin_scatter_plot(df,"input","chip")+coord_fixed() + ggtitle( chip %>% colnames)
dev.off()

peaks = mosaicsPeak( fit, signalModel="2S", FDR= opt$fdr / 100,
                        maxgap= opt$maxgap, thres= opt$thresh )

peaks = print(peaks) %>% as.tbl


write_tsv(peaks,path = opt$peakfile)

