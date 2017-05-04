
rm(list = ls())

library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFiles)
library(tidyverse)


library(parallel)
options(mc.cores = 20)

indr = "data/BAM/hg19/bowtie_dip"
files = list.files(indr,full.names = TRUE,pattern = "sort")

files = files[grep("bai",files,invert = TRUE)]
files = files[grep("rep0",files)]

sizes = read_tsv("/p/keles/SOFTWARE/hg19.chrom.sizes",col_names = FALSE) %>%
    filter(X1 != "chrM")

binsize = 200
fraglen = 300
shift = 0

count_bins <- function(file,bins,fraglen)
{
    rr = readGAlignments(file,param = NULL)
    rr = as(rr,"GRanges")
    rr = resize(rr,fraglen)
    countOverlaps(bins,rr)   
}


devtools::load_all("~/Desktop/Docs/Code/ChIPUtils")
bins = ChIPUtils::create_bins(binsize,sizes )

count_vectors = mclapply(files,count_bins,bins,fraglen)
names(count_vectors) = gsub(".sort.bam","",basename(files))

mindepth = 1

depths = readr::read_tsv("data/metadata/MeDIPseq_sequencingDepth.tsv") 
depths = depths %>% 
    filter(!grepl("rep0",file)) %>%
    filter(depth > mindepth) %>%
    mutate(
        cell = ifelse(grepl("akata",file),"EBV","NOKS"),
        treat= ifelse(grepl("CaFBS",file),"CaFBS","NoTr"),
        meth = ifelse(grepl("Input",file),"Input",
               ifelse(grepl("hmC",file),"hmC","mC"))) %>%
    group_by(cell,treat,meth) %>%
    summarize(
        total_depth = sum(depth),
        nrep = length(file)
    ) %>% ungroup() %>% 
    mutate(file = basename(files)) %>%
    mutate(file = gsub(".sort.bam","",file))

rpkm_vectors = count_vectors %>%
    map2(depths$total_depth,
         .f = function(x,y)(1e9/(binsize* y)) * x)
             
pal = viridis::viridis(1e3)


hexbin_plot <- function(x_idx,y_idx,rpkm_vectors,pal,depths,l1,l2)
{
    xbins = rpkm_vectors[x_idx]
    ybins = rpkm_vectors[y_idx]

    bins = xbins %>% map2(ybins,.f = function(x,y)tibble(x,y))

    bins = bins %>%
        map2(names(xbins),
             .f = function(x,z)
             {
                 x %>% mutate(file = z) %>%
                     left_join(depths,by = "file") %>%
                     mutate(
                         file = NULL,
                         total_depth = NULL,
                         nrep = NULL,
                         meth = NULL)
             }) %>% bind_rows()


    ggplot(bins,aes_string("x","y"))+stat_binhex(bins = 80)+
        scale_x_continuous(limits = c(0,l1))+
        scale_y_continuous(limits = c(0,l2))+
        scale_fill_gradientn(colours = pal,trans = "log10",
                             labels = trans_format('log10',math_format(10^.x)),
                             guide = guide_colourbar(title = NULL,
                                                     barheight = unit(0.92,"npc"),
                                                     barwidth = unit(0.01,"npc")))+
        facet_grid(treat ~ cell)+geom_abline(slope =1,intercept = 0,linetype = 2,colour = "red")
}

library(scales)

input = c(2,5,8,11)
mc = input + 1
hmc = input - 1

figsdr = "figs/methylation/bins"

pdf(file.path(figsdr,"Hexbin_plot_mC_vs_hmC.pdf"))
hexbin_plot(hmc,mc,rpkm_vectors,pal,depths,15,15)+xlab("hmC")+ylab("mC") %>% print()
hexbin_plot(hmc,mc,rpkm_vectors,pal,depths,25,25)+xlab("hmC")+ylab("mC") %>% print()
hexbin_plot(hmc,mc,rpkm_vectors,pal,depths,50,50)+xlab("hmC")+ylab("mC") %>% print()
dev.off()

pdf(file.path(figsdr,"Hexbin_plot_Input_vs_mC.pdf"))
hexbin_plot(input,mc,rpkm_vectors,pal,depths,10,15)+xlab("Input")+ylab("mC") %>% print()
hexbin_plot(input,mc,rpkm_vectors,pal,depths,15,20)+xlab("Input")+ylab("mC") %>% print()
hexbin_plot(input,mc,rpkm_vectors,pal,depths,40,50)+xlab("Input")+ylab("mC") %>% print()
dev.off()

pdf(file.path(figsdr,"Hexbin_plot_Input_vs_hmC.pdf"))
hexbin_plot(input,hmc,rpkm_vectors,pal,depths,10,15)+xlab("Input")+ylab("hmC") %>% print()
hexbin_plot(input,hmc,rpkm_vectors,pal,depths,15,20)+xlab("Input")+ylab("hmC") %>% print()
hexbin_plot(input,hmc,rpkm_vectors,pal,depths,40,50)+xlab("Input")+ylab("hmC") %>% print()
dev.off()
