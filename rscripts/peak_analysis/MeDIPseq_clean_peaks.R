
rm(list = ls())

library(tidyverse)
library(parallel)

bindr = "data/bins/MeDIPseq"
binfiles = list.files(bindr,full.names = TRUE)

binfiles = binfiles %>% {.[grep("300",.)]}

medips_bins = read_tsv(binfiles[2],n_max = 3)
input_bins = read_tsv(binfiles[1],n_max = 3)

medips = medips_bins %>% names() %>% {.[!. %in% c("seqnames","start")]}
inputs = input_bins %>% names()  %>% {.[!. %in% c("seqnames","start")]}

rm(medips_bins,input_bins)

match_input <- function(medip,inputs)
{

    chars = medip %>%
        {strsplit(.,"_")} %>%
        {.[[1]]}
    inputs %>%
        {.[grepl(chars[2],.)]} %>%
        {.[grepl(chars[3],.)]}
               
}   

depths = read_tsv("data/metadata/MeDIPseq_sequencingDepth.tsv")
input_depths = depths %>% filter(grepl("Input",file))

depths = depths %>%
    filter(!grepl("Input",file)) 

get_depth <- function(medip,
                      depths)
{
    chars = medip %>%
        {strsplit(.,"_")} %>%
        {.[[1]]}

    if(chars[2] == "EBV"){
        depths = depths %>% filter(grepl("akata",file))
    }else{
        depths = depths %>% filter(!grepl("akata",file))
    }

    depths = depths %>%
        filter(grepl(ifelse(chars[3] == "NoTrt","mono","CaFBS"),file))


    if(chars[4] == "hmC"){
        depths = depths %>%
            filter(grepl("hmC",file))
    }else{
        depths = depths %>%
            filter(!grepl("hmC",file))
    }

    depths = depths %>%
        filter(grepl(tolower(chars[5]),file))

    depths$depth
        

}    



get_input_depth <- function(input,depths)
{

    depths = depths %>%
        filter(grepl("rep0",file))

    chars = input %>%
        {strsplit(.,"_")} %>%
        {.[[1]]}

    if(chars[2] == "EBV"){
        depths = depths %>% filter(grepl("akata",file))
    }else{
        depths = depths %>% filter(!grepl("akata",file))
    }

    depths = depths %>%
        filter(grepl(ifelse(chars[3] == "NoTrt","mono","CaFBS"),file))


    depths = depths %>%
        filter(grepl(tolower(chars[5]),file))

    depths$depth
        
}    


file_relation = tibble(medips) %>%
    mutate(
        input = map_chr(medips,match_input,inputs),
        dip_depth = map_int(medips,get_depth,depths),
        input_depth = map_int(input,get_input_depth,input_depths)
    )

peakdir = "data/peaks/hg19/MeDIPseq/fragLen300"
peakfiles = list.files(peakdir,full.names = TRUE)

peaks = peakfiles %>%
    map(read_tsv)

names(peaks) = peakfiles %>%
    basename() %>%
    {gsub("Peaks_","",.)} %>%
    {gsub(".tsv","",.)}

peak_summary = tibble(medips = names(peaks) , npeaks = peaks %>% map_int(nrow)) %>%
    left_join(file_relation ,by = "medips") %>%
    select(-contains("nput")) 
    
## peak_summary %>% mutate_if(is.numeric,funs(prettyNum(.,big.mark=","))) %>% as.data.frame()
##                         medips  npeaks  dip_depth
## 1  MeDIPseq_EBV_CaFBS_hmC_Rep1  37,996 34,316,338
## 2  MeDIPseq_EBV_CaFBS_hmC_Rep2 131,394 23,663,657
## 3  MeDIPseq_EBV_CaFBS_hmC_Rep3  87,314 27,639,538
## 4   MeDIPseq_EBV_CaFBS_mC_Rep1 289,003 21,247,394
## 5   MeDIPseq_EBV_CaFBS_mC_Rep2 277,002 26,243,419
## 6   MeDIPseq_EBV_CaFBS_mC_Rep3 196,920 24,132,875
## 7  MeDIPseq_EBV_NoTrt_hmC_Rep1  49,531 32,231,395
## 8  MeDIPseq_EBV_NoTrt_hmC_Rep2 250,295 27,179,844
## 9  MeDIPseq_EBV_NoTrt_hmC_Rep3  74,984 16,502,496
## 10  MeDIPseq_EBV_NoTrt_mC_Rep1  88,994 21,058,739
## 11  MeDIPseq_EBV_NoTrt_mC_Rep3 204,490 19,931,469
## 12 MeDIPseq_NOK_CaFBS_hmC_Rep1  47,963 51,556,777
## 13 MeDIPseq_NOK_CaFBS_hmC_Rep2 171,874 25,713,480
## 14 MeDIPseq_NOK_CaFBS_hmC_Rep3 143,982 21,986,186
## 15  MeDIPseq_NOK_CaFBS_mC_Rep1 170,736 21,570,297
## 16  MeDIPseq_NOK_CaFBS_mC_Rep2 150,315 14,503,924
## 17  MeDIPseq_NOK_CaFBS_mC_Rep3 247,377 20,330,070
## 18 MeDIPseq_NOK_NoTrt_hmC_Rep1  46,190 31,291,638
## 19 MeDIPseq_NOK_NoTrt_hmC_Rep2 204,439 34,397,949
## 20 MeDIPseq_NOK_NoTrt_hmC_Rep3 160,621 22,184,075
## 21  MeDIPseq_NOK_NoTrt_mC_Rep1 248,447 19,041,385

    

fix_details <- function(peak,name,file_relation)
{
    file_rel = file_relation %>%
        filter(medips == name)

    dip_depth = file_rel$dip_depth
    input_depth = file_rel$input_depth

    peak %>%
        mutate(
            aveInputCountScaled = (dip_depth / input_depth * aveInputCount ),            
            aveLog2Ratio = log2( ( 1 + aveChipCount ) / (1 + aveInputCountScaled))
        )
    

}    

peaks = peaks %>% map2(names(peaks),
                       fix_details,file_relation)

## summary plots:

## aveLog2Ratio                       
## min expectation, (much) greater than one for the majority of peaks and
## al

file_rel = file_relation %>% select(medips,dip_depth) %>%
    separate(medips,into = c("extra","Cell","Treat","Meth","Rep"),
             sep = "\\_",remove = FALSE)

all_together = peaks %>% map2(names(peaks),
                              .f = function(x,y)
                                  x %>%
                                  mutate(medips = y)) %>%
    bind_rows() %>%
    left_join(file_rel,by = "medips")


figsdir = "figs/peaks"

theme_set(theme_bw())


pdf(file.path(figsdir,"aveLog2Ratio_across_samples_FDR5.pdf"),width = 9)
all_together %>%
    ggplot(aes(medips,aveLog2Ratio,fill = Cell))+geom_boxplot(outlier.shape = NA)+
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank()        
    )+
    geom_abline(slope = 0,intercept = 1,linetype = 2,colour = "red")+
    facet_grid( ~ Meth,scales = "free_x")+
    scale_fill_brewer(palette = "Pastel1")+
    ylim(0,4)
dev.off()

pdf(file.path(figsdir,"aveLogP_across_samples_FDR5.pdf"),width = 9)
all_together %>%
    ggplot(aes(medips,logAveP,fill = Cell))+geom_boxplot(outlier.shape = NA)+
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank()        
    )+
    geom_abline(slope = 0,intercept = 0,linetype = 2,colour = "red")+
    facet_grid( ~ Meth,scales = "free_x")+
    scale_fill_brewer(palette = "Pastel1")+
    ylim(0,8)    
dev.off()

pdf(file.path(figsdir,"Depth_vs_medianAveLog2Ratio_FDR5.pdf"),width = 9,height = 5)
all_together %>% group_by(medips,
                          Cell,
                          Treat,
                          Meth,
                          Rep,
                          dip_depth
                          ) %>%
    summarize(
        medlog2 = median(aveLog2Ratio)
    ) %>%
    ggplot(aes(dip_depth,
               medlog2,colour = Cell))+geom_point(size = 3)+
    facet_grid(~ Meth,scales = "free_x")+
    scale_colour_brewer(palette = "Set1")+
    theme(
        legend.position = "bottom"
    )+
    xlab("Sequencing Depth")+
    ylab("Median average Log2 Ratio")
dev.off()

              
library(Vennerable)
library(GenomicRanges)
library(ChIPseeker)

groups = file_rel %>%
    filter(dip_depth > 1e7) %>%
    group_by(Cell,Treat,Meth) %>%
    summarize(
        medips = paste(medips,collapse = ","),
        n = n()
    ) %>%
    filter(n > 1) %>%
    split(seq_len(nrow(.)))


generate_venn_weights <- function( group,peaks,topK = Inf)
{
    sets = group$medips %>%
        {strsplit(.,",")} %>%
        .[[1]]

    peaks = peaks[sets]

    if(is.infinite(topK)){
        gr = peaks
    }else{
        gr = peaks %>% map(top_n,topK,aveLog2Ratio)
    }
       
    gr = gr %>% map(.f = function(x)
        GRanges(
            seqnames = x$chrID,
            ranges = IRanges(
                start = x$peakStart,
                end = x$peakStop)))
    names(gr) = NULL

    reduced = reduce(do.call(c,gr))
    
    names(gr) = sets
    
    
    mcols(reduced) = gr %>% 
        map(.f = function(x)
            as.integer(
                ifelse(countOverlaps(reduced,x) > 0, 1, 0))) %>%
        DataFrame()

    reduced %>%
        as.data.frame() %>%
        as.tbl() %>%
        mutate(
            seqnames = as.character(seqnames),
            strand = NULL,
            peakID = paste0(seqnames,":",start,"-",end)
        ) %>%
        select(peakID,everything())

    
}   

all_peaks = groups %>% map(generate_venn_weights,peaks)

groups_summary = groups %>% bind_rows() %>%
    mutate(
        name = paste(Cell,Treat,Meth,sep = "_")
        ) %>% ungroup()

names(all_peaks) = groups_summary$name


topK = 35e3
##

## topK_peaks = groups %>% map(generate_venn_weights,peaks,topK)
## names(topK_peaks) = groups_summary$name

## ## Venn diagrams

## peak_venn_plot <- function(peaks,name,topK)
## {
##     ids = peaks %>% select(peakID) %>% .[[1]]
##     peaks = peaks %>% select(contains("MeDIP")) %>%
##         as.list()

##     peaks = peaks %>% map(.f = function(x)ids[x == 1])

##     pdf(file.path(figsdir,paste0(name,"_VennDiagram_",topK,".pdf")))
##     vennplot(peaks,by = "Vennerable")
##     dev.off()

## }    

## topK_peaks %>% map2(names(topK_peaks),peak_venn_plot,topK)


diff_groups = file_rel %>%
    filter(dip_depth > 1e7) %>%
    group_by(Cell,Meth) %>%
    summarize(
        medips = paste(medips,collapse = ","),
        n = n()
    ) %>%
    filter(n > 1) %>%
    split(seq_len(nrow(.)))


topK_peaks = diff_groups %>% map(generate_venn_weights,peaks,topK)

groups_summary = diff_groups %>% bind_rows() %>%
    mutate(
        name = paste(Cell,Meth,sep = "_")
        ) %>% ungroup()


names(topK_peaks) = groups_summary$name


## take intersection

filter_intersection <- function(peaks)
{

    which.peaks = peaks %>% select(contains("MeDI")) %>%
        mutate(
            nsets = rowSums(.)
        ) %>%
        select(nsets)

    ncol = peaks %>% select(contains("MeDI")) %>% ncol()
    peaks %>% filter(which.peaks$nsets == ncol)
       

}

diff_regions = topK_peaks %>% map(filter_intersection)


diff_regions = diff_regions %>%
    map(.f = function(x)
        GRanges(seqnames = x$seqnames,
                ranges = IRanges(
                    start = x$start,
                    end = x$end)))


regions = list()
regions[["mC"]] = reduce(c(diff_regions$EBV_mC,diff_regions$NOK_mC))
regions[["hmC"]] = reduce(c(diff_regions$EBV_hmC,diff_regions$NOK_hmC))              

outdr = "data/Diff_Methyl"
dir.create(outdr)

save(regions,file = file.path(outdr,"Regions_DiffMethyl.RData"))

