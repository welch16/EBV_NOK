
rm(list = ls())

library(tidyverse)
library(parallel)
library(mosaics)

bindr = "data/bins/MeDIPseq"
binfiles = list.files(bindr,full.names = TRUE)

binfiles = binfiles %>% {.[grep("300",.)]}

medips_bins = read_tsv(binfiles[2])
input_bins = read_tsv(binfiles[1])

medips = medips_bins %>% names() %>% {.[!. %in% c("seqnames","start")]}
inputs = input_bins %>% names()  %>% {.[!. %in% c("seqnames","start")]}


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



file_relation = tibble(medips) %>%
    mutate(
        input = map_chr(medips,match_input,inputs),
        dip_depth = map_int(medips,get_depth,depths)        
    ) %>%
    filter(dip_depth > 1e7)

create_bins <- function(medip,
                        file_relation,
                        medips_bins,
                        input_bins)
{

    file_rel = file_relation %>%
        filter(medips == medip)


    new("BinData",
        chrID = medips_bins$seqnames,
        coord = as.integer(medips_bins$start),
        tagCount = medips_bins[[medip]],
        input = input_bins[[file_rel$input]],
        dataType = c("chip","input"),
        seqDepth = file_rel$dip_depth)
   
}

bins = map(file_relation$medips,
           create_bins,file_relation,medips_bins,input_bins)

fits = map(bins,mosaicsFit,analysisType = "IO",
           parallel = TRUE,nCore = 16)

peak_list = map(fits,mosaicsPeak,signalModel = "2S",FDR = .05,
                maxgap = 200,thres = 10)


export_peaks <- function(peak)
{
    peak@peakList %>%
        as.tbl() %>%
        mutate(
            peakID = paste0(chrID,":",peakStart,"-",peakStop) 
        ) %>%
        select(
            peakID,
            everything()
            )


}    

peaks = map(peak_list,
            export_peaks)

outdr = "data/peaks/hg19/MeDIPseq/fragLen300"

map2(peaks,file.path(outdr,paste0("Peaks_",file_relation$medips,".tsv")),write_tsv)



library(viridis)
library(scales)


hexbin_scatter_plot <- function(bin)
{
    counts = tibble(dip = bin@tagCount,
                    input = bin@input)

    counts = counts %>%
        mutate_all(funs( 1e9 * . / sum(.)))
       
    xvar = "dip"
    yvar = "input"
    

    pal = viridis(1e3, option = "D")

    counts %>%
        ggplot(aes_string(xvar,yvar))+stat_binhex(bins = 140) +
        scale_fill_gradientn(colours = pal,trans = "log10",
                             labels = trans_format('log10',math_format(10^.x)),
                             guide = guide_colourbar(title = NULL,
                                                     barheight = unit(0.92,"npc"),
                                                     barwidth = unit(0.01,"npc")))+
        geom_abline(slope = 1,intercept = 0,colour = "red",linetype = 2)+
        xlim(0,2e3)+ylim(0,2.5e3)
}


generate_plots <- function(fit,bin,file_rel)
{
    pdf(file.path(figsdr,paste0("diagnostic_peaks_",file_rel$medips,".pdf")))
    hexbin_scatter_plot(bin) + ylab(file_rel$medips) + xlab(file_rel$input) %>% print()
    plot(fit) 
    dev.off()
}    


figsdr = "figs/peaks"

map3(
    fits,
    bins,
    file_relation %>% {split(.,seq_len(nrow(.)))},
    generate_plots)



