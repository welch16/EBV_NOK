
rm(list = ls())

library(GenomicAlignments)
library(magrittr)
        
indr = "data/Diff_Methyl"

load(file.path(indr,"Regions_DiffMethyl.RData")) # regions

readdr = "data/BAM/hg19/bowtie_dip"
files = list.files(readdr,full.names = TRUE) %>%
    {.[grep("sort.bam",.)]} %>%
    {.[grep("bai",.,invert = TRUE)]} %>%
    {.[grep("rep0",.,invert = TRUE)]} %>%
    {.[grep("Input",.,invert = TRUE)]}


library(dplyr)

depths = readr::read_tsv("data/metadata/MeDIPseq_sequencingDepth.tsv")
depths = depths %>%
    filter(!grepl("Input",file)) %>%
    filter(depth > 1e7)

files = files[basename(files) %in% depths$file]

mC_files = files[!grepl("hmC",files)]
hmC_files = files[grepl("hmC",files)]

generate_count_matrix <- function(files,regions,border,fraglen)
{

    regions_border = regions
    start(regions_border) = start(regions) - border
    end(regions_border) = end(regions) + border

    reads = lapply(files,
                   readGAlignments,
                   param = ScanBamParam(which = regions_border))
    names(reads) = gsub(".sort.bam","",basename(files))

    reads = lapply(reads,as,"GRanges")
    reads = lapply(reads,resize,fraglen)

    mcols(regions) = lapply(reads,function(x){
        countOverlaps(regions,x) %>% DataFrame()
        })
                          

    regions

}    

regions_counts = list()
regions_counts[["mC"]] = generate_count_matrix(mC_files,regions[["mC"]],500,300)
regions_counts[["hmC"]] = generate_count_matrix(hmC_files,regions[["hmC"]],500,300)

save( regions_counts, file = file.path(indr,"Regions_DiffMethyl_counts.RData"))





