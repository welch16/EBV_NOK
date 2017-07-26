
rm(list =ls())

library(tximport)
library(DESeq2)
library(tidyverse)

## get all files at gene / isoform level
indr = "./data/RSEM/hg19"

genefiles = indr %>%
    list.files(full.names = TRUE,
               pattern = "gene")

isofiles = indr %>%
    list.files(full.names = TRUE,
               pattern = "isof")

getNames <- function(files)
{
    files %>%
        basename() %>%
        strsplit("\\.") %>%
        map_chr( ~ .[1])
}    

getRep <- function(names)
{
    names %>%
        strsplit("-") %>%
        map( ~ .[grepl("rep",.) | grepl("clone",.)]) %>%
        unlist() %>%
        gsub("rep","",.) %>%
        gsub("clone","",.) %>%
        as.numeric()

}    

all_samples =
    tibble(
        name = isofiles %>% getNames(),
        isofiles,genefiles) %>%
    mutate(
        Strain = "NOKS",
        Strain = if_else(grepl("EBV",name) | grepl("akata",name),"EBV",Strain),
        Strain = if_else(grepl("clone",name),"EBV.dRdZ",Strain),
        Treatment = "none",
        Treatment = if_else(grepl("CaFBS",name),"CaFBS",Treatment),
        Treatment = if_else(grepl("methyl",name) | grepl("MC",name),"methyl",Treatment),
        Replicate = getRep(name),
        Lab = if_else(grepl("Noks",name),"Scott","Johannsen")
    ) %>%
    select(name,Lab,Strain,Treatment,Replicate,everything())


##
my_test =
    all_samples %>%
    filter(Treatment %in% c("methyl","none")) %>%
    filter(Strain %in% c("EBV","EBV.dRdZ")) %>%
    split(.$Lab)

my_test[["Johannsen"]] = my_test[["Johannsen"]] %>% filter(Strain == "EBV.dRdZ")

my_test = my_test %>%
    bind_rows()

genedata = tximport(my_test$genefiles,type = "rsem",importer = read_tsv)

isodata = tximport(my_test$isofiles,type = "rsem",importer = read_tsv)
                   
rsem_iso = my_test$isofiles[1] %>% read_tsv()

isodata[1:3] = isodata[1:3] %>%
    map( ~ {
            rownames(.) = rsem_iso$transcript_id
        .
    })



coldata = my_test %>%
    select(Strain,Treatment) %>%
    as.data.frame()
rownames(coldata) = my_test$name

library(BiocParallel)

options(mc.cores = 12)


deseq_analysis <- function(txdata,coldata)
{
    deseq =
        DESeqDataSetFromMatrix(floor(txdata[["counts"]]),
                               colData = coldata,
                               design= ~ Strain + Treatment)

    deseq = deseq[rowSums(counts(deseq)) > 1,]

    DESeq(deseq,minReplicatesForReplace = Inf,parallel = TRUE)

}    


gene_model = deseq_analysis(genedata,coldata)
iso_model = deseq_analysis(isodata,coldata)


build_RSEM_TPM <- function(files,var = "TPM")
{

    files %>%
        map(read_tsv,progress = FALSE) %>%
        map2(getNames(files) , ~ mutate(.x,set = .y)) %>%
        bind_rows() %>%
        select(contains("id"),set,contains(var)) %>%
        select(gene_id,everything()) %>%
        spread( set, TPM)

}    



tpm_gene = build_RSEM_TPM(my_test$genefiles)
tpm_iso = build_RSEM_TPM(my_test$isofiles)


get_summary_TPM <- function( tpm,fun)
{
    summary_col = tpm %>% 
        select(-contains("id")) %>%
        split(seq_len(nrow(.))) %>%
        map(unlist) %>%
        map_dbl(fun)                    # perhaps change this to parallel
    
    names(summary_col) = NULL
    summary_col
      
}                                    

tpm_gene = tpm_gene %>%
    mutate(
        maxTPM = get_summary_TPM(.,max),
        aveTPM = get_summary_TPM(.,mean)
        )

tpm_iso = tpm_iso %>%
    mutate(
        maxTPM = get_summary_TPM(.,max),
        aveTPM = get_summary_TPM(.,mean)
        )


res_gene = results(gene_model,
                   cooksCutoff = FALSE,
                   contrast = c("Strain","EBV","EBV.dRdZ"),
                   tidy = TRUE,
                   parallel = TRUE) %>% as_tibble() 

res_iso = results(iso_model,
                   cooksCutoff = FALSE,
                   contrast = c("Strain","EBV","EBV.dRdZ"),
                   tidy = TRUE,
                   parallel = TRUE) %>% as_tibble()

th = 1e-5

res_gene = res_gene %>%
    rename(gene_id = row) %>%
    left_join(
        tpm_gene %>% select(contains("id"),contains("TPM")),
        by = "gene_id") %>%
    mutate(
        DiffExp = if_else(padj <= th,"yes","no","no"),
        reg = if_else(log2FoldChange > 0,"up","down"),
        maxRank = dense_rank(desc(maxTPM)),
        aveRank = dense_rank(desc(aveTPM)))

res_iso = res_iso %>%
    rename(transcript_id = row) %>%
    left_join(
        tpm_iso %>% select(contains("id"),contains("TPM")),
        by = "transcript_id") %>%
    mutate(
        DiffExp = if_else(padj <= th,"yes","no","no"),
        reg = if_else(log2FoldChange > 0,"up","down"),
        maxRank = dense_rank(desc(maxTPM)),
        aveRank = dense_rank(desc(aveTPM)))



