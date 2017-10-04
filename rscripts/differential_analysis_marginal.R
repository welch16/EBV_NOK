#!/usr/bin/env Rscript

library(optparse)

optList = list(
    make_option("--no_treat_files", action = "store_true",type = "character",
                help = "RSEM output for samples where no treatment was applied"),
    make_option("--treat_files", action = "store_true", type = "character",
                help = "RSEM output for samples where the treatment was applied"),
    make_option("--treats",action = "store_true",type = "character",default = "no,yes",
                help = "Names of the two treatments separated by ','.
                        The default values are 'no,yes'"),
    make_option(c("-i","--iso"),action = "store_true",type = "character",default = "gene",
                help = "Flag indicating if the analysis is made at gene or isoform levels.
                        Either 'gene' or 'isoform'"),
    make_option("--tpm_file",action = "store_true",type = "character",
                help = "File with the TPM values for each gene / isoform and replicate"),
    make_option("--out_file", action = "store_true",default = tempfile(),type = "character",
                help = "Name of the file where the output is saved.")
)

opt = parse_args(OptionParser(option_list = optList))

library(base,quietly = TRUE)
library(tidyverse,quietly = TRUE)

separate_files <- function(ff)
{
  if(grepl(",",ff)){
      ff = ff %>% strsplit(',') %>% unlist     
  }else{
      ff = Sys.glob(ff)      
  }
  ff
}

library(tximport,quietly = TRUE)

get_reps <- function(x)
{
    x %>%
        basename() %>%
        map_chr( ~ gsub(paste0(".",tools::file_ext(.)), "",.)) %>%
        map_chr( ~ gsub(paste0(".",tools::file_ext(.)),"",.))      
    

}    

treats = opt$treats %>% strsplit(",") %>% unlist()

opt$no_treat_files = opt$no_treat_files %>%
    separate_files() %>%
    set_names(
        get_reps(.))

opt$treat_files = opt$treat_files %>%
    separate_files() %>%
    set_names(
        get_reps(.))



stopifnot(
    all(file.exists(opt$no_treat_files)),
    all(file.exists(opt$treat_files)),
    file.exists(opt$tpm_file),
    opt$iso %in% c("gene","isoform"))
       

opt$iso = opt$iso == "isoform"

all_files = c(opt$no_treat_files,opt$treat_files)

## load files with tximport
txdata = tximport(all_files,type = "rsem" , importer = read_tsv)

tpm_mat = read_tsv(opt$tpm_file)

if(opt$iso){                            # isoform level
    txdata[1:3] = txdata[1:3] %>%
    map( ~ {
            rownames(.) = tpm_mat$transcript_id
        .
    })
}

coldata = data.frame(
    treat = c(
        rep_len(treats[1],length(opt$no_treat_files)),
        rep_len(treats[2],length(opt$treat_files))))

rownames(coldata) = names(all_files)


library(DESeq2,quietly = TRUE)
library(BiocParallel,quietly = TRUE)

options(mc.cores = 20)



deseq_analysis <- function(txdata,coldata)
{
    deseq =
        DESeqDataSetFromMatrix(floor(txdata[["counts"]]),
                               colData = coldata,
                               design= ~  treat)

    deseq[rowSums(counts(deseq)) > 1,]

}    

deseq = deseq_analysis(txdata,coldata)

message("Fitting main model...")
model = DESeq(deseq,minReplicatesForReplace = Inf,parallel = TRUE)

## the test is: yes vs no
results = results(model,
                  cooksCutoff = FALSE,tidy = TRUE) %>%
    as_tibble()


get_mean_tpm <- function(col,samples,tpm_mat,coldata,
                         varname = paste("rank_mean_tpm",samples,sep = ":"))
{
    var = coldata %>%
        mutate(file = rownames(.)) %>% 
        filter("["(.,col) == samples ) %>%
        select(file) %>% .[[1]]

    tpm_mat %>%
        select(contains("id")) %>%
        mutate(
            !!varname  := tpm_mat[,var] %>%
                rowMeans()
            )

}

message("Adding ranking variables...")

tpm_mat = bind_cols(
    tpm_mat[,1:2],
    tpm_mat[,names(all_files)])

no_treat_tpm = get_mean_tpm("treat","no",tpm_mat,coldata)
treat_tpm = get_mean_tpm("treat","yes",tpm_mat,coldata)

by_ids = intersect(
    no_treat_tpm %>% names(),
    treat_tpm %>% names())

rank_tpm = inner_join(no_treat_tpm,
                      treat_tpm,by = by_ids)

if(opt$iso){
    results = results %>%
        dplyr::rename(transcript_id = row) %>%
        inner_join(rank_tpm, by = "transcript_id")
}else{
    results = results %>%
        dplyr::rename(gene_id = row) %>%
        inner_join(rank_tpm, by ="gene_id")
}


## desc means that if the number in rank_mean_tpm is
## low then there is high expression

results = results %>%
    mutate_at(
        .vars = vars(contains("rank_mean_tpm")),
        .funs = funs(
            dense_rank(desc(.))
        ) ) %>% 
    mutate(
        baseMean = NULL) %>%
    dplyr::rename(
               log2FC = log2FoldChange,
               log2FC_se = lfcSE) %>%
    select(contains("id"),everything())

if(opt$iso){
    results = results %>%        
        select(transcript_id,everything())
}else{
    results = results %>%
        select(gene_id,everything())
}


message("Saving results...")
write_tsv( results,
          opt$out_file)


