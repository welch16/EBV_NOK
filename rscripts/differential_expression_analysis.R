#!/usr/bin/env Rscript

library(optparse)

opt_list = list(
    make_option("--sample_dir",action = "store_true",type = "character",default = "./",
                help = "Directory where the RSEM results are stored"),
    make_option("--samples_file",action = "store_true",type = "character",
                help = "File with the samples to used in the analysis and the
                        design of the experiment. The first column is always the files and
                        the last one are the relicates"),
    make_option("--contrast_file",action = "store_true",type = "character",
                help = "File with the contrasts to be calculated after the models is fitted"),
    make_option("--tpm_file",action = "store_true",type = "character",
                help = "File with the TPM matrix of all samples"),
    make_option(c("-i","--iso"),action = "store_true",type = "character",default = "gene",
                help = "Flag indicating if is genes or isoform"),
    make_option("--figs_dir",action = "store_true",type = "character",default = tempdir(),
                help = "Directory where the summary plots are stored"),
    make_option("--cores",action = "store_true",type = "numeric",default = 20,
                help = "Number of parallel workers to use"),
    make_option("--out_dir",action = "store_true",type = "character",default = tempdir(),
                help = "Directory where the contrast are going to be saved")
)

opt = parse_args(OptionParser(option_list = opt_list))

opt$sample_dir = "data/RSEM/hg19/Sept17"
opt$samples_file = "data/Diff.Genes/hg19/Sept17/full_model/Sept17_Genes_samples_full.tsv"
opt$contrast_file = "data/Diff.Genes/hg19/Sept17/full_model/Sept17_contrasts_full.tsv"
opt$tpm_file = "data/TPM_matrices/Sept17/Genes_TPM_matrix_newBatch.tsv"

library(base,quietly = TRUE)
library(tidyverse,quietly = TRUE)
library(tximport,quietly = TRUE)
library(scales)

source("rfuns/full_model_gene_expression.R")

options(mc.cores = opt$cores)

stopifnot(file.exists(opt$samples_file),
          file.exists(opt$contrast_file),
          file.exists(opt$tpm_file),
          opt$iso %in% c("gene","iso","isoform"))

samples = opt$samples_file %>%
    read_tsv()

contrasts = opt$contrast_file %>%
    read_tsv()

sample_files = file.path(opt$sample_dir,samples$Files) %>%
    set_names(get_reps(.))

stopifnot(all(file.exists(sample_files)))

opt$iso = opt$iso == "isoform"


## load files with tximport
txdata = tximport(sample_files,type = "rsem" , importer = read_tsv)

tpm_mat = read_tsv(opt$tpm_file)

if(opt$iso){                            # isoform level
    txdata[1:3] = txdata[1:3] %>%
    map( ~ {
            rownames(.) = tpm_mat$transcript_id
        .
    })
}

coldata = samples %>%
    select(-Files,-Replicate) %>%
    mutate(
        interac = paste(Treatment,Cell,sep = "_and_")
        ) %>% 
    as.data.frame()
rownames(coldata) = names(sample_files)


library(DESeq2,quietly = TRUE)
library(BiocParallel,quietly = TRUE)

factors = samples %>%
    select( -c(1,ncol(samples))) %>%
    names()

if(length(factors) == 2){
    full_formula = paste0("~",
                          paste(c(factors,
                                  paste(factors,collapse = ":")),
                                collapse= "+"))
}else{
    stop("need to figure out this ammount of factors")
}


deseq_full = DESEQ2_full_model(txdata,coldata,full_formula )
deseq_interac = DESEQ2_full_model(txdata,coldata,"~interac")

message("Fitting main model...")
model_full = DESeq(deseq_full,minReplicatesForReplace = Inf,parallel = TRUE)
model_interac = DESeq(deseq_interac,minReplicatesForReplace = Inf , parallel = TRUE)

stopifnot(
    identical(
        colData(model_full)$sizeFactor,
        colData(model_interac)$sizeFactor))

message("Performing contrasts...")
 ## Note: The ranking is decreasing, hence highest expressed genes
 ## will have low ranked values
results = contrasts %>%
    split(.$Contrast_name) %>%
    map(evaluate_contrast,
        model_full,
        model_interac,tpm_mat,opt$iso)

theme_set(theme_bw())

pdf(file.path(opt$figs_dir,paste0("Pvalue_histograms_",if_else(opt$iso,"Isoforms","Genes"),".pdf")))
u = results %>%
    map2(names(.),pvalue_histogram) %>%
    map(print)
dev.off()


pdf(file.path(opt$figs_dir,paste0("MA_plots_",if_else(opt$iso,"Isoforms","Genes"),".pdf")))
u = results %>%
    map2(names(.),MA_plot) %>%
    map(print)
dev.off()

pdf(file.path(opt$figs_dir,paste0("Volcano_plots_",if_else(opt$iso,"Isoforms","Genes"),".pdf")))
u = results %>%
    map2(names(.),volcano_plot) %>%
    map(print)
dev.off()

out_files = file.path(opt$out_dir,
                      paste0(if_else(opt$iso,"Isoforms","Genes"),"_",
                             names(results),".tsv"))

message("Saving results...")
map2(results,
     out_files,
     write_tsv)

