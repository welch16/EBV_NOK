#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--A_noTr", action = "store_true",type = "character",
                help = "Files corresponding to A cell line and no treatment was applied"),
    make_option("--B_noTr", action = "store_true", type = "character",
                help = "Files corresponding to B cell line and no treatment was applied"),
    make_option("--A_Tr", action = "store_true", type = "character",
                help = "Files corresponding to A cell line and the treatment was applied"),
    make_option("--B_Tr", action = "store_true", type = "character",
                help = "Files corresponding to B cell line and the treatment was applied"),
    make_option("--cells",action = "store_true",default = "A,B",type = "character",
                help = "Names of the two cell lines separated by ','.
                        The default value is 'A,B', the test that is going to perform are of
                        the type B vs A."),
    make_option("--treats",action = "store_true",type = "character",default = "no,yes",
                help = "Names of the two treatments separated by ','.
                        The default values are 'no,yes'"),
    make_option(c("-i","--iso"),action = "store_true",type = "character",default = "gene",
                help = "Flag indicating if the analysis is made at gene or isoform levels.
                        Either 'gene' or 'isoform'"),
    make_option("--tpm_file",action = "store_true",type = "character",
                default = "data/TPM_matrices/Genes_TPM_matrix.tsv",
                help = "File with the TPM values for each gene / isoform and replicate"),
    make_option("--out_file_suff", action = "store_true",default = tempfile(),type = "character",
                help = "Suffix of the files where the test results are saved.")
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
 
opt$A_noTr = separate_files(opt$A_noTr)
opt$B_noTr = separate_files(opt$B_noTr)
opt$A_Tr = separate_files(opt$A_Tr)
opt$B_Tr = separate_files(opt$B_Tr)

library(tximport,quietly = TRUE)

source("rfuns/geneExpression_analysis.R")

cells = opt$cells %>% strsplit(",") %>% unlist
treats = opt$treats %>% strsplit(",") %>% unlist

names(opt$A_noTr) = getRep(opt$A_noTr,paste(cells[1],treats[1],sep = "_"))
names(opt$B_noTr) = getRep(opt$B_noTr,paste(cells[2],treats[1],sep = "_"))
names(opt$A_Tr) = getRep(opt$A_Tr,paste(cells[1],treats[2],sep = "_"))
names(opt$B_Tr) = getRep(opt$B_Tr,paste(cells[2],treats[2],sep = "_"))

stopifnot(all(file.exists(opt$A_noTr)),
          all(file.exists(opt$B_noTr)),
          all(file.exists(opt$A_Tr)),
          all(file.exists(opt$B_Tr)),
          opt$iso %in% c("gene","isoform"))

opt$iso = opt$iso == "isoform"

all_files = c(opt$A_noTr,opt$B_noTr,opt$A_Tr,opt$B_Tr)

## load files with tximport
txdata = tximport(all_files,
                  type = "rsem" , importer = read_tsv)

tpm_mat = read_tsv(opt$tpm_file)

if(opt$iso){                            # isoform level
    txdata[1:3] = txdata[1:3] %>%
    map( ~ {
            rownames(.) = tpm_mat$transcript_id
        .
    })
}

coldata = data.frame(
    cell = c(rep(cells[1],length(opt$A_noTr)),
             rep(cells[2],length(opt$B_noTr)),
             rep(cells[1],length(opt$A_Tr)),
             rep(cells[2],length(opt$B_Tr))),
    treat = c(rep(treats[1],length(opt$A_noTr)),
             rep(treats[1],length(opt$B_noTr)),
             rep(treats[2],length(opt$A_Tr)),
             rep(treats[2],length(opt$B_Tr))))

coldata = coldata %>%
    mutate(
        interac = paste(cell,treat,sep = "_"),
        cell = factor(cell,levels = cells),
        treat = factor(treat,levels = treats),
        interac = factor(interac,
                         levels =
                             paste(
                                 rep( cells, 2),
                                 rep(treats,each = 2),
                                 sep = "_")))
        

rownames(coldata) = names(all_files)


library(DESeq2,quietly = TRUE)
library(BiocParallel,quietly = TRUE)

options(mc.cores = 20)

deseq_analysis <- function(txdata,coldata)
{
    deseq =
        DESeqDataSetFromMatrix(floor(txdata[["counts"]]),
                               colData = coldata,
                               design= ~ cell + treat)

    deseq[rowSums(counts(deseq)) > 1,]

}    

deseq = deseq_analysis(txdata,coldata)

message("Fitting main model...")
model_1 = DESeq(deseq,minReplicatesForReplace = Inf,parallel = TRUE)

design(deseq) <- ~ interac

message("Fitting interactions model ...")
model_2 = DESeq(deseq,minReplicatesForReplace = Inf ,parallel = TRUE)

test = paste(cells,collapse = "_vs_")

message("Contrasting hypothesis...")
test_results = list()
test_results[[paste(rev(cells),collapse = "_vs_")]] =
    results(model_1,
            cooksCutoff = FALSE,
            contrast = c("cell",rev(cells)),tidy = TRUE) %>%
    as_tibble()
test_results[[paste(treats[1],
                    paste(rev(cells),collapse = "_vs_"),sep = ":")]] =
    results(model_2,
            cooksCutoff = FALSE,
            contrast = c("interac",
                         paste(rev(cells), treats[1],sep = "_")),tidy = TRUE) %>%
    as_tibble()
test_results[[paste(treats[1],
                    paste(rev(cells),collapse = "_vs_"),sep = ":")]] =
    results(model_2,
            cooksCutoff = FALSE,
            contrast = c("interac",
                         paste(rev(cells), treats[1],sep = "_")),tidy = TRUE) %>%
    as_tibble()
test_results[[paste(treats[2],
                    paste(rev(cells),collapse = "_vs_"),sep = ":")]] =
    results(model_2,
            cooksCutoff = FALSE,
            contrast = c("interac",
                         paste(rev(cells), treats[2],sep = "_")),tidy = TRUE) %>%
    as_tibble()


coldata = coldata %>%
    mutate(
        cols = all_files %>%
            basename() %>%
            strsplit("\\.") %>%
            map_chr( ~.[1]))

get_mean_tpm <- function(col,samples,tpm_mat,coldata,
                         varname = paste("rank_mean_tpm",samples,sep = ":"))
{
    var = coldata %>%
        filter("["(.,col) == samples ) %>%
        select(cols) %>% .[[1]]

    tpm_mat %>%
        select(contains("id")) %>%
        mutate(
            !!varname  := tpm_mat[,var] %>%
                rowMeans()
            )

}

message("Adding ranking variables...")
summary_B = c(cells[2],paste(cells[2],treats,sep = "_"))
summary_A = c(cells[1],paste(cells[1],treats,sep = "_"))

mean_tpm_B = map2(
    c("cell",rep("interac",2)),
    summary_B,
    get_mean_tpm,
    tpm_mat,coldata)
names(mean_tpm_B) = summary_B

mean_tpm_A = map2(
    c("cell",rep("interac",2)),
    summary_A,
    get_mean_tpm,
    tpm_mat,coldata)
names(mean_tpm_A) = summary_A

clean_results <- function(res,summ_B,summ_A,iso )
{
    nms_A = summ_A %>%
        select(-contains("id")) %>%
        names()

    nms_B = summ_B %>%
        select(-contains("id")) %>%
        names()

    summ = inner_join(summ_B,summ_A,
                      by = summ_A %>%
                          select(contains("id")) %>%
                          names())    
    if(iso){
        
        out = res %>%
            dplyr::rename(transcript_id = row) %>%
            left_join(summ, by = "transcript_id") %>%
            select(contains("id"),everything()) 
            

    }else{

        out = res %>%
            dplyr::rename(gene_id = row) %>%
            left_join(summ,by = "gene_id") %>%
            select(contains("id"),everything())
   
    }

    out %>%
        mutate_at(
            .vars = vars(nms_B,nms_A),
            .funs = funs(
                dense_rank(desc(.))
            ) )  %>%
        mutate(
            baseMean = NULL) %>%
        dplyr::rename(
            log2FC = log2FoldChange,
            log2FC_se = lfcSE)
            
  
}    

out_results = pmap(
    list(
        test_results,
        mean_tpm_B,
        mean_tpm_A),
     clean_results,opt$iso)

message("Saving results...")
out_results %>%
    map2( names(.),
         ~ write_tsv( .x , path = paste0(opt$out_file_suff,.y,".tsv")))


