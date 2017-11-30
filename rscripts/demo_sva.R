
rm(list = ls())

library(edge)
library(tidyverse)
library(sva)

library(tximport)

rsemdr = "data/RSEM/hg19/Sept17"


stri_subset <- function(string,pattern)
{

    string[ ! str_detect(string,pattern)]

}
       

rsem_files = rsemdr %>%
    list.files(full.names = TRUE) %>%
    str_subset("genes") %>%
    stri_subset("clone2") %>%
    stri_subset("clone4") %>% 
    set_names(
        basename(.) %>%
        str_replace(".genes.results",""))

rsem = rsem_files %>%
    map(read_tsv)  %>%
    map(separate,gene_id,c("ensembl_id","gene_id"),sep = "_")

build_matrix <- function(rsem_list,var = "TPM")
{
    cnames = names(rsem_list)
    rnames = rsem_list %>%
        pluck(1) %>%
        select(ensembl_id) %>%
        unlist() %>%
        set_names(NULL)

    mat = rsem_list %>%
        map(select_,.dots = "TPM") %>%
        bind_cols() %>%
        as.matrix()

    colnames(mat) = cnames
    rownames(mat) = rnames

    mat
}

mat = rsem %>% build_matrix("FPKM") 

cov = tibble(
    file = names(rsem_files)) %>%
    mutate(
        cell = case_when(
            str_detect(file,"clone") ~ "EBV.dRdZ",
            str_detect(file,"akata") ~ "EBV",
            TRUE ~ "NOKS"),
        treatment = if_else(str_detect(file,"methyl"),"methyl","no_treatmnet"),
        file = NULL) %>%
    as.data.frame()
rownames(cov) = colnames(mat)

de_obj = build_models(
    data = mat,
    cov = cov,
    full.model = ~ cell + treatment,
    null.model = ~ treatment )

fit = fit_models(de_obj,stat.type = "lrt")

## works until here

de_op = odp(de_obj,fit)

##  when trying to run this part, it crashes by some subscript out of
##  bound or something

## > de_op = odp(de_obj,fit, n.mods = 10)
##  Null iteration:  100
## Error in smooth.spline(lambda, pi0, df = smooth.df) : 
##   missing or infinite values in inputs are not allowed
## > ?odp
## > de_op = odp(de_obj,fit)
## Error in `[<-`(`*tmp*`, l, , value = colMeans(de.fit@fit.full[mod.member ==  : 
##   subscript out of bounds

## in the vignette, it wasn't quite clear what gene expression
## measurement was used. Tried TPM, expected_counts and FPKM, it
## crashed with the same errors every run.


de_sva = apply_sva(de_obj, B = 10)




## quick note, this is a very simple exercise and it is crashing.
## tried without success:

## - changing the number of surrogate variable
## - changing the number of bootstrap samples  

