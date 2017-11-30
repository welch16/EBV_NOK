
rm(list = ls())


library(tidyverse)
library(ggrepel)
library(tximport)
library(DESeq2)
library(RUVSeq)
library(BiocParallel,quietly = TRUE)

stri_subset <- function(string,pattern)
{
    string[!str_detect(string,pattern)]
}
    
sample_dir = "data/RSEM/hg19/Sept17"
samples_file = "data/Diff.Genes/hg19/Sept17/full_model/Sept17_Genes_samples_full.tsv"
contrast_file = "data/Diff.Genes/hg19/Sept17/full_model/Sept17_contrasts_full.tsv"
tpm_file = "data/TPM_matrices/Sept17/Genes_TPM_matrix_newBatch.tsv"


source("rfuns/full_model_gene_expression.R")

options(mc.cores = 8)

samples = samples_file %>% read_tsv()
contrasts = contrast_file %>% read_tsv()

sample_files = file.path(sample_dir,samples$Files) %>%
    set_names(get_reps(.))

## load files with tximport
txdata = tximport(sample_files,type = "rsem" , importer = read_tsv)

tpm_mat = read_tsv(tpm_file)

coldata = samples %>%
    select(-Files,-Replicate) %>%
    mutate(
        interac = paste(Treatment,Cell,sep = "_and_")
        ) %>% 
    as.data.frame()
rownames(coldata) = names(sample_files)

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


rowVars <- function (x,na.rm = TRUE)
{
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}

theme_set( theme_bw())

PCA_plot <- function(rl,vars = rownames(coldata),ntop = 500,x = 1,y=2)
{
    my_data = assay(rl)[,vars]
    rv = rowVars(my_data)

    select = order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]

    pca = prcomp(t(my_data[select,]))

    percentVar = pca$sdev^2/sum(pca$sdev^2)

    my_pca = colData(rl)[vars,,drop = FALSE] %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(
            sample = vars,
            rep = vars %>%
                strsplit("-") %>%
                map_chr( ~ .[length(.)]),
            x = pca$x[,x],
            y = pca$x[,y])

    my_pca %>%
        ggplot(aes(x,y,colour = interaction(Cell,Treatment)))+
        geom_point(size =2)+
        theme(legend.position = "top")+
        coord_fixed()+
        xlab(paste0("PC",x,": ",round(percentVar[x] * 100), "% variance"))+
        ylab(paste0("PC",y,": ",round(percentVar[y] * 100), "% variance"))+
        scale_color_brewer(palette = "Dark2",name = "Cell.Treatment")+
        geom_label_repel(aes(x,y,label = rep),size = 3,
                         box.padding = unit(.3,"lines"),
                         point.padding = unit(.7,"lines"),
                         show.legend = FALSE)


}

all_vars = rownames(coldata)
figsdr = "figs/PCA_MDS/Sept17"


deseq_full = DESEQ2_full_model(txdata,coldata,full_formula )

## message("Fitting main model...")
model_full = DESeq(deseq_full,minReplicatesForReplace = Inf,parallel = TRUE)

cd1 = colData(model_full)
residuals_full = residuals_DESeq2(model_full)

## this process RUVr
deseq_full_corr = DESEQ2_full_model_corr(model_full,
                                          residuals_full,
                                         paste(c(full_formula , "W_1"),collapse = "+"))

model_full_corr = DESeq(deseq_full_corr,
                   minReplicatesForReplace = Inf,parallel = TRUE)

cd2 = colData(model_full_corr)


rl = rlog(model_full_corr)


pdf(file.path(figsdr,"PCA_plots_allReplicates_afterRUVr.pdf"))
plots = list(
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 1,y= 2),
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 1,y= 3),
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 1,y= 4),
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 2,y= 3),
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 2,y= 4))
u = plots %>% print()
dev.off()

