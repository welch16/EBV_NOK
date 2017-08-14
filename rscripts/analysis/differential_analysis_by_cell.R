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
    make_option("--iso",action = "store_true",type = "logical",default = FALSE,
                help = "Flag indicating if the analysis is made at gene or isoform levels"),
    make_option("--tpm_file",action = "store_true",type = "character",
                default = "data/TPM_matrices/Genes_TPM_matrix.tsv",
                help = "File with the TPM values for each gene / isoform and replicate"),
    make_option("--outfile", action = "store_true",default = tempfile(),type = "character",
                help = "Name of the outfile where the test results are saved.")
)

opt = parse_args(OptionParser(option_list = optList))

opt$A_noTr = "data/RSEM/hg19/RNAseq-Noks-mono-rep?.genes.results"
opt$B_noTr = "data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep2.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep3.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep4.genes.results"
opt$A_Tr = "data/RSEM/hg19/RNAseq-Noks-MC-rep?.genes.results"
opt$B_Tr = "data/RSEM/hg19/RNAseq-Noks_EBV-MC-rep?.genes.results"  

opt$cells = "NOK,EBV"
opt$treatments = "none,MC"  

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
treats = opt$treatment %>% strsplit(",") %>% unlist

names(opt$A_noTr) = getRep(opt$A_noTr,paste(cells[1],treats[1],sep = "_"))
names(opt$B_noTr) = getRep(opt$B_noTr,paste(cells[2],treats[1],sep = "_"))
names(opt$A_Tr) = getRep(opt$A_Tr,paste(cells[1],treats[2],sep = "_"))
names(opt$B_Tr) = getRep(opt$B_Tr,paste(cells[2],treats[2],sep = "_"))

stopifnot(all(file.exists(opt$A_noTr)),
          all(file.exists(opt$B_noTr)),
          all(file.exists(opt$A_Tr)),
          all(file.exists(opt$B_Tr)))

all_files = c(opt$A_noTr,opt$B_noTr,opt$A_Tr,opt$B_Tr)

## load files with tximport
txdata = tximport(all_files,
                  type = "rsem" , importer = read_tsv)

tpm_mat = read_tsv(opt$tpm_file)

if(opt$iso){                            # isoform level
    txtdata[1:3] = txdata[1:3] %>%
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

options(mc.cores = 12)

deseq_analysis <- function(txdata,coldata)
{
    deseq =
        DESeqDataSetFromMatrix(floor(txdata[["counts"]]),
                               colData = coldata,
                               design= ~ cell + treat)

    deseq[rowSums(counts(deseq)) > 1,]

}    

deseq = deseq_analysis(txdata,coldata)

model_1 = DESeq(deseq,minReplicatesForReplace = Inf,parallel = TRUE)

design(deseq) <- ~ interac

model_2 = DESeq(deseq,minReplicatesForReplace = Inf ,parallel = TRUE)

test = paste(cells,collapse = "_vs_")

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

get_ave_tpm_rank <- function(tpm_mat,col,samples,coldata)
{
    browser()


}

#%>% strsplit("_") %>% sapply(function(x)x[2])

clean_results <- function(res, genes)
{   
    res %>% as.data.frame %>%
        as.tbl %>% mutate(gene = genes) %>%
        select(gene,everything())
}

my_results_DF = my_results %>% map(clean_results,genes)


rename_cols <- function(x,label)
{
    nms = names(x)
    nms[nms != "gene"] = paste(label, nms[nms != "gene"],sep = ":")
    names(x) = nms
    x    
}

my_results_DF = my_results_DF %>% map2(names(my_results_DF),rename_cols)


my_results = Reduce(function(...) merge(..., by='gene', all.x=TRUE),my_results_DF) %>% as.tbl

write_delim(my_results,opt$outfile,delim = "\t")

## plot

pvalue_histogram <- function(mod,my_results)
{
    dd = my_results %>% select(gene,contains(mod))

    if(mod %in% c("cell","treat")){
        dd = dd[,1:7]
    }

    dd = dd %>%
        select(gene,contains("pvalue"))
    names(dd) = c("gene","pvalue")

    dd %>% ggplot(aes(pvalue)) +
        geom_histogram(bins = 50,colour = "blue",fill = "white",boundary = 0) +
        ggtitle(mod)

}


pdf(paste0(opt$figs,"_pvalue_histograms.pdf"),height = 5)
plots = lapply(models,pvalue_histogram,my_results)
u = lapply(plots,print)
dev.off()

volcano_plot <- function(mod,my_results,fdr = opt$fdr)
{
    nms = c("gene","log2FoldChange","pvalue","padj")
    dd = my_results %>% select(gene,contains(mod))

    if( mod %in% c("cell","treat")){
        dd = dd[,1:7]
    }

    dd = dd %>%
        select(gene,contains("Fold"),contains("pval"),contains("padj") )
    names(dd) = nms

    

    dd = dd %>% mutate(col = ifelse(padj <= fdr,"yes","no")) %>% filter(!is.na(padj))

    mm = 1.3 * min(min(dd$log2FoldChange),max(dd$log2FoldChange)) %>% abs

    dd %>% ggplot(aes(log2FoldChange,-log10(pvalue),colour = col))+geom_point()+
        scale_colour_brewer(name = paste("adjusted pval <=",fdr),palette = "Set1")+
        theme(legend.position = "top")+ggtitle(mod)+xlim(-mm,mm)+
        geom_vline(xintercept = 0,linetype = 2)


}

pdf(paste0(opt$figs,"_volcano_plots.pdf"))
plots = lapply(models,volcano_plot,my_results,fdr = 1e-5)
u = lapply(plots,print)
dev.off()

signal2noise_barplot <- function(mod,my_results)
{    

    dd = my_results %>% select(gene,contains(mod))

    if(mod %in% c("cell","treat")){
        dd = dd[,1:7]
    }

    names(dd) = gsub(paste0(mod,":"),"",names(dd) )

    dd = dd %>% arrange(stat) %>%
        mutate(gene = factor(gene,levels = gene),
               rank = seq_len(nrow(dd)))

    dd %>% ggplot(aes(rank,stat))+geom_line()+
        ggtitle(mod)+geom_abline(slope = 0,intercept = 0,linetype = 2)
}

pdf(paste0(opt$figs,"_signal_to_noise_plot.pdf"))
plots = lapply(models,signal2noise_barplot,my_results)
u = lapply(plots,print)
dev.off()


rld = rlog(deseq_model1,blind = FALSE)

pdf(paste0(opt$figs,"_PCA_analysis.pdf"))
u = plotPCA(rld,intgroup = c("cell","treat"))
print(u)
dev.off()

library(pheatmap)

rowVars <- function (x,na.rm = TRUE)
{
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}


mat = assay(rld)[,]
mat = mat - rowMeans(mat)
rownames(mat) = genes

r = viridis::viridis(1e3,option = "D")

K = 20
topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)
pheatmap(mat[topgenes,],annotation_col = coldata,color = r,
             filename = paste0(opt$figs,"_Heatmap_top",K,"genes.pdf"))

K = 50
topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)
pheatmap(mat[topgenes,],annotation_col = coldata,color = r,
             filename = paste0(opt$figs,"_Heatmap_top",K,"genes.pdf"))

K = 100
topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)
pheatmap(mat[topgenes,],annotation_col = coldata,show_rownames = FALSE,color = r,
             filename = paste0(opt$figs,"_Heatmap_top",K,"genes.pdf"))

K = 5e2
topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)
pheatmap(mat[topgenes,],annotation_col = coldata,show_rownames = FALSE,color = r,
             filename = paste0(opt$figs,"_Heatmap_top",K,"genes.pdf"))

K = 1e3
topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)
pheatmap(mat[topgenes,],annotation_col = coldata,show_rownames = FALSE,color = r,
             filename = paste0(opt$figs,"_Heatmap_top",K,"genes.pdf"))

K = 5e3
topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)
pheatmap(mat[topgenes,],annotation_col = coldata,show_rownames = FALSE,color = r,
             filename = paste0(opt$figs,"_Heatmap_top",K,"genes.pdf"))





