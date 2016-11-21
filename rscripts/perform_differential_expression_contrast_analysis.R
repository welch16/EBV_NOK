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
                help = "Names of the two cell lines separated by ','. The default value is 'A,B'"),
    make_option("--treatments",action = "store_true",type = "character",default = "no,yes",
                help = "Names of the two treatments separated by ','. The default value is 'no,yes'"),
    make_option("--treattestfile", action = "store_true",default = tempfile(),type = "character",
                help = "Name of the outfile where the output list of genes is saved"),
    make_option("--celltestfile",action = "store_true",default = tempfile(),type = "character",
                help = "Name of the outfile where the differentially expressed genes by cell are saved"),
    make_option(c("-t","--type"),action = "store_true",default = "rsem",type = "character",
                help = "Tool used to quantify transcripts. The default value is 'rsem'"),
    make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
              help = "Directory where the figures are being saved")
  ## make_option("--fdr",action = "store_true",type = "numeric",default = 0.05,
  ##             help = "FDR used to call differentially expressed genes")
  )

opt = parse_args(OptionParser(option_list = optList))

library(base,quietly = TRUE)
library(magrittr,quietly = TRUE)

separateFiles <- function(ff)
{
  if(grepl(",",ff)){
      ff = ff %>% strsplit(',') %>% unlist     
  }else{
      ff = Sys.glob(ff)      
  }
  ff
}
 
opt$A_noTr = separateFiles(opt$A_noTr)
opt$B_noTr = separateFiles(opt$B_noTr)
opt$A_Tr = separateFiles(opt$A_Tr)
opt$B_Tr = separateFiles(opt$B_Tr)

library(tximport,quietly = TRUE)
library(readr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(tidyr,quietly = TRUE)
library(ggplot2,quietly = TRUE)

source("rfuns/geneExpression_analysis.R")

names(opt$A_noTr) = getRep(opt$A_noTr,"A_noTr")
names(opt$B_noTr) = getRep(opt$B_noTr,"B_noTr")
names(opt$A_Tr) = getRep(opt$A_Tr,"A_Tr")
names(opt$B_Tr) = getRep(opt$B_Tr,"B_Tr")

stopifnot(all(file.exists(opt$A_noTr)),
          all(file.exists(opt$B_noTr)),
          all(file.exists(opt$A_Tr)),
          all(file.exists(opt$B_Tr)))

stopifnot(opt$type %in% c("none", "kallisto", "salmon", "sailfish","rsem"))

allfiles = c(opt$A_noTr,opt$B_noTr,opt$A_Tr,opt$B_Tr)

## load files with tximport
txdata = tximport(allfiles,
                  type = opt$type , importer = read_tsv)

cells = opt$cells %>% strsplit(",") %>% unlist
treat = opt$treatment %>% strsplit(",") %>% unlist


coldata = data.frame(
    cell = c(rep(cells[1],length(opt$A_noTr)),
             rep(cells[2],length(opt$B_noTr)),
             rep(cells[1],length(opt$A_Tr)),
             rep(cells[2],length(opt$B_Tr))),
    treat = c(rep(treat[1],length(opt$A_noTr)),
             rep(treat[1],length(opt$B_noTr)),
             rep(treat[2],length(opt$A_Tr)),
             rep(treat[2],length(opt$B_Tr))))

rownames(coldata) = names(allfiles)


theme_set(theme_bw())

library(DESeq2,quietly = TRUE)

deseq = DESeqDataSetFromMatrix(
    round(txdata[["counts"]]),colData = coldata,
                               design= ~ cell + treat)
deseq = deseq[rowSums(counts(deseq)) > 1,]

deseq_model = DESeq(deseq)

deseq_results = results(deseq_model,contrast = c("treat",rev(treat)))
deseq_cell = results(deseq_model,contrast = c("cell",cells))

genes = deseq_results %>% rownames %>% strsplit("_") %>%
    sapply(function(x)x[2])

results = deseq_results %>% as.data.frame %>%
    as.tbl %>% mutate(gene = genes) %>%
    select(gene,padj,log2FC = log2FoldChange,everything())

print(
    results %>% arrange(desc(log2FC)) %>% head(10)
)

print(
    results %>% arrange(log2FC) %>% head(10)
)

cell_results = deseq_cell %>% as.data.frame %>% as.tbl %>%
    mutate(gene = genes) %>%
    select(gene,padj,log2FC = log2FoldChange,everything())

print(
    cell_results %>% arrange(desc(log2FC)) %>% head(10)
)

print(
    cell_results %>% arrange(log2FC) %>% head(10)
)


rld = rlog(deseq_model,blind = FALSE)

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

K = 20
topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)
pheatmap(mat[topgenes,],annotation_col = coldata,
             filename = paste0(opt$figs,"_Heatmap_top",K,"genes.pdf"))

K = 50
topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)
pheatmap(mat[topgenes,],annotation_col = coldata,
             filename = paste0(opt$figs,"_Heatmap_top",K,"genes.pdf"))

write_delim(results %>% select(gene,baseMean,log2FC,pvalue,padj),
            opt$treattestfile , delim = "\t")
write_delim(cell_results %>% select(gene,baseMean,log2FC,pvalue,padj),
            opt$celltestfile,delim = "\t")



