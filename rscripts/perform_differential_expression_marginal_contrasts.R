#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--mono_files",action = "store_true",type = "character",
                help = "Files corresponding to untreated samples in cell"),
    make_option("--treat_files",action = "store_true", type = "character",
                help = "Files corresponding to treated samples in cell"),
    make_option("--treatment",action = "store_true",type = "character",default = "no,yes",
                help = "Names of the two treatments separated by ','. The default value is 'no,yes'"),
    make_option("--outfile", action = "store_true",default = tempfile(),type = "character",
                help = "Name of the outfile where the test results are saved."),
    make_option(c("-t","--type"),action = "store_true",default = "rsem",type = "character",
                help = "Tool used to quantify transcripts. The default value is 'rsem'"),
    make_option("--plot_title",action = "store_true",type = "character",
                help = "Title used in plots"),
    make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
                help = "Directory where the figures are being saved"),
    make_option("--fdr",action = "store_true",type = "numeric",default = 1e-3,
               help = "FDR used to call differentially expressed genes")
  )

opt = parse_args(OptionParser(option_list = optList))

## opt$mono_files = 'data/RSEM/hg19/RNAseq-noks-CaFBS-rep?.genes.results'
## opt$treat_files = 'data/RSEM/hg19/RNAseq-akata-noks-CaFBS-rep?.genes.results'
## opt$treatment = 'NOKS,EBV'

## opt$mono_files = "data/RSEM/hg19/RNAseq-Noks-mono-rep?.genes.results"
## opt$treat_files = "data/RSEM/hg19/RNAseq-Noks-MC-rep?.genes.results"
## opt$plot_title = "Scott NOKS: MC vs mono"

## opt$treat_files = "data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep2.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep3.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep4.genes.results"
## opt$B_Tr = "data/RSEM/hg19/RNAseq-Noks_EBV-MC-rep?.genes.results"  

## opt$cells = "NOK,EBV"
## opt$treatments = "none,MC"  

library(base,quietly = TRUE)
library(tidyverse,quietly = TRUE)
library(tximport,quietly = TRUE)

separateFiles <- function(ff)
{
  if(grepl(",",ff)){
      ff = ff %>% strsplit(',') %>% unlist     
  }else{
      ff = Sys.glob(ff)      
  }
  ff
}

opt$mono_files = separateFiles(opt$mono_files)
opt$treat_files = separateFiles(opt$treat_files)

## print(opt)


source("rfuns/geneExpression_analysis.R")

treat = opt$treatment %>% strsplit(",") %>% unlist


names(opt$mono_files) = strsplit(opt$mono_files,"rep") %>% map_chr(.f = function(x)x[2]) %>% {paste0(treat[1],"-",gsub(".genes.results","",.))}
names(opt$treat_files) = strsplit(opt$treat_files,"rep") %>% map_chr(.f = function(x)x[2]) %>% {paste0(treat[2],"-",gsub(".genes.results","",.))}

stopifnot(all(file.exists(opt$mono_files)),
          all(file.exists(opt$treat_files)))

stopifnot(opt$type %in% c("none", "kallisto", "salmon", "sailfish","rsem"))

allfiles = c(opt$mono_files,opt$treat_files)

## load files with tximport
txdata = tximport(allfiles,
                  type = opt$type , importer = read_tsv)

coldata = data.frame( treat = c(rep(treat[1],length(opt$mono_files)),
                                rep(treat[2],length(opt$treat_files))))
coldata$treat = factor(coldata$treat, levels = treat)
rownames(coldata) = names(allfiles)

theme_set(theme_bw())

library(DESeq2,quietly = TRUE)

deseq = DESeqDataSetFromMatrix(
    round(txdata[["counts"]]),colData = coldata,
                               design= ~ treat)

deseq = deseq[rowSums(counts(deseq)) > 1,]

## 
deseqModel = DESeq(deseq,minReplicatesForReplace = Inf)

my_results = results(deseqModel,cooksCutoff = FALSE,tidy = TRUE) %>% as.tbl %>%
    separate(row,into = c("ensembl_id","gene_id"),sep = "\\_")

write_delim(my_results,opt$outfile,delim = "\t")

## plot

pvalue_histogram <- function(mod,my_results)
{
    my_results %>% ggplot(aes(pvalue)) +
        geom_histogram(bins = 50,colour = "blue",fill = "white",
                       boundary = 0) +
        ggtitle(mod)

}

pdf(paste0(opt$figs,"_pvalue_histograms.pdf"),height = 5)
pvalue_histogram(opt$plot_title,my_results) %>% print
dev.off()

volcano_plot <- function(mod,my_results,fdr = opt$fdr)
{
    mm = my_results %>% select(log2FoldChange) %>% {1.3 * max(-min(.),max(.))}

    my_results %>% filter(!is.na(padj)) %>% mutate(col = ifelse(padj <= fdr,"yes","no")) %>%
    ggplot(aes(log2FoldChange,-log10(pvalue),colour = col))+geom_point()+
        scale_colour_brewer(name = paste("adjusted pval <=",fdr),palette = "Set1")+
        theme(legend.position = "top")+ggtitle(mod)+xlim(-mm,mm)+
        geom_vline(xintercept = 0,linetype = 2)


}

pdf(paste0(opt$figs,"_volcano_plots.pdf"))
volcano_plot(opt$plot_title,my_results,opt$fdr) %>% print
dev.off()

signal2noise_barplot <- function(mod,my_results)
{
    mm = my_results %>% select(stat) %>% {1.2 * max(-min(.),max(.))}

    my_results %>% arrange(stat) %>% mutate(rank = seq_len(nrow(my_results))) %>%
        ggplot(aes(rank,stat))+geom_line()+ylim(-mm,mm)+
        ggtitle(mod)+geom_abline(slope = 0,intercept = 0,linetype = 2)
}

pdf(paste0(opt$figs,"_signal_to_noise_plot.pdf"))
signal2noise_barplot(opt$plot_title,my_results) %>% print
dev.off()


rld = rlog(deseqModel,blind = FALSE)

pdf(paste0(opt$figs,"_PCA_analysis.pdf"))
u = plotPCA(rld,intgroup = c("treat"))+scale_color_brewer(palette = "Set1",name = "Treatment")+
    theme(legend.position = "top")
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



pheat_pval <- function(rld,my_results,K)
{
    mat = assay(rld)[,]
    mat = mat - rowMeans(mat)
    rownames(mat) = my_results$gene_id
    r = viridis::viridis(1e3,option = "D")

    topgenes = head(order(my_results$pvalue),K)
    pheatmap(mat[topgenes,],annotation_col = coldata,color = r,
             show_rownames = K <= 75,
             filename = paste0(opt$figs,"_Heatmap_bottom",K,"genes_pval.pdf"))
    
}    

pheat_var <- function(rld,my_results,K)
{
    mat = assay(rld)[,]
    mat = mat - rowMeans(mat)
    rownames(mat) = my_results$gene_id
    r = viridis::viridis(1e3,option = "D")

    topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)

    pheatmap(mat[topgenes,],annotation_col = coldata,color = r,
             show_rownames = K <= 75,
             filename = paste0(opt$figs,"_Heatmap_top",K,"genes_var.pdf"))
    
}    


Ks = c(20,50,100,5e2,1e3,5e3,1e4)

a = lapply(Ks,function(x)pheat_pval(rld,my_results,x))
b = lapply(Ks,function(x)pheat_var(rld,my_results,x))



