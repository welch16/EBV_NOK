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
    make_option("--outfile", action = "store_true",default = tempfile(),type = "character",
                help = "Name of the outfile where the test results are saved."),
    make_option(c("-t","--type"),action = "store_true",default = "rsem",type = "character",
                help = "Tool used to quantify transcripts. The default value is 'rsem'"),
    make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
                help = "Directory where the figures are being saved"),
    make_option("--fdr",action = "store_true",type = "numeric",default = 1e-3,
               help = "FDR used to call differentially expressed genes")
  )

opt = parse_args(OptionParser(option_list = optList))

## opt$A_noTr = "data/RSEM/hg19/RNAseq-Noks-mono-rep?.genes.results"
## opt$B_noTr = "data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep2.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep3.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep4.genes.results"
## opt$A_Tr = "data/RSEM/hg19/RNAseq-Noks-MC-rep?.genes.results"
## opt$B_Tr = "data/RSEM/hg19/RNAseq-Noks_EBV-MC-rep?.genes.results"  

## opt$cells = "NOK,EBV"
## opt$treatments = "none,MC"  

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

cells = opt$cells %>% strsplit(",") %>% unlist
treat = opt$treatment %>% strsplit(",") %>% unlist

names(opt$A_noTr) = getRep(opt$A_noTr,paste(cells[1],treat[1],sep = "_"))
names(opt$B_noTr) = getRep(opt$B_noTr,paste(cells[2],treat[1],sep = "_"))
names(opt$A_Tr) = getRep(opt$A_Tr,paste(cells[1],treat[2],sep = "_"))
names(opt$B_Tr) = getRep(opt$B_Tr,paste(cells[2],treat[2],sep = "_"))

stopifnot(all(file.exists(opt$A_noTr)),
          all(file.exists(opt$B_noTr)),
          all(file.exists(opt$A_Tr)),
          all(file.exists(opt$B_Tr)))

stopifnot(opt$type %in% c("none", "kallisto", "salmon", "sailfish","rsem"))

allfiles = c(opt$A_noTr,opt$B_noTr,opt$A_Tr,opt$B_Tr)

## load files with tximport
txdata = tximport(allfiles,
                  type = opt$type , importer = read_tsv)

coldata = data.frame(
    cell = c(rep(cells[1],length(opt$A_noTr)),
             rep(cells[2],length(opt$B_noTr)),
             rep(cells[1],length(opt$A_Tr)),
             rep(cells[2],length(opt$B_Tr))),
    treat = c(rep(treat[1],length(opt$A_noTr)),
             rep(treat[1],length(opt$B_noTr)),
             rep(treat[2],length(opt$A_Tr)),
             rep(treat[2],length(opt$B_Tr))))

coldata$interac = paste(coldata$cell,coldata$treat,sep = "_")

coldata$cell = factor(coldata$cell)
coldata$treat = factor(coldata$treat)
coldata$interac = factor(coldata$interac)

rownames(coldata) = names(allfiles)

theme_set(theme_bw())

library(DESeq2,quietly = TRUE)

deseq = DESeqDataSetFromMatrix(
    round(txdata[["counts"]]),colData = coldata,
                               design= ~ cell + treat + cell:treat)

deseq = deseq[rowSums(counts(deseq)) > 1,]

## 
deseq_model1 = DESeq(deseq,minReplicatesForReplace = Inf)

my_results = list()
my_results[["full"]] = results(deseq_model1,cooksCutoff = FALSE)
my_results[["cell"]] = results(deseq_model1,contrast = c("cell",cells),cooksCutoff = FALSE)
my_results[["treat"]] = results(deseq_model1,contrast = c("treat",rev(treat)),cooksCutoff = FALSE)

get_columns <- function(cond,columns)
{
    if(is.null(cond)){
        out = rep(TRUE,length(columns))
    }else{
        if(grepl(",",cond,fixed = TRUE)){
            cond = cond %>% strsplit(",") %>% unlist
            out = lapply(cond,function(x)grepl(x,columns))
            out = do.call(and,out)
        }else{
            out = grepl(cond,columns)
        }
    }

    columns[out]
}

signal_to_noise <- function(results,model,condA,condB,cells,treats)
{
    
    ## Under their model rhe counts for gene i and condition j are distributed
    ## K_ij ~ NegBin(mu_ij , alpha_i)
    ## where mu_ij is the mean and alpha_i is the gene specific dispersion
    ## We define
    ## signaltonoise_i = log2FC_i / 2* alpha_i

    ## cond1 & cond2 contains the restrictions of the hypothesis test
    ## ',' means 'and'
    ## '_' means 'when'


    num = results$log2FoldChange
    
    ## mu = assays(model)[["mu"]]
    ## columns = mu %>% colnames

    ## columnsA = get_columns(condA,columns)
    ## columnsB = get_columns(condB,columns)

    ## muA = rowMeans(mu[,columns %in% columnsA])
    ## muB = rowMeans(mu[,columns %in% columnsB])
    
    alpha = dispersions(model)
   
    num / (2 * alpha)
}

signal2noise = list()
signal2noise[["full"]] = signal_to_noise(my_results[["full"]],deseq_model1,paste(cells[1],treat[1],sep = ","),
                                         paste(cells[2],treat[2],sep = ","),cells,treat)

signal2noise[["cell"]] = signal_to_noise(my_results[["cell"]],deseq_model1,cells[1],cells[2],cells,treat)
signal2noise[["treat"]] = signal_to_noise(my_results[["treat"]],deseq_model1,treat[2],treat[1],cells,treat)

design(deseq) <- ~ interac
deseq_model2 = DESeq(deseq, minReplicatesForReplace = Inf)

interac = expand.grid(cell = cells,treat = treat) %>% mutate(interac = paste(cell,treat,sep = "_")) %>%
    select(interac)

grepv <- function(pattern, x, ignore.case = FALSE, perl = FALSE, value = FALSE, 
                  fixed = FALSE, useBytes = FALSE, invert = FALSE)
{
  x[grep(pattern,x,ignore.case,perl,value,fixed,useBytes,invert)]
}

## print(cells[1])
## print(interac)
## print(grepv(cells[1],interac[[1]]) %>% rev)
## print(coldata)


my_results[[paste0("cell_",cells[1])]] = results(deseq_model2,
    contrast = c("interac",grepv(cells[1],interac[[1]]) %>% rev),cooksCutoff = FALSE)
my_results[[paste0("cell_",cells[2])]] = results(deseq_model2,
    contrast = c("interac",grepv(cells[2],interac[[1]]) %>% rev),cooksCutoff = FALSE)

signal2noise[[paste0("cell_",cells[1])]] = signal_to_noise(my_results[[paste0("cell_",cells[1])]],
    deseq_model2,paste(cells[1],treat[2],sep = ","),paste(cells[1],treat[1],sep = ","),cells,treat)
signal2noise[[paste0("cell_",cells[2])]] = signal_to_noise(my_results[[paste0("cell_",cells[2])]],
    deseq_model2,paste(cells[2],treat[2],sep = ","),paste(cells[2],treat[1],sep = ","),cells,treat)

my_results[[paste0("treat_",treat[1])]] = results(deseq_model2,
    contrast = c("interac",grepv(treat[1],interac[[1]]) %>% rev), cooksCutoff = FALSE)
my_results[[paste0("treat_",treat[2])]] = results(deseq_model2,
    contrast = c("interac",grepv(treat[2],interac[[1]]) %>% rev), cooksCutoff = FALSE)

signal2noise[[paste0("treat_",treat[1])]] = signal_to_noise(my_results[[paste0("treat_",treat[1])]],
    deseq_model2,paste(treat[1],cells[1],sep = ","),paste(treat[1],cells[2],sep = ","),cells,treat)
signal2noise[[paste0("treat_",treat[2])]] = signal_to_noise(my_results[[paste0("treat_",treat[2])]],
    deseq_model2,paste(treat[2],cells[1],sep = ","),paste(treat[2],cells[2],sep = ","),cells,treat)

models = names(my_results)

genes = my_results[[1]] %>% rownames

#%>% strsplit("_") %>% sapply(function(x)x[2])

clean_results <- function(res,s2n, genes)
{   
    res %>% as.data.frame %>%
        as.tbl %>% mutate(gene = genes) %>%
        select(gene,baseMean,log2FC = log2FoldChange,pvalue,padj) %>%
        mutate(signal2noise = s2n)
}

my_results_DF = mapply(clean_results,my_results,signal2noise,MoreArgs = list(genes),SIMPLIFY = FALSE)

rename_cols <- function(x,label)
{
    nms = names(x)
    nms[nms != "gene"] = paste(label, nms[nms != "gene"],sep = ":")
    names(x) = nms
    x    
}

my_results_DF = mapply(rename_cols,my_results_DF,names(my_results_DF),
                       SIMPLIFY = FALSE)

my_results = Reduce(function(...) merge(..., by='gene', all.x=TRUE),my_results_DF) %>% as.tbl

write_delim(my_results,opt$outfile,delim = "\t")

## plot

pvalue_histogram <- function(mod,my_results)
{
    dd = my_results %>% select(gene,contains(mod))

    if(mod %in% c("cell","treat")){
        dd = dd[,1:6]
    }

    dd = dd %>%
        select(gene,contains("pvalue"))
    names(dd) = c("gene","pvalue")

    dd %>% ggplot(aes(pvalue)) + geom_histogram(bins = 150) +
        ggtitle(mod)

}


pdf(paste0(opt$figs,"_pvalue_histograms.pdf"),height = 5)
plots = lapply(models,pvalue_histogram,my_results)
u = lapply(plots,print)
dev.off()

volcano_plot <- function(mod,my_results,fdr = opt$fdr)
{
    
    nms = c("gene","log2FC","pvalue","padj")
    dd = my_results %>% select(gene,contains(mod))

    if( mod %in% c("cell","treat")){
        dd = dd[,1:6]
    }

    dd = dd %>%
        select(gene,contains("FC"),contains("pval"),contains("padj") )
    names(dd) = nms

    dd = dd %>% mutate(col = ifelse(padj <= fdr,"yes","no")) %>% filter(!is.na(padj))

    mm = 1.3 * min(min(dd$log2FC),max(dd$log2FC)) %>% abs

    dd %>% ggplot(aes(log2FC,-log10(pvalue),colour = col))+geom_point()+
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
        dd = dd[,1:6]
    }

    names(dd) = gsub(paste0(mod,":"),"",names(dd) )

    dd = dd %>% arrange(signal2noise) %>%
        mutate(gene = factor(gene,levels = gene),
               rank = seq_len(nrow(dd)))

    dd %>% ggplot(aes(rank,signal2noise))+geom_line()+
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




