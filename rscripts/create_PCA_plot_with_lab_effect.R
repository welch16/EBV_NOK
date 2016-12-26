#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--metadatafile",action = "store_true",type = "character",
                help = "Input tsv file consisting of the metadata used to generate the PCA plot:
                          file | cell | treatment | lab | rep "),
    make_option("--ntop",action = "store_true",type = "numeric", default = 500,
                help = "Number of top genes (in variance) used to plot the principal components"),
    make_option("--figsfile", action = "store_true",type = "character",default = "./Rplots",
                help = "File where the plot is going to be stored")                        
    )

opt = parse_args(OptionParser(option_list = optList))

##opt$metadatafile = "data/metadata/PCA_definition.tsv"

stopifnot(file.exists(opt$metadatafile))

library(base,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(tximport,quietly = TRUE)
library(readr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(tidyr,quietly = TRUE)
library(ggplot2,quietly = TRUE)
library(ggrepel,quietly = TRUE)

metadata = read_tsv(opt$metadatafile)

txdata = tximport(metadata$file,
                  type = "rsem" , importer = read_tsv)

library(DESeq2,quietly = TRUE)



deseq = DESeqDataSetFromMatrix(
    round(txdata[["counts"]]),colData = as.data.frame(metadata),
                               design= ~ cell + treatment  + lab)

deseq = deseq[rowSums(counts(deseq)) > 1,]

rl = rlog(deseq)

rowVars <- function (x,na.rm = TRUE)
{
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}


my_PCA_plot <- function(rl,ntop = 500,x = 1,y = 2)
{
    
    rv = rowVars(assay(rl))

    select = order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]

    pca = prcomp(t(assay(rl)[select, ]))
    percentVar = pca$sdev^2/sum(pca$sdev^2)

    intgroup.df =  as.data.frame(colData(rl)[, ,drop = FALSE]) %>% as.tbl

    dt = intgroup.df %>% mutate(file = basename(file)) %>%
        mutate(x = pca$x[,x], y = pca$x[,y])

    dt %>%
        ggplot(aes(x,y,colour = interaction(cell,treatment),shape = lab))+
        geom_point(size = 2)+coord_fixed()+
        xlab(paste0("PC",x,": ",round(percentVar[x] * 100), "% variance"))+
        ylab(paste0("PC",y,": ",round(percentVar[y] * 100), "% variance"))+
        scale_color_brewer(palette = "Dark2",name = "Cell.Treatment")+        
        scale_shape_manual(name = "Lab",values = 1:2)+
        geom_label_repel(aes(x,y,label = rep),size = 2,
                         box.padding = unit(.3,"lines"),
                         point.padding = unit(.7,"lines"),
                         show.legend = FALSE)+theme_bw()

}


pdf(opt$figsfile,width = 9 ,height = 9)
plots = list()
plots[[1]] = my_PCA_plot(rl,500,x = 1, y = 2)
plots[[2]] = my_PCA_plot(rl,500,x = 1, y = 3)
plots[[3]] = my_PCA_plot(rl,500,x = 1, y = 4)
plots[[4]] = my_PCA_plot(rl,500,x = 2, y = 3)
plots[[5]] = my_PCA_plot(rl,500,x = 2, y = 4)
u = lapply(plots,print)
dev.off()
