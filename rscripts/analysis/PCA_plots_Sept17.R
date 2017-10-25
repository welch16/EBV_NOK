
rm(list = ls())

library(tidyverse)
library(ggrepel)
library(tximport)

indr = "data/RSEM/hg19/Sept17"
files = indr %>%
    list.files(full.names = TRUE,pattern = "genes.results")

txdata = tximport(files,
                  type = "rsem" , importer = read_tsv)


library(DESeq2)

coldata = tibble(files = basename(files)) %>%
    mutate(
        cell = if_else(
            grepl("clone",files),"EBV_dRdZ",
            if_else(
                grepl("akata",files),"EBV","NOKS")),
        treat = if_else(grepl("methyl",files),"MC","NoTr")) %>%
    as.data.frame() %>%
    mutate(files = NULL)

countmat = round(txdata[["counts"]])
colnames(countmat) = gsub(".genes.results","",basename(files))
rownames(coldata) = colnames(countmat)


deseq = DESeqDataSetFromMatrix(
    countmat,colData = coldata,design =  ~ cell + treat)

deseq = deseq[rowSums(counts(deseq)) > 1,]

rl = rlog(deseq)

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
        ggplot(aes(x,y,colour = interaction(cell,treat)))+
        geom_point(size = 2)+
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

only_first_clones = all_vars %>%
    {.[ !(grepl("clone3",.) | grepl("clone4",.))]}

only_new_clones = all_vars %>%
    {.[ !(grepl("clone1",.) | grepl("clone2",.))]}



figsdr = "figs/PCA_MDS/Sept17"

pdf(file.path(figsdr,"PCA_plots_allReplicates.pdf"))
plots = list(
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 1,y= 2),
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 1,y= 3),
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 1,y= 4),
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 2,y= 3),
    PCA_plot(rl,vars = all_vars,ntop = 1e3,x = 2,y= 4))
u = plots %>% print()
dev.off()

pdf(file.path(figsdr,"PCA_plots_onlyFirstClones.pdf"))
plots = list(
    PCA_plot(rl,vars = only_first_clones,ntop = 1e3,x = 1,y= 2),
    PCA_plot(rl,vars = only_first_clones,ntop = 1e3,x = 1,y= 3),
    PCA_plot(rl,vars = only_first_clones,ntop = 1e3,x = 1,y= 4),
    PCA_plot(rl,vars = only_first_clones,ntop = 1e3,x = 2,y= 3),
    PCA_plot(rl,vars = only_first_clones,ntop = 1e3,x = 2,y= 4))
u = plots %>% print()
dev.off()

pdf(file.path(figsdr,"PCA_plots_onlyNewClones.pdf"))
plots = list(
    PCA_plot(rl,vars = only_new_clones,ntop = 1e3,x = 1,y= 2),
    PCA_plot(rl,vars = only_new_clones,ntop = 1e3,x = 1,y= 3),
    PCA_plot(rl,vars = only_new_clones,ntop = 1e3,x = 1,y= 4),
    PCA_plot(rl,vars = only_new_clones,ntop = 1e3,x = 2,y= 3),
    PCA_plot(rl,vars = only_new_clones,ntop = 1e3,x = 2,y= 4))
u = plots %>% print()
dev.off()
