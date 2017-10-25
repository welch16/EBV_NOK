rm(list = ls())

library(tidyverse)
library(ggrepel)
library(tximport)

## The bad replicates to be removed are clone2 and clone4

indr <- "data/RSEM/hg19/Sept17"
files <- indr %>%
    list.files(full.names = TRUE, pattern = "genes.results") %>% {
      .[ !grepl("clone2", .) & !grepl("clone4", .)]
    }

txdata <- tximport(files,type = "rsem", importer = read_tsv)

library(DESeq2)

coldata <- tibble(files = basename(files)) %>%
    mutate(
        cell = if_else(
            grepl("clone", files), "EBV_dRdZ",
            if_else(
                grepl("akata", files), "EBV", "NOKS")),
        treat = if_else(grepl("methyl", files), "MC", "NoTr")) %>%
    as.data.frame() %>%
    mutate(files = NULL)

countmat <- round(txdata[["counts"]])
colnames(countmat) <- gsub(".genes.results", "", basename(files))
rownames(coldata) <- colnames(countmat)


deseq <- DESeqDataSetFromMatrix(
    countmat, colData = coldata, design =  ~ cell + treat)

deseq <- deseq[rowSums(counts(deseq)) > 1, ]

rl <- rlog(deseq)

row_vars <- function (x, na.rm = TRUE){
    sqr <- function(x) x * x
    n <- rowSums(!is.na(x))
    n[n <= 1] <- NA
    rowSums(sqr(x - rowMeans(x, na.rm = na.rm)), na.rm = na.rm) / (n - 1)
}

theme_set( theme_bw())

PCA_plot <- function(rl,
  vars = rownames(coldata),
  ntop = 500, x = 1, y=2) {
    my_data <- assay(rl)[, vars]
    rv <- row_vars(my_data)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(my_data[select, ]))
    percent_var <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    my_pca <- colData(rl)[vars, , drop = FALSE] %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(
            sample = vars,
            rep = vars %>%
                strsplit("-") %>%
                map_chr( ~ .[length(.)]),
            x = pca$x[, x],
            y = pca$x[, y])
    my_pca %>%
        ggplot(aes(x, y, colour = interaction(cell, treat))) +
        geom_point(size = 2) +
        theme(legend.position = "top") +
        coord_fixed() +
        xlab(paste0("PC", x, ": ", round(percent_var[x] * 100), "% variance")) +
        ylab(paste0("PC", y, ": ", round(percent_var[y] * 100), "% variance")) +
        scale_color_brewer(palette = "Dark2", name = "Cell.Treatment") +
        geom_label_repel(aes(x, y, label = rep), size = 3,
                         box.padding = unit(.3, "lines"),
                         point.padding = unit(.7, "lines"),
                         show.legend = FALSE)
}

all_vars <- rownames(coldata)
figsdr <- "figs/PCA_MDS/Sept17"

pdf(file.path(figsdr, "PCA_plots_allReplicates_woContamination.pdf"))
plots <- list(
    PCA_plot(rl, vars = all_vars, ntop = 1e3, x = 1, y = 2),
    PCA_plot(rl, vars = all_vars, ntop = 1e3, x = 1, y = 3),
    PCA_plot(rl, vars = all_vars, ntop = 1e3, x = 1, y = 4),
    PCA_plot(rl, vars = all_vars, ntop = 1e3, x = 2, y = 3),
    PCA_plot(rl, vars = all_vars, ntop = 1e3, x = 2, y = 4))
u <- plots %>% print()
dev.off()
