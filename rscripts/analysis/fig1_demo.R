

library(patchwork)
library(scales)
library(tidyverse)
library(DESeq2)
library(hexbin)
library(viridis)

## From Mark's description

## Compare NOKS vs NOKS+EBV (undifferentiated) with Scott's data
##
## Volcano plot
## Pathway analysis with GSEA
## Possibly other pathway analysis

datadr = "data/RSEM/hg19"
figsdr = "figs/fig2"

files = datadr %>%
    list.files(full.names = TRUE) %>%
    str_subset("genes.results") %>%
    str_subset("Nok")

rsem = files %>%
    set_names(
        basename(.) %>%
        str_replace(".genes.results","")) %>%
    map(read_tsv) %>%
    tibble(
        file = names(.),
        rsem = .)  %>%
    mutate(
        cell = if_else(file %>% str_detect("EBV"),"NOKS_EBV","NOKS"),
        treat = if_else(file %>% str_detect("mono"),"None","MC"))

## quick check of count histograms
## 

library(ggjoy)

theme_set(theme_bw())

log2_count_histogram = rsem %>%
    unnest() %>%
    ggplot(aes(expected_count))+
    geom_histogram(aes(fill = interaction(cell,treat)),
                   bins = 45,colour = "black")+
    scale_fill_brewer(palette = "Pastel1",name = "")+
    facet_grid( file ~ .)+
    theme(
        legend.position = "top",
        strip.text.y = element_text(angle = 0))+
    scale_x_log10(labels = comma)+
    geom_vline(xintercept = 100,colour = "red",
               linetype = 2)

ggsave(file = file.path(figsdr,"Histogram_log_expected_counts_scott.pdf"),
       plot = log2_count_histogram,
       height = unit(8,"in"),
       width = unit(6,"in"))

## RLE analysis

rle_plot <- function(rsem)
{
    medians = rsem %>%
        unnest() %>%
        group_by(gene_id) %>%
        summarize(
            median = median(expected_count)) %>%
        ungroup()
    
    rsem = rsem %>%
        unnest() %>%
        inner_join(medians,by = "gene_id") %>%
        mutate(
            rle = log2( expected_count / median)) %>%
        mutate(median = NULL)
    
    rsem %>% 
        ggplot(aes(rle))+
        geom_histogram(aes(fill = interaction(cell,treat)),
                       bins = 45,colour = "black")+
    scale_fill_brewer(palette = "Pastel1",name = "")+
        facet_grid( file ~ .)+
        theme(
            legend.position = "top",
            strip.text.y = element_text(angle = 0))+
        geom_vline(xintercept = 0,colour = "red",
                   linetype = 2)+
    scale_x_continuous(limits = c(-5,5))
}

rle_histogram = rle_plot(rsem)

ggsave(file = file.path(figsdr,"Histogram_RLE_scott.pdf"),
       plot = rle_histogram,
       height = unit(8,"in"),
       width = unit(6,"in"))

logs = read_tsv(file.path(figsdr,"Scott_aligned_reads_table.tsv"))

library(ggrepel)

summary_plot =
    logs %>%
    ggplot(
        aes(aligned,aligned_percent))+
    geom_point()+
    geom_text_repel(aes(label = replicate))+
    facet_grid(treat ~ cell)+
    theme(axis.text.x = element_text(angle = 25,
                                     hjust = 1))+
    scale_x_continuous(
        labels = comma)+
    scale_y_continuous(
        labels = percent)+
    xlab("Number of aligned reads")+
    ylab("Percentage of aligned reads")

ggsave(file = file.path(figsdr,"Scatter_aligned_reads_scott.pdf"),
       plot = summary_plot,
       height = unit(5,"in"),
       width = unit(5.5,"in"))

to_remove = logs %>%
    arrange(aligned,aligned_percent) %>%
    head(3) %>% pluck("file")

ggsave(file = file.path(figsdr,"Histogram_RLE_scott_woLowAlignment.pdf"),
       plot = rle_plot(rsem %>%
                       filter(!file %in% to_remove)),
       height = unit(8,"in"),
       width = unit(6,"in"))

## DESeq2 analysis

pre_count = rsem %>%
    mutate(
        rsem = map(rsem,select,gene_id,expected_count)) %>% 
    select(file,rsem) %>%
    unnest() %>%
    mutate(
        expected_count = round(expected_count,0))

count_matrix = pre_count %>% 
    spread(file, expected_count)

## this step is to remove genes with low counts,
summary_count_matrix =
    pre_count %>%
    group_by(
        gene_id) %>%
    summarize(
        total = sum(expected_count),
        mean_counts = mean(expected_count),
        median_counts = median(expected_count),
        q25_counts = quantile(expected_count,.25),
        q75_counts = quantile(expected_count,.75))


coldata = logs %>%
    select(file,cell,treat) %>%
    as.data.frame() %>% 
    tibble::column_to_rownames("file") %>%
    mutate(
        interac = paste(cell,treat,sep = "."))

genes_to_remove = summary_count_matrix %>%
    filter(
        total <= 100)

count_matrix = count_matrix %>%
    anti_join(
        genes_to_remove %>%
        select(gene_id),by ="gene_id") %>%
    as.data.frame() %>%
    tibble::column_to_rownames("gene_id") %>%
    as.matrix()

tpm_summary = rsem %>%
    mutate(
        rsem = map(rsem,select,gene_id,TPM)) %>% 
    select(file,rsem) %>%
    unnest()  %>%
    anti_join(
        genes_to_remove %>%
        select(gene_id),by ="gene_id")  %>%
    inner_join(
        logs %>%
        select(file,cell,treat),by = "file") %>%
    filter(treat == "None") %>%
    group_by(cell,gene_id) %>%
    summarize(
        mean_TPM = mean(TPM)) %>%
    ungroup() %>%
    spread(cell,mean_TPM)

deseq = DESeqDataSetFromMatrix(
    count_matrix,colData = coldata,
    design = ~ interac)

## fitting model
deseq = DESeq(deseq)

rl = rlog(deseq)

no_treat_noks_vs_ebv =
    results(deseq,
        contrast = c("interac","NOKS.None","NOKS_EBV.None"),
        tidy  = TRUE) %>%
    as_tibble() %>%
    dplyr::rename(gene_id = row) %>%
    inner_join(
        tpm_summary, by = "gene_id") 

pal = viridis(1e3,option = "D")

diff_genes = no_treat_noks_vs_ebv %>%
    filter(padj <= .05) %>%
    pluck("gene_id")

library(pheatmap)
cc = assay(rl)

pheatmap(cc[rownames(cc) %in% diff_genes,],
         color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100))


volcano_plot <- function(res,name)
{

    pal = viridis::viridis(1e3, option = "D")

    res = res %>%
        mutate(
            log10pval = -log10(pvalue))
    
    yul = quantile( res$log10pval, .99)

    xl = quantile(abs(res$log2FoldChange),.9999)
    
    res %>%
        ggplot(aes(log2FoldChange,log10pval))+stat_binhex(bins = 51)+
        ggtitle(name)+
        scale_y_continuous(limits = c(0,yul))+
        scale_x_continuous(limits = c(-xl,xl))+
        ylab(expression(-log[10](p.value)))+
        scale_fill_gradientn(
            colours = pal,trans = "log10",
                         labels = trans_format('log10',math_format(10^.x)),
                         guide = guide_colourbar(title = NULL,
                           barheight = unit(0.92,"npc"),
                           barwidth = unit(0.01,"npc")))+
        geom_vline(xintercept = 0,linetype = 2)

}


## MA plot


hexbin_maplot =
    no_treat_noks_vs_ebv %>%
    ggplot(aes(baseMean,log2FoldChange))+
    stat_binhex(bins = 81)+
    scale_x_continuous(labels = comma,limits = c(1,5e3))+
    scale_fill_gradientn(
        colours = pal,trans = "log10",
        labels = trans_format('log10',math_format(10^.x)),
        guide = guide_colourbar(title = NULL,
                                barheight = unit(0.92,"npc"),
                                barwidth = unit(0.01,"npc")))+
    scale_y_continuous(limits = c(-5,5))+
    ylab("log2 NOKS / NOKS_EBV")+
    geom_smooth(
        method = "loess",
        colour = "orange",
        se = FALSE)+
    xlab("Weighted average of number of reads per gene")


ggsave(
    file = file.path(figsdr,"MA_plot_Hexbin_Undiff_NOKS_vs_NOKS_EBV.pdf"),
    plot = hexbin_maplot,
    dpi = 400)

## this is just to check that there is no bias towards one cell-type,
## from visually exploring it appears it doesn't

## "quick" PCA analysis


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

    intgroup.df =  colData(rl)[, ,drop = FALSE] %>%
        as.data.frame() %>% 
        tibble::rownames_to_column() %>%
        as_tibble() %>%
        dplyr::rename(
                   file = rowname) %>%
        mutate(
            replicate = map_chr(
                strsplit(file,"-"), ~ .[length(.)]))
        
    dt = intgroup.df %>%
        mutate(file = basename(file)) %>%
        mutate(x = pca$x[,x], y = pca$x[,y])

    dt %>%
        ggplot(aes(x,y,colour = interaction(cell,treat)))+
        geom_point(size = 2)+coord_fixed()+
        xlab(paste0("PC",x,": ",round(percentVar[x] * 100), "% variance"))+
        ylab(paste0("PC",y,": ",round(percentVar[y] * 100), "% variance"))+
        scale_color_brewer(palette = "Set1",name = "Cell.Treatment")+        
        geom_label_repel(aes(x,y,label = replicate),size = 2,
                         box.padding = unit(.3,"lines"),
                         point.padding = unit(.7,"lines"),
                         show.legend = FALSE)+theme_bw()

}

ggsave(
    file = file.path(figsdr,"PCA_plot_scott_Comp1_vs_Comp2.pdf"),
    plot = my_PCA_plot(rl,ntop = 500),dpi = 400)

ggsave(
    file = file.path(figsdr,"PCA_plot_scott_Comp1_vs_Comp3.pdf"),
    plot = my_PCA_plot(rl,ntop = 500,x = 1,y = 3),dpi = 400)


## only first two component, ~ 85% of the variation is occurring due
## to the treatment. The third components shows everything very
## shuffled and without any clear interpretation.




## Proportion of variance explained:

## rv = rowVars(assay(rl))
## select = order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
## pca = prcomp(t(assay(rl)[select, ]))
## percentVar = pca$sdev^2/sum(pca$sdev^2)

## > summary(pca)
## Importance of components:
##                            PC1     PC2     PC3     PC4     PC5     PC6     PC7
## Standard deviation     36.0815 8.79209 7.44406 4.69532 4.28396 4.00685 2.85140
## Proportion of Variance  0.8553 0.05079 0.03641 0.01448 0.01206 0.01055 0.00534
## Cumulative Proportion   0.8553 0.90612 0.94252 0.95701 0.96907 0.97961 0.98496
##                            PC8     PC9   PC10    PC11    PC12   PC13    PC14
## Standard deviation     2.32472 2.16687 1.7887 1.72152 1.60296 1.3507 1.25679
## Proportion of Variance 0.00355 0.00308 0.0021 0.00195 0.00169 0.0012 0.00104
## Cumulative Proportion  0.98851 0.99159 0.9937 0.99564 0.99733 0.9985 0.99957
##                           PC15      PC16
## Standard deviation     0.81357 9.605e-15
## Proportion of Variance 0.00043 0.000e+00
## Cumulative Proportion  1.00000 1.000e+00

