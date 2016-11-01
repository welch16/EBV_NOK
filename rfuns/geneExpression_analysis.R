

## Returns a 'rep:i' tag for each element in files vector
getRep <- function(files,name = "rep")
{
  paste0(name,seq_along(files))
}

## using geometric mean, we average per replicate and create a tibble
gemean <- function(x,...)exp(mean(log(x),...))
rowGeoMeans <- function(mat)apply(mat,1,gemean)

## creates a data table using the tximport outputs
createDataTable <- function(A,A_diff,B,B_diff,what = c("abundance","counts"))
{
  require(dplyr)
  require(magrittr)

  genes = A[[what]] %>% rownames %>%
    strsplit("_",fixed = TRUE) %>%
    sapply(function(x)x[2])
  
  DT = tibble(
    genes = genes,
    A = rowGeoMeans(A[[what]]),
    A_diff = rowGeoMeans(A_diff[[what]]),
    B = rowGeoMeans(B[[what]]),
    B_diff = rowGeoMeans(B_diff[[what]]))

  DT = DT %>%
    mutate(log2FC_A = log2(1 + A_diff) - log2(1 + A),
           log2FC_B = log2(1 + B_diff) - log2(1 + B))

  DT
}

geneExpression_count_plot <- function(DT,xvar,yvar,opt,sc)
{
  require(viridis)
  require(scales)
  require(ggplot2)
  pal = viridis(1e3, option = "D")
  
  ggplot(DT,aes_string(xvar,yvar))+stat_binhex(bins = 140) +
    scale_fill_gradientn(colours = pal,trans = "log10",
                         labels = trans_format('log10',math_format(10^.x)),
                         guide = guide_colourbar(title = NULL,
                           barheight = unit(0.92,"npc"),
                           barwidth = unit(0.01,"npc")))+
    xlim(-sc,sc)+ylim(-sc,sc)+
    xlab(opt$xlab)+ylab(opt$ylab)
                          
}

DESeq_plots <- function(countDataSet, results,fdr)
{
  require(viridis,quietly = TRUE)
  require(ggplot2,quietly = TRUE)
  require(scales,quietly = TRUE)
  require(base,quietly = TRUE)
  require(dplyr,quietly = TRUE)
  pal = viridis(1e3, option = "D")
  plots = list()
  plots[[1]] = DESeq_mean_vs_dispersion(countDataSet,results,pal,fdr)
  plots[[2]] = DESeq_MA_plot(results,pal,fdr)
  plots[[3]] = DESeq_pval_hist(results,fdr)
  plots
}

DESeq_mean_vs_dispersion <- function(countDataSet,results,pal,fdr)
{
  info = fitInfo(countDataSet)
  px = rowMeans(counts(countDataSet, normalized = TRUE))

  DT = tibble(x = px , y = info$perGeneDispEsts,
    fit = info$fittedDispEsts , padj = results$padj)  
  miny = DT %>% select(y) %>% filter(y > 0) %>% min
  DT = DT %>% filter(x > 0)  %>% filter(y >= miny)
  
  p = ggplot(DT,aes(x,y))+stat_binhex(bins = 140) +
    scale_fill_gradientn(colours = pal,trans = "log10",
                         labels = trans_format('log10',math_format(10^.x)),
                         guide = guide_colourbar(title = NULL,
                           barheight = unit(0.92,"npc"),
                           barwidth = unit(0.01,"npc")))+
    scale_x_log10()+scale_y_log10()+
    geom_line(data = DT ,aes(x,fit),linetype = 2,colour = "orange",size = .8)+
    xlab("mean of normalized counts")+ylab("dispersion")+
    geom_density2d(data = DT %>% filter(padj <= fdr),colour = "red",size = .5)
                                
}

DESeq_MA_plot <- function(results,pal,fdr)
{  
  results = results %>% filter(baseMean != 0)%>%
    filter(foldChange < quantile(foldChange,prob = .99))
  yl = results %>% select(foldChange) %>%
    mutate(aa = log2(foldChange)) %>% range %>% abs %>% max %>% round
  results = results %>% mutate(log2 = log2(foldChange)) %>%
    filter( between(log2,-yl,yl))
  ggplot(results,
         aes(baseMean,log2(foldChange)))+stat_binhex(bins = 140)+
    scale_fill_gradientn(colours = pal,trans = "log10",
                         labels = trans_format('log10',math_format(10^.x)),
                         guide = guide_colourbar(title = NULL,
                           barheight = unit(0.92,"npc"),
                           barwidth = unit(0.01,"npc")))+
    scale_x_log10()+ylim(-yl ,yl)+
    xlab("mean of normalized counts")+ylab("log2 fold change")+
    geom_abline(slope = 0,intercept = 0,linetype = 2,colour = "orange",size = .8)+
    geom_density2d(data = results %>% filter(padj <= fdr),colour = "red",size = .5)

}

DESeq_pval_hist <- function(results,fdr)
{

  results = results %>% filter(baseMean != 0 ) %>%
    mutate( dd = ifelse(padj <= fdr,"yes","no"))

  ggplot(results,
         aes(pval,fill = dd))+geom_histogram(binwidth = .025)+
    xlab("p.value")+ylab("number of genes")+
    theme(legend.position = "top")+
    scale_fill_manual(values = c("black","red") , name  = paste0("fdr <= ",fdr))

}
