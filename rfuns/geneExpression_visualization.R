

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
