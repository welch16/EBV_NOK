

separateFiles <- function(ff)
{
  if(grepl(",",ff)){
    ff = ff %>% strsplit(',') %>% unlist     
  }else{
    ff = Sys.glob(ff)      
  }
  ff
}



create_bins <- function(.,bin_size)
{
  GRanges(seqnames = .$seqnames,
          ranges = IRanges(
            start = seq(1,.$size,by = bin_size),
            width = bin_size
          ))
}

hexbin_scatter_plot <- function(counts,xvar,yvar)
{
  require(viridis)
  require(scales)
  require(ggplot2)
  pal = viridis(1e3, option = "D")

  counts %>%
      ggplot(aes_string(xvar,yvar))+stat_binhex(bins = 140) +
      scale_fill_gradientn(colours = pal,trans = "log10",
                           labels = trans_format('log10',math_format(10^.x)),
                           guide = guide_colourbar(title = NULL,
                                                   barheight = unit(0.92,"npc"),
                                                   barwidth = unit(0.01,"npc")))
}
