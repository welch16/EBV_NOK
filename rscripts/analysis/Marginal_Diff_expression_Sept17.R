
rm(list = ls())

library(tidyverse)
library(scales)

indr = "data/Diff.Genes/hg19/Sept17"
figsdr = "figs/diff_expression/Sept17"
rsemdr = "data/RSEM/hg19/Sept17"

files = list.files(indr,full.names = TRUE)

meta = tibble(
  file = files %>% basename()
) %>% 
  mutate(
    cell = if_else(grepl("dRdZ",file),"EBV_dRdZ",
                   if_else(grepl("NOKS",file),"NOKS","EBV")),
    unit = if_else(grepl("genes",file),"gene","isoform"),
    file = map_chr(file , ~ gsub(paste0(".",tools::file_ext(.)),"",.))
  ) 
  

diff_genes = files %>% 
  set_names(meta$file) %>% 
  map(read_tsv) %>% 
  map2( names(.), ~ mutate(.x,file = .y))


theme_set(theme_bw())

## p.value histograms
pdf(file.path(figsdr,"Pvalue_histogram_all.pdf"),width = 7 , height = 6)
diff_genes %>% 
  map(select,file,pvalue) %>% 
  bind_rows() %>% 
  inner_join(meta,by = "file") %>% 
  ggplot(aes(pvalue))+
  geom_histogram(boundary = 0,colour = "blue",fill = "white",bins = 31)+
  facet_grid( unit ~ cell , scales = "free_y")
dev.off()

fdr = 1e-10
## volcano plots
pdf(file.path(figsdr,"Volcano_plot_all.pdf"),width = 7, height = 6)
diff_genes %>% 
  map(select,file,log2FC,padj,pvalue) %>% 
  bind_rows() %>% 
  inner_join(meta, by = "file") %>%
  mutate(log10pval = -log10(pvalue),
         DE = if_else(padj <= fdr,"yes","no")) %>% 
  filter(!is.na(DE)) %>% 
  ggplot(aes(log2FC,log10pval,colour = DE))+
  geom_point(alpha = 1/2)+
  geom_vline(xintercept = 0,linetype = 2)+
  facet_grid(unit  ~ cell)+ylim(0,30)+
  theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = paste0("DE: padj <=",fdr))+
  ylab(expression(-log[10](p.value)))
dev.off()
  
## pheatmaps for each test top K
options(mc.cores = 16)

rsem_files = rsemdr %>%
    list.files(full.names = TRUE,pattern = "results")

rsem_meta = tibble(
    file = rsem_files %>% basename()) %>%
    mutate(
    cell = if_else(grepl("clone",file),"EBV_dRdZ",
                   if_else(grepl("akata",file),"EBV","NOKS")),
    unit = if_else(grepl("genes",file),"gene","isoform")) %>% 
    split(.$cell)

library(tximport)

gene_files = rsem_meta %>%
    map(filter,unit == "gene") %>%
    map( ~ file.path(rsemdr,.$file))

gene_data = gene_files %>%
    map( tximport,type = "rsem",importer =read_tsv)

iso_files = rsem_meta %>%
    map(filter,unit == "isoform") %>%
    map(~ file.path(rsemdr,.$file))

iso_data = iso_files %>%
    map(tximport,type = "rsem",importer = read_tsv) 

gene_coldata = gene_files %>%
    map(basename) %>%
    map( ~ {
        df = data.frame(file = . ) %>%
            mutate(treat = if_else(grepl("methyl",file),"Methyl_Cell","No_Treatment"))
        rownames(df) = gsub(".genes.results","",df$file)
        df %>%
            mutate(file = NULL)
    })

iso_coldata = iso_files %>%
    map(basename) %>%
    map( ~ {
        df = data.frame(file = . ) %>%
            mutate(treat = if_else(grepl("methyl",file),"Methyl_Cell","No_Treatment"))
        rownames(df) = gsub(".genes.results","",df$file)
        df %>%
            mutate(file = NULL)
    })

library(DESeq2)

deseq_rl <- function(tpm_data,col_data)
{
    deseq = DESeqDataSetFromMatrix(
        floor(tpm_data[["counts"]]),
        colData = col_data,
        design = ~ treat)
    deseq = deseq[rowSums(counts(deseq)) > 1,]

    rlog(deseq,blind = FALSE)

}

gene_rlog = map2(gene_data,
                 gene_coldata,
                 deseq_rl)

iso_rlog = map2(iso_data,
                iso_coldata,
                deseq_rl)

library(pheatmap)
pal = viridis::viridis(1e3)

rowVars <- function (x,na.rm = TRUE)
{
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}

generate_heatmaps <- function(rld,name,Ks = c(20,50))
{

    coldata = colData(rld)
    coldata$sizeFactor = NULL

    coldata = as.data.frame(coldata)
    
    aux_heatmap <- function(K,rld,coldata,name){
        mat = assay(rld)[,]
        mat = mat - rowMeans(mat)

        topgenes = head(order(rowVars(assay(rld)),decreasing = TRUE),K)
        pheatmap(mat[topgenes,],
                 annotation_col = coldata,color = pal,
                 show_rownames = FALSE,
                 filename = file.path(figsdr,paste0(name,"_Heatmap_top",K,".pdf")))
    }

    Ks %>% map(aux_heatmap,
               rld,coldata,name)


    
}

Ks = c(20,50,100,500,1e3,5e3)

gene_rlog %>%
    map2(paste(names(.),"Genes",sep = "_"),generate_heatmaps,K = Ks)

iso_rlog %>%
    map2(paste(names(.),"Isoforms",sep = "_"),generate_heatmaps,K = Ks)




