
rm(list = ls())

library(tidyverse)
library(viridis)
library(GGally)
library(hexbin)
library(grid)
library(gridExtra)
library(scales)
library(ggrepel)
library(RColorBrewer)

figsdr = "./figs/dRdZ_summary"
datadr = "./data"


all_files = list.files(datadr,recursive = TRUE,pattern = "tsv",full.names = TRUE)

## diff. expressed genes
deseq2_files = file.path(datadr,"Diff.Genes/hg19/DESeq2_strain") %>%
    list.files(full.names = TRUE,recursive = TRUE)

## tpm files
tpm_mat_files = all_files %>%
    {.[grep("TPM",.)]}
tpm_matrices = tpm_mat_files %>%
    map2( c(" ","\t"),.f = function(x,y)read_delim(x,y))


## clean both tpm matrices
des1 = tibble(col = colnames(tpm_matrices[[1]][,-(1:2)])) %>% 
  mutate(
    cell = ifelse(grepl("akata",col),"EBV_dRdZ","EBV-"),
    treat = ifelse(grepl("meth",col),"MC","no"),
    rep = ifelse(grepl(1,col),1,2),
    lab = "Johannsen"
  )

des2 = tibble(col = colnames(tpm_matrices[[2]][,-c(1:2)])) %>% 
  separate(col , into = c("lab","cell","treat","rep"),sep = "\\.",remove = FALSE) %>% 
  mutate(
    cell = ifelse(cell == "EBV","EBV","EBV-"),
    treat = ifelse(treat == "NoTr","no",treat),
    rep = gsub("Rep","",rep) %>% as.numeric(),
    lab = ifelse(lab == "Ben","Johannsen","Scott")
  ) %>% 
  select(col,cell,treat,rep,lab)

tpm_matrices[[2]] = tpm_matrices[[2]] %>% 
  select(-contains("transcrip")) %>% 
  separate(gene_id,into = c("ensembl_id","gene_id"),sep = "\\_")

create_matrix <- function(base,key_treat,key_no)
{
  base %>% 
    mutate(
      treat_avg = base %>% 
        select(contains(key_treat)) %>% 
        rowMeans(),
      notreat_avg = base %>% 
        select(contains(key_no)) %>% 
        rowMeans(),
      log2FC = log2((1 + treat_avg) / (1 + notreat_avg))
    ) %>% 
    select(
      ensembl_id,
      gene_id,
      treat_avg,
      notreat_avg,
      log2FC
    )
}


summary_matrices = list()
summary_matrices[["EBV.dRdZ"]] = create_matrix(
  tpm_matrices[[1]],"meth","no_tr"
)
summary_matrices[["EBV_Ben"]] = create_matrix(
  tpm_matrices[[2]] %>% 
    select(ensembl_id,gene_id,contains("Ben.EBV")),"MC","NoTr"
)
summary_matrices[["NOKS_Ben"]] = create_matrix(
  tpm_matrices[[2]] %>% 
    select(ensembl_id,gene_id,contains("Ben.NOKS")),"MC","NoTr"
)
summary_matrices[["EBV"]] = create_matrix(
  tpm_matrices[[2]] %>% 
    select(ensembl_id,gene_id,contains("Scott.EBV")),"MC","NoTr"
)
summary_matrices[["NOKS"]] = create_matrix(
  tpm_matrices[[2]] %>% 
    select(ensembl_id,gene_id,contains("Scott.NOKS")),"MC","NoTr"
)


create_scatter_component <- function(summary_matrices,x_axis,y_axis,var,minL = 0,maxL = -minL,nbins = 40,
                                     distance_thr = 3)
{
    
    pal = viridis(1e3, option = "D")
    nms = names(summary_matrices)
    names(nms) = nms
    my_data = summary_matrices[c(x_axis,y_axis)] %>%
        map(select_,.dots = c("ensembl_id","gene_id",var)) %>%
        {inner_join(.[[1]],.[[2]],by = c("ensembl_id","gene_id"))}

    cpal = brewer.pal(3,"Set1")
    
    my_genes = locate_genes(my_data,nbins,distance_thr) %>%
        mutate(uplow = factor(
                   if_else(log2FC.x < log2FC.y,"Low","Up"),
                   levels = c("Up","Low")))
    my_data %>%
        ggplot(aes_string(x = paste0(var,".x"),y = paste0(var,".y")))+
        stat_binhex(bins = nbins) +
        scale_fill_gradientn(colours = pal,trans = "log10",
                             labels = trans_format('log10',math_format(10^.x)),
                             guide = guide_colourbar(title = NULL,
                                                     barwidth = unit(0.4,"npc"),
                                                     barheight = unit(0.01,"npc")))+      
        xlim(minL,maxL)+ ylim(minL,maxL)+
        geom_abline(slope = 1,intercept = 0,linetype = 2,colour = "red")+
        xlab(nms[x_axis])+ylab(nms[y_axis])+coord_fixed()+
        theme(legend.position = "top")+
        geom_text_repel(data = my_genes,
                        size = 4,
                        point.padding = unit(1,"lines"),
                        aes(log2FC.x,log2FC.y,label = gene_id,colour = uplow,
                            angle = if_else(my_genes$uplow == "Up", 90,0)),
                        force = 4)+
        scale_color_brewer(palette = "Set1",guide = FALSE)
        
    
}


locate_genes <- function(my_data,nbins,dthr)
{
    
    bindata = hexbin(my_data$log2FC.x , my_data$log2FC.y,xbins = nbins) %>%
        {tibble(
            xcoord = .@xcm,
            ycoord = .@ycm,
            counts = .@count)} %>%
        filter(counts == 1) %>%
        mutate(
            dd = abs(xcoord - ycoord)) %>%
        filter(dd > dthr) %>%
        split( seq_len(nrow(.)))
          
    my_distance <- function(v1,v2)norm(v1  - v2)

    find_nearest <- function(v1, vtib)
    {
        

        tovector <- function(v1,col1,col2)
            v1 %>%
                {c(.[,col1], .[,col2])} %>%
                unlist() %>%
                as.matrix()

        v1 = tovector(v1,1,2)

        vtib = vtib %>%
            mutate(dd = abs(log2FC.x - log2FC.y)) %>%
            filter( dd > dthr )

        distances = vtib %>%
            split(seq_len(nrow(.))) %>%
            map( tovector, 3,4) %>%
            map_dbl( ~ my_distance( v1 , .))

        vtib = vtib %>% mutate(distances)

        vtib %>% arrange(distances) %>% head(1)        
        

    }

    bindata %>%
        map( ~ find_nearest( . , my_data)) %>%
        bind_rows()
    

}    


ll = 10
nbins = 40
dist_thr = 3.4

theme_set(theme_bw())

plots = list()
plots[[1]] = create_scatter_component(summary_matrices,5,4,"log2FC",
                                      minL = -ll,nbins = nbins,distance_thr = 3)+
    xlab("NOKS log2FC")+ylab("EBV log2FC")

plots[[2]] = create_scatter_component(summary_matrices,1,4,"log2FC",
                                      minL = -ll,nbins = nbins,distance_thr = 3)+
    xlab("EBV.dRdZ log2FC")+ylab("EBV log2FC")
plots[[3]] = create_scatter_component(summary_matrices,1,5,"log2FC",
                                      minL = -ll,nbins = nbins,distance_thr = 4)+
    xlab("EBV.dRdZ log2FC")+ylab("EBV- log2FC")


hexbin_scatter = arrangeGrob(
    plots[[1]],
    plots[[2]],
    ggplot()+theme(panel.border = element_blank()),
    plots[[3]],nrow = 2)


ggsave(file = file.path(figsdr,"Log2FC_comparison_EBV_NOKS_EBVdRdZ.pdf"),hexbin_scatter,
       width = 12,
       height = 12,
       units = "in")

strain_test_files = all_files %>%
    {.[grep("strain:",.)]}

strain_tests = strain_test_files %>%
    map(read_tsv)
names(strain_tests) = strain_test_files %>%
    basename() %>% {gsub(".tsv","",.)}

rm(strain_test_files)

drdz_reps = des1 %>%
    split(.$treat) %>%
    map( ~ .$col)

ebv_reps = des2 %>%
    filter(lab == "Scott") %>%
    filter(cell == "EBV") %>%
    split(.$treat) %>%
    map( ~.$col)

noks_reps = des2 %>%
    filter(lab == "Scott") %>%
    filter(cell == "EBV-") %>%
    split(.$treat) %>%
    map( ~ .$col)



average_TPM_boxplot <- function(genes,dRdZ_reps,
                                EBV_reps,tpmmat1,tpmmat2,ll = 200)
{

    tpmmat1 = tpmmat1 %>%
        filter(gene_id %in% genes)
    tpmmat2 = tpmmat2 %>%
        filter(gene_id %in% genes)

    calculate_average <- function(tpm,reps){
        tpm[,reps] %>%
            split(seq_len(nrow(.))) %>%
            map(unlist) %>% map_dbl(mean)}
    dRdZ = calculate_average(tpmmat1,dRdZ_reps)
    EBV = calculate_average(tpmmat2,EBV_reps)

    dt = tpmmat1 %>%
        select(gene_id) %>%
        mutate( dRdZ,EBV)

    dt = dt %>% gather( strain, tpm, -gene_id) %>%
        mutate(strain = factor(strain,levels = c("dRdZ","EBV")))

    dt %>%
        ggplot(aes(strain,tpm,fill = strain))+
        geom_boxplot(outlier.shape = NA) +
        coord_cartesian(ylim = c(0,ll))+
        theme(axis.title.x = element_blank(),
              legend.position = "top",
              axis.text.x = element_text(angle = 25,hjust =1))+
        scale_fill_brewer(palette = "Pastel1")
       
}


pdf(file.path(figsdr,"Boxplot_aveTPM_MC_all_genes.pdf"))
average_TPM_boxplot( strain_tests[[1]]$GENE_ID ,
                    drdz_reps$MC,
                    ebv_reps$MC,
                    tpm_matrices[[1]],
                    tpm_matrices[[2]],ll = 220)+
    ylab("All genes")
dev.off()                    
                    

pdf(file.path(figsdr,"Boxplot_aveTPM_MC_all_diff_exp_genes_EBV_vs_dRdZ.pdf"))
average_TPM_boxplot( strain_tests[[1]] %>%
                     filter(padj <= 1e-6) %>%
                      {.$GENE_ID},
                    drdz_reps$MC,
                    ebv_reps$MC,
                    tpm_matrices[[1]],
                    tpm_matrices[[2]],ll = 220)+
    ylab("Diff. Expressed EBV vs dRdZ genes ")
dev.off()                    


pdf(file.path(figsdr,"Boxplot_aveTPM_MC_all_upreg_diff_exp_genes_EBV_vs_dRdZ.pdf"))
average_TPM_boxplot( strain_tests[[1]] %>%
                     filter(padj <= 1e-6 & log2FoldChange > 0) %>%
                      {.$GENE_ID},
                    drdz_reps$MC,
                    ebv_reps$MC,
                    tpm_matrices[[1]],
                    tpm_matrices[[2]],ll = 220)+
    ylab("Upregulated EBV vs dRdZ genes ")
dev.off()                    


pdf(file.path(figsdr,"Boxplot_aveTPM_MC_all_downreg_diff_exp_genes_EBV_vs_dRdZ.pdf"))
average_TPM_boxplot( strain_tests[[1]] %>%
                     filter(padj <= 1e-6 & log2FoldChange < 0) %>%
                      {.$GENE_ID},
                    drdz_reps$MC,
                    ebv_reps$MC,
                    tpm_matrices[[1]],
                    tpm_matrices[[2]],ll = 220)+
    ylab("Downregulated EBV vs dRdZ genes ")
dev.off()                    
