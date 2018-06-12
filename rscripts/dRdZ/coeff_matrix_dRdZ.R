
## Extract the coefficient matrices for the ratio specific tests


library(magrittr)
library(tidyverse)
library(here)
library(rlang)
library(rwlib)
library(DESeq2)


thr = 20

rsem_dr = here("data/RSEM/hg19/Sept17")
align_dr = here("manuscript/logs/Sept17")

alignment = align_dr %>% 
  list.files(full.names = TRUE) %>% 
  stri_subset(".sh") %>% 
  tibble( file = .) %>% 
  mutate(
    stats = map(file,read_tsv,col_names = FALSE,progress = FALSE))

parse_stats = function(stats)
{
  stats %>% 
    filter( 
      str_detect(X1,"total") | str_detect(X1,"mapped")) %>% 
    filter(
      negate(str_detect)(X1,"mate")
      ) %>% 
    mutate(
      type = if_else(str_detect(X1,"total"),"total","aligned"),
      reads = str_split(X1," ") %>% 
        map_chr( ~ .[1]) %>% 
        as.numeric()) %>% 
    dplyr::select(-X1) 
}

alignment %<>% 
  mutate(
    stats = map(stats,parse_stats),
    file = basename(file)   ) %>% 
  unnest() %>% 
  spread(type,reads) %>% 
  mutate(
    perc = aligned / total,
    file = str_replace(file,".logs","")
    ) %>% 
  mutate(
    treatment = file %>% 
      str_detect("methyl") %>% 
      if_else(
        "methyl","none" ),
    cell = case_when(
      str_detect(file,"clone") ~ "dRdZ",
      str_detect(file,"akata") ~ "EBV" ,
      TRUE ~ "NOKS"),
    rep = file %>% 
      str_split("-") %>% 
      map_chr( ~ .[length(.)]) %>% 
      str_replace("clone","") %>% 
      str_replace("rep","") %>% 
      as.integer())

rsem_data = tibble(
  file = list.files(rsem_dr,full.names = TRUE)) %>% 
  mutate(
    rsem = map(file,read_tsv)
    ) %>% 
  mutate(
    file = basename(file) %>% 
      str_replace(".genes.results","")) %>% 
  inner_join(
    select_if(alignment,negate(is.list)), by ="file") %>% 
  dplyr::select(-aligned,-total,-perc)

## remove contaminated samples clone2 and clone4
alignment %<>% filter(!( cell == "dRdZ" & rep %in% c(2,4) ))

rsem_data %<>% 
  filter(!( cell == "dRdZ" & rep %in% c(2,4) ))

as_matrix <- function(x)
{
  x %>% 
  as.data.frame() %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("gene_id") %>% 
  as.matrix()
}

count_matrix = rsem_data %>% 
  dplyr::select(file,rsem) %>% 
  unnest() %>% 
  dplyr::select(file,gene_id,expected_count) %>% 
  mutate(
    expected_count = floor(expected_count)
    ) %>% 
  spread(file,expected_count) %>% 
  as_matrix()

coldata = rsem_data %>% 
  dplyr::select(file,cell,treatment) %>% 
  mutate(interac = paste(cell,treatment, sep = ".")) %>% 
  as.data.frame() %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("file")


deseq = DESeqDataSetFromMatrix(
  count_matrix,colData = coldata,design = ~ interac)

deseq = deseq[ rowSums(assay(deseq) ) > thr,]
deseq = DESeq(deseq)


do_contrast = function(deseq,contrast)
{
  results(deseq,
  cooksCutoff = FALSE,
  contrast = contrast,
  tidy = TRUE) %>% 
  as_tibble() %>% 
  dplyr::rename(
  gene_id = row
  ) %>% 
  mutate(
  log10pval = -log10(pvalue)
  )
}

## The contrast is EBV.none vs NOKS.none

diff_genes = tribble(
  ~ cell, ~ test, ~ contrast,
  "dRdZ", "MC_vs_none", c("interac","dRdZ.methyl","dRdZ.none"),
  "EBV", "MC_vs_none",c("interac","EBV.methyl","EBV.none"),
  "NOKS","MC_vs_none",c("interac","NOKS.methyl","NOKS.none")
  ) %>% 
  mutate(
    results = map(contrast, ~ do_contrast(deseq,.))
    )

common_genes = diff_genes %>% 
  dplyr::select(results) %>% 
  unnest() %>% 
  bind_rows() %>% 
  group_by(gene_id) %>% 
  summarize(
    n = n()
    ) %>% 
  filter(n == 3) %>% 
  ungroup() %>% 
  dplyr::select(gene_id)

rsem_data = rsem_data %>% 
  mutate(
    rsem = map(rsem , inner_join,common_genes,by = "gene_id")
    )

rlog = rlog(deseq)  

rlogmat = rlog %>% 
  assay() %>% 
  as.matrix()

rm(coldata,common_genes,count_matrix,rlog)


ratio_of_ratios_deseq <- function(rsem_data,thr = 20)
{
  count_matrix =  rsem_data %>% 
  dplyr::select(file,rsem) %>% 
  unnest() %>% 
  dplyr::select(file,gene_id,expected_count) %>% 
  mutate(
  expected_count = floor(expected_count)
  ) %>% 
  spread(file,expected_count) %>% 
  as_matrix()
  
  coldata = rsem_data %>% 
  dplyr::select(file,cell,treatment) %>% 
  mutate(interac = paste(cell,treatment, sep = ".")) %>% 
  as.data.frame() %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("file")
  
  deseq = DESeqDataSetFromMatrix(
  count_matrix,colData = coldata,
  design = ~ cell + treatment + cell:treatment)
  
  deseq = deseq[ rowSums(assay(deseq) ) > thr,]
  deseq = DESeq(deseq, test = "LRT", reduced =  ~ cell + treatment)
  
  deseq
  
}

ratio_of_ratios = tribble(
  ~ test , ~ deseq,
  "EBV_vs_NOKS" , rsem_data %>% 
    filter(cell != "dRdZ") %>% 
    ratio_of_ratios_deseq(thr = thr),
  "EBV_vs_dRdZ" , 
  rsem_data %>% 
    filter(cell != "NOKS") %>% 
    ratio_of_ratios_deseq(thr = thr),
  "NOKS_vs_dRdZ",
  rsem_data %>% 
    filter(cell != "EBV") %>% 
    ratio_of_ratios_deseq(thr = thr))

clean_results <- function(deseq)
{
  results(deseq,
          cooksCutoff = FALSE,
          tidy = TRUE) %>% 
    as_tibble() %>% 
    dplyr::rename(
      gene_id = row  ) %>% 
    mutate(
      log10pval = -log10(pvalue)
      )
  
}


outdr = here("Data/Ratio_tests")

ratio_of_ratios %<>% 
  mutate(
    results = map(deseq, clean_results) %>% 
      map(separate,gene_id , into = c("ensembl","symbol"),
          sep = "\\_",extra = "drop" ),
    coeff_mats = map(deseq,coef) %>% map(as.data.frame) %>% 
      map(rownames_to_column,"gene_id") %>% 
      map(dplyr::select,gene_id,everything()) %>% 
      map(as_data_frame),
    outfiles = file.path(outdr, paste0(test,"_cofficient_matrix.tsv")),
    coeff_mats = map2(coeff_mats,outfiles,write_tsv))
    











