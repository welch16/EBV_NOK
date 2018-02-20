
library(tidyverse)
library(DESeq2)
library(here)

rsem_dr = here("data/RSEM/hg19/Sept17")
align_dr = here("manuscript/logs/Sept17")

## alignment = read_tsv("figs/fig2/Scott_aligned_reads_table.tsv")

stri_subset = function(string,pattern)string[negate(str_detect)(string,pattern)]

alignment = align_dr %>% 
  list.files(full.names = TRUE) %>% 
  stri_subset(".sh") %>% 
  tibble(
    file = .
  ) %>% 
  mutate(
    stats = map(file,read_tsv,col_names = FALSE))

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
    select(-X1) 
}

alignment = alignment %>% 
  mutate(
    stats = map(stats,parse_stats),
    file = basename(file)
  ) %>% 
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
        "methyl","none"
      ),
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
      str_replace(".genes.results","")
  ) %>% 
  inner_join(
    select_if(alignment,negate(is.list)), by ="file"
  ) %>% 
  select(-aligned,-total,-perc)

## remove contaminated samples clone2 and clone4
alignment = alignment %>% 
  filter(!( cell == "dRdZ" & rep %in% c(2,4) ))

rsem_data = rsem_data %>% 
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
  select(file,rsem) %>% 
  unnest() %>% 
  select(file,gene_id,expected_count) %>% 
  mutate(
    expected_count = floor(expected_count)
  ) %>% 
  spread(file,expected_count) %>% 
  as_matrix()

coldata = rsem_data %>% 
  select(file,cell,treatment) %>% 
  mutate(interac = paste(cell,treatment, sep = ".")) %>% 
  as.data.frame() %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("file")


## remove genes with low read counts, filtered genes with no reads but probably can 
## filter genes with low expression

deseq = DESeqDataSetFromMatrix(
  count_matrix,colData = coldata,
  design = ~ interac
)

thr = 20

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
  ~ cell , ~ contrast ,
  "EBV",  c("interac","EBV.methyl","EBV.none"),
  "NOKS", c("interac","NOKS.methyl","NOKS.none"),
  "dRdZ", c("interac","dRdZ.methyl", "dRdZ.none")
) %>% 
  mutate(
    results = map(contrast, ~ do_contrast(deseq,.))
  )

common_genes = diff_genes %>% 
  select(results) %>% 
  unnest() %>% 
  bind_rows() %>% 
  group_by(gene_id) %>% 
  summarize(
    n = n()
  ) %>% 
  filter(n == 3) %>% 
  ungroup() %>% 
  select(gene_id)

library(patchwork)

pval_histogram = function(result,cell)
{
  result %>% 
    ggplot(aes(pvalue))+
    geom_histogram(bins = 51,fill = "white",color = "black")+
    ggtitle(
      cell
    )
}

histograms = diff_genes %>% 
  mutate(
    plots = map2(results,cell,pval_histogram)
  ) %>% 
  pull(plots)

histograms[[1]] + histograms[[2]] + histograms[[3]] + plot_layout(nrow = 1)

rsem_data = rsem_data %>% 
  mutate(
    rsem = map(rsem , inner_join,common_genes,by = "gene_id")
  )

message("generating rlog object")
rlog = rlog(deseq)  

rlogmat = rlog %>% 
  assay() %>% 
  as.matrix()

message("saving...")
save(
  alignment,
  rsem_data,
  diff_genes,
  rlogmat,
  file = here("apps/Fig5_dashboard/fig5data.RData"))

