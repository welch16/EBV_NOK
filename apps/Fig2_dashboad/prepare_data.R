
library(tidyverse)
library(DESeq2)

dr = "data/RSEM/hg19"
alignment = read_tsv("figs/fig2/Scott_aligned_reads_table.tsv")


rsem_data = tibble(
  file = list.files(dr,full.names = TRUE)) %>% 
  mutate(
    genes = if_else(str_detect(file,"genes"),"genes","isoforms"),
    lab = if_else(str_detect(file,"Noks"),"Scott","Johansenn")) %>% 
  filter( lab == "Scott") %>% 
  mutate(
    cell = if_else(str_detect(file,"EBV"),"EBV","NOKS"),
    treatment = if_else(str_detect(file,"MC"),"MC","none"),
    rsem = map(file,read_tsv),
    genes = NULL,
    lab = NULL,
    file = basename(file) %>% 
      str_replace(".genes.results","")
  ) %>% 
  filter( ! (str_detect(file , "rep1")  & cell == "EBV" & treatment == "none")) %>% 
  mutate(
    rsem = map(rsem, ~ select(.,-contains("transcript"))) ## didn't use transcript_id(s)
  )
  ## removed the EBV-no_treatment-rep1 sample due to bad quality, 
  ## in the alignment plot, it can be see that both the number of 
  ## aligned reads and the percentage of aligned reada is quite bad

  ## a second note: NOKS-MC-rep3 could be treated as another file that looks as that

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

deseq = DESeqDataSetFromMatrix(
  count_matrix,colData = coldata,
  design = ~ interac
)

## remove genes with low read counts, filtered genes with no reads but probably can 
## filter genes with low expression

thr = 20

deseq = deseq[ rowSums(assay(deseq) ) > thr,]
deseq = DESeq(deseq)

## The contrast is EBV.none vs NOKS.none
diff_genes = results(
  deseq , 
  cooksCutoff = FALSE,
  contrast = c("interac","EBV.none","NOKS.none"),tidy = TRUE) %>% ## 
  as_tibble() %>% 
  dplyr::rename(
    gene_id = row
  )

diff_genes %>%
  ggplot(aes(pvalue))+
  ggtitle(thr)+
  geom_histogram(bins = 51,fill = "white",color = "black")

rsem_data = rsem_data %>% 
  mutate(
    rsem = map(rsem , ~ filter(., gene_id %in% pluck(diff_genes,"gene_id")) )
  )

message("generating rlog object")
rlog = rlog(deseq)  
rlogmat = rlog %>% 
  assay() %>% 
  as.matrix()

message("saving...")
save(alignment,
     rsem_data,
     diff_genes,
     rlogmat,
     file = "./apps/Fig2_dashboad/fig2data.RData")
