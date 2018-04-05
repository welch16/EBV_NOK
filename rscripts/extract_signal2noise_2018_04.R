
library(tidyverse)
library(here)
library(glue)

outdr = here("data/signal2noise/")

load2env = function(file,env = new.env()){load(file,envir = env);env}

app_data = here("apps/") %>% 
  list.files(full.names = TRUE,recursive = TRUE) %>% 
  str_subset("RData") %>% 
  map(load2env) %>% 
  set_names(c("fig2","fig3","fig5"))

process_rsem = function(data)
{
  data %>% 
    select(file,cell,treatment,rep,rsem) %>% 
    unnest() %>% 
    select(-rep) %>% 
    group_by(cell,treatment,gene_id) %>% 
    summarize_if(is.numeric,funs(mean)) %>% 
    ungroup() %>% 
    select(-contains("length")) %>% 
    mutate(expected_count = as.integer(floor(expected_count))) %>% 
    nest(-cell,-treatment,.key = "rsem")
}

join_rsem = function(fig,fig_rsem)
{
  fig %>% 
    inner_join(fig_rsem,by = "cell") %>% 
    spread(treatment,rsem) %>% 
    mutate(
      data = map2(data,none,inner_join,by = "gene_id") %>% 
        map2(MC,inner_join,by = "gene_id",suffix = c(".none",".MC"))
    ) %>% 
    select(cell,data) %>% 
    mutate(
      data = map(data,separate,gene_id,into = c("ensembl","symbol"),sep = "\\_") %>% 
        map(select,-ensembl)
    )

}
## Fig 3

fig3 = tribble( ~ cell ,"EBV","NOKS") %>% 
  mutate(
    data = map(cell, ~ pluck(app_data$fig3,paste0("diff_genes_",.))) %>% 
      map(select,gene_id,log2FoldChange,stat,padj)
  )

fig3_rsem = app_data$fig3 %>% 
  pluck("rsem_data") %>% 
  mutate(
    rep = map(file,str_split,"rep") %>% 
      map(unlist) %>% 
      map_chr( ~ .[length(.)]) %>% 
      as.integer()
  ) %>% process_rsem()


fig3 = fig3 %>% join_rsem(fig3_rsem) %>% 
  mutate(
    file = glue("{outdr}{fig}_{cell}_signal2noise.tsv",
                outdr = outdr,
                cell = cell,
                fig = "fig3")
  )

fig3 %>% 
{
  map2(.$data,.$file,write_tsv)
}


## Fig 5

fig5 = diff_genes = pluck(app_data$fig5,"diff_genes") %>% 
  select(-contrast) %>% 
  rename(data = results) %>% 
  mutate(
    data = data %>% 
      map(select,gene_id,log2FoldChange,stat,padj)
  )

fig5_rsem = app_data$fig5 %>% 
  pluck("rsem_data") %>% 
  mutate(treatment = if_else(treatment == "methyl","MC","none")) %>% 
  process_rsem()

fig5 = fig5 %>% join_rsem(fig5_rsem) %>% 
  mutate(
    file = glue("{outdr}{fig}_{cell}_signal2noise.tsv",
                outdr = outdr,
                cell = cell,
                fig = "fig5")
  ) 

fig5 %>% 
  {
    map2(.$data,.$file,write_tsv)
  }





