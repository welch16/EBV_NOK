
rm(list = ls())

library(tidyverse)
library(scales)
library(ggbeeswarm)
library(ggjoy)


stri_subset = function(string,pattern)
{
  string[! str_detect(string,pattern)]
  
}

files = list.files("data/RSEM/hg19/Sept17",
                   pattern = "genes",
                   full.names = TRUE) %>% 
  stri_subset("clone[2|4]") %>% 
  set_names(
    basename(.) %>% 
      str_replace(".genes.results","")
  )  


rsem = files %>% 
  map(read_tsv,progress = FALSE) %>% 
  map(select,gene_id,expected_count) %>% 
  tibble(expected_counts = ., 
         file = names(.)) %>% 
  select(file,everything()) %>% 
  unnest()

summaries = rsem %>% 
  group_by(gene_id) %>% 
  summarize(
    median = median(expected_count),
    mean = mean(expected_count)
  )
  
rle = rsem %>% 
  left_join(
    summaries %>% 
      select(gene_id,median),by = "gene_id") %>% 
  mutate(
    rle = log2( expected_count / median)
  ) %>% 
  select(file,gene_id,rle) %>% 
  mutate(
    cell = case_when(
      str_detect(file,"clone") ~ "EBV.dRdZ",
      str_detect(file,"akata") ~ "EBV",
      TRUE ~ "NOKS"),
    treatment = case_when(
      str_detect(file,"meth") ~ "MC",
      TRUE ~ "No treat"
    ),
    file = str_replace(file,"RNAseq-","")) %>% 
  filter(!is.infinite(rle))

figsdr = "figs/Batch"

pdf(file.path(figsdr,"RLE_vs_sample.pdf"),width = 9 ,height = 7)
ggplot(rle,aes(file,rle))+
  geom_boxplot(aes(colour = interaction(cell,treatment)), 
               outlier.size = .3)+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "top",
        axis.title.x = element_blank())+
  geom_hline(yintercept = 0,linetype = 2,colour = "red")+
  scale_color_brewer(palette = "Dark2")
ggplot(rle,aes(file,rle))+
  geom_quasirandom(aes(colour = interaction(cell,treatment)))+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.position = "top")+
  geom_hline(yintercept = 0,linetype = 2,colour = "red")+
  scale_color_brewer(palette = "Dark2")
dev.off()




rle %>% 
  ggplot(aes(x = rle, y = file))+
  geom_joy(aes(fill = interaction(cell,treatment)))+
  geom_vline(xintercept = 0,linetype = 2,colour = "black")+
  scale_fill_brewer(palette = "Pastel1",
                    name = "")+
  scale_x_continuous(limits = c(-4,4))+
  theme_joy()+
  theme(axis.title.y = element_blank(),
        legend.position = "top")+
  scale_y_discrete(
    labels = sort(unique(rle$file)) %>% 
      rev()
  )+
  xlab(
    expression(log[2](Counts[ig] / median[g]))
  )


ggsave(file.path(figsdr,"RLE_density_by_sample.pdf"),
       width = unit(8,"in"),
       height = unit(9,"in"))

