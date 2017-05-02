  
library(shiny)
library(tidyverse)
library(grid)
library(gridExtra)
library(ggrepel)
  
## load gene expression data
tpmmat = read_tsv("data/metadata/RSEM_gene_TPM_matrix.tsv")
  
tpmmat = tpmmat %>% 
  separate(gene_id ,into = c("ensembl","gene"),sep = "\\_") %>%
  dplyr::rename(transcript = `transcript_id(s)`)
  
dipdr = "data/MeDIPseq_results"
dipmat = read_tsv(file.path(dipdr,"MeDIPseq_PromotersCounts_upstr500_downstr1000fraglen300.tsv"))
colnames(dipmat)[1] = "ensembl"
  
dipdepths = read_tsv("data/metadata/MeDIPseq_sequencingDepth.tsv")
  
  
## clean MeDIP-seq, removing input1,...,input3 because we
## pooled them into input0
  
dipdepths = dipdepths %>%
  filter(!grepl("Input-rep1",file),
         !grepl("Input-rep2",file),
         !grepl("Input-rep3",file))

dipmat = dipmat %>% dplyr::select(-contains("Input-rep1"),
                           -contains("Input-rep2"),
                           -contains("Input-rep3"))
  
## remove samples with very low sequencing depths
minDepth = 10e6

base = dipmat %>% dplyr::select(-contains("MeDIP"))

w = dipdepths %>% filter(depth > minDepth) %>% {.$file} %>% 
  {gsub(".sort.bam","",.)}

dipmat = bind_cols(base,dipmat[w]) %>% right_join(tpmmat[,1:2],by = "ensembl") %>%
  dplyr::select(ensembl,gene,seqnames,start,end,description,everything())

genes = intersect(dipmat$gene,tpmmat$gene)
  
tpmmat = tpmmat %>% filter(gene %in% genes) %>% dplyr::select(gene,ensembl,everything())
dipmat = dipmat %>% filter(gene %in% genes) %>% dplyr::select(gene,ensembl,everything())
  
create_boxplot <- function(mat)
{
  ggplot(mat,aes(cell,tpm))+geom_boxplot() +
    geom_point()+geom_text_repel(aes(label = sample),colour = "blue",show.legend = FALSE,
                                 size = 5)+
    theme_minimal()+
    theme(axis.title.x = element_blank(),axis.text = element_text(size = 15),
          axis.title.y = element_text(size = 20))
}

    
ui = fluidPage(
  
  # Application title
  titlePanel("Gene Expression and Gene Methylation"),
  
  # Sidebar with controls to select the random distribution type
  # and number of observations to generate. Note the use of the
  # br() element to introduce extra vertical spacing
  selectInput("searchgene",
              label = "Genes",
              choices = genes,
              selected = genes[1]),
  selectInput("methyltre",
              label = "mC or hmC",
              choices = c("mC","hmC"),
              selected = "hmC"),
  fixedRow(
    tabsetPanel(type = "tabs", 
                tabPanel("Comparison", plotOutput("plot")), 
                tabPanel("Gene Expression", tableOutput("geneExp_table")),
                tabPanel("Methylation",tableOutput("meth_table"))
             )
    )
  )
  # sidebarLayout(
  #   sidebarPanel(
  #     selectInput("searchgene",
  #                 label = "Genes",
  #                 choices = genes,
  #                 selected = genes[1]),
  #   ),
      
    # Show a tabset that includes a plot, summary, and table view
    # of the generated distribution
#     mainPanel(
#       tabsetPanel(type = "tabs", 
#                   tabPanel("Comparison", plotOutput("plot")), 
#                   tabPanel("Gene Expression", tableOutput("geneExp_table")),
#                   tabPanel("Methylation",tableOutput("meth_table"))
#         )
#     )
#   )
# )
#   
  
# Define server logic for random distribution application
server = function(input, output) {
  
  # Reactive expression to generate the requested distribution.
  # This is called whenever the inputs change. The output
  # functions defined below then all use the value computed from
  # this expression
  geneExp <- reactive({
      
    mat = tpmmat %>% filter(gene == input$searchgene) %>%
      dplyr::select(contains("Scott"))
    mat = mat %>% t 
    nms = mat %>% rownames
    mat = tibble(sample = gsub("Scott.","",nms), tpm = mat[,1] ) %>% 
      mutate(cell = ifelse(grepl("EBV",sample),"EBV","NOKS"))
    mat
  })
    
  methyl <- reactive({
    mat  = dipmat %>% filter(gene == input$searchgene) %>%
      dplyr::select(-contains("Input")) %>% dplyr::select(contains("MeDIP"))
    if(input$methyltre == "hmC"){
      mat = mat  %>% dplyr::select(contains("hmC"))
    }else{
      mat = mat %>% dplyr::select(-contains("hmC"))
    }
    
    mat = mat %>% t 
    nms = mat %>% rownames
    mat = tibble(sample = gsub("MeDIPseq-","",nms), tpm = mat[,1] ) %>% 
        mutate(cell = ifelse(grepl("akata",sample),"EBV","NOKS"))
    mat
  })
    
  output$plot <- renderPlot({
    
    grid.arrange(create_boxplot(geneExp()) + ylab("Gene Expression"),
                 create_boxplot(methyl()) + ylab("Methylation"),
                 nrow = 1,top = input$searchgene)
  })
    
  # Generate an HTML table view of the data
  output$geneExp_table <- renderTable({
    geneExp()
  })
    
  output$meth_table <- renderTable({
    methyl()
  })
    
}
  
  
# Run the application 
shinyApp(ui = ui, server = server)
  
