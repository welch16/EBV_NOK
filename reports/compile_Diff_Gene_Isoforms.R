#!/usr/bin/env Rscript


## USAGE: TBD
## EXAMPLE: TBD

library(optparse)
library(tools)

option_list = list(
    make_option("--report_name",type = "character",action = "store_true",
                help = "A different name touse for the report file (excluding file extension)"),
    make_option("--gene_file",type = "character",action = "store_true",
                help = "File with the results of the diff. exp. test at gene level to be considered"),
    make_option("--iso_file",type = "character",action = "store_true",
                help = "File with the results of the diff. exp. test at isoform level to be considered"),
    make_option("--test_name",type = "character",action = "store_true",
                help = "String with the name of the test used")
)


opt = parse_args(OptionParser(option_list = option_list))

report_name = opt$report_name
gene_file = opt$gene_file
iso_file = opt$iso_file
test_name = opt$test_name

# report template file
Rmdfile = "reports/Results_Diff_Gene_Isoforms.Rmd"

ext = file_ext(Rmdfile)
new_file = file.path(dirname(Rmdfile),
                     paste(report_name,ext,sep = "."))

file.copy(from = Rmdfile , to = new_file,overwrite = TRUE)

Rmdfile = new_file

message("new report base " ,Rmdfile)

rmarkdown::render(
   input = Rmdfile,
   output_format = "pdf_document",
   params = list(
     gene_file = gene_file,
     iso_file = iso_file,
     test_name = test_name))

