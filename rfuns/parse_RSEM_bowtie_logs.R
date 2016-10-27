

##' Parses an RSEM log generated when aligning
##' using bowtie
##'
##' @param file File where the logs are stored
##' @return A tibble with the columns
##' processed | aligned | failed | supressed
##'
##' @name parse_RSEM_bowtie_log
##' @export
parse_RSEM_bowtie_log <- function(file)
{
  require(readr)
  require(dplyr)
  require(magrittr)

  logs = read_file(file)
  logs = strsplit(logs,"\n")
  logs = logs[[1]]
  reads = logs[grep('#',logs)]
  reads = strsplit(reads,":",fixed = TRUE)

  reads = lapply(reads,strsplit," ")

  dt = data.frame(t( as.numeric(sapply(reads,function(x)x[[2]][2]))))
  colnames(dt) = c("processed","aligned","failed","supressed")
  
  dt %>% as.tbl
}
