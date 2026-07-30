library(tidyverse)
library(magrittr)

arguments <- optparse::parse_args(OptionParser(
  usage = "", 
  description="Script for splitting annotated junctions lists into smaller chunks for parallel processing.
  Junction list path in position 1. 
  Junction list must contain start, end, strand and chr columns.",
  
  positional_arguments = 1))

annotated_junction_list<-read_csv(arguments$args[1])

chunks <- split(junction_list, ceiling(seq_len(nrow(junction_list)) / 10000))

str(chunks)

names(chunks)

# if novel donor or acceptor

if (grepl(pattern = "acceptor", x = arguments$args[1])) {
  

for (i in names(chunks)) {
  write_csv(x = chunks[i], file=paste0("annotated_junctions_novel_acceptor_chunk_",i, ".csv"))  

  }
  
} else if (grepl(pattern = "donor", x = arguments$args[1])) {
  
  for (i in names(chunks)) {
    write_csv(x = chunks[i], file=paste0("annotated_junctions_novel_donor_chunk_",i, ".csv"))  
    
  }  
  
} else{
  cat(paste0("error with this file", arguments$args[1]), file = "error_annotated_junction_chunk_", arguments$args[1], "_", i, ".txt")
}
