# Adding DASPER annotation to exon-exon junctionsw
shhh <- suppressPackageStartupMessages
shhh(library(tidyverse))
shhh(library(magrittr))
shhh(library(GenomeInfoDb))
shhh(library(GenomicRanges))
shhh(library(GenomicFeatures))
shhh(library(optparse))
shhh(library(dasper))


arguments <- optparse::parse_args(OptionParser(
  usage = "", 
  description="Script for adding junction annotation to list of junctions of interest. 
  GTF required in position 1. 
  Junction list in position 2. 
  Junction list must contain start, end, strand and chr columns.",
 
  positional_arguments = 2))

  gtf_path<-arguments$args[1]
  junction_list_path<-arguments$args[2]
  
# load junction list based on file format and delim  
  if (grepl(pattern = "\\.csv", x = junction_list_path)) {
    junction_list<-vroom::vroom(file = junction_list_path, delim = ",")
    print("Junction list loaded from csv format")
  } else if(grepl(pattern = "\\.tsv", x = junction_list_path, ignore.case = T)) {
    junction_list<-vroom::vroom(file = junction_list_path, delim = "\t")
    print("Junction list loaded from tsv format")
  } else if(grepl(pattern = "\\.txt", x = junction_list_path, ignore.case = T)) {
    junction_list<-vroom::vroom(file = junction_list_path)
    print("Junction list loaded from guessed delim txt format")
  } else if(grepl(pattern = "\\.rds", x = junction_list_path, ignore.case = T)) {
    junction_list<-readRDS(junction_list_path)
    print("Junction list loaded from RDS/rds format")
  } else {
    print("Junction list not in tsv, csv, txt or RDS/rds format")
  }
  
  
  # TxDb_ref <- dasper:::ref_load(gtf_path)
  TxDb_ref <-GenomicFeatures::makeTxDbFromGFF(gtf_path,format = 'gtf')
  
  # change exon - exon junction to intron - intron for DASPER functions.
  junction_list_gr<-junction_list %>% 
    GenomicRanges::makeGRangesFromDataFrame(start.field = "start",
                                            end.field = "end",
                                            seqnames.field = "chr",
                                            strand.field = "strand",
                                            keep.extra.columns = T)
  
  # convert back into intron-intron junction format for DASPER
  GenomicRanges::end(junction_list_gr) <- GenomicRanges::end(junction_list_gr) - 1
  # 
  GenomicRanges::start(junction_list_gr) <- GenomicRanges::start(junction_list_gr) + 1
  
  GenomeInfoDb::seqlevelsStyle(junction_list_gr) <- seqlevelsStyle(TxDb_ref)
  
  
  annotated_introns <- junction_list_gr %>% dasper::junction_annot(ref = TxDb_ref)
  
  # remove any rows with multiple genes assigned per junction
  annotated_introns <- annotated_introns[which(sapply(annotated_introns$genes, length) < 2),]
  
  # convert back into exon-exon junction format for downstream work
  GenomicRanges::end(annotated_introns) <- GenomicRanges::end(annotated_introns) + 1
  # 
  GenomicRanges::start(annotated_introns) <- GenomicRanges::start(annotated_introns) - 1
  
  annotated_junctions<-annotated_introns %>%  GenomicRanges::as.data.frame() %>% 
    # as_tibble() %>% 
    dplyr::rename('chr'='seqnames') 
  
  # split the annotated junctions into separate output files containing novel donor and novel acceptor types for better parallel processing
  annotated_junctions_donor<-annotated_introns %>% dplyr::filter(type=="novel_donor")
  annotated_junctions_acceptor<-annotated_introns %>%  dplyr::filter(type=="novel_acceptor")
  
  vroom::vroom_write(x = annotated_junctions_donor, file = "annotated_junctions_novel_donors.csv", delim = ",")
  vroom::vroom_write(x = annotated_junctions_acceptor, file = "annotated_junctions_novel_acceptors.csv", delim = ",")
  
  