# Creating a .gtf for novel donor junctions

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Setup ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load libraries
# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(readr)
  library(GenomicRanges)
  library(Biostrings)
  library(optparse)
  library(tidyverse)
  library(rtracklayer)
  library(purrr)
})


arguments <- optparse::parse_args(
  OptionParser(
    usage = "", 
    description="Script for adding junction annotation to list of junctions of interest. 
  GTF required in position 1. 
  Fasta required in position 2. 
  Annotated junction list in position 3, that must contain start, end, strand and chr columns.",
    positional_arguments = 3)
)

gtf_path<-arguments$args[1]
fasta_path<-arguments$args[2]
junction_list_path<-arguments$args[3]

conflicted::conflict_prefer_all(winner = "dplyr", quiet = T)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 2.2 novel donor junctions ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load relevant gtf
gtf<-rtracklayer::import(gtf_path)

gtf_df<- GenomicRanges::as.data.frame(gtf)

novel_donor_list<-read_csv(junction_list_path)

novel_donor_gtf<-c()

annotated_donor_gtf<-c()

for(i in 1:length(novel_donor_list$start)){  
  # select row i only
  filt_novel_donor_list<-novel_donor_list[i,]

  junction_start<-as.numeric(novel_donor_list$start)
  junction_end<-as.numeric(novel_donor_list$end)
  
  if (filt_novel_donor_list$strand=="-") { # for negative strand
    # take transcripts the junction acceptor can belong to. Is end since on negative strand donor is exon end and exon start is acceptor
    transcript_acceptors<-purrr::flatten(filt_novel_donor_list$tx_name_start)
    
    # make vector to put results into
    new_transcript_list_donors <- vector("list", length(transcript_acceptors))
    
    for (j in 1:length(transcript_acceptors)) {
      # filter for exons and arrange in order of exon number within transcript
      filt_gtf_df<-gtf_df[gtf_df$transcript_id %in% transcript_acceptors[[j]],] %>%
        filter(type == "exon") %>% dplyr::arrange(as.numeric(exon_number))
      
      # troubleshooting:
      print(paste("j is: ", j))
      print(filt_gtf_df$transcript_id)
      
      filt_gtf_df_ann_ts<-filt_gtf_df
      
      pos_ex<-which(filt_gtf_df$end %in% junction_start)  # find the end of the exon which the acceptor position matches
      
      # if the new donor site goes before the first exon - will move on as won't know where that exon starts.
      if (pos_ex-1<min(as.numeric(filt_gtf_df$exon_number))) { 
        next
      }
      
      if (as.numeric(junction_end)>filt_gtf_df$end[pos_ex-1]) { 
        # if the junction end is before the end of the previous exon then the exon will be 
        #negative width/won't make sense in terms of using the other exon coordinate.
        next
      }
      # go one row/exon back to change the start of the exon to the new donor position
      filt_gtf_df$start[pos_ex-1]<-junction_end
      
      
      new_transcript_gtf<-makeGRangesFromDataFrame(filt_gtf_df, keep.extra.columns = T)
      
      new_transcript_gtf$Transcript<-str_c(new_transcript_gtf$Transcript, "_ND")
      
      annotated_transcript_gtf<-makeGRangesFromDataFrame(filt_gtf_df_ann_ts, keep.extra.columns = T)
      
    }
    
    novel_donor_gtf[i]<-new_transcript_gtf
    
    annotated_donor_gtf[i]<-annotated_transcript_gtf
    
    } else { # for positive strand
    
    
    # take transcripts the junction acceptor can belong to. Is end since on positive strand donor is exon end and exon start is acceptor
    transcript_acceptors<-purrr::flatten(filt_novel_donor_list$tx_name_end)
    
    # make vector to put results into
    new_transcript_list_donors <- vector("list", length(transcript_acceptors))
    
    
    for (j in 1:length(transcript_acceptors)) {
      
      filt_gtf_df<-gtf_df[gtf_df$transcript_id %in% transcript_acceptors[[j]],] %>% filter(type == "exon")
      
      # troubleshooting:
      print(paste("j is: ", j))
      print(filt_gtf_df$transcript_id)
      
      
      filt_gtf_df_ann_ts<-filt_gtf_df
      
      pos_ex<-which(filt_gtf_df$start %in% filt_novel_d_sig_leaf_tdp$end)  # find the end of the junction in the start of a new exon - which row that is
      
      # if the new donor site goes before the first exon - will move on as won't know where that exon starts.
      if (pos_ex-1<min(as.numeric(filt_gtf_df$exon_number))) { 
        next
      }
      
      if (as.numeric(junction_start)<filt_gtf_df$start[pos_ex-1]) { # skips ahead if the donor site is before the start of the previous exon - would cause negative new exon length
        next
      }
      
      filt_gtf_df$end[pos_ex-1]<-junction_start
      
      
      new_transcript_gtf<-makeGRangesFromDataFrame(filt_gtf_df, keep.extra.columns = T)
      
      new_transcript_gtf$Transcript<-str_c(new_transcript_gtf$Transcript, "_ND")
      
      annotated_transcript_gtf<-makeGRangesFromDataFrame(filt_gtf_df_ann_ts, keep.extra.columns = T)
    
  }
  
    novel_donor_gtf[i]<-new_transcript_gtf
    
    annotated_donor_gtf[i]<-annotated_transcript_gtf
    
  }
  
  }

# Combine all results into one dataframe

chunk_num<-gsub(x = basename(junction_list_path), pattern = ".*_(chunk_.).*", replacement="\\1")

write_rds(annotated_donor_gtf, paste0(chunk_num,"_annotated_donor_transcripts_gtf.rds"))

write_rds(novel_donor_gtf, paste0(chunk_num,"_novel_donor_transcripts_gtf.rds"))