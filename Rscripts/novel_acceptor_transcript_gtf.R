# Creating a .gtf for novel acceptor junctions

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Setup ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
## 2.1 novel acceptor junctions ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load relevant gtf
gtf<-rtracklayer::import(gtf_path)

gtf_df<- GenomicRanges::as.data.frame(gtf)

novel_acceptor_list<-read_csv(junction_list_path)

novel_acceptor_gtf<-c()

annotated_acceptor_gtf<-c()


for(i in 1:length(novel_acceptor_list$start)){
  # select row i only
  filt_novel_acceptor_list<-novel_acceptor_list[i,]
  
  junction_start<-filt_novel_acceptor_list$start
  junction_end<-filt_novel_acceptor_list$end
  
if (filt_novel_acceptor_list$strand=="-") { # for negative strand
    
    # take transcripts the junction donor can belong to. Is end since on negative strand donor is exon end and exon start is acceptor  
    transcript_donors<-purrr::flatten(filt_novel_acceptor_list$tx_name_end)
    
    new_transcript_list_acceptors <- vector("list", length(transcript_donors))
    
    # troubleshooting:
    print(paste("i is: ", i))
    
    # loop across each transcript the donor site is found in  
    for (j in 1:length(transcript_donors)) {
      
      # get gtf of transcript that contains donor site that matches annotation
      filt_gtf_df<-gtf_df[gtf_df$transcript_id %in% transcript_donors[[j]],] %>% filter(type == "exon") %>% 
        dplyr::arrange(as.numeric(exon_number))
      
      print(paste("j is ", j))
      print(filt_gtf_df$transcript_id)
      
      # keep a gtf granges for canonical/annotated transcript
      filt_gtf_df_ann_ts<-filt_gtf_df
      
      # find the end of the junction in the start of a new exon - which row that is - find donor site that matches annotation
      pos_ex<-which(filt_gtf_df$start %in%  junction_end)  
      
      # if the new acceptor site goes beyond the last exon - will move on as won't know where that exon ends.
      if (pos_ex+1>max(as.numeric(filt_gtf_df$exon_number))) { 
        next
      }
      
      # if junction start is before the start of the next exon then junction will be negative so avoid and move on as this would not work as replacing the next exon end.
      # If so suggests more complicated splicing that just replacing the next exon end needed.
      if (as.numeric(junction_start)<filt_gtf_df$start[pos_ex+1]) {
        next
      }
      
      # change the next exon end to the start of the novel junction / acceptor site
      filt_gtf_df$end[pos_ex+1]<-as.numeric(junction_start)
      
      # run get sequences on new transcript with novel acceptor site, to get NT, CDS and AA sequences (AA sequence is longest possible one based on possible start sites)
      new_transcript_gtf<-makeGRangesFromDataFrame(filt_gtf_df, keep.extra.columns = T)
      
      # change transcript name to reflect additional of novel acceptor
      new_transcript_gtf$transcript_name<-str_c(new_transcript$transcript_name, "_NA")
      # run get sequences on normal/annotated transcript to get NT, CDS and AA sequences (AA sequence is longest possible one based on possible start sites)
      annotated_transcript_gtf<-makeGRangesFromDataFrame(filt_gtf_df_ann_ts, keep.extra.columns = T)
      
      
  }
    
    novel_acceptor_gtf[i]<-new_transcript_gtf
    
    annotated_acceptor_gtf[i]<-annotated_transcript_gtf    
    
    
    } else { # for positive strand
    
    # take transcripts the junction donor can belong to. Is end since on positive strand donor is exon end and exon start is acceptor  
    transcript_donors<-purrr::flatten(filt_novel_acceptor_list$tx_name_start)
    
    new_transcript_list_acceptors <- vector("list", length(transcript_donors))
    
    for (j in 1:length(transcript_donors)) {
      
      # get gtf of transcript that contains donor site that matches annotation  
      filt_gtf_df<-gtf_df[gtf_df$transcript_id %in% transcript_donors[[j]],] %>% filter(type == "exon")
      
      # troubleshooting:
      print(paste("j is ", j))
      print(filt_gtf_df$transcript_id)
      
      # keep a gtf granges for canonical/annotated transcript
      filt_gtf_df_ann_ts<-filt_gtf_df
      
      # find the start of the junction in the end of an exon - which row that is - find donor site that matches annotation
      pos_ex<-which(filt_gtf_df$end %in% junction_start)  
      
      # if the new acceptor site goes beyond the last exon - will move on as won't know where that exon ends.
      if (pos_ex+1>max(as.numeric(filt_gtf_df$exon_number))) { 
        next
      }
      
      # if junction end / novel exon start goes beyond the end of the next exon, move on as this would not work as replacing the next exon start.
      # If so suggests more complicated splicing that just replacing the next exon start needed.
      if (as.numeric(junction_end)>filt_gtf_df$end[pos_ex+1]) { 
        next
      }
      
      # change the start of the next exon to the end of the junction / new acceptor
      filt_gtf_df$start[pos_ex+1]<-as.numeric(junction_end)
      
      new_transcript<-makeGRangesFromDataFrame(filt_gtf_df, keep.extra.columns = T)
      

      new_transcript$transcript_name<-str_c(new_transcript$transcript_name, "_NA")
      
      annotated_transcript<-makeGRangesFromDataFrame(filt_gtf_df_ann_ts, keep.extra.columns = T)
    
  
    }
    
    novel_acceptor_gtf[i]<-new_transcript_gtf
    
    annotated_acceptor_gtf[i]<-annotated_transcript_gtf    
    
  }  
    
    
}

chunk_num<-gsub(x = basename(junction_list_path), pattern = ".*_(chunk_.).*", replacement="\\1")

write_rds(annotated_acceptor_gtf, paste0(chunk_num, "_annotated_acceptor_transcripts_gtf.rds"))

write_rds(novel_acceptor_gtf, paste0(chunk_num,"_novel_acceptor_transcripts_gtf.rds"))

