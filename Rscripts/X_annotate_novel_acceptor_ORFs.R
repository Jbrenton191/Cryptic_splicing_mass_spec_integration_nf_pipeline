#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Setup ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(readr)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(ORFik)
  library(optparse)
  library(tidyverse)
  library(rtracklayer)
  library(purrr)
  library(parallel)
  library(ORFik)
  library(factR)
  library(optparse)
})


arguments <- optparse::parse_args(OptionParser(
  usage = "", 
  description="Script for adding junction annotation to list of junctions of interest. 
  GTF required in position 1. 
  Fasta required in position 2. 
  Annotated junction list in position 3, that must contain start, end, strand and chr columns.",
  
  positional_arguments = 3))

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

# importing factR way for GTF and FASTA
ref.gtf <- factR::importGTF(gtf_path)

fasta <- factR::importFASTA(fasta_path)

# Get sequences functions - from Emil - now with strand selection to reverse seqeunce - added in lines to remove N sequences.
get_sequences <- function(gtf){
  NT_seq <-
    subset(gtf, type == "exon") %>%
    split(., mcols(.)$transcript_id) %>%
    lapply(., function(x) {
      if (unique(strand(x)) == "-") {
        x %>%
          sort(., decreasing = TRUE) %>%  # Sort in decreasing order, which changes start and end of Iranges to make it correct for converting to NTs
          getSeq(Hsapiens, .) %>% paste(., collapse = "")
      } else {
        x %>%
          sort(.) %>%  # Sort in increasing order
          getSeq(Hsapiens, .) %>% paste(., collapse = "")
      }
    }) %>%
    utils::stack() %>%
    setNames(c("NT_seq", "Transcript")) %>%
    subset(., select=c("Transcript", "NT_seq"))
  
  if (str_detect(string = NT_seq$NT_seq,
                 pattern="N")) {
    NT_seq$CDS_seq <- NA
    NT_seq$AA_seq <- NA
    return(NT_seq)
  }
  
  NT_seq$CDS_seq <- NA
  NT_seq$AA_seq <- NA
  
  for (i in 1:NROW(NT_seq)) {
    ORF_range <- ORFik::findORFs(
      NT_seq$NT_seq[i], longestORF = TRUE, startCodon = "ATG") %>%
      unlist() %>%
      data.frame()
    
    # Check if the resulting data frame is empty and adjust accordingly
    if (nrow(ORF_range) == 0) {
      ORF_range <- data.frame(width = NA)
    } else {
      ORF_range <- ORF_range %>%
        dplyr::filter(width == max(width))
    }
    
    # Use ifelse to handle the condition and assign values accordingly
    NT_seq$CDS_seq[i] <- ifelse(is.na(ORF_range$width),
                                NA,
                                substr(as.character(NT_seq$NT_seq[i]),
                                       start = ORF_range$start, stop = ORF_range$end))
    
    # Get AA
    NT_seq$AA_seq[i] <- ifelse(is.na(ORF_range$width),
                               NA,
                               NT_seq$CDS_seq[i] %>% 
                                 DNAString() %>% 
                                 Biostrings::translate() %>% 
                                 as.character())
  }
  
  
  return(NT_seq)
}



# make blank df to combine rows into
# run across each row of the novel sig tdp dataframe
# for (i in 1:length(novel_a_sig_leaf_tdp$start)) {
results_list <- mclapply(1:length(novel_a_sig_leaf_tdp$start), function(i) {  
  # select row i only
  filt_novel_a_sig_leaf_tdp<-novel_a_sig_leaf_tdp[i,]
  
  junction_start<-filt_novel_a_sig_leaf_tdp$start
  junction_end<-filt_novel_a_sig_leaf_tdp$end
  
  if (filt_novel_a_sig_leaf_tdp$strand=="-") { # for negative strand
    
    # take transcripts the junction donor can belong to. Is end since on negative strand donor is exon end and exon start is acceptor  
    transcript_donors<-purrr::flatten(filt_novel_a_sig_leaf_tdp$tx_name_end)
    
    new_transcript_list_acceptors <- vector("list", length(transcript_donors))
    
    # troubleshooting:
    print(paste("i is: ", i))
    
    # loop across each transcript the donor site is found in  
    for (j in 1:length(transcript_donors)) {
      
      # get gtf of transcript that contains donor site that matches annotation
      filt_gtf_df<-gtf_df[gtf_df$transcript_id %in% transcript_donors[[j]],] %>% filter(type == "exon") %>% dplyr::arrange(as.numeric(exon_number))
      
      print(paste("j is ", j))
      print(filt_gtf_df$transcript_id)
      
      # keep a gtf granges for canonical/annotated transcript
      filt_gtf_df_norm_ts<-filt_gtf_df
      
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
      new_transcript<-get_sequences(makeGRangesFromDataFrame(filt_gtf_df, keep.extra.columns = T))
      
      # if there is no open reading frame possible then move to next iteration of loop and also skips any N containing NT sequences
      if (is.na(new_transcript$AA_seq)) {
        next
      }
      
      # change transcript name to reflect additional of novel acceptor
      new_transcript$Transcript<-str_c(new_transcript$Transcript, "_NA")
      # run get sequences on normal/annotated transcript to get NT, CDS and AA sequences (AA sequence is longest possible one based on possible start sites)
      normal_transcript<-get_sequences(makeGRangesFromDataFrame(filt_gtf_df_norm_ts, keep.extra.columns = T))
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      ### NMD prediction! ----  
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      # get transcript ID to add back in to make GTF - transcript id needed to build CDS
      transcript_df<-gtf_df[gtf_df$transcript_id %in% transcript_donors[[j]],] %>% filter(type == "transcript") 
      # turn new df with updated junction into a GRanges object/GTF
      filt_gtf_df_gr<-makeGRangesFromDataFrame(rbind(filt_gtf_df, transcript_df), keep.extra.columns = T) %>% 
        subset(type == "exon" | type == "transcript") 
      # match gene indo
      custom_matched_1.gtf <- matchGeneInfo(filt_gtf_df_gr, ref.gtf)  
      # Build CDS and predict NMD
      custom_new.gtf <- subsetNewTranscripts(custom_matched_1.gtf, ref.gtf)  
      custom_new_CDS.gtf <- buildCDS(query = custom_new.gtf, 
                                     ref = ref.gtf, fasta = fasta)  
      
      if ("CDS" %in% custom_new_CDS.gtf$type) {
        NMDresult1 <- predictNMD(custom_new_CDS.gtf)
        NMDresult <- NMDresult1$is_NMD
      } else {
        message("Skipping NMD prediction: no CDS found.")
        NMDresult <- NA
      }
      
      # NMDprediction.out <- predictNMD(custom_new_CDS.gtf)
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      ### Add in NT seq difference! ----  
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      difference_new_to_old_NT<-!str_split(new_transcript$NT_seq, "")[[1]] == 
        str_split(normal_transcript$NT_seq, "")[[1]]
      
      # find first difference in AA sequence between new and old - if there is a difference, add that to the result, if not state the longest AA seq would produce the same protein
      if (any(difference_new_to_old_NT=="TRUE")) {
        # find first mismatch / TRUE value
        first_difference_new_to_old_NT<-which(difference_new_to_old_NT)[1]
        # take sequence from first mismatch and take all AAs from that start point
        NT_seq_from_first_diff<-str_sub(string = new_transcript$NT_seq, start = first_difference_new_to_old_NT, end=length(difference_new_to_old_NT))
        # add this to result
        new_transcript$NT_seq_from_first_diff<-NT_seq_from_first_diff
        
      } else {
        # if no difference - state that
        new_transcript$NT_seq_from_first_diff<-"No difference in NT sequence between normal and novel ones"
        
      }
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      
      # if there is no open reading frame for normal transcript then set to first value in new_transcript being TRUE
      if (is.na(normal_transcript$AA_seq)) {
        difference_new_to_old<-TRUE
      } else {
        
        # compare longest possible AA sequences to see if there is a deviation from cannonical / annotated protein produced
        # TRUE indicates a difference
        difference_new_to_old<-!str_split(new_transcript$AA_seq, "")[[1]] == 
          str_split(normal_transcript$AA_seq, "")[[1]]
      }
      
      
      # find first difference in AA sequence between new and old - if there is a difference, add that to the result, if not state the longest AA seq would produce the same protein
      if (any(difference_new_to_old=="TRUE")) {
        # find first mismatch / TRUE value
        first_difference_new_to_old<-which(difference_new_to_old)[1]
        # take sequence from first mismatch
        seq_aa_from_first_diff<-str_sub(string = new_transcript$AA_seq, start = first_difference_new_to_old, 
                                        end=length(difference_new_to_old))
        # add this to result
        new_transcript$seq_aa_from_first_diff<-seq_aa_from_first_diff
        
      } else {
        # if no difference - state that
        new_transcript$seq_aa_from_first_diff<-"No difference in longest transcript between normal and novel ones"
        
      }
      
      # add annotated transcript NT seq
      new_transcript$old_transcript_NT_seq<-normal_transcript$NT_seq
      # add annotated transcript AA seq
      new_transcript$old_transcript_CDS_seq<-normal_transcript$CDS_seq
      # add annotated transcript AA seq
      new_transcript$old_transcript_AA_seq<-normal_transcript$AA_seq
      
      
      # add novel junction type to result
      new_transcript$novel_junction_type<-'novel_acceptor'
      # add gene novel junction comes from to result
      new_transcript$genes<-filt_novel_a_sig_leaf_tdp$genes
      # add leafcutter LRRK2 analysis deltapsi of junction to result
      new_transcript$novel_junction_dpsi<-filt_novel_a_sig_leaf_tdp$deltapsi
      # add tdp kd deltapsi of junction to result
      new_transcript$novel_junction_dpsi_tdp_kd_avg<-filt_novel_a_sig_leaf_tdp$tdp_kd_average_deltapsi
      
      
      #add in NMD prediction
      new_transcript$NMD<-NMDresult
      
      # add in TSL support level
      filt_gtf_df_gr2<-makeGRangesFromDataFrame(rbind(filt_gtf_df, transcript_df), keep.extra.columns = T) %>% 
        subset( type == "transcript") 
      
      new_transcript$TSL<-unique(filt_gtf_df_gr2@elementMetadata$transcript_support_level)
      
      
      new_transcript_list_acceptors[[j]]<-new_transcript
      
    }
    
    return(new_transcript_list_acceptors)
    
    
  } else { # for positive strand
    
    
    # take transcripts the junction donor can belong to. Is end since on positive strand donor is exon end and exon start is acceptor  
    transcript_donors<-purrr::flatten(filt_novel_a_sig_leaf_tdp$tx_name_start)
    
    new_transcript_list_acceptors <- vector("list", length(transcript_donors))
    
    for (j in 1:length(transcript_donors)) {
      
      # get gtf of transcript that contains donor site that matches annotation  
      filt_gtf_df<-gtf_df[gtf_df$transcript_id %in% transcript_donors[[j]],] %>% filter(type == "exon")
      
      # troubleshooting:
      print(paste("j is ", j))
      print(filt_gtf_df$transcript_id)
      
      
      # keep a gtf granges for canonical/annotated transcript
      filt_gtf_df_norm_ts<-filt_gtf_df
      
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
      
      new_transcript<-get_sequences(makeGRangesFromDataFrame(filt_gtf_df, keep.extra.columns = T))
      
      # if there is no open reading frame possible then move to next iteration of loop
      if (is.na(new_transcript$AA_seq)) {
        next
      }
      
      new_transcript$Transcript<-str_c(new_transcript$Transcript, "_NA")
      
      normal_transcript<-get_sequences(makeGRangesFromDataFrame(filt_gtf_df_norm_ts, keep.extra.columns = T))
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      ### NMD prediction! ----  
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      # get transcript ID to add back in to make GTF - transcript id needed to build CDS
      transcript_df<-gtf_df[gtf_df$transcript_id %in% transcript_donors[[j]],] %>% filter(type == "transcript") 
      # turn new df with updated junction into a GRanges object/GTF
      filt_gtf_df_gr<-makeGRangesFromDataFrame(rbind(filt_gtf_df, transcript_df), keep.extra.columns = T) %>% 
        subset(type == "exon" | type == "transcript") 
      # match gene indo
      custom_matched_1.gtf <- matchGeneInfo(filt_gtf_df_gr, ref.gtf)  
      # Build CDS and predict NMD
      custom_new.gtf <- subsetNewTranscripts(custom_matched_1.gtf, ref.gtf)  
      custom_new_CDS.gtf <- buildCDS(query = custom_new.gtf, 
                                     ref = ref.gtf, fasta = fasta)  
      
      if ("CDS" %in% custom_new_CDS.gtf$type) {
        NMDresult1 <- predictNMD(custom_new_CDS.gtf)
        NMDresult <- NMDresult1$is_NMD
      } else {
        message("Skipping NMD prediction: no CDS found.")
        NMDresult <- NA
      }
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      ### Add in NT seq difference! ----  
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      difference_new_to_old_NT<-!str_split(new_transcript$NT_seq, "")[[1]] == 
        str_split(normal_transcript$NT_seq, "")[[1]]
      
      # find first difference in AA sequence between new and old - if there is a difference, add that to the result, if not state the longest AA seq would produce the same protein
      if (any(difference_new_to_old_NT=="TRUE")) {
        # find first mismatch / TRUE value
        first_difference_new_to_old_NT<-which(difference_new_to_old_NT)[1]
        # take sequence from first mismatch and take all AAs from that start point
        NT_seq_from_first_diff<-str_sub(string = new_transcript$NT_seq, start = first_difference_new_to_old_NT, end=length(difference_new_to_old_NT))
        # add this to result
        new_transcript$NT_seq_from_first_diff<-NT_seq_from_first_diff
        
      } else {
        # if no difference - state that
        new_transcript$NT_seq_from_first_diff<-"No difference in NT sequence between normal and novel ones"
        
      }
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      
      
      # if there is no open reading frame for normal transcript then set to first value in new_transcript being TRUE
      if (is.na(normal_transcript$AA_seq)) {
        difference_new_to_old<-TRUE
      } else {
        
        # compare longest possible AA sequences to see if there is a deviation from cannonical / annotated protein produced
        # TRUE indicates a difference
        difference_new_to_old<-!str_split(new_transcript$AA_seq, "")[[1]] == 
          str_split(normal_transcript$AA_seq, "")[[1]]
      }
      
      
      # find first difference in AA sequence between new and old - if there is a difference, add that to the result, if not state the longest AA seq would produce the same protein
      if (any(difference_new_to_old=="TRUE")) {
        # find first mismatch / TRUE value
        first_difference_new_to_old<-which(difference_new_to_old)[1]
        # take sequence from first mismatch
        seq_aa_from_first_diff<-str_sub(string = new_transcript$AA_seq, 
                                        start = first_difference_new_to_old, 
                                        end=length(difference_new_to_old))
        # add this to result
        new_transcript$seq_aa_from_first_diff<-seq_aa_from_first_diff
        
      } else {
        # if no difference - state that
        new_transcript$seq_aa_from_first_diff<-"No difference in longest transcript between normal and novel ones"
        
      }
      
      # add annotated transcript NT seq
      new_transcript$old_transcript_NT_seq<-normal_transcript$NT_seq
      # add annotated transcript AA seq
      new_transcript$old_transcript_CDS_seq<-normal_transcript$CDS_seq
      # add annotated transcript AA seq
      new_transcript$old_transcript_AA_seq<-normal_transcript$AA_seq
      
      
      new_transcript$novel_junction_type<-'novel_acceptor'
      
      new_transcript$genes<-filt_novel_a_sig_leaf_tdp$genes
      
      # add leafcutter LRRK2 analysis deltapsi of junction to result
      new_transcript$novel_junction_dpsi<-filt_novel_a_sig_leaf_tdp$deltapsi
      # add tdp kd deltapsi of junction to result
      new_transcript$novel_junction_dpsi_tdp_kd_avg<-filt_novel_a_sig_leaf_tdp$tdp_kd_average_deltapsi
      
      #add in NMD prediction
      new_transcript$NMD<-NMDresult
      
      # add in TSL support level
      filt_gtf_df_gr2<-makeGRangesFromDataFrame(rbind(filt_gtf_df, transcript_df), keep.extra.columns = T) %>% 
        subset( type == "transcript") 
      
      new_transcript$TSL<-unique(filt_gtf_df_gr2@elementMetadata$transcript_support_level)
      
      new_transcript_list_acceptors[[j]]<-new_transcript
      
    }
    
    return(new_transcript_list_acceptors)
    
    
  }
  
}, mc.cores = num_cores)

# Combine all results into one dataframe
new_transcript_list_acceptors <- dplyr::bind_rows(results_list)

# remove rows that will cause excel to break - 32,767 characters per cell in excel - also likely a multi-junction case that would require complicated knowing of which ones work together 
# so cut off at 12k since largent known exon is in dystrophin at 10kb and unlikely bigger than that
new_transcript_list_acceptors<-new_transcript_list_acceptors[!nchar(new_transcript_list_acceptors$NT_seq)>12000,]

# write results
# write_csv(x = new_transcript_list_acceptors, 
#           file = file.path(results_path, 
#                            "Upreg_tdp_tsl_nmd_included_Rio_LRRK2_novel_acceptor_junction_NT_CDS_AA_sequences_full.csv"))
# 
# 
# write_csv(x = new_transcript_list_acceptors %>% 
#             dplyr::filter(!seq_aa_from_first_diff=="No difference in longest transcript between normal and novel ones") %>%
#             arrange(desc(novel_junction_dpsi)), file = file.path(results_path, "Upreg_tdp_tsl_nmd_included_Rio_LRRK2_novel_acceptor_junction_NT_CDS_AA_sequences_filtered_no_annotation_AA_matches.csv"))

write_csv(x = new_transcript_list_acceptors %>% 
            dplyr::filter(!seq_aa_from_first_diff=="No difference in longest transcript between normal and novel ones") %>%
            select(Transcript, AA_seq, seq_aa_from_first_diff, novel_junction_type,genes, novel_junction_dpsi) 
          %>% arrange(desc(novel_junction_dpsi)),
          file = file.path(results_path, "Novel_acceptor_junction_ORF_AA_sequences.csv"))




