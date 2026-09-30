
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Setup ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(readxl)
  library(Biostrings)
})


conflicted::conflict_prefer_all(winner = "dplyr", quiet = T)



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


# load predicted protein sequences

## novel

### load the summary and aa sequences and combine

novel_prot<-Biostrings::readAAStringSet(orfannotate_novel_prot_fasta)

novel_summary<-read_tsv(orfannotate_novel_summary)


# load MS ions - they need a transcript name (separated by commas for multiple), 
# the amino acid sequence and counts in columns afterwards.

# Column 1 - protein names with E

ms_sequences<-read_csv()

# protein id

## aa seq

## transcript id - if multiple and comma separated change.

ms_sequences %>%
  separate_rows(Protein.Group, sep = ";") %>%
  mutate(transcript = str_remove_all(pattern = "\\_.*", string = Protein.Group)) %>%
  mutate(`Gene names`=str_remove_all(pattern = ";.*", string = `Gene names`))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 2.1 matching ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


for (i in 1:length(novel_summary$transcript_id)) {
  
  filt_ms_sequences<-ms_sequences[which(ms_sequences$transcript_id %in% novel_summary$transcript_id[i]),]
  
  for (j in 1:length(filt_ms_sequences$transcript_id)) {
    
    # novel alignment
    seq1<-filt_ms_sequences$aa_seq[j]
    seq2<-novel_prot[[i]]
    
    alignment <- pairwiseAlignment(pattern = seq1, subject = seq2,
                                    substitutionMatrix = "BLOSUM62",
                                    gapOpening = -10, gapExtension = -0.5,
                                    type = "local")
    
    annotated_AA_seq<-as.character(seq3)
    annotated_alignment_score<-alignment2@score
    annotated_full_match_within<-alignedPattern(alignment2) == alignedSubject(alignment2)
    
    annotated_number_matches<-nmatch(alignment2)     # Number of matches
    annotated_percent_matching<-(annotated_number_matches/ms_ion_seq_length)*100
    annotated_num_mismatch<-nmismatch(alignment2)  # Number of mismatches
    
    annotated_indel<-nindel(alignment2)     # Number of inserted/deleted residues (gaps)
    
    annotated_inserted_residues_number<-annotated_indel@insertion[1]
    annotated_deleted_residues_number<-annotated_indel@deletion[1]
    
    annotated_aligned_query<-as.character(alignedPattern(alignment2)) 
    
    alignment_length = nchar(aligned_query)
    # create normalised alignment score for 
    annotated_alignment_length = nchar(annotated_aligned_query)
    
    normalized_score_by_alignment_length = alignment_score / alignment_length
    normalized_annotated_score_by_alignment_length = annotated_alignment_score / annotated_alignment_length
    
    normalized_score_by_full_query_length = alignment_score / ms_ion_seq_length
    normalized_annotated_score_by_full_query_length = annotated_alignment_score / ms_ion_seq_length
    
    delta_alignment_score_by_alignment_length = normalized_score_by_alignment_length - normalized_annotated_score_by_alignment_length
    
    delta_alignment_score_by_full_query_length = normalized_score_by_full_query_length - normalized_annotated_score_by_full_query_length
    
    ratio_alignment_score_by_alignment_length = normalized_score_by_alignment_length / normalized_annotated_score_by_alignment_length
    
    ratio_alignment_score_by_full_query_length = normalized_score_by_full_query_length / normalized_annotated_score_by_full_query_length
    
    stat_results <- matching_results %>% 
      filter(
        number_matches>=7,
        number_matches<=30
      ) %>% 
      filter(delta_alignment_score_by_full_query_length>=0.5)
    
    # DO T-Test on only these proteins to see if they are differentially upregulated.
    
    # Make graphing dataframe with only Significant Upregulated novel peptides
    # If multiple possible significant novel peptides per gene - average of these is taken for Log2FC
    top_proteins <- stat_results %>%
      filter(
        t.test_fdr<0.05,
        Log2FC>0
      ) %>%
      group_by(gene, transcript_status) %>%
      mutate(avg_Log2FC=mean(Log2FC)) %>%
      arrange(desc(avg_Log2FC)) %>%
      distinct(gene,.keep_all = T) %>%
      mutate(gene = factor(gene, levels = gene)) %>%
      mutate(junction_type=case_when(junction_type=="novel_donor" ~ "Novel donor",
                                     junction_type=="novel_acceptor" ~ "Novel acceptor",
      )) %>%
      mutate(junction_type=factor(junction_type, levels=c("Novel donor", "Novel acceptor"))) 
    
    
    
    # annotated alignment
    seq1<-filt_ms_sequences$aa_seq[j]
    seq3<-annotated_prot[[i]]
    
    alignment <- pairwiseAlignment(pattern = seq1, subject = seq3,
                                   substitutionMatrix = "BLOSUM62",
                                   gapOpening = -10, gapExtension = -0.5,
                                   type = "local")
    
  }
  
  
}

ms_aa_sequences<-Biostrings::AAStringSet(x = novel_sequences_donors$Stripped.Sequence)
















