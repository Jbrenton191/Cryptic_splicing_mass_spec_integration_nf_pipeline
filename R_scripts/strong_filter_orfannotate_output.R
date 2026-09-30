# Creating a .gtf for novel donor junctions

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Setup ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load libraries
suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(tidyverse)
  
})

conflicted::conflict_prefer_all(winner = "dplyr", quiet = T)

arguments <- optparse::parse_args(
  OptionParser(
    usage = "", 
    description="Script for extracting the junctions with ORF peptide predictions and those predicted to undergo NMD from ORFannotate output. 
  novel protein summary required in position 1. 
  annotated protein summary required in position 2. 
  novel protein fasta required in position 3. 
  annotated protein fasta required in position 4.
  name of folder with chunk name needed in position 5.  
    ",
    positional_arguments = 5)
)

orfannotate_novel_summary_in<-arguments$args[1]
orfannotate_annotated_summary_in<-arguments$args[2]

orfannotate_novel_prot_fasta_in<-arguments$args[3]
orfannotate_annotated_prot_fasta_in<-arguments$args[4]

results_folder<-arguments$args[5]

orfannotate_novel_summary<-file.path(novel_results_folder, "ORFannotate_summary.tsv")
orfannotate_annotated_summary<-file.path(annotated_results_folder, "ORFannotate_summary.tsv")

orfannotate_novel_prot_fasta<-file.path(novel_results_folder, "protein.fa")
orfannotate_annotated_prot_fasta<-file.path(annotated_results_folder, "protein.fa")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 2: Load and filter ORFannotate output  ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 2.1
novel_summary<-read_tsv(orfannotate_novel_summary)
annotated_summary<-read_tsv(orfannotate_annotated_summary)


# keep open reading frames that aren't NMD sensitive
pass<-which(novel_summary$has_orf == TRUE & 
        novel_summary$NMD_sensitive == FALSE & 
        novel_summary$has_uORF == FALSE & 
        novel_summary$kozak_strength %in% c('strong'))

# filter novel summary for those passing:
novel_out<-novel_summary[pass,]

# filter annotated summary for the matching related transcripts:
annotated_out<-annotated_summary[pass,]

# 2.2 - keep NMD predicted sequences to compare with gene level results
nmd<-which(summ$has_orf==TRUE & summ$NMD_sensitive==TRUE)
# filter for transcripts predicted to undergo NMD
novel_nmd<-novel_summary[nmd,]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 3: Load and filter ORFannotate protein sequences for the same transcripts  ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# filter the protein fastas for the same transcripts/peptides
## novel
novel_prot<-Biostrings::readAAStringSet(orfannotate_novel_prot_fasta)

novel_prot_out<-novel_prot[pass,]

## annotated
annotated_prot<-Biostrings::readAAStringSet(orfannotate_annotated_prot_fasta)

annotated_prot_out<-annotated_prot[pass,]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 4: Output all filtered summaries and protein sequences  ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Summaries #
## novel protein/ORF predictions
write_csv(novel_out, paste0(results_folder, "_novel_transcript_summary_out.csv"))
## associated annotated protein predictions that match the novel ones
write_csv(annotated_out, paste0(results_folder, "_annotated_transcript_summary_out.csv"))
## NMD sensitive transcripts
write_csv(novel_nmd, paste0(results_folder, "_novel_NMD_predicted_transcripts.csv"))

# Protein sequences #
# novel ORF protein seqs
writeXStringSet(novel_prot_out, filepath = paste0(results_folder, "_novel_prot_out.fa"))
# annotated ORF protein seqs
writeXStringSet(annotated_prot_out, filepath = paste0(results_folder, "_annotated_prot_out.fa"))
