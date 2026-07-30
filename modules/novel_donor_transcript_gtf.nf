process novel_donor_transcript_gtf {
  
  publishDir "${projectDir}/output/gtf_transcript", mode: 'copy', overwrite: true
  
  input:
  path(gtf_path)
  path(fasta_path)
  path(junction_list)
  
  output:
  path("*novel_donor_transcripts_gtf.rds"), emit: novel_donor_gtf
  path("*annotated_donor_transcripts_gtf.rds"), emit: annotated_donor_gtf
  
  script:
    """
    Rscript ${projectDir}/Rscripts/novel_donor_transcript_gtf.R $gtf_path $fasta_path $junction_list
    """
}