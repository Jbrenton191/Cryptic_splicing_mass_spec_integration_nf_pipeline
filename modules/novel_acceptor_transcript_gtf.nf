process novel_acceptor_transcript_gtf {
  
  publishDir "${projectDir}/output/gtf_transcript", mode: 'copy', overwrite: true
  
  input:
  path(gtf_path)
  path(fasta_path)
  path(junction_list)
  
  output:
  path("*novel_acceptor_transcripts_gtf.rds"), emit: novel_acceptor_gtf
  path("*annotated_acceptor_transcripts_gtf.rds"), emit: annotated_acceptor_gtf
  
  script:
    """
    Rscript ${projectDir}/Rscripts/novel_acceptor_transcript_gtf.R $gtf_path $fasta_path $junction_list
    """
}