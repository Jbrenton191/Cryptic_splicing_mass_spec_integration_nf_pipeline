process junction_chunk_split {
  
  publishDir "${projectDir}/output/junction_chunk_split", mode: 'copy', overwrite: true
  
  input:
  path(junction_list)

  output:
  path("*novel_acceptor_chunk*.csv"), emit: novel_acceptor_chunk
  path("*novel_donor_chunk*.csv"), emit: novel_donor_chunk
  path("*.txt"), emit: error_out
  

  script:
    """
    Rscript ${projectDir}/R_scripts/Junction_chunk_split.R $junction_list
    """
}