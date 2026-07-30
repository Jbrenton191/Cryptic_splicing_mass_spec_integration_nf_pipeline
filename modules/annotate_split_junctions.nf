process annotate_junctions {
  
  publishDir "${projectDir}/output/junction_annotation", mode: 'copy', overwrite: true
  
  input:
    val(Rpackages_done)
    path(gtf_path)
    path(junction_list_path)
    
  output:
    path("*donors.csv"), emit: junction_list_donors
    path("*acceptors.csv"), emit: junction_list_acceptors
    
  script:
    """
    Rscript ${projectDir}/Rscripts/Junction_annotation.R $gtf_path $junction_list_path
    """
}