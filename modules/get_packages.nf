process get_packages {
  
  myDir = file("${params.output}")
  myDir.mkdirs()
  
  publishDir "${myDir}", mode: 'copy', overwrite: true
  
  output:
    val(fin_val), emit: pack_done_val
  
  script:
  fin_val="Package download didn't fail"
  """
	Rscript ${projectDir}/R_scripts/Rpackage_download.R
	"""

}