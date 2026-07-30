process genome_download {
  
  
  myDir = file("${params.output}/reference_downloads")
  myDir.mkdir()
  
  publishDir "${params.output}/reference_downloads", mode: 'copy', overwrite: true
  
  output:
  path("*primary_assembly*.fa"), emit: fasta
  path("*.gtf"), emit: gtf

  script:
  """
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v[0-9]?.annotation.gtf.gz
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.primary_assembly.genome.fa.gz

	gunzip GRCh38.primary_assembly.genome.fa.gz
	gunzip gencode.*.annotation.gtf.gz
	"""
}