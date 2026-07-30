process ORFannotate {
  
  conda 'orfannotate'
  
  publishDir "${projectDir}/output/ORFannotate", mode: 'copy', overwrite: true
  
  input:
  // path(gtf_path)
  path(fasta_path)
  path(novel_junction_list_gtf)
  path(annotated_junction_list_gtf)
  
  output:
  path("*ORFannotate_annotated.gtf"), emit: gtf
  path("*protein.fa"), emit: protein_fasta
  path("*ORFannotate_summary.tsv"), emit: summary
  path("*_annotated_*/**"),  emit: annotated_files
  path("*_novel_*/**"),      emit: novel_files
  val("*novel*.txt"), emit: novel_folder_name_text
  val("*annotated*.txt"), emit: annotated_folder_name_text
  
  
  script:
  
    """
    name_novel=`echo "$novel_junction_list_gtf" | sed -E 's/_transcripts_gtf.rds//g'`
    
    name_annotated=`echo "$annotated_junction_list_gtf" | sed -E 's/_transcripts_gtf.rds//g'`
    
    python ${projectDir}/orfannotate/ORFannotate.py --gtf $novel_junction_list_gtf --fa $fasta_path --outdir \$name_novel
    
    python ${projectDir}/orfannotate/ORFannotate.py --gtf $annotated_junction_list_gtf --fa $fasta_path --outdir \$name_annotated

    touch \${name_novel}/\${name_novel}.txt
    touch \${name_annotated}/\${name_annotated}.txt
    
    """
}