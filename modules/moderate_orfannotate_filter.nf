process moderate_orfannotate_filter {
  
  publishDir "${projectDir}/output/ORFannotate_filtered", mode: 'copy', overwrite: true
  
  input:
  path(summaries)
  path(protein_fastas)
  val(novel_name_file)
  val(annotated_name_file)
  
  output:
  path("*novel_acceptor_transcripts_gtf.rds"), emit: novel_acceptor_gtf
  path("*annotated_acceptor_transcripts_gtf.rds"), emit: annotated_acceptor_gtf
  
  path("*novel_transcript_summary_out.csv"), emit: filtered_novel_summary
  path("*annotate_transcript_summary_out.csv"), emit: filtered_annotate_summary
  
  path("*novel_NMD_predicted_transcripts.csv"), emit: NMD_predictions
  
  path("*novel_prot_out.fa"), emit: filtered_novel_proteins
  path("*annotated_prot_out.fa"), emit: filtered_annotated_proteins
  
  script:
    """
    annotated_folder_name=`echo ${annotated_name_file} | sed 's/\\.txt//'`
    annotated_folder=${projectDir}/output/ORFannotate/\${annotated_folder_name}
    
    novel_folder_name=`echo ${novel_name_file} | sed 's/\\.txt//'`
    novel_folder=${projectDir}/output/ORFannotate/\${novel_folder_name}
    
    results_name=`echo "\${novel_folder_name}" | sed -E 's/_novel//g'`

    Rscript ${projectDir}/R_scripts/moderate_filter_orfannotate_output.R \${novel_folder}/ORFannotate_summary.tsv \${annotated_folder}/ORFannotate_summary.tsv \${novel_folder}/protein.fa \${annotated_folder}/protein.fa \$results_name
    """
}
