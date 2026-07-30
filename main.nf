nextflow.enable.dsl=2

     params.data_mass_spec=""  
     params.data_junctions="test.csv" 
    // params.metadata_csv= "${projectDir}/ASAP_samples_master_spreadsheet_25.8.21.csv"
  //   params.metadata_key= "${projectDir}/key_for_metadata.txt"
     params.output= "${projectDir}/output"
     
     params.filter_mode = 'moderate' // default value


include { get_packages } from './modules/get_packages'
     
include { genome_download } from './modules/gencode_genome_download'
     
include { annotate_junctions } from './modules/annotate_split_junctions'
include { junction_chunk_split as chunk_split_donors } from './modules/junction_chunk_split'
include { junction_chunk_split as chunk_split_acceptors } from './modules/junction_chunk_split'

include { novel_donor_transcript_gtf } from './modules/novel_donor_transcript_gtf'
include { novel_acceptor_transcript_gtf } from './modules/novel_acceptor_transcript_gtf'

include { ORFannotate as ORFannotate_donors } from './modules/ORFannotate'
include { ORFannotate as ORFannotate_acceptors } from './modules/ORFannotate'

include { moderate_orfannotate_filter as moderate_orfannotate_filter_donor } from './modules/moderate_orfannotate_filter'
include { moderate_orfannotate_filter as moderate_orfannotate_filter_acceptor } from './modules/moderate_orfannotate_filter'
include { strong_orfannotate_filter as strong_orfannotate_filter_donor } from './modules/strong_orfannotate_filter'
include { strong_orfannotate_filter as strong_orfannotate_filter_acceptor } from './modules/strong_orfannotate_filter'


// Run with:  nextflow run pipeline.nf --filter_mode strong

workflow {
  
// rpackage_download!!!
// gtf and fasta download!!

// Get input data and view:
junction_files=Channel.fromPath("$params.data_junctions")
junction_files.view()

get_packages()

genome_download()

annotate_junctions(get_packages.out.pack_done_val, genome_download.out.gtf, junction_files)

chunk_split_donors(annotate_junctions.out.junction_list_donors)
chunk_split_acceptors(annotate_junctions.out.junction_list_acceptors)

novel_donor_transcript_gtf(genome_download.out.gtf, genome_download.out.fasta, chunk_split_donors.out.novel_donor_chunk.flatten())
novel_acceptor_transcript_gtf(genome_download.out.gtf, genome_download.out.fasta, chunk_split_acceptors.out.novel_acceptor_chunk.flatten())

ORFannotate_donors(genome_download.out.fasta, novel_donor_transcript_gtf.out.novel_donor_gtf, novel_donor_transcript_gtf.out.annotated_donor_gtf)

ORFannotate_acceptors(genome_download.out.fasta, novel_acceptor_transcript_gtf.out.novel_acceptor_gtf, novel_acceptor_transcript_gtf.out.annotated_acceptor_gtf)

if (params.filter_mode == 'moderate') {
   moderate_orfannotate_filter_donor(ORFannotate_donors.out.summary,ORFannotate_donors.out.protein_fasta,ORFannotate_donors.out.novel_folder_name_text, ORFannotate_donors.out.annotated_folder_name_text)
   moderate_orfannotate_filter_acceptor(ORFannotate_acceptors.out.summary,ORFannotate_acceptors.out.protein_fasta,ORFannotate_acceptors.out.novel_folder_name_text, ORFannotate_acceptors.out.annotated_folder_name_text)
    
} else if (params.filter_mode == 'strong') {
    strong_orfannotate_filter_donor(ORFannotate_donors.out.summary,ORFannotate_donors.out.protein_fasta,ORFannotate_donors.out.novel_folder_name_text, ORFannotate_donors.out.annotated_folder_name_text)
    strong_orfannotate_filter_acceptor(ORFannotate_acceptors.out.summary,ORFannotate_acceptors.out.protein_fasta,ORFannotate_acceptors.out.novel_folder_name_text, ORFannotate_acceptors.out.annotated_folder_name_text)
    
} else {
    error "Invalid filter_mode: '${params.filter_mode}'. Must be 'moderate' or 'strong'."
} // filtering end




} // pipeline end
