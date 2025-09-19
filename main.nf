#!/usr/bin/env nextflow

/*
========================================================================================
   CRISPRi-seq Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-crispriseq
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
pipeline_title = """\
   S C B I R   C R I S P R i   P O O L E D   F I T N E S S   P I P E L I N E
   =========================================================================
   Nextflow pipeline to count guides from SRA files and calculate fitness changes.
   """
   .stripIndent()

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println pipeline_title + """\

         Usage:
            nextflow run sbcirlab/nf-crispriseq --help
            nextflow run sbcirlab/nf-crispriseq --sample_sheet <csv> --fastq-dir <dir>
            nextflow run sbcirlab/nf-crispriseq -c <config-file>

         Required parameters:
            sample_sheet      Path to SRA run table identifying the SRA IDs to download, and with columns 
                              corresponding to conditions and replicates.
            ---
            conditions
            sequencing_group
            expansion_group
            reference


         Optional parameters (with defaults):
            name_column = ${params.name_column}         Which column from the guide table to use as the guide name.
            sequence_column = ${params.sequence_column}         Which column from the guide table to use as the guide sequence.
            trim_qual = ${params.trim_qual}             For `cutadapt`, the minimum Phred score for trimming 3' calls
            min_length = ${params.min_length}           For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded.

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   System.exit(0)
}

/*
========================================================================================
   Check parameters
========================================================================================
*/

if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}
if ( !params.from_sra ) {
   if ( !params.fastq_dir ) {
      throw new Exception("!!! PARAMETER MISSING: Please provide a path to fastq_dir")
   }
}

log.info pipeline_title + """\
   inputs
      input dir.       : ${params.inputs}
      FASTQ dir.       : ${params.fastq_dir}
      sample sheet     : ${params.sample_sheet}
      guides provided  : ${params.guides}
   UMI mode            : ${params.use_umis}
   Clone BC mode       : ${params.use_clone_bc}
   SRA options
      SRA mode         : ${params.from_sra}
   table to FASTA
      name columns     : ${params.name_column}
      seq columns      : ${params.sequence_column}
   trimming 
      quality          : ${params.trim_qual}
      minimum length   : ${params.min_length}
   output              : ${params.outputs}
   """
   .stripIndent()

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

include { 
   gff2table;
   table2fasta;
} from './modules/bioino.nf'
include { 
   design_guides_with_crispio;
   map_guides_to_genome_features;
} from './modules/crispio.nf'
include { 
   reads_per_guide;
   reads_per_umi;
} from './modules/count_utils.nf'
include { 
   anchor_sequences; 
   count_guides_with_cutadapt;
} from './modules/demux.nf'
include { 
   download_eggnog_databases;
   get_functional_sets_with_eggnog;
} from './modules/eggnog.nf'
include { 
   calculate_relative_fitness;
   stack_tables;
} from './modules/fitness.nf'
include { multiQC } from './modules/multiqc.nf'
include { 
   fetch_genome_from_NCBI; 
   prefetch_from_SRA;
   download_FASTQ_from_SRA;
} from './modules/ncbi.nf'
include { 
   plot_count_distributions;
} from './modules/plots.nf'
include { fastQC } from './modules/qc.nf'
include { 
   reverse_complement;
} from './modules/seqtk.nf'
include { 
   trim_using_cutadapt; 
} from './modules/trimming.nf'
include { 
   fastq2tab;
   count_tab;
   join_umi_and_read_counts;
   UMItools_count_tab;
   UMItools_extract;
} from './modules/umitools.nf'

workflow {

   Channel.fromPath( 
         params.sample_sheet, 
         checkIfExists: true 
      )
      .splitCsv( header: true )
      .set { csv_ch }

   csv_ch
      .map { tuple(
         it.sample_id, //[expt: it.expt_id, sample: it.sample_id],
         tuple( it.adapter_read1_5prime, it.adapter_read2_5prime ),
         tuple( it.adapter_read1_3prime, it.adapter_read2_3prime ),
      ) }
      .unique()
      .set { adapter_ch }  // sample_name, [adapt5], [adapt3]

   csv_ch
      .map { tuple(
         it.sample_id,
         it.expt_id, 
      ) }
      .unique()
      .set { sample2expt }  // sample_name, expt_id

   csv_ch
      .map { tuple( 
         it.sample_id, //[expt: it.expt_id, sample: it.sample_id],
         tuple( it.umi_read1, it.umi_read2 ),
      ) }
      .unique()
      .set { umi_ch }  // sample_id, [umis]

   csv_ch
      .map { tuple( 
         it.sample_id, //[expt: it.expt_id, sample: it.sample_id],
         it.genome,
         it.pam,
         it.scaffold,
      ) }
      .unique()
      .set { genome_pam_ch }

   if ( params.from_sra ) {

      csv_ch
         .map { tuple( 
            it.sample_id, //[expt: it.expt_id, sample: it.sample_id], 
            it.Run,
         ) }
         | prefetch_from_SRA
         | download_FASTQ_from_SRA
      download_FASTQ_from_SRA.out
         .transpose()
         .groupTuple( by: 0 )
         .set { reads_ch }  // sample_id, [reads]

   } 
   
   else {

      csv_ch
         .map { tuple( 
            it.sample_id, //[expt: it.expt_id, sample: it.sample_id],
            file( 
               "${params.fastq_dir}/*${it.reads}*",
               checkIfExists: true
            ).sort()
         ) }
         .unique()
         .transpose()
         .groupTuple( by: 0 )
         .set { reads_ch }  // sample_id, [reads]

   }

   /*
   ========================================================================================
      Processing
   ========================================================================================
   */


   reads_ch | fastQC

   fetch_genome_from_NCBI(
      genome_pam_ch
         .map { it[1] }  // genome_acc
         .unique(),
      Channel.value( true ), // include protein FASTA for eggNOG
   )

   Channel.of( params.eggnog_url ) 
      | download_eggnog_databases
   fetch_genome_from_NCBI.out
      .map { tuple( it[0], it[1][1], it[2] ) }
      .combine( download_eggnog_databases.out ) 
      | get_functional_sets_with_eggnog

   fetch_genome_from_NCBI.out
      .map { tuple( it[0], it[1][0] ) }
      .combine( get_functional_sets_with_eggnog.out.gff, by: 0 )
      .set { genome_info }
   
   trim_using_cutadapt(
      reads_ch.combine( adapter_ch, by: 0 ),  // sample_id, [reads], [adapt5], [adapt3]
      Channel.value( params.trim_qual ), 
      Channel.value( params.min_length ),
      Channel.value( params.retain_5prime ),
   )  // sample_id, [reads]
   trim_using_cutadapt.out.main.set { trimmed }
   trim_using_cutadapt.out.logs.set { trim_logs }
   
   if ( params.guides ) {

      csv_ch
         .map { tuple( 
            it.sample_id, 
            it.guides_filename,
            file( 
               "${params.inputs}/${it.guides_filename}", 
               checkIfExists: true,
            ),
         ) }
         .unique()
         .set { guide_csv }  // sample_id, guide_filename, guide_file

      table2fasta(
         guide_csv
            .map { it[1..2] }  // guide_filename, guide_file
            .unique(),
         Channel.value( params.sequence_column ),
         Channel.value( params.name_column ),
      )

      guide_csv
         .map { it[1..0] }  // guide_filename, sample_id
         .unique()
         .combine( table2fasta.out, by: 0 )  // guide_filename, sample_id, guide_fasta
         .map { it[1..-1] }  // sample_id, guide_fasta
         .unique()     
         .set { guide_fasta0 }

      guide_fasta0
         .combine( genome_pam_ch, by: 0 )  // sample_id, guide_fasta, genome_acc, pam, scaffold
         .map { it[2..4] + [ it[1] ] }  // genome_acc, pam, scaffold, guide_fasta
         .unique()
         .combine( 
            genome_info, 
            by: 0,
         )  // genome_acc, pam, scaffold, guide_fasta, genome_fasta, genome_gff
         | map_guides_to_genome_features
      map_guides_to_genome_features.out.main
         .set { guide_gff }  // genome_acc, pam, guide_gff

   } else {

      genome_pam_ch
         .map { it[1..-1] }  // genome_acc, pam, scaffold
         .unique()
         .combine( 
            genome_info, 
            by: 0,
         )  // genome_acc, pam, genome_fasta, genome_gff
         .combine( Channel.of( params.guide_length ) )  // genome_acc, pam, genome_fasta, genome_gff, guide_length
         | design_guides_with_crispio
      design_guides_with_crispio.out.main
         .set { guide_gff }  // genome_acc, guide_gff

   }

   guide_gff
      .map { tuple( [id: it[0], pam: it[1]], it[2] ) }
      | gff2table
      | map { tuple( it[0].id, it[0].pam, it[1] ) }
      | set { gff_table }

   if ( ! params.guides ) {

      gff_table
         .map { tuple( 
            it[-1].simpleName, 
            it[-1], 
            "guide_sequence", 
            params.guide_name,
         ) } 
         | table2fasta  // guide_id, guide_fasta

      gff_table
         .map { [ it[-1].simpleName ] + it[0..1] }  // guide_id, genome_acc, pam
         .combine( table2fasta.out, by: 0 )  // guide_id, genome_acc, pam, guide_fasta
         .map { it[1..-1] }  // genome_acc, pam, guide_fasta
         .combine( 
            genome_pam_ch
               .map { it[1..2] + [ it[0] ] },  // genome_acc, pam, sample_id
            by: [0, 1],
         )  // genome_acc, pam, guide_fasta, sample_id
         .map { it[-1..-2] }  // sample_id, guide_fasta
         .unique()
         .set { guide_fasta0 }
   }

   if ( params.rc ) {  // reverse complement

      guide_fasta0
         .map { tuple( it[1], it[1] ) }
         .unique()
         | reverse_complement
      guide_fasta0
         .map { it[-1..0] }
         .combine( reverse_complement.out, by: 0 )
         .map { it[1..-1] }
         .unique()
         | set { guide_fasta }

   }
   else {  // pass through

      guide_fasta0.set { guide_fasta }
      
   }

   if ( params.use_umis ) {

      if ( params.use_clone_bc ) {

         UMItools_whitelist(
            trimmed.combine( umi_ch, by: 0 ),
            Channel.value( params.umi_error_cutoff ),
         )
            | set { whitelist }
         
      } else {

         trimmed
            .map { it[0] }
            .unique()
            .combine( Channel.of( file( "placeholder" ) ) )
            .set { whitelist }

      }

      UMItools_extract(
         trimmed
            .combine( umi_ch, by: 0 )
            .combine( whitelist, by: 0 ),
         Channel.value( params.trim_qual ),
         Channel.value( params.use_clone_bc ),
      )
      UMItools_extract.out.main
         .set { pre_demux }

   } else {

      trimmed.set { pre_demux }

   }

   count_guides_with_cutadapt(
      pre_demux  // sample_id, reads
         .combine( guide_fasta, by: 0 ).unique(),   // sample_id, reads, guide_fasta
      Channel.value( params.allow_guide_errors ),
   )

   count_guides_with_cutadapt.out.main 
      | fastq2tab  // sample_id, tab
      | count_tab
   
   if ( params.use_umis ) {

      UMItools_count_tab(
         fastq2tab.out,
         Channel.value( params.umi_method ),
         Channel.value( params.use_clone_bc ),
      )
      UMItools_count_tab.out.main
         .combine( count_tab.out, by: 0 )
         | join_umi_and_read_counts
         | set { guide_counts }

   } else {

      count_tab.out
         .set { guide_counts }
         
   }

   guide_counts | plot_count_distributions
   ANNOTATE_COUNTS_WITH_GENOME_FEATURES(
      gff_table  // genome_acc, pam, guide_tsv
         .map { tuple( it[0], it[2] ) }  // genome_acc, guide_tsv
         .unique()
         .combine( 
            genome_pam_ch
               .map { it[1..0] }
               .unique(), 
            by: 0,
         )  // genome_acc, guide_tsv, sample_id
         .map { it[2..1] }  // sample_id, guide_tsv
         .unique()
         .combine( guide_counts, by: 0 ),
      Channel.value( params.guide_name ),
   )
   

   if ( params.do_fitness ) {
      
      guide_counts
         .combine( sample2expt, by: 0 )
         .map { tuple( it[-1], it[1] ) }
         .groupTuple( by: 0 )
         | stack_tables
      calculate_relative_fitness(
         stack_tables.out,
         Channel.value( file( "${params.sample_sheet}" ) ),
         Channel.value( params.growth ? file( "${params.growth}" ) : file( "placeholder" ) ),
         Channel.value( params.guide_name ),
         Channel.value( params.reference_guide ),
         Channel.value( params.use_umis ),
         Channel.value( !params.growth && params.use_spike ),
         Channel.value( params.timepoint_column ),
         Channel.value( params.concentration_column ),
      )
      // JOIN_GFF(calculate_relative_fitness.out.table, gff_table)
      // PLOT_FITNESS(JOIN_GFF.out, essential_ch)
   }

   trim_logs
      .concat(
         fastQC.out.logs,
         // count_guides_with_cutadapt.out.logs,
      )
      .flatten()
      .unique()
      .collect()
      | multiQC

}


// stack count TSV files and merge the condition table
process STACK_JOIN_CONDITIONS {

   tag "${conditions}"

   label 'med_mem'
   time '24h'

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   path '??/*'
   path conditions

   output:
   path "counts.tsv" 

   script:
   """
   files=(??/*)
   for f in "\${files[@]}"
   do
      python ${projectDir}/bin/join.py "${conditions}" "," \
      < "\$f" \
      > "\$f.joined.tsv"
   done

   outputs=(*.joined.tsv)

   head -n1 "\${outputs[0]}" \
   | cat - <(tail -q -n+2 "\${outputs[@]}") \
   > counts.tsv

   """
}

// merge the guide table
process ANNOTATE_COUNTS_WITH_GENOME_FEATURES {
   
   tag "${id}"
   label 'med_mem'

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( guide_tsv ), path( counts )
   val guide_name

   output:
   tuple val( id ), path( "annotated.tsv" )

   script:
   """
   GUIDE_NAME="${guide_name}"
   GUIDE_NAME_COL=\$(
      head -n1 "${guide_tsv}" \
      | tr \$'\\t' \$'\\n' \
      | grep -n "\$GUIDE_NAME"  \
      | cut -d: -f 1
   )

   cut -f10-21,\$GUIDE_NAME_COL < "${guide_tsv}" \
   | sed 's/'"\$GUIDE_NAME"'/guide_name/' \
   > mini

   sed 's/'"\$GUIDE_NAME"'/guide_name/' < "${counts}" \
   | python ${projectDir}/bin/join.py mini \
   > annotated.tsv

   """
}

// merge the guide table
process JOIN_GFF {
   tag "${gfftable}"

   label 'med_mem'

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   path fit_params 
   path fitted 
   path gfftable 

   output:
   path( "*params*.tsv", includeInputs: true )
   path "*-fit-annotated.tsv"

   script:
   """
   python ${projectDir}/bin/join.py ${gfftable} \
   < fitness_params-guide_name.tsv \
   > fitness_params-guide_name-annotated.tsv
   python ${projectDir}/bin/join.py ${gfftable} \
   < ${fitted} \
   > ${fitted.baseName}-annotated.tsv

   """
}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/