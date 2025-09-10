#!/usr/bin/env nextflow

/*
========================================================================================
   INSPECT Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-crispriseq
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
pipeline_title = """\
   S C B I R   C R I S P R i   P O O L E D   F I T N E S S   P I P E L I N E
   =========================================================================
   Nextflow pipeline to count guides from SRA files and calculate
         fitness changes.
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
   anchor_sequences; 
   count_guides_with_cutadapt;
} from './modules/demux.nf'
include { multiQC } from './modules/multiqc.nf'
include { 
   fetch_genome_from_NCBI; 
   prefetch_from_SRA;
   download_FASTQ_from_SRA;
} from './modules/ncbi.nf'
include { 
   plot_UMI_distributions;
} from './modules/plots.nf'
include { fastQC } from './modules/qc.nf'
include { 
   reverse_complement;
} from './modules/seqtk.nf'
include { 
   trim_using_cutadapt; 
} from './modules/trimming.nf'
include { 
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
         it.sample_id,
         tuple( it.adapter_read1_5prime, it.adapter_read2_5prime ),
         tuple( it.adapter_read1_3prime, it.adapter_read2_3prime ),
      ) }
      .unique()
      .set { adapter_ch }  // sample_name, [adapt5], [adapt3]

   csv_ch
      .map { tuple( 
         it.sample_id,
         tuple( it.umi_pattern )
      ) }
      .unique()
      .set { umi_ch }  // sample_id, [umis]

   csv_ch
      .map { tuple( 
         it.sample_id,
         it.genome,
         it.pam,
         it.scaffold,
      ) }
      .unique()
      .set { genome_pam_ch }

   if ( params.from_sra ) {
      csv_ch
         .map { tuple( it.sample_id, it.Run ) }
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
            it.sample_id,
            file( 
               "${params.fastq_dir}/*${it.fastq_pattern}*",
               checkIfExists: true
            ).sort()
         ) }
         .set { reads_ch }  // sample_id, [reads]
   }

   /*
   ========================================================================================
      Processing
   ========================================================================================
   */


   reads_ch | fastQC

   genome_pam_ch
      .map { it[1] }  // genome_acc
      .unique()
      | fetch_genome_from_NCBI   // genome_acc, genome, gff
   
   trim_using_cutadapt(
      reads_ch.combine( adapter_ch, by: 0 ),  // sample_id, [reads], [adapt5], [adapt3]
      Channel.value( params.trim_qual ), 
      Channel.value( params.min_length ),
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
         .combine( table2fasta.out, by: 0 )  // guide_filename, sample_id, guide_fasta
         .map { it[1..-1] }  // sample_id, guide_fasta      
         .set { guide_fasta0 }

      guide_fasta0
         .combine( genome_pam_ch, by: 0 )  // sample_id, guide_fasta, genome_acc, pam, scaffold
         .map { it[2..4] + [ it[1] ] }  // genome_acc, pam, scaffold, guide_fasta
         .unique()
         .combine( fetch_genome_from_NCBI.out, by: 0)  // genome_acc, pam, scaffold, guide_fasta, genome_fasta, genome_gff
         | map_guides_to_genome_features
      map_guides_to_genome_features.out.main
         .set { guide_gff }  // genome_acc, pam, guide_gff

   } else {

      genome_pam_ch
         .map { it[1..-1] }  // genome_acc, pam, scaffold
         .unique()
         .combine( fetch_genome_from_NCBI.out, by: 0 )  // genome_acc, pam, genome_fasta, genome_gff
         .combine( Channel.of( params.guide_length ) )  // genome_acc, pam, genome_fasta, genome_gff, guide_length
         | design_guides_with_crispio
      design_guides_with_crispio.out.main
         .set { guide_gff }  // genome_acc, guide_gff

   }

   gff2table(
      guide_gff
         .map { tuple( [id: it[0], pam: it[1]], it[2] ) },
      // Channel.value( "guide_sequence" ),
      // Channel.value( params.guide_name ),
   )
   gff2table.out
      .map { tuple( it[0].id, it[0].pam, it[1] ) }
      .set { gff_table }

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
            by: [0, 1] 
         )  // genome_acc, pam, guide_fasta, sample_id
         .map { it[-1..-2] }  // sample_id, guide_fasta
         .set { guide_fasta0 }

   }

   if ( params.rc ) {  // reverse complement

      guide_fasta0
         .map { tuple( it[1].name, it[1] ) }
         .unique()
         | reverse_complement
      guide_fasta0
         .map { tuple( it[1].name, it[0] ) }
         .combine( reverse_complement.out, by:0 )
         .map { it[1..-1] }
         .set { guide_fasta }

   }
   else {  // pass through

      guide_fasta0.set { guide_fasta }
      
   }

   if ( params.use_umis ) {

      umi_ch
         .combine( 
            trimmed, 
            by: 0 
         )  // sample_id, umi_pattern, reads 
         .set { pre_umi }
      pre_umi 
         | UMITOOLS_EXTRACT  // sample_id, reads
      UMITOOLS_EXTRACT.out.main.set { pre_demux }

   } else {

      trimmed.set { pre_demux }

   }

   pre_demux  // sample_id, reads
      .combine( guide_fasta, by: 0 )   // sample_id, reads, guide_fasta
      | count_guides_with_cutadapt  // sample_id, reads
   
   if ( params.use_umis ) {

      count_guides_with_cutadapt.out.main 
         | FASTQ2TAB  // sample_id, tab
         | UMItools_count_tab  // sample_id, counts
      UMItools_count_tab.out.main
         | plot_UMI_distributions
      FASTQ2TAB.out 
         | READS_PER_UMI_AND_PER_GUIDE  // sample_id, per_umi, per_guide
      UMItools_count_tab.out.main
         .set { guide_counts }

   } else {

      COUNTS_PER_GUIDE(
         count_guides_with_cutadapt.out.main
         .combine( guide_fasta, by: 0 ),  // sample_id, reads, guide_fasta
         Channel.value( params.guide_name ),
      )
         | set { guide_counts }
         
   }

   ANNOTATE_COUNTS_WITH_GENOME_FEATURES(
      gff_table  // genome_acc, pam, guide_tsv
         .map { tuple( it[0], it[2] ) }  // genome_acc, guide_tsv
         .combine( 
            genome_pam_ch.map { it[1..0] }, 
            by: 0,
         )  // genome_acc, guide_tsv, sample_id
         .map { it[2..1] }  // sample_id, guide_tsv
         .combine( guide_counts, by: 0 ),
      Channel.value( params.guide_name ),
   )
   

   if ( params.do_fitness ) {
      csv_ch
         .map { 
            [ it.sample_id, it.ref_guide, it.ref_timepoint ] +
            it.findAll { k, v -> k.startsWith( "condition_" ) }
         }
         .set { conditions_ch }
      STACK_JOIN_CONDITIONS( 
         guide_counts
            .map { [ it[1] ] }
            .collect(),
         conditions_ch,
      )  
      | FITNESS
      JOIN_GFF(FITNESS.out, gff_table)
      PLOT_FITNESS(JOIN_GFF.out, essential_ch)
   }

   trim_logs
      .concat(
         fastQC.out.logs,
         count_guides_with_cutadapt.out.logs,
      )
      .map { it[-1] }
      .flatten()
      .unique()
      .collect()
      | multiQC

}

process FASTQ2TAB {

   tag "${id}"

   publishDir( 
        "${params.outputs}/counts", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

   input:
   tuple val( id ), path( fastqs )

   output:
   tuple val( id ), path( "tab.tsv" )

   script:
   """
   zcat ${fastqs[0]} \
      | awk '(NR + 3) % 4 == 0' \
      | tr ' ' \$'\t' \
      | cut -f1-2 \
      | sort -k2 \
      > tab.tsv

   """
}

process PLOT_READS_VS_UMIS {

   tag "${id}"

   publishDir( 
        "${params.outputs}/plot-counts", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

   input:
   tuple val( id ), path( guide_umi_counts )

   output:
   tuple val( id ), path( "umi-vs-reads.png" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd
   
   df = pd.read_csv(
      "${guide_umi_counts}", 
      sep='\\t',
   ) 

   fig, axes = scattergrid(
      df,
      grid_columns=["read_count", "umi_count"],
      log=["read_count", "umi_count"]
   )
   figsaver()(
      fig=fig,
      name='umi-vs-reads',
   )
   
   """
}


process READS_PER_UMI_AND_PER_GUIDE {

   tag "${id}"

   publishDir( 
        "${params.outputs}/counts", 
        mode: 'copy',
      //   saveAs: { "${id}.${it}" },
    )

   input:
   tuple val( id ), path( tabfile )

   output:
   tuple val( id ), path( "*.umi.tsv" ), path( "*.guide_name.tsv" )

   script:
   """
   cut -f1 ${tabfile} | cut -d _ -f2 > umi.tsv
   cut -f2 ${tabfile} > guide_name.tsv

   for f in umi.tsv guide_name.tsv
   do
      BASENAME=\$(basename \$f .tsv)
      NLINES=\$(cat \$f | wc -l)

      printf 'sample_id\\t'\$BASENAME'\\t'\$BASENAME'_read_count\\n' \
         > ${sample_id}.\$BASENAME.tsv

      sort \$f | uniq -c \
         | awk -F' ' -v OFS=\$'\\t' '{ print "${sample_id}", \$2, \$1 }' \
         | sort -k3 -n \
         >> ${sample_id}.\$BASENAME.tsv
   done

   """
}


process COUNTS_PER_GUIDE {

   tag "${id}"

   publishDir( 
        "${params.outputs}/counts", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

   input:
   tuple val( id ), path( fastqs ), path( fastas )
   val guide_name

   output:
   tuple val( id ), path( "counts.tsv" )

   script:
   """
   printf '${guide_name}\\tsample_id\\tguide_count\\n' \
      > counts.tsv

   zcat ${fastqs} \
      | awk '(NR + 3) % 4 == 0' \
      | tr ' ' \$'\\t' \
      | cut -f2 \
      | sort -k1 \
      | uniq -c \
      | awk -F' ' -v OFS='\\t' -v id="${id}" \
         '{ print \$2, id, \$1 }' \
      | sort -k1 \
      > counts0.tsv

   grep '^>' ${fastas} \
      | cut -d'>' -f2 \
      | tr -d ' ' \
      | sort -u \
      > guide-names.txt

   join -j 1 -a 1 -t\$'\\t' \
      guide-names.txt counts0.tsv  \
      | tr ' ' \$'\\t' \
      | awk -F '\\t' -v OFS='\\t' -v id="${id}" \
         'NF == 1 { print \$0, id, 0 }; NF > 1' \
      | sort -k4 -n \
      >> counts.tsv

   """
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
   path 'counts*.in.tsv'
   path conditions

   output:
   path "counts.tsv" 

   script:
   """
   for f in counts*.in.tsv
   do
      cat \$f \
      | python ${projectDir}/bin/join.py ${conditions} "," \
      > \$f.joined.tsv
   done

   head -n1 counts1.in.tsv.joined.tsv \
   | cat - <(tail -q -n+2 counts*.joined.tsv) \
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