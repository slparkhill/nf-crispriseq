process stack_tables {

   tag "${id}"

   publishDir(
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( 'counts-??/*' )

   output:
   tuple val( id ), path( "all-expt-counts.tsv" )

   script:
   """
   #!/usr/bin/env python

   from glob import glob

   import pandas as pd

   (
      pd.concat(
         [
            pd.read_csv(f, sep="\\t")
            for f in glob("counts-??/*.tsv")
         ], 
         axis=0,
      )
      .to_csv(
         "all-expt-counts.tsv", 
         sep="\\t", 
         index=False,
      )
   )
   
   """

}

process calculate_relative_fitness {

   tag "${id}"
   label "big_time"

   publishDir(
      "${params.outputs}/fitness", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( counts )
   path sample_sheet
   path growth
   val guide_name
   val reference_guide
   val use_umis
   val use_spike
   val timepoint_column
   val concentration_column

   output:
   tuple val( id ), path( "fitness.tsv" ), emit: table
   tuple val( id ), path( "predictions.tsv" ), emit: predictions
   tuple val( id ), path( "regression-input.tsv" ), emit: reg_inputs
   tuple val( id ), path( "fitness/" ), emit: plots
   path "*.log", emit: logs

   script:
   """
   python "${projectDir}"/bin/fitness.py fitness \
      "${counts}" \
      --count-column ${use_umis ? "umi_count" : "read_count"} \
      --guide-column "guide_name" \
      --conditions "${sample_sheet}" \
      --reference "${reference_guide}" \
      --timepoint-column "${timepoint_column}" ${concentration_column ? "--concentration-column ${concentration_column}" : ""} \
      ${use_spike ? "--spike spike" : "--growth ${growth}"} \
      --sample-column "sample_id" \
      --format TSV \
      --plot fitness \
      --output fitness.tsv \
   2> fitness.log

   """
}


// Use `crispin` to plot fitness 
// process PLOT_FITNESS {
//    tag{"${fitted}"}

//    label 'big_mem'

//    publishDir( model_o, 
//                mode: 'copy' )

//    input:
//    path fit_params 
//    path fitted 
//    path essentials 

//    output:
//    path "*.png"

//    script:
//    """
//    guideplot \
//       --fitness fitness_params-guide_name-annotated.tsv \
//       --expansion fitness_params-exp_group.tsv \
//       --essentials ${essentials} \
//       --essential_calls ${params.essential_call} \
//       --essential_scores ${params.essential_score} \
//       --fitted ${fitted} \
//       --reference "${params.reference}" \
//       --initial "${params.initial}" \
//       --control_column "${params.control_column}" \
//       --negative ${params.negative} \
//       --count guide_count \
//       --format TSV \
//       -o fitness
//    """
}