// Design guides from scratch.
process design_guides_with_crispio {
   
   tag "${id}:${pam}"
   label 'big_time'

   publishDir( 
      "${params.outputs}/guides", 
      mode: 'copy',
      saveAs: { "${id}.${pam}-l=${guide_length}.${it}" },
   )

   input:
   tuple val( id ), val( pam ), path( genome ), path( gff ), val( guide_length )

   output:
   tuple val( id ), val( pam ), path( "guide-design.gff" ), emit: main
   path "guide-design.log", emit: logs

   script:
   """
   crispio generate "${fasta}" \
      --genome "${genome}" \
      --annotations "${gff}" \
      --pam ${pam} \
      -o "guide-design.gff"
      2> "guide-design.log"

   """
}


// Map a FASTA of guides to a genome and annotate.
process map_guides_to_genome_features {
   
   tag "${id}:${pam}:${scaffold}"
   label 'big_time'

   publishDir( 
      "${params.outputs}/guides", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), val( pam ), val( scaffold ), path( guide_fasta ), path( genome_fasta ), path( gff )

   output:
   tuple val( id ), val( pam ), path( "mapped.gff" ), emit: main
   path "map.log", emit: logs

   script:
   """
   crispio map "${guide_fasta}" \
      --genome "${genome_fasta}" \
      --annotations "${gff}" \
      --pam "${pam}" \
   2> map.log \
   | crispio featurize \
      --scaffold "${scaffold}" \
   > mapped.gff

   """
}


// Use `crispio` to calculate fitness 
process FITNESS {

   tag "${counts}" 

   label 'big_gpu'
   time '24h'

   input:
   path counts 

   output:
   path "fitness_params-*.tsv"
   path "fitness_*-fit.tsv"

   script:
   """
   guidefitness ${counts} \
      --sequencing_group "${params.sequencing_group}" \
      --expansion_group "${params.expansion_group}" \
      --culture "${params.culture_group}" \
      --reference "${params.reference}" \
      --initial "${params.initial}" \
      --name ${guide_name} \
      --count guide_count \
      --format TSV \
      -o fitness
   """
}


// Use `crispin` to plot fitness 
process PLOT_FITNESS {
   tag{"${fitted}"}

   label 'big_mem'

   publishDir( model_o, 
               mode: 'copy' )

   input:
   path fit_params 
   path fitted 
   path essentials 

   output:
   path "*.png"

   script:
   """
   guideplot \
      --fitness fitness_params-guide_name-annotated.tsv \
      --expansion fitness_params-exp_group.tsv \
      --essentials ${essentials} \
      --essential_calls ${params.essential_call} \
      --essential_scores ${params.essential_score} \
      --fitted ${fitted} \
      --reference "${params.reference}" \
      --initial "${params.initial}" \
      --control_column "${params.control_column}" \
      --negative ${params.negative} \
      --count guide_count \
      --format TSV \
      -o fitness
   """
}
