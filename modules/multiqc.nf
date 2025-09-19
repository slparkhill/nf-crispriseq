/*
 * Make log report
 */
process multiQC {

   publishDir( 
      "${params.outputs}/multiqc", 
      mode: 'copy',
   )

   input:
   path( '*', stageAs: '??/*' )

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}