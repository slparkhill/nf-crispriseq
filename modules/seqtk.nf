// Reverse complement the FASTA file.
process reverse_complement {

   tag "${id}"

   publishDir( 
        "${params.outputs}/references", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

   input:
   tuple val( id ), path( fasta )

   output:
   tuple val( id ), path( "rc.fasta" )

   script:
   """
   seqtk seq -r "${fasta}" > rc.fasta
   
   """
}