// Convert a table of barcodes to a FASTA file for mapping.
process table2fasta {

   tag "${id}"

   publishDir( 
        "${params.outputs}/references", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

   input:
   tuple val( id ), path( table )
   val sequence_column
   val name_column

   output:
   tuple val( id ), path( "to-map.fasta" )

   script:
   if ( table.getExtension() == "csv" )
      """
      bioino table2fasta "${table}" \
         --sequence "${sequence_column}" \
         --format CSV \
         --name "${name_column}" \
         --output to-map.fasta
      """
   else
      """
      cp "${table}" to-map.fasta
      """

   stub:
   if ( table.getExtension() == "csv" )
      """
      head -n2000 "${table}" > in.csv
      bioino table2fasta in.csv \
         --sequence "${sequence_column}" \
         --format CSV \
         --name "${name_column}" \
         --output to-map.fasta
      """
   else
      """
      head -n2000 "${table}" > to-map.fasta
      """
}


// Convert a GFF to a TSV table
process gff2table {
   
   tag "${id}"

   publishDir( 
        "${params.outputs}/references", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

   input:
   tuple val( id ), path( gff )

   output:
   tuple val( id ), path( "annotations.tsv" )

   script:
   """
   bioino gff2table "${gff}" > annotations.tsv
   """
}

