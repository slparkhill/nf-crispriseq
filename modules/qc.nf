// Do quality control checks
process fastQC {

   tag "${id}"
   label 'big_mem'
   
   input:
   tuple val( id ), path ( reads )

   output:
   tuple val( id ), path ( "*.zip" ), emit: main
   path "*.zip", emit: logs

   script:
   """
   zcat ${reads} > "${id}.fastq"
   fastqc --noextract --memory 10000 --threads ${task.cpus} "${id}.fastq"
   rm "${id}.fastq"
   """
   stub:
   """
   zcat ${reads} | head -n1000 > "${id}.fastq"
   fastqc --noextract --memory 10000 --threads ${task.cpus} "${id}.fastq"
   rm "${id}.fastq"
   """

}