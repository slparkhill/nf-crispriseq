/*
 * Extract cell barcodes and UMIs
 */
process UMItools_extract {

   tag "${id}"
   label "big_mem"
   time "2d"  // in case of very large whitelists

   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( 
      "${params.outputs}/extracted", 
      mode: 'copy',
      pattern: "*.extract.log", 
   )

   input:
   tuple val( id ), path( reads ), val( umis ), path( whitelist )

   output:
   tuple val( id ), path( "*.extracted.fastq.gz" ), emit: main
   path "*.log", emit: logs

   script:
   def error_correct = ( params.allow_cell_errors ? "--error-correct-cell" : "" )
   def second_reads  = ( 
      reads[1]
      ? "--bc-pattern '${umis[0]}' --bc-pattern2 '${umis[1]}' --read2-in '${reads[1]}' --read2-out '${id}'_R2.extracted0.fastq.gz"
      : "--bc-pattern '${umis}'" 
   )
   """
   zcat "${whitelist}" > wl.txt
   umi_tools extract \
      --extract-method regex ${error_correct} \
      --quality-encoding phred33 \
      --whitelist "wl.txt" \
      --log "${id}.${whitelist.baseName}.extract.log" \
      --stdin "${reads[0]}" \
      --stdout "${id}"_R1.extracted0.fastq.gz ${second_reads}

   # make sure no 0-length reads
   for f in "${id}"_R?.extracted0.fastq.gz
   do
      zcat \$f \
      | awk \
         '
         NR % 4 == 1 { name = \$0 } 
         NR % 4 == 2 { if (length(\$0) > 0) { seq = \$0 } else { seq = "N" } } 
         NR % 4 == 0 { if (length(\$0) > 0) { qual = \$0 } else { qual = "?" } } 
         ( NR > 1 && NR%4 == 1 ) { print name ORS seq ORS "+" ORS qual } 
         END { print name ORS seq ORS "+" ORS qual }
         ' \
      | pigz -v --stdout -p ${task.cpus} \
      > \$(basename \$f .fastq.gz).extracted.fastq.gz
   done

   rm "${id}"_R?.extracted0.fastq.gz

   """
}

process concat_extractions {

   tag "${id}"

   // errorStrategy 'retry'
   // maxRetries 2

   // publishDir( 
   //    "${params.outputs}/extracted", 
   //    mode: 'copy',
   //    pattern: "*.extract.log", 
   // )

   input:
   tuple val( id ), path( reads, stageAs: '??/*' )

   output:
   tuple val( id ), path( "*.extracted-concat.fastq.gz", arity: 1..2 )

   script:
   """
   cat "??/*_R1.extracted.fastq.gz" > "${id}"_R1.extracted-concat.fastq.gz
   if ls ??/*_R2.extracted.fastq.gz 1> /dev/null 2>&1
   then
      cat "??/*_R2.extracted.fastq.gz" > "${id}"_R2.extracted-concat.fastq.gz
   fi

   """

}

/*
 * Count unique UMIs per cell per gene
 */
process UMItools_count {

   tag "${sample_id}"

   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}-${it}" }, 
   )

   input:
   tuple val( sample_id ), path( bamfile )
   val paired

   output:
   tuple val( sample_id ), path( "*.tsv" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   samtools sort ${bamfile} -o ${bamfile.getBaseName()}.sorted.bam
   samtools index ${bamfile.getBaseName()}.sorted.bam
   umi_tools count \
		--per-gene ${paired ? '--paired --chimeric-pairs discard --unpaired-reads discard' : ''} \
      --per-cell \
		--gene-tag XT \
      --log ${sample_id}.count.log \
		--stdin ${bamfile.getBaseName()}.sorted.bam \
      --stdout ${sample_id}.umitools_count.tsv \
   && rm ${bamfile.getBaseName()}.sorted.bam

   """
}

// Count unique UMIs per cell per gene
process UMItools_count_tab {

   tag "${id}"

   label 'big_mem'
   time '48h'

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}-${it}" },
   )

   input:
   tuple val( id ), path( tabfile )

   output:
   tuple val( id ), path( "${id}.umitools_count.tsv" ), emit: main
   tuple val( id ), path( "${id}.umitools_count.log" ), emit: logs

   script:
   """
   umi_tools count_tab \
      --method unique \
		--stdin "${tabfile}" \
      --stdout umitools_count0.tsv \
      --log "${id}.umitools_count.log"

   awk -v OFS='\\t' -v id="${id}" \
      BEGIN { print "guide_name", "umi_count", "sample_id" }
      NR > 1 { \$1 = \$1; print \$0, id }
   ' umitools_count0.tsv \
   | sort -k1 \
   > umitools_count-a.tsv

   cut -f2 "${tabfile}" > read_count0.tsv

   printf 'guide_name\\tread_count\\n' \
      > read_count.tsv
   sort read_count0.tsv \
      | uniq -c \
      | awk -F' ' -v OFS='\\t' \
         '{ print \$2, \$1 }' \
      | sort -k1 \
      >> read_count.tsv

   join --header umitools_count-a.tsv read_count.tsv \
      | awk -F' ' -v OFS='\\t' \
         '{ print \$3, \$1, \$2, \$4 }' \
      > "${id}.umitools_count.tsv"

   rm read_count0.tsv

   """
}
