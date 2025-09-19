// Identify clone barcodes
process UMItools_whitelist {

   label 'big_mem'
   time '24h'

   tag "${sample_id}"

   publishDir( 
      "${params.outputs}/clone-barcodes", 
      mode: 'copy',
      // pattern: "*.whitelist.log", 
   )

   input:
   tuple val( id ), path( reads ), val( umis )
   val error_threshold

   output:
   tuple val( id ), path( "*.txt" ), emit: main
   tuple val( id ), path( "*.png" ), emit: plots
   path "*.log", emit: logs

   script:
   def second_reads  = ( 
      reads[1]
      ? "--bc-pattern2 '${umis[1]}' --read2-in '${reads[1]}'"
      : "" 
   )
   """
   umi_tools whitelist \
      --knee-method distance \
      --method umis \
		--bc-pattern '${umis[0]}' ${second_reads} \
      --extract-method regex \
      --error-correct-threshold ${error_threshold} \
      --ed-above-threshold correct \
      --plot-prefix "${id}.whitelist" \
		--log "${id}.whitelist.log" \
		--stdin "${reads[0]}" \
      --stdout "${id}.whitelist.txt"
   """
}


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
   val trim_qual
   val use_whitelist

   output:
   tuple val( id ), path( "*.extracted.fastq.gz" ), emit: main
   path "*.log", emit: logs

   script:
   def second_reads  = ( 
      reads[1]
      ? "--bc-pattern2 '${umis[1]}' --read2-in '${reads[1]}' --read2-out '${id}'_R2.extracted0.fastq.gz"
      : "" 
   )
   def wl_flag = ( use_whitelist ? "--whitelist ${whitelist}" : "" )
   """
   umi_tools extract \
      --extract-method regex \
      --quality-filter-mask ${trim_qual} \
      --quality-encoding phred33 \
      --bc-pattern '${umis[0]}' \
      --log "${id}.extract.log" \
      --stdin "${reads[0]}" \
      --stdout "${id}_R1.extracted0.fastq.gz" ${second_reads} ${wl_flag}

   # make sure no 0-length reads
   for f in "${id}"_R?.extracted0.fastq.gz
   do
      zcat \$f \
      | awk \
         '
         NR % 4 == 1 { name = \$0; next } 
         NR % 4 == 2 { if (length(\$0) > 0) { seq = \$0 } else { seq = "N" }; next } 
         NR % 4 == 0 { 
            if (length(\$0) > 0) { qual = \$0 } else { qual = "?" }; 
            print name ORS seq ORS "+" ORS qual 
         } 
         ' \
      | pigz -v --stdout -p ${task.cpus} \
      > \$(basename \$f .extracted0.fastq.gz).extracted.fastq.gz
   done

   #rm "${id}"_R?.extracted0.fastq.gz

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

process fastq2tab {

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
   zcat "${fastqs[0]}" \
   | awk -F' ' -v OFS='\\t' \
      '(NR + 3) % 4 == 0 { print \$1, \$2 }' \
   | sort -k2 \
   > tab.tsv

   """
}

process count_tab {

   tag "${id}"

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}-${it}" },
   )

   input:
   tuple val( id ), path( tabfile )

   output:
   tuple val( id ), path( "read_count.tsv" )

   script:
   """
   cut -f2 "${tabfile}" \
   | sort \
   | uniq -c \
   | awk -F' ' -v OFS='\\t' -v id="${id}" '
      BEGIN { print "sample_id", "guide_name", "read_count" }
      { print id, \$2, \$1 }
   ' \
   > read_count.tsv

   """
}


// Count unique UMIs per cell per gene
process UMItools_count_tab {

   tag "${id}"

   label 'big_mem'

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}-${it}" },
   )

   input:
   tuple val( id ), path( tabfile )
   val umi_method
   val per_clone

   output:
   tuple val( id ), path( "${id}.umitools_count.tsv" ), emit: main
   path "${id}.umitools_count.log", emit: logs

   script:
   """
   export MPLCONFIGDIR=tmp
   mkdir \$MPLCONFIGDIR

   umi_tools count_tab \
      --method ${umi_method} ${per_clone ? "--per-cell" : ""} \
		--stdin "${tabfile}" \
      --stdout umitools_count0.tsv \
      --log "${id}.umitools_count.log"

   awk -v OFS='\\t' -v id="${id}" '
      BEGIN { print "sample_id", "guide_name", "umi_count" }
      NR > 1 { \$1 = \$1; print id, \$0 }
   ' umitools_count0.tsv \
   > "${id}.umitools_count.tsv"

   """
}


process join_umi_and_read_counts {

   tag "${id}"

   publishDir( 
      "${params.outputs}/counts", 
      mode: 'copy',
      saveAs: { "${id}-${it}" },
   )

   input:
   tuple val( id ), path( umi_count ), path( read_count )

   output:
   tuple val( id ), path( "${id}.umi+read_count.tsv" )

   script:
   """
   #!/usr/bin/env python

   import pandas as pd

   (
      pd.read_csv("${read_count}", sep="\\t")
      .merge(
         pd.read_csv("${umi_count}", sep="\\t"),
         how="outer",
      )
      .fillna(0)
      .to_csv("${id}.umi+read_count.tsv", sep="\\t", index=False)
   )

   """
}
