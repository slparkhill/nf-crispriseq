process anchor_sequences {

    tag "${id}" 

    label 'big_cpu'

    publishDir( 
        "${params.outputs}/references", 
        mode: 'copy',
        // saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( reads ), path( fastas )

    output:
    tuple val( id ), path( "${id}.appended.fastq.gz" ), path( "appended-dedup.fasta" )

    script:
    def seq_to_append = "TGGG"  // to anchor guides of different lengths
    def qual_to_append = "FFFF" * 2
    """
    APPEND="${seq_to_append}"
    RC_APPEND=\$(echo \$APPEND | tr ACGTacgt TGCAtgca | rev)
    PAL_APP="\$APPEND\$RC_APPEND"

    awk -v extras="\$PAL_APP" \
        '/^>/; /^[ATCGatcg]/ { print extras toupper(\$0) extras }' \
        ${fastas} \
    > "appended.fasta"

    # de-duplicate
    grep '^[ATCG]' "appended.fasta" \
    | sort \
    | uniq -c \
    | awk -F' ' '\$1 > 1 { print \$2 }' \
    > duplicate-seqs.txt

    if [ \$(cat duplicate-seqs.txt | wc -l) -gt 0 ] 
    then
        grep -Fx -B1 \
            --no-group-separator \
            -f duplicate-seqs.txt \
            "appended.fasta" \
        > duplicate-seqs2.txt

        grep -Fvx \
            --no-group-separator \
            -f duplicate-seqs2.txt \
            "appended.fasta" \
        > appended-dedup.fasta
    else
        ln -s appended.fasta appended-dedup.fasta
    fi

    zcat ${reads} \
    | awk -v extras="\$PAL_APP" -v quals="${qual_to_append}" '
        (((NR + 3) % 4 == 0) || ((NR + 3) % 4 == 2)); 
        (NR + 3) % 4 == 1 { print extras \$0 extras }; 
        (NR + 3) % 4 == 3 { print quals \$0 quals }
    ' \
    | pigz -v --threads ${task.cpus} --stdout \
    > "${id}.appended.fastq.gz"
    
   """
}


process count_guides_with_cutadapt {

    tag "${id}" 

    label 'big_cpu'

    publishDir( 
        "${params.outputs}/demultiplexed", 
        mode: 'copy',
        // saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( reads ), path( fasta )

    output:
    tuple val( id ), path( "*.matched.fastq.gz" ), emit: main
    tuple val( id ), path( "*.log" ), emit: logs
    tuple val( id ), path( "*.unmatched.fastq.gz" ), emit: unmatched

    script:
    """
    cutadapt \
        -g '^file:${fasta}' \
        -e 0 \
        -j ${task.cpus} \
        --no-indels \
        --report full \
        --action lowercase \
        --rename '{id} {adapter_name} {comment}' \
        --untrimmed-output "${id}.unmatched.fastq.gz" \
        -o "${id}.matched.fastq.gz" \
        ${reads} \
    > "${id}.matched.cutadapt.log"

   """
}
