/* Local 16S sequence identification, classic DSL2 / Nextflow 24.10.6.
 * Input: tuple(meta, nucleotide FASTA). No database downloads occur here.
 */

process VERIFY_16S_DATABASE {
    tag 'local-16s-database'
    label 'python_container'
    label 'process_single'
    cache 'deep'
    publishDir "${params.outdir}/rrna16s/database", mode: 'copy', pattern: 'database.json'

    input:
    path(db_dir, stageAs: 'db')
    path(helper, stageAs: 'rrna16s.py')

    output:
    tuple path('db'), path('database.json'), emit: verified

    script:
    """
    set -euo pipefail
    python3 rrna16s.py verify --db-dir db --output database.json
    """
}

process PREPARE_16S_FASTA {
    tag { meta.sample_id }
    label 'python_container'
    label 'process_single'
    publishDir "${params.outdir}/rrna16s/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(fasta, stageAs: 'source_sequence')
    path(helper, stageAs: 'rrna16s.py')

    output:
    tuple val(meta), path('queries.fasta'), path('query_map.tsv'), emit: queries

    script:
    if (!(meta.sample_id?.toString() ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error "Invalid 16S sample ID '${meta.sample_id}'."
    }
    """
    set -euo pipefail
    python3 rrna16s.py prepare --input source_sequence --fasta queries.fasta --map query_map.tsv
    """
}

process BLAST_16S {
    tag { meta.sample_id }
    label 'blast_container'
    label 'process_low'
    publishDir "${params.outdir}/rrna16s/${meta.sample_id}", mode: 'copy',
        pattern: '*.{blast.tsv,blast.log,versions.txt,database-info.txt}'

    input:
    tuple val(meta), path(query, stageAs: 'queries.fasta'), path(query_map, stageAs: 'query_map.tsv')
    tuple path(db_dir, stageAs: 'db'), path(db_metadata, stageAs: 'database.json')

    output:
    tuple val(meta), path("${meta.sample_id}.blast.tsv"), path('query_map.tsv'), path('database.json'), emit: hits
    tuple val(meta), path("${meta.sample_id}.versions.txt"), emit: versions
    tuple val(meta), path("${meta.sample_id}.database-info.txt"), emit: database_info
    tuple val(meta), path("${meta.sample_id}.blast.log"), emit: log

    script:
    def maxHits = params.rrna16s_max_hits as int
    """
    set -euo pipefail
    export BLASTDB="\$PWD/db"
    export BLAST_USAGE_REPORT=false
    export NCBI_DONT_USE_NCBIRC=1
    export NCBI_CONFIG__BLAST__DATA_LOADERS=none
    blastn -version > "${meta.sample_id}.versions.txt"
    blastdbcmd -db db/16S_ribosomal_RNA -info > "${meta.sample_id}.database-info.txt"
    blastn \\
        -query queries.fasta \\
        -db db/16S_ribosomal_RNA \\
        -task blastn -dust yes -evalue 1e-10 \\
        -perc_identity ${params.rrna16s_min_identity as double} \\
        -qcov_hsp_perc ${params.rrna16s_min_coverage as double} \\
        -max_hsps 1 -max_target_seqs ${maxHits} \\
        -num_threads ${task.cpus} \\
        -outfmt '6 qseqid saccver pident length qstart qend qlen sstart send slen evalue bitscore staxids sscinames stitle' \\
        -out "${meta.sample_id}.blast.tsv" \\
        2> "${meta.sample_id}.blast.log" || { cat "${meta.sample_id}.blast.log" >&2; exit 1; }
    cat "${meta.sample_id}.blast.log" >&2
    # BLAST can print a database error yet exit zero; do not accept that run.
    if grep -q '^Error:' "${meta.sample_id}.blast.log"; then exit 1; fi
    """
}

process SUMMARISE_16S {
    tag { meta.sample_id }
    label 'python_container'
    label 'process_single'
    publishDir "${params.outdir}/rrna16s/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(hits, stageAs: 'hits.tsv'), path(query_map, stageAs: 'query_map.tsv'), path(db_metadata, stageAs: 'database.json')
    path(helper, stageAs: 'rrna16s.py')

    output:
    tuple val(meta), path("${meta.sample_id}.16s.tsv"), emit: tsv
    tuple val(meta), path("${meta.sample_id}.16s.json"), emit: json

    script:
    """
    set -euo pipefail
    python3 rrna16s.py summarise \\
        --sample-id '${meta.sample_id}' --hits hits.tsv --map query_map.tsv \\
        --database database.json --output '${meta.sample_id}.16s' \\
        --min-identity ${params.rrna16s_min_identity as double} \\
        --min-coverage ${params.rrna16s_min_coverage as double} \\
        --min-length ${params.rrna16s_min_length as int} \\
        --max-ambiguous ${params.rrna16s_max_ambiguous as double} \\
        --score-window ${params.rrna16s_score_window as double} \\
        --max-hits ${params.rrna16s_max_hits as int}
    """
}
process MERGE_16S {
    tag '16s-cohort-summary'
    label 'python_container'
    label 'process_single'
    publishDir "${params.outdir}/rrna16s/summary", mode: 'copy'

    input:
    path(summaries)
    path(helper, stageAs: 'rrna16s.py')

    output:
    path('rrna16s_summary.tsv'), emit: tsv
    path('rrna16s_summary.json'), emit: json

    script:
    """
    set -euo pipefail
    python3 rrna16s.py merge --output rrna16s_summary
    """
}
