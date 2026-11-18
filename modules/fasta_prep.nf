process PREPARE_FASTA {
    tag { meta.sample_id }
    label 'python_container'
    label 'process_single'

    publishDir "${params.outdir}/assemblies", mode: 'copy'

    input:
    tuple val(meta), path(source_fasta, stageAs: 'input_sequence')
    path(preparation_script, stageAs: 'prepare_fasta.py')

    output:
    tuple val(meta), path("${meta.sample_id}.fasta"), emit: assemblies

    script:
    def sampleId = meta.sample_id?.toString()
    if (!sampleId || !(sampleId ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error "Invalid sample ID '${sampleId}' for FASTA input."
    }
    """
    set -euo pipefail
    python3 prepare_fasta.py "${source_fasta}" "${sampleId}.fasta"
    """
}
