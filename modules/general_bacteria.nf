/*
 * General bacterial characterisation processes.
 * Input convention throughout this pipeline: tuple(meta, assembly_fasta).
 * All tools run against assemblies, not raw reads.
 */

process AMRFINDERPLUS {
    tag "${meta.sample_id}"
    label 'amrfinder_container'
    label 'process_medium'
    cpus 8

    publishDir "${params.outdir}/amr_reports/${meta.species_slug}", mode: 'copy'

    input:
    tuple val(meta), path(assembly, stageAs: 'assembly.fasta')

    output:
    tuple val(meta), path("${meta.sample_id}.amrfinder.tsv"), emit: results
    tuple val(meta), path("${meta.sample_id}.versions.txt"), emit: versions

    script:
    def sampleId = meta.sample_id.toString()
    if (!(sampleId ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error "Invalid sample ID '${sampleId}': use letters, numbers, underscores, dots or hyphens; start with a letter or number."
    }
    """
    set -euo pipefail

    amrfinder --database_version > "${sampleId}.versions.txt" 2>&1
    amrfinder \\
        --nucleotide assembly.fasta \\
        --threads ${task.cpus} \\
        --quiet \\
        --plus \\
        --name "${sampleId}" \\
        > "${sampleId}.amrfinder.tsv"
    """
}

process MLST {
    tag "${meta.sample_id}"
    label 'mlst_container'
    label 'process_single'

    publishDir "${params.outdir}/mlst/${meta.species_slug}", mode: 'copy'

    input:
    tuple val(meta), path(assembly, stageAs: 'assembly.fasta')

    output:
    tuple val(meta), path("${meta.sample_id}_mlst.txt"), emit: results

    script:
    """
    set -euo pipefail
    mlst assembly.fasta > "${meta.sample_id}_mlst.txt"
    """
}

process VIRULENCEFINDER {
    tag "${meta.sample_id}"
    label 'virulence_container'
    label 'process_medium'

    publishDir "${params.outdir}/virulencefinder/${meta.species_slug}", mode: 'copy'

    input:
    tuple val(meta), path(assembly, stageAs: 'assembly.fasta')
    path(virulence_converter, stageAs: 'virulence_json_to_tsv.py')

    output:
    tuple val(meta), path("${meta.sample_id}.virulence.json"), emit: json
    tuple val(meta), path("${meta.sample_id}_virulence"), emit: raw_dir
    tuple val(meta), path("${meta.sample_id}_virulence.tsv"), emit: tsv

    script:
    """
    set -euo pipefail

    outdir="${meta.sample_id}_virulence"
    mkdir -p "\$outdir"

    python3 -m virulencefinder \\
        -ifa assembly.fasta \\
        -o "\$outdir" \\
        -j "${meta.sample_id}.virulence.json"

    test -s "${meta.sample_id}.virulence.json"

    python3 virulence_json_to_tsv.py \\
        "${meta.sample_id}.virulence.json" \\
        "${meta.sample_id}_virulence.tsv"
    """
}

process PLASMIDFINDER {
    tag "${meta.sample_id}"
    label 'plasmidfinder_container'
    label 'process_medium'

    publishDir "${params.outdir}/plasmidfinder/${meta.species_slug}", mode: 'copy'

    input:
    tuple val(meta), path(assembly, stageAs: 'assembly.fasta')
    path(plasmid_converter, stageAs: 'plasmid_json_to_tsv.py')

    output:
    tuple val(meta), path("${meta.sample_id}.plasmid.json"), emit: json
    tuple val(meta), path("${meta.sample_id}_plasmidfinder"), emit: raw_dir
    tuple val(meta), path("${meta.sample_id}_plasmid.tsv"), emit: tsv

    script:
    """
    set -euo pipefail

    outdir="${meta.sample_id}_plasmidfinder"
    mkdir -p "\$outdir"

    python3 -m plasmidfinder \\
        -i assembly.fasta \\
        -o "\$outdir" \\
        -j "${meta.sample_id}.plasmid.json"

    test -s "${meta.sample_id}.plasmid.json"

    python3 plasmid_json_to_tsv.py \\
        "${meta.sample_id}.plasmid.json" \\
        "${meta.sample_id}_plasmid.tsv"
    """
}
