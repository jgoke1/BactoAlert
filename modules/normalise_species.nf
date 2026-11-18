process NORMALISE_SPECIES {
    tag "${meta.sample_id}"
    label 'process_single'
    label 'python_container'

    input:
    tuple val(meta), path(assembly, stageAs: 'assembly.fasta'), path(speciator_report, stageAs: 'speciator_report')
    path(normalisation_script, stageAs: 'normalise_species.py')

    output:
    tuple val(meta), path('classified_assembly.fasta'), path("${meta.sample_id}.taxonomy.json"), emit: classified

    script:
    """
    set -euo pipefail
    cp --reflink=auto assembly.fasta classified_assembly.fasta
    python3 normalise_species.py \\
        --input speciator_report \\
        --sample-id "${meta.sample_id}" \\
        --output "${meta.sample_id}.taxonomy.json"
    """
}
