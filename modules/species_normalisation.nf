/* Publish each per-sample normalisation result into a species-specific folder. */
process GROUP_SPECIES_NORMALISATION {
    tag "${meta.species_slug}/${meta.sample_id}"
    label 'process_single'
    label 'python_container'

    publishDir "${params.outdir}/species_groups/${meta.species_slug}", mode: 'copy'

    input:
    tuple val(meta), path(taxonomy_file, stageAs: 'input_taxonomy.json')

    output:
    tuple val(meta), path("${meta.sample_id}.taxonomy.json"), emit: grouped

    script:
    """
    set -euo pipefail
    cp "${taxonomy_file}" "${meta.sample_id}.taxonomy.json"
    """
}

/* Combine every sample's JSON into one cohort TSV and one cohort JSON file. */
process MERGE_SPECIES_NORMALISATION {
    tag 'cohort-summary'
    label 'process_single'
    label 'python_container'

    publishDir "${params.outdir}/species_normalisation/summary", mode: 'copy'

    input:
    path(taxonomy_files)
    path(merge_script, stageAs: 'merge_species_normalisation.py')

    output:
    path('species_normalisation_summary.tsv'), emit: tsv
    path('species_normalisation_summary.json'), emit: json

    script:
    """
    set -euo pipefail
    python3 merge_species_normalisation.py
    """
}
