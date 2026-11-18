process COLLECT_RESULT_FILES {
    tag 'report-ready-results'
    label 'python_container'
    publishDir "${params.outdir}/collected_results", mode: 'copy'

    input:
    val(records)
    path(artifacts, stageAs: 'inputs/????/*')
    path(helper, stageAs: 'collect_results.py')

    output:
    path('report_data'), emit: results

    script:
    def stagedArtifacts = artifacts instanceof List ? artifacts : [artifacts]
    def entries = records.withIndex().collect { row, i ->
        [analysis: row[0], sample_id: row[1], path: stagedArtifacts[i].toString()]
    }
    def manifest = groovy.json.JsonOutput.toJson(entries)
    """
    set -euo pipefail
    cat > inputs.json <<'MANIFEST'
    ${manifest}
    MANIFEST
    python3 collect_results.py --manifest inputs.json --output report_data
    """
}

workflow COLLECT_RESULTS {
    take:
    artifacts
    helper
    main:
    // Collect complete channel events, not the asynchronously published output folder.
    batch = artifacts.toList()
    COLLECT_RESULT_FILES(batch.map { rows -> rows.collect { [it[0], it[1]] } },
                        batch.map { rows -> rows.collect { it[2] } }, helper)
    emit:
    results = COLLECT_RESULT_FILES.out.results
}
