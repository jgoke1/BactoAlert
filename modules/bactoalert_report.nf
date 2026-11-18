process PREPARE_REPORT_INPUT {
    tag 'cohort report tables'
    label 'python_container'
    input:
    tuple val(rows), path(artifacts, stageAs: 'items/????????/*')
    path helper, stageAs: 'prepare_report_input.py'
    output:
    path 'pipeline_input', emit: bundle
    script:
    def paths = artifacts instanceof List ? artifacts : [artifacts]
    def records = rows.withIndex().collect { r,i -> [r[0], r[1], paths[i].toString()] }
    def payload = groovy.json.JsonOutput.toJson(records)
    """
    cat > records.json <<'RECORDS'
    ${payload}
    RECORDS
    python3 prepare_report_input.py records.json pipeline_input
    """
}

process BACTOALERT_REPORT {
    tag "${run_id}"
    container 'jgoke1/bactoalert-report:0.1.0'
    publishDir "${params.outdir}/reports", mode: 'copy', saveAs: { name -> run_id.toString() }
    input:
    val run_id
    path bundle, stageAs: 'pipeline_input'
    path snapshot, stageAs: 'snapshot'
    path references, stageAs: 'reference'
    output:
    path 'report', emit: bundle
    script:
    def context = references ? '--reference reference' : ''
    """
    bactoalert-report --input pipeline_input --output report --run-id '${run_id}' --qualibact snapshot ${context}
    """
}
