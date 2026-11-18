process KLEBORATE {
    tag "${meta.sample_id} (${meta.kleborate_preset})"
    label 'kleborate_container'
    label 'process_medium'

    publishDir "${params.outdir}/kleborate", mode: 'copy', saveAs: { filename ->
        // Preserve the raw directory channel; publish separate file copies below.
        if (filename.toString() == 'kleborate_results') return null

        def folder = meta.kleborate_preset.toString() == 'escherichia'
            ? 'Escherichia_output'
            : 'Klebsiella_output'
        def name = filename.toString().tokenize('/').last()
        def prefix = "${meta.sample_id}_"
        if (!name.startsWith(prefix)) name = "${prefix}${name}"

        "${folder}/${name}"
    }

    input:
    tuple val(meta), path(assembly, stageAs: 'assembly.fasta')

    output:
    tuple val(meta), path('kleborate_results'), emit: results
    tuple val(meta), path('kleborate.version.txt'), emit: versions
    tuple val(meta), path('kleborate_publish/*', type: 'file'), emit: report_files

    script:
    def sampleId = meta.sample_id?.toString()
    if (!sampleId || !(sampleId ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error "Invalid sample ID '${sampleId}' for Kleborate."
    }
    def preset = meta.kleborate_preset?.toString()
    if (!(preset in ['kpsc', 'kosc', 'escherichia'])) {
        error "Unsupported Kleborate preset '${preset}' for ${meta.sample_id}"
    }
    """
    set -euo pipefail
    # Keep Matplotlib configuration/cache writable for Docker's task user.
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$MPLCONFIGDIR"

    kleborate --version > kleborate.version.txt 2>&1
    kleborate \\
        --assemblies assembly.fasta \\
        --outdir kleborate_results \\
        --preset ${preset} \\
        --trim_headers

    # Name the per-sample reports consistently across Kleborate presets.
    shopt -s nullglob
    reports=()
    for report in kleborate_results/*_output.txt; do
        case "\$report" in
            *_hAMRonization_output.txt) continue ;;
            *) reports+=("\$report") ;;
        esac
    done

    if [[ \${#reports[@]} -ne 1 ]]; then
        echo "ERROR: expected one main Kleborate report for ${sampleId}; found \${#reports[@]}" >&2
        exit 1
    fi

    main_report="\${reports[0]}"
    amr_report="\${main_report%_output.txt}_hAMRonization_output.txt"

    mv -- "\$main_report" "kleborate_results/${sampleId}_kleborate.txt"
    if [[ -f "\$amr_report" ]]; then
        mv -- "\$amr_report" "kleborate_results/${sampleId}_hAMRonization.txt"
    fi

    # Separate copies avoid overlapping directory/file outputs in publishDir.
    mkdir -p kleborate_publish
    for report in kleborate_results/*; do
        [[ -f "\$report" ]] || continue
        cp -- "\$report" kleborate_publish/
    done
    """
}
