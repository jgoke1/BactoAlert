#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
        IMPORT PLUGINS
========================================================================================
*/

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { startMessage } from './modules/messages'
/*
========================================================================================
    IMPORT SUBWORKFLOWS
========================================================================================
*/

include { SR_ASSEMBLY          } from './workflows/sr_assembly'
include { LR_ASSEMBLY          } from './workflows/lr_assembly'
include { HY_ASSEMBLY          } from './workflows/hybrid_assembly'
include { FASTA_INPUT          } from './workflows/fasta_input'
include { RRNA16S              } from './workflows/rrna16s'
include { CHARACTERISATION     } from './workflows/characterisation'
include { SPECIES_COMPARISONS } from './workflows/species_comparisons'
include { COMPARATIVE_ANALYSIS } from './workflows/comparative_analysis'
include { PREPARE_REPORT_INPUT; BACTOALERT_REPORT } from './modules/bactoalert_report'
include { COLLECT_RESULTS } from './workflows/collect_results'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Classic Groovy helpers; paths are resolved before any sample tasks launch.
def checkedSampleId(value) {
    def id = value?.toString()
    if (!id || !(id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error "Invalid sample ID '${id}': start with a letter or number; use letters, numbers, underscores, dots or hyphens."
    }
    return id
}

def checkUniqueSamples(rows, label) {
    def ids = rows.collect { it[0].sample_id }
    if (ids.size() != ids.toSet().size()) error "Duplicate sample_id values in ${label} inputs."
}

def fastaPathRows(pattern, kind) {
    def resolved = file(pattern, checkIfExists: true)
    def paths = resolved instanceof List ? resolved : [resolved]
    return paths.sort { it.toString() }.collect { fasta ->
        def filename = fasta.getName()
        if (!(filename ==~ /(?i).+\.(fa|fna|fasta)(\.gz)?/)) {
            error "Expected .fa, .fna or .fasta (optionally .gz), got '${filename}'."
        }
        def sampleId = checkedSampleId(filename.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, ''))
        tuple([sample_id: sampleId, type: kind == 'wgs' ? 'fasta' : '16s', input_type: kind, input_format: 'fasta'], fasta)
    }
}

workflow {
    startMessage(workflow.manifest.version)
    // Validate input parameters
    validateParameters()
    log.info paramsSummaryLog(workflow)

    if (!(params.input_type in ['wgs', '16s', 'both'])) {
        error "Unsupported --input_type '${params.input_type}'; choose wgs, 16s or both."
    }
    if (!(params.input_format in ['fastq', 'fasta'])) {
        error "Unsupported --input_format '${params.input_format}'; choose fastq or fasta."
    }
    if (!params.run_id) error 'Supply --run_id with the sequencing run identifier.'
    def reportId = checkedSampleId(params.run_id)
    def snapshot = file(params.qualibact_dir, checkIfExists: true)
    def contextDir = params.report_reference_dir ? file(params.report_reference_dir, checkIfExists: true) : []
    def runWgs = params.input_type in ['wgs', 'both']
    def run16s = params.input_type in ['16s', 'both']
    if (!runWgs && (params.samplesheet || params.fasta_samplesheet || params.fasta)) {
        error 'For 16S-only runs use --rrna16s_samplesheet or --rrna16s_fasta. WGS inputs need --input_type both.'
    }
    if (!run16s && (params.rrna16s_samplesheet || params.rrna16s_fasta)) {
        error '16S inputs were supplied; choose --input_type both or --input_type 16s.'
    }

    def readRows = []
    def assemblyRows = []
    if (runWgs) {
        // Keep the legacy fallback only when no explicit WGS assembly input is given.
        def sheet = params.samplesheet ?: ((params.fasta_samplesheet || params.fasta) ? null : "${launchDir}/samplesheet.csv")
        if (sheet && params.input_format == 'fastq') {
            def rows = samplesheetToList(sheet, 'assets/schema_input.json')
            if (!rows) error 'The FASTQ samplesheet contains no samples.'
            readRows = rows.collect { meta, short_r1, short_r2, long_reads ->
                def sampleId = checkedSampleId(meta.sample_id)
                if ((short_r1 as boolean) != (short_r2 as boolean)) {
                    error "Sample '${sampleId}': provide both short_reads1 and short_reads2, or leave both empty."
                }
                if (!short_r1 && !long_reads) {
                    error "Sample '${sampleId}': no FASTQ input; provide paired short reads or long_reads. For assemblies use --input_format fasta."
                }
                tuple(meta + [sample_id: sampleId, input_type: 'wgs', input_format: 'fastq'], short_r1, short_r2, long_reads)
            }
        } else if (sheet) {
            assemblyRows = samplesheetToList(sheet, 'assets/schema_input_fasta.json').collect { meta, fasta ->
                tuple(meta + [sample_id: checkedSampleId(meta.sample_id), type: 'fasta', input_type: 'wgs', input_format: 'fasta'], fasta)
            }
        }
        // A separate assembly sheet can accompany the FASTQ sheet in the same run.
        // --input_format applies only to --samplesheet; this sheet is always FASTA.
        if (params.fasta_samplesheet) {
            def rows = samplesheetToList(params.fasta_samplesheet, 'assets/schema_input_fasta.json')
            if (!rows) error 'The FASTA samplesheet contains no samples.'
            assemblyRows.addAll(rows.collect { meta, fasta ->
                tuple(meta + [sample_id: checkedSampleId(meta.sample_id), type: 'fasta', input_type: 'wgs', input_format: 'fasta'], fasta)
            })
        }
        if (params.fasta) assemblyRows.addAll(fastaPathRows(params.fasta, 'wgs'))
        if (!readRows && !assemblyRows) error 'No WGS samples supplied.'
        checkUniqueSamples(readRows + assemblyRows, 'WGS')
    }

    collected_artifacts = Channel.empty()
    def rrnaRows = []
    if (run16s) {
        if (params.rrna16s_samplesheet && params.rrna16s_fasta) {
            error 'Choose one 16S input source: --rrna16s_samplesheet OR --rrna16s_fasta.'
        }
        if (params.rrna16s_samplesheet) {
            rrnaRows = samplesheetToList(params.rrna16s_samplesheet, 'assets/schema_input_16s.json').collect { meta, fasta ->
                tuple(meta + [sample_id: checkedSampleId(meta.sample_id), type: '16s', input_type: '16s', input_format: 'fasta'], fasta)
            }
        } else if (params.rrna16s_fasta) {
            rrnaRows = fastaPathRows(params.rrna16s_fasta, '16s')
        } else {
            error 'Supply --rrna16s_samplesheet or --rrna16s_fasta for the 16S workflow.'
        }
        if (!rrnaRows) error 'The 16S input contains no samples.'
        checkUniqueSamples(rrnaRows, '16S')
        // No dependency on assembly, Speciator or characterisation outputs.
        RRNA16S(
            Channel.fromList(rrnaRows),
            file(params.rrna16s_db, checkIfExists: true),
            file("${projectDir}/bin/rrna16s.py", checkIfExists: true)
        )
    }

    if (run16s) {
        collected_artifacts = collected_artifacts.mix(
            RRNA16S.out.summary_tsv.map { p -> tuple('rrna16s_summary','cohort',p) },
            RRNA16S.out.results.map { meta, result -> tuple('16s', meta.sample_id, result) },
            RRNA16S.out.versions.map { meta, result -> tuple('16s_versions', meta.sample_id, result) },
            RRNA16S.out.database_info.map { meta, result -> tuple('16s_database', meta.sample_id, result) }
        )
    }

    if (runWgs) {
        // These are validated Lists, so selection is known before channel creation.
        def shortRows = readRows.findAll { meta, r1, r2, long_reads ->
            r1 && r2 && !long_reads
        }.collect { meta, r1, r2, long_reads ->
            tuple(meta + [type: 'short'], r1, r2)
        }
        def longRows = readRows.findAll { meta, r1, r2, long_reads ->
            long_reads && !r1 && !r2
        }.collect { meta, r1, r2, long_reads ->
            tuple(meta + [type: 'long'], long_reads)
        }
        def hybridRows = readRows.findAll { meta, r1, r2, long_reads ->
            r1 && r2 && long_reads
        }.collect { meta, r1, r2, long_reads ->
            tuple(meta + [type: 'hybrid'], r1, r2, long_reads)
        }

        log.info "WGS inputs: short=${shortRows.size()}, long=${longRows.size()}, hybrid=${hybridRows.size()}, fasta=${assemblyRows.size()}"

        // Only invoked workflows register their processes in the progress display.
        // Keep empty placeholders so skipped workflows need no .out references.
        all_assemblies = Channel.empty()
        all_species_reports = Channel.empty()

        if (shortRows) {
            SR_ASSEMBLY(Channel.fromList(shortRows))
            all_assemblies = all_assemblies.mix(SR_ASSEMBLY.out.assemblies)
            all_species_reports = all_species_reports.mix(SR_ASSEMBLY.out.species_reports)
                collected_artifacts = collected_artifacts.mix(
                    SR_ASSEMBLY.out.quast_reports.map { m, p -> tuple('quast', m.sample_id, p) },
                    SR_ASSEMBLY.out.checkm_reports.map { m, p -> tuple('checkm', m.sample_id, p) })
        }

        if (longRows) {
            LR_ASSEMBLY(Channel.fromList(longRows))
            all_assemblies = all_assemblies.mix(LR_ASSEMBLY.out.assemblies)
            all_species_reports = all_species_reports.mix(LR_ASSEMBLY.out.species_reports)
                collected_artifacts = collected_artifacts.mix(
                    LR_ASSEMBLY.out.quast_reports.map { m, p -> tuple('quast', m.sample_id, p) },
                    LR_ASSEMBLY.out.checkm_reports.map { m, p -> tuple('checkm', m.sample_id, p) })
        }

        if (hybridRows) {
            hybSplit = Channel.fromList(hybridRows).multiMap { meta, r1, r2, long_reads ->
                shortreads: tuple(meta, r1, r2)
                longreads: tuple(meta, long_reads)
            }
            HY_ASSEMBLY(hybSplit.shortreads, hybSplit.longreads)
            all_assemblies = all_assemblies.mix(HY_ASSEMBLY.out.assemblies)
            all_species_reports = all_species_reports.mix(HY_ASSEMBLY.out.species_reports)
                collected_artifacts = collected_artifacts.mix(
                    HY_ASSEMBLY.out.quast_reports.map { m, p -> tuple('quast', m.sample_id, p) },
                    HY_ASSEMBLY.out.checkm_reports.map { m, p -> tuple('checkm', m.sample_id, p) })
        }

        if (assemblyRows) {
            // Existing assemblies get assembly QC and a new Speciator report.
            FASTA_INPUT(
                Channel.fromList(assemblyRows),
                file("${projectDir}/bin/prepare_fasta.py", checkIfExists: true)
            )
            all_assemblies = all_assemblies.mix(FASTA_INPUT.out.assemblies)
            all_species_reports = all_species_reports.mix(FASTA_INPUT.out.species_reports)
                collected_artifacts = collected_artifacts.mix(
                    FASTA_INPUT.out.quast_reports.map { m, p -> tuple('quast', m.sample_id, p) },
                    FASTA_INPUT.out.checkm_reports.map { m, p -> tuple('checkm', m.sample_id, p) })
        }

        // Characterise the enabled branches using their paired Speciator reports.
        CHARACTERISATION(
        all_assemblies,
        all_species_reports,
        file("${projectDir}/bin/normalise_species.py", checkIfExists: true),
        file("${projectDir}/bin/virulencefinder_json_tsv.py", checkIfExists: true),
        file("${projectDir}/bin/plasmid_json_tsv.py", checkIfExists: true),
        file("${projectDir}/bin/merge_species_normalisation.py", checkIfExists: true)
        )
        if (params.run_genome_similarity || params.run_ska_snp || params.run_snp_phylogeny) {
            SPECIES_COMPARISONS(CHARACTERISATION.out.classified_assemblies,
                               file("${projectDir}/bin/comparison_io.py",checkIfExists:true))
            COMPARATIVE_ANALYSIS(SPECIES_COMPARISONS.out.groups,
                                 file("${projectDir}/bin/comparison_io.py",checkIfExists:true))
            collected_artifacts = collected_artifacts.mix(COMPARATIVE_ANALYSIS.out.results)
        }
            collected_artifacts = collected_artifacts.mix(
                all_species_reports.map { m, p -> tuple('speciator', m.sample_id, p) },
                CHARACTERISATION.out.kleborate_results.map { m, p -> tuple('kleborate', m.sample_id, p) },
                CHARACTERISATION.out.kleborate_versions.map { m, p -> tuple('kleborate_versions', m.sample_id, p) },
                CHARACTERISATION.out.general_amr_results.map { m, p -> tuple('amrfinder', m.sample_id, p) },
                CHARACTERISATION.out.general_amr_versions.map { m, p -> tuple('amrfinder_versions', m.sample_id, p) },
                CHARACTERISATION.out.general_mlst_results.map { m, p -> tuple('mlst', m.sample_id, p) },
                CHARACTERISATION.out.general_virulence_results.map { m, p -> tuple('virulencefinder', m.sample_id, p) },
                CHARACTERISATION.out.general_plasmid_results.map { m, p -> tuple('plasmidfinder', m.sample_id, p) },
                CHARACTERISATION.out.species_summary_tsv.map { p -> tuple('taxonomy', 'cohort', p) },
                CHARACTERISATION.out.species_summary_json.map { p -> tuple('taxonomy_json', 'cohort', p) }
            )
    }
    if (params.collect_results) {
        COLLECT_RESULTS(collected_artifacts, file("${projectDir}/bin/collect_results.py"))
    }
    report_batch = collected_artifacts.toList().map { rows ->
        def needed = rows.findAll { kind, id, artifact ->
            !(kind in ['genome_similarity','ska_snp','snp_phylogeny']) ||
                (kind != 'snp_phylogeny' && artifact.name in ['distance_report','not_performed'])
        }
        tuple(needed.collect { [it[0], it[1]] }, needed.collect { it[2] })
    }
    PREPARE_REPORT_INPUT(report_batch, file("${projectDir}/bin/prepare_report_input.py"))
    BACTOALERT_REPORT(reportId, PREPARE_REPORT_INPUT.out.bundle, snapshot, contextDir)

}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/
 workflow.onComplete {
        workDir = new File("${workflow.workDir}")

        println """
        BACTOALERT Execution Summary
        ---------------------------
        Pipeline Version : ${workflow.manifest.version}
        Nextflow Version : ${nextflow.version}
        Command Line     : ${workflow.commandLine}
        Resumed          : ${workflow.resume}
        Completed At     : ${workflow.complete}
        Duration         : ${workflow.duration}
        Success          : ${workflow.success}
        Exit Code        : ${workflow.exitStatus}
        Error Report     : ${workflow.errorReport ?: '-'}
        Launch Dir       : ${workflow.launchDir}
        """
    }

/*
========================================================================================
    THE END
========================================================================================
*/
