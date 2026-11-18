include { NORMALISE_SPECIES } from '../modules/normalise_species'
include { KLEBORATE } from '../modules/kleborate'
include { GENERAL_CHARACTERISATION } from '../subworkflows/general_characterisation.nf'
include { GROUP_SPECIES_NORMALISATION; MERGE_SPECIES_NORMALISATION } from '../modules/species_normalisation'

workflow CHARACTERISATION {
    take:
    assemblies
    species_reports
    normalisation_script
    virulence_converter
    plasmid_converter
    merge_script

    main:
    // Match reports and assemblies by sample ID; never rely on completion order.
    assemblies_by_id = assemblies.map { meta, assembly ->
        def sampleId = meta.sample_id?.toString()
        if (!sampleId || !(sampleId ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
            error "Invalid sample ID '${sampleId}': use letters, numbers, underscores, dots or hyphens; start with a letter or number."
        }
        tuple(sampleId, meta, assembly)
    }
    reports_by_id = species_reports.map { meta, report ->
        tuple(meta.sample_id.toString(), report)
    }

    normalisation_input = assemblies_by_id
        .join(reports_by_id, failOnDuplicate: true, failOnMismatch: true)
        .map { sample_id, meta, assembly, report ->
            tuple(meta, assembly, report)
        }

    NORMALISE_SPECIES(normalisation_input, normalisation_script)

    normalised_taxonomy = NORMALISE_SPECIES.out.classified.map { meta, assembly, taxonomy_json ->
        def taxonomy = new groovy.json.JsonSlurper().parseText(taxonomy_json.text)
        tuple(meta + taxonomy, assembly, taxonomy_json)
    }

    // Save one uniquely named JSON per sample under a folder named for its species.
    GROUP_SPECIES_NORMALISATION(
        normalised_taxonomy.map { meta, assembly, taxonomy_json -> tuple(meta, taxonomy_json) }
    )

    // Build one cohort-level TSV and JSON after every sample has normalised.
    MERGE_SPECIES_NORMALISATION(
        GROUP_SPECIES_NORMALISATION.out.grouped.map { meta, taxonomy_json -> taxonomy_json }.collect(),
        merge_script
    )

    classified = normalised_taxonomy.map { meta, assembly, taxonomy_json -> tuple(meta, assembly) }

    routes = classified.branch { meta, assembly ->
        // The normaliser maps Escherichia and Shigella calls to this preset.
        kleborate: meta.route in ['kpsc', 'kosc', 'escherichia']
        general: true
    }

    kleborate_input = routes.kleborate.map { meta, assembly ->
        def preset = meta.route.toString()
        tuple(meta + [kleborate_preset: preset], assembly)
    }

    KLEBORATE(kleborate_input)
    GENERAL_CHARACTERISATION(
        routes.general,
        virulence_converter,
        plasmid_converter
    )

    emit:
    classified_assemblies = classified
    kleborate_results = KLEBORATE.out.results
    kleborate_versions = KLEBORATE.out.versions
    general_amr_results = GENERAL_CHARACTERISATION.out.amr_results
    general_amr_versions = GENERAL_CHARACTERISATION.out.amr_versions
    general_mlst_results = GENERAL_CHARACTERISATION.out.mlst_results
    general_virulence_results = GENERAL_CHARACTERISATION.out.virulence_results
    general_virulence_json = GENERAL_CHARACTERISATION.out.virulence_json
    general_virulence_raw = GENERAL_CHARACTERISATION.out.virulence_raw
    general_plasmid_results = GENERAL_CHARACTERISATION.out.plasmid_results
    general_plasmid_json = GENERAL_CHARACTERISATION.out.plasmid_json
    general_plasmid_raw = GENERAL_CHARACTERISATION.out.plasmid_raw
    species_grouped = GROUP_SPECIES_NORMALISATION.out.grouped
    species_summary_tsv = MERGE_SPECIES_NORMALISATION.out.tsv
    species_summary_json = MERGE_SPECIES_NORMALISATION.out.json
}
