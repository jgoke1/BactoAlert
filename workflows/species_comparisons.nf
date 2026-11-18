include { PREPARE_SPECIES_GROUP } from '../modules/comparison_tools'

workflow SPECIES_COMPARISONS {

    take:
    assemblies
    helper

    main:
    selected = assemblies.filter { meta, assembly ->
        def wanted = []

        if (params.comparison_species) {
            wanted = params.comparison_species
                .toString()
                .split(',')
                .collect { value -> value.trim() }
        }

        return wanted.isEmpty() || wanted.contains(meta.species_slug)
    }

    checked = selected.ifEmpty {
        error 'No assemblies match the selected species groups.'
    }

    keyed = checked.map { meta, assembly ->
        tuple(meta.species_slug, meta, assembly)
    }

    grouped = keyed.groupTuple().map { species, metas, genomes ->
        if (!species) {
            error 'An assembly is missing its species_slug.'
        }

        def speciesSlug = species.toString()

        if (!(speciesSlug ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
            error "Unsafe species group '${speciesSlug}'"
        }

        def speciesNames = metas.collect { item ->
            item.species_normalised
        }.unique()

        if (speciesNames.size() != 1) {
            error "Species slug collision: ${speciesSlug}"
        }

        def context = [
            group_id    : speciesSlug,
            species_slug: speciesSlug,
            species     : speciesNames[0],
            n_samples   : metas.size()
        ]

        tuple(context, metas, genomes)
    }

    PREPARE_SPECIES_GROUP(grouped, helper)

    emit:
    groups = PREPARE_SPECIES_GROUP.out.groups
}
