include {
    AMRFINDERPLUS
    MLST
    VIRULENCEFINDER
    PLASMIDFINDER
} from '../modules/general_bacteria'

workflow GENERAL_CHARACTERISATION {
    take:
    classified_assemblies
    virulence_converter
    plasmid_converter

    main:
    // classified_assemblies items are tuple(meta, assembly_fasta).
    // Keep the metadata object attached to every result from every tool.
    AMRFINDERPLUS(classified_assemblies)
    MLST(classified_assemblies)
    VIRULENCEFINDER(classified_assemblies, virulence_converter)
    PLASMIDFINDER(classified_assemblies, plasmid_converter)

    emit:
    amr_results       = AMRFINDERPLUS.out.results
    amr_versions      = AMRFINDERPLUS.out.versions
    mlst_results      = MLST.out.results
    virulence_results = VIRULENCEFINDER.out.tsv
    virulence_json    = VIRULENCEFINDER.out.json
    virulence_raw     = VIRULENCEFINDER.out.raw_dir
    plasmid_results   = PLASMIDFINDER.out.tsv
    plasmid_json      = PLASMIDFINDER.out.json
    plasmid_raw       = PLASMIDFINDER.out.raw_dir
}
