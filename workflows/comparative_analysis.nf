include { GENOME_SIMILARITY } from './genome_similarity'
include { SKA_SNP } from './ska_snp'
include { SNP_PHYLOGENY } from './snp_phylogeny'
workflow COMPARATIVE_ANALYSIS {
    take:
    groups
    helper
    main:
    if (params.run_ska_snp && ((params.ska_kmer as int) < 5 || (params.ska_kmer as int) > 63 || (params.ska_kmer as int) % 2 == 0)) {
        error '--ska_kmer must be an odd integer from 5 to 63.'
    }
    artifacts = Channel.empty()
    if (params.run_genome_similarity) {
        GENOME_SIMILARITY(groups,helper)
        artifacts = artifacts.mix(GENOME_SIMILARITY.out.results.map { id,p -> tuple('genome_similarity',id,p) })
    }
    if (params.run_ska_snp) {
        SKA_SNP(groups,helper)
        artifacts = artifacts.mix(SKA_SNP.out.results.map { id,p -> tuple('ska_snp',id,p) })
    }
    if (params.run_snp_phylogeny) {
        if (params.snp_groups) error 'Use --snp_references references.csv with species_slug,reference columns; grouping now comes from species channels.'
        refs = Channel.empty()
        if (params.snp_references) {
            refs = Channel.fromPath(params.snp_references,checkIfExists:true).splitCsv(header:true,strip:true)
                .map { row ->
                    def id = row.species_slug?.toString()
                    if (!id || !(id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Reference CSV requires a safe species_slug column.'
                    if (!row.reference) error "Missing reference FASTA for ${id}"
                    tuple(id,file(row.reference,checkIfExists:true))
                }
        }
        SNP_PHYLOGENY(groups,refs,helper)
        artifacts = artifacts.mix(SNP_PHYLOGENY.out.results.map { id,p -> tuple('snp_phylogeny',id,p) })
    }
    emit:
    results = artifacts
}
