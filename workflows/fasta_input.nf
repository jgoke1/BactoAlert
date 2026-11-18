/* Existing genome assemblies enter here; no assembly from reads is performed. */
include { QUAST } from '../modules/quast'
include { CONTAMINATION_CHECKM } from '../modules/contamination'
include { SPECIATION } from '../modules/speciation'
include { PREPARE_FASTA } from '../modules/fasta_prep'


workflow FASTA_INPUT {
    take:
    fasta_samples
    preparation_script

    main:
    PREPARE_FASTA(fasta_samples, preparation_script)
    QUAST(PREPARE_FASTA.out.assemblies)
    CONTAMINATION_CHECKM(PREPARE_FASTA.out.assemblies)
    SPECIATION(PREPARE_FASTA.out.assemblies)

    emit:
    assemblies = PREPARE_FASTA.out.assemblies
    species_reports = SPECIATION.out.species_report
    quast_reports = QUAST.out.report
    checkm_reports = CONTAMINATION_CHECKM.out
}
