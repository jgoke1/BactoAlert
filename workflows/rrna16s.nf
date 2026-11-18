include {
    VERIFY_16S_DATABASE
    PREPARE_16S_FASTA
    BLAST_16S
    SUMMARISE_16S
    MERGE_16S
} from '../modules/16srna'

workflow RRNA16S {
    take:
    sequences
    database_dir
    helper_script

    main:
    VERIFY_16S_DATABASE(database_dir, helper_script)
    PREPARE_16S_FASTA(sequences, helper_script)
    // first() supplies one reusable database value to every query task.
    BLAST_16S(PREPARE_16S_FASTA.out.queries, VERIFY_16S_DATABASE.out.verified.first())
    SUMMARISE_16S(BLAST_16S.out.hits, helper_script)
    MERGE_16S(SUMMARISE_16S.out.json.map { meta, report -> report }.collect(), helper_script)

    emit:
    results = SUMMARISE_16S.out.tsv
    json = SUMMARISE_16S.out.json
    raw_hits = BLAST_16S.out.hits
    versions = BLAST_16S.out.versions
    database_info = BLAST_16S.out.database_info
    summary_tsv = MERGE_16S.out.tsv
    summary_json = MERGE_16S.out.json
}
