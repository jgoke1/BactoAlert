include {
    SNIPPY; SNIPPY_CORE; REFERENCE_PROVENANCE; GUBBINS; MASK_RECOMBINATION
    FILTER_COMPARISON_ALIGNMENT; IQTREE3; FORMAT_TREE
    COMPARISON_STATUS as REFERENCE_STATUS
    COMPARISON_STATUS as GUBBINS_STATUS
    COMPARISON_STATUS as TREE_STATUS
} from '../modules/comparison_tools'
include { DISTANCE_REPORTING } from './distance_reporting'

workflow SNP_PHYLOGENY {
    take:
    groups
    references
    helper
    main:
    // Reference table: one species_slug/group_id and one local FASTA per group.
    matched = groups.map { ctx,g -> tuple(ctx.group_id,ctx,g) }
        .join(references, remainder:true, failOnDuplicate:true)
        .map { row ->
            // An unmatched right-hand join row has only [id, null, ref].
            if (row[1] == null) error "Reference group '${row[0]}' is not among selected species groups."
            tuple(row[1],row[2],row.size() > 3 ? row[3] : null)
        }
    available = matched.filter { ctx,g,ref -> ref != null }
        .map { ctx,g,ref -> tuple(ctx + [reference_source:ref.toString()],g,ref) }
    REFERENCE_STATUS(matched.filter { ctx,g,ref -> ref == null }
        .map { ctx,g,ref -> tuple(ctx,g,'snp_phylogeny','no_local_reference_supplied') },helper)
    REFERENCE_PROVENANCE(available,helper)
    jobs = REFERENCE_PROVENANCE.out.verified
    SNIPPY(jobs.flatMap { ctx,g,ref ->
        def rows = new groovy.json.JsonSlurper().parseText(g.resolve('members.json').text)
        rows.collect { row -> tuple(ctx,row.alias,row.sample_id,g.resolve("genomes/${row.alias}.fasta"),ref) }
    })
    core_inputs = SNIPPY.out.calls.groupTuple()
        .map { ctx,calls -> tuple(ctx.group_id,ctx,calls) }
        .join(jobs.map { ctx,g,ref -> tuple(ctx.group_id,ref) })
        .map { id,ctx,calls,ref -> tuple(ctx,calls,ref) }
    // Snippy runs even for a singleton; only the cohort core step needs >=2.
    SNIPPY_CORE(core_inputs.filter { ctx,calls,ref -> ctx.n_samples >= 2 })
    core_with_group = SNIPPY_CORE.out.results.map { ctx,raw -> tuple(ctx.group_id,ctx,raw) }
        .join(groups.map { ctx,g -> tuple(ctx.group_id,g) })
        .map { id,ctx,raw,g -> tuple(ctx+[method:'snippy_unmasked',route:'snp_phylogeny'],raw.resolve('core.full.aln'),g) }
    FILTER_COMPARISON_ALIGNMENT(core_with_group,helper)
    unmasked = FILTER_COMPARISON_ALIGNMENT.out.results.map { ctx,a ->
        def stats = new groovy.json.JsonSlurper().parseText(a.resolve('status.json').text)
        tuple(ctx + stats,a)
    }
    all_alignments = unmasked
    tree_candidates = Channel.empty()
    extra_results = Channel.empty()
    if (params.snp_recombination) {
        GUBBINS(unmasked.filter { ctx,a -> ctx.gubbins_eligible })
        mask_input = GUBBINS.out.results.map { ctx,gub -> tuple(ctx.group_id,ctx,gub) }
            .join(unmasked.map { ctx,a -> tuple(ctx.group_id,a) })
            .join(groups.map { ctx,g -> tuple(ctx.group_id,g) })
            .map { id,ctx,gub,a,g -> tuple(ctx+[method:'snippy_filtered'],a,gub,g) }
        MASK_RECOMBINATION(mask_input,helper)
        filtered = MASK_RECOMBINATION.out.results.map { ctx,a ->
            def stats = new groovy.json.JsonSlurper().parseText(a.resolve('status.json').text)
            tuple(ctx + stats,a)
        }
        all_alignments = all_alignments.mix(filtered)
        tree_candidates = filtered
        skipped_gubbins = unmasked.filter { ctx,a -> !ctx.gubbins_eligible }
            .map { ctx,a -> tuple(ctx.group_id,ctx) }
            .join(groups.map { ctx,g -> tuple(ctx.group_id,g) })
            .map { id,ctx,g -> tuple(ctx+[method:'gubbins_status'],g,'snp_phylogeny','insufficient_samples_variation_or_callability_for_gubbins_and_filtered_tree') }
        GUBBINS_STATUS(skipped_gubbins,helper)
        extra_results = GUBBINS.out.results.mix(MASK_RECOMBINATION.out.results,GUBBINS_STATUS.out.results)
    } else {
        tree_candidates = unmasked
    }
    DISTANCE_REPORTING(all_alignments,groups,helper)
    IQTREE3(tree_candidates.filter { ctx,a -> ctx.tree_eligible })
    tree_input = IQTREE3.out.results.map { ctx,p -> tuple(ctx.group_id,ctx,p) }
        .join(groups.map { ctx,g -> tuple(ctx.group_id,g) })
        .map { id,ctx,p,g -> tuple(ctx,p,g) }
    FORMAT_TREE(tree_input,helper)
    skipped_trees = tree_candidates.filter { ctx,a -> !ctx.tree_eligible }
        .map { ctx,a -> tuple(ctx.group_id,ctx) }
        .join(groups.map { ctx,g -> tuple(ctx.group_id,g) })
        .map { id,ctx,g -> tuple(ctx+[method:'tree_status'],g,'snp_phylogeny','insufficient_samples_variation_or_callability_for_tree') }
    singleton = core_inputs.filter { ctx,calls,ref -> ctx.n_samples < 2 }
        .map { ctx,calls,ref -> tuple(ctx.group_id,ctx) }
        .join(groups.map { ctx,g -> tuple(ctx.group_id,g) })
        .map { id,ctx,g -> tuple(ctx+[method:'singleton_status'],g,'snp_phylogeny','snippy_calls_complete_no_between_sample_comparison') }
    TREE_STATUS(skipped_trees.mix(singleton),helper)
    emit:
    results = REFERENCE_PROVENANCE.out.results.mix(REFERENCE_STATUS.out.results,SNIPPY.out.versions,SNIPPY_CORE.out.results,
        FILTER_COMPARISON_ALIGNMENT.out.results,DISTANCE_REPORTING.out.results,extra_results,
        IQTREE3.out.results,FORMAT_TREE.out.results,TREE_STATUS.out.results)
        .map { ctx,p -> tuple("${ctx.group_id}/${ctx.method ?: 'reference'}".toString(),p) }
}
