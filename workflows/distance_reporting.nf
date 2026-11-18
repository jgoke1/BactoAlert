include {
    SNP_DISTS
    FORMAT_DISTANCE_REPORT
} from '../modules/comparison_tools'

workflow DISTANCE_REPORTING {

    take:
    alignments
    groups
    helper

    main:
    eligible = alignments.filter { ctx, alignment ->
        ctx.distance_eligible
    }

    SNP_DISTS(eligible)

    alignment_info = eligible.map { ctx, alignment ->
        def comparisonKey = "${ctx.group_id}/${ctx.method}".toString()
        tuple(comparisonKey, ctx, alignment)
    }

    distance_results = SNP_DISTS.out.results.map { ctx, result ->
        def comparisonKey = "${ctx.group_id}/${ctx.method}".toString()
        tuple(comparisonKey, result)
    }

    matched = distance_results.join(alignment_info)

    results_by_group = matched.map {
        comparisonKey, result, ctx, alignment ->
        tuple(ctx.group_id, ctx, result, alignment)
    }

    group_files = groups.map { ctx, group ->
        tuple(ctx.group_id, group)
    }

    combined = results_by_group.combine(group_files, by: 0)

    report_inputs = combined.map {
        groupId, ctx, result, alignment, group ->
        tuple(ctx, result, group, alignment)
    }

    FORMAT_DISTANCE_REPORT(report_inputs, helper)

    emit:
    results = SNP_DISTS.out.results.mix(
        FORMAT_DISTANCE_REPORT.out.results
    )
}
