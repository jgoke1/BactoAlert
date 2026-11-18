include { MASHTREE; FORMAT_DISTANCE_REPORT; COMPARISON_STATUS } from '../modules/comparison_tools'
workflow GENOME_SIMILARITY {
    take:
    groups
    helper
    main:
    eligible = groups.filter { ctx, group -> ctx.n_samples >= 2 }
    small = groups.filter { ctx, group -> ctx.n_samples < 2 }
    MASHTREE(eligible)
    by_id = groups.map { ctx, group -> tuple(ctx.group_id, group) }
    formatted_input = MASHTREE.out.results.map { ctx, raw -> tuple(ctx.group_id,ctx,raw) }
        .join(by_id)
        .map { id, ctx, raw, group -> tuple(ctx + [method:'mash',route:'genome_similarity'], raw, group, group) }
    FORMAT_DISTANCE_REPORT(formatted_input, helper)
    COMPARISON_STATUS(small.map { ctx,g -> tuple(ctx,g,'genome_similarity','fewer_than_two_samples') },helper)
    emit:
    results = MASHTREE.out.results.mix(FORMAT_DISTANCE_REPORT.out.results, COMPARISON_STATUS.out.results)
        .map { ctx, result -> tuple(ctx.group_id,result) }
}
