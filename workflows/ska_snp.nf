include { SKA_BUILD; SKA_ALIGN; FILTER_COMPARISON_ALIGNMENT; COMPARISON_STATUS } from '../modules/comparison_tools'
include { DISTANCE_REPORTING } from './distance_reporting'
workflow SKA_SNP {
    take:
    groups
    helper
    main:
    eligible = groups.filter { ctx,g -> ctx.n_samples >= 2 }
    SKA_BUILD(eligible)
    SKA_ALIGN(SKA_BUILD.out.results)
    alignment_input = SKA_ALIGN.out.results.map { ctx,raw -> tuple(ctx.group_id,ctx,raw) }
        .join(groups.map { ctx,g -> tuple(ctx.group_id,g) })
        .map { id,ctx,raw,g -> tuple(ctx + [method:'ska2',route:'ska_snp'],raw.resolve('full.fasta'),g) }
    FILTER_COMPARISON_ALIGNMENT(alignment_input,helper)
    checked = FILTER_COMPARISON_ALIGNMENT.out.results.map { ctx,a ->
        def stats = new groovy.json.JsonSlurper().parseText(a.resolve('status.json').text)
        tuple(ctx + stats,a)
    }
    DISTANCE_REPORTING(checked,groups,helper)
    COMPARISON_STATUS(groups.filter { ctx,g -> ctx.n_samples < 2 }
        .map { ctx,g -> tuple(ctx,g,'ska_snp','fewer_than_two_samples') },helper)
    emit:
    results = SKA_BUILD.out.results.mix(SKA_ALIGN.out.results,FILTER_COMPARISON_ALIGNMENT.out.results,
        DISTANCE_REPORTING.out.results,COMPARISON_STATUS.out.results)
        .map { ctx,p -> tuple(ctx.group_id,p) }
}
