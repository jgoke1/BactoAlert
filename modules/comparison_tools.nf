/* Classic DSL2. Each tool runs in its own public, versioned container. */
process PREPARE_SPECIES_GROUP {
    tag "${ctx.group_id}"
    container params.comparison_python_container
    publishDir "${params.outdir}/comparison_groups/${ctx.group_id}", mode: 'copy'
    input:
    tuple val(ctx), val(metas), path(assemblies, stageAs: 'assemblies/????/input.fasta')
    path(helper, stageAs: 'comparison_io.py')
    output:
    tuple val(ctx), path('group'), emit: groups
    script:
    def paths = assemblies instanceof List ? assemblies : [assemblies]
    def manifest = groovy.json.JsonOutput.toJson(metas.withIndex().collect { m,i ->
        [sample_id:m.sample_id, species:m.species_normalised, input_type:m.type,
         species_confidence:m.speciator_confidence, assembly:paths[i].toString()]
    })
    """
    set -euo pipefail
    cat > manifest.json <<'MANIFEST'
    ${manifest}
    MANIFEST
    python3 comparison_io.py prepare --manifest manifest.json --output group
    """
}

process COMPARISON_STATUS {
    tag "${ctx.group_id}: ${method}: ${reason}"
    container params.comparison_python_container
    publishDir "${params.outdir}/${method}/${ctx.group_id}", mode: 'copy'
    input:
    tuple val(ctx), path(group), val(method), val(reason)
    path(helper, stageAs: 'comparison_io.py')
    output:
    tuple val(ctx), path('not_performed'), emit: results
    script:
    """
    python3 comparison_io.py status --group '${group}' --method '${method}' \\
      --reason '${reason}' --output not_performed
    """
}

process MASHTREE {
    tag "${ctx.group_id}"
    container params.mashtree_container
    cpus { Math.min(params.max_cpus as int, 8) }
    publishDir "${params.outdir}/genome_similarity/${ctx.group_id}", mode: 'copy'
    input:
    tuple val(ctx), path(group, stageAs: 'group')
    output:
    tuple val(ctx), path('mash_raw'), emit: results
    script:
    """
    set -euo pipefail
    mkdir mash_raw
    mashtree --version > mash_raw/mashtree.version.txt 2>&1
    mashtree --numcpus ${task.cpus} --mindepth 0 --seed 42 \\
      --outmatrix mash_raw/matrix.tsv --outtree mash_raw/tree.nwk \\
      group/genomes/*.fasta > mash_raw/run.log 2>&1
    """
}

process SKA_BUILD {
    tag "${ctx.group_id}"
    container params.ska2_container
    cpus { Math.min(params.max_cpus as int, 8) }
    publishDir "${params.outdir}/ska_snp/${ctx.group_id}", mode: 'copy'
    input:
    tuple val(ctx), path(group, stageAs: 'group')
    output:
    tuple val(ctx), path('ska_build'), emit: results
    script:
    """
    set -euo pipefail
    mkdir ska_build
    ska --version > ska_build/ska.version.txt 2>&1
    ska build -f group/ska_input.tsv -k ${params.ska_kmer as int} \\
      --threads ${task.cpus} -o ska_build/species > ska_build/run.log 2>&1
    """
}

process SKA_ALIGN {
    tag "${ctx.group_id}"
    container params.ska2_container
    cpus 1
    publishDir "${params.outdir}/ska_snp/${ctx.group_id}", mode: 'copy'
    input:
    tuple val(ctx), path(build, stageAs: 'build')
    output:
    tuple val(ctx), path('ska_align'), emit: results
    script:
    """
    set -euo pipefail
    mkdir ska_align
    ska align build/species.skf --min-freq 1 --filter no-filter \\
      -o ska_align/full.fasta > ska_align/run.log 2>&1
    """
}

process FILTER_COMPARISON_ALIGNMENT {
    tag "${ctx.group_id}: ${ctx.method}"
    container params.comparison_python_container
    publishDir "${params.outdir}/${ctx.route}/${ctx.group_id}/${ctx.method}", mode: 'copy'
    input:
    tuple val(ctx), path(alignment, stageAs: 'input.fasta'), path(group, stageAs: 'group')
    path(helper, stageAs: 'comparison_io.py')
    output:
    tuple val(ctx), path('alignment'), emit: results
    script:
    """
    python3 comparison_io.py filter --alignment input.fasta --group group \\
      --method '${ctx.method}' --min-fraction ${params.snp_min_core_fraction as double} --output alignment
    """
}

process SNP_DISTS {
    tag "${ctx.group_id}: ${ctx.method}"
    container params.snp_dists_container
    cpus { Math.min(params.max_cpus as int, 8) }
    publishDir "${params.outdir}/${ctx.route}/${ctx.group_id}/${ctx.method}", mode: 'copy'
    input:
    tuple val(ctx), path(alignment, stageAs: 'alignment')
    output:
    tuple val(ctx), path('snp_dists'), emit: results
    script:
    """
    set -euo pipefail
    mkdir snp_dists
    snp-dists -v > snp_dists/snp-dists.version.txt 2>&1
    snp-dists -j ${task.cpus} alignment/core.fasta > snp_dists/matrix.tsv
    """
}

process FORMAT_DISTANCE_REPORT {
    tag "${ctx.group_id}: ${ctx.method}"
    container params.comparison_python_container
    publishDir "${params.outdir}/${ctx.route}/${ctx.group_id}/${ctx.method}", mode: 'copy'
    input:
    tuple val(ctx), path(raw, stageAs: 'raw'), path(group, stageAs: 'group'), path(info, stageAs: 'alignment_info')
    path(helper, stageAs: 'comparison_io.py')
    output:
    tuple val(ctx), path('distance_report'), emit: results
    script:
    def extras = ctx.method == 'mash' ? '--tree raw/tree.nwk' : '--info alignment_info/status.json'
    """
    python3 comparison_io.py report --matrix raw/matrix.tsv --group group \\
      --method '${ctx.method}' ${extras} --output distance_report
    """
}

process SNIPPY {
    tag "${ctx.group_id}: ${sample_id}"
    container params.snippy_container
    cpus { Math.min(params.max_cpus as int, 8) }
    publishDir "${params.outdir}/snp_phylogeny/${ctx.group_id}", mode: 'copy', saveAs: { filename -> filename.endsWith('.version.txt') ? "versions/${filename}" : "calls/${sample_id}" }
    input:
    tuple val(ctx), val(alias), val(sample_id), path(assembly, stageAs: 'assembly.fasta'), path(reference, stageAs: 'reference.fasta')
    output:
    tuple val(ctx), path(alias), emit: calls
    tuple val(ctx), path("${sample_id}.snippy.version.txt"), emit: versions
    script:
    """
    set -euo pipefail
    snippy --outdir '${alias}' --ref reference.fasta --ctgs assembly.fasta --cpus ${task.cpus}
    snippy --version > '${sample_id}.snippy.version.txt' 2>&1
    """
}

process SNIPPY_CORE {
    tag "${ctx.group_id}"
    container params.snippy_container
    publishDir "${params.outdir}/snp_phylogeny/${ctx.group_id}", mode: 'copy'
    input:
    tuple val(ctx), path(calls, stageAs: 'calls/*'), path(reference, stageAs: 'reference.fasta')
    output:
    tuple val(ctx), path('snippy_core'), emit: results
    script:
    """
    set -euo pipefail
    mkdir snippy_core
    snippy-core --ref reference.fasta --prefix snippy_core/core calls/* > snippy_core/run.log 2>&1
    snippy --version > snippy_core/snippy.version.txt 2>&1
    """
}

process REFERENCE_PROVENANCE {
    tag "${ctx.group_id}"
    container params.comparison_python_container
    publishDir "${params.outdir}/snp_phylogeny/${ctx.group_id}", mode: 'copy'
    input:
    tuple val(ctx), path(group), path(reference, stageAs: 'reference.fasta')
    path(helper, stageAs: 'comparison_io.py')
    output:
    tuple val(ctx), path('reference_info'), emit: results
    tuple val(ctx), path(group), path('reference.fasta'), emit: verified
    script:
    def provenance = groovy.json.JsonOutput.toJson([group_id:ctx.group_id,reference_source:ctx.reference_source])
    """
    set -euo pipefail
    mkdir reference_info
    cat > reference_info/source.json <<'MANIFEST'
    ${provenance}
    MANIFEST
    python3 comparison_io.py reference --reference reference.fasta --output reference_info/checksum.json
    """
}

process GUBBINS {
    tag "${ctx.group_id}"
    container params.gubbins_container
    cpus { Math.min(params.max_cpus as int, 8) }
    publishDir "${params.outdir}/snp_phylogeny/${ctx.group_id}", mode: 'copy'
    input:
    tuple val(ctx), path(alignment, stageAs: 'alignment')
    output:
    tuple val(ctx), path('gubbins'), emit: results
    script:
    """
    set -euo pipefail
    export NUMBA_CACHE_DIR="\$PWD/.numba"
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir gubbins
    run_gubbins.py --version > gubbins/gubbins.version.txt 2>&1
    run_gubbins.py --threads ${task.cpus} --prefix gubbins/result \\
      --filter-percentage 100 alignment/full.fasta > gubbins/run.log 2>&1
    """
}

process MASK_RECOMBINATION {
    tag "${ctx.group_id}"
    container params.comparison_python_container
    publishDir "${params.outdir}/snp_phylogeny/${ctx.group_id}/snippy_filtered", mode: 'copy'
    input:
    tuple val(ctx), path(alignment, stageAs: 'alignment'), path(gubbins, stageAs: 'gubbins'), path(group, stageAs: 'group')
    path(helper, stageAs: 'comparison_io.py')
    output:
    tuple val(ctx), path('alignment_filtered'), emit: results
    script:
    """
    python3 comparison_io.py filter --alignment alignment/full.fasta --group group \\
      --gff gubbins/result.recombination_predictions.gff --method snippy_filtered \\
      --min-fraction ${params.snp_min_core_fraction as double} --output alignment_filtered
    """
}

process IQTREE3 {
    tag "${ctx.group_id}: ${ctx.method}"
    container params.iqtree3_container
    cpus { Math.min(params.max_cpus as int, 8) }
    publishDir "${params.outdir}/snp_phylogeny/${ctx.group_id}/${ctx.method}", mode: 'copy'
    input:
    tuple val(ctx), path(alignment, stageAs: 'alignment')
    output:
    tuple val(ctx), path('iqtree3'), emit: results
    script:
    """
    set -euo pipefail
    mkdir iqtree3
    iqtree3 --version > iqtree3/iqtree3.version.txt 2>&1
    iqtree3 -s alignment/core.fasta -m '${ctx.tree_model}' -nt ${task.cpus} \\
      -seed 42 -B 1000 --prefix iqtree3/tree > iqtree3/run.log 2>&1
    """
}

process FORMAT_TREE {
    tag "${ctx.group_id}"
    container params.comparison_python_container
    publishDir "${params.outdir}/snp_phylogeny/${ctx.group_id}/${ctx.method}", mode: 'copy'
    input:
    tuple val(ctx), path(tree, stageAs: 'tree'), path(group, stageAs: 'group')
    path(helper, stageAs: 'comparison_io.py')
    output:
    tuple val(ctx), path('tree_report'), emit: results
    script:
    """
    python3 comparison_io.py tree --tree tree/tree.treefile --group group --output tree_report
    """
}
