nextflow.enable.dsl=2

// ── Parameters ────────────────────────────────────────────────────────────────

params.node_names      = null
params.ids             = "${projectDir}/ids.txt"
params.OUTDIR          = "${projectDir}/results"
params.SPECIES_TREE    = "${projectDir}/data/species_tree.full.newick"
params.REFSPECIES      = "Mmus"
params.REFNAMES        = "${projectDir}/data/Mmus_gene_names.csv"
params.gene_trees_dir  = "${params.OUTDIR}/gene_trees"

if (!params.node_names) {
    error "ERROR: --node_names is required. Provide a comma-separated list of " +
          "clade names ordered broad→narrow. Example: --node_names Metazoa,Bilateria,Vertebrata"
}

// ── Channels ─────────────────────────────────────────────────────────────────

// Ordered clade levels (broad → narrow), as a single string value channel
levels_val = Channel.value(params.node_names)

// One entry per clade
node_names_ch = Channel.from(
    params.node_names.split(',').collect { it.trim() }.findAll { it }
)

// HG IDs
hg_ids_ch = Channel
    .fromPath(params.ids)
    .splitText()
    .map  { it.trim() }
    .filter { it }

// Static files as value channels
species_tree_val = Channel.value( file(params.SPECIES_TREE) )
refnames_val     = Channel.value( file(params.REFNAMES) )

// ── Process 1 — Extract clade species lists and pruned subtree ────────────────

process EXTRACT_CLADE {

    tag "${node}"

    publishDir "${params.OUTDIR}/ancestry/${node}", mode: 'copy'

    input:
    tuple val(node), path(species_tree)

    output:
    tuple val(node),
          path("${node}.in_species.txt"),
          path("${node}.ignore_species.txt"),
          path("${node}.pruned.tree")

    script:
    """
    python ${projectDir}/workflow/extract_clade.py \
        --tree   ${species_tree} \
        --node   "${node}" \
        --out_prefix ${node}
    """
}

// ── Process 2 — POSSVM per clade × HG ────────────────────────────────────────
// Runs with the in-group species defined for each clade, so OG assignments
// reflect the diversity at that particular taxonomic level.

process PVM_CLADE {

    tag "${node}/${id}"

    publishDir "${params.OUTDIR}/ancestry/${node}/possvm", mode: 'copy'

    cpus 1
    memory { 500.MB * task.attempt }
    time   { 10.min * task.attempt }

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(node), val(id), path(treefile), path(ignore_species), path(refnames)

    output:
    // Emit id as well so it can be used as a grouping key in LINK_HOGS
    tuple val(node), val(id), path("${node}.${id}.ortholog_groups.csv"), optional: false

    script:
    """
    python ${projectDir}/phylogeny/main.py possvm \
        -t          ${treefile} \
        --refsps    ${params.REFSPECIES} \
        -r          ${refnames} \
        --outgroup  ${ignore_species} \
        -o          ${node}.${id}. \
        -p          ${node}.${id} \
    """
}

// ── Process 3 — Link OGs across clade levels for one HG ──────────────────────
// Groups all per-clade POSSVM outputs for a given HG and builds the
// parent→child OG relationship table based on shared gene membership.
// A Metazoa OG splitting into 4 Vertebrata OGs (WGD) will appear as
// one parent row with four children.

process LINK_HOGS {

    tag "${id}"

    publishDir "${params.OUTDIR}/ancestry/hog_links", mode: 'copy'

    input:
    tuple val(id), val(nodes), path(csvs)
    val(levels_str)

    output:
    tuple path("${id}.og_links.tsv"), path("${id}.og_stats.tsv")

    script:
    """
    python ${projectDir}/workflow/link_hog_levels.py \
        --hg     '${id}' \
        --levels '${levels_str}' \
        --csvs   ${csvs} \
        --output_links ${id}.og_links.tsv \
        --output_stats ${id}.og_stats.tsv
    """
}

// ── Process 4 — Generate interactive Sankey hierarchy visualisation ───────────
// Collects all per-HG link/stats TSVs and produces a single self-contained
// HTML report showing OG flow from the broadest to the narrowest clade level.

process VISUALIZE_HIERARCHY {

    publishDir "${params.OUTDIR}/ancestry", mode: 'copy'

    input:
    path(all_links)
    path(all_stats)
    val(levels_str)

    output:
    path("hog_hierarchy.html")

    script:
    """
    python ${projectDir}/workflow/visualize_hog_hierarchy.py \
        --links  ${all_links} \
        --stats  ${all_stats} \
        --levels '${levels_str}' \
        --output hog_hierarchy.html
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow {

    // Step 1: extract per-clade species lists and pruned trees
    clade_info = node_names_ch
        .combine(species_tree_val)
        | EXTRACT_CLADE
    // clade_info: (node, in_species, ign_species, pruned_tree)

    // Step 2: run POSSVM for each clade × HG combination
    pvm_input = clade_info
        .map   { node, in_sp, ign_sp, pruned -> tuple(node, ign_sp) }
        .combine(hg_ids_ch)
        .combine(refnames_val)
        .map { node, ign_sp, id, refnames ->
            def treefile = file("${params.gene_trees_dir}/${id}.treefile")
            tuple(node, id, treefile, ign_sp, refnames)
        }
        .filter { node, id, treefile, ign_sp, refnames -> treefile.exists() }

    pvm_out = pvm_input | PVM_CLADE
    // pvm_out: (node, id, og_csv)

    // Step 3: group all clade-level results for the same HG and link OGs
    link_in = pvm_out
        .map   { node, id, csv -> tuple(id, node, csv) }
        .groupTuple(by: 0)
    // link_in: (id, [node1, node2, ...], [csv1, csv2, ...])

    link_out = LINK_HOGS(link_in, levels_val)
    // link_out: (links_tsv, stats_tsv)

    // Step 4: collect everything and build the visualisation
    all_links = link_out.map { links, stats -> links }.collect()
    all_stats = link_out.map { links, stats -> stats }.collect()

    VISUALIZE_HIERARCHY(all_links, all_stats, levels_val)
}
