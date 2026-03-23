nextflow.enable.dsl=2

// ── Parameters ────────────────────────────────────────────────────────────────

params.node_names      = null
params.ids             = "${projectDir}/ids.txt"
params.OUTDIR          = "${projectDir}/results"
params.SPECIES_TREE    = "${projectDir}/data/species_tree.full.newick"
params.REFSPECIES      = "Mmus"
params.REFNAMES        = "${projectDir}/data/Mmus_gene_names.csv"
params.gene_trees_dir  = "${params.OUTDIR}/gene_trees"
params.min_presence    = 2

if (!params.node_names) {
    error "ERROR: --node_names is required. Example: --node_names Bilateria"
}

// ── Channels ─────────────────────────────────────────────────────────────────

// One value per clade
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

// ── Process 1 — Extract clade species lists and pruned tree ──────────────────

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

// ── Process 2 — POSSVM per clade × HG ─────────────────────────────────────────

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
    tuple val(node), path("${node}.${id}.ortholog_groups.csv"), optional: false

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

// ── Process 3 — Build presence/absence matrix ────────────────────────────────

process BUILD_PAM {

    tag "${node}"

    publishDir "${params.OUTDIR}/ancestry/${node}", mode: 'copy'

    input:
    tuple val(node), path(csvs), path(in_species)

    output:
    tuple val(node), path("${node}.pam.tsv")

    script:
    """
    python ${projectDir}/workflow/build_pam.py \
        --csvs ${csvs} \
        --species ${in_species} \
        --min_presence ${params.min_presence} \
        --output ${node}.pam.tsv
    """
}

// ── Process 4 — Ancestral reconstruction ─────────────────────────────────────

process ANCESTRAL_RECON {

    tag "${node}"

    publishDir "${params.OUTDIR}/ancestry/${node}", mode: 'copy'

    memory '4 GB'
    time   '4 h'

    input:
    tuple val(node), path(pam), path(pruned_tree)

    output:
    tuple val(node),
          path("${node}.ancestral_states.tsv"),
          path("${node}.node_probs.tsv"),
          path(pruned_tree)

    script:
    """
    python ${projectDir}/workflow/ancestral_reconstruction.py \
        --pam        ${pam} \
        --tree       ${pruned_tree} \
        --output     ${node}.ancestral_states.tsv \
        --node_probs ${node}.node_probs.tsv \
        --node       ${node}
    """
}

// ── Process 5 — Visualisation ────────────────────────────────────────────────

process VISUALIZE {

    tag "${node}"

    publishDir "${params.OUTDIR}/ancestry/${node}", mode: 'copy'

    input:
    tuple val(node),
          path(states),
          path(node_probs),
          path(pam),
          path(pruned_tree)

    output:
    path("${node}.html")

    script:
    """
    python ${projectDir}/workflow/visualize_ancestry.py \
        --tree       ${pruned_tree} \
        --node_probs ${node_probs} \
        --states     ${states} \
        --pam        ${pam} \
        --output     ${node}.html \
        --node       ${node}
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow {

    // Step 1
    clade_info = node_names_ch
        .combine(species_tree_val)
        | EXTRACT_CLADE

    // Step 2
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

    // Step 3
    pam_input = pvm_out
        .groupTuple(by: 0)
        .join(
            clade_info.map { node, in_sp, ign_sp, pruned -> tuple(node, in_sp) }
        )

    pam_out = pam_input | BUILD_PAM

    // Step 4
    asr_input = pam_out
        .join(
            clade_info.map { node, in_sp, ign_sp, pruned -> tuple(node, pruned) }
        )

    asr_out = asr_input | ANCESTRAL_RECON

    // Step 5 — pruned_tree flows through from ANCESTRAL_RECON output
    viz_input = asr_out
        .join(pam_out)
        .map { node, states, node_probs, pruned_tree, pam ->
            tuple(node, states, node_probs, pam, pruned_tree)
        }

    viz_input | VISUALIZE
}