nextflow.enable.dsl=2

// ── Parameters ────────────────────────────────────────────────────────────────
//
// Required:
//   --node_names   Comma-separated clade node name(s) matching named internal
//                  nodes in the species tree, e.g. "Bilateria" or
//                  "Bilateria,Metazoa".
//
// Optional (all have defaults matching the rest of the pipeline):
//   --ids               Path to ids.txt (list of HG IDs to process)
//   --SPECIES_TREE      Full species tree with named internal nodes
//   --REFSPECIES        Reference species prefix for POSSVM
//   --REFNAMES          Reference species gene-name CSV for POSSVM
//   --gene_trees_dir    Directory containing *.treefile outputs from step2
//   --min_presence      Minimum number of in-clade species for an OG to be
//                       kept in the PAM (default: 2)
//   --OUTDIR            Root output directory
//
// Example run (local, fast):
//   nextflow run step4.ancestry.nf \
//       -profile local \
//       --node_names Bilateria \
//       --ids ids.txt

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

species_tree_val = Channel.value( file(params.SPECIES_TREE) )

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

// ── Process 2 — POSSVM per clade × HG (ignoring out-of-clade species) ────────
//
// The --ignoretips flag tells POSSVM to disregard those species when building
// OGs, so the resulting ortholog groups are scoped to the target clade.
// We use the raw IQ-TREE / FastTree gene tree (not GeneRax) because
// GeneRax has already been pulled towards the species tree topology,
// which would bias the ancestral reconstruction.

process PVM_CLADE {

    tag "${node}/${id}"

    publishDir "${params.OUTDIR}/ancestry/${node}/possvm", mode: 'copy'

    cpus 1
    memory { 500.MB * task.attempt }
    time   { 10.min * task.attempt }

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(node), val(id), path(treefile), path(ignore_species)

    output:
    // POSSVM may produce 0 or more CSV files per run depending on tree size
    tuple val(node), path("${node}.${id}.*.ortholog_groups.csv"), optional: true

    script:
    """
    python ${projectDir}/phylogeny/main.py possvm \
        -t          ${treefile} \
        --refsps    ${params.REFSPECIES} \
        -r          ${params.REFNAMES} \
        --ignoretips ${ignore_species} \
        -o          ${node}.${id}.
    """
}

// ── Process 3 — Build presence/absence matrix across all HGs for a clade ─────

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

// ── Process 4 — Mk-model ancestral state reconstruction ──────────────────────
//
// Fits a two-rate (q01 = gain, q10 = loss) continuous-time Markov model
// by maximum likelihood for each OG column in the PAM, then computes the
// marginal probability P(present) at every internal node of the pruned
// species tree using the two-pass (Felsenstein pruning + peeling) algorithm.

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
          path("${node}.node_probs.tsv")

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

// ── Workflow ──────────────────────────────────────────────────────────────────

workflow {

    // Step 1: extract clade info for each named node
    clade_info = node_names_ch
        .combine(species_tree_val)
        | EXTRACT_CLADE
    // clade_info: tuple(node, in_species.txt, ignore_species.txt, pruned.tree)

    // Step 2: cross-product of clades × HGs, filtered to existing trees
    pvm_input = clade_info
        .map   { node, in_sp, ign_sp, pruned -> tuple(node, ign_sp) }
        .combine(hg_ids_ch)   // every (node, ignore_file) × every HG id
        .map { node, ign_sp, id ->
            def treefile = file("${params.gene_trees_dir}/${id}.treefile")
            tuple(node, id, treefile, ign_sp)
        }
        .filter { node, id, treefile, ign_sp -> treefile.exists() }

    pvm_out = pvm_input | PVM_CLADE
    // pvm_out: tuple(node, csv_file)  (one or more CSV per node/HG)

    // Step 3: collect all CSVs per node and build PAM
    // groupTuple aggregates all CSV path(s) emitted for the same node
    pam_input = pvm_out
        .groupTuple(by: 0)
        // Join with in_species file for each node
        .join(
            clade_info.map { node, in_sp, ign_sp, pruned -> tuple(node, in_sp) }
        )
    // pam_input: tuple(node, [csvs...], in_species.txt)

    pam_out = pam_input | BUILD_PAM
    // pam_out: tuple(node, pam.tsv)

    // Step 4: ancestral reconstruction
    asr_input = pam_out
        .join(
            clade_info.map { node, in_sp, ign_sp, pruned -> tuple(node, pruned) }
        )
    // asr_input: tuple(node, pam.tsv, pruned.tree)

    asr_input | ANCESTRAL_RECON
}
