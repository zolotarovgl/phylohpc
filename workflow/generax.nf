nextflow.enable.dsl=2

params.resources_tsv = null
params.ALIGN_DIR     = "${projectDir}/../results/align"
params.TREE_DIR      = "${projectDir}/../results/gene_trees"
params.OUTDIR        = "${projectDir}/../results"

// ── species_tree resolution ─────────────────────────────────────────────────
// CRITICAL 4 (2026-08-19 whole-branch review): this file is a second, independent
// GeneRax implementation with its own workflow{} block (a documented live entry point,
// docs/manual.md:416-436) and it was still reading params.SPECIES_TREE directly --
// exactly the parameter that caused the 2026-08-14 incident in step2.nf, and since
// nextflow.config no longer defaults ANY of the uppercase keys (round 1 of this
// refactor), that channel now resolves a NULL tree here with no local default either.
// Same resolution + guard as step2.nf: canonical() makes the uppercase key a hard
// error rather than a silent null, and the binary-tree check runs before the DAG is
// built. Reuses the SAME count_multifurcations.py and message template step2.nf uses
// -- not a third copy of either.
def canonical(String name, String upper, Object fallback = null) {
    if( params.containsKey(upper) )
        error "Parameter '${upper}' was renamed to '${name}'. Rename it in your params file / profile and re-run."
    if( params.containsKey(name) && params[name] != null )
        return params[name]
    return fallback
}
def _species_tree_val = canonical('species_tree', 'SPECIES_TREE')
if( !_species_tree_val )
    error "'species_tree' is not set -- pass --species_tree <path> (no default in this file, unlike step2.nf)"
def species_tree_file = file(_species_tree_val)
if( !species_tree_file.exists() )
    error "species tree not found: ${species_tree_file} -- pass --species_tree <path>"
def _cm = ["python3", "${projectDir}/count_multifurcations.py", species_tree_file.toString()].execute()
_cm.waitFor()
def _bad = _cm.text.trim() as Integer
if( _bad > 0 ) {
    def _msgTpl = file("${projectDir}/messages/species_tree_multifurcating.txt").text
    error _msgTpl.replace('__TREE__', species_tree_file.toString()).replace('__N__', _bad.toString()).trim()
}
// REFSPECIES/REFNAMES are left uppercase in this file deliberately (out of scope for
// this fix -- see the CRITICAL 4 comment above) but must NOT silently resolve null now
// that nextflow.config no longer defaults them either: fail loudly rather than let
// PVM/GR_watcher fail downstream with an opaque error far from the real cause.
if( !params.REFSPECIES )
    error "'REFSPECIES' is not set -- pass --REFSPECIES <code>, e.g. Mmus"
if( !params.REFNAMES || !file(params.REFNAMES.toString()).exists() )
    error "'REFNAMES' points at a missing file: ${params.REFNAMES} -- pass --REFNAMES <path>"
log.info "species tree: ${species_tree_file} -- strictly binary"

species_tree_ch = Channel.fromPath(species_tree_file)
// -----------------------------
// Load per-family resources
// -----------------------------
def res = [:]
if( params.resources_tsv && file(params.resources_tsv).exists() ) {
    file(params.resources_tsv).eachLine { line, n ->
        line = line.trim()
        if( !line || line.startsWith('#') )
            return
        def cols = line.split('\t', -1).collect { it.trim() }
        if( n == 1 && cols[0].toLowerCase() == 'id' )
            return
        def rid = cols[0]
        def m = [:]
        if( cols.size() > 1 && cols[1] ) m.mem = cols[1] as nextflow.util.MemoryUnit
        if( cols.size() > 2 && cols[2] ) m.time = cols[2] as nextflow.util.Duration
        res[rid] = m
    }
}

// -----------------------------
// Channel of IDs
// -----------------------------
Channel
    .fromPath(params.ids)
    .splitText()
    .map { it.trim() }
    .filter { it }
    .map { id ->
        tuple(
            id,
            file("${params.ALIGN_DIR}/${id}.aln.fasta"),
            file("${params.TREE_DIR}/${id}.treefile")
        )
    }
    .set { hg_inputs }

// -----------------------------
// GeneRax process
// -----------------------------
process GR_watcher {

    tag "${id}"

    publishDir "${params.OUTDIR}/generax", mode: 'copy'

    cpus params.NCPU_GENERAX
    maxForks 30

    memory {
        def base = res[id]?.mem ?: 500.MB
        return base * Math.pow(2, task.attempt-1)
    }

    time {
        def base = res[id]?.time ?: 30.min
        def scaled = base * Math.pow(2, task.attempt-1)
        return scaled > 24.h ? 24.h : scaled
    }

    errorStrategy {
        if( task.exitStatus == 10 ) {
            log.warn "GeneRax | ${id} | Exit 10 | Family parsing error — ignored"
            return 'ignore'
        }
        else if( task.exitStatus == 137 ) {
            log.warn "GeneRax | ${id} | Exit 137 | Likely OOM — retrying (attempt ${task.attempt})"
            return 'retry'
        }
        else {
            log.warn "GeneRax | ${id} | Exit ${task.exitStatus} | Retrying"
            return 'retry'
        }
    }

    maxRetries 2
    maxErrors -1

    input:
    tuple val(id), path(aln), path(tree), path(species_tree)

    output:
    tuple val(id),
          path("${id}.generax.tree"),
          path("${id}.generax.log"),
          path("${id}.progress.tree"),
          path(aln)

	script:

	def existing = file("${params.OUTDIR}/generax/${id}.generax.tree")

	if (existing.exists()) {
		"""
		echo "Using existing GeneRax result for ${id}"
		ln -sf ${existing} ${id}.generax.tree
		touch ${id}.generax.log
		touch ${id}.progress.tree
		"""
	}
	else {
		"""
		set -euo pipefail

		export OMP_NUM_THREADS=${task.cpus}
		export OPENBLAS_NUM_THREADS=${task.cpus}
		export MKL_NUM_THREADS=${task.cpus}
		export NUMEXPR_NUM_THREADS=${task.cpus}

		touch ${id}.progress.tree

		progress_watcher() {
			while kill -0 \$MAIN_PID 2>/dev/null; do
				if [[ -f ${id}_generax/results/${id}/geneTree.newick ]]; then
					cp ${id}_generax/results/${id}/geneTree.newick \
					${id}.progress.tree 2>/dev/null || true
				fi
				sleep 60
			done
		}

		python ${projectDir}/../phylogeny/main.py generax \
			--name ${id} \
			--alignment ${aln} \
			--gene_tree ${tree} \
			--species_tree ${species_tree} \
			--output_dir ${id}_generax \
			--subs_model ${params.SUBS_MODEL} \
			--max_spr ${params.MAX_SPR} \
			--cpus ${task.cpus} \
			--logfile ${id}.generax.log \
			--outfile ${id}.generax.tree &

		MAIN_PID=\$!

		progress_watcher &
		WATCH_PID=\$!

		cleanup() {
			kill \$WATCH_PID 2>/dev/null || true
			wait \$WATCH_PID 2>/dev/null || true
		}
		trap cleanup EXIT INT TERM

		wait \$MAIN_PID
		EXIT_CODE=\$?

		echo "GeneRax exit code: \$EXIT_CODE"

		if [[ -f ${id}_generax/results/${id}/geneTree.newick ]]; then
			cp ${id}_generax/results/${id}/geneTree.newick ${id}.progress.tree || true
		fi

		exit \$EXIT_CODE
		"""
	}
}

process GR {

    tag "${id}"

    publishDir "${params.OUTDIR}/generax", mode: 'copy'

    cpus params.NCPU_GENERAX
	maxForks 10

    memory {
        def base = res[id]?.mem ?: 500.MB
        return base * Math.pow(2, task.attempt-1)
    }

    time {
        def base = res[id]?.time ?: 30.min
        def scaled = base * Math.pow(2, task.attempt-1)
        return scaled > 24.h ? 24.h : scaled
    }

    errorStrategy {
        if( task.exitStatus == 10 ) {
            log.warn "GeneRax | ${id} | Exit 10 | Family parsing error — ignored"
            return 'ignore'
        }
        else if( task.exitStatus == 137 ) {
            log.warn "GeneRax | ${id} | Exit 137 | Likely OOM — retrying with more memory (attempt ${task.attempt})"
            return 'retry'
        }
        else {
            log.warn "GeneRax | ${id} | Exit ${task.exitStatus} | Retrying"
            return 'retry'
        }
    }

    maxRetries 2
    maxErrors -1

    input:
    tuple val(id), path(aln), path(tree), path(species_tree)

	output:
	tuple val(id), path("${id}.generax.tree"), path("${id}.generax.log"), path(aln)

	script:
	"""
	set -euo pipefail
	touch ${id}.progress.tree
	python ${projectDir}/../phylogeny/main.py generax \
		--name ${id} \
		--alignment ${aln} \
		--gene_tree ${tree} \
		--species_tree ${species_tree} \
		--output_dir ${id}_generax \
		--subs_model ${params.SUBS_MODEL} \
		--max_spr ${params.MAX_SPR} \
		--cpus ${task.cpus} \
		--logfile ${id}.generax.log \
		--outfile ${id}.generax.tree

	EXIT_CODE=\$?

	echo "GeneRax exit code: \$EXIT_CODE"

	exit \$EXIT_CODE
	"""
}
process PVM {

    tag "${id}"

    publishDir "${params.OUTDIR}/generax", mode: 'copy'

    cpus 1

    memory {
        def base = res[id]?.pvm_mem ?: 500.MB
        return base + (task.attempt - 1) * 500.MB
    }

    time {
        def base = res[id]?.pvm_time ?: 5.min
        return base * task.attempt
    }

    errorStrategy {
        return task.attempt <= 3 ? 'retry' : 'ignore'
    }

	maxRetries 3

	input:
	tuple val(id), path(tree), path(log), path(aln), path(refnames_file)

    output:
    tuple val(id),
          path("${id}.generax.tree.ortholog_groups.newick"),
          path("${id}.generax.tree.ortholog_groups.csv"),
          path("${id}.generax.tree.pairs_orthologs.csv")

    script:
    """
    python ${projectDir}/../phylogeny/main.py possvm \
        -t ${tree} \
        --refsps ${params.REFSPECIES} \
        -r ${refnames_file} \
        -o ${id}.
    """
}
// -----------------------------
// Workflow
// -----------------------------
refnames_ch = Channel.value( file(params.REFNAMES) )

workflow {

	gr_out = hg_inputs
		.combine(species_tree_ch)
		| GR_watcher

	pvm_in = gr_out
		.map { id, tree, log, progress, aln ->
			tuple(id, tree, log, aln)
		}
		.combine(refnames_ch)

    pvm_in | PVM
}