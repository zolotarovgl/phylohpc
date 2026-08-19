nextflow.enable.dsl=2

// ── Canonical parameter names ─────────────────────────────────────────────────
// One canonical name per parameter, lower_snake_case. The uppercase twin
// (SPECIES_TREE, OUTDIR, REFNAMES, REFSPECIES, SUBS_MODEL, MAX_SPR, NCPU_GENERAX,
// TREE_METHOD, MAFFT_OPT, IQTREE2_MODEL) is a HARD ERROR if present at all -- see below.
//
// WHY THIS EXISTS: the uppercase SPECIES_TREE param (nextflow.config only) and
// params.species_tree (lowercase, the documented --species_tree flag) both existed and
// the GeneRax channel read the uppercase one. The documented flag therefore could not
// affect GeneRax, and on 2026-08-14 that handed it a multifurcating tree and killed all
// 475 tasks. Read every parameter through canonical() so a name can never fork again.
//
// WHY ERROR, NOT WARN (round 2, 2026-08-18): a warning here can only be reliable if
// `params.containsKey(upper)` is true exactly when a USER set the uppercase key --
// which requires nextflow.config to never default it either. But the profiles
// (`fast`/`precise`) legitimately need to override tree_method/mafft_opt/max_spr/
// ncpu_generax with real values, and doing that in lowercase makes
// `params.containsKey(lower)` true from profile load alone, defeating the OTHER branch
// of the shim the same way the top-level config defaults defeated this one (round 1).
// There is no config-only-lowercase-there / config-only-uppercase-nowhere split that
// keeps both directions honest AND lets profiles tune method params. So: nothing may
// ever define the uppercase keys again (not in nextflow.config, not in a profile) --
// `containsKey(upper)` is therefore true if and only if a user actually passed one, on
// the CLI or in a params-file, and a legacy key is refused outright rather than risking
// a silently-wrong resolution.
def canonical(String name, String upper, Object fallback = null) {
    if( params.containsKey(upper) )
        error "Parameter '${upper}' was renamed to '${name}'. Rename it in your params file / profile and re-run."
    if( params.containsKey(name) && params[name] != null )
        return params[name]
    return fallback
}

params.ids           = "${projectDir}/ids.txt"
params.resources_tsv = "${projectDir}/resources.tsv"
// MEASURED bug, fixed here while adding preflightChecks() below (config-layer only, no
// DAG/process change): nextflow.config's top-level params{} ALWAYS sets genefam_info =
// "genefam.csv" (a relative path, no projectDir prefix, and the file does not exist in
// the repo), so the `containsKey('genefam_info')` branch fired unconditionally and the
// exists-based fallbacks below it were dead code -- family_info silently resolved to a
// nonexistent path on every run with no --family_info override. Invisible until
// preflightChecks() actually checks existence. Now genefam_info is honoured only if it
// resolves to a real file.
params.family_info   = params.containsKey('family_info')
    ? params.family_info
    : (params.containsKey('genefam_info') && file(params.genefam_info.toString()).exists()
        ? params.genefam_info
        : (file("${projectDir}/genefam.csv").exists()
            ? "${projectDir}/genefam.csv"
            : "${projectDir}/data/gene_families_searchinfo.csv"))
params.species_tree_report = params.containsKey('species_tree_report')
    ? params.species_tree_report
    : "${projectDir}/data/species_tree.full.newick"
params.refnames       = canonical('refnames', 'REFNAMES', "${projectDir}/data/Mmus_gene_names.csv")
params.refsps         = canonical('refsps', 'REFSPECIES', 'Mmus')
params.run_generax    = params.containsKey('run_generax') ? params.run_generax : false
params.outdir         = canonical('outdir', 'OUTDIR', "${projectDir}/results")
params.mafft_opt      = canonical('mafft_opt', 'MAFFT_OPT', '')
params.tree_method    = canonical('tree_method', 'TREE_METHOD', 'iqtree2')
params.iqtree2_model  = canonical('iqtree2_model', 'IQTREE2_MODEL', 'TEST')
params.subs_model     = canonical('subs_model', 'SUBS_MODEL', 'LG')
params.max_spr        = canonical('max_spr', 'MAX_SPR', 3)
params.ncpu_generax   = canonical('ncpu_generax', 'NCPU_GENERAX', 2)


params.tag_prefix = ''

// Print what actually got resolved, BEFORE the DAG is built. The 2026-08-14 failure was
// invisible partly because nothing ever stated which file configured the run. Absorbs
// the old standalone "tree_method: ..." log.info line (still asserted by
// tests/test_params.sh's profile checks, via the "tree_method: <value>" line below).
def printResolvedConfig(String step) {
    log.info """
== resolved configuration (${step})
   params file   : ${workflow.configFiles.join(', ')}
   revision      : ${workflow.revision ?: 'local checkout'} (${workflow.commitId ?: 'no commit id'})
   outdir        : ${canonical('outdir','OUTDIR')}
   work dir      : ${workflow.workDir}
   species tree  : ${canonical('species_tree','SPECIES_TREE')}
   run_generax   : ${params.run_generax}
   tree_method: ${params.tree_method}
   family info   : ${params.family_info}
   ids restrict  : ${params.ids ?: '(all clusters)'}
"""
}
printResolvedConfig('step2')

// ── Pre-flight ─────────────────────────────────────────────────────────────────
// One named-cause check before the DAG is built -- the in-workflow twin of
// workflow/preflight.sh. The multifurcating-tree message is NOT hand-duplicated here:
// both this function and preflight.sh read the SAME template,
// workflow/messages/species_tree_multifurcating.txt, and substitute __TREE__/__N__ the
// same way, so the failure reads identically whichever way the pipeline was launched.
// (Round 2, 2026-08-19: an earlier version of this comment claimed the two messages
// were "kept word-for-word compatible" while in fact each had its own hand-written
// copy of the sentence, and they had already drifted -- caught in review, not by any
// test, which is why tests/test_preflight.sh now diffs the two outputs directly.)
// Exists because the 2026-08-14 run staged a multifurcating species tree, GeneRax said
// only "Cannot parse species tree file", the errorStrategy relabelled that as a family
// parsing error, and the run reported success with an empty output directory.
def species_tree_file = file(canonical('species_tree', 'SPECIES_TREE', "${projectDir}/data/species_tree.newick"))

def preflightChecks(treeFile) {
    ['infasta', 'family_info', 'refnames'].each { key ->
        def v = params[key]
        if( !v )
            error "'${key}' is not set"
        if( !file(v.toString()).exists() )
            error "'${key}' points at a missing file: ${v}"
    }

    if( params.run_generax ) {
        if( !treeFile.exists() )
            error "species tree not found: ${treeFile} (--species_tree)"
        // Delegate to the ONE tested implementation. An inline Groovy re-implementation
        // was tried first and disagreed with it (20 vs 4 multifurcations on the tree
        // that killed the 2026-08-14 run, 1 vs 0 on data/species_tree.newick) -- caught
        // by running both.
        def _cm = ["python3", "${projectDir}/workflow/count_multifurcations.py", treeFile.toString()].execute()
        _cm.waitFor()
        def bad = _cm.text.trim() as Integer
        if( bad > 0 ) {
            def msgTpl = file("${projectDir}/workflow/messages/species_tree_multifurcating.txt").text
            def msg = msgTpl.replace('__TREE__', treeFile.toString()).replace('__N__', bad.toString()).trim()
            error msg
        }
        log.info "species tree: ${treeFile} -- strictly binary"
    }
}
preflightChecks(species_tree_file)

// ── Failure diagnostics ───────────────────────────────────────────────────────
// Log the work dir + tail of .command.err whenever a task fails (any attempt).
def errReport(task) {
    def errf = task.workDir?.resolve('.command.err')
    def tail = (errf && errf.exists()) ? errf.text.readLines().takeRight(15).join('\n  ')
                                       : '(.command.err not found)'
    log.warn "✗ ${task.process} [${task.tag ?: task.index}] exit=${task.exitStatus} attempt=${task.attempt}\n" +
             "  workDir: ${task.workDir}\n  ${tail}"
}

workflow.onError {
    log.error "Pipeline stopped: ${workflow.errorMessage ?: 'see the failed task above'}"
}

workflow.onComplete {
    log.info "step2 ${workflow.success ? 'completed OK' : 'FAILED'} after ${workflow.duration}"
    if( !workflow.success )
        log.info "List failed tasks:  nextflow log ${workflow.runName} -f name,status,exit,workdir -F \"status=='FAILED'\""
}

// ── Environment preflight ─────────────────────────────────────────────────────
// Runs once before alignment/tree building; fails fast with one clear message if
// Python deps or core tools are missing (e.g. the conda env wasn't activated).
process PREFLIGHT {

    tag 'env_check'
    cpus 1
    memory '500.MB'
    time '5.min'
    cache false

    output:
    path 'preflight.ok'

    script:
    """
    export PYTHONNOUSERSITE=1
    missing=0
    python -c 'import Bio, ete3, yaml' 2>/dev/null || { echo "PREFLIGHT: missing Python deps (Bio/ete3/yaml)." >&2; missing=1; }
    for tool in mafft; do
        command -v \$tool >/dev/null 2>&1 || { echo "PREFLIGHT: '\$tool' not on PATH." >&2; missing=1; }
    done
    if [ "\$missing" -ne 0 ]; then
        echo "PREFLIGHT FAILED: environment not ready -- run 'mamba activate phylo' (and 'module load ...' on the HPC) before launching." >&2
        exit 1
    fi
    echo "PREFLIGHT OK"
    touch preflight.ok
    """
}

def countFastaSeqs(path) {
    def n = 0
    path.eachLine { line ->
        if( line.startsWith('>') )
            n++
    }
    return n
}

def res = [:]

if( params.resources_tsv && file(params.resources_tsv).exists() ) {

    def header = null

    file(params.resources_tsv).eachLine { line, n ->

        line = line.trim()
        if( !line || line.startsWith('#') ) return

        def cols = line.split('\t').collect{ it.trim() }

        if( n == 1 ) {
            header = cols
            return
        }

        def row = [:]
        header.eachWithIndex { h, i -> row[h] = cols[i] }

        def rid = row.id
        def m = [:]

        if( row.aln_mem ) m.aln_mem = row.aln_mem as nextflow.util.MemoryUnit
        if( row.aln_time ) m.aln_time = row.aln_time as nextflow.util.Duration

        if( row.phy_mem ) m.phy_mem = row.phy_mem as nextflow.util.MemoryUnit
        if( row.phy_time ) m.phy_time = row.phy_time as nextflow.util.Duration

        if( row.pvm_mem ) m.pvm_mem = row.pvm_mem as nextflow.util.MemoryUnit
        if( row.pvm_time ) m.pvm_time = row.pvm_time as nextflow.util.Duration

        if( row.gr_watcher_mem ) m.mem = row.gr_watcher_mem as nextflow.util.MemoryUnit
        if( row.gr_watcher_time ) m.time = row.gr_watcher_time as nextflow.util.Duration

        if( row.gr_mem ) m.mem = row.gr_mem as nextflow.util.MemoryUnit
        if( row.gr_time ) m.time = row.gr_time as nextflow.util.Duration

        res[rid] = m
    }
}

if( params.ids && file(params.ids).exists() ) {
    Channel
        .fromPath(params.ids)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { id -> tuple(id, file("${params.outdir}/clusters/${id}.fasta")) }
        .filter { id, fasta -> countFastaSeqs(fasta) >= 2 }
        .set { hg_fastas }
}
else {
    Channel
        .fromPath("${params.outdir}/clusters/*.fasta")
        .map { fasta -> tuple(fasta.baseName, fasta) }
        .filter { id, fasta -> countFastaSeqs(fasta) >= 2 }
        .set { hg_fastas }
}

process ALN {

    tag "${id}"

    publishDir "${params.outdir}/align", mode: 'copy'

    cpus 4

    errorStrategy = { errReport(task); task.attempt <= 10 ? 'retry' : 'ignore' }
    maxRetries 10
    maxErrors -1

    input:
    tuple val(id), path(fasta)

    output:
    tuple val(id), path("${id}.aln.fasta")

    script:
    def existing = file("${params.outdir}/align/${id}.aln.fasta")

    if (existing.exists()) {
        """
        ln -s ${existing} ${id}.aln.fasta
        """
    }
	    else {
	        """
	        export PYTHONNOUSERSITE=1
	        NSEQ=\$(grep -c '^>' ${fasta} || true)
	        if [[ "\$NSEQ" -lt 2 ]]; then
	            echo "Alignment input for ${id} contains fewer than 2 sequences; cannot build a trimmed alignment." >&2
	            exit 1
	        fi
	        if ! clipkit --help >/dev/null 2>&1; then
	            echo "clipkit is required for Step 2 trimming but is not runnable in the current environment." >&2
	            exit 1
	        fi
	        python ${projectDir}/phylogeny/main.py align -f ${fasta} -o ${id}.aln.fasta -c ${task.cpus} -m "${params.mafft_opt}"
	        python ${projectDir}/workflow/remove_gaponly.py ${id}.aln.fasta ${id}.aln.fasta_tmp
	        mv ${id}.aln.fasta_tmp ${id}.aln.fasta
	        """
	    }
}
process PHY {

    maxForks 50
    tag "${params.tag_prefix ? params.tag_prefix + '_' : ''}${id}"

    publishDir "${params.outdir}/gene_trees", mode: 'copy'

    cpus 4

    memory {
        def base = res[id]?.phy_mem ?: 300.MB
        return base * task.attempt
    }

    time {
        def base = res[id]?.phy_time ?: 30.min
        return base + (task.attempt - 1) * 6.h
    }

    errorStrategy = { errReport(task); task.attempt <= 10 ? 'retry' : 'ignore' }
    maxRetries 10
    maxErrors -1

    input:
    tuple val(id), path(aln)

    output:
    tuple val(id), path("${id}.treefile"), path(aln), path("${id}.log"), emit: trees
    path("${id}.ckp.gz"), optional: true, emit: ckp

    script:

    def existing     = file("${params.outdir}/gene_trees/${id}.treefile")
    def existing_ckp = file("${params.outdir}/gene_trees/${id}.ckp.gz")

    if (existing.exists()) {
        """
        echo "Using existing tree for ${id}"
        ln -sf ${existing} ${id}.treefile
        if [[ -e ${params.outdir}/gene_trees/${id}.log ]]; then
            ln -sf ${params.outdir}/gene_trees/${id}.log ${id}.log
        else
            printf "Using existing tree for %s\\nOriginal phylogeny log unavailable in %s\\n" "${id}" "${params.outdir}/gene_trees" > ${id}.log
        fi
        """
    }
    else if (params.tree_method == "iqtree2" && existing_ckp.exists()) {
        """
        echo "Resuming IQ-TREE2 from checkpoint for ${id}"
        cp ${existing_ckp} ${id}.ckp.gz
        python ${projectDir}/phylogeny/main.py phylogeny \
            -f ${aln} \
            --outprefix ${id} \
            -c ${task.cpus} \
            --method ${params.tree_method} \
            --iqtree2_model ${params.iqtree2_model} \
            > ${id}.phy_run.log 2>&1
        # iqtree2 writes ${id}.log natively (-pre); fasttree does not — fall back
        [ -s ${id}.log ] || cp ${id}.phy_run.log ${id}.log
        """
    }
	    else {
	        """
	        export PYTHONNOUSERSITE=1
	        python ${projectDir}/phylogeny/main.py phylogeny \
	            -f ${aln} \
	            --outprefix ${id} \
            -c ${task.cpus} \
            --method ${params.tree_method} \
            --iqtree2_model ${params.iqtree2_model} \
            > ${id}.phy_run.log 2>&1
        # iqtree2 writes ${id}.log natively (-pre); fasttree does not — fall back
        [ -s ${id}.log ] || cp ${id}.phy_run.log ${id}.log
        """
    }
}
process PVM {

	tag "${params.tag_prefix ? params.tag_prefix + '_' : ''}${id}"
	publishDir "${params.outdir}/possvm", mode: 'copy'

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
		errReport(task); return task.attempt <= 3 ? 'retry' : 'ignore'
	}

	maxRetries 3

	input:
	tuple val(id), path(tree), path(aln), path(refnames_file)

	output:
	tuple val(id),
		      path(tree),
		      path("${id}.*.ortholog_groups.newick"),
		      path("${id}.*.ortholog_groups.csv"),
		      path("${id}.*.pairs_orthologs.csv")

		script:
		"""
		export PYTHONNOUSERSITE=1
		python ${projectDir}/phylogeny/main.py possvm \
		    -t ${tree} \
		    --refsps ${params.refsps} \
	    -r ${refnames_file} \
	    -o ${id}.
	"""
}

process PVM_PREV {

	tag "${params.tag_prefix ? params.tag_prefix + '_' : ''}${id}"
	publishDir "${params.outdir}/possvm_prev", mode: 'copy'

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
		errReport(task); return task.attempt <= 3 ? 'retry' : 'ignore'
	}

	maxRetries 3
	maxErrors -1

	input:
	tuple val(id), path(tree), path(aln), path(refnames_file)

	output:
	tuple val(id),
	      path(tree),
	      path("${id}.*.ortholog_groups.newick"),
	      path("${id}.*.ortholog_groups.csv"),
	      path("${id}.*.pairs_orthologs.csv")

		script:
		"""
		export PYTHONNOUSERSITE=1
		python ${projectDir}/phylogeny/main.py possvm \
		    -t ${tree} \
		    --refsps ${params.refsps} \
	    -r ${refnames_file} \
	    -o ${id}.
	"""
}

// -----------------------------
// Step-2 HTML report
// -----------------------------
process REPORT {

    publishDir "${params.outdir}", mode: 'copy'

    cpus   1
    memory 2.GB
    time   30.min

    input:
    path(newicks)   // collected PVM newick outputs — used as a completion barrier

    output:
    path("report_step2.html")

	    script:
        def reportRefArgs = []
        if (params.refnames) reportRefArgs << "--refnames ${params.refnames}"
        if (params.refsps)   reportRefArgs << "--refsps ${params.refsps}"
        def refArgs = reportRefArgs.join(' ')
	    """
		    export PYTHONNOUSERSITE=1
		    python ${projectDir}/workflow/report_step2.py \
		        --results_dir     ${params.outdir} \
		        --family_info     ${params.family_info} \
		        --species_tree    ${params.species_tree_report} \
		        --species_info    ${projectDir}/data/species_info.tsv \
	            ${refArgs} \
	        --output          report_step2.html
    """
}

// -----------------------------
// GeneRax process
// -----------------------------
process GR_watcher {

    tag "${id}"

    publishDir "${params.outdir}/generax", mode: 'copy'

    cpus params.ncpu_generax
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
        errReport(task)
        def max_attempts = 5
        if( task.exitStatus == 10 ) {
            log.warn "GeneRax | ${id} | Exit 10 | Family parsing error — ignored"
            return 'ignore'
        }
        else if( task.attempt <= max_attempts ) {
            if( task.exitStatus == 137 ) {
                log.warn "GeneRax | ${id} | Exit 137 | Likely OOM — retrying (attempt ${task.attempt}/${max_attempts})"
            }
            else {
                log.warn "GeneRax | ${id} | Exit ${task.exitStatus} | Retrying (attempt ${task.attempt}/${max_attempts})"
            }
            return 'retry'
        }
        else {
            log.warn "GeneRax | ${id} | Exit ${task.exitStatus} | Retries exhausted — ignored"
            return 'ignore'
        }
    }

    maxRetries 5
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

    def existing = file("${params.outdir}/generax/${id}.generax.tree")

    if (existing.exists()) {
        """
        echo "Using existing GeneRax result for ${id}"

        ln -sf ${existing} ${id}.generax.tree

        # dummy files so outputs exist
        touch ${id}.generax.log
        touch ${id}.progress.tree
        """
    }
    else {
        """
	        set -euo pipefail
	        export PYTHONNOUSERSITE=1

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
                sleep 10
            done
        }

        python ${projectDir}/phylogeny/main.py generax \
            --name ${id} \
            --alignment ${aln} \
            --gene_tree ${tree} \
            --species_tree ${species_tree} \
            --output_dir ${id}_generax \
            --subs_model ${params.subs_model} \
            --max_spr ${params.max_spr} \
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

//workflow {
//	hg_fastas|ALN|PHY|map { id, tree, aln -> tuple(id, tree, aln, refnames_file) } | PVM
//}

// ---------------------------------------------------------------------------------
// SPECIES TREE -- one parameter, validated before anything runs.
//
// BUG FIXED 2026-08-18. This channel used to read the uppercase SPECIES_TREE param,
// which is set ONLY in nextflow.config and defaults to data/species_tree.newick. The
// documented switch, `--species_tree` (lowercase, resolved at the top of this file and
// used by the ancestral step), therefore had NO EFFECT ON GENERAX. On 2026-08-14 that
// silently handed GeneRax the UNRESOLVED tree and all 475 GR_watcher tasks died with
//     [Error] Cannot parse species tree file (species_tree.newick)
// which the errorStrategy relabelled "Family parsing error", sending the diagnosis at
// the family files for four days. GeneRax requires a STRICTLY BINARY species tree.
//
// Precedence now: --species_tree (or a params-file / YAML key) wins; SPECIES_TREE is
// honoured only as a legacy fallback via canonical(), and the tree is checked by
// preflightChecks() (called right after printResolvedConfig, near the top of this
// file) before the DAG is built. species_tree_file itself is also defined up there so
// preflightChecks() can see it; this section just wires it into a channel.
species_tree_ch = Channel.value( species_tree_file )
refnames_ch     = Channel.value( file(params.refnames) )

workflow {

    // Fail fast if the environment is not ready.
    def ready = PREFLIGHT()
    def hg_in = hg_fastas.combine(ready).map { id, fasta, ok -> tuple(id, fasta) }

    PHY(hg_in | ALN)
    phy_out = PHY.out.trees

    if (params.run_generax) {

        // ---------- PVM on original trees ----------
        pvm_prev_out = phy_out
            .map { id, tree, aln, log -> tuple(id, tree, aln) }
            .combine(refnames_ch)
            | PVM_PREV


        // ---------- GeneRax ----------
        gr_input = phy_out
            .map { id, tree, aln, log ->
                tuple(id, aln, tree)
            }
            .combine(species_tree_ch)

        gr_out = gr_input | GR_watcher


        // ---------- PVM on GeneRax trees ----------
        pvm_out = gr_out
            .map { id, generax_tree, log, progress, aln ->
                tuple(id, generax_tree, aln)
            }
            .combine(refnames_ch)
            | PVM

        // ---------- Report (wait for both PVM and PVM_PREV) ----------
        pvm_out.map { id, tree, nwk, csv, pairs -> nwk }
            .mix(pvm_prev_out.map { id, tree, nwk, csv, pairs -> nwk })
            .collect()
            | REPORT

    }
    else {

        pvm_out = phy_out
            .map { id, tree, aln, log ->
                tuple(id, tree, aln)
            }
            .combine(refnames_ch)
            | PVM

        // ---------- Report ----------
        pvm_out
            .map { id, tree, nwk, csv, pairs -> nwk }
            .collect()
            | REPORT
    }
}
