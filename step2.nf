nextflow.enable.dsl=2
params.ids           = "${projectDir}/ids.txt"
params.resources_tsv = "${projectDir}/resources.tsv"
params.family_info   = params.containsKey('family_info')
    ? params.family_info
    : (params.containsKey('genefam_info')
        ? params.genefam_info
        : (file("${projectDir}/genefam.csv").exists()
            ? "${projectDir}/genefam.csv"
            : "${projectDir}/data/gene_families_searchinfo.csv"))
params.species_tree  = params.containsKey('species_tree')
    ? params.species_tree
    : "${projectDir}/data/species_tree.full.newick"
params.refnames      = params.containsKey('refnames')
    ? params.refnames
    : (params.containsKey('REFNAMES') ? params.REFNAMES : null)
params.refsps        = params.containsKey('refsps')
    ? params.refsps
    : (params.containsKey('REFSPECIES') ? params.REFSPECIES : null)
params.run_generax   = params.containsKey('run_generax') ? params.run_generax : false
// Compute SH-aLRT support on the reconciled topology before POSSVM (see process GXSUP).
// Set --gxsup false to hand POSSVM the raw GeneRax tree, whose "support" is a constant 1.0.
params.gxsup         = params.containsKey('gxsup') ? params.gxsup : true
params.OUTDIR        = params.containsKey('OUTDIR')
    ? params.OUTDIR
    : (params.containsKey('outdir') ? params.outdir : "${projectDir}/results")
//params.MAFFT_OPT = "--maxiterate 1000 --localpair"


params.tag_prefix = ''

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
        .map { id -> tuple(id, file("${params.OUTDIR}/clusters/${id}.fasta")) }
        // nseq is carried so ALN/PHY can size resources by family, in the process body or
        // from a -c config (conf/resources_bvw.config). countFastaSeqs is called once here
        // rather than once per directive evaluation.
        .map { id, fasta -> tuple(id, fasta, countFastaSeqs(fasta)) }
        .filter { id, fasta, nseq -> nseq >= 2 }
        .set { hg_fastas }
}
else {
    Channel
        .fromPath("${params.OUTDIR}/clusters/*.fasta")
        .map { fasta -> tuple(fasta.baseName, fasta) }
        .map { id, fasta -> tuple(id, fasta, countFastaSeqs(fasta)) }
        .filter { id, fasta, nseq -> nseq >= 2 }
        .set { hg_fastas }
}

process ALN {

    tag "${id}"

    publishDir "${params.OUTDIR}/align", mode: 'copy'

    cpus 4

    // ALN had NO memory directive, so it took the profile default of 4 GB and every one of
    // its 10 retries asked for the same 4 GB -- ten identical OOM kills. Measured 2026-08-25:
    // mafft L-INS-i on adh.fibronectin.HG9 (153 seqs) died as
    //   mafft: line 2842: ... tbfast ... Killed
    // i.e. SIGKILL, not a mafft error. Ladder the memory the way PHY and GR_watcher already
    // do. Attempt 1 stays at 4 GB so nothing that works today changes; the floor is 4 GB
    // because resources.tsv's aln_mem is often SMALLER (800 MB - 2.1 GB) and would make
    // the first attempt worse.
    memory {
        def base = res[id]?.aln_mem
        base = (base && base > 4.GB) ? base : 4.GB
        def scaled = base * Math.pow(2, task.attempt-1)
        return scaled > 64.GB ? 64.GB : scaled
    }

    errorStrategy = { task.attempt <= 10 ? 'retry' : 'ignore' }
    maxRetries 10
    maxErrors -1

    input:
    tuple val(id), path(fasta), val(nseq)

    output:
    tuple val(id), path("${id}.aln.fasta"), val(nseq)

    script:
    def existing = file("${params.OUTDIR}/align/${id}.aln.fasta")

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
	        python ${projectDir}/phylogeny/main.py align -f ${fasta} -o ${id}.aln.fasta -c ${task.cpus} -m "${params.MAFFT_OPT}"
	        python ${projectDir}/workflow/remove_gaponly.py ${id}.aln.fasta ${id}.aln.fasta_tmp
	        mv ${id}.aln.fasta_tmp ${id}.aln.fasta
	        """
	    }
}
process PHY {

    maxForks 50
    tag "${params.tag_prefix ? params.tag_prefix + '_' : ''}${id}"

    publishDir "${params.OUTDIR}/gene_trees", mode: 'copy'

    cpus 4

    memory {
        def base = res[id]?.phy_mem ?: 300.MB
        return base * task.attempt
    }

    time {
        def base = res[id]?.phy_time ?: 30.min
        return base + (task.attempt - 1) * 6.h
    }

    // --signal=B:USR2@180 asks SLURM to warn the batch shell 3 min before the wall-clock
    // kill. That window is what lets the checkpoint stash in the script body run; without
    // it the task is killed outright and the IQ-TREE2 checkpoint dies with the work dir.
    // nseq is carried on the channel so the qos can follow family size.
    clusterOptions {
        def qos = (nseq > 500 || task.attempt > 1) ? '--qos=long' : '--qos=normal'
        return "${qos} --signal=B:USR2@180"
    }

    // maxRetries must EXCEED the closure's threshold: if the closure still returns 'retry'
    // on the attempt where maxRetries is reached, the 'ignore' branch is never evaluated
    // and the whole run terminates instead of skipping the family. (BvW)
    errorStrategy { task.attempt <= 10 ? 'retry' : 'ignore' }
    maxRetries 11
    maxErrors -1

    input:
    tuple val(id), path(aln), val(nseq)

    output:
    tuple val(id), path("${id}.treefile"), path(aln), path("${id}.log"), emit: trees
    path("${id}.ckp.gz"), optional: true, emit: ckp

    script:

    def existing = file("${params.OUTDIR}/gene_trees/${id}.treefile")
    def ckp_dir  = "${params.OUTDIR}/phy_ckp"

    if (existing.exists()) {
        """
        echo "Using existing tree for ${id}"
        ln -sf ${existing} ${id}.treefile
        if [[ -e ${params.OUTDIR}/gene_trees/${id}.log ]]; then
            ln -sf ${params.OUTDIR}/gene_trees/${id}.log ${id}.log
        else
            printf "Using existing tree for %s\\nOriginal phylogeny log unavailable in %s\\n" "${id}" "${params.OUTDIR}/gene_trees" > ${id}.log
        fi
        """
    }
    else {
        """
        export PYTHONNOUSERSITE=1

        CKP_DIR="${ckp_dir}"
        CKP_STASH="\$CKP_DIR/${id}.ckp.gz"
        mkdir -p "\$CKP_DIR"

        # A task that times out or OOMs never reaches publishDir, so its IQ-TREE2
        # checkpoint dies with the work directory and the next attempt redoes
        # ModelFinder from scratch. Keep it in a stable per-family location instead,
        # so retries resume where the previous attempt was killed. This is what the
        # giant TF families needed: exit 140 (128+SIGUSR2) on a 24 h limit currently
        # throws away every hour of work done so far.
        for src in "\$CKP_STASH" "${params.OUTDIR}/gene_trees/${id}.ckp.gz"; do
            if [[ -s "\$src" ]] && gzip -t "\$src" 2>/dev/null; then
                cp "\$src" ${id}.ckp.gz
                echo "Resuming ${id} from checkpoint \$src" >> ${id}.log
                break
            fi
        done

        # gzip -t before installing: SLURM's warning signal can arrive while IQ-TREE2
        # is mid-write, and a truncated checkpoint makes the next attempt fail instantly
        # rather than resume. A bad copy is discarded and the previous good stash kept.
        stash_ckp() {
            if [[ -s ${id}.ckp.gz ]] && gzip -t ${id}.ckp.gz 2>/dev/null; then
                cp -f ${id}.ckp.gz "\$CKP_STASH.tmp" && mv -f "\$CKP_STASH.tmp" "\$CKP_STASH"
            fi
            return 0
        }
        trap stash_ckp EXIT
        trap 'stash_ckp; exit 140' USR2
        trap 'stash_ckp; exit 143' TERM

        # Background + wait so the shell can service USR2 promptly; a foreground child
        # would defer the trap until after it exits.
        set +e
        python ${projectDir}/phylogeny/main.py phylogeny -f ${aln} --outprefix ${id} -c ${task.cpus} --method ${params.TREE_METHOD} --iqtree2_model ${params.IQTREE2_MODEL} >> ${id}.log 2>&1 &
        PY_PID=\$!
        wait "\$PY_PID"
        RC=\$?
        set -e

        trap - EXIT USR2 TERM
        if [[ "\$RC" -eq 0 ]]; then
            # Succeeded: drop the stash so a later rerun on changed input can never
            # resume from a stale checkpoint.
            rm -f "\$CKP_STASH"
        else
            stash_ckp
        fi
        exit \$RC
        """
    }
}

process PVM {

	tag "${params.tag_prefix ? params.tag_prefix + '_' : ''}${id}"
	publishDir "${params.OUTDIR}/possvm", mode: 'copy'

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
		# --skiproot: this leg reads a GeneRax tree, which the reconciliation ALREADY rooted.
		# POSSVM roots by default (iterative midpoint), and re-rooting a reconciled tree
		# collapsed 38 of 146 TF families to one orthogroup spanning the whole tree
		# (tfs.Forkhead.HG2: 1 group rooted vs 12 with the flag). PVM_PREV deliberately does
		# NOT get this -- its IQ-TREE trees are genuinely unrooted and need the rooting.
		python ${projectDir}/phylogeny/main.py possvm \
		    -t ${tree} \
		    --refsps ${params.REFSPECIES} \
		    --skiproot \
	    -r ${refnames_file} \
	    -o ${id}.

		# POSSVM reports a PLACEHOLDER as if it were support: a singleton group has no ancestor
		# to test and a root-spanning group gets the root's own default, and ete3 returns 1.0 for
		# both. Measured over 400 families: 15.4 % singleton, 4.8 % root-MRCA, and 19 groups whose
		# 1.0 is REAL -- which is why this recomputes the MRCA instead of matching on "1.0".
		# Drops singleton rows; sets support to -1 where the MRCA is the root.
		for f in ${id}.*.ortholog_groups.csv; do
		    if [ -s "\$f" ]; then python ${projectDir}/workflow/possvm_postprocess.py "\$f"; fi
		done
	"""
}

process PVM_PREV {

	tag "${params.tag_prefix ? params.tag_prefix + '_' : ''}${id}"
	publishDir "${params.OUTDIR}/possvm_prev", mode: 'copy'

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
		    --refsps ${params.REFSPECIES} \
	    -r ${refnames_file} \
	    -o ${id}.

		# POSSVM reports a PLACEHOLDER as if it were support: a singleton group has no ancestor
		# to test and a root-spanning group gets the root's own default, and ete3 returns 1.0 for
		# both. Measured over 400 families: 15.4 % singleton, 4.8 % root-MRCA, and 19 groups whose
		# 1.0 is REAL -- which is why this recomputes the MRCA instead of matching on "1.0".
		# Drops singleton rows; sets support to -1 where the MRCA is the root.
		for f in ${id}.*.ortholog_groups.csv; do
		    if [ -s "\$f" ]; then python ${projectDir}/workflow/possvm_postprocess.py "\$f"; fi
		done
	"""
}

// -----------------------------
// Step-2 HTML report
// -----------------------------
process REPORT {

    publishDir "${params.OUTDIR}", mode: 'copy'

    cpus   1
    memory 2.GB
    time   30.min

    input:
    path(newicks)   // collected PVM newick outputs — used as a completion barrier

    output:
    path("report_step2.html")

	    script:
        def reportRefArgs = []
        if (params.refnames) reportRefArgs << "--refnames ${file(params.refnames)}"
        if (params.refsps)   reportRefArgs << "--refsps ${params.refsps}"
        def refArgs = reportRefArgs.join(' ')
	    """
		    export PYTHONNOUSERSITE=1
		    python ${projectDir}/workflow/report_step2.py \
		        --results_dir     ${params.OUTDIR} \
		        --family_info     ${file(params.family_info)} \
		        --species_tree    ${file(params.species_tree)} \
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

    // --signal=B:USR2@300 gives the batch shell 5 min of warning before the wall-clock
    // kill, which is what lets the trap stash GeneRax's latest progress tree. GeneRax has
    // no checkpoint of its own, so without this a killed attempt discards everything.
    clusterOptions { '--signal=B:USR2@300' }

    // one above the errorStrategy closure's threshold, or 'ignore' is unreachable
    maxRetries 6
    maxErrors -1

    input:
    tuple val(id), path(aln), path(tree), path(species_tree)

    output:
    tuple val(id),
          path("${id}.generax.tree"),
          path("${id}.generax.log"),
          path("${id}.progress.tree"),
          path(aln), emit: trees
    // GeneRax writes a reconciliations/ directory next to the gene tree holding
    // _orthogroups.txt, _orthogroups_all.txt, _events.newick, _reconciliated.nhx/.xml,
    // _transfers.txt and the event counts. None of it was declared, so all of it died with
    // the work directory. Emitted as its own channel, and OPTIONAL because the
    // "existing result" branch below re-uses a published tree and never re-runs GeneRax.
    path("${id}.reconciliations"), optional: true, emit: recs

    script:

    def existing = file("${params.OUTDIR}/generax/${id}.generax.tree")

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

        GR_STASH_DIR="${params.OUTDIR}/generax_ckp"
        GR_STASH="\$GR_STASH_DIR/${id}.progress.tree"
        mkdir -p "\$GR_STASH_DIR"

        touch ${id}.progress.tree

        # GeneRax has NO checkpoint of its own: a killed run loses everything and the retry
        # restarts from the original IQ-TREE topology. Its SPR search improves the tree
        # iteratively though, and the watcher below already snapshots that. So stash the
        # snapshot somewhere stable and start the next attempt from it. Measured on the TF
        # run: the giants were being killed at 16-18 h having discarded every previous hour.
        #
        # Validation before use: same tip count as the alignment and a trailing ';'. The
        # snapshot is a plain cp of a file GeneRax may be mid-write on, so a truncated copy
        # is expected occasionally and must not be fed back in.
        GENE_TREE="${tree}"
        if [[ -s "\$GR_STASH" ]]; then
            want=\$(grep -c '^>' ${aln})
            got=\$(tr -cd ',' < "\$GR_STASH" | wc -c)
            got=\$(( got + 1 ))
            if [[ "\$got" -eq "\$want" ]] && tail -c 2 "\$GR_STASH" | grep -q ';'; then
                cp "\$GR_STASH" ${id}.resume.tree
                GENE_TREE="${id}.resume.tree"
                echo "Resuming ${id} from stashed progress tree (\$got tips)" >&2
            else
                echo "Discarding stashed progress tree for ${id}: \$got/\$want tips, or no terminal ';'" >&2
            fi
        fi

        stash_progress() {
            if [[ -s ${id}.progress.tree ]] && tail -c 2 ${id}.progress.tree | grep -q ';'; then
                cp -f ${id}.progress.tree "\$GR_STASH.tmp" && mv -f "\$GR_STASH.tmp" "\$GR_STASH"
            fi
            return 0
        }

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
            --gene_tree \$GENE_TREE \
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
            stash_progress
        }
        trap cleanup EXIT INT TERM
        trap 'cleanup; exit 140' USR2

        wait \$MAIN_PID
        EXIT_CODE=\$?

        echo "GeneRax exit code: \$EXIT_CODE"

        if [[ -f ${id}_generax/results/${id}/geneTree.newick ]]; then
            cp ${id}_generax/results/${id}/geneTree.newick ${id}.progress.tree || true
        fi

        # Keep GeneRax's reconciliation products (orthogroups, events, transfers). They are
        # written whether or not we ask, and were being discarded with the work directory.
        if [[ -d ${id}_generax/reconciliations ]]; then
            cp -r ${id}_generax/reconciliations ${id}.reconciliations || true
        fi

        if [[ "\$EXIT_CODE" -eq 0 ]]; then
            rm -f "\$GR_STASH"
        else
            stash_progress
        fi

        exit \$EXIT_CODE
        """
    }
}

// -----------------------------
// Branch support for the reconciled trees
// -----------------------------
// GeneRax computes no branch support. Its geneTree.newick carries a constant placeholder
// (one empty label, one "0", N "1"s -- identical across families of every size), which
// POSSVM reads as support 1.0, so --min_support_transfer was being compared against a
// constant and the label-transfer step was silently disabled on every reconciled family.
//
// The original UFBoot values cannot be reused: PHY runs `-bb 1000` without --wbtl, so no
// replicate trees were written, and GeneRax moves a lot of the topology anyway (measured
// over 15 random families: median 68.7% of splits shared, range 37.5-86.7%).
//
// So compute support FOR the reconciled topology instead: fix the GeneRax tree, re-optimise
// branch lengths on the same alignment under the family's own best-fit model, and run
// SH-aLRT. Verified 2026-08-27 on 3 families -- with --keep-ident the output topology is
// identical to GeneRax's (RF 0/88, 0/88, 0/204); without it, IQ-TREE drops and reinserts
// duplicate sequences and moves up to 5 of 44 splits. ~7 s per family at 2 threads.
// Median SH-aLRT 62.8-78.3, and 52-69% of nodes clear 50, so the existing threshold works.
process GXSUP {

    tag "${id}"

    // pattern: the alignment is in the output tuple because PVM needs it downstream, but it
    // must not be published -- without this every alignment is duplicated into
    // generax_support/ (~1500 files, ~200 MB of copies). The .iqtree report is kept because
    // it is the only record of the model and likelihood the support was computed under.
    publishDir "${params.OUTDIR}/generax_support", mode: 'copy',
               pattern: "*.{generax.support.tree,gxsup.iqtree}"

    cpus 2
    memory { 2.GB * task.attempt }
    time   { 30.min * task.attempt }

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(id), path(gtree), path(aln), path(phylog)

    output:
    tuple val(id), path("${id}.generax.support.tree"), path(aln), emit: trees
    path("${id}.gxsup.iqtree"), emit: report

    script:
    """
    set -euo pipefail
    export PYTHONNOUSERSITE=1

    # Reuse the model PHY selected for this family; they differ per family
    # (LG+I+G4, WAG+F+G4, JTT+G4, ...), so a single hard-coded model would be wrong.
    MODEL=\$(sed -n 's/.*Best-fit model: \\([^ ]*\\).*/\\1/p' ${phylog} | head -1)
    if [[ -z "\$MODEL" ]]; then MODEL="${params.SUBS_MODEL}+G4"; fi
    echo "${id}: fixed-topology support under \$MODEL"

    iqtree2 -s ${aln} -te ${gtree} -m "\$MODEL" \\
        --alrt 1000 --keep-ident -nt ${task.cpus} \\
        -pre ${id}.gxsup --quiet -redo

    # --alrt writes TWO values per node, "<parametric aLRT>/<SH-aLRT>"; POSSVM needs one
    # float, so keep the second. (Checked against a --alrt 0 run, which emits the
    # parametric value in the first field.)
    sed -E 's/\\)([0-9.]*)\\/([0-9.]*):/)\\2:/g' \\
        ${id}.gxsup.treefile > ${id}.generax.support.tree

    grep -q ';' ${id}.generax.support.tree || { echo "empty support tree" >&2; exit 1; }
    """
}

//workflow {
//	hg_fastas|ALN|PHY|map { id, tree, aln -> tuple(id, tree, aln, refnames_file) } | PVM
//}

species_tree_ch = Channel.value( file(params.SPECIES_TREE) )
refnames_ch     = Channel.value( file(params.REFNAMES) )

workflow {

    PHY(hg_fastas | ALN)
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

        GR_watcher(gr_input)
        gr_out = GR_watcher.out.trees


        // ---------- Branch support for the reconciled trees ----------
        // GeneRax emits no support (see the GXSUP comment). Off => PVM reads the raw
        // GeneRax tree and every orthogroup_support is the placeholder 1.0.
        gr_trees = gr_out.map { id, generax_tree, log, progress, aln ->
            tuple(id, generax_tree, aln)
        }

        if (params.gxsup) {
            GXSUP( gr_trees.join( phy_out.map { id, tree, aln, log -> tuple(id, log) } ) )
            pvm_trees = GXSUP.out.trees
        }
        else {
            pvm_trees = gr_trees
        }


        // ---------- PVM on GeneRax trees ----------
        pvm_out = pvm_trees
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
