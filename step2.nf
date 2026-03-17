nextflow.enable.dsl=2
params.ids           = "${projectDir}/ids.txt"
params.resources_tsv = "${projectDir}/resources.tsv"
params.run_generax   = params.run_generax ?: false
params.OUTDIR        = "${projectDir}/results"
//params.MAFFT_OPT = "--maxiterate 1000 --localpair"


params.tag_prefix = ''

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
Channel.fromPath(params.ids).splitText().map { it.trim() }.filter { it }.map { id -> tuple(id, file("${projectDir}/results/clusters/${id}.fasta")) }.set { hg_fastas }
process ALN {

    tag "${id}"

    publishDir "${params.OUTDIR}/align", mode: 'copy'

    errorStrategy = { task.attempt <= 10 ? 'retry' : 'ignore' }
    maxRetries 10
    maxErrors -1

    input:
    tuple val(id), path(fasta)

    output:
    tuple val(id), path("${id}.aln.fasta")

    script:
    def existing = file("${params.OUTDIR}/align/${id}.aln.fasta")

    if (existing.exists()) {
        """
        ln -s ${existing} ${id}.aln.fasta
        """
    }
    else {
        """
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

    errorStrategy = { task.attempt <= 10 ? 'retry' : 'ignore' }
    maxRetries 10
    maxErrors -1

    input:
    tuple val(id), path(aln)

    output:
    tuple val(id), path("${id}.treefile"), path(aln)
    path("${id}.ckp.gz"), optional: true

    script:

    def existing     = file("${params.OUTDIR}/gene_trees/${id}.treefile")
    def existing_ckp = file("${params.OUTDIR}/gene_trees/${id}.ckp.gz")

    if (existing.exists()) {
        """
        echo "Using existing tree for ${id}"
        ln -sf ${existing} ${id}.treefile
        """
    }
    else if (params.TREE_METHOD == "iqtree2" && existing_ckp.exists()) {
        """
        echo "Resuming IQ-TREE2 from checkpoint for ${id}"
        cp ${existing_ckp} ${id}.ckp.gz
        python ${projectDir}/phylogeny/main.py phylogeny \
            -f ${aln} \
            --outprefix ${id} \
            -c ${task.cpus} \
            --method ${params.TREE_METHOD} \
            --iqtree2_model ${params.IQTREE2_MODEL}
        """
    }
    else {
        """
        python ${projectDir}/phylogeny/main.py phylogeny \
            -f ${aln} \
            --outprefix ${id} \
            -c ${task.cpus} \
            --method ${params.TREE_METHOD} \
            --iqtree2_model ${params.IQTREE2_MODEL}
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
	python ${projectDir}/phylogeny/main.py possvm \
	    -t ${tree} \
	    --refsps ${params.REFSPECIES} \
	    -r ${refnames_file} \
	    -o ${id}.
	"""
}

process PVM_PREV {

	tag "${params.tag_prefix ? params.tag_prefix + '_' : ''}${id}"
	publishDir "${params.OUTDIR}/possvm_prev", mode: 'copy'

	cpus 1

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
	python ${projectDir}/phylogeny/main.py possvm \
	    -t ${tree} \
	    --refsps ${params.REFSPECIES} \
	    -r ${refnames_file} \
	    -o ${id}.
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
    """
    python ${projectDir}/workflow/report_step2.py \
        --possvm_dir      ${params.OUTDIR}/possvm \
        --possvm_prev_dir ${params.OUTDIR}/possvm_prev \
        --search_dir      ${params.OUTDIR}/search \
        --cluster_dir     ${params.OUTDIR}/clusters \
        --family_info     ${projectDir}/data/gene_families_searchinfo.csv \
        --species_tree    ${projectDir}/data/species_tree.full.newick \
        --output report_step2.html
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

        # dummy files so outputs exist
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
                sleep 10
            done
        }

        python ${projectDir}/phylogeny/main.py generax \
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

//workflow {
//	hg_fastas|ALN|PHY|map { id, tree, aln -> tuple(id, tree, aln, refnames_file) } | PVM
//}

species_tree_ch = Channel.value( file(params.SPECIES_TREE) )
refnames_ch     = Channel.value( file(params.REFNAMES) )

workflow {

    phy_out = hg_fastas | ALN | PHY

    if (params.run_generax) {

        // ---------- PVM on original trees ----------
        pvm_prev_out = phy_out
            .map { id, tree, aln -> tuple(id, tree, aln) }
            .combine(refnames_ch)
            | PVM_PREV


        // ---------- GeneRax ----------
        gr_input = phy_out
            .map { id, tree, aln ->
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
            .map { id, tree, aln ->
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