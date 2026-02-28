nextflow.enable.dsl=2

params.resources_tsv = "${projectDir}/resources.tsv"
params.ALIGN_DIR     = "${projectDir}/results/align"
params.TREE_DIR      = "${projectDir}/results/gene_trees"
params.OUTDIR        = "${projectDir}/results"


species_tree_ch = Channel.value( file(params.SPECIES_TREE) )

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
process GR{

    tag "${id}"

    publishDir "${params.OUTDIR}/generax", mode: 'copy'

    cpus params.NCPU_GENERAX

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
    tuple val(id), path("${id}.generax.tree"), path("${id}.generax.log")

	script:
	"""
	set +e

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
		--outfile ${id}.generax.tree


	# These lines are to ignore the generax failure - may not finish going through the whole range of SPRs
	EXIT_CODE=\$?
	echo "GeneRax exit code: \$EXIT_CODE"
	
	if [[ \$EXIT_CODE -eq 10 ]]; then
    echo "GeneRax family ${id} exited with code 10 (Family parsing error)"
	fi
	exit \$EXIT_CODE
	"""
}

// -----------------------------
// Workflow
// -----------------------------
workflow {
    hg_inputs
        .combine(species_tree_ch)
        .map { id, aln, tree, species_tree ->
            tuple(id, aln, tree, species_tree)
        }
        | GR
}