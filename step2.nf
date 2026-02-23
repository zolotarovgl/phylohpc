nextflow.enable.dsl=2

params.ids    = "${projectDir}/ids.txt"
params.config = "${projectDir}/configs/config.txt"

def cfg = [:]
file(params.config).eachLine { line ->
	line = line.trim()
	if( line && !line.startsWith('#') && line.contains('=') ) {
		def (k,v) = line.split('=',2)
		cfg[k.trim()] = v.trim()
	}
}

params.ALIGN_DIR    = cfg.ALIGN_DIR
params.TREE_DIR     = cfg.TREE_DIR
params.TREE_METHOD  = 'iqtree2'
params.SPECIES_TREE = cfg.SPECIES_TREE
params.REFSPECIES   = cfg.REFSPECIES
params.REFNAMES     = cfg.REFNAMES

Channel
	.fromPath(params.ids)
	.splitText()
	.map { it.trim() }
	.filter { it }
	.map { id ->
		tuple(id, file("${projectDir}/results/clusters/${id}.fasta"))
	}
	.set { hg_fastas }

process ALN {

	maxForks 50

	tag "$id"

	publishDir "${projectDir}/results/align", mode: 'copy'

	cpus 8

	memory {
		def base = 500.MB
		def prev_exit = task.attempt > 1 ? task.previousTrace?.exit : null

		if( prev_exit in [1,137] )
			return base + (task.attempt - 1) * 2.GB
		else
			return base
	}

	time {
		def base = 5.min
		def prev_exit = task.attempt > 1 ? task.previousTrace?.exit : null

		if( prev_exit == 143 )
			return base + (task.attempt - 1) * 30.min
		else
			return base
	}

	errorStrategy = {
		if (task.attempt <= 10) {
			return 'retry'
		} else {
			return 'ignore' // Continue even after retries
		}
	}

	maxRetries 10
	maxErrors -1

	//Inputs 
	input:
	tuple val(id), path(fasta)

	output:
	tuple val(id), path("${id}.aln.fasta")

	script:
	"""
	python ${projectDir}/phylogeny/main.py align \
		-f ${fasta} \
		-o ${id}.aln.fasta \
		-c ${task.cpus} \
		-m "--maxiterate 1000 --localpair"
	python ${projectDir}/workflow/remove_gaponly.py ${id}.aln.fasta ${id}.aln.fasta_tmp
	mv ${id}.aln.fasta_tmp ${id}.aln.fasta
	"""
}

process PHY {

	maxForks 50

	tag "$id"

	publishDir "${projectDir}/results/gene_trees", mode: 'copy'
	
	cpus 4

	memory {
		def base = 300.MB
		if( task.attempt > 1 && task.previousTrace?.exitStatus == 137 )
			return base + (task.attempt - 1) * 500.MB
		else
			return base
	}

	time {
		def prev = task.attempt > 1 ? task.previousTrace?.exitStatus : null
		if( task.attempt == 1 )
			return 5.min
		else if( task.attempt == 2 && prev == 143 )
			return 30.min
		else if( task.attempt == 3 && prev == 143 )
			return 1.h
		else
			return 6.h
	}

	errorStrategy = {
		if (task.attempt <= 5) {
			return 'retry'
		} else {
			return 'ignore' // Continue even after retries
		}
	}
	maxRetries 5
	maxErrors -1

	input:
	tuple val(id), path(aln)

	output:
	tuple val(id), path("${id}.treefile"), path(aln)

	script:
	"""
	python ${projectDir}/phylogeny/main.py phylogeny \
		-f ${aln} \
		--outprefix ${id} \
		-c ${task.cpus} \
		--method ${params.TREE_METHOD}
	"""
}

process PVM {

	tag "$id"

	publishDir "${projectDir}/results/possvm", mode: 'copy'

	cpus 1
	memory {
		def base = 300.MB
		if( task.attempt > 1 && task.previousTrace?.exitStatus == 137 )
			return base + (task.attempt - 1) * 500.MB
		else
			return base
	}
	time {
		def base = 1.min
		if( task.attempt > 1 && task.previousTrace?.exitStatus == 143 )
			return base + (task.attempt - 1) * 5.min
		else
			return base
	}


	errorStrategy 'ignore'

	maxRetries 3

	input:
	tuple val(id), path(tree), path(aln), path(refnames_file)

	output:
		tuple val(id),
			path(tree),
			path("${id}.treefile.ortholog_groups.newick"),
			path("${id}.treefile.ortholog_groups.csv"),
			path("${id}.treefile.pairs_orthologs.csv")

	script:
	"""
	python ${projectDir}/phylogeny/main.py possvm \
		-t ${tree} \
		--refsps ${params.REFSPECIES} \
		-r ${refnames_file} \
		-o ${id}.
	"""
}

refnames_file = file(params.REFNAMES)

workflow {
	hg_fastas \
		| ALN \
		| PHY \
		| map { id, tree, aln -> tuple(id, tree, aln, refnames_file) } \
		| PVM
}