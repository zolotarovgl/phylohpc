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

process ALIGN {

	tag "$id"

	publishDir "${projectDir}/results/align", mode: 'copy'

	cpus 8
	memory { 500.MB + (task.attempt - 1) * 1.GB }
	time { 30.min + (task.attempt - 1) * 1.h }

	errorStrategy 'retry'
	maxRetries 5

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
		-m ""
	python ${projectDir}/workflow/remove_gaponly.py ${id}.aln.fasta ${id}.aln.fasta_tmp
	mv ${id}.aln.fasta_tmp ${id}.aln.fasta
	"""
}

process PHYLOGENY {

	tag "$id"

	publishDir "${projectDir}/results/gene_trees", mode: 'copy'

	cpus 4
	memory { 500.MB + (task.attempt - 1) * 1.GB }
	time { 1.h + (task.attempt - 1) * 2.h }

	errorStrategy 'retry'
	maxRetries 3

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

process POSSVM {

	tag "$id"

	publishDir "${projectDir}/results/possvm", mode: 'copy'

	cpus 1
	memory { 500.MB * task.attempt }
	time   { 5.min  * task.attempt }

	errorStrategy 'retry'
	maxRetries 5

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
		| ALIGN \
		| PHYLOGENY \
		| map { id, tree, aln -> tuple(id, tree, aln, refnames_file) } \
		| POSSVM
}