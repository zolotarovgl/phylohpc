nextflow.enable.dsl=2
params.ids = "${projectDir}/ids.txt"
params.config = "${projectDir}/configs/config.txt"
params.resources_tsv = "${projectDir}/resources.tsv"
def cfg = [:]
file(params.config).eachLine { line ->
	line = line.trim()
	if( line && !line.startsWith('#') && line.contains('=') ) {
		def (k,v) = line.split('=',2)
		cfg[k.trim()] = v.trim()
	}
}
params.ALIGN_DIR = cfg.ALIGN_DIR
params.TREE_DIR = cfg.TREE_DIR
params.TREE_METHOD = 'iqtree2'
params.SPECIES_TREE = cfg.SPECIES_TREE
params.REFSPECIES = cfg.REFSPECIES
params.REFNAMES = cfg.REFNAMES
params.OUTDIR = "${projectDir}/results"
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
		if( cols.size() > 1 && cols[1] ) m.aln_mem = cols[1] as nextflow.util.MemoryUnit
		if( cols.size() > 2 && cols[2] ) m.aln_time = cols[2] as nextflow.util.Duration
		if( cols.size() > 3 && cols[3] ) m.phy_mem = cols[3] as nextflow.util.MemoryUnit
		if( cols.size() > 4 && cols[4] ) m.phy_time = cols[4] as nextflow.util.Duration
		if( cols.size() > 5 && cols[5] ) m.pvm_mem = cols[5] as nextflow.util.MemoryUnit
		if( cols.size() > 6 && cols[6] ) m.pvm_time = cols[6] as nextflow.util.Duration
		res[rid] = m
	}
}
Channel.fromPath(params.ids).splitText().map { it.trim() }.filter { it }.map { id -> tuple(id, file("${projectDir}/results/clusters/${id}.fasta")) }.set { hg_fastas }
process ALN {
	maxForks 50
	tag "test_${id}"
	publishDir "${params.OUTDIR}/align", mode: 'copy'
	cpus 8
	memory {
		def base = res[id]?.aln_mem ?: 500.MB
		return base + (task.attempt - 1) * 1.GB
	}
	time {
		def base = res[id]?.aln_time ?: 30.min
		return base + (task.attempt - 1) * 30.min
	}
	errorStrategy = { task.attempt <= 10 ? 'retry' : 'ignore' }
	maxRetries 10
	maxErrors -1
	input:
	tuple val(id), path(fasta)
	output:
	tuple val(id), path("${id}.aln.fasta")
	script:
	"""
	python ${projectDir}/phylogeny/main.py align -f ${fasta} -o ${id}.aln.fasta -c ${task.cpus} -m "--maxiterate 1000 --localpair"
	python ${projectDir}/workflow/remove_gaponly.py ${id}.aln.fasta ${id}.aln.fasta_tmp
	mv ${id}.aln.fasta_tmp ${id}.aln.fasta
	"""
}
process PHY {
	maxForks 50
	tag "test_${id}"
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
	script:
	"""
	python ${projectDir}/phylogeny/main.py phylogeny -f ${aln} --outprefix ${id} -c ${task.cpus} --method ${params.TREE_METHOD}
	"""
}
process PVM {
	tag "test_${id}"
	publishDir "${params.OUTDIR}/possvm", mode: 'copy'
	cpus 1
	memory {
		def base = res[id]?.pvm_mem ?: 300.MB
		return base * task.attempt
	}
	time {
		def base = res[id]?.pvm_time ?: 5.min
		return base * task.attempt
	}
	errorStrategy 'ignore'
	maxRetries 3
	input:
	tuple val(id), path(tree), path(aln), path(refnames_file)
	output:
	tuple val(id), path(tree), path("${id}.treefile.ortholog_groups.newick"), path("${id}.treefile.ortholog_groups.csv"), path("${id}.treefile.pairs_orthologs.csv")
	script:
	"""
	python ${projectDir}/phylogeny/main.py possvm -t ${tree} --refsps ${params.REFSPECIES} -r ${refnames_file} -o ${id}.
	"""
}
refnames_file = file(params.REFNAMES)
workflow {
	hg_fastas|ALN|PHY|map { id, tree, aln -> tuple(id, tree, aln, refnames_file) }|PVM
}