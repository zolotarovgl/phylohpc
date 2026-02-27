nextflow.enable.dsl=2

params.pref_family   = null
params.genefam_info  = null
params.infasta       = null
params.search_dir    = "results/search"
params.cluster_dir   = "results/clusters"
params.max_n         = ""
params.pfam_db       = null
params.domain_expand = 30
params.s2_ncpu      = 4
params.s2_inflation = 2.0

workflow {

    if( !params.genefam_info )
        error "Please provide --genefam_info"

    if( !params.infasta )
        error "Please provide --infasta"

    // Channel of prefix + family parsed from file
	families_ch = Channel
		.fromPath(params.genefam_info)
		.splitText()
		.filter { it.trim() }
		.map { line ->
			def cols = line.trim().split('\t')
			def family = cols[0].trim()
			def pref   = cols[-1].trim()
			tuple(pref, family)
		}

	SEARCH(
		families_ch,
		file(params.genefam_info),
		file(params.infasta)
	)
	.filter { pref, family, fasta -> fasta != null && fasta.size() > 0 }
	| CLUSTER
}


process SEARCH {

    tag "${pref}.${family}"

    input:
    tuple val(pref), val(family)
    path genefam_info
    path infasta

	output:
	tuple val(pref), val(family),
		path("${pref}.${family}.domains.fasta", optional: true)
    publishDir "${params.search_dir}", mode: 'copy'

    script:
    """
	echo "Running hmmsearch for ${family}"

	python ${projectDir}/phylogeny/main.py hmmsearch \
		-f ${infasta} \
		-g ${genefam_info} \
		${family} \
		-o . \
		--pfam_db ${params.pfam_db} \
		--domain_expand ${params.domain_expand} \
		--ncpu ${task.cpus}

	OUTFILE=${pref}.${family}.domains.fasta

	# Ensure file exists so Nextflow does not crash
	if [[ ! -f "\$OUTFILE" ]]; then
		echo "No domains found for ${family} — creating empty file"
		touch "\$OUTFILE"
	fi
	"""
}
process CLUSTER {

    tag "${pref}.${family}"

    input:
    tuple val(pref), val(family), path(domains_fasta)

    output:
    path("${pref}.${family}_cluster.tsv")
    path("${pref}.${family}.*.fasta")

    publishDir "${params.cluster_dir}", mode: 'copy'

    script:
    """
    echo "Clustering: ${domains_fasta}"

    python ${projectDir}/phylogeny/main.py cluster \
        -f ${domains_fasta} \
        --out_file ${pref}.${family}_cluster.tsv \
        -c ${params.s2_ncpu} \
        -m ${params.max_n} \
        -i ${params.s2_inflation}

    echo "Indexing fasta"
    samtools faidx ${domains_fasta}

    OUTPREFIX=\$(echo ${pref}.${family}_cluster.tsv | sed -E 's/_cluster.tsv\$//g')

    for ID in \$(cut -f 1 ${pref}.${family}_cluster.tsv | sort | uniq); do
        awk -v ID=\${ID} '\$1==ID { print \$2 }' ${pref}.${family}_cluster.tsv \
        | xargs samtools faidx ${domains_fasta} \
        > \${OUTPREFIX}.\${ID}.fasta

        echo "Created \${OUTPREFIX}.\${ID}.fasta"
    done
    """
}