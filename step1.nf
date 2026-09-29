nextflow.enable.dsl=2


// Single source of truth for these params is nextflow.config.
// Only pref_family (not in the config) keeps an inline default here.
params.pref_family   = null
// Rebuild every family even if its outputs already exist (default: skip complete ones).
params.redo          = false

workflow {

    if( !params.genefam_info )
        error "Please provide --genefam_info"

    if( !params.infasta )
        error "Please provide --infasta"

    def families_ch = Channel
        .fromPath(params.genefam_info)
        .splitText()
        .filter { it.trim() }
        .map { line ->
            def cols   = line.trim().split('\t')
            def family = cols[0].trim()
            def pref   = cols[-1].trim()
            tuple(pref, family)
        }

    def genefam_ch = Channel.value(file(params.genefam_info))
    def infasta_ch = Channel.value(file(params.infasta))
    // Custom HMMs: staged as a path input so an edited .hmm invalidates the cache.
    // [] is Nextflow's placeholder for an absent optional path input.
    def hmmdir_ch  = Channel.value(params.hmm_dir
        ? file(params.hmm_dir, type: 'dir', checkIfExists: true)
        : [])

    // ---- skip families that are already complete -------------------------------
    // The genefam table is a per-task INPUT to SEARCH (see the process below), so
    // Nextflow hashes its CONTENT into every task's cache key: appending one family
    // row invalidates -resume for ALL of them and re-searches the whole table.
    // Adding a family must cost one family's compute, so completeness is decided
    // from the published outputs instead of from the resume cache.
    //
    // A family is DONE in one of two ways, and both must count or the second kind
    // re-runs forever:
    //   1. it produced domain hits and was clustered  -> cluster tsv exists
    //   2. it produced NO domain hits, so there was nothing to cluster -> the
    //      domains fasta exists and is empty. 31 families are legitimately in this
    //      state (plant/fungal TF families absent from animals: WRKY, GRAS, YABBY,
    //      NAC, AP2, TCP, B3, zf-Dof, EIN3, FLO_LFY, SBP, APSES, ...).
    // Pass --redo true to rebuild everything regardless.
    def redo = params.redo as boolean

    def familyIsDone = { String pref, String family ->
        def dom = file("${params.search_dir}/${pref}.${family}.domains.fasta")
        def clu = file("${params.cluster_dir}/${pref}.${family}_cluster.tsv")
        dom.exists() && ( dom.size() == 0 || clu.exists() )
    }

    def todo_ch = families_ch.filter { pref, family ->
        def done = familyIsDone(pref, family)
        if( done && !redo ) log.info "skip (already complete): ${pref}.${family}"
        redo || !done
    }

    def search = SEARCH(todo_ch, genefam_ch, infasta_ch, hmmdir_ch)

    def nonempty = search.main
        .filter { pref, family, fasta -> fasta && fasta.size() > 0 }

    CLUSTER(nonempty)
}

process SEARCH {



  stageInMode 'copy'
  tag "${pref}.${family}"

  cpus   { params.s1_ncpu as int }
  memory { 500.MB + (task.attempt - 1) * 500.MB }
  time   { 5.min + (task.attempt - 1) * 10.min }

  errorStrategy = { task.attempt <= 5 ? 'retry' : 'terminate' }
  maxRetries 5

  input:
  tuple val(pref), val(family)
  path(genefam_info, stageAs: 'genefam.csv')
  path(infasta,      stageAs: 'input.fasta')
  path(hmm_dir,      stageAs: 'custom_hmms')

  output:
  tuple val(pref), val(family),
        path("${pref}.${family}.domains.fasta"),
        emit: main

  tuple val(pref), val(family),
        path("${pref}.${family}.domains.csv", optional: true),
        path("${pref}.${family}.domains_ummerged.csv", optional: true),
        path("${pref}.${family}.genes.list",  optional: true),
        emit: aux

  publishDir "${params.search_dir}", mode: 'copy'

  script:
  """
	set -e
	export PYTHONNOUSERSITE=1
	echo "Running hmmsearch for ${family}"

	python ${projectDir}/phylogeny/main.py hmmsearch \
		-f input.fasta \
		-g genefam.csv \
		${family} \
		-o . \
		--pfam_db ${params.pfam_db} \
		${params.hmm_dir ? '--hmm_dir custom_hmms' : ''} \
		--domain_expand ${params.domain_expand} \
		--ncpu ${task.cpus}

	touch ${pref}.${family}.domains.fasta 
	touch ${pref}.${family}.domains.csv 
	touch ${pref}.${family}.domains_ummerged.csv
	touch ${pref}.${family}.genes.list
  """
}
process CLUSTER {

    tag "${pref}.${family}"

    cpus   { params.s2_ncpu as int }
    memory { 300.MB + (task.attempt - 1) * 1.GB }
    time   { 10.min + (task.attempt - 1) * 1.h }

    errorStrategy = { task.attempt <= 5 ? 'retry' : 'terminate' }
    maxRetries 5

    input:
    tuple val(pref), val(family), path(domains_fasta)

    output:
    path("${pref}.${family}_cluster.tsv")
    path("${pref}.${family}.*.fasta")

    publishDir "${params.cluster_dir}", mode: 'copy'

	script:
	"""
	echo "Clustering: ${domains_fasta}"
	export PYTHONNOUSERSITE=1

	python ${projectDir}/phylogeny/main.py cluster \
		-f ${domains_fasta} \
		--out_file ${pref}.${family}_cluster.tsv \
		-c ${task.cpus} \
		-m ${params.max_n} \
		-i ${params.s2_inflation}

	samtools faidx ${domains_fasta}
	while read -r ID; do
		[ -n "\$ID" ] || continue
		xargs samtools faidx ${domains_fasta} < <(awk -v ID="\$ID" '\$1==ID { print \$2 }' ${pref}.${family}_cluster.tsv) > ${pref}.${family}.\${ID}.fasta
	done < <(cut -f 1 ${pref}.${family}_cluster.tsv | sort -u)

	# Guarantee structural outputs
	touch ${pref}.${family}_cluster.tsv
	"""
}
