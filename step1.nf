nextflow.enable.dsl=2


// Defaults 
params.pref_family   = null
params.genefam_info  = "genefam.csv"
params.infasta       = "data/input.fasta"
params.search_dir    = "results/search"
params.cluster_dir   = "results/clusters"
params.max_n         = 3000
params.pfam_db       = "/users/asebe/xgraubove/data/pfam/Pfam-A.hmm"
params.domain_expand = 30
params.s1_ncpu       = 4
params.s2_ncpu       = 4
params.s2_inflation  = 1.1

// ── Failure diagnostics ───────────────────────────────────────────────────────
// Log the work dir + tail of .command.err whenever a task fails (any attempt),
// so a failure is readable in the run log without digging into work dirs.
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
    log.info "step1 ${workflow.success ? 'completed OK' : 'FAILED'} after ${workflow.duration}"
    if( !workflow.success )
        log.info "List failed tasks:  nextflow log ${workflow.runName} -f name,status,exit,workdir -F \"status=='FAILED'\""
}

// ── Environment preflight ─────────────────────────────────────────────────────
// Runs once before the search swarm and fails the whole run with one clear
// message if Python deps or tools are missing (e.g. the conda env wasn't activated).
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
    for tool in hmmsearch mafft mcl samtools; do
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

workflow {

    if( !params.genefam_info )
        error "Please provide --genefam_info"

    if( !params.infasta )
        error "Please provide --infasta"

    // Fail fast with one clear message if the environment is not ready.
    def ready = PREFLIGHT()

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
        .combine(ready)
        .map { pref, family, ok -> tuple(pref, family) }

    def genefam_ch = Channel.value(file(params.genefam_info))
    def infasta_ch = Channel.value(file(params.infasta))

    def search = SEARCH(families_ch, genefam_ch, infasta_ch)

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

  errorStrategy = { errReport(task); task.attempt <= 5 ? 'retry' : 'terminate' }
  maxRetries 5

  input:
  tuple val(pref), val(family)
  path(genefam_info, stageAs: 'genefam.csv')
  path(infasta,      stageAs: 'input.fasta')

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

    errorStrategy = { errReport(task); task.attempt <= 5 ? 'retry' : 'terminate' }
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
