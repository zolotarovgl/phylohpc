configfile: "config/step1.yaml"

shell.executable("/bin/bash")

################################
# Parse genefam table
################################
families = []
prefs = {}

with open(config["genefam_info"]) as f:
	for raw in f:
		line = raw.strip()
		if not line or line.startswith("#"):
			continue
		cols = [c.strip() for c in raw.rstrip("\n").split("\t")]
		cols = [c for c in cols if c]
		family = cols[0]
		pref = cols[-1]
		families.append(family)
		prefs[family] = pref


cluster_targets = [
	f"{config['cluster_dir']}/{prefs[f]}.{f}_cluster.tsv"
	for f in families
]


################################
# Targets
################################
rule all_clusters:
	input:
		expand(
			f"{config['cluster_dir']}/{{pref}}.{{family}}_cluster.tsv",
			zip,
			pref=[prefs[f] for f in families],
			family=families,
		)


USE_PREPARE = bool(int(config.get("use_prepare", 0)))

if USE_PREPARE:
	rule prepare:
		input:
			config["species_list"]
		output:
			config["infasta"]
		params:
			db_dir=config.get("prepare_db_dir", "")
		shell:
			"""
			mkdir -p $(dirname {output})
			bash workflow/prepare_fasta.sh {input} {output} {params.db_dir}
			"""

SEARCH_FASTA = rules.prepare.output[0] if USE_PREPARE else config["infasta"]


################################
# SEARCH
################################
rule search:
	input:
		fasta=SEARCH_FASTA
	output:
		fasta=f"{config['search_dir']}/{{pref}}.{{family}}.domains.fasta",
		domains=f"{config['search_dir']}/{{pref}}.{{family}}.domains.csv",
		domains_ummerged=f"{config['search_dir']}/{{pref}}.{{family}}.domains_ummerged.csv",
		genes=f"{config['search_dir']}/{{pref}}.{{family}}.genes.list"
	threads: int(config["s1_ncpu"])
	resources:
		mem_mb=1000,
		runtime=10
	params:
		pfam=config["pfam_db"],
		genefam=config["genefam_info"],
		expand=config["domain_expand"],
		search_dir=config["search_dir"]
	shell:
		"""
		mkdir -p {params.search_dir}
		export PYTHONNOUSERSITE=1

		echo "Running hmmsearch for {wildcards.family}"

		python phylogeny/main.py hmmsearch \
			-f {input.fasta} \
			-g {params.genefam} \
			{wildcards.family} \
			-o {params.search_dir} \
			--pfam_db {params.pfam} \
			--domain_expand {params.expand} \
			--ncpu {threads}

		touch {output.fasta} {output.domains} {output.genes}
		"""


################################
# CLUSTER
################################
rule cluster:
	input:
		fasta=lambda wc: f"{config['search_dir']}/{wc.pref}.{wc.family}.domains.fasta"
	output:
		cluster=f"{config['cluster_dir']}/{{pref}}.{{family}}_cluster.tsv"
	threads: int(config["s2_ncpu"])
	resources:
		mem_mb=1000,
		runtime=10
	params:
		cluster_dir=config["cluster_dir"],
		max_n=config["max_n"],
		inflation=config["s2_inflation"]
	run:
		import os
		from pathlib import Path

		Path(params.cluster_dir).mkdir(parents=True, exist_ok=True)

		if os.path.getsize(input.fasta) == 0:
			print(f"Skipping clustering: {input.fasta} is empty")
			shell("touch {output.cluster}")
		else:
			shell(
				"""
				export PYTHONNOUSERSITE=1
				echo "Clustering {input.fasta}"

				python phylogeny/main.py cluster \
					-f {input.fasta} \
					--out_file {output.cluster} \
					-c {threads} \
					-m {params.max_n} \
					-i {params.inflation}

				samtools faidx {input.fasta}
				while read -r ID; do
					[ -n "$ID" ] || continue
					xargs samtools faidx {input.fasta} < <(awk -v ID="$ID" '$1==ID {{ print $2 }}' {output.cluster}) > {params.cluster_dir}/{wildcards.pref}.{wildcards.family}.${{ID}}.fasta
				done < <(cut -f 1 {output.cluster} | sort -u)
				"""
			)
