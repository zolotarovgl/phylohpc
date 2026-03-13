#!/usr/bin/env python3

import argparse
import sys
import tempfile
from pathlib import Path


def load_species_ids(id_arg):
	"""
	Accept either a species prefix or a file containing prefixes.
	Returns a list of species IDs.
	"""
	p = Path(id_arg)

	if p.exists() and p.is_file():
		with open(p) as f:
			return [line.strip() for line in f if line.strip()]
	else:
		return [id_arg]


def parse_args():
	parser = argparse.ArgumentParser(description="Gather species annotation")

	parser.add_argument("--id", required=True, help="Species prefix or file with prefixes")
	parser.add_argument("--outfile", help="Optional output file (single species only)")
	parser.add_argument("--outdir", help="Output directory for multiple species")
	parser.add_argument("--search-dir", required=True, help="Search directory")

	parser.add_argument(
		"--tree-dir",
		nargs="+",
		required=True,
		help="Directories with POSSVM outputs (searched in order of priority)",
	)

	parser.add_argument("--prefix", default=None, help="Optional family prefix")

	parser.add_argument(
		"--split-prefix",
		action="store_true",
		help="Split annotations by prefix (e.g. Nvec.tfs.tsv)",
	)

	return parser.parse_args()


def collect_tmp_anno(tree_dirs, prefix, species_id):
	result = {}

	pattern = "*groups.csv" if not prefix else f"{prefix}*groups.csv"

	for tree_dir in tree_dirs:
		tree_dir = Path(tree_dir)

		if not tree_dir.exists():
			continue

		for file in tree_dir.glob(pattern):
			with open(file) as f:
				for line in f:
					if species_id in line:
						parts = line.rstrip("\n").split("\t")

						if len(parts) >= 4:
							gene = parts[0]

							if gene not in result:
								result[gene] = parts[1:4]

	return result


def collect_ids_todo(search_dir, prefix, species_id):
	ids = set()

	pattern = "*genes.list" if not prefix else f"{prefix}*.genes.list"

	for file in Path(search_dir).glob(pattern):
		with open(file) as f:
			for line in f:
				if species_id in line:
					ids.add(line.strip().split()[0])

	return sorted(ids)


def collect_gene2class(search_dir, prefix, species_id):
	gene2class = {}

	pattern = "*genes.list" if not prefix else f"{prefix}*.genes.list"

	for file in Path(search_dir).glob(pattern):

		pref = file.name.replace(".genes.list", "")

		with open(file) as f:
			for line in f:
				if species_id in line:
					gene = line.strip().split()[0]
					gene2class[gene] = pref

	return gene2class


def collect_pep2hg(cluster_dir, species_id):
	cluster_dir = Path(cluster_dir)

	if not cluster_dir.exists():
		raise FileNotFoundError(f"Cluster directory not found: {cluster_dir}")

	mapping = {}

	for file in cluster_dir.glob("*_cluster.tsv"):

		pref = file.name.replace("_cluster.tsv", "")

		with open(file) as f:
			for line in f:
				parts = line.rstrip("\n").split("\t")

				if len(parts) < 2:
					continue

				cluster_id = parts[0]
				gene = parts[1]

				if gene.startswith(species_id):
					mapping[gene] = f"{pref}.{cluster_id}"

	return mapping


def build_result(tmp_anno, ids_todo, gene2class, pep2hg):
	results = []

	for gene in ids_todo:

		if gene in tmp_anno:
			class_info = tmp_anno[gene]
		else:
			class_info = ["Unclassified"]

		hg = pep2hg.get(gene)

		if class_info[0] == "Unclassified" and hg:
			class_info = [f"{hg}:Unclassified"]

		gene_class = gene2class.get(gene)

		if class_info[0] == "Unclassified" and gene_class:
			class_info = [f"{gene_class}:Unclassified"]

		row = [gene] + class_info
		results.append(row)

	return results


def write_split_outputs(rows, species_id, outdir):
	files = {}

	for row in rows:
		gene = row[0]
		class_info = row[1]

		if ":" in class_info:
			prefix = class_info.split(":", 1)[0].split('.')[0]
		
		else:
			prefix = "Unclassified"

		if prefix not in files:
			files[prefix] = []
		
		files[prefix].append("\t".join(row))

	for prefix, lines in files.items():
		outfile = outdir / f"{species_id}.{prefix}.tsv"
		with open(outfile, "w") as out:
			out.write("\n".join(lines) + "\n")


def main():

	args = parse_args()

	species_ids = load_species_ids(args.id)

	search_dir = args.search_dir
	tree_dirs = args.tree_dir
	prefix = args.prefix
	cluster_dir = "results/clusters"

	if len(species_ids) > 1 or args.split_prefix:
		if not args.outdir:
			sys.exit("When multiple species or --split-prefix is used, --outdir is required")

		outdir = Path(args.outdir)
		outdir.mkdir(parents=True, exist_ok=True)
	else:
		outdir = None

	with tempfile.TemporaryDirectory(prefix="gather_anno_"):

		for species_id in species_ids:
			print(species_id)
			tmp_anno = collect_tmp_anno(tree_dirs, prefix, species_id)

			ids_todo = collect_ids_todo(search_dir, prefix, species_id)

			gene2class = collect_gene2class(search_dir, prefix, species_id)

			pep2hg = collect_pep2hg(cluster_dir, species_id)

			rows = build_result(tmp_anno, ids_todo, gene2class, pep2hg)

			if args.split_prefix:
				write_split_outputs(rows, species_id, outdir)

			else:
				result = ["\t".join(r) for r in rows]

				if outdir and len(species_ids) > 1:
					outfile = outdir / f"{species_id}.annotation.tsv"
					with open(outfile, "w") as out:
						out.write("\n".join(result) + "\n")

				elif args.outfile:
					with open(args.outfile, "w") as out:
						out.write("\n".join(result) + "\n")

				else:
					print("\n".join(result))


if __name__ == "__main__":
	main()