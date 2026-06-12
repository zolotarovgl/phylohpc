#!/usr/bin/env python3
import os, sys, argparse, subprocess

def parse_bash_config(path):
	cfg = {}
	with open(path) as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			if "=" in line:
				key, val = line.split("=", 1)
				cfg[key.strip()] = val.strip()
	return cfg

def main():
	parser = argparse.ArgumentParser(description="Run HG step without scheduler")
	parser.add_argument("configfile")
	parser.add_argument("--input")
	parser.add_argument("--output")
	parser.add_argument("--pref", required=True)
	parser.add_argument("--family", required=True)
	parser.add_argument("--hg", required=True)
	parser.add_argument("--mode", choices=["align","phylogeny","possvm","generax"], required=True)
	parser.add_argument("--cpus", type=int, default=4)
	parser.add_argument("--mafft", default="--maxiterate 1000 --localpair")
	args = parser.parse_args()

	config = parse_bash_config(args.configfile)

	ALIGN_DIR = config.get("ALIGN_DIR")
	TREE_DIR = config.get("TREE_DIR")
	TREE_METHOD_BIG = config.get("TREE_METHOD_BIG","iqtree")
	REFSPECIES = config.get("REFSPECIES")
	REFNAMES = config.get("REFNAMES")
	SPECIES_TREE = config.get("SPECIES_TREE")

	PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
	PHYLO_MAIN = os.path.join(PROJECT_DIR, "phylogeny", "main.py")

	PREF, FAMILY, HG = args.pref, args.family, args.hg

	if args.mode == "align":
		if not args.input or not args.output:
			sys.exit("Align mode requires --input and --output")
		if not os.path.isfile(args.input):
			sys.exit(f"Error: file {args.input} does not exist")
		cmd = [
			"python", PHYLO_MAIN, "align",
			"-f", args.input,
			"-o", args.output,
			"-c", str(args.cpus),
			"-m", args.mafft
		]

	elif args.mode == "phylogeny":
		if not args.input or not args.output:
			sys.exit("Phylogeny mode requires --input and --output")
		if not os.path.isfile(args.input):
			sys.exit(f"Error: file {args.input} does not exist")
		cmd = [
			"python", PHYLO_MAIN, "phylogeny",
			"-f", args.input,
			"--outprefix", args.output,
			"-c", str(args.cpus),
			"--method", TREE_METHOD_BIG
		]

	elif args.mode == "possvm":
		if not args.input:
			sys.exit("Possvm mode requires --input (tree file)")
		if not os.path.isfile(args.input):
			sys.exit(f"Error: file {args.input} does not exist")
		if not REFSPECIES or not REFNAMES:
			sys.exit("Error: REFSPECIES and REFNAMES must be set in config")
		cmd = [
			"python", PHYLO_MAIN, "possvm",
			"-t", args.input,
			"--refsps", REFSPECIES,
			"-r", REFNAMES,
			"-o", args.output if args.output else f"{PREF}.{FAMILY}.{HG}."
		]

	elif args.mode == "generax":
		if not args.input or not args.output:
			sys.exit("Generax mode requires --input (alignment) and --output (tree)")
		if not os.path.isfile(SPECIES_TREE):
			sys.exit(f"Error: species tree not found: {SPECIES_TREE}")
		cmd = [
			"python", PHYLO_MAIN, "generax",
			"--alignment", args.input,
			"--gene_tree", args.output + ".gene_tree",
			"--species_tree", SPECIES_TREE,
			"--output_dir", args.output + "_dir",
			"--subs_model", "LG",
			"--name", f"{PREF}.{FAMILY}.{HG}",
			"--outfile", args.output,
			"-c", str(args.cpus)
		]

	else:
		sys.exit("Unknown mode")

	print("Running:", " ".join(cmd))
	res = subprocess.run(cmd)
	sys.exit(res.returncode)

if __name__ == "__main__":
	main()