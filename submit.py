#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess

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
    parser = argparse.ArgumentParser()
    parser.add_argument("configfile", help="Bash-style config file with KEY=VALUE lines")
    parser.add_argument("--pref", required=True, help="Prefix (PREF)")
    parser.add_argument("--family", required=True, help="Gene family (FAMILY)")
    parser.add_argument("--hg", required=True, help="HG identifier")
    parser.add_argument("--mode", choices=["align", "phylogeny", "possvm"], required=True, help="Which step to run")
    parser.add_argument("--cpus", type=int, default=4, help="CPUs per task")
    parser.add_argument("--mem", default="1G", help="Memory")
    parser.add_argument("--time", default="1:00:00", help="Walltime")
    parser.add_argument("-n", "--dry-run", action="store_true", help="Dry run: only print the sbatch command")
    parser.add_argument("--mafft",default = "--maxiterate 1000 --localpair", help = "MAFFT settings")
    args = parser.parse_args()

    config = parse_bash_config(args.configfile)

    ALIGN_DIR = config["ALIGN_DIR"]
    TREE_DIR = config["TREE_DIR"]
    TREE_METHOD_BIG = config.get("TREE_METHOD_BIG", "iqtree")
    REFSPECIES = config.get("REFSPECIES")
    REFNAMES = config.get("REFNAMES")

    PREF = args.pref
    FAMILY = args.family
    HG = args.hg
    jobname = ".".join([args.mode, PREF, FAMILY, HG])

    if args.mode == "phylogeny":
        in_fasta = f"{ALIGN_DIR}/{PREF}.{FAMILY}.{HG}.aln.fasta"
        out_prefix = f"{TREE_DIR}/{PREF}.{FAMILY}.{HG}"
        if not os.path.isfile(in_fasta):
            sys.exit(f"Error: file {in_fasta} does not exist")
        wrap = f"python phylogeny/main.py phylogeny -f {in_fasta} --outprefix {out_prefix} -c {args.cpus} --method {TREE_METHOD_BIG}"

    elif args.mode == "align":
        in_fasta = f"results/clusters/{PREF}.{FAMILY}.{HG}.fasta"
        out_fasta = f"{ALIGN_DIR}/{PREF}.{FAMILY}.{HG}.aln.fasta"
        if not os.path.isfile(in_fasta):
            sys.exit(f"Error: file {in_fasta} does not exist")
        wrap = f'python phylogeny/main.py align -f {in_fasta} -o {out_fasta} -c {args.cpus} -m "{args.mafft}"'

    elif args.mode == "possvm":
        tree_file = f"{TREE_DIR}/{PREF}.{FAMILY}.{HG}.treefile"
        if not os.path.isfile(tree_file):
            sys.exit(f"Error: file {tree_file} does not exist")
        if not REFSPECIES or not REFNAMES:
            sys.exit("Error: REFSPECIES and REFNAMES must be set in the config file for poss mode")
        out_prefix = f"{PREF}.{FAMILY}.{HG}."
        wrap = f"python phylogeny/main.py possvm -t {tree_file} --refsps {REFSPECIES} -r {REFNAMES} -o {out_prefix}"

    cmd = [
        "sbatch",
        f"--job-name={jobname}",
        f"--cpus-per-task={args.cpus}",
        f"--mem={args.mem}",
        f"--time={args.time}",
        "--wrap", wrap
    ]
    print(" ".join(cmd))
    if not args.dry_run:
        subprocess.run(cmd)

if __name__ == "__main__":
    main()

