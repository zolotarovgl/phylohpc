#!/usr/bin/env python3
import argparse, subprocess, os, sys
from workflow.helper import parse_bash_config

def read_json(json_fn):
    if not os.path.isfile(json_fn):
        print(f"WARNING: JSON file {json_fn} doesn't exist!")
        data = {}
    elif os.path.getsize(json_fn) == 0:
        print(f"WARNING: JSON file {json_fn} is empty!")
        data = {}
    else:
        with open(json_fn) as f:
            data = json.load(f)
    return(data)

def run_cmd(cmd, dry=False):
	print(" ".join(cmd))
	if not dry:
		out = subprocess.run(cmd, capture_output=True, text=True, check=True)
		return out.stdout.strip()
	return ""

def submit_family(config_fn, pref, family, time_s1=None, time_s2=None, mem_s1=None, mem_s2=None, dry=False,verbose = True):
	config = parse_bash_config(config_fn)

	PROJECT = config.get("PROJECT", "project")
	QOS_S1 = config.get("QOS_S1", "normal")
	QOS_S2 = config.get("QOS_S2", "normal")
	LOG_DIR = config.get("LOG_DIR", "logs")
	SEARCH_DIR = config.get("SEARCH_DIR", "results/search")
	CLUSTER_DIR = config.get("CLUSTER_DIR", "results/cluster")
	GENEFAM_INFO = config.get("GENEFAM_INFO")
	INFASTA = config.get("INFASTA")
	MAX_N = config.get("MAX_N", "")
	TIME_S1 = time_s1 or config.get("TIME_S1", "01:00:00")
	TIME_S2 = time_s2 or config.get("TIME_S2", "01:00:00")
	MEM_S1 = mem_s1 or config.get("MEM_S1", "4G")
	MEM_S2 = mem_s2 or config.get("MEM_S2", "8G")

	#print(f"Prefix: {pref}; Family: {family}")
	#print(f"MEM_S1: {MEM_S1}\nTIME_S1: {TIME_S1}")
	#print(f"MEM_S2: {MEM_S2}\nTIME_S2: {TIME_S2}")



	#if dry:
	#	print("Dry run mode: no jobs will be submitted.")
	#	return

	job_name1 = f"{PROJECT}.s1.{pref}.{family}"
	cmd1 = [
		"sbatch", f"--time={TIME_S1}", f"--qos={QOS_S1}", f"--mem={MEM_S1}",
		f"--job-name={job_name1}",
		f"--output={LOG_DIR}/s1/%x.%j.out",
		f"--error={LOG_DIR}/s1/%x.%j.out",
		"workflow/s01_search.sh", family, GENEFAM_INFO, INFASTA, SEARCH_DIR
	]
	if verbose:
		print(' '.join(cmd1))
	
	if not dry:
		out1 = run_cmd(cmd1, dry)
	else:
		out1 = None
	jid1 = out1.split()[-1] if out1 else "DRYRUN"
	if not dry:
		print(f"Submitted {jid1}: {job_name1}")

	job_name2 = f"{PROJECT}.s2.{pref}.{family}"
	cmd2 = [
		"sbatch", f"--qos={QOS_S2}", f"--time={TIME_S2}",
		f"--dependency=afterok:{jid1}", f"--mem={MEM_S2}",
		f"--job-name={job_name2}",
		f"--output={LOG_DIR}/s2/%x.%j.out",
		f"--error={LOG_DIR}/s2/%x.%j.out",
		"workflow/s02_cluster.sh",
		f"{SEARCH_DIR}/{pref}.{family}.domains.fasta",
		f"{CLUSTER_DIR}/{pref}.{family}_cluster.tsv",
		str(MAX_N)
	]
	if verbose:
		print(' '.join(cmd2))
	if not dry:
		out2 = run_cmd(cmd2, dry)
	else:
		out2 = None
	jid2 = out2.split()[-1] if out2 else "DRYRUN"
	if not dry:
		print(f"Submitted {jid2}: {job_name2}")
	job_ids = {'search' : jid1, 'cluster': jid2}
	return(job_ids)


def main():
	parser = argparse.ArgumentParser(description="Submit hmmsearch and/or clustering for a family.")
	parser.add_argument("config", help="Path to config file")
	parser.add_argument("pref_family", help="Prefix.Family or separate values")
	parser.add_argument("--time_s1")
	parser.add_argument("--time_s2")
	parser.add_argument("--mem_s1")
	parser.add_argument("--mem_s2")
	parser.add_argument("--dry", "-n", action="store_true", help="Dry run")
	args = parser.parse_args()

	verbose = True
	if "." in args.pref_family and not args.pref_family.startswith("--"):
		pref, family = args.pref_family.split(".", 1)
	else:
		sys.exit("Error: please use PREFIX.FAMILY format")

	submit_family(
		config_fn=args.config,
		pref=pref,
		family=family,
		time_s1=args.time_s1,
		time_s2=args.time_s2,
		mem_s1=args.mem_s1,
		mem_s2=args.mem_s2,
		dry=args.dry
	)

if __name__ == "__main__":
	main()
