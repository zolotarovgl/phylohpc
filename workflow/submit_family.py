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

def run_cmd(cmd, dry=False,verbose = True):
	if verbose:
		print(" ".join(cmd))
	if not dry:
		out = subprocess.run(cmd, capture_output=True, text=True, check=True)
		return out.stdout.strip()
	return ""



def submit_search(pref,family,config,time = None,mem = None,verbose = True, dry = False):

	PROJECT = config.get("PROJECT", "project")
	QOS_S1 = config.get("QOS_S1", "normal")
	QOS_S2 = config.get("QOS_S2", "normal")
	LOG_DIR = config.get("LOG_DIR", "logs")
	SEARCH_DIR = config.get("SEARCH_DIR", "results/search")
	CLUSTER_DIR = config.get("CLUSTER_DIR", "results/cluster")
	GENEFAM_INFO = config.get("GENEFAM_INFO")
	INFASTA = config.get("INFASTA")
	MAX_N = config.get("MAX_N", "")
	TIME = time or config.get("TIME_S1", "01:00:00")
	MEM = mem or config.get("MEM_S1", "4G")
	

	job_name = f"{PROJECT}.s1.{pref}.{family}"
	cmd = [
		"sbatch", f"--time={TIME}", f"--qos={QOS_S1}", f"--mem={MEM}",
		f"--job-name={job_name}",
		f"--output={LOG_DIR}/s1/%x.%j.out",
		f"--error={LOG_DIR}/s1/%x.%j.out",
		"workflow/s01_search.sh", family, GENEFAM_INFO, INFASTA, SEARCH_DIR
	]
	if verbose:
		print(' '.join(cmd))
	
	if not dry:
		out = run_cmd(cmd, dry, verbose)
	else:
		out = None
	jid = out.split()[-1] if out else "DRYRUN"
	if not dry and verbose:
		print(f"Submitted {jid}: {job_name}")
	return(jid)

def submit_cluster(pref,family,config,after_jid = None,time = None,mem = None,verbose = True, dry = False):

	PROJECT = config.get("PROJECT", "project")
	QOS_S1 = config.get("QOS_S1", "normal")
	QOS_S2 = config.get("QOS_S2", "normal")
	LOG_DIR = config.get("LOG_DIR", "logs")
	SEARCH_DIR = config.get("SEARCH_DIR", "results/search")
	CLUSTER_DIR = config.get("CLUSTER_DIR", "results/cluster")
	GENEFAM_INFO = config.get("GENEFAM_INFO")
	INFASTA = config.get("INFASTA")
	MAX_N = config.get("MAX_N", "")
	TIME = time or config.get("TIME_S2", "01:00:00")
	MEM = mem or config.get("MEM_S2", "4G")

	job_name = f"{PROJECT}.s2.{pref}.{family}"
	cmd = [
		"sbatch", f"--qos={QOS_S2}", f"--time={TIME}",
		f"--mem={MEM}",
		f"--job-name={job_name}",
		f"--output={LOG_DIR}/s2/%x.%j.out",
		f"--error={LOG_DIR}/s2/%x.%j.out",
		"workflow/s02_cluster.sh",
		f"{SEARCH_DIR}/{pref}.{family}.domains.fasta",
		f"{CLUSTER_DIR}/{pref}.{family}_cluster.tsv",
		str(MAX_N)
	]

	if after_jid:
		cmd.insert(3, f"--dependency=afterok:{after_jid}")
	if verbose:
		print(' '.join(cmd))
	if not dry:
		out = run_cmd(cmd, dry, verbose)
	else:
		out = None
	jid = out.split()[-1] if out else None
	if not dry and verbose:
		print(f"Submitted {jid}: {job_name}")
	return(jid)


def submit_family(config_fn, pref, family, mode = 'both', time_dict=None, mem_dict = None, dry=False,verbose = True):
	# submit the jobs for families 
	config = parse_bash_config(config_fn)
	#print(f'Submit family mode: {mode}')

	# parse the config file 
	PROJECT = config.get("PROJECT", "project")
	QOS_S1 = config.get("QOS_S1", "normal")
	QOS_S2 = config.get("QOS_S2", "normal")
	LOG_DIR = config.get("LOG_DIR", "logs")
	SEARCH_DIR = config.get("SEARCH_DIR", "results/search")
	CLUSTER_DIR = config.get("CLUSTER_DIR", "results/cluster")
	GENEFAM_INFO = config.get("GENEFAM_INFO")
	INFASTA = config.get("INFASTA")
	MAX_N = config.get("MAX_N", "")

	mem_s1 = mem_s2 = time_s1 = time_s2 = None
	if time_dict:
		time_s1 = time_dict['search']
		time_s2 = time_dict['cluster']

	if mem_dict:
		mem_s1 = mem_dict['search']
		mem_s2 = mem_dict['cluster']
			 

	TIME_S1 = time_s1 or config.get("TIME_S1", "01:00:00")
	TIME_S2 = time_s2 or config.get("TIME_S2", "01:00:00")
	MEM_S1 = mem_s1 or config.get("MEM_S1", "4G")
	MEM_S2 = mem_s2 or config.get("MEM_S2", "8G")

	if mode == 'both':
		# Submit as dependent jobs
		#print(f"Search: {pref}.{family}")
		jid1 = submit_search(pref,family,config,time = TIME_S1,mem = MEM_S1,verbose = verbose, dry = dry)
		#print(f"Submitted {jid1}")
		#print(f"Cluster: {pref}.{family}")
		jid2 = submit_cluster(pref = pref, family = family,config = config,time = TIME_S2,mem = MEM_S2,after_jid = jid1,verbose = False,dry = dry)
		#print(f"Submitted {jid2}")
		return({'search': jid1,'cluster': jid2})
		
	elif mode == 'search':
		print(f"Search: {pref}.{family}")
		jid = submit_search(pref,family,config,time = TIME_S1,mem = MEM_S1,verbose = verbose, dry = dry)
		#print(f"Submitted {jid}")
		return({'search': jid})
	
	elif mode == 'cluster':
		print(f"Cluster: {pref}.{family}")
		jid = submit_cluster(pref = pref, family = family,config = config,after_jid = None,time = TIME_S2,mem = MEM_S2,verbose = verbose,dry = dry)
		#print(f"Submitted {jid}")
		return({'cluster': jid})
		
	else:
		print(f'ERROR: unknown family submit type: {mode}')
		quit()


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
