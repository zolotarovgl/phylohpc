#!/usr/bin/env python3
import os, sys, argparse, subprocess, json
import subprocess
def increase_mem(mem_str,step = 1024):
    # step - increase memory by this amount in Mb
    s = mem_str.strip().upper().rstrip("B")  # handle Mb, MB, etc.
    if s.endswith("G"):
        val = int(s[:-1]) + 1
        return f"{val}G"
    elif s.endswith("M"):
        val = int(s[:-1]) + step
        return f"{val}M"
    else:
        # assume in MB if no suffix
        val = int(s) + step
        return f"{val}M"


def check_job(jobid):
    def parse_mem(val):
        if val.endswith("K"): return int(val[:-1]) / 1024
        if val.endswith("M"): return int(val[:-1])
        if val.endswith("G"): return int(val[:-1]) * 1024
        if val.endswith("T"): return int(val[:-1]) * 1024 * 1024
        try: return float(val)
        except: return 0

    status = subprocess.run(["squeue","-j",jobid,"-h","-o","%T"],capture_output=True,text=True)
    if status.stdout.strip():
        state = status.stdout.strip()
        return {
            "jobid": jobid,
            "state": state,
            "exitcode": None,
            "elapsed": None,
            "timelimit": None,
            "reqmem": None,
            "maxrss": None
        }
    else:
        out = subprocess.run(
            ["sacct","-j",jobid,"-o","JobID,State,ExitCode,Elapsed,Timelimit,ReqMem,MaxRSS","-n","-P"],
            capture_output=True,text=True
        )
        lines = out.stdout.splitlines()
        main = None
        max_mem = 0
        for line in lines:
            parts = line.strip().split("|")
            if not parts or len(parts) < 7: continue
            if parts[0] == jobid:
                main = parts
            if parts[6].strip() not in ("","None"):
                mem_val = parse_mem(parts[6].strip())
                if mem_val > max_mem: max_mem = mem_val
        if main:
            return {
                "jobid": main[0],
                "state": main[1],
                "exitcode": main[2],
                "elapsed": main[3],
                "timelimit": main[4],
                "reqmem": main[5],
                "maxrss": f"{round(max_mem)}M" if max_mem else None
            }
    return {
        "jobid": jobid,
        "state": "UNKNOWN",
        "exitcode": None,
        "elapsed": None,
        "timelimit": None,
        "reqmem": None,
        "maxrss": None
    }



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
    parser = argparse.ArgumentParser(description="A python wrapper to submit HG jobs")
    parser.add_argument("configfile", help="Bash-style config file with KEY=VALUE lines")
    parser.add_argument("--pref", required=True, help="Prefix (PREF)")
    parser.add_argument("--family", required=True, help="Gene family (FAMILY)")
    parser.add_argument("--hg", required=True, help="HG identifier")
    parser.add_argument("--mode", choices=["align","phylogeny","possvm","generax"], required=True, help="Which step to run")
    parser.add_argument("--cpus", type=int, default=4, help="CPUs per task")
    parser.add_argument("--mem", default="1G", help="Memory")
    parser.add_argument("--mem_increase", default=10240, help="Memory increase step, Mb. Default: 10240")
    parser.add_argument("--time", default="1:00:00", help="Walltime")
    parser.add_argument("-n", "--dry_run", action="store_true", help="Dry run: only print the sbatch command")
    parser.add_argument("--mafft", default="--maxiterate 1000 --localpair", help="MAFFT settings")
    parser.add_argument("--json", help="Optional JSON file to store job IDs")
    args = parser.parse_args()

    config = parse_bash_config(args.configfile)
    ALIGN_DIR = config["ALIGN_DIR"]
    TREE_DIR = config["TREE_DIR"]
    TREE_METHOD_BIG = config.get("TREE_METHOD_BIG","iqtree")
    REFSPECIES = config.get("REFSPECIES")
    REFNAMES = config.get("REFNAMES")
    SPECIES_TREE = config.get("SPECIES_TREE")
    
    PREF, FAMILY, HG = args.pref, args.family, args.hg
    jobname = ".".join([args.mode,PREF,FAMILY,HG])

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
    elif args.mode == "generax":
        time = config.get("TIME_S5")
        mem = config.get("MEM_S5")
        print(f'Generax:\nTime: {time}\nMem: {mem}')
        args.time = time
        args.mem = mem
        if not os.path.isfile(SPECIES_TREE):
            sys.exit(f"Error: species tree not found! {SPECIES_TREE}")
        out_tree = f"results/generax/{PREF}.{FAMILY}.{HG}.treefile"
        aln_file = f"{ALIGN_DIR}/{PREF}.{FAMILY}.{HG}.aln.fasta"
        gene_tree = f"{TREE_DIR}/{PREF}.{FAMILY}.{HG}.treefile"
        output_dir = f"results/generax/{PREF}.{FAMILY}.{HG}"
        wrap = f"python phylogeny/main.py generax --alignment {aln_file} --gene_tree {gene_tree} --species_tree {SPECIES_TREE} --output_dir {output_dir} --subs_model LG -c 1 --name {PREF}.{FAMILY}.{HG} --outfile {out_tree} -c {args.cpus}"
   
    # check if the pref.fam.hg has been ran before and increase the resources
    if args.json:
        if os.path.exists(args.json):
            with open(args.json) as jf:
                data = json.load(jf)
        else:
            data = {}
        key = f"{PREF}.{FAMILY}.{HG}"
        if key in data:
            # job handling rules 
            #print(f"WARNING: Entry for {key} already exists in {args.json}: ")
            #print(f"{args.json}: ")
            #print(data[key])
            #print("Job info:")
            jobid = data[key]['jobid']
            info = check_job(jobid)
            print(info)
            job_state = info.get("state")
            if "CANCELLED" in job_state:
                job_state = "CANCELLED"
            print(f"{key} {jobid} {job_state}")
            # check if job died due to out of memory
            if info.get("exitcode") == "0:125" or info.get("state") == "OUT_OF_MEMORY":
                prev_mem = data[key]['mem']
                new_mem = increase_mem(data[key].get("mem", args.mem),int(args.mem_increase))
                print(f"OOM, increasing memory: {prev_mem}  -> {new_mem}")
                args.mem = new_mem
            elif job_state == "CANCELLED":
                print(f"Job was cancelled -> re-submitting")
            elif job_state == "FAILED":
                print(f"Job has failed -> re-submitted")
            elif job_state == "TIMEOUT":
                print(info)
                prev_time = info['elapsed']
                new_time = args.time
                new_time = "7-00:00:00"
                args.time = new_time
                print(f"Timeout, increasing time: {prev_time} -> {new_time} ")
            else:
                sys.exit()

    cmd = ["sbatch", f"--job-name={jobname}", f"--cpus-per-task={args.cpus}", f"--mem={args.mem}", f"--time={args.time}", "--wrap", wrap]
    #print(" ".join(cmd))
    
    if args.dry_run:
        print('dry run is set!')
        print(args.dry_run)
    else:
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            sys.exit(f"sbatch failed:\n{res.stderr}")
        output = res.stdout.strip()
        jobid = None
        for token in output.split():
            if token.isdigit():
                jobid = token
                break
        if not jobid:
            sys.exit(f"Could not parse job ID from sbatch output:\n{output}")

    # Update the json
    if args.json and not args.dry_run:
        data[key] = {
            "jobid": jobid,
            "mode": args.mode,
            "command": wrap,
            "cpus": args.cpus,
            "mem": args.mem,
            "time": args.time
        }
        with open(args.json,"w") as jf:
            json.dump(data, jf, indent=2)


    print(f"Submitted batch job {jobid}")

if __name__ == "__main__":
    main()

