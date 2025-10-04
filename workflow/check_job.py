#!/usr/bin/env python3
import argparse, subprocess, math

def parse_mem(val):
    if val.endswith("K"): return int(val[:-1]) / 1024
    if val.endswith("M"): return int(val[:-1])
    if val.endswith("G"): return int(val[:-1]) * 1024
    if val.endswith("T"): return int(val[:-1]) * 1024 * 1024
    try: return float(val)
    except: return 0

parser = argparse.ArgumentParser(description="Check SLURM job status")
parser.add_argument("jobid", help="SLURM job ID")
args = parser.parse_args()

status = subprocess.run(["squeue", "-j", args.jobid, "-h", "-o", "%T"], capture_output=True, text=True)
if status.stdout.strip():
    state = status.stdout.strip()
    print(f"{args.jobid}\t{state}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN")
else:
    out = subprocess.run(
        ["sacct", "-j", args.jobid, "-o", "JobID,State,ExitCode,Elapsed,Timelimit,ReqMem,MaxRSS", "-n", "-P"],
        capture_output=True, text=True
    )
    lines = out.stdout.splitlines()
    main = None
    max_mem = 0
    for line in lines:
        parts = line.strip().split("|")
        if not parts or len(parts) < 7: continue
        if parts[0] == args.jobid:
            main = parts
        if parts[6].strip() not in ("", "None"):
            mem_val = parse_mem(parts[6].strip())
            if mem_val > max_mem: max_mem = mem_val
    if main:
        main[6] = f"{round(max_mem)}M" if max_mem else "NaN"
        print("\t".join(main))

