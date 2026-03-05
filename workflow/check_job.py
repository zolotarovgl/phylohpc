#!/usr/bin/env python3
import argparse, subprocess, math

def parse_mem(val):
    if val.endswith("K"): return int(val[:-1]) / 1024
    if val.endswith("M"): return int(val[:-1])
    if val.endswith("G"): return int(val[:-1]) * 1024
    if val.endswith("T"): return int(val[:-1]) * 1024 * 1024
    try: return float(val)
    except: return 0

parser = argparse.ArgumentParser(
    description="""Check SLURM job status
JobID    State    ExitCode    Elapsed    Timelimit    ReqMem    MaxRSS""",
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("jobid", nargs="+", help="SLURM job ID(s), space or comma separated")
args = parser.parse_args()

# Flatten comma-separated values
job_ids = []
for item in args.jobid:
    job_ids.extend(item.split(","))
for jobid in job_ids:
    jobid = jobid.strip()
    if not jobid:
        continue

    status = subprocess.run(
        ["squeue", "-j", jobid, "-h", "-o", "%T"],
        capture_output=True, text=True
    )

    if status.stdout.strip():
        state = status.stdout.strip()
        print(f"{jobid}\t{state}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN")
        continue

    out = subprocess.run(
        ["sacct", "-j", jobid,
         "-o", "JobID,State,ExitCode,Elapsed,Timelimit,ReqMem,MaxRSS",
         "-n", "-P"],
        capture_output=True, text=True
    )

    lines = out.stdout.splitlines()
    main = None
    max_mem = 0

    for line in lines:
        parts = line.strip().split("|")
        if not parts or len(parts) < 7:
            continue

        if parts[0] == jobid:
            main = parts

        if parts[6].strip() not in ("", "None"):
            mem_val = parse_mem(parts[6].strip())
            if mem_val > max_mem:
                max_mem = mem_val

    if main:
        main[6] = f"{round(max_mem)}M" if max_mem else "NaN"
        print("\t".join(main))
    else:
        print(f"{jobid}\tUNKNOWN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN")