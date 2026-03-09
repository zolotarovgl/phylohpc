#!/usr/bin/env python3
import argparse
import subprocess
import sys


def parse_mem(val):
    try:
        if val.endswith("K"): return float(val[:-1]) / 1024
        if val.endswith("M"): return float(val[:-1])
        if val.endswith("G"): return float(val[:-1]) * 1024
        if val.endswith("T"): return float(val[:-1]) * 1024 * 1024
        return float(val)
    except:
        return 0


parser = argparse.ArgumentParser(
    description="""Check SLURM job status
JobID    State    ExitCode    Elapsed    Timelimit    ReqMem    MaxRSS""",
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument(
    "jobid",
    nargs="*",
    help="SLURM job ID(s), space or comma separated"
)

parser.add_argument(
    "-f", "--file",
    help="Text file containing job IDs (one per line)"
)

args = parser.parse_args()

# -----------------------------
# Collect job IDs
# -----------------------------

job_ids = []

# from command line
for item in args.jobid:
    job_ids.extend(item.split(","))

# from file
if args.file:
    try:
        with open(args.file) as f:
            for line in f:
                job_ids.extend(line.strip().split(","))
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

# clean list
job_ids = [j.strip() for j in job_ids if j.strip()]

if not job_ids:
    print("No job IDs provided.", file=sys.stderr)
    sys.exit(1)

job_str = ",".join(job_ids)

# -----------------------------
# Check running jobs (squeue)
# -----------------------------

running = {}

sq = subprocess.run(
    ["squeue", "-j", job_str, "-h", "-o", "%i|%T"],
    capture_output=True,
    text=True
)

for line in (sq.stdout or "").splitlines():
    try:
        jid, state = line.split("|")
        running[jid] = state
    except:
        continue

# -----------------------------
# Query sacct once
# -----------------------------

out = subprocess.run(
    [
        "sacct",
        "-j", job_str,
        "-o", "JobID,State,ExitCode,Elapsed,Timelimit,ReqMem,MaxRSS",
        "-n",
        "-P"
    ],
    capture_output=True,
    text=True
)

if out.returncode != 0:
    print("Error running sacct:", out.stderr, file=sys.stderr)
    sys.exit(1)

records = {}
maxrss = {}

for line in (out.stdout or "").splitlines():

    parts = line.strip().split("|")
    if len(parts) < 7:
        continue

    jid = parts[0].split(".")[0]

    if jid not in records:
        records[jid] = parts

    mem = parts[6].strip()

    if mem not in ("", "None"):
        val = parse_mem(mem)
        if jid not in maxrss or val > maxrss[jid]:
            maxrss[jid] = val

# -----------------------------
# Print results
# -----------------------------

for jid in job_ids:

    # running jobs
    if jid in running:
        print(f"{jid}\t{running[jid]}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN")
        continue

    # finished jobs
    if jid in records:
        rec = records[jid]
        mem = round(maxrss.get(jid, 0))
        rec[6] = f"{mem}M" if mem else "NaN"
        print("\t".join(rec))
    else:
        print(f"{jid}\tUNKNOWN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN")
