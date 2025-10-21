import os, sys, argparse, subprocess, json
import subprocess
import json

def check_job(jobid):
    def parse_mem(val):
        try:
            val = val.strip()
            if val.endswith("K"): return int(val[:-1]) / 1024
            if val.endswith("M"): return int(val[:-1])
            if val.endswith("G"): return int(val[:-1]) * 1024
            if val.endswith("T"): return int(val[:-1]) * 1024 * 1024
            return float(val)
        except Exception:
            return 0

    # try squeue first
    try:
        status = subprocess.run(
            ["squeue", "-j", jobid, "-h", "-o", "%T"],
            capture_output=True,
            text=True,
            check=False
        )
    except Exception as e:
        status = None

    state = status.stdout.strip() if status and status.stdout else ""
    if state:
        return {
            "jobid": jobid,
            "state": state,
            "exitcode": None,
            "elapsed": None,
            "timelimit": None,
            "reqmem": None,
            "maxrss": None
        }

    # fallback to sacct if squeue gave nothing
    try:
        out = subprocess.run(
            ["sacct", "-j", jobid, "-o", "JobID,State,ExitCode,Elapsed,Timelimit,ReqMem,MaxRSS", "-n", "-P"],
            capture_output=True,
            text=True,
            check=False
        )
        lines = [l.strip() for l in out.stdout.splitlines() if l.strip()]
    except Exception as e:
        lines = []

    main = None
    max_mem = 0

    for line in lines:
        parts = line.split("|")
        if len(parts) < 7:
            continue

        job_main = parts[0].split(".")[0]
        if job_main != jobid:
            continue  # skip batch/extern

        # record the main entry
        main = parts

        # track MaxRSS from all parts (to include batch steps)
        if parts[6].strip() not in ("", "None"):
            mem_val = parse_mem(parts[6])
            max_mem = max(max_mem, mem_val)

    if main:
        return {
            "jobid": main[0].replace('.extern',''),
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


parser = argparse.ArgumentParser()
parser.add_argument('--json', required=True, help='Path to JSON file')
args = parser.parse_args()

json_fn = args.json
update = False
update = True

with open(json_fn) as f:
    data = json.load(f)

# gather infos 
n_update = 0
infos = {}
states = {}
for key, value in data.items():
    jobid = value['jobid']
    info = check_job(jobid)
    new_state = info['state']
    infos[key] = info
    print(f"{key}\t{new_state}")
    if 'status' in data[key]:
        old_state = data[key]['status']
        if old_state == "UNKNOWN" and old_state != new_state:
            print(f"Updated {key}: {old_state} -> {new_state}")
    else:
        old_state = None
    # decide how to update the status 
    if not old_state and new_state or old_state and old_state != new_state:
        #print("found a new state!")
        states.update({key: new_state})
        n_update += 1
        print(f"{key}: {old_state} - {new_state}")
    else:
        #print('keeping the old state')
        states.update({key: old_state})
for key, value in infos.items():
    data[key]['status'] = states[key]


print(f"Updated {n_update} / {len(data)} entries.")




# update the file 
if update:
    with open(json_fn, 'w') as f:
        json.dump(data, f, indent=2)
