import os, sys, argparse, subprocess, json
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


def parse_genefam(fn,append_prefix = True):
    genefam = {}
    with open(fn) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) == 7:
                family = fields[0]
                hmms = fields[1].split(',')
                inflation = float(fields[2])
                min_seq = int(fields[3])
                threshold = str(fields[4])
                group = str(fields[5])
                prefix = str(fields[6])
                if append_prefix:
                    key = prefix + '.' + family
                else:
                    key = family
                genefam.update({
                    key: 
                        {
                        'group': group,
                        'prefix': prefix,
                        'family': family,
                        'hmms': hmms,
                        'min_seq' : min_seq,
                        'threshold': threshold
                        }
                        
                    })
            else:
                print(f'Unknown genefam format!')
    return(genefam)


def check_job(jobid):
    def parse_mem(val):
        if val.endswith("K"): return float(val[:-1]) / 1024
        if val.endswith("M"): return float(val[:-1])
        if val.endswith("G"): return float(val[:-1]) * 1024
        if val.endswith("T"): return float(val[:-1]) * 1024 * 1024
        try:
            return float(val)
        except:
            return 0

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