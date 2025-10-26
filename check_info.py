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

    try:
        cmd = [
            "sacct", "-j", jobid,
            "-o", "JobID,State,ExitCode,Elapsed,Timelimit,ReqMem,MaxRSS",
            "-n", "-P"
        ]
        out = subprocess.run(cmd, capture_output=True, text=True, check=False)
        lines = [l.strip() for l in out.stdout.splitlines() if l.strip()]
    except Exception:
        lines = []

    main = None
    max_mem = 0

    for line in lines:
        parts = line.split("|")
        if len(parts) < 7:
            continue

        jid = parts[0]
        root = jid.split(".")[0]

        if root == jobid and "." not in jid:
            main = parts  # main job only

        # still track memory from any substep
        if parts[6].strip() not in ("", "None"):
            mem_val = parse_mem(parts[6])
            max_mem = max(max_mem, mem_val)


    if main:
        state = main[1]
        if "CANCEL" in state:
            state = "CANCELLED"
        return {
            "jobid": jobid,
            "status": state, 
            "exitcode": main[2],
            "elapsed": main[3],
            "timelimit": main[4],
            "reqmem": main[5],
            "maxrss": f"{round(max_mem)}M" if max_mem else None
        }

    return {
        "jobid": jobid,
        "status": "UNKNOWN",
        "exitcode": None,
        "elapsed": None,
        "timelimit": None,
        "reqmem": None,
        "maxrss": None
    }

def update_hg(hg,data,verbose =  True):
    if hg in data.keys():
        info = data[hg]
        jobid = info['jobid']
        job_info = check_job(jobid)
        if 'status' in info.keys():
            status_prev = info['status']
        else:
            status_prev = None
        status_new = job_info['status']
        if verbose:
            print(status_prev,status_new) 
        # update the value 
        if job_info["status"] != "UNKNOWN":
            for key,value in job_info.items():
                if key == 'state':
                    key = 'status'
                if key not in info.keys():
                    if verbose:
                        print(f'Adding {key}: {value}')
                else:
                    if job_info[key] != info[key]:
                        if verbose:
                            print(f'Updating {key}: {info[key]} -> {job_info[key]}')
                    if key == 'state':
                        key = 'status'        
                info.update({key:value})
            data.update({hg:info})
            return(status_prev,status_new)
        else:
            print(f'ERROR: unknown job status for {jobid}!')
    else:
        print(f'ERROR: {hg} not found in {json_fn}!')


parser = argparse.ArgumentParser(description = 'Given a JSON file, check and update HG status.')
parser.add_argument('--json', required=True, help='Path to JSON file')
parser.add_argument('--hg', required=False,default = None, help='Optional PREF.FAMILY.HG string')
parser.add_argument('--output', required=False,default = None, help='Optional tab-file to write status for each job')
args = parser.parse_args()

json_fn = args.json
hg = args.hg

out_fn = args.output

update = True
verbose = True

if not os.path.isfile(json_fn):
    print(f"WARNING: JSON file {json_fn} doesn't exist!")
    data = {}
else:
    with open(json_fn) as f:
        data = json.load(f)


if hg:
    print("Provided HG string. Will check the job and update the status")
    #if verbose:
    #    print(f'{hg}')
    if hg in data.keys():
        info = data[hg]
        jobid = info['jobid']
        job_info = check_job(jobid)
        if 'status' in info.keys():
            status_prev = info['status']
        if 'state' in info.keys():
            info.pop('state')
        else:
            status_prev = None
        status_new = job_info['status']
        
        # update the value 
        if job_info["status"] != "UNKNOWN":
            for key,value in job_info.items():
                if key == 'state':
                    key = 'status'
                if key not in info.keys():
                    print(f'Adding {key}: {value}')
                else:
                    if job_info[key] != info[key]:
                        print(f'Updating {key}: {info[key]} -> {job_info[key]}')
                    if key == 'state':
                        key = 'status'        
                info.update({key:value})
            data.update({hg:info})
            if verbose:
                for k,v in data[hg].items():
                    print('\t',k,v)
            if update:
                with open(json_fn, 'w') as f:
                    json.dump(data, f, indent=2)
                if verbose:
                    print(f'Updated {hg} in {json_fn}') 
        else:
            print(f'ERROR: unknown job status for {jobid}!')
    else:
        print(f'WARNING: {hg} not found in {json_fn}!')
    quit()




print(f'{json_fn}: {len([x for x in data.keys()])} entries')


# gather infos 
n_update = 0
infos = {}
states = {}
ntot = len([x for x in data.keys()])
print('Checking JSON entries...\n')
for i,(key, value) in enumerate(data.items()):
    if (i+1)<= ntot:
        if (i+1) % 100 == 0:
            print(f'{i+1}/{ntot}')
        jobid = data[key]['jobid']
        #print(i,key,jobid)
        #print('------')
        old_state, new_state = update_hg(key,data,verbose = False)
        if not old_state and new_state or old_state and old_state != new_state:
            #print("found a new state!")
            states.update({key: new_state})
            n_update += 1
            #print(f"{key}: {jobid}: {old_state} - {new_state}")
        else:
            #print('keeping the old state')
            states.update({key: old_state})
    else:
        break
for key, value in infos.items():
    data[key]['status'] = states[key]


from collections import Counter

counts = Counter(states.values())
print('\nJob status summary:')
for state, n in counts.items():
    print(f"\t{state}: {n}")
print(f"Total: {sum(counts.values())}")



# update the file 
if update:
    with open(json_fn, 'w') as f:
        json.dump(data, f, indent=2)
print(f"Updated {n_update} / {len(data)} entries.")




# Report jobs of a specific status 
def get_jobs(status,states):
    jobs = [k for k,v in states.items() if v == status]
    return(jobs)

if out_fn:
    with open(out_fn, "w") as f:
        for k, v in states.items():
            f.write(f"{k}\t{v}\n")
    print(f'HG states written to {out_fn}')


