import subprocess
import os
import json
import argparse 
import glob
import sys


from workflow.helper import parse_bash_config
from workflow.helper import parse_genefam
from workflow.helper import check_job
from workflow.helper import read_json

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


def increase_time(base, step):
    from datetime import timedelta 
    def parse(t):
        if '-' in t:
            d, hms = t.split('-')
            days = int(d)
        else:
            days = 0
            hms = t
        h, m, s = map(int, hms.split(':'))
        return timedelta(days=days, hours=h, minutes=m, seconds=s)
    def fmt(td):
        d = td.days
        h, rem = divmod(td.seconds, 3600)
        m, s = divmod(rem, 60)
        return f"{d}-{h:02}:{m:02}:{s:02}" if d else f"{h:02}:{m:02}:{s:02}"
    return fmt(parse(base) + parse(step))


def convert_state(job_state):
    if not job_state:
        return("None")
    elif "CANCEL" in job_state:
        return("CANCELLED")
    else:
        return(job_state)


def get_family_job_info(fam,json_info,job_name = None,verbose = True):
    # output: dict job_id, mem, time 
    out = {}
    if not fam in json_info.keys():
        if verbose:
            print(f'{fam}: WARNING not found in JSON')
        return(out)
    jobs = [k for k in json_info[fam].keys()]
    if verbose:
        print(f'{fam}: {",".join(jobs)}')
    if not job_name in jobs:
        if verbose:
            print(f'{fam}: WARNING: job {job_name} not found in the job log! Available jobs: {",".join(jobs)}')
    else:
        info = json_info[fam][job_name]
        if 'job' in info.keys():
            if 'jobid' in info['job'].keys():
                jobid = info['job']['jobid']
                if verbose:
                    print(f'Found job_id ({jobid}). Checking the info ...')
                check_out = check_job(jobid)
                out.update({'jobid' : jobid})
                if 'timelimit' in check_out.keys():
                    out.update({'time' : check_out['timelimit']})
                if 'reqmem' in check_out.keys():
                    out.update({'mem' : check_out['reqmem']})
                if 'state' in check_out.keys():
                    job_state = check_out['state']
                    job_state = convert_state(job_state)
                    out.update({'state' : job_state})
            else:
                if verbose:
                    print(f'WARNING: No job id found')
        else:
            if verbose:
                print(f'WARNING: job info not found!')
    return(out)

def check_outfile(family,config,job_type = None):
    if job_type == 'search':
        SEARCH_DIR = config['SEARCH_DIR']
        file_glob = f"{SEARCH_DIR}/*.domains.fasta"
        files = glob.glob(file_glob)
        toreplace = file_glob.split('*')[1]
        fams = [".".join(os.path.basename(f).replace(toreplace,'').split('.')[:2]) for f in files]
        state = int(family in fams)
        return(state)
    elif job_type == 'cluster':
        CLUSTER_DIR = config['CLUSTER_DIR']
        file_glob = f"{CLUSTER_DIR}/*_cluster.tsv"
        files = glob.glob(file_glob)
        toreplace = file_glob.split('*')[1]
        fams = [".".join(os.path.basename(f).replace(toreplace,'').split('.')[:2]) for f in files]
        state = int(family in fams)
        return(state)
    else:
        print(f'ERROR: Unknown job type: {job_type}')
        quit()

def decide(job_state,file_state,verbose = True):
    # This functoin will decide which job to run for a particular job 
    # decides whether to run a particular job based on the file state and the job status
    # Output: job_state, file_state, dict torun, increase_mem, increase time
    increase_time = False
    increase_mem = False
    if file_state == 1:
        if verbose:
            print(f'Output file found! Not running the job!')
        torun = 0
    else:
        if not job_state:
            job_state = "None"
            torun = 1
        elif 'CANCEL' in job_state:
            if verbose:
                print('Job was cancelled -> rerunning ...')
            job_state = "CANCELLED"
            torun = 1
        elif job_state in ["RUNNING","PENDING"]:
            if verbose:
                print(f'WARNING: {fam:} {k}: job state is {job_state} - Skipping')
            torun = 0
        elif job_state in ["OUT_OF_MEMORY"]:
            # increase the memory!
            if verbose:
                print(f'OOM: Increasing memory')
            increase_mem = True
            torun = 1
        elif job_state in ["TIMEOUT"]:
            if verbose:
                print(f'TIMEOUT: Increasing time')
            increase_time = True
            torun = 1
        elif job_state in [ "COMPLETING",'COMPLETED']:
            if verbose:
                print(f'Completed job')
            torun = 0
        elif job_state in ["FAILED","CANCELLED"]:
            if verbose:
                print(f'Cancelled job')
            torun = 1
        else:
            print(f'ERROR: unsure what to do with this state: {job_state}')
            quit()
    if verbose:
        print('job_state\tfile_state\ttorun\tincrease_time\tincrease_mem')
        print(f'{job_state}\t{file_state}\t{torun}\t{increase_time}\t{increase_mem}')
    return({'torun' : torun, 'increase_mem' : increase_mem, 'increase_time' : increase_time})


def get_default_res(config,job_name,verbose = True):
    # get the default memory and time settings from the config
    default_res = {}
    config_dict = {'search' : {'mem' : "MEM_S1",'time': "TIME_S1"},'cluster': {'mem' : "MEM_S2",'time': "TIME_S2"}}
    for job_name,params in config_dict.items():
        job_res = {}
        for k,v in params.items():
            if v in config.keys():
                job_res.update({k:config[v]})
        default_res.update({job_name : job_res})
    return(default_res)
        

def get_job_res(fam,json_info,config,job_name,add_mem = False,add_time = False,time_step = '1:00:00',mem_step = 100,verbose = True):
    # collect the default and the last resources per job, increase if needed 
    job_info = get_family_job_info(fam = fam,json_info = json_info,job_name = job_name,verbose = verbose)
    
    res = get_default_res(config,job_name,verbose = verbose)
    # Decide on resources
    if not job_info:
        return(res[job_name])
    if 'mem' in job_info.keys():
        if verbose:
            print(f'{fam}: {job_name}: Picking old memory from JSON')
        mem = job_info['mem']
    else:
        if verbose:
            print(f'{fam}: {job_name}: Picking from config')
        mem = res[job_name]['mem']
    if 'time' in job_info.keys():
        if verbose:
            print(f'{fam}: {job_name}: Picking old time from JSON')
        time = job_info['time']
    else:
        if verbose:
            print('Picking from config')
        time = res[job_name]['time']
    if add_mem:
        mem_old = mem
        mem_new = increase_mem(mem_old,mem_step)
        if verbose:
            print(f'{fam}: {job_name}: Increasing memory by {mem_step}: {mem_old} -> {mem_new}')
        mem = mem_new
    if add_time:
        time_old = time
        time_new = increase_time(time_old,time_step)
        if verbose:
            print(f'{fam}: {job_name}: Increasing memory by {time_step}: {time_old} -> {time_new}')
        time = time_new
    res = {'time' : time, 'mem' : mem}
    return(res)

def check_family(fam,json_info,config,time_step_dict,mem_step_dict,verbose = True):
    # The list of available jobs, get their execution params and decide on whether to relaunch them
    output = {} 
    jobs = ['search','cluster']

    if not fam in json_info.keys():
        if verbose:
            print(f'WARNING: {fam} not found in json_info!')
        for job_name in jobs:
            res = get_job_res(fam = fam,json_info = json_info,config = config,job_name = job_name,add_mem = False,add_time = False,time_step = None,mem_step = None,verbose = verbose)
            output.update({job_name : {'res' : res}})
    else:
        for job_name in jobs:
            file_state = check_outfile(fam,config,job_type = job_name)
            job_info = get_family_job_info(fam = fam,json_info = json_info,job_name = job_name,verbose = verbose)
            if job_info:
                job_state = job_info['state']
            else:
                job_state = None
            decision = decide(file_state = file_state,job_state = job_state)
            print(decision)
            time_step = time_step_dict[job_name]
            mem_step = mem_step_dict[job_name]
            if decision['torun']:
                res = get_job_res(fam = fam,json_info = json_info,config = config,job_name = job_name,add_mem = decision['increase_mem'],add_time = decision['increase_time'],time_step = time_step,mem_step = mem_step,verbose = verbose)
                output.update({job_name : {'res' : res}})

    return(output)

#########################################################################
# A simple script to deal with a single family, trigger re-submission etc
########################################################################
# based on the job_state and a file state, decide whether to relaunch the jobs 
# what am I doing here? 
verbose = True
json_fn = 'info/families.json'
config_fn = 'configs/config.txt'
# Load family info json
json_info = read_json(json_fn)
config = parse_bash_config(config_fn)

mem_step_dict = {'search' : 100,'cluster' : 1024}
time_step_dict = {'search' : '00:01:00','cluster' : '00:01:00'}

resubmit = False
dry = True
from workflow.submit_family import submit_family
fam = 'tfs.Filamin'
pref,family = fam.split('.')
res_dict = check_family(fam = fam,json_info = json_info,config = config,time_step_dict = time_step_dict,mem_step_dict = mem_step_dict,verbose = verbose)

if(len(res_dict)>0):
    submit_ids = {}
    print(f'{len(res_dict)} jobs to run')
    for job_name in res_dict.keys():
    
        pref = 'tfs'
        family = 'Filamin'
        mode = 'search'
        # wtf is this really? 

        job_ids = submit_family(
            config_fn=config_fn,
            pref=pref,
            family=family,
            mode=mode,
            time_dict = None,
            mem_dict = None,
            dry=dry,
            verbose=verbose
        )
        # gather the info on the submitted jobs:
        
        quit()
        submit_ids.update({job_name : job_ids})
    print(submit_ids)

quit()
# what the fuck grisha. really. 

