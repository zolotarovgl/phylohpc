# This script will check the families 
import subprocess
import os
import json
import argparse 
import glob


from workflow.helper import parse_bash_config
from workflow.helper import parse_genefam
from workflow.helper import check_job
from workflow.helper import read_json

def count_lines(path):
	with open(path) as f:
		return sum(1 for _ in f)


def get_family_job_info(fam,json_info,job_name,verbose = True):
    # output: dict job_id, mem, time 
    out = {}
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
            else:
                if verbose:
                    print(f'WARNING: No job id found')
        else:
            if verbose:
                print(f'WARNING: job info not found!')
    return(out)

##################
# Check family

# 1. search
# f"{config['SEARCH_DIR']}/*.domains.fasta"
# 2. clustering 
# f"{config['CLUSTER_DIR']}/*_cluster.tsv"
# Explicitly define which files to check 
def check_outfile(family,config,job_type = None):
    if job_type == 'search':
        SEARCH_DIR = config['SEARCH_DIR']
        file_glob = f"{config['SEARCH_DIR']}/*.domains.fasta"
        files = glob.glob(file_glob)
        toreplace = file_glob.split('*')[1]
        fams = [".".join(os.path.basename(f).replace(toreplace,'').split('.')[:2]) for f in files]
        state = int(family in fams)
        return(state)
    elif job_type == 'cluster':
        CLUSTER_DIR = config['CLUSTER_DIR']
        file_glob = f"{config['CLUSTER_DIR']}/*_cluster.tsv"
        files = glob.glob(file_glob)
        toreplace = file_glob.split('*')[1]
        fams = [".".join(os.path.basename(f).replace(toreplace,'').split('.')[:2]) for f in files]
        state = int(family in fams)
        return(state)
    else:
        print(f'ERROR: Unknown job type: {job_type}')
        quit()



def check_families(family = None,json_info = None,json_fn = None,genefam_dict = None,config = None,update_json = True,verbose = True):
    update_counter = 0
    new_counter = 0

    if family:
        fams = family.split(',')
    else:
        # If the fams haven't been provided, check all from the gene family file
        if verbose:
            print('No family provided, checking all of them in the genefam') 
        fams = [family for family in genefam_dict.keys()]

    # Check if the families are defined 
    for family in fams:
        if not family in genefam_dict.keys():
            print(f'ERROR: family "{family}" is not present in {g_fn}')
            quit()

    # check output files
    for family in fams:        
        status = {'search': None, 'cluster': None}
        # 1. search output 
        status['search'] = check_outfile(family = family,config = config,job_type = 'search')
        # 2. clustering output 
        status['cluster'] = check_outfile(family = family,config = config,job_type = 'cluster')
        if not family in json_info.keys():
            status = {'search': {'status': status['search']}, 'cluster' : {'status': status['cluster']}}
            json_info.update({family: status})
            new_counter += 1
        else:
            for k,v in status.items():
                json_info[family][k]['status'] = v
            
                # check the job status 
                if 'job' in json_info[family][k].keys():
                    jobid = json_info[family][k]['job']['jobid']
                    if 'state' in json_info[family][k].keys():
                        state_prev = json_info[family][k]['state']
                    else:
                        state_prev = None
                    json_info[family][k]['job'] = check_job(jobid)
                    state_new = json_info[family][k]['job']['state']  
                    if state_new != state_prev:
                        update_counter += 1
    if verbose:
        print(f'Added {new_counter} new familes')
        print(f'Updated {update_counter} job states')
    
    if update_json and json_fn:
        with open(json_fn, 'w') as f:
            json.dump(json_info, f, indent=2)
        if verbose:
            print(f'Updated: {json_fn}')

    return(json_info)

########################################
# Search job results   
########################################
parser = argparse.ArgumentParser(description="Check the status of families")
parser.add_argument("configfile", help="Bash-style config file with KEY=VALUE lines")
parser.add_argument("--genefam", required=True, help="Path to gene family definition file")
parser.add_argument("--json", required = True, help="JSON file to store the info")
parser.add_argument("--output",default = None, help="Optional output tab file with family states")
parser.add_argument("--family",default = None, help="Optional family name(s) to check. Comma-separated")
parser.add_argument("--stdout", action="store_true", help="Print family statuses to stdout")
parser.add_argument("--resubmit", action="store_true", help="Whether to re-submit failed family jobs")
parser.add_argument("--dry", action="store_true", help="Re-submission dry mode: print the commands and exit")
parser.add_argument("--max_iter", default = int(3), help="Maximum number of times to try to submit a job")
parser.add_argument("--ignore-iter", action="store_true", help="Use to ignore the maximum number of iterations")
args = parser.parse_args()


# Given a list of families or the whole genefam file, check their status

g_fn = args.genefam
c_fn = args.configfile
json_fn = args.json
out_fn = args.output
family = args.family 
stdout = args.stdout
dry = args.dry
max_iter = int(args.max_iter) # maximum number of times to try to submit a job 
if args.ignore_iter:
    max_iter = 100000

from workflow.submit_family import submit_family
resubmit = args.resubmit

update_json = True
verbose = True

json_info = read_json(json_fn)
genefam_dict = parse_genefam(g_fn)
config = parse_bash_config(c_fn)

# checks the status of the families by checking the output files, also update the job states 
json_info = check_families(family = args.family,json_info = json_info,json_fn = json_fn,genefam_dict = genefam_dict,config = config,update_json = update_json,verbose = verbose)


states = {}
for k,v in json_info.items():
    if k in genefam_dict.keys():
        states.update({k : "".join([str(v['search']['status']),str(v['cluster']['status'])])})


if out_fn:
    with open(out_fn, "w") as f:
        for k, v in states.items():
            f.write(f"{k}\t{v}\n")
    if not stdout and verbose:
        print(f'Family states written to: {out_fn}')

# Report job states
from collections import Counter
counts = Counter([x for x in [v for v in states.values()]])
if not stdout:
    print(f'Family status [search,clustering]:')
    for state, n in counts.items():    
        print(f"-  {state}: {n}")
    print(f"Total: {sum(counts.values())}")

if stdout:
    for k,v in states.items():
        print(f'{k}\t{v}')

################################################################################
# Resubmit jobs
# Add an increase  in memory!
################################################################################
# the submit family will submit 2 jobs at the same time 

if resubmit:
    print("-------------------------------")
    print("""# Family jobs re-submission""")
    print("-------------------------------")

    # just make the dictionary with resourcnes
    res_dict = {'search' : {'mem' : config['MEM_S1'], 'time' : config['TIME_S1']},'cluster' : {'mem' : config['MEM_S2'], 'time' : config['TIME_S2']}}
    
    fams = [k for k,v in states.items() if v in ['00','01','10']]

    if args.family:
        fams = [x for x in fams if x in args.family.split(',')]
    if len(fams)>0:
        print(f'{len(fams)} families to re-submit')
    else:
        print(f'No families to re-submit!')
        quit()


    # the function to pull the family job states
    def check_job_states(fam,json_info):
        job_states = {}
        for job_type, v  in json_info[fam].items():
            job_state = None
            if 'job' in v.keys():
                if 'jobid' in v['job'].keys():
                    job_state = check_job(v['job']['jobid'])['state']
                
            job_states.update({job_type : job_state})
        return(job_states)

    def check_file_states(fam,json_info):
        file_states = {}
        for job_type, v  in json_info[fam].items():
            file_state = None
            if 'status' in v.keys():
                file_state = v['status']
            else:
                file_state = 0
            file_states.update({job_type : file_state})
        return(file_states)

    submit_counter = 0

    for fam in fams:
        # What is the best logic to re-submit the family?
        # the file states are stored somewhere
        # which script checks for the file states? 

        pref = fam.split(".")[0]
        family = fam.split(".")[1]
        state = str(states[fam])
        mode = None
        job_states = check_job_states(fam,json_info) # gets the state of the last job 
        file_states = check_file_states(fam,json_info) # gets the file status of the last job 

       
        pref = fam.split(".")[0]
        family = fam.split(".")[1]
        state = str(states[fam])
        mode = None
        job_states = check_job_states(fam,json_info)
        file_states = check_file_states(fam,json_info)
        # the file state should be always above the job states 
        d = {}
        for k,v in file_states.items():
            d.update({k : {'job' : job_states[k], 'file' : file_states[k]}})       
        
        # decide on a final state (what to do):
        dd = {}
        increase_mem_dict = {}
        increase_time_dict = {}
        for k,v in d.items():
            final_state = 0
            increase_mem = False
            increase_time = False
            job_state = v['job']
            file_state = v['file']
            if file_state == 1:
                final_state = 1
            else:
                final_state = 0
                if job_state in ["RUNNING","PENDING"]:
                    js = check_job(v['job'])
                    print(f'WARNING: {fam:} {k}: job state is {job_state} - Skipping')
                    final_state = 1
            if False:
                if not job_state:
                    final_state = 0
                elif "CANCEL" in job_state:
                    job_state = "CANCELLED"
                    final_state = 0
                elif job_state in ["RUNNING","PENDING"]:
                    js = check_job(v['job'])
                    print(f'WARNING: {fam:} {k}: job state is {job_state} - Skipping')
                    final_state = 1
                elif job_state in ["OUT_OF_MEMORY"]:
                    # increase the memory!
                    print(f'OOM: Increasing memory')
                    increase_mem = True
                    final_state = 0
                elif job_state in ["TIMEOUT"]:
                    # increase the memory!
                    print(f'TIMEOUT: Increasing time')
                    increase_time = True
                    final_state = 0
                elif job_state in [ "COMPLETING",'COMPLETED']:
                    final_state = 1
                elif job_state in ["FAILED","CANCELLED"]:
                    final_state = 0
                else:
                    print(f'ERROR: what to do with this state? {job_state}')
                    quit()
            
            dd.update({k : str(final_state)})
            increase_mem_dict.update({k : increase_mem})
            increase_time_dict.update({k : increase_time})

        state = "".join([str(x) for x in dd.values()])
        
        if state == '00':
            mode = 'both'
        elif state == '10':
            mode = 'cluster'
        elif state == '01':
            mode = 'search'
        elif state == '11':
            mode = 'DONE'
        else:
            print(f'Unknown state: {state}')
            quit()
        if mode != 'DONE':
            # Check if the number of iterations has been exceeded:
            n_iter = {}
            if 'iter' in json_info[fam]['search'].keys():
                search_iter = json_info[fam]['search']['iter']
            else:
                search_iter = 0
            
            if 'iter' in json_info[fam]['cluster'].keys():
                cluster_iter = json_info[fam]['cluster']['iter']
            else: 
                cluster_iter  = 0
            iters = [search_iter,cluster_iter]
            if max(iters) >= max_iter:
                print("-----------------------------------")         
                print(f'{fam}: ERROR: Can not submit the job: maximum number of iterations reached (iter: {max(iters)}, max_iter: {max_iter})\nExplore the log files to figure out why the job has been failing!\n(Disable using --ignore-iter)')
                print("-----------------------------------")    
            else:
            
                for job_type in res_dict.keys():
                    if increase_mem_dict[job_type]:
                        if json_info[fam][job_type]['job']['reqmem']:
                            mem_prev = json_info[fam][job_type]['job']['reqmem']
                        else:
                            mem_prev = '100M'
                        mem_new = increase_mem(mem_prev)
                        res_dict[job_type].update(mem = mem_new)
                    if increase_time_dict[job_type]:
                        time_step = '01:00:00'
                        def increase_time(time,time_step):
                            from datetime import timedelta 
                            def increase_time(base, step):
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
                        
                        prev_time = json_info[fam][job_type]['job']['timelimit']
                        new_time = increase_time(prev_time,time_step)
                        #print(f"{fam}: Timeout, increasing time: {prev_time} by {time_step} -> {new_time} ")
                        res_dict[job_type].update(mem = new_time)


                mem_dict = {k : v['mem'] for k,v in res_dict.items()}
                time_dict = {k : v['time'] for k,v in res_dict.items()}
                
                job_ids = submit_family(
                    config_fn=c_fn,
                    pref=pref,
                    family=family,
                    mode=mode,
                    mem_dict = mem_dict,
                    time_dict = time_dict,
                    dry=dry,
                    verbose=verbose
                )

                submit_counter += 1  

                if mode == 'both':
                    search_iter += 1
                    cluster_iter += 1
                elif mode == 'search':
                    search_iter += 1
                elif mode == 'cluster':
                    cluster_iter += 1

                iters = {'search' : search_iter, 'cluster' : cluster_iter}
                
                for k,v in job_ids.items():
                    job_state = check_job(v)
                    json_info[fam][k].update({'iter' : iters[k] })
                    json_info[str(fam)][k]['job'] = check_job(v)
                    
                        

    print(f'Submitted {submit_counter} new jobs!')

    if update_json and json_fn:
        with open(json_fn, 'w') as f:
            json.dump(json_info, f, indent=2)
        print(f'Updated: {json_fn}')
    else:
        print("Skipping job re-submission")
        quit()



