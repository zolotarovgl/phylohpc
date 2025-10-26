# This script will check the families 
import subprocess
import os
import json
import argparse 
import glob


from workflow.helper import parse_bash_config
from workflow.helper import parse_genefam
from workflow.helper import check_job

def count_lines(path):
	with open(path) as f:
		return sum(1 for _ in f)


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
        # update job info if present 
        #if 'job' in json_info['family']
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
parser.add_argument("--max_iter", default = int(2), help="Maximum number of times to try to submit a job")
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


update_json = True
verbose = True

json_info = read_json(json_fn)
genefam_dict = parse_genefam(g_fn)
config = parse_bash_config(c_fn)

# checks the status of the families by checking the output files, also update the job states 
json_info = check_families(family = args.family,json_info = json_info,json_fn = json_fn,genefam_dict = genefam_dict,config = config,update_json = update_json,verbose = verbose)


states = {}
for k,v in json_info.items():
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


########################################
# Resubmit jobs
# Make sure to always update json file 
# do not re-submit the jobs that are running 
########################################
# the submit family will submit 2 jobs at the same time 
from workflow.submit_family import submit_family
resubmit = args.resubmit


if resubmit:
    print("""# Job re-submission""")
    fams = [k for k,v in states.items() if v in ['00','01','10']]

    if len(fams)>0:
        print(f'{len(fams)} families to re-submit')
    else:
        print(f'No families to re-submit!')
        quit()



    def check_job_states(fam,json_info):
        job_states = {}
        for job_type, v  in json_info[fam].items():
            job_state = None
            if 'job' in v.keys():
                if 'state' in v['job'].keys():
                    job_state = v['job']['state']
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
        # check and decide what to do 
        # Took care of not resubmitting when there are RUNNING  / PENDING jobs 
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

        # decide on a final state:
        dd = {}
        for k,v in d.items():
            if v['file'] == 1:
                file_state = 1
            elif v['file'] == 0 and v['job'] in ["RUNNING","PENDING"]:
                js = v['job']
                print(f'WARNING: {fam:} {k}: job sat is {js} - Skipping')
                file_state = 1
            else:
                file_state = 0
            dd.update({k : str(file_state)})

        state = "".join([str(x) for x in dd.values()])
        
        
        if state == '00':
            mode = 'both'
        elif state == '10':
            mode = 'cluster'
        elif state == '01':
            mode = 'search'
        elif state == '11':
            mode = 'both'
        else:
            print(f'Unknown state: {state}')
            quit()

        # IMPORTANT: check if there jobs that are already running for this family!

        # Check if the number of iterations has been exceeded:
        if 'iter' in json_info[fam]['search'].keys():
            iter = json_info[fam]['search']['iter']
        elif 'iter' in json_info[fam]['cluster'].keys():
            iter = json_info[fam]['cluster']['iter']
        else: 
            iter = 0
        if iter > max_iter:
            print(f'{fam}: Can not submit the job: iter {iter} > max_iter {max_iter}')
        else:
            job_ids = submit_family(config_fn = c_fn,pref = pref,family = family,mode = mode,time_s1=None, time_s2=None, mem_s1=None, mem_s2=None,dry = dry,verbose = verbose)  
            submit_counter += 1  
            for k,v in job_ids.items():
                job_state = check_job(v)
                if not 'iter' in json_info[fam][k].keys():
                    json_info[fam][k].update({'iter' : 1})
                else:
                    json_info[fam][k].update({'iter' : int(json_info[fam][k]['iter']+1)})
                # update job info
                json_info[str(fam)][str(k)]['job'] = check_job(v)

    print(f'Submitted {submit_counter} new jobs!')

    if update_json and json_fn:
        with open(json_fn, 'w') as f:
            json.dump(json_info, f, indent=2)
        print(f'Updated: {json_fn}')
    else:
        print("Skipping job re-submission")
        quit()



