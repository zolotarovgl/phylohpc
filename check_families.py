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

# Example
#family = 'tfs.Forkhead'
# search output 
#out = check_outfile(family = family,config = config,file_glob = f"{config['SEARCH_DIR']}/*.domains.fasta")
#print(out)
#out = check_outfile(family = family,config = config,file_glob = f"{config['CLUSTER_DIR']}/*_cluster.tsv")
#print(out)
#output = check_families(genefam_dict = genefam_dict,config = config,json_fn = json_fn,verbose = False)
#quit()

def check_families(family = None,json_info = None,json_fn = None,genefam_dict = None,config = None,update_json = True,verbose = True):
    if family:
        fams = args.family.split(',')
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
        else:
            for k,v in status.items():
                json_info[family][k]['status'] = v
            
                # check the job status 
                if 'job' in json_info[family][k].keys():
                    jobid = json_info[family][k]['job']['jobid']
                    json_info[family][k]['job'] = check_job(jobid)        
        # update job info if present 
        #if 'job' in json_info['family']
        
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
args = parser.parse_args()


# Given a list of families or the whole genefam file, check their status

g_fn = args.genefam
c_fn = args.configfile
json_fn = args.json
out_fn = args.output
family = args.family 
stdout = args.stdout

update_json = True
verbose = False


json_info = read_json(json_fn)
genefam_dict = parse_genefam(g_fn)
config = parse_bash_config(c_fn)

# checks the status of the families 
json_info = check_families(family = args.family,json_info = json_info,json_fn = json_fn,genefam_dict = genefam_dict,config = config,update_json = update_json,verbose = verbose)

# If stdout, report their statuses
states = {}
for k,v in json_info.items():
    states.update({k : "".join([str(v['search']['status']),str(v['cluster']['status'])])})

if out_fn:
    with open(out_fn, "w") as f:
        for k, v in states.items():
            f.write(f"{k}\t{v}\n")
    if not stdout:
        print(f'Output: {out_fn}')

# Report job states
from collections import Counter
counts = Counter([x for x in [v for v in states.values()]])
if not stdout:
    print(f'Family status [search,clustering]:\n')
    for state, n in counts.items():    
        print(f"-  {state}: {n}")
    print(f"Total: {sum(counts.values())}")

if stdout:
    for k,v in states.items():
        print(f'{k}\t{v}')


########################################
# Resubmit jobs? 
# - use submit family functionality 
########################################
# the submit family will submit 2 jobs at the same time 
from submit_family import submit_family
resubmit = True

print(
"""
########################################
# Job re-submission
########################################
"""
)
if resubmit:  
    print('resubmit=True: checking failed families and re-submitting')

fams = [k for k,v in states.items() if v in ['00','01']]

if len(fams)>0:
    print(f'{len(fams)} families to re-submit')
else:
    print(f'No families to re-submit!')


quit()
#################################
for fam in fams:
    job_ids = submit_family(config_fn = c_fn,pref = 'vil',family = 'Fascin',dry = False,verbose = False)    
    for k,v in job_ids.items():
        job_state = check_job(v)
        json_info[str(fam)][str(k)]['job'] = check_job(v)

if update_json and json_fn:
    with open(json_fn, 'w') as f:
        json.dump(json_info, f, indent=2)

