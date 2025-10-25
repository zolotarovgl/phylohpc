import subprocess
import os
import json
import argparse 


def parse_genefam(fn):
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
                genefam.update({
                    family: 
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


def check_family(file_glob,config,status = None,variable = None,verbose = False):
    import glob
    import os
    SEARCH_DIR = config['SEARCH_DIR']
    files = glob.glob(file_glob)
    toreplace = file_glob.split('*')[1]
    fams = [os.path.basename(f).replace(toreplace,'').split('.')[1] for f in files]
    if verbose:
        print(f'Found {len(fams)} family search results.')
    fams = [x for x in fams if x in genefam.keys()]
    if verbose:
        print(f'{len(fams)} from genefam.')

    if not status:
        status = {}
    for fam in genefam.keys():
        if fam in fams:
            state = 1
        else:
            state = 0
        if not fam in status.keys():
            status.update({fam: {}})
        status[fam].update({variable: state})
    return(status)

def count_lines(path):
	with open(path) as f:
		return sum(1 for _ in f)


########################################
# Search job results   
########################################
parser = argparse.ArgumentParser(description="A python wrapper to check families")
parser.add_argument("configfile", help="Bash-style config file with KEY=VALUE lines")
parser.add_argument("--genefam", required=True, help="Path to gene family definition file")
parser.add_argument("--json", help="Optional JSON file to store the info")
parser.add_argument("--output",default = None, help="Optional output tab file with family states")
args = parser.parse_args()



g_fn = 'genefam.csv'
c_fn = 'configs/config.txt'
json_fn = 'fam.info.json'
out_fn = 'boo.tab'

g_fn = args.genefam
c_fn = args.configfile
json_fn = args.json
out_fn = args.output

update = True


genefam = parse_genefam(g_fn)
config = parse_bash_config(c_fn)

# check family status 
file_glob = f"{config['SEARCH_DIR']}/*.domains.fasta"
status = check_family(file_glob,config,status = None,variable = 'search',verbose = False)
file_glob = f"{config['CLUSTER_DIR']}/*_cluster.tsv"
status = check_family(file_glob,config,status = status,variable = 'cluster',verbose = False)

######################
# Number of sequences   
######################
nseq = {}
for fam in genefam.keys():
    pat = genefam[fam]['prefix'] + '.' + fam
    fn = f"{config['SEARCH_DIR']}/{pat}.domains.csv"
    if os.path.isfile(fn):
        n = count_lines(fn)
    else:
        n = 0
    nseq.update({fam: n})


######################
# Update the JSON file  
######################
# update the file 
if update and json_fn:
    with open(json_fn, 'w') as f:
        json.dump(status, f, indent=2)

######################
# Report 
######################
states = [str("".join([str(v) for v in status[fam].values()])) for fam in status.keys()]
from collections import Counter
counts = Counter(states)
print(f'\nFamily status [search,clustering]:')
for state, n in counts.items():    
    print(f"{state}: {n}")
print(f"Total: {sum(counts.values())}")


states = {fam: str("".join([str(v) for v in status[fam].values()])) for fam in status.keys()}


if out_fn:
    with open(out_fn, "w") as f:
        for k, v in states.items():
            f.write(f"{k}\t{v}\n")
    print(f'Family states written to {out_fn}')


