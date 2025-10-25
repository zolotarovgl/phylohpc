#!/usr/bin/env python3
import argparse, subprocess, os, tempfile

parser = argparse.ArgumentParser(description="Gather species annotation")
parser.add_argument("--config",required = True, help = 'Bash-formatted config file configs/config.txt')
parser.add_argument("--id", required=True, help = 'Species prefox')
parser.add_argument("--outfile", help = 'Optional output tab file. If not specified, the annotations are printed to stdout')
parser.add_argument("--search-dir", default=None)
parser.add_argument("--tree-dir", default=None)
parser.add_argument("--prefix", default=None, help = 'Optional prefix of the families, e.g. tfs')
args = parser.parse_args()

ID = args.id
SEARCH_DIR = args.search_dir
TREE_DIR = args.tree_dir
OUTFILE = args.outfile
verbose = False


from workflow.helper import parse_bash_config
config_fn = 'configs/config.txt'
config = parse_bash_config(config_fn)
if not SEARCH_DIR:
	SEARCH_DIR = config['SEARCH_DIR']
if not TREE_DIR:
	TREE_DIR = config['TREE_DIR']

def run(cmd):
	r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	return r.stdout.strip()

prefix = args.prefix
if not prefix:
	prefix = ''



#with tempfile.TemporaryDirectory() as tmp:
tmp = tempfile.mkdtemp(prefix="gather_anno_")
if OUTFILE and verbose:
    print(f"Temporary directory: {tmp}")
tmp_anno = os.path.join(tmp, "tmp_anno")
ids_todo = os.path.join(tmp, "ids_todo")
gene2class = os.path.join(tmp, "gene2class")
pep2hg = os.path.join(tmp, "pep2hg")

# 1. tmp_anno
cmd = f"cat {TREE_DIR}/{prefix}*groups.csv | grep '{ID}' > {tmp_anno}"
if OUTFILE and verbose:
    print(cmd)
subprocess.run(cmd, shell=True)

# 2. ids_todo
cmd = f"cat {SEARCH_DIR}/{prefix}*.genes.list | grep '{ID}' | sort | uniq > {ids_todo}"
if OUTFILE and verbose:
    print(cmd)
subprocess.run(cmd, shell=True)

# 3. gene2class
with open(gene2class, "w") as out:
    files = files = [f for f in os.listdir(SEARCH_DIR) if f.startswith(prefix)]
    for f in files:
        if f.endswith(".genes.list"):
            pref = f.replace(".genes.list", "")
            cmd = f"grep '{ID}' {os.path.join(SEARCH_DIR,f)} | awk -v PREF={pref!r} '{{print $1\"\\t\"PREF}}'"
        #if OUTFILE:
            #print(cmd)
        subprocess.run(cmd, shell=True, stdout=out)

# 4. pep2hg
cmd = "bash workflow/get_gene2cluster.sh configs/config.txt pep2hg"
if OUTFILE and verbose:
    print(cmd)
subprocess.run(cmd, shell=True)
cmd = f"grep -E '^{ID}' pep2hg > {pep2hg}"
if OUTFILE and verbose:
    print(cmd)
subprocess.run(cmd, shell=True)

# 5. RESULT (awk logic)
awk_cmd = (
	f"awk 'NR==FNR{{dict[$1]=$2\"\\t\"$3\"\\t\"$4;next}}{{if($1 in dict) print $1\"\\t\"dict[$1]; "
	f"else print $1\"\\tUnclassified\"}}' {tmp_anno} {ids_todo} | "
	f"awk 'BEGIN{{OFS=\"\\t\"}}FNR==NR{{d[$1]=$2;next}}{{if($2==\":Unclassified\"){{$2=d[$1]\":Unclassified\"}};print $0}}' {pep2hg} - | "
	f"awk 'FNR==NR{{d[$1]=$2;next}}{{if($2 ~ /Unclass/){{$2=d[$1]\":Unclassified\"}};print $0}}' {pep2hg} - | "
	f"awk 'FNR==NR{{d[$1]=$2;next}}{{if($2==\":Unclassified\"){{$2=d[$1]\":Unclassified\"}};print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4}}' {gene2class} -"
)

if OUTFILE and verbose:
    print(awk_cmd)
RESULT = run(awk_cmd)

if OUTFILE:
    with open(OUTFILE, "w") as out:
        out.write(RESULT + "\n")
else:
    print(RESULT)
