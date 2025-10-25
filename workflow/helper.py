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