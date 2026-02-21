import sys

inp=sys.argv[1]
out=sys.argv[2]

with open(inp) as f, open(out,"w") as o:
	name=None
	seq=[]
	for line in f:
		line=line.rstrip()
		if line.startswith(">"):
			if name:
				s="".join(seq)
				if s.replace("-","").replace(".","")!="":
					o.write(name+"\n"+s+"\n")
			name=line
			seq=[]
		else:
			seq.append(line)
	if name:
		s="".join(seq)
		if s.replace("-","").replace(".","")!="":
			o.write(name+"\n"+s+"\n")
