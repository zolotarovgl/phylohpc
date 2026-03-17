#!/usr/bin/env python3
"""Generate an interactive HTML report for step1 (gene family search + clustering) outputs.

Reads
-----
  results/search/*.genes.list   one gene per line
  results/clusters/*.fasta      per-HG FASTA files (headers give gene IDs)
  data/gene_families_searchinfo.csv   family → class mapping (7-column TSV)
  data/species_tree.full.newick (or any newick) for species ordering

Outputs
-------
  A self-contained HTML file with:
    • A phylogenetic tree cladogram (left panel) aligned to species rows
    • A scrollable heatmap of gene counts (species × families or HGs)
    • Filter tabs by gene family class
    • Hover tooltips showing counts and z-scores
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


def get_species_prefix(gene_id: str) -> str:
    """Return species prefix from a gene ID such as 'Mmus_ENSG123' → 'Mmus'."""
    for sep in ("_", "."):
        if sep in gene_id:
            return gene_id.split(sep)[0]
    return gene_id


def parse_genes_list(path) -> list:
    """Return list of gene IDs from a *.genes.list file."""
    genes = []
    try:
        with open(path) as fh:
            for line in fh:
                g = line.strip()
                if g and not g.startswith("#"):
                    genes.append(g)
    except (FileNotFoundError, OSError):
        pass
    return genes


def parse_fasta_species(path) -> dict:
    """Return {species: count} by reading gene-ID headers of a FASTA file."""
    counts = defaultdict(int)
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    gene = line[1:].strip().split()[0]
                    counts[get_species_prefix(gene)] += 1
    except (FileNotFoundError, OSError):
        pass
    return dict(counts)


def load_family_info(csv_path) -> dict:
    """Return {family_name: class_label} from gene_families_searchinfo.csv."""
    info = {}
    try:
        with open(csv_path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 7:
                    family = cols[0].strip()
                    cls = cols[6].strip()
                    if family and cls:
                        info[family] = cls
    except (FileNotFoundError, OSError):
        pass
    return info


def tree_to_dict(node) -> dict:
    """Recursively convert ete3 node to JSON-serialisable dict."""
    d = {"name": node.name or "", "dist": round(float(node.dist), 6)}
    if node.is_leaf():
        d["leaf"] = True
    else:
        d["children"] = [tree_to_dict(c) for c in node.children]
    return d


def load_tree_data(tree_path: str):
    """Return (species_order, tree_dict) using ete3, or ([], {}) on failure."""
    try:
        from ete3 import Tree  # type: ignore
        t = Tree(tree_path, format=1)
        species_order = [n.name for n in t.get_leaves()]
        tree_dict = tree_to_dict(t)
        return species_order, tree_dict
    except Exception:
        return [], {}


def build_family_records(search_dir: Path, family_info: dict) -> list:
    """Return list of family records from *.genes.list files."""
    records = []
    for genes_file in sorted(search_dir.glob("*.genes.list")):
        name = genes_file.name
        if not name.endswith(".genes.list"):
            continue
        stem = name[: -len(".genes.list")]
        parts = stem.split(".", 1)
        if len(parts) < 2:
            continue
        pref, family = parts[0], parts[1]
        genes = parse_genes_list(genes_file)
        if not genes:
            continue
        sp_counts = defaultdict(int)
        for g in genes:
            sp_counts[get_species_prefix(g)] += 1
        cls = family_info.get(family, pref)
        records.append({
            "id": f"{pref}.{family}",
            "family": family,
            "pref": pref,
            "class": cls,
            "species_counts": dict(sp_counts),
            "total": len(genes),
        })
    return records


def build_hg_records(cluster_dir: Path, family_info: dict) -> list:
    """Return list of HG records from *.fasta files in the cluster directory."""
    records = []
    for fasta_file in sorted(cluster_dir.glob("*.fasta")):
        stem = fasta_file.stem
        parts = stem.split(".", 2)
        if len(parts) < 3:
            continue
        pref, family, hg_id = parts
        sp_counts = parse_fasta_species(fasta_file)
        if not sp_counts:
            continue
        cls = family_info.get(family, pref)
        records.append({
            "id": stem,
            "family": family,
            "pref": pref,
            "hg": hg_id,
            "class": cls,
            "species_counts": sp_counts,
            "total": sum(sp_counts.values()),
        })
    return records


#!/usr/bin/env python3

# (Python part unchanged — omitted for brevity, keep yours exactly the same)


#!/usr/bin/env python3

# (Python part unchanged — keep your existing code)

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Step1 Report</title>
<script src="https://d3js.org/d3.v7.min.js"></script>

<style>
body { font-family: sans-serif; }
#header { background:#2c3e50;color:white;padding:10px; display:flex; gap:10px; align-items:center; }
#layout { display:flex; height:90vh; }
#tree-panel { width:280px; border-right:1px solid #ccc; }
#heatmap-panel { flex:1; overflow:auto; }
.cell:hover { stroke:black; stroke-width:1px; }
.col { font-size:9px; }

#tooltip {
  position:fixed;
  display:none;
  background:rgba(0,0,0,0.85);
  color:white;
  padding:6px 10px;
  border-radius:4px;
  font-size:11px;
  pointer-events:none;
}
</style>
</head>

<body>

<div id="header">
  Step1 Report
  <button id="backBtn" style="display:none;">⬅ Back</button>

  Prefix:
  <select id="prefixSelect"></select>
</div>

<div id="layout">
  <div id="tree-panel"></div>
  <div id="heatmap-panel"></div>
</div>

<div id="tooltip"></div>

<script>

// ── CONSTANTS ────────────────────────────────
const TOP_MARGIN = 110;

// ── DATA ────────────────────────────────────
const SPECIES_ORDER = %%SPECIES_ORDER%%;
const TREE_DATA     = %%TREE_DATA%%;
const FAMILY_DATA   = %%FAMILY_DATA%%;
const HG_DATA       = %%HG_DATA%%;

// ── STATE ───────────────────────────────────
let viewMode = "family";
let activeFamily = null;
let activePrefix = "all";

// ── TOOLTIP ─────────────────────────────────
const tooltip = d3.select("#tooltip");

function showTip(html, event) {
  tooltip.style("display","block")
    .html(html)
    .style("left",(event.clientX+12)+"px")
    .style("top",(event.clientY-10)+"px");
}
function hideTip(){ tooltip.style("display","none"); }

// ── PREFIX SELECTOR ─────────────────────────
const allPrefixes = ["all", ...new Set(FAMILY_DATA.map(d => d.pref).filter(Boolean).sort())];

const prefixSelect = d3.select("#prefixSelect");

prefixSelect.selectAll("option")
  .data(allPrefixes)
  .enter()
  .append("option")
  .attr("value",d=>d)
  .text(d=>d);

prefixSelect.on("change", function(){
  activePrefix=this.value;
  drawHeatmap();
});

// ── SPECIES ORDER ───────────────────────────
const dataSpecies = new Set([...FAMILY_DATA, ...HG_DATA]
  .flatMap(d => Object.keys(d.species_counts)));

let speciesOrder = SPECIES_ORDER.filter(s => dataSpecies.has(s));
for (const s of dataSpecies) {
  if (!speciesOrder.includes(s)) speciesOrder.push(s);
}

// ── TREE ────────────────────────────────────
function drawTree(){

  const svg = d3.select("#tree-panel").html("")
    .append("svg")
    .attr("width",280)
    .attr("height", speciesOrder.length*14 + TOP_MARGIN + 40);

  const leafY = {};
  speciesOrder.forEach((s,i)=> leafY[s]=TOP_MARGIN+i*14);

  function clone(n){ return JSON.parse(JSON.stringify(n)); }

  function prune(n){
    if(!n.children) return speciesOrder.includes(n.name)?n:null;
    const kids=n.children.map(prune).filter(x=>x);
    if(kids.length===0) return null;
    if(kids.length===1) return kids[0];
    n.children=kids;
    return n;
  }

  let tree = prune(clone(TREE_DATA));
  if(!tree) return;

  function assignY(n){
    if(!n.children){ n._y=leafY[n.name]; return n._y; }
    const ys=n.children.map(assignY);
    n._y=d3.mean(ys); return n._y;
  }

  function assignX(n,d=0){
    n._d=d;
    if(n.children) n.children.forEach(c=>assignX(c,d+1));
  }

  function flat(n){
    return [n].concat(n.children? n.children.flatMap(flat):[]);
  }

  assignY(tree); assignX(tree);

  const maxD=d3.max(flat(tree),d=>d._d);
  const sx=d=>10+(d/maxD)*150;

  function draw(n){
    if(!n.children) return;

    const ys=n.children.map(c=>c._y);

    svg.append("line")
      .attr("x1",sx(n._d)).attr("x2",sx(n._d))
      .attr("y1",d3.min(ys)).attr("y2",d3.max(ys))
      .attr("stroke","#aaa");

    n.children.forEach(c=>{
      svg.append("line")
        .attr("x1",sx(n._d)).attr("x2",sx(c._d))
        .attr("y1",c._y).attr("y2",c._y)
        .attr("stroke","#aaa");
      draw(c);
    });
  }

  draw(tree);

  speciesOrder.forEach((s,i)=>{
    svg.append("text")
      .attr("x",270)
      .attr("y",leafY[s]+4)
      .attr("text-anchor","end")
      .text(s);
  });
}

// ── HEATMAP ─────────────────────────────────
function drawHeatmap(){

  let data;

  if(viewMode==="family"){
    data = FAMILY_DATA
      .filter(d => activePrefix==="all" || d.pref===activePrefix)
      .slice()
      .sort((a,b)=>b.total-a.total);
  } else {
    data = HG_DATA
      .filter(d => d.family===activeFamily)
      .filter(d => activePrefix==="all" || d.pref===activePrefix)
      .sort((a,b)=>b.total-a.total);
  }

  const svg = d3.select("#heatmap-panel").html("")
    .append("svg")
    .attr("width", data.length*18+100)
    .attr("height", speciesOrder.length*14 + TOP_MARGIN + 40);

  // z-score
  const zMat = data.map(rec=>{
    const vals=speciesOrder.map(s=>rec.species_counts[s]||0);
    const m=d3.mean(vals);
    const sd=d3.deviation(vals)||1;
    return vals.map(v=>(v-m)/sd);
  });

  const zAll=zMat.flat();
  const zMax=Math.max(Math.abs(d3.min(zAll)), Math.abs(d3.max(zAll)));

  const color=d3.scaleDiverging()
    .domain([zMax,0,-zMax])
    .interpolator(d3.interpolateRdBu);

  data.forEach((rec,ci)=>{
    speciesOrder.forEach((sp,ri)=>{

      const count=rec.species_counts[sp]||0;
      const z=zMat[ci][ri];

      svg.append("rect")
        .attr("class","cell")
        .attr("x",ci*18+60)
        .attr("y",ri*14 + TOP_MARGIN)
        .attr("width",16)
        .attr("height",12)
        .attr("fill",color(z))

        .on("mouseover",(event)=>{
          showTip(
            `<b>${sp}</b><br>` +
            `${viewMode==="family" ? rec.family : rec.id}<br>` +
            `class: ${rec.class}<br>` +
            `count: <b>${count}</b><br>` +
            `z-score: <b>${z.toFixed(2)}</b>`,
            event
          );
        })
        .on("mousemove",(event)=>{
          tooltip
            .style("left",(event.clientX+12)+"px")
            .style("top",(event.clientY-10)+"px");
        })
        .on("mouseout",hideTip);
    });
  });

  // column labels
  svg.selectAll("text.col")
    .data(data)
    .enter()
    .append("text")
    .attr("class","col")
    .attr("transform",(d,i)=>
      `translate(${i*18+60},${TOP_MARGIN-10}) rotate(-65)`
    )
    .text(d=> viewMode==="family" ? d.family : d.id)
    .style("cursor","pointer")
    .on("click",(event,d)=>{
      if(viewMode==="family"){
        viewMode="hg";
        activeFamily=d.family;
        document.getElementById("backBtn").style.display="inline";
        drawHeatmap();
      }
    });
}

// ── BACK BUTTON ─────────────────────────────
document.getElementById("backBtn").onclick=()=>{
  viewMode="family";
  activeFamily=null;
  document.getElementById("backBtn").style.display="none";
  drawHeatmap();
};

// ── INIT ────────────────────────────────────
drawTree();
drawHeatmap();

</script>
</body>
</html>
"""


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate interactive HTML report for step1 outputs."
    )
    p.add_argument("--search_dir", default="results/search",
                   help="Directory containing *.genes.list files")
    p.add_argument("--cluster_dir", default="results/clusters",
                   help="Directory containing *.fasta cluster files")
    p.add_argument("--family_info", default="data/gene_families_searchinfo.csv",
                   help="gene_families_searchinfo.csv (TSV, 7 columns)")
    p.add_argument("--tree", default="data/species_tree.full.newick",
                   help="Newick species tree for species ordering")
    p.add_argument("--output", required=True,
                   help="Output HTML file path")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    search_dir = Path(args.search_dir)
    cluster_dir = Path(args.cluster_dir)

    family_info = load_family_info(args.family_info)
    species_order, tree_dict = load_tree_data(args.tree)

    family_records = build_family_records(search_dir, family_info)
    hg_records = build_hg_records(cluster_dir, family_info)

    species_json = json.dumps(species_order)
    tree_json = json.dumps(tree_dict)
    family_json = json.dumps(family_records)
    hg_json = json.dumps(hg_records)

    html = (HTML_TEMPLATE
            .replace("%%SPECIES_ORDER%%", species_json)
            .replace("%%TREE_DATA%%", tree_json)
            .replace("%%FAMILY_DATA%%", family_json)
            .replace("%%HG_DATA%%", hg_json))

    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()