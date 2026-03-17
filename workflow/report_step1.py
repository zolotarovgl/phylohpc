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
    • Hover tooltips showing counts
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


# ── Species prefix extraction ─────────────────────────────────────────────────

def get_species_prefix(gene_id: str) -> str:
    """Return species prefix from a gene ID such as 'Mmus_ENSG123' → 'Mmus'."""
    for sep in ("_", "."):
        if sep in gene_id:
            return gene_id.split(sep)[0]
    return gene_id


# ── File parsers ──────────────────────────────────────────────────────────────

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
    counts: dict = defaultdict(int)
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
    info: dict = {}
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


# ── Tree helpers ──────────────────────────────────────────────────────────────

def get_species_order(tree_path: str) -> list:
    """Return leaf names in depth-first (left-to-right) order."""
    try:
        from ete3 import Tree  # type: ignore
        t = Tree(tree_path, format=1)
        return [n.name for n in t.get_leaves()]
    except Exception:
        return []


def tree_to_dict(node) -> dict:
    """Recursively convert ete3 node to JSON-serialisable dict."""
    d: dict = {"name": node.name or "", "dist": round(float(node.dist), 6)}
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


# ── Data assembly ─────────────────────────────────────────────────────────────

def build_family_records(search_dir: Path, family_info: dict) -> list:
    """Return list of family records from *.genes.list files.

    Each record::

        {id, family, pref, class, species_counts: {sp: n}, total}
    """
    records = []
    for genes_file in sorted(search_dir.glob("*.genes.list")):
        name = genes_file.name
        if not name.endswith(".genes.list"):
            continue
        stem = name[: -len(".genes.list")]   # "{pref}.{family}"
        parts = stem.split(".", 1)
        if len(parts) < 2:
            continue
        pref, family = parts[0], parts[1]
        genes = parse_genes_list(genes_file)
        if not genes:
            continue
        sp_counts: dict = defaultdict(int)
        for g in genes:
            sp_counts[get_species_prefix(g)] += 1
        cls = family_info.get(family, pref)
        records.append({
            "id":            f"{pref}.{family}",
            "family":        family,
            "pref":          pref,
            "class":         cls,
            "species_counts": dict(sp_counts),
            "total":         len(genes),
        })
    return records


def build_hg_records(cluster_dir: Path, family_info: dict) -> list:
    """Return list of HG records from *.fasta files in the cluster directory.

    Each record::

        {id, family, pref, class, species_counts: {sp: n}, total}
    """
    records = []
    for fasta_file in sorted(cluster_dir.glob("*.fasta")):
        stem = fasta_file.stem          # "{pref}.{family}.{hg_id}"
        parts = stem.split(".", 2)
        if len(parts) < 3:
            continue
        pref, family, hg_id = parts
        sp_counts = parse_fasta_species(fasta_file)
        if not sp_counts:
            continue
        cls = family_info.get(family, pref)
        records.append({
            "id":            stem,
            "family":        family,
            "pref":          pref,
            "hg":            hg_id,
            "class":         cls,
            "species_counts": sp_counts,
            "total":         sum(sp_counts.values()),
        })
    return records


# ── HTML template ─────────────────────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Step1 Report – Gene Family Counts</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
         font-size: 12px; background: #f7f7f7; color: #333; }
  #header { padding: 12px 16px; background: #2c3e50; color: #ecf0f1; }
  #header h2 { font-size: 16px; font-weight: 600; margin-bottom: 6px; }
  #tabs { display: flex; flex-wrap: wrap; gap: 4px; }
  .tab { padding: 3px 10px; border-radius: 12px; cursor: pointer;
         background: #4a6278; color: #dde; font-size: 11px;
         border: 1px solid transparent; user-select: none; }
  .tab:hover { background: #607d8b; }
  .tab.active { background: #1abc9c; color: #fff; border-color: #16a085; }
  #view-toggle { margin-left: auto; display: flex; align-items: center; gap: 6px; }
  .vtog { padding: 3px 10px; border-radius: 12px; cursor: pointer;
          background: #4a6278; color: #dde; font-size: 11px; }
  .vtog.active { background: #e67e22; color: #fff; }
  #layout { display: flex; height: calc(100vh - 56px); overflow: hidden; }
  #tree-panel { flex: 0 0 220px; overflow: hidden; border-right: 1px solid #ccc;
                background: #fff; }
  #heatmap-panel { flex: 1 1 auto; overflow: auto; background: #fff;
                   position: relative; }
  #heatmap-panel svg { display: block; }
  .node circle { fill: #2c3e50; r: 2; }
  .node text { font-size: 10px; }
  .link { fill: none; stroke: #aaa; stroke-width: 1px; }
  #tooltip { position: fixed; padding: 6px 10px; background: rgba(0,0,0,.78);
             color: #fff; border-radius: 4px; font-size: 11px; pointer-events: none;
             display: none; z-index: 100; max-width: 260px; line-height: 1.5; }
  .col-label { font-size: 9px; }
  .cell:hover { stroke: #e74c3c; stroke-width: 1px; cursor: default; }
  #summary { padding: 4px 16px; background: #ecf0f1; border-bottom: 1px solid #ccc;
             font-size: 11px; color: #666; }
</style>
</head>
<body>
<div id="header">
  <h2>Step 1 – Gene Family Report</h2>
  <div id="tabs"></div>
</div>
<div id="summary" id="summary-bar">Loading…</div>
<div id="layout">
  <div id="tree-panel"></div>
  <div id="heatmap-panel"></div>
</div>
<div id="tooltip"></div>

<script>
// ── Embedded data ─────────────────────────────────────────────────────────────
const SPECIES_ORDER = %%SPECIES_ORDER%%;
const TREE_DATA     = %%TREE_DATA%%;
const FAMILY_DATA   = %%FAMILY_DATA%%;
const HG_DATA       = %%HG_DATA%%;

// ── Constants ─────────────────────────────────────────────────────────────────
const ROW_H   = 14;          // px per species row
const COL_W   = 18;          // px per family/HG column
const TREE_W  = 215;         // px for the tree panel
const MARGIN  = { top: 80, right: 10, bottom: 10, left: 4 };
const MAX_LABEL_W = 100;     // species label area in tree panel

// ── State ─────────────────────────────────────────────────────────────────────
let activeClass  = "all";
let activeView   = "family";   // "family" | "hg"
let currentData  = [];

// ── Helpers ───────────────────────────────────────────────────────────────────
const tooltip = d3.select("#tooltip");
function showTip(html, event) {
  tooltip.style("display", "block").html(html)
    .style("left", (event.clientX + 14) + "px")
    .style("top",  (event.clientY - 10) + "px");
}
function hideTip() { tooltip.style("display", "none"); }

// ── Class tabs ────────────────────────────────────────────────────────────────
const allClasses = ["all", ...new Set(FAMILY_DATA.map(d => d.class).filter(Boolean).sort())];
const tabsEl = d3.select("#tabs");
tabsEl.selectAll(".tab").data(allClasses).join("div")
  .attr("class", d => "tab" + (d === "all" ? " active" : ""))
  .text(d => d)
  .on("click", function(_, cls) {
    tabsEl.selectAll(".tab").classed("active", d => d === cls);
    activeClass = cls;
    updateHeatmap();
  });

// ── Species ordering ──────────────────────────────────────────────────────────
// Use SPECIES_ORDER when available; fall back to union of all species in data.
let speciesOrder = SPECIES_ORDER.length
  ? SPECIES_ORDER
  : [...new Set([...FAMILY_DATA, ...HG_DATA].flatMap(d => Object.keys(d.species_counts)))];

// Keep only species that appear in the data
const knownSpecies = new Set([...FAMILY_DATA, ...HG_DATA]
  .flatMap(d => Object.keys(d.species_counts)));
speciesOrder = speciesOrder.filter(s => knownSpecies.has(s));
// Append any species in data but not in tree
for (const s of knownSpecies) { if (!speciesOrder.includes(s)) speciesOrder.push(s); }

const N = speciesOrder.length;

// ── Tree cladogram ─────────────────────────────────────────────────────────────
function drawTree() {
  const el = document.getElementById("tree-panel");
  const H  = MARGIN.top + N * ROW_H + MARGIN.bottom;
  const W  = TREE_W;

  const svg = d3.select("#tree-panel").append("svg")
    .attr("width", W).attr("height", H);

  if (!TREE_DATA || !TREE_DATA.name) {
    // No tree: just draw species labels
    svg.selectAll("text").data(speciesOrder).join("text")
      .attr("x", W - 4).attr("y", (_, i) => MARGIN.top + i * ROW_H + ROW_H / 2 + 4)
      .attr("text-anchor", "end").attr("font-size", 10)
      .text(d => d);
    return;
  }

  // Build leaf-to-row-index lookup
  const leafY = {};
  speciesOrder.forEach((s, i) => { leafY[s] = MARGIN.top + i * ROW_H + ROW_H / 2; });

  // Assign Y to every node (internal nodes get average of children leaf Y)
  function assignY(node) {
    if (node.leaf || !node.children) {
      node._y = leafY[node.name] ?? 0;
      return node._y;
    }
    const ys = node.children.map(assignY);
    node._y = (Math.min(...ys) + Math.max(...ys)) / 2;
    return node._y;
  }
  assignY(TREE_DATA);

  // Assign X by depth (proportion of total max depth)
  function maxDepth(n) {
    if (!n.children) return 0;
    return 1 + Math.max(...n.children.map(maxDepth));
  }
  const MD = maxDepth(TREE_DATA);

  function assignX(node, depth) {
    node._x = (depth / MD) * (W - MAX_LABEL_W - 10);
    if (node.children) node.children.forEach(c => assignX(c, depth + 1));
  }
  assignX(TREE_DATA, 0);

  // Draw links
  function drawLinks(node) {
    if (!node.children) return;
    node.children.forEach(child => {
      // Elbow connector: horizontal then vertical then horizontal
      svg.append("path")
        .attr("class", "link")
        .attr("d", `M${node._x},${node._y} H${child._x} V${child._y} H${child._x}`);
      drawLinks(child);
    });
  }
  drawLinks(TREE_DATA);

  // Draw leaf labels
  function drawLabels(node) {
    if (node.leaf || !node.children) {
      svg.append("text")
        .attr("x", W - 4)
        .attr("y", node._y + 4)
        .attr("text-anchor", "end")
        .attr("font-size", Math.min(10, ROW_H - 2))
        .text(node.name);
    } else {
      node.children.forEach(drawLabels);
    }
  }
  drawLabels(TREE_DATA);
}

drawTree();

// ── Heatmap ───────────────────────────────────────────────────────────────────
let hmSvg = null;

function filteredData() {
  const src = activeView === "family" ? FAMILY_DATA : HG_DATA;
  return activeClass === "all" ? src : src.filter(d => d.class === activeClass);
}

function updateHeatmap() {
  currentData = filteredData();
  const M = currentData.length;

  // Update summary bar
  const totalGenes = currentData.reduce((a, d) => a + (d.total || 0), 0);
  d3.select("#summary").text(
    `Showing ${M} ${activeView === "family" ? "families" : "HGs"}` +
    (activeClass !== "all" ? ` [class: ${activeClass}]` : " [all classes]") +
    ` · ${totalGenes.toLocaleString()} genes total`
  );

  const W = MARGIN.left + M * COL_W + 10;
  const H = MARGIN.top  + N * ROW_H + MARGIN.bottom;

  const panel = d3.select("#heatmap-panel");
  panel.selectAll("*").remove();

  if (M === 0) {
    panel.append("p").style("padding", "20px").style("color", "#999")
      .text("No data for this selection.");
    return;
  }

  hmSvg = panel.append("svg").attr("width", W).attr("height", H);

  // Color scale
  const maxCount = d3.max(currentData, d => d3.max(Object.values(d.species_counts))) || 1;
  const color = d3.scaleSequential([0, maxCount], d3.interpolateBlues);

  // X scale (column index → pixel)
  const xScale = d3.scaleBand()
    .domain(d3.range(M))
    .range([MARGIN.left, MARGIN.left + M * COL_W])
    .padding(0.05);

  // Y scale (species order)
  const yScale = d3.scaleBand()
    .domain(speciesOrder)
    .range([MARGIN.top, MARGIN.top + N * ROW_H])
    .padding(0.05);

  // Column labels (rotated)
  hmSvg.selectAll(".col-label").data(currentData).join("text")
    .attr("class", "col-label")
    .attr("transform", (_, i) =>
      `translate(${xScale(i) + xScale.bandwidth() / 2},${MARGIN.top - 4}) rotate(-65)`)
    .attr("text-anchor", "start")
    .attr("font-size", Math.min(9, COL_W - 1))
    .text(d => activeView === "family" ? d.family : (d.hg || d.id))
    .append("title").text(d => `${d.id} [${d.class}] total=${d.total}`);

  // Cells
  currentData.forEach((rec, ci) => {
    speciesOrder.forEach(sp => {
      const n = rec.species_counts[sp] || 0;
      hmSvg.append("rect")
        .attr("class", "cell")
        .attr("x", xScale(ci))
        .attr("y", yScale(sp))
        .attr("width",  xScale.bandwidth())
        .attr("height", yScale.bandwidth())
        .attr("fill", n === 0 ? "#f5f5f5" : color(n))
        .on("mouseover", (event) => {
          showTip(
            `<b>${sp}</b><br>${rec.id}<br>` +
            `class: ${rec.class}<br>genes: <b>${n}</b>`,
            event
          );
        })
        .on("mouseout", hideTip);
    });
  });

  // Colour legend
  const legW = 120, legH = 10;
  const legX = MARGIN.left, legY = H - 20;
  const legDefs = hmSvg.append("defs");
  const grad = legDefs.append("linearGradient").attr("id", "legGrad");
  grad.selectAll("stop").data(d3.range(0, 1.01, 0.1)).join("stop")
    .attr("offset", d => d)
    .attr("stop-color", d => color(d * maxCount));
  hmSvg.append("rect").attr("x", legX).attr("y", legY)
    .attr("width", legW).attr("height", legH)
    .attr("fill", "url(#legGrad)").attr("stroke", "#aaa").attr("stroke-width", 0.5);
  hmSvg.append("text").attr("x", legX).attr("y", legY - 3)
    .attr("font-size", 9).text("0");
  hmSvg.append("text").attr("x", legX + legW).attr("y", legY - 3)
    .attr("text-anchor", "end").attr("font-size", 9).text(maxCount + " genes");
}

updateHeatmap();
</script>
</body>
</html>
"""


# ── Argument parsing and main ─────────────────────────────────────────────────

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate interactive HTML report for step1 outputs.")
    p.add_argument("--search_dir",  default="results/search",
                   help="Directory containing *.genes.list files")
    p.add_argument("--cluster_dir", default="results/clusters",
                   help="Directory containing *.fasta cluster files")
    p.add_argument("--family_info", default="data/gene_families_searchinfo.csv",
                   help="gene_families_searchinfo.csv (TSV, 7 columns)")
    p.add_argument("--tree",        default="data/species_tree.full.newick",
                   help="Newick species tree for species ordering")
    p.add_argument("--output",      required=True,
                   help="Output HTML file path")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    search_dir  = Path(args.search_dir)
    cluster_dir = Path(args.cluster_dir)

    # Load family → class mapping
    family_info = load_family_info(args.family_info)

    # Load species order from tree
    species_order, tree_dict = load_tree_data(args.tree)

    # Build data records
    family_records = build_family_records(search_dir, family_info)
    hg_records     = build_hg_records(cluster_dir, family_info)

    # Serialise to JSON
    species_json  = json.dumps(species_order)
    tree_json     = json.dumps(tree_dict)
    family_json   = json.dumps(family_records)
    hg_json       = json.dumps(hg_records)

    html = (HTML_TEMPLATE
            .replace("%%SPECIES_ORDER%%", species_json)
            .replace("%%TREE_DATA%%",     tree_json)
            .replace("%%FAMILY_DATA%%",   family_json)
            .replace("%%HG_DATA%%",       hg_json))

    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
