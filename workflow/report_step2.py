#!/usr/bin/env python3
"""Generate a self-contained interactive HTML report for step2 (POSSVM) outputs.

Reads
-----
  results/possvm/*.ortholog_groups.newick   gene trees with OG labels on internal nodes
  results/possvm/*.ortholog_groups.csv      gene → OG mapping (tab-sep, no header)
    Columns: gene_id  og_name  reference_og  reference_gene_name

Memory design
-------------
The HTML embeds two kinds of data:

1. TREE_INDEX  (JS constant, parsed at startup)
   Lightweight metadata only — ~200 bytes per HG.
   No recursive tree dicts.

2. Per-HG <script type="application/json" id="treedata-{id}"> tags
   Inert strings never evaluated by JS until the user selects that HG.
   Parsed on demand via JSON.parse().

This keeps startup memory well under 1 MB even for 5 000+ HGs.

Outputs
-------
  A self-contained HTML file with:
    • Sidebar: family-class accordion → HG items (search collapses / expands groups)
    • Main panel: collapsible D3.js gene tree
        – Leaf nodes coloured by species prefix
        – OG nodes labelled and highlighted
        – Click to collapse / expand subtrees
        – Zoom / pan
        – Tooltip with gene / OG details
"""

import argparse
import html as _html
import json
import sys
from collections import defaultdict
from pathlib import Path


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_species_prefix(gene_id: str) -> str:
    """Return species prefix from a gene ID such as 'Mmus_ENSG123' → 'Mmus'."""
    for sep in ("_", "."):
        if sep in gene_id:
            return gene_id.split(sep)[0]
    return gene_id


# ── Parsers ───────────────────────────────────────────────────────────────────

def gene_tree_to_dict(node) -> dict:
    """Recursively convert an ete3 node to a JSON-serialisable dict.

    Leaf nodes include a ``species`` field.
    Internal nodes carry their name (OG label from POSSVM) if present.
    """
    name = node.name or ""
    d: dict = {"name": name, "dist": round(float(node.dist), 6)}
    if node.is_leaf():
        d["leaf"] = True
        d["species"] = get_species_prefix(name) if name else ""
    else:
        d["children"] = [gene_tree_to_dict(c) for c in node.children]
    return d


def load_og_csv(csv_path: Path) -> dict:
    """Return {og_name: [gene_id, ...]} from a POSSVM ortholog_groups.csv file."""
    og_members: dict = defaultdict(list)
    try:
        with open(csv_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                gene, og = parts[0], parts[1]
                og_members[og].append(gene)
    except OSError:
        pass
    return dict(og_members)


def load_possvm_trees(possvm_dir: Path) -> tuple[list, list]:
    """Return (tree_records, all_species).

    Each record::

        {id, hg, family, prefix, n_leaves, species, og_names, og_member_counts,
         tree_dict, ogs}

    ``tree_dict`` and ``ogs`` are the heavy data kept separate for lazy loading.
    """
    try:
        from ete3 import Tree  # type: ignore
    except ImportError:
        print("ERROR: ete3 not installed. Run: pip install ete3", file=sys.stderr)
        return [], []

    records = []
    all_species: set = set()

    for nwk in sorted(possvm_dir.glob("*.ortholog_groups.newick")):
        # Derive HG/family from filename
        stem = nwk.stem  # strips last .newick
        for suffix in (".treefile.ortholog_groups", ".ortholog_groups"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        parts = stem.split(".", 2)
        if len(parts) >= 3:
            prefix, family, hg = parts[0], parts[1], parts[2]
        elif len(parts) == 2:
            prefix, family, hg = parts[0], parts[1], stem
        else:
            prefix, family, hg = stem, "", stem

        try:
            t = Tree(str(nwk), format=1)
        except Exception as exc:
            print(f"WARN: cannot parse {nwk.name}: {exc}", file=sys.stderr)
            continue

        tree_dict = gene_tree_to_dict(t)
        leaves = t.get_leaves()
        species = sorted({get_species_prefix(l.name) for l in leaves if l.name})
        all_species.update(species)

        csv_path = nwk.parent / (stem + ".ortholog_groups.csv")
        if not csv_path.exists():
            csv_path = nwk.parent / nwk.name.replace(
                ".ortholog_groups.newick", ".ortholog_groups.csv"
            )
        ogs = load_og_csv(csv_path) if csv_path.exists() else {}

        records.append({
            "id":               stem,
            "hg":               hg,
            "family":           family,
            "prefix":           prefix,
            "n_leaves":         len(leaves),
            "species":          species,
            # lightweight index fields
            "og_names":         sorted(ogs.keys()),
            "n_ogs":            len(ogs),
            # heavy fields (for lazy loading)
            "tree_dict":        tree_dict,
            "ogs":              ogs,
        })

    return records, sorted(all_species)


# ── HTML template ─────────────────────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Step 2 Report – Gene Tree Explorer</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
       font-size: 12px; background: #f7f7f7; color: #333; overflow: hidden; }

/* ── Header ── */
#header { padding: 10px 16px; background: #2c3e50; color: #ecf0f1; display: flex;
          align-items: center; gap: 12px; flex-wrap: wrap; min-height: 46px; }
#header h2 { font-size: 15px; font-weight: 600; white-space: nowrap; }
#tree-count { font-size: 11px; color: #95a5a6; }

/* ── Layout ── */
#app { display: flex; height: calc(100vh - 46px); overflow: hidden; }

/* ── Sidebar ── */
#sidebar { flex: 0 0 220px; display: flex; flex-direction: column;
           border-right: 1px solid #ccc; background: #fff; }
#sidebar-top { padding: 8px; border-bottom: 1px solid #eee; }
#hg-search { width: 100%; padding: 5px 7px; font-size: 11px;
             border: 1px solid #ccc; border-radius: 3px; }
#hg-count { font-size: 10px; color: #999; margin-top: 4px; }
#hg-list { flex: 1; overflow-y: auto; }

/* Family accordion */
.fam-group {}
.fam-header { padding: 5px 8px; background: #ecf0f1; border-bottom: 1px solid #ddd;
              cursor: pointer; display: flex; align-items: center; gap: 5px;
              font-size: 11px; font-weight: 600; color: #444; user-select: none; }
.fam-header:hover { background: #dce4ec; }
.fam-arrow { font-size: 9px; transition: transform .15s; }
.fam-header.open .fam-arrow { transform: rotate(90deg); }
.fam-body { display: none; }
.fam-body.open { display: block; }

/* HG items */
.hg-item { padding: 5px 10px 5px 16px; cursor: pointer;
           border-bottom: 1px solid #f0f0f0; line-height: 1.35; }
.hg-item:hover { background: #f5f5f5; }
.hg-item.selected { background: #d5f5e3; border-left: 3px solid #1abc9c; padding-left: 13px; }
.hg-item .hg-name { font-weight: 600; font-size: 11px; }
.hg-item .hg-meta { font-size: 10px; color: #888; }

/* ── Main ── */
#main { flex: 1; display: flex; flex-direction: column; overflow: hidden; }

/* Controls bar */
#controls { padding: 5px 10px; background: #f5f5f5; border-bottom: 1px solid #ddd;
            display: flex; align-items: center; gap: 8px; flex-wrap: wrap; }
.ctrl-btn { padding: 3px 9px; border: 1px solid #bbb; border-radius: 3px;
            cursor: pointer; background: #fff; font-size: 11px; }
.ctrl-btn:hover { background: #eee; }
#tree-title { font-size: 11px; color: #555; margin-left: 4px; }
#n-ogs-label { font-size: 10px; color: #888; margin-left: auto; }

/* Species legend */
#species-legend { padding: 4px 10px; background: #fafafa;
                  border-bottom: 1px solid #eee;
                  display: flex; flex-wrap: wrap; gap: 8px; align-items: center;
                  font-size: 10px; min-height: 24px; }
.leg-item { display: flex; align-items: center; gap: 3px; }
.leg-dot { width: 9px; height: 9px; border-radius: 50%; flex-shrink: 0; }

/* Tree area */
#tree-wrap { flex: 1; overflow: hidden; position: relative; background: #fff; }
#tree-svg { width: 100%; height: 100%; cursor: grab; display: block; }
#tree-svg:active { cursor: grabbing; }

/* SVG elements */
.link { fill: none; stroke: #d0d0d0; stroke-width: 1px; }
.node-g circle { cursor: pointer; transition: r .12s; }
.node-g circle:hover { stroke-width: 2px !important; }
.leaf-label { font-size: 10px; fill: #333; pointer-events: none; }
.og-label { font-size: 9px; fill: #c0392b; pointer-events: none; font-style: italic; }

/* Tooltip */
#tooltip { position: fixed; display: none; pointer-events: none;
           background: rgba(255,255,255,.97); border: 1px solid #bbb;
           border-radius: 5px; padding: 8px 10px; font-size: 11px;
           box-shadow: 0 2px 8px rgba(0,0,0,.15); z-index: 100; max-width: 240px; }
#tt-name { font-weight: 700; margin-bottom: 4px; font-size: 12px; }
.tt-row { display: flex; justify-content: space-between; gap: 10px; color: #555;
          margin-top: 2px; }

#empty-msg { padding: 40px; color: #999; text-align: center; }
</style>
</head>
<body>

<div id="header">
  <h2>Step 2 – Gene Tree Explorer</h2>
  <span id="tree-count"></span>
</div>

<div id="app">
  <div id="sidebar">
    <div id="sidebar-top">
      <input id="hg-search" type="text" placeholder="Search HG / family…">
      <div id="hg-count"></div>
    </div>
    <div id="hg-list"></div>
  </div>

  <div id="main">
    <div id="controls">
      <button class="ctrl-btn" onclick="expandAll()">Expand all</button>
      <button class="ctrl-btn" onclick="collapseToOGs()">Collapse to OGs</button>
      <button class="ctrl-btn" onclick="collapseAll()">Collapse all</button>
      <span id="tree-title"></span>
      <span id="n-ogs-label"></span>
    </div>
    <div id="species-legend"></div>
    <div id="tree-wrap">
      <svg id="tree-svg"></svg>
    </div>
  </div>
</div>

<div id="tooltip">
  <div id="tt-name"></div>
  <div id="tt-body"></div>
</div>

<!-- ── Per-HG lazy data (inert JSON strings, parsed only on selection) ── -->
<div id="lazy-data" style="display:none">
%%LAZY_SCRIPTS%%
</div>

<script>
// ── Lightweight index (metadata only, no tree dicts) ─────────────────────────
const TREE_INDEX  = %%TREE_INDEX_JSON%%;
const ALL_SPECIES = %%SPECIES_JSON%%;

// ── Species colour scale ──────────────────────────────────────────────────────
const SP_COLORS = {};
const palette = [
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
  "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
  "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5"
];
ALL_SPECIES.forEach((sp, i) => { SP_COLORS[sp] = palette[i % palette.length]; });
function spColor(sp) { return SP_COLORS[sp] || "#aaa"; }

// ── State ─────────────────────────────────────────────────────────────────────
let currentIndex  = null;   // lightweight index record
let currentDetail = null;   // {tree, ogs} — loaded on demand
let rootNode = null;
let gMain    = null;
let _uid     = 0;

// ── Lazy loader ───────────────────────────────────────────────────────────────
function loadDetail(id) {
  const el = document.getElementById("treedata-" + id);
  if (!el) return null;
  try {
    return JSON.parse(el.textContent);
  } catch (e) {
    console.error("Failed to parse tree data for", id, e);
    return null;
  }
}

// ── Sidebar ───────────────────────────────────────────────────────────────────
document.getElementById("tree-count").textContent =
  TREE_INDEX.length + " gene tree" + (TREE_INDEX.length !== 1 ? "s" : "");

// Group index by family
function groupByFamily(records) {
  const groups = {};
  for (const r of records) {
    const f = r.family || "(other)";
    if (!groups[f]) groups[f] = [];
    groups[f].push(r);
  }
  return groups;
}

function renderSidebar(filter) {
  const lc = (filter || "").toLowerCase().trim();
  const list = document.getElementById("hg-list");
  list.innerHTML = "";

  const groups = groupByFamily(TREE_INDEX);
  let totalVisible = 0;

  for (const [fam, recs] of Object.entries(groups).sort()) {
    const matching = lc
      ? recs.filter(r => r.hg.toLowerCase().includes(lc) || r.family.toLowerCase().includes(lc))
      : recs;
    if (matching.length === 0) continue;

    totalVisible += matching.length;
    const forceOpen = lc.length > 0;

    const groupDiv = document.createElement("div");
    groupDiv.className = "fam-group";

    const hdr = document.createElement("div");
    hdr.className = "fam-header" + (forceOpen ? " open" : "");
    hdr.dataset.fam = fam;
    hdr.innerHTML =
      `<span class="fam-arrow">▶</span>` +
      `<span>${fam}</span>` +
      `<span style="font-weight:400;color:#888;margin-left:auto">${matching.length}</span>`;
    hdr.addEventListener("click", () => toggleFamily(hdr));
    groupDiv.appendChild(hdr);

    const body = document.createElement("div");
    body.className = "fam-body" + (forceOpen ? " open" : "");

    for (const rec of matching) {
      const item = document.createElement("div");
      item.className = "hg-item" +
        (currentIndex && currentIndex.id === rec.id ? " selected" : "");
      item.innerHTML =
        `<div class="hg-name">${rec.hg}</div>` +
        `<div class="hg-meta">${rec.n_leaves} genes · ${rec.n_ogs} OGs</div>`;
      item.addEventListener("click", () => selectTree(rec));
      body.appendChild(item);
    }

    groupDiv.appendChild(body);
    list.appendChild(groupDiv);
  }

  document.getElementById("hg-count").textContent = totalVisible + " shown";
}

function toggleFamily(hdr) {
  hdr.classList.toggle("open");
  const body = hdr.nextElementSibling;
  body.classList.toggle("open", hdr.classList.contains("open"));
}

document.getElementById("hg-search").addEventListener("input", function () {
  renderSidebar(this.value);
});

// ── Species legend ────────────────────────────────────────────────────────────
function buildLegend(species) {
  const el = document.getElementById("species-legend");
  el.innerHTML = "";
  species.forEach(sp => {
    const item = document.createElement("div");
    item.className = "leg-item";
    item.innerHTML =
      `<div class="leg-dot" style="background:${spColor(sp)}"></div>${sp}`;
    el.appendChild(item);
  });
}

// ── Tree selection (lazy load) ────────────────────────────────────────────────
function selectTree(rec) {
  currentIndex  = rec;
  currentDetail = loadDetail(rec.id);
  if (!currentDetail) {
    console.warn("No detail data for", rec.id);
    return;
  }
  renderSidebar(document.getElementById("hg-search").value);
  buildLegend(rec.species);
  document.getElementById("tree-title").textContent =
    rec.id + " · " + rec.n_leaves + " genes";
  document.getElementById("n-ogs-label").textContent =
    rec.n_ogs + " orthogroups";
  drawTree(currentDetail.tree);
}

// ── D3 tree rendering ─────────────────────────────────────────────────────────
const svg     = d3.select("#tree-svg");
const tooltip = document.getElementById("tooltip");

function isOGNode(d) {
  return !d.data.leaf && d.data.name && d.data.name !== "";
}

function countDescLeaves(children) {
  if (!children) return 0;
  let n = 0;
  for (const c of children) {
    if (c.data.leaf) n++;
    else n += countDescLeaves(c.children) + countDescLeaves(c._children);
  }
  return n;
}

function drawTree(treeData) {
  svg.selectAll("*").remove();
  _uid = 0;

  const wrap = document.getElementById("tree-wrap");
  const W = wrap.clientWidth  || 800;
  const H = wrap.clientHeight || 600;
  svg.attr("width", W).attr("height", H);

  const zoom = d3.zoom().scaleExtent([0.05, 20]).on("zoom", e => {
    gMain.attr("transform", e.transform);
  });
  svg.call(zoom).on("dblclick.zoom", null);
  gMain = svg.append("g");

  rootNode = d3.hierarchy(treeData, d => d.children);
  rootNode.each(d => { d._uid = ++_uid; });

  renderTree();
  fitView(W, H, zoom);
}

function renderTree() {
  const wrap = document.getElementById("tree-wrap");
  const W = wrap.clientWidth  || 800;
  const H = wrap.clientHeight || 600;
  const mg = { top: 20, right: 200, bottom: 20, left: 30 };
  const iW = W - mg.left - mg.right;
  const iH = H - mg.top  - mg.bottom;

  const nVis = rootNode.leaves().length;
  const rowH = Math.max(13, Math.min(24, Math.floor(iH / Math.max(nVis, 1))));
  const tH   = Math.max(iH, nVis * rowH);

  const layout = d3.cluster().size([tH, iW]);
  layout(rootNode);

  // Links
  gMain.selectAll(".link")
    .data(rootNode.links(), d => d.target._uid)
    .join("path")
    .attr("class", "link")
    .attr("d", d =>
      `M${d.source.y + mg.left},${d.source.x + mg.top}` +
      `H${d.target.y + mg.left}` +
      `V${d.target.x + mg.top}`);

  // Nodes
  const nodeSel = gMain.selectAll(".node-g")
    .data(rootNode.descendants(), d => d._uid);

  const entered = nodeSel.enter().append("g").attr("class", "node-g");
  entered.append("circle");
  entered.append("text").attr("class", "leaf-label");
  entered.append("text").attr("class", "og-label");
  nodeSel.exit().remove();

  const merged = nodeSel.merge(entered)
    .attr("transform", d => `translate(${d.y + mg.left},${d.x + mg.top})`);

  merged.select("circle")
    .attr("r", d => d.data.leaf ? 3.5 : d._children ? 7 : isOGNode(d) ? 5 : 3)
    .attr("fill", d => {
      if (d.data.leaf)  return spColor(d.data.species || "");
      if (d._children)  return "#e67e22";
      if (isOGNode(d))  return "#e74c3c";
      return "#bbb";
    })
    .attr("stroke", d =>
      (d._children || isOGNode(d)) ? "#c0392b" : "#888")
    .attr("stroke-width", d =>
      (d._children || isOGNode(d)) ? 1.5 : 0.8)
    .on("click", (event, d) => {
      if (d.data.leaf) return;
      event.stopPropagation();
      if (d.children) { d._children = d.children; d.children = null; }
      else if (d._children) { d.children = d._children; d._children = null; }
      renderTree();
    })
    .on("mouseover", showTip)
    .on("mousemove", moveTip)
    .on("mouseout",  hideTip);

  // Leaf labels
  merged.select(".leaf-label")
    .attr("x", d => d.data.leaf ? 6 : 0)
    .attr("dy", "0.32em")
    .attr("text-anchor", "start")
    .attr("font-size", d => d.data.leaf ? Math.min(11, rowH - 1) : 0)
    .attr("fill", d => spColor(d.data.species || ""))
    .attr("display", d => d.data.leaf ? null : "none")
    .text(d => d.data.leaf ? d.data.name : "");

  // OG / collapsed labels
  merged.select(".og-label")
    .attr("x", d => d._children ? 9 : -7)
    .attr("dy", "0.32em")
    .attr("text-anchor", d => d._children ? "start" : "end")
    .attr("font-size", 9)
    .attr("display", d =>
      (!d.data.leaf && (isOGNode(d) || d._children)) ? null : "none")
    .text(d => {
      if (d._children) {
        const n = countDescLeaves(d._children);
        return (d.data.name ? d.data.name + " " : "") + `[${n}]`;
      }
      return isOGNode(d) ? d.data.name : "";
    });
}

function fitView(W, H, zoom) {
  const b = gMain.node().getBBox();
  if (!b.width || !b.height) return;
  const scale = Math.min((W - 10) / b.width, (H - 10) / b.height, 1.5);
  const tx = (W - b.width  * scale) / 2 - b.x * scale;
  const ty = (H - b.height * scale) / 2 - b.y * scale;
  svg.call(zoom.transform, d3.zoomIdentity.translate(tx, ty).scale(scale));
}

// ── Tree controls ─────────────────────────────────────────────────────────────
function expandAll() {
  if (!rootNode) return;
  rootNode.each(d => {
    if (d._children) { d.children = d._children; d._children = null; }
  });
  renderTree();
}

function collapseToOGs() {
  if (!rootNode) return;
  rootNode.each(d => {
    if (d.data.leaf || !isOGNode(d)) return;
    (d.children || []).forEach(c => {
      if (!c.data.leaf && c.children) { c._children = c.children; c.children = null; }
    });
  });
  renderTree();
}

function collapseAll() {
  if (!rootNode) return;
  rootNode.each(d => {
    if (d.depth > 0 && !d.data.leaf && d.children) {
      d._children = d.children; d.children = null;
    }
  });
  if (rootNode._children) { rootNode.children = rootNode._children; rootNode._children = null; }
  renderTree();
}

// ── Tooltip ───────────────────────────────────────────────────────────────────
function showTip(event, d) {
  const name = d.data.name || (d.data.leaf ? "leaf" : "internal");
  document.getElementById("tt-name").textContent = name;
  let body = "";
  if (d.data.leaf) {
    body += `<div class="tt-row"><span>Species</span>` +
      `<strong style="color:${spColor(d.data.species)}">${d.data.species || "?"}</strong></div>`;
    if (currentDetail) {
      for (const [og, genes] of Object.entries(currentDetail.ogs)) {
        if (genes.includes(d.data.name)) {
          body += `<div class="tt-row"><span>OG</span><strong>${og}</strong></div>`;
          break;
        }
      }
    }
  } else {
    const nL = d._children ? countDescLeaves(d._children) : d.leaves().length;
    body += `<div class="tt-row"><span>Subtree leaves</span><strong>${nL}</strong></div>`;
    body += `<div class="tt-row" style="color:#888"><span>${d._children ? "collapsed" : "expanded"}</span>` +
      `<span>click to ${d._children ? "expand" : "collapse"}</span></div>`;
    if (isOGNode(d) && currentDetail && currentDetail.ogs[d.data.name]) {
      body += `<div class="tt-row"><span>OG members</span>` +
        `<strong>${currentDetail.ogs[d.data.name].length}</strong></div>`;
    }
  }
  document.getElementById("tt-body").innerHTML = body;
  tooltip.style.display = "block";
  moveTip(event);
}

function moveTip(event) {
  const x = event.clientX + 14, y = event.clientY - 10;
  tooltip.style.left = Math.min(x, window.innerWidth  - tooltip.offsetWidth  - 5) + "px";
  tooltip.style.top  = Math.min(y, window.innerHeight - tooltip.offsetHeight - 5) + "px";
}
function hideTip() { tooltip.style.display = "none"; }

// ── Init ──────────────────────────────────────────────────────────────────────
renderSidebar("");
if (TREE_INDEX.length > 0) {
  selectTree(TREE_INDEX[0]);
} else {
  document.getElementById("tree-wrap").innerHTML =
    '<div id="empty-msg">No gene trees found.<br>' +
    'Run POSSVM (step 2) and pass <code>--possvm_dir</code>.</div>';
}
</script>
</body>
</html>
"""


# ── Argument parsing and main ─────────────────────────────────────────────────

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate interactive HTML gene-tree report for step2 (POSSVM) outputs.")
    p.add_argument("--possvm_dir", default="results/possvm",
                   help="Directory with POSSVM *.ortholog_groups.newick files")
    p.add_argument("--output", required=True, help="Output HTML file path")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    possvm_dir = Path(args.possvm_dir)
    if not possvm_dir.exists():
        print(f"WARN: {possvm_dir} does not exist – generating empty report.", file=sys.stderr)

    records, all_species = load_possvm_trees(possvm_dir)
    print(f"Loaded {len(records)} gene trees, {len(all_species)} species.", file=sys.stderr)

    # Build lightweight index (no tree/ogs dicts)
    index_records = [
        {k: v for k, v in rec.items() if k not in ("tree_dict", "ogs")}
        for rec in records
    ]

    # Build per-HG lazy <script> tags
    lazy_parts = []
    for rec in records:
        detail = {"tree": rec["tree_dict"], "ogs": rec["ogs"]}
        tag_id = _html.escape(rec["id"], quote=True)
        lazy_parts.append(
            f'<script type="application/json" id="treedata-{tag_id}">'
            + json.dumps(detail, separators=(",", ":"))
            + "</script>"
        )
    lazy_scripts = "\n".join(lazy_parts)

    html = (HTML_TEMPLATE
            .replace("%%LAZY_SCRIPTS%%",    lazy_scripts)
            .replace("%%TREE_INDEX_JSON%%", json.dumps(index_records))
            .replace("%%SPECIES_JSON%%",    json.dumps(all_species)))

    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
