#!/usr/bin/env python3
"""
Generate a self-contained interactive HTML tree visualization of ancestral
state reconstruction results.

Features:
  - Collapsible D3.js species tree (click internal nodes to collapse/expand)
  - Nodes coloured by P(present) for the selected orthogroup
    (diverging red → white → blue colour scale)
  - Leaf nodes coloured by observed PAM state (if --pam is provided)
  - Left sidebar: searchable OG list grouped by gene-family class and HG
  - Hover tooltip: node name, P(present), support classification
  - Zoom / pan on the tree panel
  - Colour legend

Usage:
    python visualize_ancestry.py \\
        --tree       Bilateria.pruned.tree \\
        --node_probs Bilateria.node_probs.tsv \\
        --states     Bilateria.ancestral_states.tsv \\
        --pam        Bilateria.pam.tsv \\
        --output     Bilateria.html \\
        --node       Bilateria
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from ete3 import Tree


# ── Tree conversion ───────────────────────────────────────────────────────────

def tree_to_dict(node: Tree) -> dict:
    """Recursively convert an ete3 node to a nested dict for D3.hierarchy."""
    d: dict = {"name": node.name or "", "dist": float(node.dist)}
    if node.is_leaf():
        d["leaf"] = True
    else:
        d["children"] = [tree_to_dict(c) for c in node.children]
    return d


# ── OG name parsing ───────────────────────────────────────────────────────────

def parse_og_name(og: str) -> dict:
    """Split 'tfs.HG001.OG1' → {prefix, hg, og_num}."""
    parts = og.split(".")
    return {
        "prefix": parts[0] if len(parts) > 0 else og,
        "hg":     parts[1] if len(parts) > 1 else "",
        "og_num": parts[2] if len(parts) > 2 else "",
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate interactive HTML ancestral-reconstruction tree."
    )
    parser.add_argument("--tree",       required=True, help="Pruned species tree (newick)")
    parser.add_argument("--node_probs", required=True, help="Per-OG per-node probabilities TSV")
    parser.add_argument("--states",     required=True, help="Per-OG ancestral states TSV")
    parser.add_argument("--pam",        default=None,  help="PAM TSV for observed leaf states (optional)")
    parser.add_argument("--output",     required=True, help="Output HTML file")
    parser.add_argument("--node",       default="root", help="Clade name for the title")
    args = parser.parse_args()

    # ── Load tree ─────────────────────────────────────────────────────────────
    t = Tree(args.tree, format=1)
    tree_dict = tree_to_dict(t)

    # ── Load node_probs → {node_name: {og: p}} ───────────────────────────────
    np_df = pd.read_csv(args.node_probs, sep="\t")
    node_probs: dict[str, dict] = {}
    for _, row in np_df.iterrows():
        node_probs.setdefault(str(row["node"]), {})[str(row["og"])] = \
            None if pd.isna(row["P_present"]) else round(float(row["P_present"]), 4)

    # ── Load ancestral states → {og: meta} ───────────────────────────────────
    states_df = pd.read_csv(args.states, sep="\t")
    og_meta: dict[str, dict] = {}
    for _, row in states_df.iterrows():
        og = str(row["og"])
        p  = None if pd.isna(row["P_at_root"]) else round(float(row["P_at_root"]), 4)
        og_meta[og] = {
            "n": int(row["n_present"]),
            "t": int(row["n_total"]),
            "s": str(row["support"]),
            "p": p,
            **parse_og_name(og),
        }

    # ── Load PAM for observed leaf states ─────────────────────────────────────
    leaf_states: dict[str, dict] = {}
    if args.pam and Path(args.pam).exists():
        pam_df = pd.read_csv(args.pam, sep="\t", index_col=0)
        for sps in pam_df.index:
            leaf_states[str(sps)] = {
                og: int(pam_df.loc[sps, og]) for og in pam_df.columns
            }

    # ── Inject into HTML template ─────────────────────────────────────────────
    html = HTML_TEMPLATE \
        .replace("%%NODE_NAME%%",       args.node) \
        .replace("%%TREE_JSON%%",        json.dumps(tree_dict, separators=(",", ":"))) \
        .replace("%%NODE_PROBS_JSON%%",  json.dumps(node_probs, separators=(",", ":"))) \
        .replace("%%OG_META_JSON%%",     json.dumps(og_meta,    separators=(",", ":"))) \
        .replace("%%LEAF_STATES_JSON%%", json.dumps(leaf_states, separators=(",", ":")))

    Path(args.output).write_text(html)
    print(f"Written: {args.output}", file=sys.stderr)


# ── HTML / D3 template ────────────────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Ancestral Reconstruction — %%NODE_NAME%%</title>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body { display: flex; height: 100vh; font-family: "Helvetica Neue", Arial, sans-serif;
       font-size: 13px; background: #f8f8f8; color: #333; overflow: hidden; }

/* ── Sidebar ── */
#sidebar {
  width: 300px; min-width: 260px; display: flex; flex-direction: column;
  background: #fff; border-right: 1px solid #ddd; z-index: 10;
}
#sidebar-header { padding: 14px 14px 10px; background: #1e3a5f; color: #fff; }
#sidebar-header h2 { font-size: 15px; font-weight: 600; }
#sidebar-header .sub { font-size: 11px; opacity: 0.75; margin-top: 3px; }
#controls { padding: 8px 10px; border-bottom: 1px solid #eee; }
#search { width: 100%; padding: 6px 8px; border: 1px solid #ccc; border-radius: 4px;
          font-size: 12px; }
#og-count { font-size: 11px; color: #888; margin-top: 5px; }
#og-list { flex: 1; overflow-y: auto; padding: 4px 0; }

.og-group { }
.og-group-header {
  display: flex; align-items: center; padding: 5px 10px;
  cursor: pointer; font-weight: 600; font-size: 12px;
  background: #f4f4f4; border-bottom: 1px solid #eee;
  user-select: none;
}
.og-group-header:hover { background: #ececec; }
.og-group-header .arrow { margin-right: 6px; font-size: 10px; transition: transform .15s; }
.og-group-header.collapsed .arrow { transform: rotate(-90deg); }
.og-group-body { }

.hg-group-header {
  padding: 3px 18px; font-size: 11px; color: #666; cursor: pointer;
  background: #fafafa; border-bottom: 1px solid #f0f0f0; user-select: none;
}
.hg-group-header:hover { background: #f0f0f0; }
.hg-group-body { }

.og-item {
  display: flex; align-items: center; padding: 4px 24px 4px 28px;
  cursor: pointer; border-bottom: 1px solid #f5f5f5;
}
.og-item:hover { background: #f0f0f0; }
.og-item.selected { background: #dbeaff; }
.og-label { flex: 1; font-size: 11px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.og-frac { font-size: 10px; color: #999; margin-right: 5px; white-space: nowrap; }
.badge {
  font-size: 9px; padding: 1px 5px; border-radius: 8px; color: #fff; white-space: nowrap;
}
.badge-present       { background: #2e7d32; }
.badge-likely_present{ background: #66bb6a; }
.badge-likely_absent { background: #f57c00; }
.badge-absent        { background: #c62828; }
.badge-uncertain     { background: #9e9e9e; }
.badge-failed        { background: #9e9e9e; }

/* ── Main tree panel ── */
#main { flex: 1; display: flex; flex-direction: column; overflow: hidden; }
#top-bar {
  display: flex; align-items: center; padding: 6px 14px;
  background: #fff; border-bottom: 1px solid #ddd; gap: 12px;
}
#selected-og-label { font-size: 12px; color: #555; flex: 1; }
#legend-bar { display: flex; align-items: center; gap: 8px; font-size: 11px; color: #555; }
#legend-gradient {
  width: 120px; height: 12px; border-radius: 2px;
  background: linear-gradient(to right, #b2182b, #f7f7f7, #2166ac);
  border: 1px solid #ccc;
}
#tree-wrap { flex: 1; overflow: hidden; position: relative; }
#tree-svg { width: 100%; height: 100%; cursor: grab; }
#tree-svg:active { cursor: grabbing; }

/* ── Tooltip ── */
#tooltip {
  position: fixed; display: none; pointer-events: none;
  background: rgba(255,255,255,0.97); border: 1px solid #bbb;
  border-radius: 5px; padding: 8px 10px; font-size: 12px;
  box-shadow: 0 2px 8px rgba(0,0,0,.15); z-index: 100; max-width: 220px;
}
#tooltip .tt-name { font-weight: 700; margin-bottom: 4px; font-size: 13px; }
#tooltip .tt-row  { display: flex; justify-content: space-between; gap: 12px; color: #555; }

/* ── Tree SVG elements ── */
.link { fill: none; stroke: #bbb; stroke-width: 1.2px; }
.node-circle { stroke: #555; stroke-width: 1px; cursor: pointer; transition: r .1s; }
.node-circle.leaf { stroke-width: 0.8px; }
.node-circle:hover { stroke: #000; stroke-width: 2px; }
.node-label { font-size: 10px; fill: #444; pointer-events: none; }
.node-label.internal { fill: #888; font-style: italic; }
.node-collapsed circle { stroke-dasharray: 3,2; }
</style>
</head>
<body>

<!-- Sidebar -->
<div id="sidebar">
  <div id="sidebar-header">
    <h2>Ancestral Reconstruction</h2>
    <div class="sub">Clade: <strong>%%NODE_NAME%%</strong></div>
  </div>
  <div id="controls">
    <input id="search" type="text" placeholder="Search orthogroups…" oninput="filterOGs()">
    <div id="og-count"></div>
  </div>
  <div id="og-list"></div>
</div>

<!-- Main -->
<div id="main">
  <div id="top-bar">
    <span id="selected-og-label">← Select an orthogroup</span>
    <div id="legend-bar">
      <span>absent</span>
      <div id="legend-gradient"></div>
      <span>present</span>
    </div>
  </div>
  <div id="tree-wrap">
    <svg id="tree-svg"></svg>
  </div>
</div>

<div id="tooltip">
  <div class="tt-name" id="tt-name"></div>
  <div id="tt-body"></div>
</div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
// ── Embedded data ─────────────────────────────────────────────────────────────
const TREE_DATA   = %%TREE_JSON%%;
const NODE_PROBS  = %%NODE_PROBS_JSON%%;
const OG_META     = %%OG_META_JSON%%;
const LEAF_STATES = %%LEAF_STATES_JSON%%;

// ── Colour scale: red (absent=0) → white (0.5) → blue (present=1) ─────────────
const colorScale = d3.scaleSequential(d3.interpolateRdBu).domain([0, 1]);
const GREY = "#ccc";

function supportBadgeClass(s) {
  return "badge badge-" + (s || "uncertain");
}

// ── State ─────────────────────────────────────────────────────────────────────
let selectedOG = null;

// ── Sidebar ───────────────────────────────────────────────────────────────────
function buildSidebar() {
  // Group OGs: prefix → hg → [og]
  const groups = {};
  for (const [og, meta] of Object.entries(OG_META)) {
    const pref = meta.prefix || "other";
    const hg   = meta.hg    || "—";
    if (!groups[pref])     groups[pref] = {};
    if (!groups[pref][hg]) groups[pref][hg] = [];
    groups[pref][hg].push(og);
  }

  const container = document.getElementById("og-list");
  container.innerHTML = "";

  for (const [pref, hgs] of Object.entries(groups).sort()) {
    const groupDiv = document.createElement("div");
    groupDiv.className = "og-group";
    groupDiv.dataset.prefix = pref;

    const hdr = document.createElement("div");
    hdr.className = "og-group-header";
    const totalOGs = Object.values(hgs).flat().length;
    hdr.innerHTML = `<span class="arrow">▾</span> ${pref} <span style="font-weight:400;color:#999;margin-left:6px">(${totalOGs})</span>`;
    hdr.onclick = () => {
      const body = groupDiv.querySelector(".og-group-body");
      const collapsed = hdr.classList.toggle("collapsed");
      body.style.display = collapsed ? "none" : "";
    };
    groupDiv.appendChild(hdr);

    const body = document.createElement("div");
    body.className = "og-group-body";

    for (const [hg, ogs] of Object.entries(hgs).sort()) {
      const hgDiv = document.createElement("div");
      const hgHdr = document.createElement("div");
      hgHdr.className = "hg-group-header";
      hgHdr.textContent = `  ${hg} (${ogs.length})`;
      hgHdr.onclick = () => {
        const b = hgDiv.querySelector(".hg-group-body");
        b.style.display = b.style.display === "none" ? "" : "none";
      };
      hgDiv.appendChild(hgHdr);

      const hgBody = document.createElement("div");
      hgBody.className = "hg-group-body";

      for (const og of ogs.sort()) {
        const meta = OG_META[og];
        const item = document.createElement("div");
        item.className = "og-item";
        item.dataset.og = og;
        item.innerHTML = `
          <span class="og-label" title="${og}">${meta.og_num || og}</span>
          <span class="og-frac">${meta.n}/${meta.t}</span>
          <span class="${supportBadgeClass(meta.s)}">${(meta.p !== null && meta.p !== undefined) ? meta.p.toFixed(2) : "—"}</span>`;
        item.onclick = () => selectOG(og);
        hgBody.appendChild(item);
      }
      hgDiv.appendChild(hgBody);
      body.appendChild(hgDiv);
    }

    groupDiv.appendChild(body);
    container.appendChild(groupDiv);
  }
  updateOGCount();
}

function updateOGCount() {
  const visible = document.querySelectorAll(".og-item:not([style*='display: none'])").length;
  document.getElementById("og-count").textContent = `${visible} orthogroups`;
}

function filterOGs() {
  const q = document.getElementById("search").value.toLowerCase().trim();
  document.querySelectorAll(".og-item").forEach(el => {
    const match = !q || el.dataset.og.toLowerCase().includes(q);
    el.style.display = match ? "" : "none";
  });
  // Show/hide HG headers if all children hidden
  document.querySelectorAll(".hg-group-body").forEach(body => {
    const anyVisible = [...body.querySelectorAll(".og-item")].some(e => e.style.display !== "none");
    body.parentElement.style.display = anyVisible ? "" : "none";
  });
  // Show/hide prefix group headers
  document.querySelectorAll(".og-group").forEach(grp => {
    const anyVisible = [...grp.querySelectorAll(".og-item")].some(e => e.style.display !== "none");
    grp.style.display = anyVisible ? "" : "none";
  });
  updateOGCount();
}

function selectOG(og) {
  selectedOG = og;
  // Highlight selected item
  document.querySelectorAll(".og-item").forEach(el => {
    el.classList.toggle("selected", el.dataset.og === og);
  });
  const meta = OG_META[og];
  document.getElementById("selected-og-label").innerHTML =
    `<strong>${og}</strong> &nbsp;·&nbsp; ${meta.n}/${meta.t} species &nbsp;·&nbsp;
     <span class="${supportBadgeClass(meta.s)}" style="font-size:11px">${meta.s} (${meta.p !== null ? meta.p : "—"})</span>`;
  updateColors();
}

// ── D3 Tree ───────────────────────────────────────────────────────────────────
const svg   = d3.select("#tree-svg");
const wrap  = document.getElementById("tree-wrap");
let gMain; // main <g> for zoom

function getNodeColor(d) {
  if (!selectedOG) return GREY;
  if (d.data.leaf) {
    // Observed leaf state from PAM
    const states = LEAF_STATES[d.data.name];
    if (!states || !(selectedOG in states)) return "#e0e0e0";
    return colorScale(states[selectedOG]);        // 0→red, 1→blue
  } else {
    // Internal node: P(present) from PastML
    const nodeMap = NODE_PROBS[d.data.name];
    if (!nodeMap || !(selectedOG in nodeMap) || nodeMap[selectedOG] === null) return GREY;
    return colorScale(nodeMap[selectedOG]);
  }
}

function getNodeRadius(d) { return d.data.leaf ? 4 : 7; }

let nodeIdCounter = 0;
let root;

function initTree() {
  const W = wrap.clientWidth;
  const H = wrap.clientHeight;
  const margin = { top: 20, right: 160, bottom: 20, left: 40 };
  const innerW = W - margin.left - margin.right;
  const innerH = H - margin.top  - margin.bottom;

  svg.attr("width", W).attr("height", H);

  // Zoom
  const zoom = d3.zoom().scaleExtent([0.1, 10]).on("zoom", e => {
    gMain.attr("transform", e.transform);
  });
  svg.call(zoom);
  svg.on("dblclick.zoom", null);

  gMain = svg.append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  // Build hierarchy
  root = d3.hierarchy(TREE_DATA);
  root.each(d => { d.id = ++nodeIdCounter; });

  // Collapse all internal nodes except root by default
  // (leave top 2 levels open)
  root.each(d => {
    if (d.depth > 1 && !d.data.leaf) {
      d._children = d.children;
      d.children  = null;
    }
  });

  const treeLayout = d3.tree().size([innerH, innerW]);
  renderTree(treeLayout, gMain, W, H);
}

function renderTree(treeLayout, g, W, H) {
  const treeData = treeLayout(root);

  // ── Links ─────────────────────────────────────────────────────────────────
  const linkSel = g.selectAll(".link")
    .data(treeData.links(), d => d.target.id);

  linkSel.enter().append("path")
    .attr("class", "link")
    .merge(linkSel)
    .attr("d", d3.linkHorizontal().x(d => d.y).y(d => d.x));

  linkSel.exit().remove();

  // ── Nodes ─────────────────────────────────────────────────────────────────
  const nodeSel = g.selectAll(".node-g")
    .data(treeData.descendants(), d => d.id);

  const nodeEnter = nodeSel.enter().append("g")
    .attr("class", "node-g")
    .attr("transform", d => `translate(${d.y},${d.x})`);

  nodeEnter.append("circle")
    .attr("class", d => "node-circle" + (d.data.leaf ? " leaf" : ""))
    .attr("r", getNodeRadius)
    .attr("fill", getNodeColor)
    .on("click", (event, d) => {
      if (!d.data.leaf) toggleNode(d, g, treeLayout, W, H);
    })
    .on("mouseover", showTooltip)
    .on("mousemove", moveTooltip)
    .on("mouseout",  hideTooltip);

  // Leaf labels
  nodeEnter.filter(d => d.data.leaf).append("text")
    .attr("class", "node-label")
    .attr("dy", "0.31em")
    .attr("x", 9)
    .text(d => d.data.name);

  // Internal node label (shown on non-collapsed named nodes near root)
  nodeEnter.filter(d => !d.data.leaf && d.data.name && d.depth <= 3).append("text")
    .attr("class", "node-label internal")
    .attr("dy", "-0.6em")
    .attr("x", -4)
    .attr("text-anchor", "end")
    .text(d => d.data.name);

  const nodeUpdate = nodeSel.merge(nodeEnter)
    .attr("transform", d => `translate(${d.y},${d.x})`);

  nodeUpdate.select("circle")
    .attr("r", getNodeRadius)
    .attr("fill", getNodeColor);

  nodeSel.exit().remove();
}

function toggleNode(d, g, treeLayout, W, H) {
  if (d.children) {
    d._children = d.children;
    d.children  = null;
  } else if (d._children) {
    d.children  = d._children;
    d._children = null;
  }
  renderTree(treeLayout, g, W, H);
}

function updateColors() {
  svg.selectAll(".node-circle").attr("fill", getNodeColor);
}

// ── Tooltip ───────────────────────────────────────────────────────────────────
const tooltip = document.getElementById("tooltip");

function showTooltip(event, d) {
  const name = d.data.name || (d.data.leaf ? "leaf" : "unnamed node");
  document.getElementById("tt-name").textContent = name;

  let body = "";
  if (selectedOG) {
    if (d.data.leaf) {
      const states = LEAF_STATES[d.data.name];
      const s = states && (selectedOG in states) ? states[selectedOG] : null;
      body += `<div class="tt-row"><span>Observed state</span><span>${s !== null ? (s ? "present" : "absent") : "n/a"}</span></div>`;
    } else {
      const nodeMap = NODE_PROBS[d.data.name];
      const p = nodeMap && (selectedOG in nodeMap) ? nodeMap[selectedOG] : null;
      body += `<div class="tt-row"><span>P(present)</span><strong>${p !== null && p !== undefined ? p.toFixed(3) : "—"}</strong></div>`;
    }
    body += `<div class="tt-row" style="color:#aaa;font-size:10px">${selectedOG}</div>`;
  } else {
    body = `<div style="color:#aaa;font-size:11px">Select an OG in the sidebar</div>`;
  }

  // Add clade size for internal nodes
  if (!d.data.leaf) {
    const leaves = d.leaves ? d.leaves().length : "?";
    body += `<div class="tt-row" style="margin-top:4px"><span>Clade size</span><span>${typeof leaves === "number" ? leaves + " spp." : leaves}</span></div>`;
  }

  document.getElementById("tt-body").innerHTML = body;
  tooltip.style.display = "block";
  moveTooltip(event);
}

function moveTooltip(event) {
  const x = event.clientX + 14;
  const y = event.clientY - 10;
  const w = tooltip.offsetWidth;
  const h = tooltip.offsetHeight;
  tooltip.style.left = (x + w > window.innerWidth  ? x - w - 28 : x) + "px";
  tooltip.style.top  = (y + h > window.innerHeight ? y - h      : y) + "px";
}

function hideTooltip() { tooltip.style.display = "none"; }

// ── Init ──────────────────────────────────────────────────────────────────────
buildSidebar();
initTree();

// Auto-select first OG
const firstOG = Object.keys(OG_META)[0];
if (firstOG) selectOG(firstOG);
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
