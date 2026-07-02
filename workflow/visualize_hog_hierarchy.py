#!/usr/bin/env python3
"""
Generate a self-contained interactive HTML Sankey visualization of hierarchical
orthogroup (hOG) relationships across multiple taxonomic levels.

For each homology group (HG), the Sankey shows how orthogroups (OGs) defined
at a broad taxonomic level (e.g. Metazoa) split or expand into OGs at narrower
levels (e.g. Bilateria → Vertebrata). Splitting/expansion indicates gene family
evolution events such as whole-genome duplications.

Inputs:
  --links  One or more *.og_links.tsv files (one per HG, from link_hog_levels.py)
  --stats  One or more *.og_stats.tsv files (one per HG, from link_hog_levels.py)
  --trees  One or more pruned species tree files ({level}.pruned.tree)
  --levels Ordered clade names, comma-separated (broad → narrow)
  --output Output HTML file

This module is also imported by build_hog_report.py, which reuses HTML_TEMPLATE,
build_data() and parse_newick_tree() so that the two entry points share one UI.
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


# ── Newick parsing (stdlib only) ───────────────────────────────────────────────

class _NewickReader:
    """Minimal recursive-descent Newick parser → JSON-serialisable nested dict.

    Each node is ``{"name": str, "dist": float, "leaf": True}`` for leaves or
    ``{"name": str, "dist": float, "children": [...]}`` for internal nodes.
    The ``dist`` key is omitted when no branch length is present.
    """

    def __init__(self, text: str):
        self.s = text
        self.i = 0
        self.n = len(text)

    def _peek(self) -> str:
        return self.s[self.i] if self.i < self.n else ""

    def read_node(self) -> dict:
        node: dict = {}
        if self._peek() == "(":
            self.i += 1  # consume '('
            children = [self.read_node()]
            while self._peek() == ",":
                self.i += 1  # consume ','
                children.append(self.read_node())
            if self._peek() == ")":
                self.i += 1  # consume ')'
            node["children"] = children
        name = self._read_label()
        if name:
            node["name"] = name
        else:
            node["name"] = ""
        if self._peek() == ":":
            self.i += 1  # consume ':'
            dist = self._read_number()
            if dist is not None:
                node["dist"] = dist
        if "children" not in node:
            node["leaf"] = True
        return node

    def _read_label(self) -> str:
        # Quoted label
        if self._peek() == "'":
            self.i += 1
            start = self.i
            while self.i < self.n and self.s[self.i] != "'":
                self.i += 1
            label = self.s[start:self.i]
            if self.i < self.n:
                self.i += 1  # consume closing quote
            return label.strip()
        start = self.i
        while self.i < self.n and self.s[self.i] not in ",():;":
            self.i += 1
        return self.s[start:self.i].strip()

    def _read_number(self):
        start = self.i
        while self.i < self.n and (self.s[self.i].isdigit()
                                   or self.s[self.i] in ".+-eE"):
            self.i += 1
        chunk = self.s[start:self.i]
        try:
            return round(float(chunk), 6)
        except ValueError:
            return None


def parse_newick_tree(path) -> dict | None:
    """Parse a Newick tree *file* into a nested dict, or return None on failure.

    Returns None (falsy) when the file is missing, empty or unparseable, so
    callers can guard with ``if tree:``.
    """
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    try:
        text = p.read_text(encoding="utf-8").strip()
    except Exception as exc:
        print(f"  WARNING: cannot read tree {p}: {exc}", file=sys.stderr)
        return None
    if not text:
        return None
    # Strip NHX / comment blocks and keep only the first tree statement.
    text = re.sub(r"\[[^\]]*\]", "", text)
    if ";" in text:
        text = text[: text.index(";") + 1]
    try:
        tree = _NewickReader(text).read_node()
    except Exception as exc:
        print(f"  WARNING: cannot parse tree {p}: {exc}", file=sys.stderr)
        return None
    return tree or None


def _level_from_tree_path(path) -> str:
    """Infer the clade level name from a pruned-tree filename.

    Filename convention: ``{level}.pruned.tree`` (falls back to the stem).
    """
    name = Path(path).name
    for suffix in (".pruned.tree", ".tree", ".newick", ".nwk"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


# ── Data loading ───────────────────────────────────────────────────────────────

def load_tsvs(paths: list[str]) -> pd.DataFrame:
    frames = []
    for p in paths:
        path = Path(p)
        if not path.exists() or path.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(path, sep="\t")
            if not df.empty:
                frames.append(df)
        except Exception as exc:
            print(f"  WARNING: Cannot read {path}: {exc}", file=sys.stderr)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# ── JSON data builder ──────────────────────────────────────────────────────────

def build_data(stats: pd.DataFrame,
               links: pd.DataFrame,
               levels: list[str],
               trees: dict | None = None,
               gene_trees: dict | None = None) -> dict:
    """Build the per-HG data structure embedded in the HTML.

    Parameters
    ----------
    stats, links : pandas.DataFrame
        Per-OG stats and cross-level link rows (from link_hog_levels.py or
        build_hog_report.py).
    levels : list[str]
        Ordered clade level names (broad → narrow).
    trees : dict | None
        Mapping ``{level: parsed_tree_dict}`` of pruned species trees, where
        each value is the nested dict returned by ``parse_newick_tree``.
    gene_trees : dict | None
        Mapping ``{hg: newick_string}`` of POSSVM-annotated gene trees.
    """
    trees = trees or {}
    gene_trees = gene_trees or {}
    hgs: dict[str, dict] = {}

    if stats.empty:
        return {"levels": levels, "hgs": hgs,
                "trees": trees, "gene_trees": gene_trees}

    for hg, hg_stats in stats.groupby("hg", sort=True):
        hg = str(hg)
        parts = hg.split(".")
        prefix = parts[0] if parts else ""
        family = parts[1] if len(parts) > 1 else ""
        hg_num = parts[2] if len(parts) > 2 else ""

        nodes = []
        og_counts: dict[str, int] = {}
        for _, row in hg_stats.iterrows():
            lv = str(row["level"])
            nodes.append({
                "id":        str(row["og"]),
                "level":     lv,
                "n_genes":   int(row["n_genes"]),
                "n_species": int(row["n_species"]),
            })
            og_counts[lv] = og_counts.get(lv, 0) + 1

        hg_links: list[dict] = []
        if not links.empty and "hg" in links.columns:
            subset = links[links["hg"].astype(str) == hg]
            for _, row in subset.iterrows():
                hg_links.append({
                    "source":  str(row["parent_og"]),
                    "target":  str(row["child_og"]),
                    "n_genes": int(row["n_genes"]),
                })

        # Summary string: "3 → 5 → 8" counts per level in declared order
        counts_str = " → ".join(
            str(og_counts.get(lv, 0)) for lv in levels if lv in og_counts
        )

        hgs[hg] = {
            "prefix":     prefix,
            "family":     family,
            "hg_num":     hg_num,
            "og_counts":  og_counts,
            "counts_str": counts_str,
            "nodes":      nodes,
            "links":      hg_links,
            "has_gene_tree": hg in gene_trees,
        }

    return {"levels": levels, "hgs": hgs,
            "trees": trees, "gene_trees": gene_trees}


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate interactive hOG hierarchy Sankey HTML."
    )
    parser.add_argument("--links",  nargs="+", required=True,
                        help="*.og_links.tsv files from link_hog_levels.py")
    parser.add_argument("--stats",  nargs="+", required=True,
                        help="*.og_stats.tsv files from link_hog_levels.py")
    parser.add_argument("--trees",  nargs="+", required=True,
                        help="Pruned species tree files ({level}.pruned.tree)")
    parser.add_argument("--levels", required=True,
                        help="Ordered clade names, comma-separated (broad→narrow)")
    parser.add_argument("--output", required=True, help="Output HTML file")
    args = parser.parse_args()

    levels = [lv.strip() for lv in args.levels.split(",") if lv.strip()]

    stats = load_tsvs(args.stats)
    links = load_tsvs(args.links)

    # ── Parse pruned species trees, keyed by inferred clade level ──────────────
    trees: dict[str, dict] = {}
    for tp in (args.trees or []):
        level = _level_from_tree_path(tp)
        tree = parse_newick_tree(tp)
        if tree:
            trees[level] = tree
        else:
            print(f"  WARNING: could not parse tree '{tp}'", file=sys.stderr)

    n_hgs = stats["hg"].nunique() if not stats.empty else 0
    print(f"Loaded {n_hgs} HGs across {len(levels)} levels, "
          f"{len(trees)} species trees", file=sys.stderr)

    data = build_data(stats, links, levels, trees)
    data_json = json.dumps(data, separators=(",", ":"))

    html = HTML_TEMPLATE.replace("%%DATA_JSON%%", data_json)
    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Written: {args.output}", file=sys.stderr)


# ── HTML template ──────────────────────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>hOG Hierarchy</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{display:flex;height:100vh;font-family:"Helvetica Neue",Arial,sans-serif;
     font-size:13px;background:#f5f5f5;color:#333;overflow:hidden}

/* ── Sidebar ── */
#sidebar{width:290px;min-width:260px;display:flex;flex-direction:column;
  background:#fff;border-right:1px solid #ddd;z-index:10}
#sidebar-header{padding:14px 14px 10px;background:#1e3a5f;color:#fff}
#sidebar-header h2{font-size:15px;font-weight:600}
#sidebar-header .sub{font-size:11px;opacity:.75;margin-top:3px}
#controls{padding:8px 10px;border-bottom:1px solid #eee}
#search{width:100%;padding:6px 8px;border:1px solid #ccc;border-radius:4px;font-size:12px}
#hg-count{font-size:11px;color:#888;margin-top:5px}
#hg-list{flex:1;overflow-y:auto;padding:4px 0}

.hg-group{}
.hg-group-header{display:flex;align-items:center;padding:5px 10px;cursor:pointer;
  font-weight:600;font-size:12px;background:#f4f4f4;border-bottom:1px solid #eee;
  user-select:none}
.hg-group-header:hover{background:#ececec}
.hg-group-header .arrow{margin-right:6px;font-size:10px;transition:transform .15s}
.hg-group-header.collapsed .arrow{transform:rotate(-90deg)}
.hg-group-body{}

.fam-header{padding:3px 20px;font-size:11px;color:#555;cursor:pointer;
  background:#fafafa;border-bottom:1px solid #f0f0f0;user-select:none}
.fam-header:hover{background:#f0f0f0}
.fam-body{}

.hg-item{display:flex;align-items:center;padding:4px 28px;cursor:pointer;
  border-bottom:1px solid #f5f5f5;gap:6px}
.hg-item:hover{background:#f0f0f0}
.hg-item.selected{background:#dbeaff}
.hg-label{flex:1;font-size:11px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.hg-counts{font-size:10px;color:#888;white-space:nowrap}

/* ── Main ── */
#main{flex:1;display:flex;flex-direction:column;overflow:hidden}
#top-bar{display:flex;align-items:center;padding:8px 16px;
  background:#fff;border-bottom:1px solid #ddd;gap:12px;min-height:40px}
#hg-title{font-size:13px;font-weight:600;color:#1e3a5f}
#hg-subtitle{font-size:11px;color:#888;flex:1}
#tree-btn{font-size:11px;padding:4px 10px;border:1px solid #ccc;border-radius:4px;
  background:#f7f7f7;cursor:pointer;color:#444}
#tree-btn:hover{background:#eee}
#tree-btn.active{background:#dbeaff;border-color:#9cc2f0;color:#1e3a5f}
#canvas-wrap{flex:1;overflow:hidden;position:relative;background:#fff}
#sankey-svg{width:100%;height:100%}

/* ── Tree overlay ── */
#tree-overlay{position:absolute;top:0;right:0;bottom:0;width:340px;display:none;
  background:rgba(255,255,255,.98);border-left:1px solid #ddd;overflow-y:auto;
  padding:10px 12px;box-shadow:-2px 0 8px rgba(0,0,0,.06);z-index:20}
#tree-overlay h3{font-size:12px;color:#1e3a5f;margin:10px 0 4px;font-weight:600}
#tree-overlay h3:first-child{margin-top:0}
#tree-overlay .tree-svg{width:100%;display:block}
#tree-overlay pre{font-size:10px;white-space:pre-wrap;word-break:break-all;
  background:#f7f7f7;border:1px solid #eee;border-radius:4px;padding:6px;color:#555}
#tree-overlay .none{font-size:11px;color:#aaa;font-style:italic}

/* ── Tooltip ── */
#tooltip{position:fixed;display:none;pointer-events:none;
  background:rgba(255,255,255,.97);border:1px solid #bbb;
  border-radius:5px;padding:8px 10px;font-size:12px;
  box-shadow:0 2px 8px rgba(0,0,0,.15);z-index:100;max-width:260px}
#tooltip .tt-name{font-weight:700;margin-bottom:4px;font-size:13px}
#tooltip .tt-row{display:flex;justify-content:space-between;gap:12px;color:#555}

/* ── Empty state ── */
#empty-msg{position:absolute;inset:0;display:flex;align-items:center;
  justify-content:center;color:#aaa;font-size:14px;pointer-events:none}
</style>
</head>
<body>

<div id="sidebar">
  <div id="sidebar-header">
    <h2>hOG Hierarchy</h2>
    <div class="sub">Hierarchical orthogroup expansion across clades</div>
  </div>
  <div id="controls">
    <input id="search" type="text" placeholder="Search gene families…" oninput="filterHGs()">
    <div id="hg-count"></div>
  </div>
  <div id="hg-list"></div>
</div>

<div id="main">
  <div id="top-bar">
    <span id="hg-title">← Select a gene family</span>
    <span id="hg-subtitle"></span>
    <button id="tree-btn" onclick="toggleTrees()">Trees</button>
  </div>
  <div id="canvas-wrap">
    <svg id="sankey-svg"></svg>
    <div id="empty-msg">No data for this HG</div>
    <div id="tree-overlay"></div>
  </div>
</div>

<div id="tooltip">
  <div class="tt-name" id="tt-name"></div>
  <div id="tt-body"></div>
</div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
// ── Embedded data ─────────────────────────────────────────────────────────────
const APP = %%DATA_JSON%%;
const LEVELS     = APP.levels;
const HGS        = APP.hgs;
const TREES      = APP.trees || {};
const GENE_TREES = APP.gene_trees || {};

// ── Colour palette for top-level OG lineages ──────────────────────────────────
const PALETTE = [
  "#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f",
  "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac",
  "#86bcb6","#f1ce63","#a0cbe8","#ffbe7d","#8cd17d",
];

// ── OG short display name ─────────────────────────────────────────────────────
function ogShort(ogId, hg, level) {
  const prefix = level + "." + hg + ".";
  const rest = ogId.startsWith(prefix) ? ogId.slice(prefix.length) : ogId;
  const colon = rest.indexOf(":");
  const num = colon >= 0 ? rest.slice(0, colon) : rest;
  const ref = colon >= 0 ? rest.slice(colon + 1) : "";
  return ref ? ref : ("OG" + num);
}

// ── State ─────────────────────────────────────────────────────────────────────
let selectedHG = null;
let treesOpen  = false;

// ── Sidebar ───────────────────────────────────────────────────────────────────
function buildSidebar() {
  // Group: prefix → family → [hg_id]
  const groups = {};
  for (const [hgId, meta] of Object.entries(HGS)) {
    const p = meta.prefix || "other";
    const f = meta.family || "—";
    if (!groups[p]) groups[p] = {};
    if (!groups[p][f]) groups[p][f] = [];
    groups[p][f].push(hgId);
  }

  const container = document.getElementById("hg-list");
  container.innerHTML = "";

  for (const [prefix, fams] of Object.entries(groups).sort()) {
    const groupDiv = document.createElement("div");
    groupDiv.className = "hg-group";
    groupDiv.dataset.prefix = prefix;

    const totalHGs = Object.values(fams).flat().length;
    const hdr = document.createElement("div");
    hdr.className = "hg-group-header";
    hdr.innerHTML = `<span class="arrow">▾</span>${prefix}<span style="font-weight:400;color:#999;margin-left:6px">(${totalHGs})</span>`;
    hdr.onclick = () => {
      const body = groupDiv.querySelector(".hg-group-body");
      const c = hdr.classList.toggle("collapsed");
      body.style.display = c ? "none" : "";
    };
    groupDiv.appendChild(hdr);

    const groupBody = document.createElement("div");
    groupBody.className = "hg-group-body";

    for (const [fam, hgIds] of Object.entries(fams).sort()) {
      const famDiv = document.createElement("div");
      const famHdr = document.createElement("div");
      famHdr.className = "fam-header";
      famHdr.textContent = `${fam} (${hgIds.length})`;
      famHdr.onclick = () => {
        const b = famDiv.querySelector(".fam-body");
        b.style.display = b.style.display === "none" ? "" : "none";
      };
      famDiv.appendChild(famHdr);

      const famBody = document.createElement("div");
      famBody.className = "fam-body";
      for (const hgId of hgIds.sort()) {
        const meta = HGS[hgId];
        const item = document.createElement("div");
        item.className = "hg-item";
        item.dataset.hg = hgId;
        item.innerHTML = `<span class="hg-label" title="${hgId}">${meta.hg_num || hgId}</span>
                          <span class="hg-counts">${meta.counts_str}</span>`;
        item.onclick = () => selectHG(hgId);
        famBody.appendChild(item);
      }
      famDiv.appendChild(famBody);
      groupBody.appendChild(famDiv);
    }
    groupDiv.appendChild(groupBody);
    container.appendChild(groupDiv);
  }
  updateHGCount();
}

function updateHGCount() {
  const n = document.querySelectorAll(".hg-item:not([style*='display: none'])").length;
  document.getElementById("hg-count").textContent = `${n} gene families`;
}

function filterHGs() {
  const q = document.getElementById("search").value.toLowerCase().trim();
  document.querySelectorAll(".hg-item").forEach(el => {
    const match = !q || el.dataset.hg.toLowerCase().includes(q);
    el.style.display = match ? "" : "none";
  });
  // Show/hide family headers
  document.querySelectorAll(".fam-body").forEach(body => {
    const anyVisible = [...body.querySelectorAll(".hg-item")].some(el => el.style.display !== "none");
    body.closest("div").style.display = anyVisible ? "" : "none";
  });
  document.querySelectorAll(".hg-group").forEach(g => {
    const anyVisible = [...g.querySelectorAll(".hg-item")].some(el => el.style.display !== "none");
    g.style.display = anyVisible ? "" : "none";
  });
  updateHGCount();
}

// ── HG selection ──────────────────────────────────────────────────────────────
function selectHG(hgId) {
  document.querySelectorAll(".hg-item").forEach(el => {
    el.classList.toggle("selected", el.dataset.hg === hgId);
  });
  selectedHG = hgId;
  const meta = HGS[hgId];
  document.getElementById("hg-title").textContent = hgId;
  const bits = [];
  if (meta.counts_str) bits.push("OGs per level: " + meta.counts_str);
  if (meta.has_gene_tree || GENE_TREES[hgId]) bits.push("gene tree available");
  document.getElementById("hg-subtitle").textContent = bits.join("  ·  ");
  renderSankey(hgId);
  if (treesOpen) renderTrees(hgId);
}

// ── Sankey rendering ──────────────────────────────────────────────────────────
function renderSankey(hgId) {
  const svg = document.getElementById("sankey-svg");
  svg.innerHTML = "";
  document.getElementById("empty-msg").style.display = "none";

  const meta = HGS[hgId];
  if (!meta || !meta.nodes || !meta.nodes.length) {
    document.getElementById("empty-msg").style.display = "";
    return;
  }

  const W = svg.clientWidth  || svg.parentElement.clientWidth;
  const H = svg.clientHeight || svg.parentElement.clientHeight;

  // Active levels: only those with at least one node for this HG
  const presentLevels = new Set(meta.nodes.map(n => n.level));
  const activeLevels  = LEVELS.filter(lv => presentLevels.has(lv));
  if (!activeLevels.length) {
    document.getElementById("empty-msg").style.display = "";
    return;
  }

  // ── Layout constants ───────────────────────────────────────────────────────
  const MARGIN_X  = 20;
  const MARGIN_Y  = 50;   // leave room for column header
  const COL_W     = 16;   // rectangle width
  const NODE_PAD  = 6;    // gap between nodes in same column
  const MIN_H     = 8;    // minimum node height in pixels
  const nCols     = activeLevels.length;
  const colSpan   = nCols > 1 ? (W - 2 * MARGIN_X - COL_W) / (nCols - 1) : W - 2 * MARGIN_X;

  // Column x positions (left edge of rectangle)
  const colX = {};
  activeLevels.forEach((lv, i) => { colX[lv] = MARGIN_X + i * colSpan; });

  // Group nodes by level
  const byLevel = {};
  activeLevels.forEach(lv => { byLevel[lv] = []; });
  meta.nodes.forEach(n => { if (byLevel[n.level]) byLevel[n.level].push({...n}); });

  // Primary parent lookup: child_id → link with most n_genes
  const primaryParent = {};
  (meta.links || []).forEach(lk => {
    if (!primaryParent[lk.target] || primaryParent[lk.target].n_genes < lk.n_genes)
      primaryParent[lk.target] = lk;
  });

  // Assign colours by top-level OG lineage
  const rootLevel = activeLevels[0];
  const rootNodes = byLevel[rootLevel];
  rootNodes.sort((a, b) => a.id.localeCompare(b.id));
  const lineageColor = {};
  rootNodes.forEach((n, i) => { lineageColor[n.id] = PALETTE[i % PALETTE.length]; });

  // Propagate lineage colour to descendants via links (BFS)
  const nodeColor = {...lineageColor};
  const queue = [...rootNodes.map(n => n.id)];
  const visited = new Set(queue);
  while (queue.length) {
    const pid = queue.shift();
    (meta.links || []).filter(lk => lk.source === pid).forEach(lk => {
      if (!visited.has(lk.target)) {
        nodeColor[lk.target] = nodeColor[pid] || "#aaa";
        visited.add(lk.target);
        queue.push(lk.target);
      }
    });
  }

  // pxPerGene: scale so first-level total fits usable height
  const firstTotal = byLevel[rootLevel].reduce((s, n) => s + n.n_genes, 0) || 1;
  const usableH    = H - MARGIN_Y - 20;
  const pxPerGene  = Math.max(
    (usableH - NODE_PAD * (byLevel[rootLevel].length - 1)) / firstTotal,
    MIN_H / Math.max(...byLevel[rootLevel].map(n => n.n_genes), 1)
  );

  // Compute node positions (x, y_top, h)
  const nodePos = {};

  function assignColumnY(level, nodes) {
    // Sort: level 0 alphabetically; subsequent levels by primary parent midY
    if (level === rootLevel) {
      nodes.sort((a, b) => a.id.localeCompare(b.id));
    } else {
      nodes.sort((a, b) => {
        const pa = primaryParent[a.id], pb = primaryParent[b.id];
        const ya = pa && nodePos[pa.source]
          ? nodePos[pa.source].yt + nodePos[pa.source].h / 2 : 1e9;
        const yb = pb && nodePos[pb.source]
          ? nodePos[pb.source].yt + nodePos[pb.source].h / 2 : 1e9;
        return ya - yb || a.id.localeCompare(b.id);
      });
    }
    let y = MARGIN_Y;
    nodes.forEach(n => {
      const h = Math.max(n.n_genes * pxPerGene, MIN_H);
      nodePos[n.id] = { x: colX[n.level], yt: y, h, level: n.level };
      y += h + NODE_PAD;
    });
  }

  activeLevels.forEach(lv => assignColumnY(lv, byLevel[lv]));

  // ── Draw using D3 ──────────────────────────────────────────────────────────
  const s = d3.select(svg);

  // Column headers
  activeLevels.forEach(lv => {
    s.append("text")
      .attr("x", colX[lv] + COL_W / 2)
      .attr("y", MARGIN_Y - 10)
      .attr("text-anchor", "middle")
      .attr("font-size", 11)
      .attr("fill", "#555")
      .attr("font-weight", "600")
      .text(lv);
  });

  // Ribbons (links)
  // For each node, track current outgoing / incoming y offset
  const outCur = {};  // node_id → current y for next outgoing ribbon
  const inCur  = {};  // node_id → current y for next incoming ribbon
  for (const [id, pos] of Object.entries(nodePos)) {
    outCur[id] = pos.yt;
    inCur[id]  = pos.yt;
  }

  // Sort links: by source y then target y (reduces crossings)
  const sortedLinks = [...(meta.links || [])].filter(
    lk => nodePos[lk.source] && nodePos[lk.target]
  ).sort((a, b) => {
    const ya = nodePos[a.source].yt, yb = nodePos[b.source].yt;
    return ya !== yb ? ya - yb : nodePos[a.target].yt - nodePos[b.target].yt;
  });

  const linkGroup = s.append("g").attr("class", "links");
  sortedLinks.forEach(lk => {
    const sp = nodePos[lk.source], tp = nodePos[lk.target];
    if (!sp || !tp) return;

    const srcH = Math.max(lk.n_genes * pxPerGene, 2);
    const tgtH = Math.max(lk.n_genes * pxPerGene, 2);

    const x1 = sp.x + COL_W, x2 = tp.x;
    const y1t = outCur[lk.source], y1b = y1t + srcH;
    const y2t = inCur[lk.target],  y2b = y2t + tgtH;
    outCur[lk.source] = y1b;
    inCur[lk.target]  = y2b;

    const mx = (x1 + x2) / 2;
    const d  = `M${x1},${y1t} C${mx},${y1t} ${mx},${y2t} ${x2},${y2t}
                L${x2},${y2b} C${mx},${y2b} ${mx},${y1b} ${x1},${y1b} Z`;

    const color = nodeColor[lk.source] || "#aaa";
    linkGroup.append("path")
      .attr("d", d)
      .attr("fill", color)
      .attr("opacity", 0.35)
      .attr("stroke", color)
      .attr("stroke-width", 0.3);
  });

  // Nodes (rectangles)
  const nodeGroup = s.append("g").attr("class", "nodes");
  meta.nodes.forEach(n => {
    const pos = nodePos[n.id];
    if (!pos) return;
    const color  = nodeColor[n.id] || "#aaa";
    const hg     = hgId;
    const lv     = n.level;
    const label  = ogShort(n.id, hg, lv);
    const numChildren = (meta.links || []).filter(lk => lk.source === n.id).length;
    const numParents  = (meta.links || []).filter(lk => lk.target === n.id).length;

    const g = nodeGroup.append("g")
      .attr("cursor", "pointer")
      .on("mousemove", (event) => {
        const parent = (meta.links || []).find(lk => lk.target === n.id);
        const parentStr = parent ? ogShort(parent.source, hg, primaryParent[n.id]?.level || "") : "—";
        showTooltip(event,
          label,
          [
            ["Level",    lv],
            ["OG",       n.id.split(".").slice(-1)[0]],
            ["Genes",    n.n_genes],
            ["Species",  n.n_species],
            numParents  ? ["From",   parentStr]   : null,
            numChildren ? ["Splits into", numChildren + " OG" + (numChildren>1?"s":"")] : null,
          ]
        );
      })
      .on("mouseleave", hideTooltip);

    // Rectangle
    g.append("rect")
      .attr("x", pos.x)
      .attr("y", pos.yt)
      .attr("width",  COL_W)
      .attr("height", pos.h)
      .attr("fill",   color)
      .attr("stroke", d3.color(color).darker(0.5))
      .attr("stroke-width", 0.8)
      .attr("rx", 2);

    // Label to the right of the rightmost column, or below for intermediate
    const isLast = lv === activeLevels[activeLevels.length - 1];
    if (isLast || activeLevels.length === 1) {
      g.append("text")
        .attr("x", pos.x + COL_W + 5)
        .attr("y", pos.yt + pos.h / 2 + 4)
        .attr("font-size", Math.min(11, pos.h * 0.85))
        .attr("fill", "#444")
        .text(label);
    } else if (pos.h >= 14) {
      // Label centred inside rectangle if tall enough
      g.append("text")
        .attr("x", pos.x + COL_W / 2)
        .attr("y", pos.yt + pos.h / 2 + 4)
        .attr("text-anchor", "middle")
        .attr("font-size", Math.min(10, pos.h * 0.75))
        .attr("fill", "#fff")
        .attr("pointer-events", "none")
        .text(label.length > 8 ? label.slice(0, 7) + "…" : label);
    }
  });
}

// ── Species / gene tree overlay ────────────────────────────────────────────────
function toggleTrees() {
  treesOpen = !treesOpen;
  document.getElementById("tree-btn").classList.toggle("active", treesOpen);
  document.getElementById("tree-overlay").style.display = treesOpen ? "block" : "none";
  if (treesOpen && selectedHG) renderTrees(selectedHG);
}

// Collect leaves of a parsed-tree dict in traversal order.
function treeLeaves(node, acc) {
  if (!node) return acc;
  if (!node.children || !node.children.length) { acc.push(node); return acc; }
  node.children.forEach(c => treeLeaves(c, acc));
  return acc;
}

// Draw a compact cladogram for a parsed-tree dict; returns an SVG element.
function drawSpeciesTree(root) {
  const leaves = treeLeaves(root, []);
  if (!leaves.length) return null;
  const rowH = 15, padTop = 6, padBot = 6;
  const width = 320, labelW = 120;
  const height = padTop + padBot + leaves.length * rowH;

  // Assign leaf y positions.
  leaves.forEach((l, i) => { l._y = padTop + i * rowH + rowH / 2; });
  // Depth of each node.
  let maxDepth = 0;
  (function depth(n, d) {
    n._d = d;
    if (n.children && n.children.length) n.children.forEach(c => depth(c, d + 1));
    else maxDepth = Math.max(maxDepth, d);
  })(root, 0);
  // Internal y = midpoint of children.
  (function midY(n) {
    if (n.children && n.children.length) {
      n.children.forEach(midY);
      n._y = (n.children[0]._y + n.children[n.children.length - 1]._y) / 2;
    }
  })(root);

  const xStep = maxDepth > 0 ? (width - labelW) / maxDepth : 0;
  const xOf = n => 4 + n._d * xStep;

  const svgNS = "http://www.w3.org/2000/svg";
  const svg = document.createElementNS(svgNS, "svg");
  svg.setAttribute("class", "tree-svg");
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.setAttribute("height", height);

  function line(x1, y1, x2, y2) {
    const el = document.createElementNS(svgNS, "line");
    el.setAttribute("x1", x1); el.setAttribute("y1", y1);
    el.setAttribute("x2", x2); el.setAttribute("y2", y2);
    el.setAttribute("stroke", "#7a93ad");
    el.setAttribute("stroke-width", "1");
    svg.appendChild(el);
  }

  (function draw(n) {
    const x = xOf(n);
    if (n.children && n.children.length) {
      const ys = n.children.map(c => c._y);
      // vertical connector spanning children
      line(x, Math.min(...ys), x, Math.max(...ys));
      n.children.forEach(c => {
        line(x, c._y, xOf(c), c._y);  // horizontal to each child
        draw(c);
      });
    } else {
      const t = document.createElementNS(svgNS, "text");
      t.setAttribute("x", x + 4);
      t.setAttribute("y", n._y + 3);
      t.setAttribute("font-size", "10");
      t.setAttribute("fill", "#444");
      t.textContent = n.name || "";
      svg.appendChild(t);
    }
  })(root);

  return svg;
}

function renderTrees(hgId) {
  const box = document.getElementById("tree-overlay");
  box.innerHTML = "";
  const meta = HGS[hgId];

  // Species trees for the levels present in this HG (broad → narrow).
  const presentLevels = meta ? new Set(meta.nodes.map(n => n.level)) : new Set();
  const active = LEVELS.filter(lv => presentLevels.has(lv) && TREES[lv]);
  const h = document.createElement("h3");
  h.textContent = "Pruned species trees";
  box.appendChild(h);
  if (!active.length) {
    const p = document.createElement("div");
    p.className = "none";
    p.textContent = "No species trees available.";
    box.appendChild(p);
  } else {
    active.forEach(lv => {
      const sub = document.createElement("h3");
      sub.textContent = lv;
      sub.style.fontSize = "11px";
      sub.style.color = "#777";
      box.appendChild(sub);
      const svg = drawSpeciesTree(TREES[lv]);
      if (svg) box.appendChild(svg);
    });
  }

  // Gene tree (Newick) for this HG, if available.
  const gh = document.createElement("h3");
  gh.textContent = "Gene tree (Newick)";
  box.appendChild(gh);
  const nwk = GENE_TREES[hgId];
  if (nwk) {
    const pre = document.createElement("pre");
    pre.textContent = nwk;
    box.appendChild(pre);
  } else {
    const p = document.createElement("div");
    p.className = "none";
    p.textContent = "No gene tree for this HG.";
    box.appendChild(p);
  }
}

// ── Tooltip ───────────────────────────────────────────────────────────────────
function showTooltip(event, name, rows) {
  const tt = document.getElementById("tooltip");
  document.getElementById("tt-name").textContent = name;
  const body = document.getElementById("tt-body");
  body.innerHTML = "";
  rows.forEach(r => {
    if (!r) return;
    const div = document.createElement("div");
    div.className = "tt-row";
    div.innerHTML = `<span>${r[0]}</span><span><b>${r[1]}</b></span>`;
    body.appendChild(div);
  });
  tt.style.display = "block";
  tt.style.left = (event.clientX + 14) + "px";
  tt.style.top  = (event.clientY - 10) + "px";
}
function hideTooltip() {
  document.getElementById("tooltip").style.display = "none";
}

// ── Init ──────────────────────────────────────────────────────────────────────
buildSidebar();

// Redraw on resize
let resizeTimer;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => { if (selectedHG) renderSankey(selectedHG); }, 120);
});
</script>
</body>
</html>"""


if __name__ == "__main__":
    main()
