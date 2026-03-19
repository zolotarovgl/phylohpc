#!/usr/bin/env python3
"""Generate a self-contained interactive HTML report for step2 (POSSVM) outputs.

Combines:
  • Heatmap start page (species × families / HGs) with species-tree cladogram
  • Drill-down:  family → HG columns → interactive gene tree
  • Gene tree explorer with OG labels, collapse/expand, clade colouring, highlight

Reads
-----
  results/possvm/*.ortholog_groups.newick   gene trees with OG labels
  results/possvm/*.ortholog_groups.csv      gene → OG mapping (tab-sep)
  results/search/*.genes.list               gene lists per family
  results/clusters/*.fasta                  per-HG FASTA files
  data/gene_families_searchinfo.csv         family → class mapping
  data/species_tree.full.newick             species tree (named internal nodes)

Outputs
-------
  A self-contained HTML file.
"""

import argparse
import html as _html
import json
import sys
from collections import defaultdict
from pathlib import Path


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_species_prefix(gene_id: str) -> str:
    """Return species prefix from a gene ID such as 'Mmus_ENSG123' → 'Mmus'.

    Also strips pipe-separated annotations (e.g. 'Mmus_ENSG123 | OG001 | ref')
    before extracting the prefix.
    """
    # Strip pipe-separated tip annotations produced by POSSVM
    gene_id = gene_id.split("|")[0].strip()
    for sep in ("_", "."):
        idx = gene_id.find(sep)
        if idx > 0:  # require a non-empty prefix
            return gene_id[:idx]
    return gene_id


# ── Heatmap data loaders (from report_step1) ─────────────────────────────────

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


def parse_fasta_species(path) -> dict:
    """Return {species: count} from FASTA headers."""
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


def build_family_records(search_dir: Path, family_info: dict) -> list:
    """Return list of family records from *.genes.list files."""
    records = []
    if not search_dir.is_dir():
        return records
    for genes_file in sorted(search_dir.glob("*.genes.list")):
        name = genes_file.name
        if not name.endswith(".genes.list"):
            continue
        stem = name[: -len(".genes.list")]
        parts = stem.split(".", 1)
        if len(parts) < 2:
            continue
        pref, family = parts[0], parts[1]
        genes = []
        try:
            with open(genes_file) as fh:
                for line in fh:
                    g = line.strip()
                    if g and not g.startswith("#"):
                        genes.append(g)
        except OSError:
            continue
        if not genes:
            continue
        sp_counts: dict = defaultdict(int)
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
    """Return list of HG records from *.fasta files."""
    records = []
    if not cluster_dir.is_dir():
        return records
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


def load_tree_data(tree_path: str) -> tuple:
    """Return (species_order, tree_dict) from a newick species tree."""
    try:
        from ete3 import Tree  # type: ignore
        t = Tree(tree_path, format=1)
        species_order = [n.name for n in t.get_leaves()]
        return species_order, _sp_tree_to_dict(t)
    except Exception:
        return [], {}


def _sp_tree_to_dict(node) -> dict:
    d: dict = {"name": node.name or "", "dist": round(float(node.dist), 6)}
    if node.is_leaf():
        d["leaf"] = True
    else:
        d["children"] = [_sp_tree_to_dict(c) for c in node.children]
    return d


# ── Gene-tree / POSSVM loaders ───────────────────────────────────────────────

def gene_tree_to_dict(node) -> dict:
    """Recursively convert an ete3 node to a JSON-serialisable dict.

    Supports pipe-separated POSSVM tip annotations of the form:
        gene_id | orthogroup_name | reference_ortholog
    The full name is kept for display; gene_id and og are stored separately so
    that JavaScript can perform OG-based collapsing and tooltip lookups.
    """
    name = node.name or ""
    d: dict = {"name": name, "dist": round(float(node.dist), 6)}
    if node.is_leaf():
        d["leaf"] = True
        parts = [p.strip() for p in name.split("|")]
        gene_id = parts[0]
        d["gene_id"] = gene_id
        d["species"] = get_species_prefix(gene_id) if gene_id else ""
        if len(parts) >= 2 and parts[1]:
            d["og"] = parts[1]
        if len(parts) >= 3 and parts[2] and parts[2].upper() != "NA":
            d["ref"] = parts[2]
    else:
        d["children"] = [gene_tree_to_dict(c) for c in node.children]
    return d


def load_og_csv(csv_path: Path) -> dict:
    """Return {og_name: [gene_id, ...]} from a POSSVM CSV."""
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
    """Return (tree_records, all_species)."""
    try:
        from ete3 import Tree  # type: ignore
    except ImportError:
        print("ERROR: ete3 not installed. Run: pip install ete3", file=sys.stderr)
        return [], []

    records = []
    all_species: set = set()

    for nwk in sorted(possvm_dir.glob("*.ortholog_groups.newick")):
        stem = nwk.stem
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
            "og_names":         sorted(ogs.keys()),
            "n_ogs":            len(ogs),
            "tree_dict":        tree_dict,
            "ogs":              ogs,
        })

    return records, sorted(all_species)


# ── Clade groupings from species tree ────────────────────────────────────────

def _direct_named_children(node) -> list:
    """Named internal descendants with no other named internal node between."""
    result = []
    for child in node.children:
        if not child.is_leaf() and child.name:
            result.append(child)
        elif not child.is_leaf():
            result.extend(_direct_named_children(child))
    return result


def parse_clade_groupings(species_tree_path: Path) -> list:
    """Return [{name, groups: {species: clade_name}}, ...] for colouring."""
    try:
        from ete3 import Tree  # type: ignore
    except ImportError:
        return []
    try:
        t = Tree(str(species_tree_path), format=1)
    except Exception as exc:
        print(f"WARN: species tree parse error: {exc}", file=sys.stderr)
        return []

    groupings = []
    for node in t.traverse("levelorder"):
        if node.is_leaf() or not node.name:
            continue
        named_children = _direct_named_children(node)
        if len(named_children) < 2:
            continue
        species_map: dict = {}
        for nc in named_children:
            for leaf in nc.get_leaves():
                species_map[leaf.name] = nc.name
        if len(set(species_map.values())) < 2 or len(species_map) < 4:
            continue
        groupings.append({
            "name": node.name,
            "groups": species_map,
            "_n": len(species_map),
        })

    for g in groupings:
        del g["_n"]
    return groupings  # levelorder already gives most-basal (root-proximal) splits first


# ── HTML template ─────────────────────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Step 2 Report</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
*{box-sizing:border-box;margin:0;padding:0}
html{height:100%;height:-webkit-fill-available}
body{height:100%;height:-webkit-fill-available;overflow:hidden;font-family:"Helvetica Neue",Helvetica,Arial,sans-serif;font-size:12px;background:#f7f7f7;color:#333;display:flex;flex-direction:column}

/* ── shared header bar ── */
.top-bar{padding:8px 14px;background:#2c3e50;color:#ecf0f1;display:flex;align-items:center;gap:10px;flex-wrap:wrap;min-height:42px}
.top-bar h2{font-size:14px;font-weight:600;white-space:nowrap}
.top-bar .btn{padding:3px 10px;border:1px solid #7f8c8d;border-radius:3px;background:transparent;color:#ecf0f1;cursor:pointer;font-size:11px}
.top-bar .btn:hover{background:#34495e}
.top-bar select{font-size:11px;padding:2px 5px;border:1px solid #7f8c8d;border-radius:3px;background:#34495e;color:#ecf0f1}
.top-bar label{font-size:11px;color:#bdc3c7}

/* ── Tab strip ── */
#body-wrap{display:flex;flex:1;overflow:hidden;min-height:0}
#tab-strip{width:36px;display:flex;flex-direction:column;background:#2c3e50;gap:2px;padding-top:4px;flex-shrink:0}
.tab-btn{writing-mode:vertical-rl;transform:rotate(180deg);padding:10px 6px;cursor:pointer;background:transparent;border:none;border-left:3px solid transparent;color:#95a5a6;font-size:11px;font-weight:600;white-space:nowrap;text-align:center}
.tab-btn.active{color:#ecf0f1;border-left-color:#1abc9c;background:#34495e}
.tab-btn:hover:not(.active){color:#ecf0f1;background:#3d5166}
.tab-pane{display:none;flex:1;overflow:hidden;flex-direction:column;min-height:0}
.tab-pane.active{display:flex}

/* ── Species tree pane ── */
#sp-tree-wrap{flex:1;overflow:auto;background:#fff;padding:16px}
#sp-tree-wrap svg{display:block}

/* ── Heatmap pane ── */
#pane-heatmap{flex-direction:row}
#hm-layout{display:flex;flex:1;overflow:hidden}
#tree-panel{width:200px;overflow-y:auto;border-right:1px solid #ccc;background:#fff;flex-shrink:0}
#heatmap-panel{flex:1;overflow:auto;background:#fff}
.hm-col{font-size:9px;cursor:pointer}
.hm-col:hover{font-weight:700}
.hm-cell:hover{stroke:#000;stroke-width:1px}

/* ── Tree pane ── */
#pane-trees{flex-direction:column}
#app{display:flex;flex:1;overflow:hidden}

/* sidebar */
#sidebar{flex:0 0 220px;display:flex;flex-direction:column;border-right:1px solid #ccc;background:#fff}
#sidebar-top{padding:8px;border-bottom:1px solid #eee}
#hg-search{width:100%;padding:5px 7px;font-size:11px;border:1px solid #ccc;border-radius:3px}
#hg-count{font-size:10px;color:#999;margin-top:4px}
#hg-list{flex:1;overflow-y:auto}
.fam-header{padding:5px 8px;background:#ecf0f1;border-bottom:1px solid #ddd;cursor:pointer;display:flex;align-items:center;gap:5px;font-size:11px;font-weight:600;color:#444;user-select:none}
.fam-header:hover{background:#dce4ec}
.fam-arrow{font-size:9px;transition:transform .15s}
.fam-header.open .fam-arrow{transform:rotate(90deg)}
.fam-body{display:none}
.fam-body.open{display:block}
.hg-item{padding:5px 10px 5px 16px;cursor:pointer;border-bottom:1px solid #f0f0f0;line-height:1.35}
.hg-item:hover{background:#f5f5f5}
.hg-item.selected{background:#d5f5e3;border-left:3px solid #1abc9c;padding-left:13px}
.hg-item .hg-name{font-weight:600;font-size:11px}
.hg-item .hg-meta{font-size:10px;color:#888}

/* main tree panel */
#main{flex:1;display:flex;flex-direction:column;overflow:hidden}
#controls{padding:5px 10px;background:#f5f5f5;border-bottom:1px solid #ddd;display:flex;align-items:center;gap:8px;flex-wrap:wrap}
.ctrl-btn{padding:3px 9px;border:1px solid #bbb;border-radius:3px;cursor:pointer;background:#fff;font-size:11px}
.ctrl-btn:hover{background:#eee}
#tree-title{font-size:11px;color:#555;margin-left:4px}
#n-ogs-label{font-size:10px;color:#888}
#hl-search{font-size:11px;padding:3px 6px;border:1px solid #bbb;border-radius:3px;width:180px}
#hl-clear{padding:2px 6px;border:1px solid #bbb;border-radius:3px;cursor:pointer;background:#fff;font-size:11px}
#hl-clear:hover{background:#eee}
#hl-tags{display:flex;flex-wrap:wrap;gap:3px;align-items:center}
.hl-tag{display:inline-flex;align-items:center;gap:3px;padding:2px 8px;border-radius:10px;font-size:10px;color:#fff;white-space:nowrap;cursor:default}
.hl-tag-x{cursor:pointer;opacity:.7;font-size:12px;line-height:1}
.hl-tag-x:hover{opacity:1}

/* tree svg */
#tree-wrap{flex:1;overflow:hidden;position:relative;background:#fff}
#tree-svg{width:100%;height:100%;cursor:grab;display:block}
#tree-svg:active{cursor:grabbing}
.link{fill:none;stroke:#d5d5d5;stroke-width:1.3px}
.node-g circle{cursor:pointer;transition:r .12s,fill .12s}
.node-g circle:hover{stroke-width:2.5px !important}
.badge-bg{fill:#fff7ed;stroke:#e8913a;stroke-width:1.3px;cursor:pointer}
.badge-bg:hover{fill:#ffecd4}
.leaf-label{pointer-events:none;font-family:monospace}
.og-label{fill:#b5371f;pointer-events:none}
.ctrl-btn.active-btn{background:#d5f5e3;border-color:#1abc9c;color:#1a6b4a}
.scale-bar-g line{stroke:#999;stroke-width:1.5px}
.scale-bar-g text{font-size:9px;fill:#888;text-anchor:middle}

/* tooltip */
#tooltip{position:fixed;display:none;pointer-events:none;background:rgba(255,255,255,.97);border:1px solid #bbb;border-radius:5px;padding:8px 10px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.15);z-index:100;max-width:260px}
/* collapsed-node popup */
#collapsed-popup{position:fixed;display:none;background:#fff;border:1px solid #bbb;border-radius:6px;padding:10px 12px;font-size:11px;box-shadow:0 3px 10px rgba(0,0,0,.2);z-index:200;min-width:230px}
#collapsed-popup .cp-title{font-weight:700;margin-bottom:7px;font-size:12px;color:#333}
#collapsed-popup .cp-row{display:flex;align-items:center;gap:6px;margin-bottom:5px}
#collapsed-popup .cp-row input{flex:1;font-size:11px;padding:2px 5px;border:1px solid #ccc;border-radius:3px}
#collapsed-popup .cp-genes-label{font-size:10px;color:#888;margin-bottom:2px}
#collapsed-popup textarea{width:100%;height:80px;font-size:9px;font-family:monospace;resize:vertical;border:1px solid #ddd;border-radius:3px;padding:3px;box-sizing:border-box}
#collapsed-popup .cp-actions{display:flex;gap:5px;margin-top:7px;flex-wrap:wrap}
#collapsed-popup .cp-btn{font-size:10px;padding:3px 8px;border:1px solid #bbb;border-radius:3px;background:#f8f8f8;cursor:pointer}
#collapsed-popup .cp-btn:hover{background:#e8e8e8}
.tt-name{font-weight:700;margin-bottom:4px;font-size:12px}
.tt-row{display:flex;justify-content:space-between;gap:10px;color:#555;margin-top:2px}
</style>
</head>
<body>

<!-- ════════ Single top bar ════════ -->
<div class="top-bar">
  <h2>Step 2 Report</h2>
  <button class="btn" id="hm-back" style="display:none" onclick="hmBack()">&#8592; Back</button>
  <span id="hm-breadcrumb" style="font-size:11px"></span>
  <span id="tree-count" style="font-size:11px;color:#95a5a6;display:none"></span>
  <label style="margin-left:auto">Prefix:
    <select id="prefixSelect"></select>
  </label>
</div>

<!-- ════════ Body: tab strip + panes ════════ -->
<div id="body-wrap">

  <!-- vertical tab strip -->
  <div id="tab-strip">
    <button class="tab-btn active" data-tab="heatmap" onclick="switchTab('heatmap')">Heatmap</button>
    <button class="tab-btn" data-tab="trees" onclick="switchTab('trees')">Gene Trees</button>
    <button class="tab-btn" data-tab="sptree" onclick="switchTab('sptree')">Species Tree</button>
  </div>

  <!-- ── Heatmap pane ── -->
  <div class="tab-pane active" id="pane-heatmap">
    <div id="hm-layout">
      <div id="tree-panel"></div>
      <div id="heatmap-panel"></div>
    </div>
  </div>

  <!-- ── Gene tree pane ── -->
  <div class="tab-pane" id="pane-trees">
    <div id="app">
      <div id="sidebar">
        <div id="sidebar-top">
          <select id="class-filter" style="width:100%;font-size:11px;padding:4px 6px;border:1px solid #ccc;border-radius:3px;margin-bottom:4px"><option value="">All classes</option></select>
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
          <button class="ctrl-btn" id="btn-lengths" onclick="toggleLengths()">Branch lengths</button>
          <button class="ctrl-btn" onclick="fitTree()">&#x2922; Fit</button>
          <button class="ctrl-btn" id="tree-toggle" style="display:none;background:#e8f0fe;border-color:#4a90d9" onclick="toggleTreeSource()">Showing: GeneRax</button>
          <span id="tree-title"></span>
          <span id="n-ogs-label"></span>
          <span style="flex:1"></span>
          <label style="font-size:11px;color:#555">Color:
            <select id="color-by" style="font-size:11px;padding:2px 4px;border:1px solid #bbb;border-radius:3px">
              <option value="species">by species</option>
            </select>
          </label>
          <div id="hl-tags"></div>
          <input id="hl-search" list="hl-list" placeholder="Highlight… (Enter to add)">
          <datalist id="hl-list"></datalist>
          <button id="hl-clear" onclick="clearHighlight()" title="Clear all highlights">&#10005;</button>
          <button class="ctrl-btn" id="btn-focus-hl" onclick="focusHighlighted()" style="display:none" title="Collapse all branches not leading to highlighted tips">Focus</button>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Label size:
            <input type="range" id="tip-font-slider" min="6" max="24" step="1" value="11" style="width:70px;cursor:pointer;accent-color:#4a90d9">
            <span id="tip-font-val" style="width:20px;text-align:right">11</span>px
          </label>
          <span style="font-size:11px;color:#555;display:flex;align-items:center;gap:8px;margin-left:4px">
            Show:
            <label style="display:flex;align-items:center;gap:2px;cursor:pointer"><input type="checkbox" id="chk-geneid" checked> gene&nbsp;ID</label>
            <label style="display:flex;align-items:center;gap:2px;cursor:pointer"><input type="checkbox" id="chk-og" checked> OG</label>
            <label style="display:flex;align-items:center;gap:2px;cursor:pointer"><input type="checkbox" id="chk-ref" checked> ref&nbsp;ortholog</label>
          </span>
        </div>
        <div id="tree-wrap">
          <svg id="tree-svg"></svg>
        </div>
      </div>
    </div>
  </div>

  <!-- ── Species tree pane ── -->
  <div class="tab-pane" id="pane-sptree">
    <div id="sp-tree-wrap"></div>
  </div>

</div>

<div id="tooltip"></div>
<div id="collapsed-popup">
  <div class="cp-title">Collapsed node</div>
  <div class="cp-row"><span>Name:</span><input id="cp-name" type="text"></div>
  <div class="cp-genes-label" id="cp-genes-label"></div>
  <textarea id="cp-genes" readonly></textarea>
  <div class="cp-actions">
    <button class="cp-btn" onclick="cpRename()">Rename</button>
    <button class="cp-btn" onclick="cpCopy()">Copy genes</button>
    <button class="cp-btn" onclick="cpExpand()">Expand</button>
    <button class="cp-btn" onclick="cpClose()">Close</button>
  </div>
</div>

<!-- ── Per-HG lazy data ── -->
<div id="lazy-data" style="display:none">
%%LAZY_SCRIPTS%%
</div>

<script>
// ═══════════════════════════════════════════════════════════════════════════════
// DATA (injected by Python)
// ═══════════════════════════════════════════════════════════════════════════════
const SPECIES_ORDER = %%SPECIES_ORDER%%;
const SP_TREE_DATA  = %%TREE_DATA%%;
const FAMILY_DATA   = %%FAMILY_DATA%%;
const HG_DATA       = %%HG_DATA%%;
const TREE_INDEX    = %%TREE_INDEX_JSON%%;
const ALL_SPECIES   = %%SPECIES_JSON%%;
const CLADE_DATA    = %%CLADE_DATA_JSON%%;

// ═══════════════════════════════════════════════════════════════════════════════
// COLOUR SYSTEM
// ═══════════════════════════════════════════════════════════════════════════════
const palette = [
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
  "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
  "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5"
];
const SP_COLORS = {};
ALL_SPECIES.forEach((sp,i) => { SP_COLORS[sp] = palette[i % palette.length]; });
function spColor(sp) { return SP_COLORS[sp] || "#aaa"; }

let colorMode     = "species";   // "species" | "clade" | "og"
let cladeSp2Color = {};
let cladeSp2Group = {};
let cladeGrpColor = {};
let ogLeaf2Color  = {};   // gene_id → color
let ogName2Color  = {};   // og_name → color
let ogGene2Name   = {};   // gene_id → og_name
let hlSet         = null;        // null = off; union Set<species> when active
let hlQueries     = [];          // committed query strings (tags)
let hlGroupIndex  = new Map();   // species → group index (for per-group color)
let tipFontSize   = null;        // null = auto; number = user override (px)
let showGeneId    = true;        // tip label parts
let showOGName    = true;
let showRefOrtho  = true;
const spCollapsed = new Set();   // node _ids collapsed in species tree
const hlTagColors = ["#e74c3c","#3498db","#27ae60","#f39c12","#8e44ad","#16a085","#e67e22","#c0392b"];

function leafColor(sp) {
  if (hlSet !== null) {
    if (!hlSet.has(sp)) return "#ccc";
    const gi = hlGroupIndex.get(sp);
    return gi !== undefined ? hlTagColors[gi % hlTagColors.length]
                            : (colorMode === "species" ? spColor(sp) : (cladeSp2Color[sp] || "#ccc"));
  }
  return colorMode === "species" ? spColor(sp) : (cladeSp2Color[sp] || "#ccc");
}
function ogLeafColor(geneId, species) {
  if (hlSet !== null) {
    if (!hlSet.has(species||"")) return "#ccc";
    const gi = hlGroupIndex.get(species||"");
    return gi !== undefined ? hlTagColors[gi % hlTagColors.length] : (ogLeaf2Color[geneId] || "#ccc");
  }
  return ogLeaf2Color[geneId] || "#ccc";
}

// ═══════════════════════════════════════════════════════════════════════════════
// NAVIGATION STATE
// ═══════════════════════════════════════════════════════════════════════════════
let hmViewMode    = "class";   // "class" | "family" | "hg" | "og"
let hmActiveClass  = null;
let hmActiveFamily = null;
let hmActiveHG     = null;   // id string from HG_DATA
let hmActiveHGRec  = null;   // full HG_DATA record (for fallback lookup)

function getSpeciesPfx(geneId){
  const g=(geneId||"").split("|")[0].trim();
  const ui=g.indexOf("_"), di=g.indexOf(".");
  const idx=[ui,di].filter(x=>x>0).reduce((a,b)=>Math.min(a,b),Infinity);
  return idx===Infinity?g:g.slice(0,idx);
}

function switchTab(name) {
  document.querySelectorAll(".tab-btn").forEach(b => b.classList.toggle("active", b.dataset.tab===name));
  document.querySelectorAll(".tab-pane").forEach(p => p.classList.toggle("active", p.id==="pane-"+name));
  const tc = document.getElementById("tree-count");
  const hb = document.getElementById("hm-back");
  const cr = document.getElementById("hm-breadcrumb");
  if (name==="trees") {
    tc.style.display = "inline";
    hb.style.display = "none"; cr.textContent = "";
    if (!currentIndex && TREE_INDEX.length) { renderSidebar(""); selectTree(TREE_INDEX[0]); }
  } else if (name==="sptree") {
    tc.style.display = "none";
    drawSpeciesTree();
  } else {
    tc.style.display = "none";
    drawCladogram(); drawHeatmap();
  }
}

// ═══════════════════════════════════════════════════════════════════════════════
// TOOLTIP (shared)
// ═══════════════════════════════════════════════════════════════════════════════
const tooltipEl = document.getElementById("tooltip");

function showTip(event, d) {
  let html;
  if (typeof d === "string") { html = d; }
  else {
    const name = d.data.name || (d.data.leaf ? "leaf" : "internal");
    html = '<div class="tt-name">' + name + '</div>';
    if (d.data.leaf) {
      const sp = d.data.species || "?";
      html += '<div class="tt-row"><span>Species</span><strong style="color:'+spColor(sp)+'">'+sp+'</strong></div>';
      if (currentDetail) {
        for (const [og,genes] of Object.entries(activeOgs())) {
          if (genes.includes(d.data.gene_id||d.data.name)) {
            html += '<div class="tt-row"><span>OG</span><strong>'+og+'</strong></div>'; break;
          }
        }
      }
    } else {
      const nL = d._children ? countDescLeaves(d._children) : d.leaves().length;
      html += '<div class="tt-row"><span>Subtree leaves</span><strong>'+nL+'</strong></div>';
      html += '<div class="tt-row" style="color:#888"><span>'+(d._children?"collapsed":"expanded")+'</span><span>click to '+(d._children?"expand":"collapse")+'</span></div>';
      if (isOGNode(d) && currentDetail && activeOgs()[d.data.name])
        html += '<div class="tt-row"><span>OG members</span><strong>'+activeOgs()[d.data.name].length+'</strong></div>';
    }
  }
  tooltipEl.innerHTML = html;
  tooltipEl.style.display = "block";
  moveTip(event);
}
function moveTip(event) {
  const x = event.clientX+14, y = event.clientY-10;
  tooltipEl.style.left = Math.min(x, window.innerWidth -tooltipEl.offsetWidth -5)+"px";
  tooltipEl.style.top  = Math.min(y, window.innerHeight-tooltipEl.offsetHeight-5)+"px";
}
function hideTip() { tooltipEl.style.display = "none"; }

// ═══════════════════════════════════════════════════════════════════════════════
// COLLAPSED NODE POPUP
// ═══════════════════════════════════════════════════════════════════════════════
const cpEl = document.getElementById("collapsed-popup");
const customNodeNames = {};
let cpActiveNode = null;

function collectLeafGenes(children) {
  const genes = [];
  (function walk(ch) {
    if (!ch) return;
    for (const c of ch) {
      if (c.data && c.data.leaf) genes.push(c.data.gene_id || c.data.name);
      else { walk(c.children); walk(c._children); }
    }
  })(children);
  return genes;
}

function showCollapsedPopup(event, d) {
  event.stopPropagation();
  cpActiveNode = d;
  const currentName = customNodeNames[d._uid] || d.data.name || collapsedLabel(d);
  const genes = collectLeafGenes(d._children);
  document.getElementById("cp-name").value = currentName;
  document.getElementById("cp-genes-label").textContent = genes.length + " genes:";
  document.getElementById("cp-genes").value = genes.join("\n");
  cpEl.style.display = "block";
  const x = Math.min(event.clientX + 12, window.innerWidth  - cpEl.offsetWidth  - 10);
  const y = Math.min(event.clientY + 12, window.innerHeight - cpEl.offsetHeight - 10);
  cpEl.style.left = Math.max(4, x) + "px";
  cpEl.style.top  = Math.max(4, y) + "px";
}

function cpClose() { cpEl.style.display = "none"; cpActiveNode = null; }

function cpRename() {
  if (!cpActiveNode) return;
  const newName = document.getElementById("cp-name").value.trim();
  if (newName) customNodeNames[cpActiveNode._uid] = newName;
  cpClose();
  renderTree(false);
}

function cpCopy() {
  const txt = document.getElementById("cp-genes").value;
  navigator.clipboard ? navigator.clipboard.writeText(txt).catch(()=>{}) : (() => {
    const ta = document.getElementById("cp-genes"); ta.select(); document.execCommand("copy");
  })();
}

function cpExpand() {
  const d = cpActiveNode;
  cpClose();
  if (d && d._children) { d.children = d._children; d._children = null; renderTree(true); }
}

// ═══════════════════════════════════════════════════════════════════════════════
// HEATMAP VIEW
// ═══════════════════════════════════════════════════════════════════════════════
const TOP_MARGIN = 110;

// species order: tree order filtered to species present in data
const dataSpecies = new Set([...FAMILY_DATA,...HG_DATA].flatMap(d=>Object.keys(d.species_counts)));
let speciesOrder = SPECIES_ORDER.filter(s=>dataSpecies.has(s));
for (const s of dataSpecies) { if (!speciesOrder.includes(s)) speciesOrder.push(s); }

// If no heatmap data but we have POSSVM trees, fall back to tree view
const hasHeatmapData = FAMILY_DATA.length > 0 || HG_DATA.length > 0;

// prefix selector
(function(){
  const sel = document.getElementById("prefixSelect");
  const prefs = ["all",...new Set(FAMILY_DATA.map(d=>d.pref).filter(Boolean).sort())];
  prefs.forEach(p => { const o=document.createElement("option"); o.value=p; o.textContent=p; sel.appendChild(o); });
  sel.addEventListener("change", drawHeatmap);
})();

function drawCladogram() {
  const tp = document.getElementById("tree-panel");
  tp.innerHTML = "";
  if (!SP_TREE_DATA || !SP_TREE_DATA.children || !speciesOrder.length) return;

  const W = 270, H = speciesOrder.length*14+TOP_MARGIN+40;
  const svg = d3.select(tp).append("svg").attr("width",W).attr("height",H);
  const leafY = {}; speciesOrder.forEach((s,i)=>{ leafY[s]=TOP_MARGIN+i*14; });

  function clone(n){ return JSON.parse(JSON.stringify(n)); }
  function prune(n){
    if(!n.children) return speciesOrder.includes(n.name)?n:null;
    const k=n.children.map(prune).filter(Boolean);
    if(!k.length) return null;
    if(k.length===1) return k[0];
    n.children=k; return n;
  }
  let tree=prune(clone(SP_TREE_DATA));
  if(!tree) return;

  function assignY(n){ if(!n.children){n._y=leafY[n.name]||0;return n._y;} const ys=n.children.map(assignY); n._y=d3.mean(ys);return n._y; }
  function assignX(n,d=0){ n._d=d; if(n.children) n.children.forEach(c=>assignX(c,d+1)); }
  function flat(n){ return [n].concat(n.children?n.children.flatMap(flat):[]); }
  assignY(tree); assignX(tree);
  const maxD=d3.max(flat(tree),d=>d._d)||1;
  flat(tree).forEach(n=>{ if(!n.children) n._d=maxD; }); // align all tips
  const sx=d=>10+(d/maxD)*150;

  function drawB(n){
    if(!n.children)return;
    const ys=n.children.map(c=>c._y);
    svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(n._d)).attr("y1",d3.min(ys)).attr("y2",d3.max(ys)).attr("stroke","#aaa");
    n.children.forEach(c=>{
      svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(c._d)).attr("y1",c._y).attr("y2",c._y).attr("stroke","#aaa");
      drawB(c);
    });
  }
  drawB(tree);
}

function drawSpeciesTree() {
  const wrap = document.getElementById("sp-tree-wrap");
  wrap.innerHTML = "";
  if (!SP_TREE_DATA || !SP_TREE_DATA.children || !SPECIES_ORDER.length) {
    wrap.innerHTML = '<p style="padding:20px;color:#888">No species tree provided (pass --species_tree).</p>';
    return;
  }

  const rowH = 22, topM = 16, leftM = 16, rightM = 260;
  const W = Math.max(Math.floor((wrap.clientWidth || 900) * 0.5), 420);
  const treeW = W - leftM - rightM;
  const inPhylo = new Set(ALL_SPECIES);
  const tipColor = n => inPhylo.has(n.name) ? "#111" : "#bbb";

  // ── helpers ──────────────────────────────────────────────────────────
  function clone(n){ return JSON.parse(JSON.stringify(n)); }
  function prune(n){
    if(!n.children) return SPECIES_ORDER.includes(n.name)?n:null;
    const k=n.children.map(prune).filter(Boolean);
    if(!k.length) return null;
    if(k.length===1) return k[0];
    n.children=k; return n;
  }
  function countLeaves(n){ return n.children?n.children.reduce((s,c)=>s+countLeaves(c),0):1; }
  function getLeaves(n){ return n.children?n.children.flatMap(getLeaves):[n]; }
  function flat(n){ return [n].concat(n.children?n.children.flatMap(flat):[]); }

  const tree = prune(clone(SP_TREE_DATA));
  if(!tree){ wrap.innerHTML='<p style="padding:20px;color:#888">Species tree is empty.</p>'; return; }

  // ── assign stable IDs and depths ──────────────────────────────────────
  let _uid=0;
  function assignId(n,d=0){ n._id=String(_uid++); n._d=d; if(n.children) n.children.forEach(c=>assignId(c,d+1)); }
  assignId(tree);
  const allNodes = flat(tree);
  const maxD = d3.max(allNodes,n=>n._d)||1;
  allNodes.forEach(n=>{ if(!n.children) n._d=maxD; }); // align tips
  const sx = d => leftM + (d/maxD)*treeW;

  // ── layout: assign _y accounting for collapsed blocks ─────────────────
  function layout(n, yStart){
    if(spCollapsed.has(n._id)){
      n._colH = countLeaves(n)*rowH;
      n._colY0 = yStart;
      n._y = yStart + n._colH/2;
      n._isCol = true;
      return yStart + n._colH;
    }
    n._isCol = false;
    if(!n.children){ n._y = yStart + rowH/2; return yStart + rowH; }
    let y = yStart;
    for(const c of n.children) y = layout(c, y);
    n._y = d3.mean(n.children.map(c=>c._y));
    return y;
  }
  const totalH = layout(tree, topM);
  const svg = d3.select(wrap).append("svg").attr("width",W).attr("height",totalH+16);

  // ── branches (skip into collapsed subtrees) ───────────────────────────
  function drawBranches(n){
    if(n._isCol||!n.children) return;
    const ys=n.children.map(c=>c._y);
    svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(n._d))
      .attr("y1",d3.min(ys)).attr("y2",d3.max(ys)).attr("stroke","#999").attr("stroke-width",1.5);
    n.children.forEach(c=>{
      svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(c._d))
        .attr("y1",c._y).attr("y2",c._y).attr("stroke","#999").attr("stroke-width",1.5);
      drawBranches(c);
    });
  }
  drawBranches(tree);

  // ── collapsed triangles ───────────────────────────────────────────────
  function drawCollapsed(n){
    if(!n.children) return;
    if(n._isCol){
      const x0=sx(n._d), x1=sx(maxD);
      const y0=n._colY0, y1=n._colY0+n._colH, yM=n._y;
      const leaves=getLeaves(n);
      // horizontal branch from parent edge to apex already drawn; draw triangle
      svg.append("polygon")
        .attr("points",`${x0},${yM} ${x1},${y0} ${x1},${y1}`)
        .attr("fill","rgba(90,130,170,0.18)").attr("stroke","#7a9ab8").attr("stroke-width",1)
        .style("cursor","pointer")
        .on("click",()=>{ spCollapsed.delete(n._id); drawSpeciesTree(); });
      // node name badge at apex
      if(n.name){
        svg.append("text").attr("x",x0+4).attr("y",yM-3)
          .attr("font-size",8).attr("fill","#5a7090").attr("font-style","italic")
          .text(n.name+" ["+leaves.length+"]");
      }
      // species names at triangle base
      const nameH=n._colH/leaves.length;
      leaves.forEach((l,i)=>{
        const ly=n._colY0+i*nameH+nameH/2;
        const tc=tipColor(l);
        svg.append("circle").attr("cx",x1+4).attr("cy",ly).attr("r",3)
          .attr("fill",tc);
        svg.append("text").attr("x",x1+10).attr("y",ly).attr("dy","0.35em")
          .attr("font-size",Math.max(7,Math.min(11,nameH-2)))
          .attr("fill",tc).attr("font-family","monospace").text(l.name);
      });
      return;
    }
    n.children.forEach(drawCollapsed);
  }
  drawCollapsed(tree);

  // ── internal nodes (click to collapse) ───────────────────────────────
  function drawInternals(n){
    if(n._isCol||!n.children) return;
    if(n.name){
      svg.append("text").attr("x",sx(n._d)+6).attr("y",n._y-4)
        .attr("font-size",9).attr("fill","#777").attr("font-style","italic").text(n.name);
    }
    // clickable node dot — only non-root internal nodes
    if(n._id!==tree._id){
      svg.append("circle").attr("cx",sx(n._d)).attr("cy",n._y).attr("r",5)
        .attr("fill","#fff").attr("stroke","#999").attr("stroke-width",1.2)
        .style("cursor","pointer")
        .attr("title","Click to collapse")
        .on("click",()=>{ spCollapsed.add(n._id); drawSpeciesTree(); });
    }
    n.children.forEach(drawInternals);
  }
  drawInternals(tree);

  // ── leaf tips ─────────────────────────────────────────────────────────
  function drawLeaves(n){
    if(n._isCol) return;
    if(!n.children){
      const tc=tipColor(n);
      svg.append("circle").attr("cx",sx(n._d)).attr("cy",n._y).attr("r",4)
        .attr("fill",tc).attr("stroke","#fff").attr("stroke-width",0.5);
      svg.append("text").attr("x",sx(n._d)+8).attr("y",n._y).attr("dy","0.35em")
        .attr("font-size",12).attr("fill",tc).attr("font-family","monospace").text(n.name);
      return;
    }
    n.children.forEach(drawLeaves);
  }
  drawLeaves(tree);
}

function drawHeatmap() {
  const prefix = document.getElementById("prefixSelect").value;
  const back   = document.getElementById("hm-back");
  const crumb  = document.getElementById("hm-breadcrumb");
  let data, colLabel, clickHandler;

  if (hmViewMode==="class") {
    // aggregate FAMILY_DATA by class
    const map={};
    FAMILY_DATA.filter(d=>prefix==="all"||d.pref===prefix).forEach(d=>{
      const cls=d.class||d.pref||"other";
      if(!map[cls]) map[cls]={id:cls,class:cls,species_counts:{}};
      for(const [sp,n] of Object.entries(d.species_counts))
        map[cls].species_counts[sp]=(map[cls].species_counts[sp]||0)+n;
    });
    data=Object.values(map).sort((a,b)=>
      Object.values(b.species_counts).reduce((x,y)=>x+y,0)-
      Object.values(a.species_counts).reduce((x,y)=>x+y,0));
    back.style.display="none"; crumb.textContent="";
    colLabel=d=>d.id;
    clickHandler=(_ev,d)=>{ hmViewMode="family"; hmActiveClass=d.id; drawHeatmap(); };

  } else if (hmViewMode==="family") {
    data=FAMILY_DATA.filter(d=>(d.class||d.pref)===hmActiveClass)
                    .filter(d=>prefix==="all"||d.pref===prefix)
                    .sort((a,b)=>b.total-a.total);
    back.style.display="inline"; crumb.textContent=hmActiveClass;
    colLabel=d=>d.family;
    clickHandler=(_ev,d)=>{ hmViewMode="hg"; hmActiveFamily=d.family; drawHeatmap(); };

  } else if (hmViewMode==="hg") {
    data=HG_DATA.filter(d=>d.family===hmActiveFamily)
                .filter(d=>prefix==="all"||d.pref===prefix)
                .sort((a,b)=>b.total-a.total);
    back.style.display="inline"; crumb.textContent=hmActiveClass+" \u203a "+hmActiveFamily;
    colLabel=d=>d.hg||d.id;
    clickHandler=(_ev,d)=>{ hmViewMode="og"; hmActiveHG=d.id; hmActiveHGRec=d; drawHeatmap(); };

  } else { // og
    // robust lookup: try id first, fall back to family+hg match
    const r=hmActiveHGRec;
    const treeRec=TREE_INDEX.find(t=>t.id===hmActiveHG)
      ||(r&&TREE_INDEX.find(t=>t.family===r.family&&t.hg===r.hg));
    const detail=treeRec?loadDetail(treeRec.id):null;
    const ogs=detail&&detail.ogs||{};
    data=Object.entries(ogs).map(([og,gids])=>{
      const sc={};
      gids.forEach(g=>{ const sp=getSpeciesPfx(g); if(sp) sc[sp]=(sc[sp]||0)+1; });
      return {id:og,species_counts:sc,total:gids.length};
    }).sort((a,b)=>b.total-a.total);
    back.style.display="inline"; crumb.textContent=hmActiveClass+" \u203a "+hmActiveFamily+" \u203a "+(r&&r.hg||hmActiveHG);
    colLabel=d=>d.id;
    clickHandler=(_ev,_d)=>{
      if(treeRec){ switchTab("trees"); selectTree(treeRec); renderSidebar(""); }
    };
  }

  const cW=18, cH=12;
  const maxNameLen=speciesOrder.reduce((m,s)=>Math.max(m,s.length),0);
  const ROW_LABEL_W=Math.max(110,Math.min(200,maxNameLen*7+14));
  const svgW=data.length*cW+ROW_LABEL_W+20, svgH=speciesOrder.length*14+TOP_MARGIN+40;
  const panel=document.getElementById("heatmap-panel");
  const svg=d3.select(panel).html("").append("svg").attr("width",svgW).attr("height",svgH);

  if(!data.length){ svg.append("text").attr("x",40).attr("y",TOP_MARGIN+30).attr("fill","#999").text("No data for this selection."); return; }

  const zMat=data.map(rec=>{
    const vals=speciesOrder.map(s=>rec.species_counts[s]||0);
    const m=d3.mean(vals), sd=d3.deviation(vals)||1;
    return vals.map(v=>(v-m)/sd);
  });
  const zAll=zMat.flat();
  const zMax=Math.max(Math.abs(d3.min(zAll)||0),Math.abs(d3.max(zAll)||0),0.01);
  const color=d3.scaleDiverging().domain([zMax,0,-zMax]).interpolator(d3.interpolateRdBu);

  data.forEach((rec,ci)=>{
    speciesOrder.forEach((sp,ri)=>{
      const count=rec.species_counts[sp]||0;
      const z=zMat[ci][ri];
      svg.append("rect").attr("class","hm-cell")
        .attr("x",ci*cW+ROW_LABEL_W).attr("y",ri*14+TOP_MARGIN).attr("width",cW-2).attr("height",cH)
        .attr("fill",color(z)).style("cursor","pointer")
        .on("mouseover",ev=>{
          showTip(ev,'<b>'+sp+'</b><br>'+rec.id+'<br>count: <b>'+count+'</b><br>z: <b>'+z.toFixed(2)+'</b>');
        })
        .on("mousemove",moveTip).on("mouseout",hideTip)
        .on("click",ev=>clickHandler(ev,rec));
    });
  });

  // row labels
  speciesOrder.forEach((sp,ri)=>{
    svg.append("text")
      .attr("x",ROW_LABEL_W-4).attr("y",ri*14+TOP_MARGIN+9)
      .attr("text-anchor","end").attr("font-size",11).attr("fill","#333")
      .text(sp);
  });

  // column headers
  svg.selectAll("text.hm-col").data(data).enter().append("text").attr("class","hm-col")
    .attr("transform",(d,i)=>`translate(${i*cW+ROW_LABEL_W},${TOP_MARGIN-10}) rotate(-65)`)
    .attr("font-size",9).style("cursor","pointer")
    .text(colLabel)
    .on("click",clickHandler);

  // hint
  svg.append("text")
    .attr("x",ROW_LABEL_W).attr("y",TOP_MARGIN-92)
    .attr("font-size",9).attr("fill","#aaa").attr("font-style","italic")
    .text("Click a column label or cell to drill down \u2193");
}

function hmBack(){
  if(hmViewMode==="og")     { hmViewMode="hg";     hmActiveHG=null; hmActiveHGRec=null; }
  else if(hmViewMode==="hg"){ hmViewMode="family";  hmActiveFamily=null; hmActiveHG=null; hmActiveHGRec=null; }
  else                      { hmViewMode="class";   hmActiveClass=null; hmActiveFamily=null; hmActiveHG=null; hmActiveHGRec=null; }
  drawHeatmap();
}

// ═══════════════════════════════════════════════════════════════════════════════
// TREE VIEW – SIDEBAR
// ═══════════════════════════════════════════════════════════════════════════════
let currentIndex  = null;
let currentDetail = null;
let treeSource    = "generax";   // "generax" | "original"

function loadDetail(id) {
  const el=document.getElementById("treedata-"+id);
  if(!el)return null;
  try{ return JSON.parse(el.textContent); }catch(e){ return null; }
}

function groupByFamily(records){
  const g={};
  for(const r of records){ const f=r.family||"(other)"; (g[f]=g[f]||[]).push(r); }
  return g;
}

// populate class filter once
(function(){
  const sel=document.getElementById("class-filter");
  const classes=[...new Set(TREE_INDEX.map(r=>r.class||"").filter(Boolean))].sort();
  classes.forEach(c=>{ const o=document.createElement("option"); o.value=c; o.textContent=c; sel.appendChild(o); });
  sel.addEventListener("change",()=>renderSidebar(document.getElementById("hg-search").value));
})();

function renderSidebar(filter){
  const lc=(filter||"").toLowerCase().trim();
  const cls=document.getElementById("class-filter").value;
  const list=document.getElementById("hg-list"); list.innerHTML="";
  const subset=cls?TREE_INDEX.filter(r=>(r.class||"")=== cls):TREE_INDEX;
  const groups=groupByFamily(subset);
  let total=0;
  for(const [fam,recs] of Object.entries(groups).sort()){
    const matching=lc?recs.filter(r=>r.hg.toLowerCase().includes(lc)||r.family.toLowerCase().includes(lc)||(r.og_names||[]).some(o=>o.toLowerCase().includes(lc))):recs;
    if(!matching.length)continue;
    total+=matching.length;
    const forceOpen=lc.length>0;
    const gDiv=document.createElement("div");
    const hdr=document.createElement("div"); hdr.className="fam-header"+(forceOpen?" open":"");
    hdr.innerHTML='<span class="fam-arrow">\u25b6</span><span>'+fam+'</span><span style="font-weight:400;color:#888;margin-left:auto">'+matching.length+'</span>';
    const body=document.createElement("div"); body.className="fam-body"+(forceOpen?" open":"");
    hdr.addEventListener("click",()=>{ hdr.classList.toggle("open"); body.classList.toggle("open",hdr.classList.contains("open")); });
    gDiv.appendChild(hdr);
    for(const rec of matching){
      const item=document.createElement("div");
      item.className="hg-item"+(currentIndex&&currentIndex.id===rec.id?" selected":"");
      item.innerHTML='<div class="hg-name">'+rec.hg+'</div><div class="hg-meta">'+rec.n_leaves+' genes \u00b7 '+rec.n_ogs+' OGs</div>';
      item.addEventListener("click",()=>selectTree(rec));
      body.appendChild(item);
    }
    gDiv.appendChild(body); list.appendChild(gDiv);
  }
  document.getElementById("hg-count").textContent=total+" shown";
}

document.getElementById("hg-search").addEventListener("input",function(){ renderSidebar(this.value); });

// ═══════════════════════════════════════════════════════════════════════════════
// TREE VIEW – COLOUR-BY & HIGHLIGHT
// ═══════════════════════════════════════════════════════════════════════════════
function populateColorBy(){
  const sel=document.getElementById("color-by");
  while(sel.options.length>1)sel.remove(1);
  if(!currentIndex)return;
  // orthogroup option (always available when a tree is loaded)
  const oOg=document.createElement("option"); oOg.value="og"; oOg.textContent="by orthogroup"; sel.appendChild(oOg);
  const tSp=new Set(currentIndex.species);
  CLADE_DATA.forEach((cd,i)=>{
    const grps=new Set();
    for(const sp of tSp){ const g=cd.groups[sp]; if(g)grps.add(g); }
    if(grps.size<2)return;
    const o=document.createElement("option"); o.value=i; o.textContent=cd.name+" ("+grps.size+" groups)";
    sel.appendChild(o);
  });
}

document.getElementById("color-by").addEventListener("change",function(){
  const v=this.value;
  if(v==="species"){
    colorMode="species"; cladeSp2Color={}; cladeSp2Group={}; cladeGrpColor={}; ogLeaf2Color={}; ogName2Color={};
    // expand so leaf colours are visible
    if(rootNode) rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;}});
  } else if(v==="og"){
    colorMode="og"; cladeSp2Color={}; cladeSp2Group={}; cladeGrpColor={};
    rebuildOgColors();
    // re-collapse to OG level
    collapseToOGs(); return;
  } else {
    colorMode="clade"; ogLeaf2Color={}; ogName2Color={};
    const cd=CLADE_DATA[+v];
    const tSp=currentIndex?currentIndex.species:[];
    const grps=new Set(); for(const sp of tSp){const g=cd.groups[sp];if(g)grps.add(g);}
    const sorted=[...grps].sort(); cladeGrpColor={};
    sorted.forEach((g,i)=>{cladeGrpColor[g]=palette[i%palette.length];});
    cladeSp2Color={}; cladeSp2Group={};
    for(const [sp,grp] of Object.entries(cd.groups)){
      cladeSp2Color[sp]=cladeGrpColor[grp]||"#ccc";
      cladeSp2Group[sp]=grp;
    }
    // expand so leaf colours are visible
    if(rootNode) rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;}});
  }
  renderTree(true);
});

function populateDatalist(){
  const dl=document.getElementById("hl-list"); dl.innerHTML="";
  for(const cd of CLADE_DATA){
    const o=document.createElement("option"); o.value=cd.name; dl.appendChild(o);
  }
  if(currentIndex){
    for(const sp of currentIndex.species){
      const o=document.createElement("option"); o.value=sp; dl.appendChild(o);
    }
  }
}

function resolveQuery(query){
  const lq=(query||"").toLowerCase().trim();
  if(!lq) return new Set();
  const clade=CLADE_DATA.find(c=>c.name.toLowerCase()===lq)||CLADE_DATA.find(c=>c.name.toLowerCase().includes(lq));
  if(clade) return new Set(Object.keys(clade.groups));
  const sp=new Set();
  if(currentIndex) for(const s of currentIndex.species){ if(s.toLowerCase().includes(lq)) sp.add(s); }
  ALL_SPECIES.forEach(s=>{ if(s.toLowerCase().includes(lq)) sp.add(s); });
  return sp;
}

function rebuildHlSet(){
  hlGroupIndex=new Map();
  if(!hlQueries.length){ hlSet=null; }
  else {
    const union=new Set();
    hlQueries.forEach((q,i)=>{ resolveQuery(q).forEach(s=>{ union.add(s); if(!hlGroupIndex.has(s)) hlGroupIndex.set(s,i); }); });
    hlSet=union.size?union:null;
  }
  renderHlTags();
  document.getElementById("btn-focus-hl").style.display=hlSet?"inline":"none";
  if(currentIndex) renderTree();
}

function renderHlTags(){
  const el=document.getElementById("hl-tags"); el.innerHTML="";
  hlQueries.forEach((q,i)=>{
    const chip=document.createElement("span"); chip.className="hl-tag";
    chip.style.background=hlTagColors[i%hlTagColors.length];
    const lbl=document.createTextNode(q+" ");
    const x=document.createElement("span"); x.className="hl-tag-x"; x.textContent="\u00d7";
    x.onclick=()=>removeHlTag(i);
    chip.appendChild(lbl); chip.appendChild(x);
    el.appendChild(chip);
  });
}

function addHlTag(query){
  query=(query||"").trim();
  if(!query||hlQueries.includes(query)) return;
  hlQueries.push(query);
  document.getElementById("hl-search").value="";
  rebuildHlSet();
}

function removeHlTag(i){ hlQueries.splice(i,1); rebuildHlSet(); }

function clearHighlight(){ hlQueries=[]; document.getElementById("hl-search").value=""; rebuildHlSet(); }

document.getElementById("hl-search").addEventListener("keydown",function(e){
  if(e.key==="Enter"&&this.value.trim()){ addHlTag(this.value); e.preventDefault(); }
});
document.getElementById("hl-search").addEventListener("change",function(){
  if(this.value.trim()) addHlTag(this.value);
});

document.getElementById("tip-font-slider").addEventListener("input",function(){
  tipFontSize=+this.value;
  document.getElementById("tip-font-val").textContent=this.value;
  if(currentIndex){ renderTree(); requestAnimationFrame(fitTree); }
});

document.getElementById("chk-geneid").addEventListener("change",function(){ showGeneId=this.checked; if(currentIndex) renderTree(); });
document.getElementById("chk-og").addEventListener("change",function(){ showOGName=this.checked; if(currentIndex) renderTree(); });
document.getElementById("chk-ref").addEventListener("change",function(){ showRefOrtho=this.checked; if(currentIndex) renderTree(); });


// ═══════════════════════════════════════════════════════════════════════════════
// TREE VIEW – SELECTION & D3 RENDERING
// ═══════════════════════════════════════════════════════════════════════════════
function selectTree(rec){
  currentIndex=rec;
  currentDetail=loadDetail(rec.id);
  if(!currentDetail){ console.warn("No detail data for",rec.id); return; }
  // reset colour state – default to OG colouring
  hlSet=null; hlQueries=[]; hlGroupIndex=new Map();
  document.getElementById("hl-search").value="";
  renderHlTags();
  cladeSp2Color={}; cladeSp2Group={}; cladeGrpColor={}; ogLeaf2Color={}; ogName2Color={}; ogGene2Name={};
  colorMode="og";

  // reset tree source toggle
  treeSource="generax";
  const toggle=document.getElementById("tree-toggle");
  if(currentDetail.prev_tree){
    toggle.style.display="inline";
    toggle.textContent="Showing: GeneRax";
  } else {
    toggle.style.display="none";
  }

  renderSidebar(document.getElementById("hg-search").value);
  populateColorBy();
  document.getElementById("color-by").value="og"; // set after options are added
  populateDatalist();
  document.getElementById("tree-title").textContent=rec.id+" \u00b7 "+rec.n_leaves+" genes";
  document.getElementById("n-ogs-label").textContent=rec.n_ogs+" orthogroups";

  // build OG colour maps (currentDetail already set above)
  const _ogs=currentDetail.ogs||{};
  Object.keys(_ogs).sort().forEach((og,i)=>{
    const col=palette[i%palette.length]; ogName2Color[og]=col;
    for(const gid of _ogs[og]){ ogLeaf2Color[gid]=col; ogGene2Name[gid]=og; }
  });
  drawGeneTree(currentDetail.tree);
  collapseToOGs();
  setTimeout(fitTree, 260);
}

function rebuildOgColors(){
  if(colorMode!=="og")return;
  ogLeaf2Color={}; ogName2Color={}; ogGene2Name={};
  const ogs=activeOgs();
  Object.keys(ogs).sort().forEach((og,i)=>{
    const col=palette[i%palette.length]; ogName2Color[og]=col;
    for(const gid of ogs[og]){ ogLeaf2Color[gid]=col; ogGene2Name[gid]=og; }
  });
}
function toggleTreeSource(){
  if(!currentDetail||!currentDetail.prev_tree) return;
  const toggle=document.getElementById("tree-toggle");
  if(treeSource==="generax"){
    treeSource="original";
    toggle.textContent="Showing: Original (bootstrap)";
    rebuildOgColors();
    drawGeneTree(currentDetail.prev_tree);
  } else {
    treeSource="generax";
    toggle.textContent="Showing: GeneRax";
    rebuildOgColors();
    drawGeneTree(currentDetail.tree);
  }
}

function activeOgs(){
  return (treeSource==="original"&&currentDetail&&currentDetail.prev_ogs)
    ? currentDetail.prev_ogs : (currentDetail?currentDetail.ogs:{});
}

const treeSvg = d3.select("#tree-svg");
let rootNode=null, gMain=null, _uid=0, _zoom=null;
let useBranchLen=false, _phyloScale=1;

function isOGNode(d){ return !d.data.leaf && d.data.name && activeOgs()[d.data.name]!==undefined; }

/** Return the MRCA node for a list of d3-hierarchy leaf nodes. */
function findMRCA(leaves){
  if(!leaves.length) return null;
  if(leaves.length===1) return leaves[0].parent||leaves[0];
  const pathsToRoot=leaves.map(l=>{const p=[];let n=l;while(n){p.push(n);n=n.parent;}return p;});
  const sets=pathsToRoot.map(p=>new Set(p));
  for(const anc of pathsToRoot[0]){if(sets.every(s=>s.has(anc)))return anc;}
  return null;
}

function countDescLeaves(children){
  if(!children)return 0; let n=0;
  for(const c of children){
    if(c.data.leaf)n++;
    else n+=countDescLeaves(c.children)+countDescLeaves(c._children);
  }
  return n;
}

function collapsedLabel(d){
  const n=countDescLeaves(d._children);
  if(d._uid&&customNodeNames[d._uid]) return customNodeNames[d._uid]+" ["+n+"]";
  const lbl=d.data._og_label||d.data.name||"";
  if(lbl) return lbl+" ["+n+"]";
  const sc={};
  (function cnt(ch){ if(!ch)return; for(const c of ch){
    if(c.data.leaf){const sp=c.data.species||"?";sc[sp]=(sc[sp]||0)+1;}
    else{cnt(c.children);cnt(c._children);}
  }})(d._children);
  const top=Object.entries(sc).sort((a,b)=>b[1]-a[1]).slice(0,3)
    .map(([sp,c])=>c>1?sp+"\u00d7"+c:sp);
  return top.join(",")+"\u00a0["+n+"]";
}

/** Assign cumulative branch-length positions (_by) to all nodes. */
function assignBranchLenPos(node, cum){
  node._by=cum;
  (node.children||[]).forEach(c=>assignBranchLenPos(c, cum+(c.data.dist||0)));
  (node._children||[]).forEach(c=>assignBranchLenPos(c, cum+(c.data.dist||0)));
}

/** Horizontal (x) pixel position of a node, honouring branch-length mode. */
function nodeX(d, mg){ return (useBranchLen?(d._by||0)*_phyloScale:d.y)+mg.left; }

/** Rounded-corner elbow path for a cladogram/phylogram link.
 *  Correct cladogram shape: vertical stem first, then horizontal to child. */
function elbowPath(s, t, mg, r){
  const sx=nodeX(s,mg), sy=s.x+mg.top;
  const tx=nodeX(t,mg), ty=t.x+mg.top;
  const dy=ty-sy;
  if(Math.abs(dy)<r*2) return `M${sx},${sy}V${ty}H${tx}`;
  return dy>0
    ? `M${sx},${sy}V${ty-r}A${r},${r} 0 0,1 ${sx+r},${ty}H${tx}`
    : `M${sx},${sy}V${ty+r}A${r},${r} 0 0,0 ${sx+r},${ty}H${tx}`;
}

function toggleLengths(){
  useBranchLen=!useBranchLen;
  const btn=document.getElementById("btn-lengths");
  btn.classList.toggle("active-btn", useBranchLen);
  if(rootNode){ assignBranchLenPos(rootNode,0); renderTree(false); }
}

function fitTree(){
  if(!gMain||!_zoom) return;
  const wrap=document.getElementById("tree-wrap");
  const W=wrap.clientWidth||800, H=wrap.clientHeight||600;
  const b=gMain.node().getBBox();
  if(!b.width||!b.height) return;
  const pad=24;
  const sc=Math.min((W-pad*2)/b.width,(H-pad*2)/b.height,2);
  const tx=(W-b.width*sc)/2-b.x*sc, ty=(H-b.height*sc)/2-b.y*sc;
  treeSvg.transition().duration(350)
    .call(_zoom.transform, d3.zoomIdentity.translate(tx,ty).scale(sc));
}

function drawScaleBar(mg, iW, tH){
  const maxBL=iW/_phyloScale;
  const raw=maxBL*0.15;
  const mag=Math.pow(10,Math.floor(Math.log10(raw||1)));
  const nice=[1,2,5].map(x=>x*mag).find(x=>x>=raw*0.5)||mag;
  const barPx=nice*_phyloScale;
  const x0=mg.left+10, y0=mg.top+tH+14;
  const g=gMain.append("g").attr("class","scale-bar-g");
  g.append("line").attr("x1",x0).attr("y1",y0).attr("x2",x0+barPx).attr("y2",y0);
  g.append("line").attr("x1",x0).attr("y1",y0-4).attr("x2",x0).attr("y2",y0+4);
  g.append("line").attr("x1",x0+barPx).attr("y1",y0-4).attr("x2",x0+barPx).attr("y2",y0+4);
  g.append("text").attr("x",x0+barPx/2).attr("y",y0+13).text(nice.toPrecision(2));
}

function drawGeneTree(treeData){
  treeSvg.selectAll("*").remove(); _uid=0;
  const wrap=document.getElementById("tree-wrap");
  const W=wrap.clientWidth||800, H=wrap.clientHeight||600;
  treeSvg.attr("width",W).attr("height",H);
  _zoom=d3.zoom().scaleExtent([0.03,30]).on("zoom",e=>{gMain.attr("transform",e.transform);});
  treeSvg.call(_zoom).on("dblclick.zoom",null);
  gMain=treeSvg.append("g");
  rootNode=d3.hierarchy(treeData,d=>d.children);
  rootNode.each(d=>{d._uid=++_uid;});
  if(useBranchLen) assignBranchLenPos(rootNode,0);
  renderTree(false);
  // auto-fit after layout
  setTimeout(fitTree, 260);
}

const BADGE_W=96, BADGE_H=17;

function renderTree(animate){
  if(!rootNode) return;
  const wrap=document.getElementById("tree-wrap");
  const W=wrap.clientWidth||800, H=wrap.clientHeight||600;
  const mg={top:20,right:240,bottom:36,left:36};
  const iW=W-mg.left-mg.right, iH=H-mg.top-mg.bottom;
  const nVis=rootNode.leaves().length;
  const rowH=Math.max(18,Math.min(32,Math.floor(iH/Math.max(nVis,1))));
  const tH=Math.max(iH,nVis*rowH);
  d3.cluster().size([tH,iW])(rootNode);

  if(useBranchLen){
    assignBranchLenPos(rootNode,0);
    const maxBL=d3.max(rootNode.leaves(),d=>d._by)||1;
    _phyloScale=iW/maxBL;
  }

  const dur=animate?240:0;

  // ── Links ──
  gMain.selectAll(".link")
    .data(rootNode.links(),d=>d.target._uid)
    .join(
      enter=>enter.append("path").attr("class","link")
        .attr("d",d=>elbowPath(d.source,d.source,mg,3)),
      update=>update,
      exit=>exit.transition().duration(dur*0.5).style("opacity",0).remove()
    )
    .transition().duration(dur)
    .attr("d",d=>elbowPath(d.source,d.target,mg,3));

  // ── Scale bar ──
  gMain.selectAll(".scale-bar-g").remove();
  if(useBranchLen) drawScaleBar(mg,iW,tH);

  // ── Nodes ──
  const nodeSel=gMain.selectAll(".node-g")
    .data(rootNode.descendants(),d=>d._uid)
    .join(
      enter=>{
        const g=enter.append("g").attr("class","node-g")
          .attr("transform",d=>{const p=d.parent||d;return `translate(${nodeX(p,mg)},${p.x+mg.top})`;})
          .style("opacity",0);
        g.append("circle");
        g.append("rect").attr("class","badge-bg")
          .attr("y",-BADGE_H/2).attr("height",BADGE_H).attr("rx",BADGE_H/2).attr("ry",BADGE_H/2);
        g.append("text").attr("class","leaf-label");
        g.append("text").attr("class","og-label");
        return g;
      },
      update=>update,
      exit=>exit.transition().duration(dur*0.5)
        .attr("transform",d=>{const p=d.parent||d;return `translate(${nodeX(p,mg)},${p.x+mg.top})`;})
        .style("opacity",0).remove()
    );

  nodeSel.transition().duration(dur)
    .attr("transform",d=>`translate(${nodeX(d,mg)},${d.x+mg.top})`)
    .style("opacity",1);

  // circle (leaves and expanded internals)
  nodeSel.select("circle")
    .attr("display",d=>d._children?"none":null)
    .attr("r",d=>d.data.leaf?4:isOGNode(d)?5.5:2.8)
    .attr("fill",d=>{
      if(d.data.leaf){
        if(colorMode==="og") return ogLeafColor(d.data.gene_id||d.data.name, d.data.species);
        const c=leafColor(d.data.species||"");
        return hlSet&&!hlSet.has(d.data.species||"")?"#e8e8e8":c;
      }
      if(isOGNode(d))return "#e74c3c";
      return "#ccc";
    })
    .attr("stroke",d=>d.data.leaf?"none":isOGNode(d)?"#b03a2e":"#aaa")
    .attr("stroke-width",d=>isOGNode(d)?1.8:0.8)
    .on("click",(event,d)=>{
      if(d.data.leaf)return; event.stopPropagation();
      if(d._children){ showCollapsedPopup(event,d); }
      else if(d.children){d._children=d.children;d.children=null; renderTree(true);}
    })
    .on("mouseover",showTip).on("mousemove",moveTip).on("mouseout",hideTip);

  // badge rectangle (collapsed nodes)
  nodeSel.select(".badge-bg")
    .attr("display",d=>d._children?null:"none")
    .attr("x",5).attr("width",BADGE_W)
    .on("click",(event,d)=>{
      if(d.data.leaf)return; event.stopPropagation();
      if(d._children) showCollapsedPopup(event,d);
    })
    .on("mouseover",showTip).on("mousemove",moveTip).on("mouseout",hideTip);

  // leaf labels: gene_id + OG name
  nodeSel.select(".leaf-label")
    .attr("x",7).attr("dy","0.32em").attr("text-anchor","start")
    .attr("font-size",d=>d.data.leaf?(tipFontSize!==null?tipFontSize:Math.min(11,rowH-2)):0)
    .attr("fill",d=>colorMode==="og"?ogLeafColor(d.data.gene_id||d.data.name,d.data.species):leafColor(d.data.species||""))
    .attr("display",d=>d.data.leaf?null:"none")
    .text(d=>{
      if(!d.data.leaf) return "";
      const gid=d.data.gene_id||d.data.name;
      const og=d.data.og||ogGene2Name[gid]||"";
      const ref=d.data.ref||"";
      const parts=[];
      if(showGeneId) parts.push(gid);
      if(showOGName && og) parts.push(og);
      if(showRefOrtho && ref) parts.push(ref);
      return parts.join(" \u00b7 ");
    });

  // OG labels (inside badge when collapsed, beside node when OG-named internal)
  nodeSel.select(".og-label")
    .attr("x",d=>d._children?BADGE_W/2+5:-7)
    .attr("dy","0.35em")
    .attr("text-anchor",d=>d._children?"middle":"end")
    .attr("font-size",9)
    .attr("display",d=>(!d.data.leaf&&(isOGNode(d)||d._children))?null:"none")
    .text(d=>d._children?collapsedLabel(d):(isOGNode(d)?d.data.name:""));
}

// ── tree controls ──
function expandAll(){
  if(!rootNode)return;
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;}});
  renderTree(true);
}
function collapseToOGs(){
  if(!rootNode)return;
  // pass 1: expand all
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;}});
  // pass 2a: collapse named OG internal nodes (POSSVM trees with annotated internals)
  let found=false;
  rootNode.each(d=>{
    if(!d.data.leaf&&isOGNode(d)&&d.children){d._children=d.children;d.children=null;found=true;}
  });
  // pass 2b: fallback – derive OGs from pipe-separated tip names (gene_id|og|ref)
  if(!found){
    const ogGroups={};
    rootNode.leaves().forEach(l=>{
      const og=l.data.og;
      if(og)(ogGroups[og]=ogGroups[og]||[]).push(l);
    });
    for(const [og,leaves] of Object.entries(ogGroups)){
      const mrca=findMRCA(leaves);
      if(mrca&&!mrca.data.leaf&&mrca.children){
        mrca.data._og_label=og;
        mrca._children=mrca.children; mrca.children=null;
      }
    }
  }
  renderTree(true);
  setTimeout(fitTree, 260);
}
function collapseAll(){
  if(!rootNode)return;
  rootNode.each(d=>{
    if(d.depth>0&&!d.data.leaf&&d.children){d._children=d.children;d.children=null;}
  });
  if(rootNode._children){rootNode.children=rootNode._children;rootNode._children=null;}
  renderTree(true);
  setTimeout(fitTree, 260);
}

function focusHighlighted(){
  if(!hlSet||!rootNode)return;
  // expand all
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;}});
  // mark nodes that have ≥1 highlighted descendant (post-order)
  const hasHl=new Map();
  rootNode.eachAfter(d=>{
    if(d.data.leaf){ hasHl.set(d,hlSet.has(d.data.species||"")); }
    else { hasHl.set(d,(d.children||[]).some(c=>hasHl.get(c))); }
  });
  // collapse subtrees with no highlighted leaf
  rootNode.each(d=>{
    if(!d.data.leaf&&d.children&&!hasHl.get(d)){d._children=d.children;d.children=null;}
  });
  renderTree(true);
  setTimeout(fitTree, 260);
}

// ═══════════════════════════════════════════════════════════════════════════════
// INIT
// ═══════════════════════════════════════════════════════════════════════════════
// scroll sync between species-tree panel and heatmap panel
(function(){
  const tp=document.getElementById("tree-panel");
  const hp=document.getElementById("heatmap-panel");
  let syncing=false;
  tp.addEventListener("scroll",()=>{ if(!syncing){syncing=true;hp.scrollTop=tp.scrollTop;syncing=false;} });
  hp.addEventListener("scroll",()=>{ if(!syncing){syncing=true;tp.scrollTop=hp.scrollTop;syncing=false;} });
})();

document.getElementById("tree-count").textContent =
  TREE_INDEX.length+" gene tree"+(TREE_INDEX.length!==1?"s":"");

if (hasHeatmapData || TREE_INDEX.length > 0) {
  switchTab(hasHeatmapData ? "heatmap" : "trees");
} else {
  document.getElementById("pane-heatmap").innerHTML =
    '<div style="padding:40px;color:#999;text-align:center">No data found.<br>'+
    'Pass <code>--possvm_dir</code> and/or <code>--search_dir</code> / <code>--cluster_dir</code>.</div>';
}
</script>
</body>
</html>
"""


# ── Argument parsing and main ─────────────────────────────────────────────────

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate interactive HTML report for step2 outputs.")
    p.add_argument("--possvm_dir", default="results/possvm",
                   help="Directory with POSSVM *.ortholog_groups.newick files")
    p.add_argument("--possvm_prev_dir", default=None,
                   help="Directory with POSSVM results on original IQ-TREE2 trees (pre-GeneRax); "
                        "enables toggle between GeneRax and original trees in the report")
    p.add_argument("--search_dir", default="results/search",
                   help="Directory with *.genes.list files (for heatmap)")
    p.add_argument("--cluster_dir", default="results/clusters",
                   help="Directory with *.fasta HG files (for heatmap)")
    p.add_argument("--family_info", default="data/gene_families_searchinfo.csv",
                   help="TSV with family→class mapping (7 columns)")
    p.add_argument("--species_tree", default=None,
                   help="Newick species tree (named internal nodes for clade colouring)")
    p.add_argument("--output", required=True, help="Output HTML file path")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    possvm_dir  = Path(args.possvm_dir)
    search_dir  = Path(args.search_dir)
    cluster_dir = Path(args.cluster_dir)

    # Heatmap data
    family_info    = load_family_info(args.family_info)
    family_records = build_family_records(search_dir, family_info)
    hg_records     = build_hg_records(cluster_dir, family_info)

    # Species tree (for ordering + cladogram)
    species_order: list = []
    tree_dict: dict = {}
    if args.species_tree and Path(args.species_tree).exists():
        species_order, tree_dict = load_tree_data(args.species_tree)

    # Gene tree data (POSSVM)
    if not possvm_dir.exists():
        print(f"WARN: {possvm_dir} does not exist – no gene trees.", file=sys.stderr)
    records, all_species = load_possvm_trees(possvm_dir) if possvm_dir.exists() else ([], [])
    print(f"Loaded {len(records)} gene trees, {len(all_species)} species.", file=sys.stderr)

    # Prev trees (IQ-TREE2 original, pre-GeneRax) — optional
    prev_records: dict = {}
    if args.possvm_prev_dir:
        prev_dir = Path(args.possvm_prev_dir)
        if prev_dir.is_dir():
            prev_list, prev_sp = load_possvm_trees(prev_dir)
            prev_records = {r["id"]: r for r in prev_list}
            all_species = sorted(set(all_species) | set(prev_sp))
            print(f"Loaded {len(prev_records)} prev gene trees (original IQ-TREE2).",
                  file=sys.stderr)
    print(f"Loaded {len(family_records)} families, {len(hg_records)} HGs for heatmap.",
          file=sys.stderr)

    # Clade groupings
    clade_groupings: list = []
    if args.species_tree and Path(args.species_tree).exists():
        clade_groupings = parse_clade_groupings(Path(args.species_tree))
        print(f"Extracted {len(clade_groupings)} clade groupings.", file=sys.stderr)

    # Build lightweight index (no tree/ogs dicts)
    index_records = []
    for rec in records:
        idx = {k: v for k, v in rec.items() if k not in ("tree_dict", "ogs")}
        idx["has_prev"] = rec["id"] in prev_records
        idx["class"] = family_info.get(rec["family"], rec.get("prefix", ""))
        index_records.append(idx)

    # Build per-HG lazy <script> tags
    lazy_parts = []
    for rec in records:
        detail = {"tree": rec["tree_dict"], "ogs": rec["ogs"]}
        prev = prev_records.get(rec["id"])
        if prev:
            detail["prev_tree"] = prev["tree_dict"]
            detail["prev_ogs"]  = prev["ogs"]
        tag_id = _html.escape(rec["id"], quote=True)
        lazy_parts.append(
            f'<script type="application/json" id="treedata-{tag_id}">'
            + json.dumps(detail, separators=(",", ":"))
            + "</script>"
        )
    lazy_scripts = "\n".join(lazy_parts)

    html = (HTML_TEMPLATE
            .replace("%%LAZY_SCRIPTS%%",    lazy_scripts)
            .replace("%%SPECIES_ORDER%%",   json.dumps(species_order))
            .replace("%%TREE_DATA%%",       json.dumps(tree_dict))
            .replace("%%FAMILY_DATA%%",     json.dumps(family_records))
            .replace("%%HG_DATA%%",         json.dumps(hg_records))
            .replace("%%TREE_INDEX_JSON%%", json.dumps(index_records))
            .replace("%%SPECIES_JSON%%",    json.dumps(all_species))
            .replace("%%CLADE_DATA_JSON%%", json.dumps(clade_groupings)))

    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
