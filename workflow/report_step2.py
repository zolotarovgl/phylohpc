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


def load_family_details(csv_path) -> dict:
    """Return {family_name: {pfam: [str,...], category: str, cls: str}} from gene_families_searchinfo.csv."""
    details = {}
    try:
        with open(csv_path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 7:
                    family   = cols[0].strip()
                    pfam_raw = cols[1].strip()
                    category = cols[5].strip() if len(cols) > 5 else ""
                    cls      = cols[6].strip()
                    if family:
                        pfam_list = [p.strip() for p in pfam_raw.split(",") if p.strip()] if pfam_raw else []
                        details[family] = {"pfam": pfam_list, "category": category, "cls": cls}
    except (FileNotFoundError, OSError):
        pass
    return details


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
        try:
            sv = float(node.support)
            if sv != 1.0:  # ete3 default is 1.0 when not set
                d["support"] = round(sv, 2)
        except (TypeError, ValueError, AttributeError):
            pass
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
                gene, og = parts[0], parts[1].replace(':', '_')
                if og.lower() == "orthogroup":   # POSSVM catch-all; skip
                    continue
                og_members[og].append(gene)
    except OSError:
        pass
    return dict(og_members)


def load_possvm_trees(possvm_dir: Path, source: str = "generax") -> tuple[list, list]:
    """Return (tree_records, all_species).  source='generax' or 'prev'."""
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
        # Strip GeneRax-specific suffixes so id == "Family.HG"
        for suffix in (".generax.tree", ".generax"):
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
            t = Tree(str(nwk), format=0)
        except Exception as exc:
            print(f"WARN: cannot parse {nwk.name}: {exc}", file=sys.stderr)
            continue

        # Re-root unrooted trees (trifurcating root) via midpoint rooting;
        # already-rooted trees (bifurcating root) are left untouched.
        if len(t.children) > 2:
            try:
                t.set_outgroup(t.get_midpoint_outgroup())
            except Exception as exc:
                print(f"WARN: midpoint rooting failed for {nwk.name}: {exc}", file=sys.stderr)

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
            "source":           source,
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
#pane-sptree{flex-direction:column}
#sptree-controls{display:flex;align-items:center;gap:12px;padding:6px 14px;background:#f5f5f5;border-bottom:1px solid #ddd;flex-shrink:0}
#sp-tree-wrap{flex:1;overflow:auto;background:#fff;padding:16px}
#sp-tree-wrap svg{display:block}

/* ── Heatmap pane ── */
#pane-heatmap{flex-direction:row}
#hm-layout{display:flex;flex:1;overflow:hidden}
#hm-right{display:flex;flex-direction:column;flex:1;overflow:hidden;min-width:0}
#tree-panel{overflow-y:auto;background:#fff}
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
.hg-cov{height:3px;background:#e8e8e8;border-radius:2px;margin:3px 0 1px}
.hg-cov-bar{height:3px;background:#1abc9c;border-radius:2px}
.src-badge{font-size:8px;padding:1px 4px;border-radius:3px;background:#dceeff;color:#2a6fa8;font-weight:600;margin-left:5px;vertical-align:middle;letter-spacing:0.02em}

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
.link{fill:none;stroke:#d5d5d5}
.node-g circle{cursor:pointer;transition:r .12s,fill .12s}
.node-g circle:hover{stroke-width:2.5px !important}
.col-tri{stroke:#222;cursor:pointer}
.col-tri:hover{filter:brightness(0.88)}
.leaf-label{font-family:monospace}
.og-label{fill:#b5371f}
.ctrl-btn.active-btn{background:#d5f5e3;border-color:#1abc9c;color:#1a6b4a}
.scale-bar-g line{stroke:#999;stroke-width:1.5px}
.scale-bar-g text{font-size:9px;fill:#888;text-anchor:middle}

/* mini species tree floating panel */
#mini-sp-panel{position:fixed;display:none;background:#fff;border:1px solid #bbb;border-radius:6px;box-shadow:0 3px 12px rgba(0,0,0,.18);z-index:300;overflow:auto;max-width:560px;max-height:70vh;padding:6px}
#mini-sp-panel .msp-title{font-size:10px;font-weight:700;color:#555;margin-bottom:4px;padding:0 2px}
#mini-sp-panel svg text.msp-node-lbl{cursor:pointer;fill:#2980b9;font-size:10px}
#mini-sp-panel svg text.msp-node-lbl:hover{fill:#e74c3c}
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
  <button id="hm-expand-og" style="display:none;padding:2px 8px;font-size:11px;border:1px solid #27ae60;color:#27ae60;border-radius:3px;background:#fff;cursor:pointer" title="Switch to Custom view with all OGs from all visible HGs" onclick="hmExpandToOGs()">&#43; Expand all to OGs</button>
  <span id="tree-count" style="font-size:11px;color:#95a5a6;display:none"></span>
</div>

<!-- ════════ Body: tab strip + panes ════════ -->
<div id="body-wrap">

  <!-- vertical tab strip -->
  <div id="tab-strip">
    <button class="tab-btn" data-tab="families" onclick="switchTab('families')">&#9783; Families</button>
    <button class="tab-btn active" data-tab="sptree" onclick="switchTab('sptree')">&#10022; Species Tree</button>
    <button class="tab-btn" data-tab="heatmap" onclick="switchTab('heatmap')">&#9639; Counts</button>
    <button class="tab-btn" data-tab="trees" onclick="switchTab('trees')">&#11044; Gene Trees</button>
  </div>

  <!-- ── Heatmap pane ── -->
  <div class="tab-pane" id="pane-heatmap">
    <div id="hm-custom-bar" style="display:none;position:fixed;z-index:500;border:1px solid #c8b96e;border-top:none;background:#fffbf0;padding:6px 12px;gap:6px;flex-wrap:wrap;align-items:center;box-shadow:0 4px 12px rgba(0,0,0,.15);border-radius:0 0 6px 6px">
      <span style="font-size:11px;color:#888;font-weight:600">Custom selection:</span>
      <span id="hm-custom-chips" style="display:flex;gap:4px;flex-wrap:wrap;align-items:center"></span>
      <button onclick="hmCustomClear()" style="margin-left:auto;padding:1px 8px;font-size:10px;border:1px solid #ccc;border-radius:3px;background:#fff;cursor:pointer">Clear all</button>
      <button onclick="hmExitCustom()" style="padding:1px 8px;font-size:10px;border:1px solid #4a90d9;color:#4a90d9;border-radius:3px;background:#fff;cursor:pointer">&#8592; Browse mode</button>
    </div>
    <div id="hm-split-bar" style="display:none;align-items:center;gap:6px;padding:4px 10px;font-size:11px;color:#555;border-bottom:1px solid #eee;background:#fafafa">
      <span style="font-weight:600">Row groups:</span>
      <span id="hm-split-tags"></span>
      <button onclick="clearHmSplits()" style="margin-left:4px;padding:1px 8px;font-size:10px;border:1px solid #ccc;border-radius:3px;background:#fff;cursor:pointer">&#10005; Clear</button>
      <span style="font-size:10px;color:#aaa">(shift+click a node in the Species Tree tab to add a group)</span>
    </div>
    <div id="hm-layout">
      <div style="display:flex;flex-direction:column;width:200px;flex-shrink:0;border-right:1px solid #ccc">
        <div style="padding:4px 8px;background:#fafafa;border-bottom:1px solid #eee;flex-shrink:0">
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Prefix: <select id="prefixSelect" style="font-size:11px;flex:1;padding:2px 4px;border:1px solid #ccc;border-radius:3px"></select>
          </label>
        </div>
        <div id="tree-panel" style="flex:1;overflow-y:auto;background:#fff"></div>
      </div>
      <div id="hm-right" style="display:flex;flex-direction:column;flex:1;overflow:hidden;min-width:0">
        <!-- column-label controls strip: always visible above the heatmap SVG -->
        <div id="hm-col-strip" style="display:flex;align-items:center;gap:10px;padding:3px 10px;background:#fafafa;border-bottom:1px solid #eee;flex-shrink:0;flex-wrap:wrap">
          <div style="position:relative;display:inline-block">
            <input id="hm-text-search" type="text" placeholder="&#128269; Search classes, families, HGs, OGs, genes&hellip;" autocomplete="off"
              style="font-size:11px;padding:2px 8px;border:1px solid #ccc;border-radius:3px;width:210px"
              oninput="hmTextSearchInput(this.value)" onkeydown="hmTextSearchKey(event)">
            <div id="hm-search-dd"
              style="display:none;position:absolute;top:100%;left:0;z-index:520;background:#fff;
                     border:1px solid #ccc;border-radius:0 0 4px 4px;max-height:220px;
                     overflow-y:auto;min-width:290px;box-shadow:0 3px 8px rgba(0,0,0,.12);
                     font-size:11px"></div>
          </div>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Label size:
            <input type="range" id="hm-col-font-slider" min="5" max="16" step="0.5" value="9" style="width:70px;cursor:pointer;accent-color:#4a90d9">
            <span id="hm-col-font-val" style="width:20px;text-align:right">9</span>px
          </label>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Rotation:
            <input type="range" id="hm-col-rot-slider" min="0" max="90" step="5" value="65" style="width:70px;cursor:pointer;accent-color:#4a90d9">
            <span id="hm-col-rot-val" style="width:24px;text-align:right">65</span>&deg;
          </label>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Colour:
            <select id="hm-color-mode" style="font-size:11px;padding:2px 4px;border:1px solid #ccc;border-radius:3px">
              <option value="zscore">Z-score (RdBu)</option>
              <option value="absolute">Absolute count</option>
            </select>
          </label>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Sort cols:
            <select id="hm-col-sort-sp" style="font-size:11px;padding:2px 4px;border:1px solid #ccc;border-radius:3px" onchange="hmColSortSp=this.value||null;hmColOrderOverride=null;drawHeatmap()">
              <option value="">— none —</option>
            </select>
          </label>
          <button onclick="downloadHeatmapPNG()" title="Download heatmap as PNG" style="padding:2px 8px;font-size:11px;border:1px solid #888;color:#555;border-radius:3px;background:#fff;cursor:pointer">&#11015; PNG</button>
          <button onclick="downloadHeatmapSVG()" title="Download heatmap as SVG" style="padding:2px 8px;font-size:11px;border:1px solid #888;color:#555;border-radius:3px;background:#fff;cursor:pointer">&#11015; SVG</button>
          <details style="margin-left:auto;font-size:10px">
            <summary style="cursor:pointer;color:#888;list-style:none">&#9432; Help</summary>
            <div style="position:absolute;z-index:50;background:#fff;border:1px solid #ddd;border-radius:4px;padding:6px 10px;box-shadow:0 2px 8px rgba(0,0,0,.12);line-height:1.7;color:#999;min-width:300px;right:10px">
              <div><b style="color:#777">Navigate:</b> click a <i>column header</i> to drill into family &rarr; HG &rarr; OG</div>
              <div><b style="color:#777">Expand to OGs:</b> at HG level, shift+click a column or use the &ldquo;Expand all to OGs&rdquo; button</div>
              <div><b style="color:#777">Open tree:</b> click a <i>cell</i> to open the gene tree and highlight that species</div>
              <div><b style="color:#777">Sort columns:</b> use the &ldquo;Sort cols&rdquo; dropdown to rank columns by a species&rsquo; gene count</div>
              <div><b style="color:#777">Reorder:</b> drag a column header left or right to reorder columns manually</div>
              <div><b style="color:#777">Row groups:</b> shift+click a node in the Species Tree tab to split rows by clade</div>
              <div><b style="color:#777">Colour:</b> intensity = gene count; grey = absent</div>
              <div><b style="color:#777">Custom:</b> build a heatmap from any OGs, classes or HGs; search by OG name, gene ID or group name</div>
            </div>
          </details>
        </div>
        <div id="heatmap-panel" style="flex:1;overflow:auto;background:#fff"></div>
      </div>
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
          <!-- group 1: tree operations -->
          <button class="ctrl-btn" onclick="expandAll()">Expand all</button>
          <button class="ctrl-btn" onclick="collapseToOGs()">Collapse to OGs</button>
          <button class="ctrl-btn" id="btn-og-labels" onclick="toggleOGLabels()">OG labels</button>
          <button class="ctrl-btn" onclick="collapseAll()">Collapse all</button>
          <span style="border-left:1px solid #ddd;margin:0 3px;height:16px;align-self:center"></span>
          <!-- group 2: view toggles -->
          <button class="ctrl-btn" id="btn-support" onclick="toggleSupport()">Support</button>
          <button class="ctrl-btn" id="btn-lengths" onclick="toggleLengths()">Branch lengths</button>
          <button class="ctrl-btn" onclick="fitTree()">&#x2922; Fit</button>
          <button class="ctrl-btn" onclick="downloadTreeSVG()" title="Download gene tree as SVG">&#11015; SVG</button>
          <button class="ctrl-btn" onclick="downloadTreePNG()" title="Download gene tree as PNG">&#11015; PNG</button>
          <button class="ctrl-btn" id="tree-toggle" style="display:none;background:#e8f0fe;border-color:#4a90d9" onclick="toggleTreeSource()">Showing: GeneRax</button>
          <span id="tree-title"></span>
          <span id="n-ogs-label"></span>
          <span style="flex:1"></span>
          <!-- group 3: visual options -->
          <label style="font-size:11px;color:#555">Color:
            <select id="color-by" style="font-size:11px;padding:2px 4px;border:1px solid #bbb;border-radius:3px">
              <option value="species">by species</option>
            </select>
          </label>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Label:
            <input type="range" id="tip-font-slider" min="6" max="24" step="1" value="11" style="width:60px;cursor:pointer;accent-color:#4a90d9">
            <span id="tip-font-val" style="width:20px;text-align:right">11</span>px
          </label>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Lines:
            <input type="range" id="line-width-slider" min="1" max="8" step="0.5" value="1.3" style="width:60px;cursor:pointer;accent-color:#4a90d9">
            <span id="line-width-val" style="width:20px;text-align:right">1.3</span>px
          </label>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Height:
            <input type="range" id="tree-height-slider" min="0.5" max="5" step="0.1" value="1" style="width:60px;cursor:pointer;accent-color:#4a90d9">
            <span id="tree-height-val" style="width:24px;text-align:right">1</span>&times;
          </label>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Width:
            <input type="range" id="tree-width-slider" min="0.3" max="4" step="0.1" value="1" style="width:60px;cursor:pointer;accent-color:#4a90d9">
            <span id="tree-width-val" style="width:24px;text-align:right">1.0</span>&times;
          </label>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px" title="Opacity of clade highlight rectangles">
            Hl&#945;:
            <input type="range" id="clade-hl-alpha-slider" min="0.05" max="0.8" step="0.05" value="0.22" style="width:55px;cursor:pointer;accent-color:#e67e22">
            <span id="clade-hl-alpha-val" style="width:28px;text-align:right">0.22</span>
          </label>
          <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:4px">
            Collapsed:
            <input type="range" id="collapsed-frac-slider" min="0.1" max="1" step="0.05" value="1" style="width:60px;cursor:pointer;accent-color:#4a90d9">
            <span id="collapsed-frac-val" style="width:28px;text-align:right">1.0</span>
          </label>
          <span style="border-left:1px solid #ddd;margin:0 3px;height:16px;align-self:center"></span>
          <!-- group 4: show checkboxes -->
          <span style="font-size:11px;color:#555;display:flex;align-items:center;gap:6px">
            Show:
            <label style="display:flex;align-items:center;gap:2px;cursor:pointer"><input type="checkbox" id="chk-geneid" checked> ID</label>
            <label style="display:flex;align-items:center;gap:2px;cursor:pointer"><input type="checkbox" id="chk-og" checked> OG</label>
            <label style="display:flex;align-items:center;gap:2px;cursor:pointer"><input type="checkbox" id="chk-ref" checked> ref</label>
            <label style="display:flex;align-items:center;gap:2px;cursor:pointer"><input type="checkbox" id="chk-hide-nonhl"> hide non-hl</label>
            <button id="btn-focus-collapse-style" class="ctrl-btn active-btn" onclick="toggleFocusCollapseStyle()" title="Toggle collapse style used by Focus and hide non-hl: triangle (MRCA) or circle" style="padding:1px 6px;font-size:10px">&#9660; MRCA</button>
          </span>
          <span style="border-left:1px solid #ddd;margin:0 3px;height:16px;align-self:center"></span>
          <!-- group 5: species + OG highlight + species tree -->
          <div id="hl-tags"></div>
          <input id="hl-search" list="hl-list" placeholder="Species… (Enter)">
          <datalist id="hl-list"></datalist>
          <button id="hl-clear" onclick="clearHighlight()" title="Clear all species highlights">&#10005;</button>
          <button class="ctrl-btn" id="btn-mini-sp" onclick="toggleMiniSpPanel(event)" title="Show species tree — click a named node to highlight that clade">&#x1F333; Species tree</button>
          <button class="ctrl-btn" id="btn-focus-hl" onclick="focusHighlighted()" style="display:none" title="Collapse all branches not leading to highlighted tips">Focus</button>
          <span style="border-left:1px solid #ddd;margin:0 3px;height:16px;align-self:center"></span>
          <div id="og-hl-tags"></div>
          <div style="position:relative;display:inline-block">
            <input id="og-hl-search" autocomplete="off" placeholder="OG name… (Enter)" style="width:140px;font-size:11px;padding:3px 6px;border:1px solid #bbb;border-radius:3px"
              oninput="ogHlSearchInput(this.value)" onkeydown="ogHlSearchKey(event)">
            <div id="og-hl-dd" style="display:none;position:absolute;top:100%;left:0;z-index:520;background:#fff;border:1px solid #ccc;border-radius:0 0 4px 4px;max-height:200px;overflow-y:auto;min-width:200px;box-shadow:0 3px 8px rgba(0,0,0,.12);font-size:11px"></div>
          </div>
          <button id="og-hl-clear" onclick="clearOgHighlight()" title="Clear OG highlights">&#10005;</button>
        </div>
        <div id="tree-wrap">
          <svg id="tree-svg"></svg>
        </div>
      </div>
    </div>
  </div>

  <!-- ── Families pane ── -->
  <div class="tab-pane" id="pane-families" style="flex-direction:column;overflow:hidden">
    <div id="fam-controls" style="display:flex;align-items:center;gap:10px;padding:6px 12px;background:#fafafa;border-bottom:1px solid #ddd;flex-shrink:0;flex-wrap:wrap">
      <input id="fam-search" type="text" placeholder="&#128269; Search family, PFAM domain, category&#8230;" autocomplete="off"
        style="font-size:11px;padding:3px 8px;border:1px solid #ccc;border-radius:3px;width:280px"
        oninput="filterFamilyTable(this.value)">
      <span id="fam-count" style="font-size:11px;color:#888"></span>
    </div>
    <div id="fam-table-wrap" style="flex:1;overflow:auto;padding:0 12px 12px">
      <table id="fam-table" style="border-collapse:collapse;width:100%;font-size:12px;margin-top:8px">
        <thead>
          <tr style="background:#f0f0f0;position:sticky;top:0;z-index:10">
            <th class="fam-th" data-col="family"   style="text-align:left;padding:6px 8px;cursor:pointer;white-space:nowrap">Family &#8597;</th>
            <th class="fam-th" data-col="pfam"     style="text-align:left;padding:6px 8px;cursor:pointer">PFAM Domain(s) &#8597;</th>
            <th class="fam-th" data-col="category" style="text-align:left;padding:6px 8px;cursor:pointer;white-space:nowrap">Category &#8597;</th>
            <th class="fam-th" data-col="n_hgs"    style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap">HGs &#8597;</th>
            <th class="fam-th" data-col="n_trees"  style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap" title="HGs with any gene tree">Trees &#8597;</th>
            <th class="fam-th fam-th-generax" data-col="n_generax" style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap" title="HGs with a GeneRax tree">GeneRax &#8597;</th>
            <th class="fam-th" data-col="total"    style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap">Genes &#8597;</th>
            <th class="fam-th" data-col="n_species"style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap">Species &#8597;</th>
          </tr>
        </thead>
        <tbody id="fam-tbody"></tbody>
      </table>
    </div>
  </div>

  <!-- ── Species tree pane ── -->
  <div class="tab-pane active" id="pane-sptree">
    <div id="sptree-controls">
      <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:6px">
        Tree width:
        <input type="range" id="sptree-width-slider" min="20" max="100" step="5" value="50" style="width:90px;cursor:pointer;accent-color:#4a90d9">
        <span id="sptree-width-val">50</span>%
      </label>
      <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:5px">
        Triangle fill:
        <input type="color" id="col-tri-fill" value="#ffffff" style="width:28px;height:22px;cursor:pointer;border:1px solid #bbb;border-radius:3px;padding:1px">
      </label>
      <button class="ctrl-btn" id="btn-prune-sptree" onclick="toggleSpPrune()" title="Toggle pruning the tree to only species present in the gene-tree data">Prune to data</button>
      <button class="ctrl-btn" id="btn-dl-anno" onclick="downloadAnnotations()">&#11015; Annotations TSV</button>
      <button class="ctrl-btn" id="btn-dl-newick" onclick="downloadNewick()">&#11015; Newick</button>
    </div>
    <div id="sp-tree-wrap"></div>
  </div>

</div>

<div id="sp-annot-popup" style="position:fixed;display:none;background:#fff;border:1px solid #bbb;border-radius:6px;padding:8px 10px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.18);z-index:300;min-width:150px">
  <div style="font-size:10px;color:#888;margin-bottom:6px;font-weight:600" id="sp-annot-popup-title"></div>
  <div id="sp-annot-popup-btns" style="display:flex;flex-direction:column;gap:4px"></div>
</div>
<div id="mini-sp-panel">
  <div class="msp-title">Species tree — click a named node to highlight that clade</div>
  <div id="mini-sp-svg-wrap"></div>
</div>
<div id="tooltip"></div>
<div id="collapse-choice-popup" style="position:fixed;display:none;background:#fff;border:1px solid #bbb;border-radius:6px;padding:6px 8px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.18);z-index:250;display:none;gap:6px;flex-direction:column">
  <div style="font-size:10px;color:#888;margin-bottom:2px">Collapse as:</div>
  <div style="display:flex;gap:6px;flex-wrap:wrap">
    <button id="ccp-triangle" style="padding:4px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer" title="Collapse to a filled triangle (proportional size)">&#x25BD; Triangle</button>
    <button id="ccp-circle" style="padding:4px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer" title="Collapse to a circle with leaf count">&#x25EF; Circle</button>
    <button id="ccp-compare" style="padding:4px 10px;font-size:11px;border:1px solid #4a90d9;border-radius:4px;background:#f0f6ff;color:#2c6090;cursor:pointer" title="Compare species coverage with another node">&#x2316; Compare</button>
    <button id="ccp-highlight" style="padding:4px 10px;font-size:11px;border:1px solid #e67e22;border-radius:4px;background:#fff8f0;color:#c0622a;cursor:pointer" title="Highlight subtree background with a colour">&#x25A0; Highlight</button>
  </div>
</div>
<div id="tri-action-popup" style="position:fixed;display:none;background:#fff;border:1px solid #bbb;border-radius:6px;padding:6px 8px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.18);z-index:260;flex-direction:column;gap:5px">
  <div style="font-size:10px;color:#888;margin-bottom:2px;font-weight:600" id="tri-popup-title"></div>
  <div style="display:flex;gap:5px;flex-wrap:wrap">
    <button id="tap-expand" style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x25B7; Expand</button>
    <button id="tap-rename" style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x270F; Rename</button>
    <button id="tap-color"  style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x1F3A8; Color</button>
    <button id="tap-compare" style="padding:3px 10px;font-size:11px;border:1px solid #4a90d9;border-radius:4px;background:#f0f6ff;color:#2c6090;cursor:pointer">&#x2316; Compare</button>
    <input id="tap-color-input" type="color" style="display:none">
  </div>
</div>
<div id="clade-hl-popup" style="position:fixed;display:none;background:#fff;border:1px solid #e67e22;border-radius:6px;padding:6px 8px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.2);z-index:270;flex-direction:column;gap:5px">
  <div style="font-size:10px;color:#888;margin-bottom:2px;font-weight:600">Clade highlight</div>
  <div style="display:flex;gap:5px">
    <button id="chl-rename" style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x270F; Label</button>
    <button id="chl-recolor" style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x1F3A8; Colour</button>
    <button id="chl-remove" style="padding:3px 10px;font-size:11px;border:1px solid #c0392b;border-radius:4px;background:#fff0ee;color:#c0392b;cursor:pointer">&#x2715; Remove</button>
  </div>
</div>
<div id="clado-action-popup" style="position:fixed;display:none;background:#fff;border:1px solid #bbb;border-radius:6px;padding:6px 8px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.18);z-index:260;flex-direction:column;gap:5px">
  <div style="font-size:10px;color:#888;margin-bottom:2px;font-weight:600" id="clado-popup-title"></div>
  <div style="display:flex;gap:5px;flex-wrap:wrap">
    <button id="cap-expand"  style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x25B7; Expand</button>
    <button id="cap-rename"  style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x270F; Rename</button>
    <button id="cap-color"   style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x1F3A8; Color</button>
    <input id="cap-color-input" type="color" style="display:none">
  </div>
</div>
<div id="collapsed-popup">
  <div class="cp-title">Collapsed node</div>
  <div class="cp-row"><span>Name:</span><input id="cp-name" type="text"></div>
  <div class="cp-genes-label" id="cp-genes-label"></div>
  <textarea id="cp-genes" readonly></textarea>
  <div class="cp-actions">
    <button class="cp-btn" onclick="cpRename()">Rename</button>
    <button class="cp-btn" onclick="cpCopy()">Copy genes</button>
    <button class="cp-btn" onclick="cpExpand()">Expand</button>
    <button class="cp-btn" onclick="cpCompare()" title="Compare species coverage with another node">&#x2316; Compare</button>
    <button class="cp-btn" onclick="cpClose()">Close</button>
  </div>
</div>

<!-- ── Species compare banner + result panel ── -->
<div id="sp-compare-banner" style="display:none;position:fixed;z-index:400;bottom:28px;left:50%;transform:translateX(-50%);background:#fffbea;border:1px solid #e8c840;border-radius:20px;padding:5px 16px;font-size:11px;color:#7a5c00;box-shadow:0 2px 8px rgba(0,0,0,.18);align-items:center;gap:10px;white-space:nowrap">
  <span>&#x2316; Click a second node to compare species&hellip;</span>
  <button onclick="exitCompareMode()" style="padding:1px 8px;font-size:10px;border:1px solid #c8a030;border-radius:10px;background:#fff;cursor:pointer;color:#7a5c00">Cancel</button>
</div>
<div id="sp-compare-panel" style="display:none;position:fixed;z-index:410;top:30%;left:50%;transform:translate(-50%,0);background:#fff;border:1px solid #bbb;border-radius:8px;padding:12px 16px;font-size:11px;box-shadow:0 6px 20px rgba(0,0,0,.22);min-width:240px;max-width:380px">
  <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px">
    <b style="font-size:12px;color:#333">Species overlap</b>
    <span onclick="document.getElementById('sp-compare-panel').style.display='none'" style="cursor:pointer;color:#aaa;font-size:16px;line-height:1;padding:0 2px" title="Close">&times;</span>
  </div>
  <div id="scp-legend" style="font-size:10px;color:#666;margin-bottom:8px;line-height:1.6"></div>
  <div id="scp-content"></div>
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
const NEWICK_RAW    = %%NEWICK_RAW%%;
const FAMILY_INFO   = %%FAMILY_INFO_JSON%%;
const HAVE_GENERAX  = %%HAVE_GENERAX_JSON%%;

// ═══════════════════════════════════════════════════════════════════════════════
// COLOUR SYSTEM
// ═══════════════════════════════════════════════════════════════════════════════
const palette = [
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
  "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
  "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5"
];
const SP_COLORS = {};
(function(){
  const n=SPECIES_ORDER.length;
  SPECIES_ORDER.forEach((sp,i)=>{ SP_COLORS[sp]=d3.interpolateTurbo(n>1?i/(n-1):0.5); });
  ALL_SPECIES.forEach(sp=>{ if(!SP_COLORS[sp]) SP_COLORS[sp]="#aaa"; });
})();
function spColor(sp) { return SP_COLORS[sp] || "#aaa"; }

// ── Stable IDs for SP_TREE_DATA nodes (allows editing original by reference) ──
const spNodeById = new Map();   // _spId → SP_TREE_DATA node reference
(function tagSpNodes(n, i={v:0}){
  n._spId = String(i.v++);
  spNodeById.set(n._spId, n);
  (n.children||[]).forEach(c=>tagSpNodes(c,i));
})(SP_TREE_DATA || {});

function treeToNewick(n){
  const name=(n.name||"").replace(/[(),:;]/g,"_");
  const dist=n.dist!=null?":"+n.dist:"";
  if(!n.children||!n.children.length) return name+dist;
  return "("+n.children.map(treeToNewick).join(",")+")"+(name||"")+dist;
}

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
let hlQueryColors = {};          // query string → custom hex color override
let ogHlSet       = null;        // null = off; Set<og_name> when active
let ogHlQueries   = [];          // committed OG query strings
let ogHlQueryColors = {};        // og query string → custom hex color override
let ogHlGroupIndex= new Map();   // og_name → group index
let showOGLabels   = true;        // toggle OG-node labels on expanded internals
let tipFontSize    = null;        // null = auto; number = user override (px)
let treeLinkWidth  = 1.3;        // branch stroke-width in screen px
let treeHeightMult = 1.0;        // vertical stretch multiplier for the tree
let treeWidthMult  = 1.0;        // horizontal stretch multiplier for the tree
let showGeneId    = true;        // tip label parts
let showOGName    = true;
let showRefOrtho  = true;
let showSupport        = false;   // show internal node support values
let hideNonHl          = false;   // hide non-highlighted tip labels when hlSet active
let focusCollapseAsTri = true;    // true=triangle (MRCA), false=circle for focus/hide-nonhl collapses
let hmFocusGids        = null;    // Set<gene_id> when navigating from heatmap cell; null otherwise
let collapsedFraction  = 1.0;     // fraction of proportional space for OG-collapsed tips (0.5–1.0)
let hmSplitSets   = [];      // array of Set<species> for heatmap row splitting
let hmSplitLabels = [];      // display label for each split group
let hmColFontSize = 9;       // column label font size (px)
let hmColRotation = 65;      // column label rotation angle (degrees)
let hmColorMode   = "zscore"; // "zscore" | "absolute"
let hmTextFilter  = "";       // free-text filter for heatmap columns
const spCollapsed = new Set();   // node _ids collapsed in species tree
let colTriFill = "#ffffff";      // fill colour for collapsed triangles
// pre-computed metadata per species: {genes, families, hgs}
const spMeta = (()=>{
  const m={};
  FAMILY_DATA.forEach(d=>{ for(const [sp,n] of Object.entries(d.species_counts)){ if(!m[sp]) m[sp]={genes:0,families:0,hgs:0}; m[sp].genes+=n; m[sp].families++; } });
  HG_DATA.forEach(d=>{ for(const sp of Object.keys(d.species_counts)){ if(!m[sp]) m[sp]={genes:0,families:0,hgs:0}; m[sp].hgs++; } });
  return m;
})();
let spTreeWidthPct = 50;         // % of pane width used by species tree SVG
let spPruneToData  = false;      // when true, drawSpeciesTree hides species absent from gene-tree data
const hlTagColors = ["#e74c3c","#3498db","#27ae60","#f39c12","#8e44ad","#16a085","#e67e22","#c0392b"];
function hlTagColor(gi){ const q=hlQueries[gi]; return (q&&hlQueryColors[q])||hlTagColors[gi%hlTagColors.length]; }
function ogHlTagColor(gi){ const q=ogHlQueries[gi]; return (q&&ogHlQueryColors[q])||hlTagColors[gi%hlTagColors.length]; }

// onPick  = called live on every input event (colour-drag preview)
// onFinal = called once on change event (picker closed); defaults to onPick
function openColorPicker(currentCol, onPick, onFinal){
  if(!onFinal) onFinal=onPick;
  const inp=document.createElement("input"); inp.type="color";
  inp.value=currentCol; inp.style.cssText="position:fixed;opacity:0;pointer-events:none";
  document.body.appendChild(inp);
  inp.addEventListener("input",()=>onPick(inp.value));
  inp.addEventListener("change",()=>{ onFinal(inp.value); document.body.removeChild(inp); });
  inp.addEventListener("blur",()=>{ if(inp.parentNode) document.body.removeChild(inp); });
  inp.click();
}

function leafColor(sp) {
  if (hlSet !== null) {
    if (!hlSet.has(sp)) return "#ccc";
    const gi = hlGroupIndex.get(sp);
    return gi !== undefined ? hlTagColor(gi)
                            : (colorMode === "species" ? spColor(sp) : (cladeSp2Color[sp] || "#ccc"));
  }
  return colorMode === "species" ? spColor(sp) : (cladeSp2Color[sp] || "#ccc");
}
function ogLeafColor(geneId, species) {
  if (hlSet !== null) {
    if (!hlSet.has(species||"")) return "#ccc";
    const gi = hlGroupIndex.get(species||"");
    return gi !== undefined ? hlTagColor(gi) : (ogLeaf2Color[geneId] || "#ccc");
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
let hmColSortSp    = null;   // sort heatmap columns by this species' count (descending)
let hmCustomOGs    = [];     // user-selected individual OG names for custom view
let hmCustomGroups = [];     // bulk selections: [{type:"class"|"family", key, label, ogs:[...]}]
let hmColOrderOverride = null; // drag-reorder: array of column ids overriding natural order
let hmOGIndex      = null;   // lazy: {og_name → {hgId,family,cls,total,species_counts}}
let hmGeneIndex    = null;   // lazy: {gene_id → og_name}
let hmSearchIndex  = null;   // lazy: [{label,type,cls,fam,hg_id}] for global search

function getEffectiveCustomOGs(){
  const all=new Set(hmCustomOGs);
  hmCustomGroups.forEach(g=>g.ogs.forEach(og=>all.add(og)));
  return [...all];
}

function getSpeciesPfx(geneId){
  const g=(geneId||"").split("|")[0].trim();
  const ui=g.indexOf("_"), di=g.indexOf(".");
  const idx=[ui,di].filter(x=>x>0).reduce((a,b)=>Math.min(a,b),Infinity);
  return idx===Infinity?g:g.slice(0,idx);
}

function switchTab(name) {
  document.querySelectorAll(".tab-btn").forEach(b => b.classList.toggle("active", b.dataset.tab===name));
  document.querySelectorAll(".tab-pane").forEach(p => p.classList.toggle("active", p.id==="pane-"+name));
  const tc  = document.getElementById("tree-count");
  const hb  = document.getElementById("hm-back");
  const cr  = document.getElementById("hm-breadcrumb");
  const pfx = document.getElementById("prefixSelect").parentElement; // the <label>
  if (name==="trees") {
    tc.style.display = "inline"; pfx.style.display = "none";
    hb.style.display = "none"; cr.textContent = "";
    if (!currentIndex && TREE_INDEX.length) { renderSidebar(""); selectTree(TREE_INDEX[0]); }
  } else if (name==="sptree") {
    tc.style.display = "none"; pfx.style.display = "none";
    drawSpeciesTree();
  } else if (name==="families") {
    tc.style.display = "none"; pfx.style.display = "none";
    hb.style.display = "none"; cr.textContent = "";
    drawFamilyTable();
  } else {
    tc.style.display = "none"; pfx.style.display = "";
    drawHeatmap(); // drawCladogram() is called at the end of drawHeatmap()
  }
}

// ═══════════════════════════════════════════════════════════════════════════════
// FAMILIES TABLE
// ═══════════════════════════════════════════════════════════════════════════════

let _famSortCol = "family", _famSortAsc = true;
let _famFilter  = "";
let _famRendered = false;

const _clsCatColors = {
  tfs:"#3498db", chr:"#9b59b6", sig:"#e67e22", neu:"#27ae60",
  rbp:"#e74c3c", ion:"#1abc9c", epi:"#f39c12", kin:"#2980b9"
};
function _clsBadge(cls){
  if(!cls) return '<span style="color:#bbb">—</span>';
  const c=_clsCatColors[cls]||"#7f8c8d";
  return `<span style="background:${c};color:#fff;font-size:10px;padding:1px 6px;border-radius:8px;white-space:nowrap;cursor:pointer"
    onclick="famGoToClass('${cls.replace(/'/g,"\\'")}')"\
    title="Show all ${cls} families in Counts tab">${cls}</span>`;
}

function famGoToClass(cls){
  switchTab("heatmap");
  hmViewMode="family"; hmActiveClass=cls; hmActiveFamily=null; hmActiveHG=null; hmActiveHGRec=null;
  drawHeatmap();
  const back=document.getElementById("hm-back");
  const cr=document.getElementById("hm-breadcrumb");
  if(back) back.style.display="inline";
  if(cr) cr.textContent=cls;
}

function drawFamilyTable(){
  const tbody = document.getElementById("fam-tbody");
  const countEl = document.getElementById("fam-count");
  const q = _famFilter.toLowerCase();

  let rows = FAMILY_INFO.filter(d=>{
    if(!q) return true;
    return d.family.toLowerCase().includes(q)
        || (d.pfam||[]).some(p=>p.toLowerCase().includes(q))
        || (d.category||"").toLowerCase().includes(q)
        || (d.cls||"").toLowerCase().includes(q);
  });

  // sort
  rows = rows.slice().sort((a,b)=>{
    let va=a[_famSortCol], vb=b[_famSortCol];
    if(_famSortCol==="pfam"){ va=(va||[]).join(","); vb=(vb||[]).join(","); }
    if(typeof va==="number") return _famSortAsc ? va-vb : vb-va;
    va=String(va||"").toLowerCase(); vb=String(vb||"").toLowerCase();
    return _famSortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
  });

  countEl.textContent = rows.length + " / " + FAMILY_INFO.length + " families";

  tbody.innerHTML = rows.map((d,i)=>{
    const pfamLinks = (d.pfam||[]).map(p=>
      `<a href="https://www.ebi.ac.uk/interpro/search/text/${encodeURIComponent(p)}/"
          target="_blank" rel="noopener"
          style="color:#2980b9;text-decoration:none;margin-right:6px;white-space:nowrap"
          title="Search InterPro for ${p}">${p}</a>`
    ).join("");
    const bg = i%2===0?"#fff":"#f9f9f9";
    // colour-code tree-coverage fractions
    function treeFrac(n, tot){
      if(!tot) return `<span style="color:#bbb">—</span>`;
      const frac = n/tot;
      const col = frac>=1 ? "#27ae60" : frac>0 ? "#e67e22" : "#c0392b";
      return `<span style="color:${col};font-weight:${frac<1?"600":"normal"}">${n}</span><span style="color:#aaa">/${tot}</span>`;
    }
    return `<tr style="background:${bg};border-bottom:1px solid #eee">
      <td style="padding:5px 8px;font-weight:600;cursor:pointer;color:#2c3e50;white-space:nowrap"
          onclick="famGoToFamily('${d.family.replace(/'/g,"\\'")}')"
          title="View in Counts tab">${d.family}</td>
      <td style="padding:5px 8px">${pfamLinks||'<span style="color:#bbb">—</span>'}</td>
      <td style="padding:5px 8px">${_clsBadge(d.cls)}<span style="color:#888;font-size:11px;margin-left:4px">${d.category||""}</span></td>
      <td style="padding:5px 8px;text-align:right;color:#555">${d.n_hgs||0}</td>
      <td style="padding:5px 8px;text-align:right">${treeFrac(d.n_trees||0, d.n_hgs||0)}</td>
      <td class="fam-td-generax" style="padding:5px 8px;text-align:right">${treeFrac(d.n_generax||0, d.n_hgs||0)}</td>
      <td style="padding:5px 8px;text-align:right;color:#555">${d.total||0}</td>
      <td style="padding:5px 8px;text-align:right;color:#555">${d.n_species||0}</td>
    </tr>`;
  }).join("");

  // wire sort headers (only once)
  if(!_famRendered){
    // hide GeneRax column when no GeneRax trees are available
    if(!HAVE_GENERAX){
      document.querySelectorAll(".fam-th-generax,.fam-td-generax").forEach(el=>el.style.display="none");
    }
    document.querySelectorAll(".fam-th").forEach(th=>{
      th.addEventListener("click",()=>{
        const col=th.dataset.col;
        if(_famSortCol===col) _famSortAsc=!_famSortAsc;
        else { _famSortCol=col; _famSortAsc=true; }
        drawFamilyTable();
      });
    });
    _famRendered=true;
  }
}

function filterFamilyTable(val){
  _famFilter=val;
  drawFamilyTable();
}

function famGoToFamily(family){
  // Navigate to the heatmap Counts tab and drill into this family
  switchTab("heatmap");
  // After heatmap draws, find matching family record and navigate to HG view
  const fRec = FAMILY_DATA.find(d=>d.family===family);
  if(!fRec) return;
  hmViewMode="hg"; hmActiveFamily=family; hmActiveClass=fRec.class||null;
  hmActiveHG=null; hmActiveHGRec=null;
  drawHeatmap();
  const back=document.getElementById("hm-back");
  const cr=document.getElementById("hm-breadcrumb");
  if(back) back.style.display="inline";
  if(cr) cr.textContent=(hmActiveClass?hmActiveClass+" › ":"")+family;
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

// ── Collapse-style choice popup ────────────────────────────────────────────
(function(){
  const pop=document.getElementById("collapse-choice-popup");
  let _ccpNode=null, _ccpTimer=null;
  function hide(){ pop.style.display="none"; _ccpNode=null; }
  function cancelHide(){ clearTimeout(_ccpTimer); }
  function schedHide(){ _ccpTimer=setTimeout(hide,200); }
  pop.addEventListener("mouseenter",cancelHide);
  pop.addEventListener("mouseleave",schedHide);
  document.getElementById("ccp-triangle").addEventListener("click",()=>{
    const node=_ccpNode; if(!node) return; hide();
    node._children=node.children; node.children=null;
    node._isOgCol=true; renderTree(true);
  });
  document.getElementById("ccp-circle").addEventListener("click",()=>{
    const node=_ccpNode; if(!node) return; hide();
    node._children=node.children; node.children=null;
    node._isOgCol=false; renderTree(true);
  });
  document.getElementById("ccp-compare").addEventListener("click",()=>{
    const node=_ccpNode; if(!node) return; hide();
    enterCompareMode(node);
  });
  document.getElementById("ccp-highlight").addEventListener("click",()=>{
    const node=_ccpNode; if(!node) return; hide();
    const existing=cladeHighlights.get(node._uid)||{};
    openColorPicker(
      existing.color||"#ffe066",
      // live: preview the colour while dragging (keep existing label)
      c=>{ cladeHighlights.set(node._uid,{color:c, label:existing.label||""}); renderTree(false); },
      // final: picker closed — now prompt for label once
      c=>{
        const labelRaw=window.prompt("Clade label (leave blank for none):", existing.label||"");
        cladeHighlights.set(node._uid,{color:c, label:(labelRaw===null ? existing.label||"" : labelRaw.trim())});
        renderTree(false);
      }
    );
  });
  window.showCollapseChoicePopup=function(event,d){
    _ccpNode=d; cancelHide();
    pop.style.display="flex";
    const x=Math.min(event.clientX+8, window.innerWidth-pop.offsetWidth-8);
    const y=Math.min(event.clientY+8, window.innerHeight-pop.offsetHeight-8);
    pop.style.left=Math.max(4,x)+"px"; pop.style.top=Math.max(4,y)+"px";
  };
  // close when clicking elsewhere
  document.addEventListener("click",(e)=>{ if(!pop.contains(e.target)) hide(); });
})();

// ── Clade highlight right-click popup ───────────────────────────────────────
(function(){
  const pop=document.getElementById("clade-hl-popup");
  let _uid=null;
  function hide(){ pop.style.display="none"; _uid=null; }
  window.showCladeHlPopup=function(event,uid){
    _uid=uid;
    pop.style.display="flex";
    const x=Math.min(event.clientX+8, window.innerWidth-pop.offsetWidth-8);
    const y=Math.min(event.clientY+8, window.innerHeight-pop.offsetHeight-8);
    pop.style.left=Math.max(4,x)+"px"; pop.style.top=Math.max(4,y)+"px";
    event.preventDefault(); event.stopPropagation();
  };
  document.getElementById("chl-rename").addEventListener("click",()=>{
    const rec=cladeHighlights.get(_uid); if(!rec) return; hide();
    const newLabel=window.prompt("Clade label:", rec.label||"");
    if(newLabel!==null){ rec.label=newLabel.trim(); renderTree(false); }
  });
  document.getElementById("chl-recolor").addEventListener("click",()=>{
    const rec=cladeHighlights.get(_uid); if(!rec) return; hide();
    openColorPicker(rec.color||"#ffe066",c=>{ rec.color=c; renderTree(false); });
  });
  document.getElementById("chl-remove").addEventListener("click",()=>{
    const uid=_uid; hide();
    cladeHighlights.delete(uid); renderTree(false);
  });
  document.addEventListener("click",(e)=>{ if(!pop.contains(e.target)) hide(); });
})();

// ── Cladogram collapsed-clade action popup ──────────────────────────────────
(function(){
  const pop=document.getElementById("clado-action-popup");
  const title=document.getElementById("clado-popup-title");
  let _key=null, _drawFn=null;
  function hide(){ pop.style.display="none"; _key=null; _drawFn=null; }
  document.getElementById("cap-expand").addEventListener("click",()=>{
    if(!_key) return; const k=_key, fn=_drawFn; hide();
    cladoCollapsed.delete(k); fn();
    if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
  });
  document.getElementById("cap-rename").addEventListener("click",()=>{
    if(!_key) return; const k=_key, name=title.textContent, fn=_drawFn; hide();
    const v=window.prompt("Rename clade:",name); if(v===null) return;
    cladoNames.set(k,v.trim()||name); fn();
    if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
  });
  const colorInput=document.getElementById("cap-color-input");
  document.getElementById("cap-color").addEventListener("click",()=>{ colorInput.click(); });
  colorInput.addEventListener("input",()=>{
    if(!_key) return; const k=_key, fn=_drawFn;
    cladoColors.set(k,colorInput.value); fn();
    if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
  });
  document.addEventListener("click",(e)=>{ if(!pop.contains(e.target)) hide(); });
  window.showCladoActionPopup=function(event,key,displayName,redrawFn){
    _key=key; _drawFn=redrawFn;
    title.textContent=displayName;
    colorInput.value=cladoColors.get(key)||"#e0e0e0";
    pop.style.display="flex";
    const x=Math.min(event.clientX+8,window.innerWidth-pop.offsetWidth-8);
    const y=Math.min(event.clientY+8,window.innerHeight-pop.offsetHeight-8);
    pop.style.left=Math.max(4,x)+"px"; pop.style.top=Math.max(4,y)+"px";
  };
})();

// ── Gene-tree collapsed-triangle action popup ───────────────────────────────
const nodeTriColors=new Map();    // _uid → custom fill override for gene-tree triangles
const cladeHighlights=new Map(); // _uid → {color, label} for subtree background highlight
let cladeHlAlpha = 0.22;         // global opacity for all clade highlight rectangles
(function(){
  const pop=document.getElementById("tri-action-popup");
  const title=document.getElementById("tri-popup-title");
  let _nd=null;
  function hide(){ pop.style.display="none"; _nd=null; }
  document.getElementById("tap-expand").addEventListener("click",()=>{
    const node=_nd; if(!node) return; hide();
    node.children=node._children; node._children=null; node._isOgCol=false; renderTree(true);
  });
  document.getElementById("tap-rename").addEventListener("click",()=>{
    const node=_nd; if(!node) return; hide();
    const cur=customNodeNames[node._uid]||collapsedLabel(node);
    const v=window.prompt("Rename collapsed node:",cur); if(v===null) return;
    customNodeNames[node._uid]=v.trim()||cur; renderTree(true);
  });
  const colorInput=document.getElementById("tap-color-input");
  document.getElementById("tap-color").addEventListener("click",()=>{ colorInput.click(); });
  colorInput.addEventListener("input",()=>{
    if(!_nd) return;
    nodeTriColors.set(_nd._uid,colorInput.value);
    // update the polygon inline without full re-render
    gMain&&gMain.selectAll(".col-tri").filter(d=>d&&d._uid===_nd._uid).style("fill",colorInput.value);
  });
  document.getElementById("tap-compare").addEventListener("click",()=>{
    const node=_nd; if(!node) return; hide();
    enterCompareMode(node);
  });
  document.addEventListener("click",(e)=>{ if(!pop.contains(e.target)) hide(); });
  window.showTriActionPopup=function(event,d){
    _nd=d; event.stopPropagation();
    title.textContent=collapsedLabel(d);
    colorInput.value=nodeTriColors.get(d._uid)||colTriFill||"#ffffff";
    pop.style.display="flex";
    const x=Math.min(event.clientX+8,window.innerWidth-pop.offsetWidth-8);
    const y=Math.min(event.clientY+8,window.innerHeight-pop.offsetHeight-8);
    pop.style.left=Math.max(4,x)+"px"; pop.style.top=Math.max(4,y)+"px";
  };
})();

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

function cpCompare() {
  const d = cpActiveNode;
  if (!d) return;
  cpClose();
  enterCompareMode(d);
}

// ═══════════════════════════════════════════════════════════════════════════════
// HEATMAP VIEW
// ═══════════════════════════════════════════════════════════════════════════════
const TOP_MARGIN = 110;
const HM_TOP     = TOP_MARGIN;   // minimum first-row y-coord (may grow with label height)
const HM_BAR_H   = 55;           // height of column-sum bar chart below the heatmap
let   hmTopActual = HM_TOP;      // actual top margin, synced between heatmap and cladogram

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
// sort-by-species selector — populate from ALL_SPECIES
(function(){
  const sel=document.getElementById("hm-col-sort-sp");
  ALL_SPECIES.forEach(sp=>{ const o=document.createElement("option"); o.value=sp; o.textContent=sp; sel.appendChild(o); });
})();

const cladoCollapsed  = new Set();   // stable node keys of collapsed cladogram nodes
const cladoNames      = new Map();   // key → custom display name (user-editable)
const cladoColors     = new Map();   // key → custom fill color override
// map: clade label → sorted array of species (updated by drawCladogram, read by drawHeatmap)
const hmCollapsedBands = [];  // [{species:[...], label}] — reset on each drawCladogram call

function drawCladogram() {
  hmCollapsedBands.length = 0;  // reset collapsed-clade shading bands
  const tp = document.getElementById("tree-panel");
  tp.innerHTML = "";
  if (!SP_TREE_DATA || !SP_TREE_DATA.children || !speciesOrder.length) return;

  // Dynamic vertical offset so cladogram rows align with heatmap rows.
  // tree-panel and heatmap-panel may start at different viewport y-positions
  // (hm-col-strip sits above heatmap-panel; prefix selector sits above tree-panel).
  const hp = document.getElementById("heatmap-panel");
  const topPad = (hp)
    ? Math.max(0, Math.round(hp.getBoundingClientRect().top - tp.getBoundingClientRect().top))
    : 0;

  const W = 270, H = speciesOrder.length*14+hmTopActual+topPad+40;
  const svg = d3.select(tp).append("svg").attr("width",W).attr("height",H);
  const leafY = {}; speciesOrder.forEach((s,i)=>{ leafY[s]=hmTopActual+topPad+i*14+7; });

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
  // stable key for a node: sorted leaf names (survives re-renders)
  function nodeKey(n){ return flat(n).filter(x=>!x.children).map(x=>x.name).sort().join(","); }
  // y-range of all leaves in a subtree
  function yRange(n){ const ys=flat(n).filter(x=>!x.children).map(x=>leafY[x.name]||0); return [d3.min(ys),d3.max(ys)]; }

  assignY(tree); assignX(tree);
  const maxD=d3.max(flat(tree),d=>d._d)||1;
  flat(tree).forEach(n=>{ if(!n.children) n._d=maxD; }); // align all tips
  const sx=d=>10+(d/maxD)*150;

  // ── branches (skip into collapsed subtrees) ──────────────────────────────
  function drawB(n){
    if(!n.children) return;
    if(cladoCollapsed.has(nodeKey(n))) return;
    const ys=n.children.map(c=>c._y);
    svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(n._d)).attr("y1",d3.min(ys)).attr("y2",d3.max(ys)).attr("stroke","#aaa");
    n.children.forEach(c=>{
      svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(c._d)).attr("y1",c._y).attr("y2",c._y).attr("stroke","#aaa");
      drawB(c);
    });
  }
  drawB(tree);

  // ── collapsed triangles ───────────────────────────────────────────────────
  function drawCollapsedTriangles(n){
    if(!n.children) return;
    const key=nodeKey(n);
    if(cladoCollapsed.has(key)){
      const [y0,y1]=yRange(n);
      const displayName=cladoNames.get(key)||(n.name||"clade");
      // collect leaf species for heatmap shading
      function getLeafSp(nd){ return nd.children?nd.children.flatMap(getLeafSp):[nd.name]; }
      hmCollapsedBands.push({species:getLeafSp(n),label:displayName});
      const triColor=cladoColors.get(key)||"#ffffff";
      const redraw=()=>{ drawCladogram(); };
      svg.append("polygon")
        .attr("points",`${sx(n._d)},${n._y} ${sx(maxD)},${y0} ${sx(maxD)},${y1}`)
        .attr("fill",triColor).attr("stroke","#222").attr("stroke-width",1.2)
        .style("cursor","pointer")
        .on("mouseover",ev=>showTip(ev,displayName+'<div style="font-size:9px;color:#aaa;margin-top:2px">click for options</div>'))
        .on("mousemove",moveTip).on("mouseout",hideTip)
        .on("click",(ev)=>{ ev.stopPropagation(); hideTip(); showCladoActionPopup(ev,key,displayName,redraw); });
      // label — normal font, clickable
      svg.append("text")
        .attr("x",sx(n._d)+6).attr("y",n._y).attr("dy","0.35em")
        .attr("font-size",9).attr("fill","#333")
        .style("cursor","pointer").text(displayName)
        .on("click",(ev)=>{ ev.stopPropagation(); hideTip(); showCladoActionPopup(ev,key,displayName,redraw); })
        .on("mouseover",ev=>showTip(ev,displayName+'<div style="font-size:9px;color:#aaa;margin-top:2px">click for options</div>'))
        .on("mousemove",moveTip).on("mouseout",hideTip);
      return;
    }
    n.children.forEach(drawCollapsedTriangles);
  }
  drawCollapsedTriangles(tree);

  // ── interactive internal nodes (click=collapse, shift+click=heatmap split) ─
  function drawCladoNodes(n){
    if(!n.children) return;
    const key=nodeKey(n);
    if(cladoCollapsed.has(key)) return;
    const nodeLabel = n.name || "(unnamed)";
    const cleavesSet = new Set(flat(n).filter(x=>!x.children).map(x=>x.name));
    const splitIdx = hmSplitSets.findIndex(s=>s.size===cleavesSet.size&&[...cleavesSet].every(v=>s.has(v)));
    const isSplit = splitIdx >= 0;
    const nodeCol = isSplit ? hmSplitLineColors[splitIdx%hmSplitLineColors.length] : "#aaa";
    const tip = nodeLabel+'<div style="font-size:9px;color:#aaa;margin-top:2px">click to collapse &nbsp;·&nbsp; shift+click to split heatmap</div>';
    svg.append("circle").attr("cx",sx(n._d)).attr("cy",n._y).attr("r",4)
      .attr("fill", isSplit ? hmSplitColors[splitIdx%hmSplitColors.length] : "#fff")
      .attr("stroke",nodeCol).attr("stroke-width", isSplit?2:1)
      .style("cursor","pointer")
      .on("mouseover",ev=>showTip(ev,tip))
      .on("mousemove",moveTip).on("mouseout",hideTip)
      .on("click",(ev)=>{
        hideTip();
        if(ev.shiftKey){
          ev.stopPropagation();
          if(isSplit){ hmSplitSets.splice(splitIdx,1); hmSplitLabels.splice(splitIdx,1); }
          else { hmSplitSets.push(cleavesSet); hmSplitLabels.push(nodeLabel); }
          updateHmSplitBar(); drawHeatmap(); // drawCladogram called at end of drawHeatmap
        } else {
          cladoCollapsed.add(key); drawCladogram();
          if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
        }
      });
    n.children.forEach(drawCladoNodes);
  }
  drawCladoNodes(tree);
}

// ── helpers for heatmap row splitting ─────────────────────────────────────
function spNodeLeaves(n) {
  if (!n.children) return n.name ? [n.name] : [];
  return n.children.flatMap(spNodeLeaves);
}
const hmSplitColors     = ["rgba(231,76,60,0.10)","rgba(52,152,219,0.10)","rgba(39,174,96,0.10)","rgba(243,156,18,0.10)"];
const hmSplitLineColors = ["#e74c3c","#3498db","#27ae60","#f39c12"];

function updateHmSplitBar() {
  const bar = document.getElementById("hm-split-bar");
  if (!hmSplitSets.length) { bar.style.display = "none"; return; }
  bar.style.display = "flex";
  document.getElementById("hm-split-tags").innerHTML = hmSplitLabels.map((lbl,i) =>
    `<span style="display:inline-flex;align-items:center;padding:1px 8px;border-radius:10px;background:${hmSplitLineColors[i%hmSplitLineColors.length]};color:#fff;font-size:10px">${lbl}</span>`
  ).join(" ");
}
function clearHmSplits() {
  hmSplitSets = []; hmSplitLabels = [];
  updateHmSplitBar();
  drawHeatmap();
}

function drawSpeciesTree() {
  const wrap = document.getElementById("sp-tree-wrap");
  wrap.innerHTML = "";
  if (!SP_TREE_DATA || !SP_TREE_DATA.children || !SPECIES_ORDER.length) {
    wrap.innerHTML = '<p style="padding:20px;color:#888">No species tree provided (pass --species_tree).</p>';
    return;
  }

  const rowH = 22, topM = 16, leftM = 16, rightM = 260;
  const W = Math.max(Math.floor((wrap.clientWidth || 900) * spTreeWidthPct / 100), 300);
  const treeW = W - leftM - rightM;
  const inPhylo = new Set(ALL_SPECIES);
  const tipColor = n => inPhylo.has(n.name) ? spColor(n.name) : "#bbb";

  // ── helpers ──────────────────────────────────────────────────────────
  function clone(n){ return JSON.parse(JSON.stringify(n)); }
  function prune(n){
    if(!n.children) return (spPruneToData?(inPhylo.has(n.name)||dataSpecies.has(n.name)):SPECIES_ORDER.includes(n.name))?n:null;
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
        .attr("fill",colTriFill).attr("stroke","#222").attr("stroke-width",1)
        .style("cursor","pointer")
        .on("click",()=>{ spCollapsed.delete(n._id); drawSpeciesTree(); })
        .on("mouseover",ev=>showTip(ev,'<div class="tt-name">'+n.name+'</div><div style="font-size:9px;color:#aaa">click to expand</div>'))
        .on("mousemove",moveTip).on("mouseout",hideTip);
      // small circle at apex for easy uncollapse
      svg.append("circle")
        .attr("cx",x0).attr("cy",yM).attr("r",5)
        .attr("fill","#f5f5f5").attr("stroke","#555").attr("stroke-width",1.2)
        .style("cursor","pointer")
        .on("click",()=>{ spCollapsed.delete(n._id); drawSpeciesTree(); })
        .on("mouseover",ev=>showTip(ev,(n.name?'<div class="tt-name">'+n.name+'</div>':'')+'<div style="font-size:9px;color:#aaa">click to expand</div>'))
        .on("mousemove",moveTip).on("mouseout",hideTip);
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
    // clickable node dot — only non-root internal nodes
    if(n._id!==tree._id){
      const nodeLabel=n.name||"(unnamed)";
      // check if this clade is already a split group
      const cleavesSet=new Set(spNodeLeaves(n));
      const splitIdx=hmSplitSets.findIndex(s=>s.size===cleavesSet.size&&[...cleavesSet].every(v=>s.has(v)));
      const isSplit=splitIdx>=0;
      const nodeCol=isSplit?hmSplitLineColors[splitIdx%hmSplitLineColors.length]:"#999";
      const tip=nodeLabel+'<div style="font-size:9px;color:#aaa;margin-top:3px">click to collapse &nbsp;·&nbsp; shift+click to split heatmap</div>';
      // always-visible node name label (italic, small, to the right of the dot)
      // double-click to edit the name; edits propagate to SP_TREE_DATA
      if(n.name)
        svg.append("text").attr("x",sx(n._d)+7).attr("y",n._y).attr("dy","0.35em")
          .attr("font-size",9).attr("fill",nodeCol).attr("font-style","italic")
          .style("cursor","text").style("user-select","none").text(n.name)
          .on("dblclick",(ev)=>{
            ev.stopPropagation();
            const orig=spNodeById.get(n._spId);
            if(!orig) return;
            const v=window.prompt("Edit node label:",orig.name||"");
            if(v===null) return;
            orig.name=v.trim();
            drawSpeciesTree(); drawCladogram();
            if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
          });
      svg.append("circle").attr("cx",sx(n._d)).attr("cy",n._y).attr("r",5)
        .attr("fill",isSplit?hmSplitColors[splitIdx%hmSplitColors.length]:"#fff")
        .attr("stroke",nodeCol).attr("stroke-width",isSplit?2:1.2)
        .style("cursor","pointer")
        .on("mouseover",ev=>showTip(ev,tip))
        .on("mousemove",moveTip)
        .on("mouseout",hideTip)
        .on("click",(ev)=>{
          if(ev.shiftKey){
            ev.stopPropagation();
            if(isSplit){ hmSplitSets.splice(splitIdx,1); hmSplitLabels.splice(splitIdx,1); }
            else{ hmSplitSets.push(cleavesSet); hmSplitLabels.push(nodeLabel); }
            updateHmSplitBar();
            hideTip(); drawSpeciesTree(); drawHeatmap();
          } else {
            spCollapsed.add(n._id); hideTip(); drawSpeciesTree();
          }
        });
    }
    n.children.forEach(drawInternals);
  }
  drawInternals(tree);

  // ── tip legend (black = in gene trees, grey = absent) ─────────────────
  const legX=W-rightM+8, legY=topM;
  svg.append("circle").attr("cx",legX).attr("cy",legY).attr("r",4).attr("fill","#111");
  svg.append("text").attr("x",legX+8).attr("y",legY).attr("dy","0.35em").attr("font-size",10).attr("fill","#444").text("in gene trees");
  svg.append("circle").attr("cx",legX).attr("cy",legY+16).attr("r",4).attr("fill","#bbb");
  svg.append("text").attr("x",legX+8).attr("y",legY+16).attr("dy","0.35em").attr("font-size",10).attr("fill","#888").text("absent");

  // ── leaf tips ─────────────────────────────────────────────────────────
  function drawLeaves(n){
    if(n._isCol) return;
    if(!n.children){
      const tc=tipColor(n);
      const inData=inPhylo.has(n.name);
      svg.append("circle").attr("cx",sx(n._d)).attr("cy",n._y).attr("r",5)
        .attr("fill",tc).attr("stroke","#fff").attr("stroke-width",0.8)
        .style("cursor",inData?"pointer":"default")
        .on("mouseover",(ev)=>{ if(inData) showTip(ev,'<div style="font-size:10px;color:#aaa">click to change color</div>'); })
        .on("mousemove",moveTip).on("mouseout",hideTip)
        .on("click",(ev)=>{
          if(!inData) return;
          ev.stopPropagation();
          openColorPicker(spColor(n.name),c=>{
            SP_COLORS[n.name]=c;
            drawSpeciesTree();
            if(rootNode) renderTree(false);
          });
        });
      svg.append("text").attr("x",sx(n._d)+10).attr("y",n._y).attr("dy","0.35em")
        .attr("font-size",12).attr("fill",inData?"#222":"#bbb").attr("font-family","monospace").text(n.name)
        .style("cursor","pointer")
        .on("mouseover",(ev)=>{
          const m=spMeta[n.name]||{genes:0,families:0,hgs:0};
          showTip(ev,'<div class="tt-name" style="color:'+tc+'">'+n.name+'</div>'
            +'<div class="tt-row"><span>Annotated genes</span><strong>'+m.genes+'</strong></div>'
            +'<div class="tt-row"><span>Families</span><strong>'+m.families+'</strong></div>'
            +'<div class="tt-row"><span>Homology groups</span><strong>'+m.hgs+'</strong></div>'
            +'<div style="font-size:9px;color:#aaa;margin-top:3px">click to export annotations</div>');
        }).on("mousemove",moveTip).on("mouseout",hideTip)
        .on("click",(ev)=>{ hideTip(); showSpAnnotPopup(ev, n.name); });
      return;
    }
    n.children.forEach(drawLeaves);
  }
  drawLeaves(tree);
}

function toggleSpPrune(){
  spPruneToData=!spPruneToData;
  const btn=document.getElementById("btn-prune-sptree");
  if(btn){ btn.style.background=spPruneToData?"#d0e8ff":""; btn.style.fontWeight=spPruneToData?"600":""; }
  drawSpeciesTree();
}

// ── helpers for sort-by-species buttons ───────────────────────────────────
function hmSortActive(sps){
  if(!hmColSortSp) return false;
  const cur=Array.isArray(hmColSortSp)?hmColSortSp:[hmColSortSp];
  const tgt=Array.isArray(sps)?sps:[sps];
  return cur.length===tgt.length&&tgt.every(s=>cur.includes(s));
}
function hmSetSort(sps){
  const tgt=Array.isArray(sps)?sps:[sps];
  hmColSortSp=hmSortActive(tgt)?null:(tgt.length===1?tgt[0]:tgt);
  hmColOrderOverride=null;
  const dd=document.getElementById("hm-col-sort-sp");
  if(dd) dd.value=(typeof hmColSortSp==="string"?hmColSortSp:"");
  drawHeatmap();
}

function drawHeatmap() {
  const prefix = document.getElementById("prefixSelect").value;
  const back   = document.getElementById("hm-back");
  const crumb  = document.getElementById("hm-breadcrumb");
  const expandBtn = document.getElementById("hm-expand-og");
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
    back.style.display="none"; crumb.textContent=""; expandBtn.style.display="none";
    colLabel=d=>d.id;
    clickHandler=(_ev,d)=>{ hmViewMode="family"; hmActiveClass=d.id; drawHeatmap(); };

  } else if (hmViewMode==="family") {
    data=FAMILY_DATA.filter(d=>(d.class||d.pref)===hmActiveClass)
                    .filter(d=>prefix==="all"||d.pref===prefix)
                    .sort((a,b)=>b.total-a.total);
    back.style.display="inline"; crumb.textContent=hmActiveClass; expandBtn.style.display="none";
    colLabel=d=>d.family;
    clickHandler=(_ev,d)=>{ hmViewMode="hg"; hmActiveFamily=d.family; drawHeatmap(); };

  } else if (hmViewMode==="hg") {
    data=HG_DATA.filter(d=>d.family===hmActiveFamily)
                .filter(d=>prefix==="all"||d.pref===prefix)
                .sort((a,b)=>b.total-a.total);
    back.style.display="inline"; crumb.textContent=hmActiveClass+" \u203a "+hmActiveFamily;
    expandBtn.style.display="inline"; expandBtn.dataset.hgData=JSON.stringify(data.map(d=>d.id));
    colLabel=d=>d.hg||d.id;
    clickHandler=(ev,d,sp)=>{
      if(sp===undefined){
        if(ev.shiftKey){ hmExpandHGToOGs([d.id]); return; }
        hmViewMode="og"; hmActiveHG=d.id; hmActiveHGRec=d; drawHeatmap(); return;
      }
      // cell click: open gene tree for this HG with species highlighted
      const tRec=TREE_INDEX.find(t=>t.id===d.id)||TREE_INDEX.find(t=>t.family===d.family&&t.hg===d.hg);
      if(!tRec) return;
      switchTab("trees"); selectTree(tRec); renderSidebar("");
      if(d.species_counts[sp]){ hlQueries=[sp]; rebuildHlSet(); setTimeout(focusHighlighted,60); }
    };

  } else if(hmViewMode==="og") { // og
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
    back.style.display="inline"; crumb.textContent=hmActiveClass+" \u203a "+hmActiveFamily+" \u203a "+(r&&r.hg||hmActiveHG); expandBtn.style.display="none";
    colLabel=d=>d.id+" ["+d.total+"]";
    // third arg 'sp': species clicked (from cell), undefined for column-header click
    clickHandler=(_ev,_d,sp)=>{
      if(!treeRec) return;
      // resolve gene IDs directly from detail.ogs before selectTree resets state
      const det=loadDetail(treeRec.id);
      const ogGids=(det&&det.ogs&&det.ogs[_d.id])||[];
      const targetGids=sp?ogGids.filter(g=>getSpeciesPfx(g)===sp):ogGids;
      switchTab("trees"); selectTree(treeRec); renderSidebar("");
      if(targetGids.length){
        hmFocusGids=new Set(targetGids);
        hideNonHl=true;
        const chk=document.getElementById("chk-hide-nonhl"); if(chk) chk.checked=true;
      }
      setTimeout(focusHighlighted,80);
    };

  } else { // custom
    buildOGIndex();
    const effOGs=getEffectiveCustomOGs();
    if(!effOGs.length){
      document.getElementById("heatmap-panel").innerHTML=
        '<p style="padding:30px 20px;color:#aaa;font-size:13px">'+
        'Use the search box above to add orthogroups to the custom view.<br>'+
        'You can search by OG name, gene name, class or homology group.</p>';
      return;
    }
    data=effOGs.map(og=>{
      const r=hmOGIndex[og]||{};
      return {id:og,species_counts:r.species_counts||{},total:r.total||0,family:r.family||"",cls:r.cls||""};
    });
    back.style.display="none"; expandBtn.style.display="none";
    crumb.textContent="Custom selection ("+effOGs.length+" OG"+(effOGs.length!==1?"s":"")+")";
    colLabel=d=>d.id;
    clickHandler=(_ev,d,sp)=>{
      const tRec=hmOGIndex[d.id]?TREE_INDEX.find(t=>t.id===hmOGIndex[d.id].hgId):null;
      if(!tRec) return;
      const det=loadDetail(tRec.id);
      const ogGids=(det&&det.ogs&&det.ogs[d.id])||[];
      const targetGids=sp?ogGids.filter(g=>getSpeciesPfx(g)===sp):ogGids;
      switchTab("trees"); selectTree(tRec); renderSidebar("");
      if(targetGids.length){
        hmFocusGids=new Set(targetGids);
        hideNonHl=true;
        const chk=document.getElementById("chk-hide-nonhl"); if(chk) chk.checked=true;
      }
      setTimeout(focusHighlighted,80);
    };
  }

  // ── apply sort-by-species (hmColSortSp may be a single string or array) ──
  if(hmColSortSp){
    const _sortSps=Array.isArray(hmColSortSp)?hmColSortSp:[hmColSortSp];
    data=[...data].sort((a,b)=>{
      const av=_sortSps.reduce((s,sp)=>s+(a.species_counts[sp]||0),0)/_sortSps.length;
      const bv=_sortSps.reduce((s,sp)=>s+(b.species_counts[sp]||0),0)/_sortSps.length;
      return bv-av;
    });
  }

  // ── apply drag-reorder override ────────────────────────────────────────
  if(hmColOrderOverride && hmColOrderOverride.length){
    const byId=new Map(data.map(d=>[d.id,d]));
    const reordered=hmColOrderOverride.filter(id=>byId.has(id)).map(id=>byId.get(id));
    const inOv=new Set(hmColOrderOverride);
    data.forEach(d=>{ if(!inOv.has(d.id)) reordered.push(d); });
    data=reordered;
  }

  // ── apply free-text search filter (browse modes only) ─────────────────
  if(hmViewMode!=="custom" && hmTextFilter){
    data=data.filter(d=>{
      const s=(d.id||"").toLowerCase()+(d.family||"").toLowerCase()+(d.hg||"").toLowerCase()+(d.class||"").toLowerCase();
      return s.includes(hmTextFilter);
    });
  }

  // ── heatmap row splitting ──────────────────────────────────────────────
  // ── row grouping: driven by species-tree splits OR active highlight queries ──
  // Priority: hmSplitSets > hlSet.  Both reorder rows and draw separator bands.
  let spSplitGroup = new Map(); // species → group index
  let bandLineColors = hmSplitLineColors;
  let bandFillColors = hmSplitColors;
  let bandLabels     = hmSplitLabels;
  let nSplitGroups   = 0;

  if (hmSplitSets.length > 0) {
    hmSplitSets.forEach((sSet,gi) => sSet.forEach(sp => { if (!spSplitGroup.has(sp)) spSplitGroup.set(sp,gi); }));
    nSplitGroups = hmSplitSets.length;
  }

  // Local copy — never mutate the global speciesOrder (it drives drawSpeciesTree too)
  let hmOrder = [...speciesOrder];
  if (nSplitGroups > 0) {
    const grouped = [];
    for (let gi=0; gi<nSplitGroups; gi++)
      hmOrder.filter(s=>spSplitGroup.get(s)===gi).forEach(s=>grouped.push(s));
    hmOrder.filter(s=>!spSplitGroup.has(s)).forEach(s=>grouped.push(s));
    hmOrder = grouped;
  }

  const cW=18, cH=12;
  const maxNameLen=hmOrder.reduce((m,s)=>Math.max(m,s.length),0);
  const ROW_LABEL_W=Math.max(110,Math.min(200,maxNameLen*7+14));
  // truncation limit scales with font size: bigger font = more characters shown
  const colLabelMaxChars = Math.round(hmColFontSize * 3.5);
  const colLabelTrunc = d => { const s=colLabel(d)||""; return s.length>colLabelMaxChars ? s.slice(0,colLabelMaxChars-1)+"\u2026" : s; };
  // top margin: enough room for rotated column labels (using truncated length)
  const maxColLen=data.reduce((m,d)=>Math.max(m,(colLabelTrunc(d)||"").length),0);
  const hmColLabelH=Math.ceil(Math.sin(hmColRotation*Math.PI/180)*maxColLen*hmColFontSize*0.6)+8;
  const hmTM=Math.max(HM_TOP, hmColLabelH+14);
  hmTopActual = hmTM; // always sync before cladogram draws
  const nRows=hmOrder.length;
  const CELLS_BOTTOM=nRows*14+hmTM;   // y just below the last cell row
  const BAR_TOP=CELLS_BOTTOM+10;      // bar chart starts here
  const svgW=data.length*cW+ROW_LABEL_W+20;
  const svgH=BAR_TOP+HM_BAR_H+50;    // extra room for colorbar legend
  const panel=document.getElementById("heatmap-panel");
  const svg=d3.select(panel).html("").append("svg").attr("width",svgW).attr("height",svgH);

  if(!data.length){ svg.append("text").attr("x",40).attr("y",hmTM+30).attr("fill","#999").text("No data for this selection."); return; }

  // draw group background bands before cells (SVG paint order: back to front)
  if (nSplitGroups > 0) {
    let gi0 = spSplitGroup.has(hmOrder[0]) ? spSplitGroup.get(hmOrder[0]) : -1;
    let bandStart = 0;
    for (let ri=1; ri<=hmOrder.length; ri++) {
      const gi1 = ri<hmOrder.length ? (spSplitGroup.has(hmOrder[ri]) ? spSplitGroup.get(hmOrder[ri]) : -1) : -2;
      if (gi1 !== gi0) {
        if (gi0 >= 0) {
          svg.append("rect")
            .attr("x",0).attr("y",bandStart*14+hmTM)
            .attr("width",svgW).attr("height",(ri-bandStart)*14)
            .attr("fill",bandFillColors[gi0%bandFillColors.length]);
          const midY=bandStart*14+hmTM+(ri-bandStart)*7;
          svg.append("text").attr("x",4).attr("y",midY).attr("dy","0.35em")
            .attr("font-size",8).attr("fill",bandLineColors[gi0%bandLineColors.length])
            .attr("font-weight","bold").text(bandLabels[gi0]||("Group "+(gi0+1)));
        }
        if (ri<hmOrder.length) {
          svg.append("line")
            .attr("x1",0).attr("x2",svgW)
            .attr("y1",ri*14+hmTM).attr("y2",ri*14+hmTM)
            .attr("stroke","#999").attr("stroke-width",1.5).attr("stroke-dasharray","4,2");
        }
        gi0 = gi1; bandStart = ri;
      }
    }
  }

  const zMat=data.map(rec=>{
    const vals=hmOrder.map(s=>rec.species_counts[s]||0);
    const m=d3.mean(vals), sd=d3.deviation(vals)||1;
    return vals.map(v=>(v-m)/sd);
  });
  const zAll=zMat.flat();
  const zMax=Math.max(Math.abs(d3.min(zAll)||0),Math.abs(d3.max(zAll)||0),0.01);
  const zColor=d3.scaleDiverging().domain([zMax,0,-zMax]).interpolator(d3.interpolateRdBu);
  // absolute count colour scale
  const allCounts=data.flatMap(rec=>hmOrder.map(s=>rec.species_counts[s]||0));
  const maxAbsCount=d3.max(allCounts)||1;
  const absColor=d3.scaleSequential().domain([0,maxAbsCount]).interpolator(d3.interpolateBlues);
  const color=hmColorMode==="absolute"?absColor:zColor;

  data.forEach((rec,ci)=>{
    hmOrder.forEach((sp,ri)=>{
      const count=rec.species_counts[sp]||0;
      const z=zMat[ci][ri];
      const cellX=ci*cW+ROW_LABEL_W, cellY=ri*14+hmTM;
      const cellFill=hmColorMode==="absolute"?absColor(count):zColor(z);
      svg.append("rect").attr("class","hm-cell")
        .attr("x",cellX).attr("y",cellY).attr("width",cW-2).attr("height",cH)
        .attr("fill",count===0?"#fff":cellFill).style("cursor","pointer")
        .on("mouseover",ev=>{
          showTip(ev,'<b>'+sp+'</b><br>'+rec.id+'<br>count: <b>'+count+'</b>'+(hmColorMode==="zscore"?'<br>z: <b>'+z.toFixed(2)+'</b>':""));
        })
        .on("mousemove",moveTip).on("mouseout",hideTip)
        .on("click",ev=>clickHandler(ev,rec,sp));
      if(count>0){
        const fc=d3.color(cellFill);
        const lum=fc?(0.2126*fc.r+0.7152*fc.g+0.0722*fc.b)/255:1;
        svg.append("text")
          .attr("x",cellX+(cW-2)/2).attr("y",cellY+cH/2).attr("dy","0.35em")
          .attr("text-anchor","middle").attr("font-size",6.5)
          .attr("fill",lum>0.45?"#222":"#eee").attr("pointer-events","none")
          .text(count>999?"≫":count);
      }
    });
  });

  // ── collapsed-clade shading bands in heatmap ─────────────────────────
  hmCollapsedBands.forEach((band,bi)=>{
    const rowsInBand=band.species.map(sp=>hmOrder.indexOf(sp)).filter(i=>i>=0).sort((a,b)=>a-b);
    if(!rowsInBand.length) return;
    const rMin=rowsInBand[0], rMax=rowsInBand[rowsInBand.length-1];
    const bandY=rMin*14+hmTM, bandH=(rMax-rMin+1)*14;
    const shadeFill=bi%2===0?"rgba(180,180,180,0.18)":"rgba(200,200,200,0.12)";
    svg.insert("rect","rect.hm-cell")  // behind cells
      .attr("x",ROW_LABEL_W).attr("y",bandY)
      .attr("width",data.length*cW).attr("height",bandH)
      .attr("fill",shadeFill).attr("pointer-events","none");
    // separator lines at top and bottom
    ["top","bottom"].forEach(edge=>{
      const lineY=edge==="top"?bandY:bandY+bandH;
      svg.append("line")
        .attr("x1",0).attr("x2",ROW_LABEL_W+data.length*cW)
        .attr("y1",lineY).attr("y2",lineY)
        .attr("stroke","#aaa").attr("stroke-width",0.8).attr("stroke-dasharray","3,2")
        .attr("pointer-events","none");
    });
    // label at left edge (offset right to make room for the group sort button)
    svg.append("text")
      .attr("x",14).attr("y",bandY+bandH/2).attr("dy","0.35em")
      .attr("font-size",8).attr("fill","#888").attr("font-style","italic")
      .text(band.label);
    // group sort button spanning the full band height
    const bandSps=band.species.filter(s=>hmOrder.includes(s));
    if(bandSps.length){
      const bActive=hmSortActive(bandSps);
      const bBtnH=Math.max(10,bandH-2);
      svg.append("rect").attr("class","hm-sort-btn")
        .attr("x",1).attr("y",bandY+1).attr("width",10).attr("height",bBtnH).attr("rx",2)
        .attr("fill",bActive?"#4a90d9":"#ddd").attr("stroke",bActive?"#2172c4":"#bbb").attr("stroke-width",0.8)
        .style("cursor","pointer")
        .on("mouseover",ev=>showTip(ev,`Sort OGs by avg count in <b>${band.label}</b> (${bandSps.length} spp.)`))
        .on("mousemove",moveTip).on("mouseout",hideTip)
        .on("click",()=>{ hideTip(); hmSetSort(bandSps); });
      svg.append("text").attr("class","hm-sort-btn-lbl")
        .attr("x",6).attr("y",bandY+bBtnH/2+1).attr("text-anchor","middle").attr("dominant-baseline","middle")
        .attr("font-size",8).attr("fill",bActive?"#fff":"#888").attr("pointer-events","none")
        .text("↓");
    }
  });

  // which species are covered by a collapsed-clade band sort button?
  const spInBand=new Map(); // sp → band index
  hmCollapsedBands.forEach((band,bi)=>{ band.species.forEach(sp=>{ if(hmOrder.includes(sp)) spInBand.set(sp,bi); }); });

  // row labels + individual sort buttons
  hmOrder.forEach((sp,ri)=>{
    const gi = spSplitGroup.get(sp);
    const labelCol = gi!==undefined
      ? bandLineColors[gi%bandLineColors.length]
      : (nSplitGroups>0 ? "#bbb" : "#333");
    svg.append("text")
      .attr("x",ROW_LABEL_W-4).attr("y",ri*14+hmTM+9)
      .attr("text-anchor","end").attr("font-size",11).attr("fill",labelCol)
      .text(sp);
    // individual sort button — skip if covered by a band button
    if(!spInBand.has(sp)){
      const active=hmSortActive([sp]);
      const btnY=ri*14+hmTM+1;
      svg.append("rect").attr("class","hm-sort-btn")
        .attr("x",1).attr("y",btnY).attr("width",10).attr("height",11).attr("rx",2)
        .attr("fill",active?"#4a90d9":"#eee").attr("stroke",active?"#2172c4":"#ccc").attr("stroke-width",0.8)
        .style("cursor","pointer")
        .on("mouseover",ev=>showTip(ev,`Sort OGs by count in <b>${sp}</b>`))
        .on("mousemove",moveTip).on("mouseout",hideTip)
        .on("click",()=>{ hideTip(); hmSetSort([sp]); });
      svg.append("text").attr("class","hm-sort-btn-lbl")
        .attr("x",6).attr("y",btnY+8).attr("text-anchor","middle")
        .attr("font-size",8).attr("fill",active?"#fff":"#888").attr("pointer-events","none")
        .text("↓");
    }
  });

  // column headers — drag to reorder, short-click to navigate
  // Drag ghost + insert marker (rendered above everything else at end of drag)
  const _dragGhost=svg.append("rect").attr("y",0).attr("height",svgH-28)
    .attr("fill","rgba(74,144,217,0.13)").attr("stroke","#4a90d9").attr("stroke-width",1)
    .attr("opacity",0).attr("pointer-events","none");
  const _dragMarker=svg.append("line").attr("y1",0).attr("y2",svgH-28)
    .attr("stroke","#4a90d9").attr("stroke-width",2.5)
    .attr("opacity",0).attr("pointer-events","none");
  let _hmDragIdx=null,_hmDragTgt=null,_hmDragMoved=false,_hmDragSX=0;
  svg.selectAll("text.hm-col").data(data).enter().append("text").attr("class","hm-col")
    .attr("transform",(d,i)=>`translate(${i*cW+ROW_LABEL_W+cW/2},${hmTM-6}) rotate(-${hmColRotation})`)
    .attr("font-size",hmColFontSize).style("cursor","grab")
    .attr("text-anchor","start")
    .text(colLabelTrunc)
    .call(d3.drag()
      .on("start",(event,d)=>{
        _hmDragIdx=data.findIndex(r=>r.id===d.id);
        _hmDragMoved=false; _hmDragSX=event.x;
        _dragGhost.attr("x",_hmDragIdx*cW+ROW_LABEL_W).attr("width",cW).attr("opacity",1);
        svg.style("cursor","grabbing");
      })
      .on("drag",(event)=>{
        if(Math.abs(event.x-_hmDragSX)>5) _hmDragMoved=true;
        if(!_hmDragMoved) return;
        const tgt=Math.max(0,Math.min(data.length,Math.round((event.x-ROW_LABEL_W)/cW)));
        _hmDragTgt=tgt;
        _dragMarker.attr("x1",tgt*cW+ROW_LABEL_W).attr("x2",tgt*cW+ROW_LABEL_W).attr("opacity",1);
      })
      .on("end",(event,d)=>{
        _dragGhost.attr("opacity",0); _dragMarker.attr("opacity",0);
        svg.style("cursor","");
        if(!_hmDragMoved){
          clickHandler(event,d,null);
        } else if(_hmDragIdx!==null&&_hmDragTgt!==null&&_hmDragTgt!==_hmDragIdx&&_hmDragTgt!==_hmDragIdx+1){
          const nd=[...data]; const [mv]=nd.splice(_hmDragIdx,1);
          const ins=_hmDragTgt>_hmDragIdx?_hmDragTgt-1:_hmDragTgt;
          nd.splice(ins,0,mv);
          hmColOrderOverride=nd.map(r=>r.id);
          drawHeatmap();
        }
        _hmDragIdx=null; _hmDragTgt=null;
      }));

  // ── column-sum bar chart (below cells) ───────────────────────────────────
  const colTotals=data.map(rec=>d3.sum(hmOrder.map(s=>rec.species_counts[s]||0)));
  const maxTotal=d3.max(colTotals)||1;
  // baseline
  svg.append("line")
    .attr("x1",ROW_LABEL_W).attr("x2",ROW_LABEL_W+data.length*cW)
    .attr("y1",BAR_TOP).attr("y2",BAR_TOP)
    .attr("stroke","#ccc").attr("stroke-width",0.8);
  svg.append("text").attr("x",ROW_LABEL_W-4).attr("y",BAR_TOP+HM_BAR_H/2)
    .attr("text-anchor","end").attr("font-size",8).attr("fill","#aaa")
    .attr("dominant-baseline","middle").text("total \u2191");
  data.forEach((_rec,ci)=>{
    const total=colTotals[ci];
    const bH=Math.max(total>0?1:0,(total/maxTotal)*HM_BAR_H);
    svg.append("rect")
      .attr("x",ci*cW+ROW_LABEL_W+1).attr("y",BAR_TOP)
      .attr("width",cW-2).attr("height",bH)
      .attr("fill","#7fb3d3").attr("opacity",0.8);
    if(total>0)
      svg.append("text")
        .attr("x",ci*cW+ROW_LABEL_W+cW/2).attr("y",BAR_TOP+bH+8)
        .attr("text-anchor","middle").attr("font-size",6.5).attr("fill","#555")
        .text(total>=10000?(total/1000).toFixed(0)+"k":total);
  });

  // colorbar legend
  const cbW=120, cbH=10, cbX=ROW_LABEL_W, cbY=svgH-28;
  const gradId="hmgrad"+Date.now();
  const defs=svg.append("defs");
  const grad=defs.append("linearGradient").attr("id",gradId).attr("x1","0%").attr("x2","100%");
  const stops=10;
  for(let i=0;i<=stops;i++){
    const t=i/stops;
    grad.append("stop").attr("offset",(t*100)+"%").attr("stop-color",color(zMax-t*2*zMax));
  }
  svg.append("rect").attr("x",cbX).attr("y",cbY).attr("width",cbW).attr("height",cbH)
    .attr("fill","url(#"+gradId+")").attr("rx",2);
  svg.append("text").attr("x",cbX).attr("y",cbY+cbH+9).attr("font-size",8).attr("fill","#888").attr("text-anchor","middle").text("high");
  svg.append("text").attr("x",cbX+cbW/2).attr("y",cbY+cbH+9).attr("font-size",8).attr("fill","#888").attr("text-anchor","middle").text("mean");
  svg.append("text").attr("x",cbX+cbW).attr("y",cbY+cbH+9).attr("font-size",8).attr("fill","#888").attr("text-anchor","middle").text("low");

  // always redraw cladogram last so it uses the updated hmTopActual
  drawCladogram();
}

function hmBack(){
  if(hmViewMode==="og")     { hmViewMode="hg";     hmActiveHG=null; hmActiveHGRec=null; }
  else if(hmViewMode==="hg"){ hmViewMode="family";  hmActiveFamily=null; hmActiveHG=null; hmActiveHGRec=null; }
  else                      { hmViewMode="class";   hmActiveClass=null; hmActiveFamily=null; hmActiveHG=null; hmActiveHGRec=null; }
  drawHeatmap();
}

// ── Heatmap unified search (navigate + custom OG/gene) ─────────────────────
let _hmSearchHits=[], _hmSearchSel=-1;
function hmTextSearchInput(val){
  buildHmSearchIndex(); buildOGIndex();
  const dd=document.getElementById("hm-search-dd");
  hmTextFilter=val.trim().toLowerCase(); drawHeatmap();
  if(!val.trim()){ dd.style.display="none"; return; }
  const q=val.toLowerCase();
  _hmSearchHits=[];

  // ── Navigate section: class / family / HG ────────────────────────────
  const navHits=hmSearchIndex.filter(r=>r.label.toLowerCase().includes(q)).slice(0,15)
    .map(h=>({...h, kind:"nav"}));

  // ── Custom section: class/family groups, individual OGs, gene matches ─
  const effOGs=new Set(getEffectiveCustomOGs());
  const clsMap=new Map(), famMap=new Map();
  for(const [og,r] of Object.entries(hmOGIndex)){
    if(og.toLowerCase().includes(q)||r.cls.toLowerCase().includes(q))
      if(r.cls){ if(!clsMap.has(r.cls)) clsMap.set(r.cls,[]); clsMap.get(r.cls).push(og); }
    if(og.toLowerCase().includes(q)||r.family.toLowerCase().includes(q))
      if(r.family){ if(!famMap.has(r.family)) famMap.set(r.family,[]); famMap.get(r.family).push(og); }
  }
  const clsGroupHits=[...clsMap.entries()].filter(([k])=>k.toLowerCase().includes(q))
    .sort((a,b)=>b[1].length-a[1].length).slice(0,5)
    .map(([k,ogs])=>({kind:"og",type:"group",groupType:"class",key:k,ogs,count:ogs.length}));
  const famGroupHits=[...famMap.entries()].filter(([k])=>k.toLowerCase().includes(q))
    .sort((a,b)=>b[1].length-a[1].length).slice(0,5)
    .map(([k,ogs])=>({kind:"og",type:"group",groupType:"family",key:k,ogs,count:ogs.length}));
  const ogHits=Object.keys(hmOGIndex).filter(k=>k.toLowerCase().includes(q))
    .sort((a,b)=>{const as=a.toLowerCase().startsWith(q),bs=b.toLowerCase().startsWith(q);
      if(as!==bs) return as?-1:1; return (hmOGIndex[b].total||0)-(hmOGIndex[a].total||0);})
    .slice(0,20).map(og=>({kind:"og",type:"og",og,matchedGene:null}));
  const seenOGs=new Set(ogHits.map(h=>h.og));
  const geneHits=Object.entries(hmGeneIndex).filter(([g])=>g.toLowerCase().includes(q))
    .filter(([,og])=>!seenOGs.has(og))
    .sort((a,b)=>(hmOGIndex[b[1]]?.total||0)-(hmOGIndex[a[1]]?.total||0)).slice(0,10)
    .map(([g,og])=>{ seenOGs.add(og); return {kind:"og",type:"og",og,matchedGene:g}; });
  const ogSection=[...clsGroupHits,...famGroupHits,...ogHits,...geneHits];

  _hmSearchHits=[
    ...navHits,
    ...(ogSection.length?[{kind:"sep",label:"Add to custom selection"}]:[]),
    ...ogSection,
  ];
  _hmSearchSel=-1;

  dd.innerHTML=_hmSearchHits.map((h,i)=>{
    if(h.kind==="sep") return `<div style="padding:3px 10px;font-size:10px;color:#aaa;background:#f5f5f5;border-top:1px solid #eee;border-bottom:1px solid #eee;font-weight:600;letter-spacing:.05em;text-transform:uppercase">${h.label}</div>`;
    if(h.kind==="nav") return `<div data-i="${i}" style="padding:4px 10px;cursor:pointer;border-bottom:1px solid #f0f0f0;color:${h.type==="hg"?"#333":h.type==="family"?"#555":"#888"}" onmousedown="hmTextSearchGo(${i})">&#8594; ${h.label}</div>`;
    if(h.type==="group"){
      const done=hmCustomGroups.some(g=>g.groupType===h.groupType&&g.key===h.key)||h.ogs.every(og=>effOGs.has(og));
      const icon=h.groupType==="class"?"\uD83D\uDCC2":"\uD83E\uDDF9";
      return `<div data-i="${i}" style="padding:4px 10px;cursor:pointer;border-bottom:1px solid #f0f0f0;background:${done?"#e8f5e9":"#f0f6ff"};color:${done?"#2e7d32":"#1a4a7a"}" onmousedown="hmTextSearchGo(${i})">${icon} <b>${h.key}</b><span style="color:#aaa;font-size:10px"> \u00b7 ${h.groupType} \u00b7 ${h.count} OGs</span>${done?`<span style="float:right">&#10003;</span>`:""}</div>`;
    }
    const r=hmOGIndex[h.og]||{family:"",total:0}; const already=effOGs.has(h.og);
    const lbl=h.matchedGene?`<span style="color:#888;font-size:10px">gene</span> <b>${h.matchedGene}</b><span style="color:#aaa;font-size:10px"> \u2192 ${h.og} \u00b7 ${r.family}</span>`:`<b>${h.og}</b><span style="color:#aaa;font-size:10px"> \u00b7 ${r.family} \u00b7 ${r.total} genes</span>`;
    return `<div data-i="${i}" style="padding:4px 10px;cursor:pointer;border-bottom:1px solid #f0f0f0;background:${already?"#e8f5e9":"#fff"};color:${already?"#2e7d32":"#222"}" onmousedown="hmTextSearchGo(${i})">${lbl}${already?`<span style="float:right;color:#27ae60">&#10003;</span>`:""}</div>`;
  }).join("")||`<div style="padding:6px 10px;color:#aaa">No matches</div>`;
  dd.style.display="block";
}
function hmTextSearchKey(ev){
  const dd=document.getElementById("hm-search-dd");
  const items=[...dd.querySelectorAll("[data-i]")];
  if(ev.key==="ArrowDown"){ ev.preventDefault();
    const ci=items.findIndex(el=>+el.dataset.i===_hmSearchSel);
    const next=items[Math.min(ci+1,items.length-1)]; if(next){ _hmSearchSel=+next.dataset.i; items.forEach(el=>el.style.background=+el.dataset.i===_hmSearchSel?"#e3f0ff":""); }
  } else if(ev.key==="ArrowUp"){ ev.preventDefault();
    const ci=items.findIndex(el=>+el.dataset.i===_hmSearchSel);
    const prev=items[Math.max(ci-1,0)]; if(prev){ _hmSearchSel=+prev.dataset.i; items.forEach(el=>el.style.background=+el.dataset.i===_hmSearchSel?"#e3f0ff":""); }
  } else if(ev.key==="Enter"){ ev.preventDefault();
    const idx=_hmSearchSel>=0?_hmSearchSel:(items.length===1?+items[0].dataset.i:-1);
    if(idx>=0) hmTextSearchGo(idx);
  } else if(ev.key==="Escape"){ dd.style.display="none"; }
}
function hmTextSearchGo(i){
  const h=_hmSearchHits[i]; if(!h||h.kind==="sep") return;
  document.getElementById("hm-search-dd").style.display="none";
  document.getElementById("hm-text-search").value=""; hmTextFilter="";
  if(h.kind==="nav"){
    if(h.type==="class"){ hmViewMode="class"; hmActiveClass=h.cls; hmActiveFamily=hmActiveHG=hmActiveHGRec=null; }
    else if(h.type==="family"){ hmViewMode="family"; hmActiveClass=h.cls; hmActiveFamily=h.fam; hmActiveHG=hmActiveHGRec=null; }
    else if(h.type==="hg"){
      const rec=TREE_INDEX.find(r=>r.id===h.hg_id); if(!rec) return;
      hmViewMode="og"; hmActiveClass=h.cls; hmActiveFamily=h.fam; hmActiveHG=h.hg_id; hmActiveHGRec=rec;
    }
  } else { // add to custom selection
    if(h.type==="group"){
      if(!hmCustomGroups.some(g=>g.groupType===h.groupType&&g.key===h.key))
        hmCustomGroups.push({groupType:h.groupType,key:h.key,label:h.key,ogs:h.ogs});
    } else {
      if(!hmCustomOGs.includes(h.og)) hmCustomOGs.push(h.og);
    }
    hmViewMode="custom";
    const bar=document.getElementById("hm-custom-bar");
    bar.style.display="flex"; _positionCustomBar();
    renderCustomChips();
  }
  drawHeatmap();
}

// ── Custom OG selection helpers ────────────────────────────────────────────
function _positionCustomBar(){
  const strip=document.getElementById("hm-col-strip");
  const bar=document.getElementById("hm-custom-bar");
  if(!strip||bar.style.display==="none") return;
  const r=strip.getBoundingClientRect();
  bar.style.top=r.bottom+"px";
  bar.style.left=r.left+"px";
  bar.style.width=(r.right-r.left)+"px";
}
function hmEnterCustom(){
  buildOGIndex();
  hmViewMode="custom";
  const bar=document.getElementById("hm-custom-bar");
  bar.style.display="flex";
  _positionCustomBar();
  document.getElementById("hm-text-search").focus();
  drawHeatmap();
}
function hmExpandHGToOGs(hgIds){
  buildOGIndex();
  const newGroups=hgIds.map(hgId=>{
    const ogs=Object.keys(hmOGIndex).filter(og=>hmOGIndex[og].hgId===hgId);
    const label=HG_DATA.find(d=>d.id===hgId)?.hg||hgId;
    return {groupType:"hg",key:hgId,label,ogs};
  }).filter(g=>g.ogs.length>0);
  if(!newGroups.length) return;
  newGroups.forEach(g=>{ if(!hmCustomGroups.some(x=>x.key===g.key)) hmCustomGroups.push(g); });
  hmViewMode="custom";
  const bar=document.getElementById("hm-custom-bar");
  bar.style.display="flex"; _positionCustomBar();
  renderCustomChips(); drawHeatmap();
}
function hmExpandToOGs(){
  const btn=document.getElementById("hm-expand-og");
  const hgIds=JSON.parse(btn.dataset.hgData||"[]");
  hmExpandHGToOGs(hgIds);
}
function hmExitCustom(){
  hmViewMode="class"; hmActiveClass=null; hmActiveFamily=null;
  hmActiveHG=null; hmActiveHGRec=null;
  hmColOrderOverride=null;
  document.getElementById("hm-custom-bar").style.display="none";
  drawHeatmap();
}
window.addEventListener("resize",_positionCustomBar);
function hmCustomClear(){
  hmCustomOGs=[]; hmCustomGroups=[]; hmColOrderOverride=null;
  renderCustomChips(); drawHeatmap();
}


function renderCustomChips(){
  const grpChips=hmCustomGroups.map((g,i)=>{
    const icon=g.groupType==="class"?"\uD83D\uDCC2":"\uD83E\uDDF9";
    return `<span style="display:inline-flex;align-items:center;gap:3px;padding:1px 7px;border-radius:10px;background:#2c6e9e;color:#fff;font-size:10px">${icon} ${g.label} <span style="font-weight:normal;opacity:.8">(${g.ogs.length})</span><span style="cursor:pointer;margin-left:2px;opacity:.8" onclick="hmCustomGroups.splice(${i},1);renderCustomChips();drawHeatmap()">&#10005;</span></span>`;
  }).join("");
  const ogChips=hmCustomOGs.map((og,i)=>`<span style="display:inline-flex;align-items:center;gap:3px;padding:1px 7px;border-radius:10px;background:#4a90d9;color:#fff;font-size:10px">${og}<span style="cursor:pointer;margin-left:2px;opacity:.8" onclick="hmCustomOGs.splice(${i},1);renderCustomChips();drawHeatmap()">&#10005;</span></span>`).join("");
  document.getElementById("hm-custom-chips").innerHTML=grpChips+ogChips;
}

// ═══════════════════════════════════════════════════════════════════════════════
// TREE VIEW – SIDEBAR
// ═══════════════════════════════════════════════════════════════════════════════
let currentIndex  = null;
let currentDetail = null;
let treeSource    = "generax";   // "generax" | "original"

function buildOGIndex(){
  if(hmOGIndex) return;
  hmOGIndex={}; hmGeneIndex={};
  for(const tRec of TREE_INDEX){
    const detail=loadDetail(tRec.id);
    if(!detail) continue;
    const ogs=detail.ogs||{};
    for(const [og,gids] of Object.entries(ogs)){
      const sc={};
      gids.forEach(g=>{
        const sp=getSpeciesPfx(g);
        if(sp) sc[sp]=(sc[sp]||0)+1;
        hmGeneIndex[g]=og;
      });
      hmOGIndex[og]={hgId:tRec.id,family:tRec.family,cls:tRec.class||tRec.prefix,
                     total:gids.length,species_counts:sc};
    }
    const UNCL="Unclassified";
    const oggedGenes=new Set(Object.values(ogs).flat());
    (function walk(n){
      if(!n) return;
      if(n.leaf){
        const gid=n.gene_id||n.name||"";
        if(gid&&!oggedGenes.has(gid)){
          hmGeneIndex[gid]=UNCL;
          if(!hmOGIndex[UNCL]) hmOGIndex[UNCL]={hgId:"",family:"",cls:UNCL,total:0,species_counts:{}};
          const sp=getSpeciesPfx(gid);
          if(sp){ hmOGIndex[UNCL].species_counts[sp]=(hmOGIndex[UNCL].species_counts[sp]||0)+1; }
          hmOGIndex[UNCL].total++;
        }
      }
      (n.children||[]).forEach(walk);
    })(detail.tree);
  }
}

function buildHmSearchIndex(){
  if(hmSearchIndex) return;
  hmSearchIndex=[];
  const classes=[...new Set(FAMILY_DATA.map(d=>d.class||"").filter(Boolean))].sort();
  classes.forEach(cls=>hmSearchIndex.push({label:cls,type:"class",cls,fam:null,hg_id:null}));
  const famSeen=new Set();
  FAMILY_DATA.forEach(d=>{
    const key=(d.class||"")+"\x00"+(d.family||"");
    if(!famSeen.has(key)){ famSeen.add(key);
      hmSearchIndex.push({label:(d.class?"["+d.class+"] ":"")+d.family,
                          type:"family",cls:d.class||null,fam:d.family,hg_id:null}); }
  });
  TREE_INDEX.forEach(r=>hmSearchIndex.push({
    label:(r.class?"["+r.class+"] ":"")+(r.family?r.family+" \u203a ":"")+r.hg,
    type:"hg",cls:r.class||null,fam:r.family||null,hg_id:r.id}));
}

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
      const badge=rec.source==="generax"?'<span class="src-badge">GeneRax</span>':'';
      const covPct=ALL_SPECIES.length?Math.round((rec.species||[]).length/ALL_SPECIES.length*100):0;
      item.innerHTML='<div class="hg-name">'+rec.hg+' '+badge+'</div>'
        +'<div class="hg-cov"><div class="hg-cov-bar" style="width:'+covPct+'%"></div></div>'
        +'<div class="hg-meta">'+rec.n_leaves+' genes \u00b7 '+rec.n_ogs+' OGs \u00b7 '+covPct+'% sp.</div>';
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
    if(hlSet) hmFocusGids=null; // user's new highlight overrides heatmap-navigation focus
  }
  renderHlTags();
  document.getElementById("btn-focus-hl").style.display=(hlSet||ogHlSet)?"inline":"none";
  if(currentIndex!==null) renderTree();
}

function renderHlTags(){
  const el=document.getElementById("hl-tags"); el.innerHTML="";
  hlQueries.forEach((q,i)=>{
    const col=hlTagColor(i);
    const chip=document.createElement("span"); chip.className="hl-tag";
    chip.style.background=col; chip.title="Click to change color";
    chip.style.cursor="pointer";
    const lbl=document.createTextNode(q+" ");
    const x=document.createElement("span"); x.className="hl-tag-x"; x.textContent="\u00d7";
    x.onclick=(e)=>{ e.stopPropagation(); removeHlTag(i); };
    chip.onclick=()=>openColorPicker(col,c=>{ hlQueryColors[q]=c; rebuildHlSet(); });
    chip.appendChild(lbl); chip.appendChild(x);
    el.appendChild(chip);
  });
}

function addHlTag(query){
  query=(query||"").trim();
  if(!query||hlQueries.includes(query)) return;
  hlQueries.push(query);
  document.getElementById("hl-search").value="";
  // auto-switch to "by species" so the highlight is visible
  if(colorMode!=="species"){
    const sel=document.getElementById("color-by");
    sel.value="species";
    sel.dispatchEvent(new Event("change"));
  }
  rebuildHlSet();
}

function removeHlTag(i){ delete hlQueryColors[hlQueries[i]]; hlQueries.splice(i,1); rebuildHlSet(); }

function clearHighlight(){ hlQueries=[]; hlQueryColors={}; document.getElementById("hl-search").value=""; rebuildHlSet(); }

// ── OG name highlight system ──────────────────────────────────────────────
function resolveOgQuery(query){
  // Returns Set<og_name> whose names contain the query string (case-insensitive)
  const lq=(query||"").toLowerCase().trim();
  if(!lq) return new Set();
  const matched=new Set();
  // Search through all known OG names in the current tree
  Object.values(ogGene2Name).forEach(og=>{ if(og.toLowerCase().includes(lq)) matched.add(og); });
  // Also search ogName2Color keys
  Object.keys(ogName2Color).forEach(og=>{ if(og.toLowerCase().includes(lq)) matched.add(og); });
  return matched;
}

function rebuildOgHlSet(){
  ogHlGroupIndex=new Map();
  if(!ogHlQueries.length){ ogHlSet=null; }
  else {
    const union=new Set();
    ogHlQueries.forEach((q,i)=>{ resolveOgQuery(q).forEach(og=>{ union.add(og); if(!ogHlGroupIndex.has(og)) ogHlGroupIndex.set(og,i); }); });
    ogHlSet=union.size?union:null;
    if(ogHlSet) hmFocusGids=null; // user's new highlight overrides heatmap-navigation focus
  }
  renderOgHlTags();
  document.getElementById("btn-focus-hl").style.display=(ogHlSet||hlSet)?"inline":"none";
  if(currentIndex!==null) renderTree();
}

function renderOgHlTags(){
  const el=document.getElementById("og-hl-tags"); el.innerHTML="";
  ogHlQueries.forEach((q,i)=>{
    const col=ogHlTagColor(i);
    const chip=document.createElement("span"); chip.className="hl-tag";
    chip.style.background=col; chip.title="Click to change color";
    chip.style.cursor="pointer";
    const lbl=document.createTextNode(q+" ");
    const btn=document.createElement("button");
    btn.textContent="\u00d7"; btn.onclick=(e)=>{ e.stopPropagation(); removeOgHlTag(i); };
    chip.onclick=()=>openColorPicker(col,c=>{ ogHlQueryColors[q]=c; rebuildOgHlSet(); });
    chip.append(lbl,btn); el.appendChild(chip);
  });
  // (datalist removed – autocomplete is now handled by the live og-hl-dd dropdown)
}

function addOgHlTag(query){
  query=(query||"").trim();
  if(!query||ogHlQueries.includes(query)) return;
  ogHlQueries.push(query);
  document.getElementById("og-hl-search").value="";
  ogHlHideDD();
  rebuildOgHlSet();
}
function removeOgHlTag(i){ delete ogHlQueryColors[ogHlQueries[i]]; ogHlQueries.splice(i,1); rebuildOgHlSet(); }
function clearOgHighlight(){ ogHlQueries=[]; ogHlQueryColors={}; document.getElementById("og-hl-search").value=""; rebuildOgHlSet(); }

// ── OG highlight search dropdown ──────────────────────────────────────────────
let _ogHlDDSel=-1;
function ogHlHideDD(){ document.getElementById("og-hl-dd").style.display="none"; _ogHlDDSel=-1; }
function ogHlSearchInput(val){
  const dd=document.getElementById("og-hl-dd");
  val=(val||"").trim().toLowerCase();
  if(!val){ ogHlHideDD(); return; }
  const allOgs=Array.from(new Set([...Object.values(ogGene2Name),...Object.keys(ogName2Color)]));
  const hits=allOgs.filter(og=>og.toLowerCase().includes(val)).slice(0,40);
  if(!hits.length){ ogHlHideDD(); return; }
  _ogHlDDSel=-1;
  dd.innerHTML=hits.map((og,i)=>`<div data-i="${i}" data-og="${og}" style="padding:4px 10px;cursor:pointer;border-bottom:1px solid #f0f0f0" onmousedown="event.preventDefault();addOgHlTag('${og.replace(/'/g,"\\'")}');document.getElementById('og-hl-search').value=''">${og}</div>`).join("");
  dd.style.display="block";
}
function ogHlSearchKey(e){
  const dd=document.getElementById("og-hl-dd");
  const items=dd.querySelectorAll("div[data-og]");
  if(e.key==="ArrowDown"){ e.preventDefault(); _ogHlDDSel=Math.min(_ogHlDDSel+1,items.length-1); items.forEach((el,i)=>el.style.background=i===_ogHlDDSel?"#e8f0fe":""); return; }
  if(e.key==="ArrowUp"){ e.preventDefault(); _ogHlDDSel=Math.max(_ogHlDDSel-1,0); items.forEach((el,i)=>el.style.background=i===_ogHlDDSel?"#e8f0fe":""); return; }
  if(e.key==="Enter"){
    e.preventDefault();
    if(_ogHlDDSel>=0&&items[_ogHlDDSel]){ addOgHlTag(items[_ogHlDDSel].dataset.og); document.getElementById("og-hl-search").value=""; }
    else if(document.getElementById("og-hl-search").value.trim()) addOgHlTag(document.getElementById("og-hl-search").value.trim());
    return;
  }
  if(e.key==="Escape"){ ogHlHideDD(); return; }
}

document.getElementById("hl-search").addEventListener("keydown",function(e){
  if(e.key==="Enter"&&this.value.trim()){ addHlTag(this.value); e.preventDefault(); }
});
document.getElementById("hl-search").addEventListener("change",function(){
  if(this.value.trim()) addHlTag(this.value);
});

document.getElementById("tip-font-slider").addEventListener("input",function(){
  tipFontSize=+this.value;
  document.getElementById("tip-font-val").textContent=this.value;
  applyTipFontSize();
});

document.getElementById("line-width-slider").addEventListener("input",function(){
  treeLinkWidth=+this.value;
  document.getElementById("line-width-val").textContent=this.value;
  if(gMain) gMain.selectAll(".link").attr("stroke-width",treeLinkWidth/_zoomScale);
  if(gMain) gMain.selectAll(".col-tri").attr("stroke-width",treeLinkWidth/_zoomScale);
});

document.getElementById("tree-height-slider").addEventListener("input",function(){
  treeHeightMult=+this.value;
  document.getElementById("tree-height-val").textContent=this.value;
  if(rootNode) renderTree(false);
});

document.getElementById("tree-width-slider").addEventListener("input",function(){
  treeWidthMult=+this.value;
  document.getElementById("tree-width-val").textContent=parseFloat(this.value).toFixed(1);
  if(rootNode) renderTree(false);
});

document.getElementById("clade-hl-alpha-slider").addEventListener("input",function(){
  cladeHlAlpha=+this.value;
  document.getElementById("clade-hl-alpha-val").textContent=parseFloat(this.value).toFixed(2);
  if(rootNode) renderTree(false);
});

document.getElementById("hm-col-font-slider").addEventListener("input",function(){
  hmColFontSize=+this.value;
  document.getElementById("hm-col-font-val").textContent=this.value;
  drawHeatmap();
});
document.getElementById("hm-col-rot-slider").addEventListener("input",function(){
  hmColRotation=+this.value;
  document.getElementById("hm-col-rot-val").textContent=this.value;
  drawHeatmap();
});
document.getElementById("hm-color-mode").addEventListener("change",function(){
  hmColorMode=this.value; drawHeatmap();
});
document.getElementById("collapsed-frac-slider").addEventListener("input",function(){
  collapsedFraction=+this.value;
  document.getElementById("collapsed-frac-val").textContent=parseFloat(this.value).toFixed(2);
  if(rootNode) renderTree(true);
});

document.getElementById("chk-geneid").addEventListener("change",function(){ showGeneId=this.checked; if(currentIndex!==null) renderTree(); });
document.getElementById("chk-og").addEventListener("change",function(){ showOGName=this.checked; if(currentIndex!==null) renderTree(); });
document.getElementById("chk-ref").addEventListener("change",function(){ showRefOrtho=this.checked; if(currentIndex!==null) renderTree(); });
document.getElementById("chk-hide-nonhl").addEventListener("change",function(){ hideNonHl=this.checked; if(currentIndex!==null) renderTree(); });

document.getElementById("sptree-width-slider").addEventListener("input",function(){
  spTreeWidthPct=+this.value;
  document.getElementById("sptree-width-val").textContent=this.value;
  drawSpeciesTree();
});
document.getElementById("col-tri-fill").addEventListener("input",function(){
  colTriFill=this.value;
  // only gene-tree .col-tri and the species-tree view use this global colour
  document.querySelectorAll("#tree-svg .col-tri").forEach(el=>{ el.style.fill=colTriFill; });
  drawSpeciesTree();
});

function _buildCombinedHeatmapSVG(){
  const cladoSvg=document.querySelector("#tree-panel svg");
  const heatSvg=document.querySelector("#heatmap-panel svg");
  if(!heatSvg) return null;
  const bb=el=>{ const v=el.viewBox.baseVal; return v.width?{w:v.width,h:v.height}:{w:el.getBoundingClientRect().width,h:el.getBoundingClientRect().height}; };
  const cBB=cladoSvg?bb(cladoSvg):{w:0,h:0};
  const hBB=bb(heatSvg);
  const W=cBB.w+hBB.w, H=Math.max(cBB.h,hBB.h);
  // collect SVG-relevant CSS from the document so class-based styles render correctly
  let css="";
  try{ for(const sh of document.styleSheets){ try{ for(const r of sh.cssRules) css+=r.cssText+"\n"; }catch(e){} } }catch(e){}
  const cladoInner=cladoSvg?cladoSvg.innerHTML:"";
  const heatInner=heatSvg.innerHTML;
  return {svgStr:
    `<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${W}" height="${H}" style="font-family:sans-serif">` +
    `<style>text{font-family:sans-serif}${css}</style>` +
    `<rect width="${W}" height="${H}" fill="#fff"/>` +
    (cladoInner?`<g>${cladoInner}</g>`:"") +
    `<g transform="translate(${cBB.w},0)">${heatInner}</g>` +
    `</svg>`,
    W, H};
}
function downloadHeatmapSVG(){
  const r=_buildCombinedHeatmapSVG();
  if(!r){ alert("No heatmap to export."); return; }
  const blob=new Blob([r.svgStr],{type:"image/svg+xml;charset=utf-8"});
  const a=document.createElement("a"); a.download="heatmap.svg"; a.href=URL.createObjectURL(blob);
  a.click(); URL.revokeObjectURL(a.href);
}
function downloadHeatmapPNG(){
  const r=_buildCombinedHeatmapSVG();
  if(!r){ alert("No heatmap to export."); return; }
  const {svgStr,W,H}=r;
  const blob=new Blob([svgStr],{type:"image/svg+xml;charset=utf-8"});
  const url=URL.createObjectURL(blob);
  const img=new Image();
  img.onload=function(){
    const canvas=document.createElement("canvas");
    canvas.width=W*2; canvas.height=H*2;
    const ctx=canvas.getContext("2d");
    ctx.scale(2,2); ctx.fillStyle="#fff"; ctx.fillRect(0,0,W,H);
    ctx.drawImage(img,0,0,W,H);
    URL.revokeObjectURL(url);
    const a=document.createElement("a"); a.download="heatmap.png";
    a.href=canvas.toDataURL("image/png"); a.click();
  };
  img.src=url;
}

function _getTreeSVGSrc(){
  const svgEl=document.getElementById("tree-svg");
  if(!svgEl||!rootNode) return null;
  const W=+svgEl.getAttribute("width")||800, H=+svgEl.getAttribute("height")||600;
  const serial=new XMLSerializer();
  const css=Array.from(document.styleSheets).flatMap(s=>{try{return Array.from(s.cssRules).map(r=>r.cssText);}catch(e){return[];}}).join("\n");
  const inner=serial.serializeToString(svgEl).replace(/^<svg[^>]*>/,"").replace(/<\/svg>$/,"");
  return {src:'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"'
    +' width="'+W+'" height="'+H+'" style="font-family:sans-serif">'
    +'<style>text{font-family:sans-serif}'+css+'</style>'
    +inner+'</svg>', W, H};
}
function downloadTreeSVG(){
  const r=_getTreeSVGSrc(); if(!r){ alert("No tree to export."); return; }
  const a=document.createElement("a"); a.download="gene_tree.svg";
  a.href=URL.createObjectURL(new Blob([r.src],{type:"image/svg+xml"})); a.click();
  URL.revokeObjectURL(a.href);
}
function downloadTreePNG(){
  const r=_getTreeSVGSrc(); if(!r){ alert("No tree to export."); return; }
  const url=URL.createObjectURL(new Blob([r.src],{type:"image/svg+xml"}));
  const img=new Image(); img.onload=()=>{
    const canvas=document.createElement("canvas");
    canvas.width=r.W*2; canvas.height=r.H*2;
    const ctx=canvas.getContext("2d"); ctx.scale(2,2);
    ctx.fillStyle="#fff"; ctx.fillRect(0,0,r.W,r.H); ctx.drawImage(img,0,0);
    const a=document.createElement("a"); a.download="gene_tree.png";
    a.href=canvas.toDataURL("image/png"); a.click();
    URL.revokeObjectURL(url);
  };
  img.src=url;
}

function downloadNewick(){
  if(!SP_TREE_DATA||!SP_TREE_DATA.children){ alert("No species tree loaded."); return; }
  // Regenerate Newick from (possibly edited) SP_TREE_DATA so label edits are preserved
  const nwk=treeToNewick(SP_TREE_DATA)+";";
  const a=document.createElement("a");
  a.href=URL.createObjectURL(new Blob([nwk],{type:"text/plain"}));
  a.download="species_tree.nwk";
  a.click();
  URL.revokeObjectURL(a.href);
}

function downloadSpeciesGenes(sp, cls){
  buildOGIndex();
  const rows=["gene_id\tog_name\thg_id\tfamily\tclass"];
  for(const [gene,og] of Object.entries(hmGeneIndex)){
    if(getSpeciesPfx(gene)!==sp) continue;
    const r=hmOGIndex[og]||{};
    if(cls && (r.cls||"other")!==cls) continue;
    rows.push([gene,og,r.hgId||"",r.family||"",r.cls||""].join("\t"));
  }
  if(rows.length===1){ alert("No gene data found for "+sp+(cls?" ("+cls+")":"")); return; }
  const blob=new Blob([rows.join("\n")],{type:"text/tab-separated-values"});
  const a=document.createElement("a");
  a.href=URL.createObjectURL(blob); a.download=sp+(cls?"_"+cls:"")+"_genes.tsv";
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(a.href);
}

function showSpAnnotPopup(event, sp){
  buildOGIndex();
  const classes=new Set();
  for(const [gene,og] of Object.entries(hmGeneIndex)){
    if(getSpeciesPfx(gene)!==sp) continue;
    classes.add((hmOGIndex[og]||{}).cls||"other");
  }
  const sorted=[...classes].sort();
  if(!sorted.length){ downloadSpeciesGenes(sp); return; }
  const pop=document.getElementById("sp-annot-popup");
  document.getElementById("sp-annot-popup-title").textContent=sp;
  const btns=document.getElementById("sp-annot-popup-btns"); btns.innerHTML="";
  const btnStyle="padding:3px 8px;font-size:11px;border:1px solid #aaa;border-radius:3px;background:#f8f8f8;cursor:pointer;text-align:left;width:100%";
  const mkBtn=(label,cls)=>{
    const b=document.createElement("button"); b.style.cssText=btnStyle; b.textContent=label;
    b.onclick=()=>{ pop.style.display="none"; downloadSpeciesGenes(sp,cls); };
    btns.appendChild(b);
  };
  mkBtn("All classes", undefined);
  sorted.forEach(c=>mkBtn(c, c));
  event.stopPropagation();
  pop.style.display="block";
  const x=Math.min(event.clientX+8,window.innerWidth-pop.offsetWidth-8);
  const y=Math.min(event.clientY+8,window.innerHeight-pop.offsetHeight-8);
  pop.style.left=Math.max(4,x)+"px"; pop.style.top=Math.max(4,y)+"px";
}
document.addEventListener("click",(e)=>{
  const pop=document.getElementById("sp-annot-popup");
  if(pop&&!pop.contains(e.target)&&e.target.id!=="sp-annot-popup") pop.style.display="none";
});

function downloadAnnotations(){
  // Build species × class annotation table from FAMILY_DATA
  const rows=[];
  const spSet=new Set(ALL_SPECIES);
  // gather unique classes
  const classes=[...new Set(FAMILY_DATA.map(d=>d.class||d.pref||"other"))].sort();
  rows.push(["species","total_genes","total_hgs",...classes].join("\t"));
  ALL_SPECIES.forEach(sp=>{
    const m=spMeta[sp]||{genes:0,hgs:0};
    const perClass={};
    classes.forEach(c=>{ perClass[c]=0; });
    FAMILY_DATA.filter(d=>(d.class||d.pref||"other")&&d.species_counts[sp])
               .forEach(d=>{ const c=d.class||d.pref||"other"; perClass[c]=(perClass[c]||0)+d.species_counts[sp]; });
    rows.push([sp,m.genes,m.hgs,...classes.map(c=>perClass[c]||0)].join("\t"));
  });
  const blob=new Blob([rows.join("\n")],{type:"text/tab-separated-values"});
  const a=document.createElement("a");
  a.href=URL.createObjectURL(blob); a.download="species_annotations.tsv";
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(a.href);
}


// ═══════════════════════════════════════════════════════════════════════════════
// TREE VIEW – SELECTION & D3 RENDERING
// ═══════════════════════════════════════════════════════════════════════════════
function selectTree(rec){
  currentIndex=rec;
  currentDetail=loadDetail(rec.id);
  if(!currentDetail){ console.warn("No detail data for",rec.id); return; }
  // reset colour state – default to OG colouring
  hlSet=null; hlQueries=[]; hlGroupIndex=new Map(); hlQueryColors={};
  document.getElementById("hl-search").value="";
  renderHlTags();
  ogHlSet=null; ogHlQueries=[]; ogHlGroupIndex=new Map(); ogHlQueryColors={};
  document.getElementById("og-hl-search").value="";
  renderOgHlTags();
  cladeSp2Color={}; cladeSp2Group={}; cladeGrpColor={}; ogLeaf2Color={}; ogName2Color={}; ogGene2Name={};
  cladeHighlights.clear();
  colorMode="og";
  hmFocusGids=null;

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
  const srcSuffix=rec.source==="generax"?" (GeneRax)":rec.source==="prev"?" (IQ-Tree)":"";
  document.getElementById("tree-title").textContent=rec.id+srcSuffix+" \u00b7 "+rec.n_leaves+" genes";
  document.getElementById("n-ogs-label").textContent=rec.n_ogs+" orthogroups";

  // build OG colour maps (currentDetail already set above)
  const _ogs=currentDetail.ogs||{};
  Object.keys(_ogs).sort().forEach((og,i)=>{
    const col=palette[i%palette.length]; ogName2Color[og]=col;
    for(const gid of _ogs[og]){ ogLeaf2Color[gid]=col; ogGene2Name[gid]=og; }
  });
  renderOgHlTags(); // refresh datalist now that ogName2Color is populated
  drawGeneTree(currentDetail.tree);
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
let rootNode=null, gMain=null, _uid=0, _zoom=null, _zoomScale=1;
let _compareNode1=null;  // first node selected for species comparison
let useBranchLen=false, _phyloScale=1;

function tipFontSVG(){ return (tipFontSize!==null?tipFontSize:11)/_zoomScale; }
function applyTipFontSize(){
  if(!gMain) return;
  const fs=tipFontSVG();
  gMain.selectAll(".leaf-label").attr("font-size",d=>d&&d.data&&d.data.leaf?fs:0);
  gMain.selectAll(".og-label").attr("font-size",fs);
  // update MRCA sub-label tspan inside og-label
  gMain.selectAll(".og-label tspan:nth-child(2)").attr("font-size",Math.max(7,fs*0.82));
  const colR2=Math.max(10,fs*0.9), leafR2=fs*0.36, colCx2=-(colR2-leafR2);
  gMain.selectAll(".count-label").attr("font-size",fs).attr("x",d=>d&&d._children&&!d._isOgCol?colCx2:null);
  gMain.selectAll("circle").filter(d=>d&&d._children&&!d._isOgCol).attr("r",colR2).attr("cx",colCx2);
  gMain.selectAll("circle").filter(d=>d&&!d._children).attr("r",d=>d.data&&d.data.leaf?fs*0.36:isOGNode(d)?fs*0.5:fs*0.26);
  gMain.selectAll(".link").attr("stroke-width",treeLinkWidth/_zoomScale);
}

function isOGNode(d){ return !d.data.leaf && d.data.name && activeOgs()[d.data.name]!==undefined; }

// Species-tree hierarchy for MRCA lookup (built once from SP_TREE_DATA)
const _spHier=(function(){
  if(!SP_TREE_DATA||!SP_TREE_DATA.children) return null;
  return d3.hierarchy(SP_TREE_DATA,d=>d.children||null);
})();

/** Given a Set of species names, return the name of the deepest named MRCA
 *  in the species tree, or null if no named ancestor is found. */
function spMRCAName(speciesSet){
  if(!_spHier||!speciesSet||!speciesSet.size) return null;
  const leaves=_spHier.leaves().filter(l=>speciesSet.has(l.data.name));
  if(!leaves.length) return null;
  // find MRCA by walking paths to root
  const paths=leaves.map(l=>{ const p=[]; let n=l; while(n){p.push(n);n=n.parent;} return p; });
  const sets=paths.map(p=>new Set(p));
  for(const anc of paths[0]){
    if(sets.every(s=>s.has(anc))){
      // anc is the MRCA — return the nearest named node at or above it
      let n=anc;
      while(n){ if(n.data.name&&!n.data.leaf) return n.data.name; n=n.parent; }
      return null;
    }
  }
  return null;
}

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

// ── Species comparison helpers ─────────────────────────────────────────────
function getSpeciesUnder(d){
  const sp=new Set();
  (function walk(n){
    const ch=n.children||n._children;
    if(!ch||!ch.length){ const s=n.data.species||getSpeciesPfx(n.data.gene_id||n.data.name||""); if(s) sp.add(s); return; }
    for(const c of ch) walk(c);
  })(d);
  return sp;
}
function getNodeCompareLabel(d){
  if(d._uid&&customNodeNames[d._uid]) return customNodeNames[d._uid];
  const lbl=d.data._og_label||(isOGNode(d)?d.data.name:"")||d.data.name||"";
  const n=countAllLeaves(d);
  return (lbl||"node")+" ("+n+" tips)";
}
function enterCompareMode(node){
  _compareNode1=node;
  document.getElementById("sp-compare-banner").style.display="flex";
  const svg=document.getElementById("tree-svg"); if(svg) svg.style.cursor="crosshair";
}
function exitCompareMode(){
  _compareNode1=null;
  document.getElementById("sp-compare-banner").style.display="none";
  const svg=document.getElementById("tree-svg"); if(svg) svg.style.cursor="";
}
function runSpeciesComparison(nodeB){
  const nodeA=_compareNode1; exitCompareMode();
  const spA=getSpeciesUnder(nodeA), spB=getSpeciesUnder(nodeB);
  const both=[...spA].filter(s=>spB.has(s)).sort();
  const onlyA=[...spA].filter(s=>!spB.has(s)).sort();
  const onlyB=[...spB].filter(s=>!spA.has(s)).sort();
  const lA=getNodeCompareLabel(nodeA), lB=getNodeCompareLabel(nodeB);
  function spPills(arr,bg,fg){
    if(!arr.length) return `<span style="color:#aaa;font-style:italic;font-size:10px">none</span>`;
    return arr.map(s=>`<span style="display:inline-block;padding:1px 7px;margin:2px 2px 0;border-radius:10px;background:${bg};color:${fg};font-size:10px">${s}</span>`).join("");
  }
  document.getElementById("scp-legend").innerHTML=
    `<span style="color:#2980b9;font-weight:600">[A]</span> ${lA}<br>`+
    `<span style="color:#c0392b;font-weight:600">[B]</span> ${lB}`;
  document.getElementById("scp-content").innerHTML=
    `<div style="margin-bottom:6px"><b style="color:#27ae60">Shared (${both.length})</b><div style="margin-top:3px;line-height:1.8">${spPills(both,"#d4efdf","#1a6b3a")}</div></div>`+
    `<div style="margin-bottom:6px"><b style="color:#2980b9">Only in A (${onlyA.length})</b><div style="margin-top:3px;line-height:1.8">${spPills(onlyA,"#d6eaf8","#1a4a6b")}</div></div>`+
    `<div><b style="color:#c0392b">Only in B (${onlyB.length})</b><div style="margin-top:3px;line-height:1.8">${spPills(onlyB,"#fadbd8","#7b241c")}</div></div>`;
  document.getElementById("sp-compare-panel").style.display="block";
}
document.addEventListener("keydown",ev=>{ if(ev.key==="Escape"&&_compareNode1) exitCompareMode(); });

function collapsedLabel(d){
  const n=countDescLeaves(d._children);
  if(d._uid&&customNodeNames[d._uid]) return customNodeNames[d._uid]+" ["+n+"]";
  // only use d.data.name as a label when it's a verified OG name (not a support value)
  const lbl=d.data._og_label||(isOGNode(d)?d.data.name:"")||"";
  if(lbl) return lbl+" ["+n+"]";
  // manually collapsed: show MRCA name from species tree + count
  const sc={};
  (function cnt(ch){ if(!ch)return; for(const c of ch){
    if(c.data.leaf){const sp=c.data.species||"?";sc[sp]=(sc[sp]||0)+1;}
    else{cnt(c.children);cnt(c._children);}
  }})(d._children);
  const mrca=spMRCAName(new Set(Object.keys(sc)));
  return (mrca||"clade")+" ["+n+"]";
}

/** Full label with species breakdown for collapsed-node tooltip. */
function collapsedTooltip(d){
  const base=collapsedLabel(d);
  const sc={};
  (function cnt(ch){ if(!ch)return; for(const c of ch){
    if(c.data.leaf){const sp=c.data.species||"?";sc[sp]=(sc[sp]||0)+1;}
    else{cnt(c.children);cnt(c._children);}
  }})(d._children);
  const rows=Object.entries(sc).sort((a,b)=>b[1]-a[1])
    .map(([sp,c])=>`<div style="display:flex;gap:6px"><span style="color:${leafColor(sp)}">${sp}</span><span>${c}</span></div>`).join("");
  return `<b>${base}</b>${rows?'<div style="margin-top:4px;font-size:10px">'+rows+'</div>':""}`;
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

/** Annotate expanded internal nodes with _og_label via MRCA from tip OG data.
 *  Used for trees where OG info is in pipe-separated tip labels (not named internals).
 *  Only sets labels on currently-visible (expanded) internal nodes. */
function annotateOGNodes(){
  if(!rootNode) return;
  // If tree already has named OG internal nodes, nothing to annotate
  if(rootNode.descendants().some(d=>isOGNode(d))) return;
  // Clear previous non-collapsed annotations
  rootNode.each(d=>{ if(!d._isOgCol) d.data._og_label=null; });
  // Build set of currently-visible nodes (only those reachable via children, not _children)
  const visibleNodes=new Set(rootNode.descendants());
  // Collect OG → ALL leaf nodes, including inside collapsed (_children) subtrees.
  // This keeps the MRCA stable when subclades are collapsed: the MRCA is computed
  // from the full leaf set, then we only annotate it if it is itself visible.
  const ogGroups={};
  (function walkAll(n){
    if(n.data.leaf){
      const og=n.data.og||ogGene2Name[n.data.gene_id||n.data.name]||"";
      if(og)(ogGroups[og]=ogGroups[og]||[]).push(n);
      return;
    }
    for(const c of (n.children||[])) walkAll(c);
    for(const c of (n._children||[])) walkAll(c);
  })(rootNode);
  for(const [og,leaves] of Object.entries(ogGroups)){
    const mrca=findMRCA(leaves);
    // Only annotate if the MRCA is a visible (expanded), non-leaf, non-collapsed node
    if(mrca&&visibleNodes.has(mrca)&&!mrca.data.leaf&&!mrca._children){ mrca.data._og_label=og; }
  }
}

function toggleOGLabels(){
  showOGLabels=!showOGLabels;
  document.getElementById("btn-og-labels").classList.toggle("active-btn", showOGLabels);
  if(rootNode) renderTree(false);
}

function toggleFocusCollapseStyle(){
  focusCollapseAsTri=!focusCollapseAsTri;
  const btn=document.getElementById("btn-focus-collapse-style");
  btn.classList.toggle("active-btn", focusCollapseAsTri);
  btn.textContent=focusCollapseAsTri?"\u25BC MRCA":"\u25CB Circle";
}

function toggleSupport(){
  showSupport=!showSupport;
  document.getElementById("btn-support").classList.toggle("active-btn", showSupport);
  if(rootNode) renderTree(false);
}

// ── Mini species tree floating panel ──────────────────────────────────────────
function toggleMiniSpPanel(ev){
  const panel=document.getElementById("mini-sp-panel");
  if(panel.style.display==="block"){ panel.style.display="none"; return; }
  // position near button
  const btn=document.getElementById("btn-mini-sp");
  const r=btn.getBoundingClientRect();
  panel.style.top=(r.bottom+6)+"px";
  panel.style.left=Math.max(0,r.left-180)+"px";
  panel.style.display="block";
  drawMiniSpTree();
}
document.addEventListener("click",ev=>{
  const panel=document.getElementById("mini-sp-panel");
  if(panel.style.display==="block"&&!panel.contains(ev.target)&&ev.target.id!=="btn-mini-sp")
    panel.style.display="none";
  // close global heatmap search dropdown when clicking outside
  if(!ev.target.closest("#hm-text-search")&&!ev.target.closest("#hm-search-dd")){
    const dd=document.getElementById("hm-search-dd");
    if(dd) dd.style.display="none";
  }
  // close OG highlight dropdown when clicking outside
  if(!ev.target.closest("#og-hl-search")&&!ev.target.closest("#og-hl-dd")) ogHlHideDD();
});

function drawMiniSpTree(){
  const wrap=document.getElementById("mini-sp-svg-wrap"); wrap.innerHTML="";
  if(!SP_TREE_DATA) return;
  // simple recursive layout
  function flat(n){ return [n].concat(n.children?n.children.flatMap(flat):[]); }
  function clone(n){ return Object.assign({},n,{children:n.children?n.children.map(clone):null}); }
  const tree=clone(SP_TREE_DATA);
  const allN=flat(tree);
  // collect only species present in current tree
  const activeSp=currentIndex?new Set(currentIndex.species):new Set(ALL_SPECIES);
  // assign y: leaves in order
  const leaves=allN.filter(n=>!n.children);
  const leafH=14, leftM=6, W=520;
  leaves.forEach((l,i)=>{ l._y=i*leafH+leafH/2; });
  // propagate y to internals
  function assignY(n){ if(n.children){ n.children.forEach(assignY); n._y=(n.children[0]._y+n.children[n.children.length-1]._y)/2; } }
  assignY(tree);
  // depth
  let maxD=0;
  function assignD(n,d){ n._d=d; maxD=Math.max(maxD,d); if(n.children) n.children.forEach(c=>assignD(c,d+1)); }
  assignD(tree,0);
  const sx=d=>leftM+(d/Math.max(1,maxD))*(W-leftM-100);
  const H=leaves.length*leafH+10;
  const svg=d3.select(wrap).append("svg").attr("width",W).attr("height",H);
  // draw branches
  function drawB(n){ if(!n.children) return; const ys=n.children.map(c=>c._y); svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(n._d)).attr("y1",d3.min(ys)).attr("y2",d3.max(ys)).attr("stroke","#bbb"); n.children.forEach(c=>{ svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(c._d)).attr("y1",c._y).attr("y2",c._y).attr("stroke","#bbb"); drawB(c); }); }
  drawB(tree);
  // leaves
  leaves.forEach(l=>{ const inTree=activeSp.has(l.name); svg.append("circle").attr("cx",sx(maxD)).attr("cy",l._y).attr("r",3).attr("fill",inTree?"#2c3e50":"#ccc"); svg.append("text").attr("x",sx(maxD)+6).attr("y",l._y).attr("dy","0.35em").attr("font-size",9).attr("fill",inTree?"#333":"#bbb").attr("font-family","monospace").text(l.name); });
  // named internal nodes — clickable to add highlight
  function drawInternals(n){ if(!n.children) return;
    if(n.name){
      svg.append("text").attr("class","msp-node-lbl").attr("x",sx(n._d)).attr("y",n._y-5)
        .attr("text-anchor","middle").text(n.name)
        .on("click",()=>{
          // collect all leaves under this clade that are in the active tree
          function sp2(nd){ return nd.children?nd.children.flatMap(sp2):[nd.name]; }
          const cladeSps=sp2(n).filter(s=>activeSp.has(s));
          if(!cladeSps.length) return;
          // try clade name first (works if in CLADE_DATA); else add species individually
          const resolved=resolveQuery(n.name);
          if(resolved.size>0){ addHlTag(n.name); }
          else { cladeSps.forEach(s=>addHlTag(s)); }
          document.getElementById("mini-sp-panel").style.display="none";
        });
    }
    n.children.forEach(drawInternals);
  }
  drawInternals(tree);
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
  _zoomScale=1;
  _zoom=d3.zoom().scaleExtent([0.03,30]).on("zoom",e=>{
    gMain.attr("transform",e.transform);
    _zoomScale=e.transform.k;
    applyTipFontSize();
  });
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
  // Annotate internal nodes with OG labels for pipe-sep tip format (no-op for named-internal trees)
  if(showOGLabels) annotateOGNodes();
  const wrap=document.getElementById("tree-wrap");
  const W=wrap.clientWidth||800, H=wrap.clientHeight||600;
  const mg={top:20,right:240,bottom:36,left:36};
  const iW=W-mg.left-mg.right, iH=H-mg.top-mg.bottom;
  const nVis=rootNode.leaves().length;
  // Count total genes (visible tips + hidden leaves inside collapsed nodes)
  // so that rowH stays proportional even when many nodes are collapsed.
  function countAllLeaves(d){ if(!d.children&&!d._children) return 1; return (d._children||d.children).reduce((s,c)=>s+countAllLeaves(c),0); }
  const nAll=countAllLeaves(rootNode);
  // Minimum row height must accommodate the current tip font size
  const fsEff = tipFontSize !== null ? tipFontSize : 11;
  const minRowH = Math.max(14, Math.ceil(fsEff * 1.3));
  const rowH = Math.max(minRowH, Math.min(Math.max(minRowH, 36), Math.floor(iH/Math.max(nAll,1))));
  const tH=Math.max(iH, nAll*rowH) * treeHeightMult;
  const effRow = rowH * treeHeightMult;  // used here and in triangle points below
  const iWeff = iW * treeWidthMult;
  d3.cluster().size([tH,iWeff])(rootNode);
  treeSvg.attr("width", Math.max(W, mg.left + iWeff + mg.right));

  // Proportional Y-spacing for OG-collapsed nodes.
  // D3 cluster treats each collapsed node as 1 leaf; we redistribute so each
  // collapsed node gets nLeaves * rowH vertical space, preventing overlap.
  {
    const vLeaves = rootNode.leaves(); // includes ._children nodes (D3 sees them as leaves)
    const wts = vLeaves.map(d => d._children ? Math.max(1, countAllLeaves(d) * collapsedFraction) : 1);
    const anyCollapsed = wts.some(w => w > 1);
    if (anyCollapsed) {
      let cy = effRow / 2;
      vLeaves.forEach((d, i) => {
        d.x = cy + wts[i] * effRow / 2;
        cy += wts[i] * effRow;
      });
      // fix internal node positions as midpoint of children (post-order)
      rootNode.eachAfter(d => {
        if (d.children && d.children.length)
          d.x = (d.children[0].x + d.children[d.children.length - 1].x) / 2;
      });
    }
  }

  if(useBranchLen){
    assignBranchLenPos(rootNode,0);
    const maxBL=d3.max(rootNode.leaves(),d=>d._by)||1;
    _phyloScale=iW/maxBL;
  }

  const dur=animate?240:0;

  // ── Clade highlight backgrounds ──
  {
    let hlLayer=gMain.select(".clade-hl-layer");
    if(hlLayer.empty()) hlLayer=gMain.insert("g",":first-child").attr("class","clade-hl-layer");
    hlLayer.lower();
    hlLayer.selectAll("*").remove();
    if(cladeHighlights.size){
      const allNodes=rootNode.descendants();
      // x2 stops just past the rightmost leaf nodes (20px padding); the label margin is NOT included
      const x2=nodeX(rootNode,mg)+iWeff+20;
      cladeHighlights.forEach((rec,uid)=>{
        const hn=allNodes.find(d=>d._uid===uid);
        if(!hn) return;
        const sub=hn.descendants(); // visible subtree only (follows d.children)
        const ys=sub.map(d=>d.x);
        const pad=effRow*0.55;
        const y1=Math.min(...ys)-pad+mg.top;
        const y2=Math.max(...ys)+pad+mg.top;
        const x1=nodeX(hn,mg)-6;
        const color=rec.color||"#ffe066";
        const label=rec.label||"";
        hlLayer.append("rect")
          .attr("x",x1).attr("y",y1)
          .attr("width",Math.max(0,x2-x1)).attr("height",Math.max(0,y2-y1))
          .attr("fill",color).attr("opacity",cladeHlAlpha).attr("rx",4)
          .style("cursor","pointer")
          .on("contextmenu",(ev)=>showCladeHlPopup(ev,uid))
          .on("dblclick",(ev)=>{
            ev.stopPropagation();
            openColorPicker(color,c=>{ rec.color=c; renderTree(false); });
          });
        if(label){
          const labelFontSize=Math.max(11, Math.min(22, (y2-y1)*0.38));
          hlLayer.append("text")
            .attr("x", x2-8)
            .attr("y", (y1+y2)/2)
            .attr("text-anchor","middle")
            .attr("dominant-baseline","central")
            .attr("transform",`rotate(-90,${x2-8},${(y1+y2)/2})`)
            .attr("font-size", labelFontSize)
            .attr("font-weight","700")
            .attr("fill", color)
            .attr("opacity", Math.min(1, cladeHlAlpha*3.5))
            .style("pointer-events","none")
            .text(label);
        }
      });
    }
  }

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
    .attr("d",d=>elbowPath(d.source,d.target,mg,3))
    .attr("stroke-width",treeLinkWidth/_zoomScale);

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
        g.append("polygon").attr("class","col-tri");
        g.append("text").attr("class","count-label")  // manual-collapse count
          .attr("dy","0.35em").attr("text-anchor","middle").attr("pointer-events","none");
        g.append("text").attr("class","leaf-label");
        g.append("text").attr("class","og-label");
        g.append("text").attr("class","support-lbl");
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
    .attr("cx",0)   // reset any offset from manual-collapse state
    .attr("r",d=>{
      if(d.data.leaf) return tipFontSVG()*0.36;
      const isOGlbl=showOGLabels&&(isOGNode(d)||d.data._og_label);
      if(isOGlbl) return tipFontSVG()*0.55;
      if(isOGNode(d)) return tipFontSVG()*0.5;
      return tipFontSVG()*0.26;
    })
    .attr("fill",d=>{
      if(d.data.leaf){
        const gid2=d.data.gene_id||d.data.name;
        const og2=d.data.og||ogGene2Name[gid2]||"";
        if(ogHlSet!==null){
          if(!ogHlSet.has(og2)) return "#e8e8e8";
          const gi=ogHlGroupIndex.get(og2)??0;
          return ogHlTagColor(gi);
        }
        if(colorMode==="og") return ogLeafColor(gid2, d.data.species);
        const c=leafColor(d.data.species||"");
        return hlSet&&!hlSet.has(d.data.species||"")?"#e8e8e8":c;
      }
      if(showOGLabels&&(isOGNode(d)||d.data._og_label)) return "#222";
      if(isOGNode(d))return "#e74c3c";
      return "#ccc";
    })
    .attr("stroke",d=>{
      if(d.data.leaf) return "none";
      if(showOGLabels&&(isOGNode(d)||d.data._og_label)) return "#000";
      if(isOGNode(d)) return "#b03a2e";
      return "#aaa";
    })
    .attr("stroke-width",d=>{
      if(d.data.leaf) return null;
      if(showOGLabels&&(isOGNode(d)||d.data._og_label)) return 2;
      if(isOGNode(d)) return 1.8;
      return 0.8;
    })
    .on("click",(event,d)=>{
      event.stopPropagation();
      // compare mode intercepts any node click as the second target; same-node click is a no-op
      if(_compareNode1!==null){ if(_compareNode1!==d) runSpeciesComparison(d); hideTip(); return; }
      if(d.data.leaf){
        const tipName=d.data.gene_id||d.data.name||"";
        if(tipName&&navigator.clipboard) navigator.clipboard.writeText(tipName).then(()=>{
          const tt=document.getElementById("tooltip");
          tt.style.display="block"; tt.innerHTML="<span style='color:#1abc9c'>&#10003; Copied!</span>";
          tt.style.left=Math.min(event.clientX+12,window.innerWidth-tt.offsetWidth-5)+"px";
          tt.style.top=Math.min(event.clientY+12,window.innerHeight-tt.offsetHeight-5)+"px";
          setTimeout(()=>{ tt.style.display="none"; },900);
        });
        return;
      }
      if(d._children){ d.children=d._children; d._children=null; d._isOgCol=false; renderTree(true); hideTip(); }
      else if(d.children){ showCollapseChoicePopup(event,d); }
    })
    .on("mouseover",(event,d)=>{
      if(!d.data.leaf&&d.children) showTip(event,'<b>'+(d.data.name||"internal")+'</b><div style="font-size:9px;color:#aaa;margin-top:3px">click to choose collapse style</div>');
      else showTip(event,d);
    })
    .on("mousemove",moveTip).on("mouseout",hideTip);

  // OG-collapse triangle (only when _isOgCol)
  nodeSel.select(".col-tri")
    .attr("display",d=>(d._children&&d._isOgCol)?null:"none")
    .attr("points",d=>{
      if(!d._children||!d._isOgCol) return "";
      const nL=countAllLeaves(d);
      // Half-height uses effRow (= rowH * treeHeightMult) so the triangle fills
      // exactly the proportional space assigned to this node; capped at 48% to
      // leave a 2px gap on each side between adjacent triangles.
      const halfH=Math.max(rowH*0.45, nL*effRow*collapsedFraction*0.48);
      return `0,0 ${BADGE_W},${-halfH} ${BADGE_W},${halfH}`;
    })
    .style("fill",d=>{
      if(!d._children||!d._isOgCol) return null;
      // per-node user color override takes priority
      if(nodeTriColors.has(d._uid)) return nodeTriColors.get(d._uid);
      const ogName=d.data._og_label||d.data.name||"";
      if(ogHlSet!==null){
        if(ogHlSet.has(ogName)){ const gi=ogHlGroupIndex.get(ogName)??0; return ogHlTagColor(gi)+"bb"; }
        return "#e8e8e8";
      }
      if(colorMode==="og"){
        const col=ogName2Color[ogName];
        return col ? col+"55" : "#ffffff";
      }
      return "#ffffff";
    })
    .attr("stroke-width",d=>(d._children&&d._isOgCol)?treeLinkWidth/_zoomScale:null)
    .on("click",(event,d)=>{
      if(!d._children) return; event.stopPropagation();
      if(_compareNode1!==null){ if(_compareNode1!==d) runSpeciesComparison(d); hideTip(); return; }
      showTriActionPopup(event,d);
    })
    .on("mouseover",(ev,d)=>{
      if(d._children) showTip(ev,collapsedTooltip(d)+'<div style="font-size:9px;color:#aaa;margin-top:4px">click for options</div>');
    }).on("mousemove",moveTip).on("mouseout",hideTip);

  // manual-collapse: larger circle + leaf count label
  // Shift circle leftward so its RIGHT edge aligns with the normal leaf column,
  // preventing it from protruding into neighbouring nodes' label zones.
  {
    const leafR = tipFontSVG()*0.36;
    const colR  = d => Math.max(10, tipFontSVG()*0.9);
    const colCx = d => { const r=colR(d); return -(r-leafR); };
    nodeSel.select("circle")
      .filter(d=>d._children&&!d._isOgCol)
      .attr("r",colR).attr("cx",colCx)
      .attr("fill","#f5f5f5").attr("stroke","#999").attr("stroke-width",1.2)
      .attr("display",null);
    nodeSel.select(".count-label")
      .attr("display",d=>(d._children&&!d._isOgCol)?null:"none")
      .attr("font-size",tipFontSVG())
      .attr("x",d=>(d._children&&!d._isOgCol)?colCx(d):0)
      .text(d=>(d._children&&!d._isOgCol)?countAllLeaves(d):"");
  }

  // leaf labels: gene_id + OG name
  nodeSel.select(".leaf-label")
    .attr("x",7).attr("dy","0.32em").attr("text-anchor","start")
    .attr("font-size",d=>d.data.leaf?tipFontSVG():0)
    .style("cursor",()=>_compareNode1?"crosshair":"default")
    .on("click",(event,d)=>{
      // In compare mode, leaf labels also act as a click target for the second node
      if(_compareNode1!==null){ event.stopPropagation(); if(_compareNode1!==d) runSpeciesComparison(d); hideTip(); }
    })
    .attr("display",d=>{
      if(!d.data.leaf) return "none";
      if(hideNonHl){
        const gid3=d.data.gene_id||d.data.name;
        if(hmFocusGids!==null){
          if(!hmFocusGids.has(gid3)) return "none";
        } else {
          const og3=d.data.og||ogGene2Name[gid3]||"";
          if(hlSet!==null&&!hlSet.has(d.data.species||"")) return "none";
          if(ogHlSet!==null&&!ogHlSet.has(og3)) return "none";
        }
      }
      return null;
    })
    .text("")
    .each(function(d){
      if(!d.data.leaf) return;
      const el=d3.select(this);
      el.selectAll("tspan").remove();
      const gid=d.data.gene_id||d.data.name;
      const og=d.data.og||ogGene2Name[gid]||"";
      const ref=d.data.ref||"";
      // dim if either species-hl or OG-hl is active and this tip doesn't match
      const notSpHl=hlSet!==null&&!hlSet.has(d.data.species||"");
      const ogHlActive=ogHlSet!==null;
      const inOgHl=ogHlActive&&ogHlSet.has(og);
      const notHl=notSpHl||(ogHlActive&&!inOgHl);
      // if OG-hl active, use that group's color; otherwise species color
      let baseCol;
      if(ogHlActive&&inOgHl){
        const gi=ogHlGroupIndex.get(og)??0;
        baseCol=ogHlTagColor(gi);
      } else {
        baseCol=notHl?"#ccc":(colorMode==="og"?ogLeafColor(d.data.gene_id||d.data.name,d.data.species):leafColor(d.data.species||""));
      }
      const sepCol=notHl?"#ddd":"#bbb";
      const ogCol=notHl?"#ccc":(inOgHl?baseCol:"#4a7aad");
      const refCol=notHl?"#ccc":"#2e8b57";
      let first=true;
      function sep(){ if(!first) el.append("tspan").attr("fill",sepCol).text(" \u00b7 "); first=false; }
      if(showGeneId){ sep(); el.append("tspan").attr("fill",baseCol).text(gid); }
      if(showOGName&&og){ sep(); el.append("tspan").attr("fill",ogCol).text(og); }
      if(showRefOrtho&&ref){ sep(); el.append("tspan").attr("fill",refCol).text(ref); }
    });

  // OG labels: beside OG-collapsed triangle, or beside expanded OG-named internal
  nodeSel.select(".og-label")
    .attr("x",d=>d._children?BADGE_W+6:-7)
    .attr("dy","0.35em")
    .attr("text-anchor",d=>d._children?"start":"end")
    .attr("font-size",tipFontSVG())
    .style("cursor",d=>{
      if(_compareNode1) return "crosshair";
      if(d._isOgCol) return "pointer";
      if(!d._children&&(isOGNode(d)||d.data._og_label)) return "pointer";
      return "default";
    })
    .attr("fill",d=>{
      // OG nodes (named or annotated) → red; all collapsed non-OG nodes → grey
      const ogName=d._isOgCol?(d.data._og_label||d.data.name||""):(isOGNode(d)?d.data.name:(d.data._og_label||""));
      if(ogName&&ogHlSet!==null){
        if(ogHlSet.has(ogName)){ const gi=ogHlGroupIndex.get(ogName)??0; return ogHlTagColor(gi); }
        return "#ccc";
      }
      if(isOGNode(d)||d.data._og_label) return "#b5371f";
      return "#444";
    })
    .attr("display",d=>(!d.data.leaf&&(d._isOgCol||(showOGLabels&&!d._children&&(isOGNode(d)||d.data._og_label))))?null:"none")
    .on("click",(event,d)=>{
      event.stopPropagation();
      if(_compareNode1!==null){ if(_compareNode1!==d) runSpeciesComparison(d); hideTip(); return; }
      // Triangle popup for collapsed OG; for expanded OG/annotated nodes toggle OG highlight
      if(d._children&&d._isOgCol) showTriActionPopup(event,d);
      else if(!d._children&&(isOGNode(d)||d.data._og_label)){
        const ogName=d.data._og_label||(isOGNode(d)?d.data.name:"");
        if(ogName){ const idx=ogHlQueries.indexOf(ogName); if(idx>=0) removeOgHlTag(idx); else addOgHlTag(ogName); }
      }
    })
    .text("")
    .each(function(d){
      const el=d3.select(this);
      el.selectAll("tspan").remove();
      let mainLbl;
      if(d._isOgCol){
        mainLbl=collapsedLabel(d);
      } else {
        mainLbl=isOGNode(d)?d.data.name:(d.data._og_label||"");
      }
      if(!mainLbl) return;
      el.append("tspan").text(mainLbl);
      // For OG-collapsed nodes (actual OG or annotated), append MRCA clade name on a second line
      // Skip for plain manual collapses: their label is already the MRCA name
      if(d._isOgCol&&(isOGNode(d)||d.data._og_label)){
        const sc={};
        (function cnt(ch){ if(!ch)return; for(const c of ch){
          if(c.data.leaf){const sp=c.data.species||"?";sc[sp]=(sc[sp]||0)+1;}
          else{cnt(c.children);cnt(c._children);}
        }})(d._children);
        const mrca=spMRCAName(new Set(Object.keys(sc)));
        if(mrca){
          el.append("tspan")
            .attr("x",BADGE_W+6).attr("dy","1.2em")
            .attr("font-size",Math.max(7,tipFontSVG()*0.82))
            .attr("fill","#888")
            .text(mrca);
        }
      }
    });

  // support labels: shown on internal nodes when showSupport is true
  nodeSel.select(".support-lbl")
    .attr("x",-4).attr("dy","-0.4em").attr("text-anchor","end")
    .attr("font-size",Math.max(7,tipFontSVG()*0.75))
    .attr("fill","#888").attr("pointer-events","none")
    .attr("display",d=>(showSupport&&!d.data.leaf&&d.data.support!=null&&!d._children)?null:"none")
    .text(d=>(d.data.support!=null)?d.data.support:"");
  applyTipFontSize();
}

// ── tree controls ──
function expandAll(){
  if(!rootNode)return;
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;d._isOgCol=false;}});
  renderTree(true);
}
function collapseToOGs(){
  if(!rootNode)return;
  // pass 1: expand all, clear flags
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;} d._isOgCol=false;});
  // pass 2a: collapse named OG internal nodes (POSSVM trees with annotated internals)
  let found=false;
  rootNode.each(d=>{
    if(!d.data.leaf&&isOGNode(d)&&d.children){d._children=d.children;d.children=null;d._isOgCol=true;found=true;}
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
        mrca.data._og_label=og; mrca._isOgCol=true;
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
  if((!hmFocusGids&&!hlSet&&!ogHlSet)||!rootNode)return;
  // expand all
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;}});
  // mark nodes that have ≥1 matching descendant (post-order)
  const hasHl=new Map();
  rootNode.eachAfter(d=>{
    if(d.data.leaf){
      let match;
      if(hmFocusGids){
        match=hmFocusGids.has(d.data.gene_id||d.data.name);
      } else {
        const sp=d.data.species||"";
        const gid=d.data.gene_id||d.data.name;
        const og=d.data.og||ogGene2Name[gid]||"";
        match=(!hlSet||hlSet.has(sp))&&(!ogHlSet||ogHlSet.has(og));
      }
      hasHl.set(d,match);
    } else {
      hasHl.set(d,(d.children||[]).some(c=>hasHl.get(c)));
    }
  });
  // collapse subtrees with no matching leaf
  // use triangle (MRCA) or circle depending on focusCollapseAsTri
  rootNode.each(d=>{
    if(!d.data.leaf&&d.children&&!hasHl.get(d)){
      d._children=d.children; d.children=null;
      d._isOgCol=focusCollapseAsTri;
    }
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

document.getElementById("btn-og-labels").classList.toggle("active-btn", showOGLabels);

if (hasHeatmapData || TREE_INDEX.length > 0) {
  switchTab("sptree");
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
    family_details = load_family_details(args.family_info)
    family_records = build_family_records(search_dir, family_info)
    hg_records     = build_hg_records(cluster_dir, family_info)

    # Species tree (for ordering + cladogram)
    species_order: list = []
    tree_dict: dict = {}
    newick_raw: str = ""
    if args.species_tree and Path(args.species_tree).exists():
        species_order, tree_dict = load_tree_data(args.species_tree)
        newick_raw = Path(args.species_tree).read_text().strip()

    # Gene tree data (POSSVM)
    if not possvm_dir.exists():
        print(f"WARN: {possvm_dir} does not exist – no gene trees.", file=sys.stderr)
    records, all_species = load_possvm_trees(possvm_dir, source="generax") if possvm_dir.exists() else ([], [])
    print(f"Loaded {len(records)} gene trees, {len(all_species)} species.", file=sys.stderr)

    # Prev trees (IQ-TREE2 original, pre-GeneRax) — optional
    prev_records: dict = {}
    if args.possvm_prev_dir:
        prev_dir = Path(args.possvm_prev_dir)
        if prev_dir.is_dir():
            prev_list, prev_sp = load_possvm_trees(prev_dir, source="prev")
            prev_records = {r["id"]: r for r in prev_list}
            all_species = sorted(set(all_species) | set(prev_sp))
            # Include HGs that have a prev tree but no GeneRax output
            generax_ids = {r["id"] for r in records}
            for r in prev_list:
                if r["id"] not in generax_ids:
                    records.append(r)
            print(f"Loaded {len(prev_records)} prev gene trees (original IQ-TREE2).",
                  file=sys.stderr)
    # Filter all data to families present in genefam.csv (when provided)
    if family_info:
        before = len(records), len(family_records), len(hg_records)
        records        = [r for r in records        if r["prefix"] in family_info or r["family"] in family_info]
        family_records = [r for r in family_records if r["family"]  in family_info]
        hg_records     = [r for r in hg_records     if r["family"]  in family_info]
        print(
            f"Filtered to genefam families: "
            f"{before[0]}→{len(records)} trees, "
            f"{before[1]}→{len(family_records)} families, "
            f"{before[2]}→{len(hg_records)} HGs.",
            file=sys.stderr,
        )
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
        idx["source"]   = rec.get("source", "generax")
        fam_key = rec["family"] if rec["family"] in family_info else rec["prefix"]
        idx["class"] = family_info.get(fam_key, rec.get("prefix", ""))
        index_records.append(idx)

    # Build Families tab data
    fam_hg_counts    = defaultdict(int)
    fam_gene_counts  = defaultdict(int)
    fam_species_sets = defaultdict(set)
    for r in hg_records:
        fam_hg_counts[r["family"]] += 1
    for r in family_records:
        fam_gene_counts[r["family"]]  = r.get("total", 0)
        fam_species_sets[r["family"]] = set(r.get("species_counts", {}).keys())

    # Count HGs-with-trees per family
    fam_generax_counts = defaultdict(int)   # GeneRax trees
    fam_tree_hg_ids    = defaultdict(set)   # union of all tree HG ids
    for r in records:                        # GeneRax (primary) trees
        fam_generax_counts[r["family"]] += 1
        fam_tree_hg_ids[r["family"]].add(r["id"])
    for r in prev_records.values():          # IQ-Tree (prev) trees
        fam_tree_hg_ids[r["family"]].add(r["id"])
    fam_tree_counts = {fam: len(ids) for fam, ids in fam_tree_hg_ids.items()}
    have_generax = bool(records)             # flag exposed to JS

    family_info_records = []
    all_families = sorted(set(family_details.keys()) | set(fam_hg_counts.keys()) | set(fam_gene_counts.keys()))
    for fam in all_families:
        det = family_details.get(fam, {})
        family_info_records.append({
            "family":    fam,
            "pfam":      det.get("pfam", []),
            "category":  det.get("category", ""),
            "cls":       det.get("cls", ""),
            "n_hgs":     fam_hg_counts.get(fam, 0),
            "total":     fam_gene_counts.get(fam, 0),
            "n_species": len(fam_species_sets.get(fam, set())),
            "n_trees":   fam_tree_counts.get(fam, 0),
            "n_generax": fam_generax_counts.get(fam, 0),
        })

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
            .replace("%%LAZY_SCRIPTS%%",       lazy_scripts)
            .replace("%%SPECIES_ORDER%%",      json.dumps(species_order))
            .replace("%%TREE_DATA%%",          json.dumps(tree_dict))
            .replace("%%FAMILY_DATA%%",        json.dumps(family_records))
            .replace("%%HG_DATA%%",            json.dumps(hg_records))
            .replace("%%TREE_INDEX_JSON%%",    json.dumps(index_records))
            .replace("%%SPECIES_JSON%%",       json.dumps(all_species))
            .replace("%%CLADE_DATA_JSON%%",    json.dumps(clade_groupings))
            .replace("%%NEWICK_RAW%%",         json.dumps(newick_raw))
            .replace("%%FAMILY_INFO_JSON%%",   json.dumps(family_info_records))
            .replace("%%HAVE_GENERAX_JSON%%",  json.dumps(have_generax)))

    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
