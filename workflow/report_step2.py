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
from typing import Optional


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

def _iter_family_info_rows(csv_path):
    """Yield normalised rows from genefam.csv or gene_families_searchinfo.csv."""
    try:
        with open(csv_path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 7:
                    continue
                family = cols[0].strip()
                if not family:
                    continue
                yield {
                    "family": family,
                    "pfam_raw": cols[1].strip(),
                    "category": cols[5].strip(),
                    "cls": cols[6].strip(),
                }
    except (FileNotFoundError, OSError):
        return


def load_family_info(csv_path) -> dict:
    """Return {family_or_alias: class_label} from supported family metadata TSVs."""
    info = {}
    for row in _iter_family_info_rows(csv_path):
        info[row["family"]] = row["cls"]
        if row["cls"]:
            info.setdefault(row["cls"], row["cls"])
    return info


def load_family_details(csv_path) -> dict:
    """Return {family_name: {pfam: [str,...], category: str, cls: str}} from supported family metadata TSVs."""
    details = {}
    for row in _iter_family_info_rows(csv_path):
        pfam_raw = row["pfam_raw"]
        pfam_list = [p.strip() for p in pfam_raw.split(",") if p.strip()] if pfam_raw else []
        details[row["family"]] = {
            "pfam": pfam_list,
            "category": row["category"],
            "cls": row["cls"],
        }
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


def parse_fasta_genes(path) -> dict:
    """Return {species: [gene_id, ...]} from FASTA headers."""
    genes: dict = defaultdict(list)
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    gene = line[1:].strip().split()[0]
                    genes[get_species_prefix(gene)].append(gene)
    except (FileNotFoundError, OSError):
        pass
    return dict(genes)


def load_domain_hits(search_dir: Path) -> dict:
    """Return {gene_id: [{name,start,end}, ...]} from *.domains.csv files."""
    hits: dict = defaultdict(list)
    if not search_dir.is_dir():
        return {}
    for domains_file in sorted(search_dir.glob("*.domains.csv")):
        try:
            with open(domains_file) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 4:
                        continue
                    gene_id = cols[0].strip()
                    dom_name = cols[3].strip()
                    if not gene_id or not dom_name:
                        continue
                    try:
                        start = int(float(cols[1]))
                        end = int(float(cols[2]))
                    except ValueError:
                        continue
                    hits[gene_id].append({
                        "name": dom_name,
                        "start": start,
                        "end": end,
                    })
        except OSError:
            continue
    for gene_id, gene_hits in hits.items():
        gene_hits.sort(key=lambda rec: (rec["start"], rec["end"], rec["name"]))
    return dict(hits)


def load_reference_names(refnames_path: Optional[str], refsps: Optional[str] = None) -> dict:
    """Return {gene_id: reference_name} from a POSSVM refnames table."""
    if not refnames_path:
        return {}
    path = Path(refnames_path)
    if not path.exists():
        return {}
    ref_species = {s.strip() for s in refsps.split(",") if s.strip()} if refsps else None
    ref_map: dict = {}
    try:
        with open(path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 2:
                    continue
                gene_id = cols[0].strip()
                ref_name = cols[1].strip()
                if not gene_id or not ref_name:
                    continue
                if ref_species and get_species_prefix(gene_id) not in ref_species:
                    continue
                ref_map[gene_id] = ref_name
    except OSError:
        return {}
    return ref_map


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
        annotated = family in family_info
        records.append({
            "id": f"{pref}.{family}",
            "family": family,
            "pref": pref,
            "class": cls,
            "annotated": annotated,
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
        annotated = family in family_info
        records.append({
            "id": stem,
            "family": family,
            "pref": pref,
            "hg": hg_id,
            "class": cls,
            "annotated": annotated,
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


def load_og_csv(csv_path: Path) -> tuple[dict, dict]:
    """Return ({og_name: [gene_id, ...]}, {gene_id: metadata}) from a POSSVM CSV."""
    og_members: dict = defaultdict(list)
    gene_meta: dict = {}
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
                if gene.lower() == "gene":
                    continue
                og_members[og].append(gene)
                meta = gene_meta.setdefault(gene, {})
                if len(parts) >= 3 and parts[2] and parts[2].upper() != "NA":
                    meta["og_support"] = parts[2]
                if len(parts) >= 4 and parts[3] and parts[3].upper() != "NA":
                    meta["ref_ortholog"] = parts[3]
                if len(parts) >= 5 and parts[4] and parts[4].upper() != "NA":
                    meta["ref_support"] = parts[4]
    except OSError:
        pass
    return dict(og_members), gene_meta


def load_gene_lengths(cluster_dir: Path) -> dict:
    """Return {gene_id: protein_length} from cluster FASTAs."""
    lengths: dict = {}
    if not cluster_dir.is_dir():
        return lengths
    for fasta_file in sorted(cluster_dir.glob("*.fasta")):
        try:
            with open(fasta_file) as fh:
                gene_id = None
                seq_chunks = []
                for line in fh:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if line.startswith(">"):
                        if gene_id is not None:
                            lengths.setdefault(gene_id, len("".join(seq_chunks)))
                        gene_id = line[1:].strip().split()[0]
                        seq_chunks = []
                    else:
                        seq_chunks.append(line.strip())
                if gene_id is not None:
                    lengths.setdefault(gene_id, len("".join(seq_chunks)))
        except OSError:
            continue
    return lengths


def load_possvm_trees(possvm_dir: Path, source: str = "generax") -> tuple[list, list, dict]:
    """Return (tree_records, all_species, gene_meta).  source='generax' or 'prev'."""
    try:
        from ete3 import Tree  # type: ignore
    except ImportError:
        print("ERROR: ete3 not installed. Run: pip install ete3", file=sys.stderr)
        return [], [], {}

    records = []
    all_species: set = set()
    gene_meta: dict = {}

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
        ogs, csv_gene_meta = load_og_csv(csv_path) if csv_path.exists() else ({}, {})
        for gene_id, meta in csv_gene_meta.items():
            gene_meta.setdefault(gene_id, {}).update(meta)

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

    return records, sorted(all_species), gene_meta


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
#main{flex:1;display:flex;overflow:hidden;min-width:0}
#controls{flex:0 0 320px;overflow:auto;padding:10px;background:#f5f5f5;border-right:1px solid #ddd;display:flex;flex-direction:column;gap:8px}
.ctrl-topline{display:flex;flex-direction:column;align-items:stretch;gap:8px}
.ctrl-primary{display:flex;align-items:center;gap:8px;flex-wrap:wrap;min-width:0}
.ctrl-drawers{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:8px}
.ctrl-panel{border:1px solid #d6dde3;border-radius:10px;background:#fff;overflow:hidden;box-shadow:0 1px 2px rgba(44,62,80,.04)}
.ctrl-panel > summary{list-style:none;cursor:pointer;padding:8px 10px;font-size:11px;font-weight:700;color:#415768;display:flex;align-items:center;justify-content:space-between;gap:10px;user-select:none}
.ctrl-panel > summary::-webkit-details-marker{display:none}
.ctrl-panel > summary::after{content:"+";font-size:14px;line-height:1;color:#7d93a3}
.ctrl-panel[open] > summary{border-bottom:1px solid #edf1f4;background:#fbfcfd}
.ctrl-panel[open] > summary::after{content:"−"}
.ctrl-panel-body{padding:10px;display:flex;flex-wrap:wrap;align-items:center;gap:8px}
.ctrl-panel-body.vert{flex-direction:column;align-items:stretch}
.ctrl-field{font-size:11px;color:#555;display:flex;align-items:center;gap:4px;flex-wrap:wrap}
.ctrl-inline-group{display:flex;flex-wrap:wrap;align-items:center;gap:6px}
.ctrl-subtle{font-size:10px;color:#8b98a3}
.ctrl-grp-label{font-size:9px;font-weight:700;color:#aaa;letter-spacing:.05em;text-transform:uppercase;white-space:nowrap;cursor:default;user-select:none}
.ctrl-sep{width:1px;height:16px;background:#ddd;align-self:center;flex-shrink:0}
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
@media (max-width: 1280px){
  #controls{flex-basis:280px}
}
@media (max-width: 980px){
  #main{flex-direction:column}
  #controls{flex:0 0 auto;max-height:42vh;border-right:none;border-bottom:1px solid #ddd}
  .ctrl-drawers{grid-template-columns:1fr}
}

/* tree svg */
#tree-wrap{flex:1;min-width:0;overflow:hidden;position:relative;background:#fff}
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
          <div class="ctrl-topline">
            <div class="ctrl-primary">
              <button class="ctrl-btn" onclick="expandAll()" title="Expand all collapsed nodes">Expand all</button>
              <button class="ctrl-btn" onclick="collapseAll()" title="Collapse all nodes to leaves">Collapse all</button>
              <button class="ctrl-btn" id="btn-reset-root" onclick="resetRoot()" title="Restore the original root (undo all rerooting)" style="display:none">&#x21BA; Reset root</button>
              <button class="ctrl-btn" id="btn-reset-focus" onclick="resetFocus()" title="Restore the full tree after focusing on a subtree" style="display:none">&#x21F1; Reset focus</button>
              <button class="ctrl-btn" onclick="fitTree()" title="Fit the tree to the current panel size">&#x2922; Fit</button>
              <button class="ctrl-btn" id="tree-toggle" style="display:none;background:#e8f0fe;border-color:#4a90d9" onclick="toggleTreeSource()" title="Switch between available tree sources (e.g. GeneRax vs FastTree)">Showing: GeneRax</button>
              <button class="ctrl-btn" id="btn-collapse-ogs" onclick="toggleCollapseToOGs()" title="Collapse each orthogroup clade into a single labelled triangle">Collapse to OGs</button>
              <button class="ctrl-btn" id="btn-highlight-ogs" onclick="toggleHighlightOGs()" title="Shade the background of each orthogroup clade with its colour">Highlight OGs</button>
              <button class="ctrl-btn" id="btn-possvm" onclick="togglePossvmPanel(event)" title="Re-call orthogroups using the POSSVM species-overlap method — choose ingroup species, then click Run">Orthogroups</button>
              <button class="ctrl-btn" onclick="downloadTreeSVG()" title="Download the current tree view as an SVG file">&#11015; SVG</button>
              <button class="ctrl-btn" onclick="downloadTreePNG()" title="Download the current tree view as a PNG image">&#11015; PNG</button>
            </div>
            <div class="ctrl-primary">
              <span id="tree-title"></span>
              <span id="n-ogs-label"></span>
            </div>
          </div>

          <div class="ctrl-drawers">
            <details class="ctrl-panel" open>
              <summary>Labels &amp; Display</summary>
              <div class="ctrl-panel-body">
                <button class="ctrl-btn" id="btn-og-labels" onclick="toggleOGLabels()" title="Show OG name at the MRCA (deepest shared ancestor) of each orthogroup">OG labels</button>
                <button class="ctrl-btn" id="btn-support" onclick="toggleSupport()" title="Show bootstrap / SH-aLRT support values at internal nodes">Support</button>
                <button class="ctrl-btn" id="btn-ds-nodes" onclick="toggleDSNodes()" title="Show duplication (D) and speciation (S) node types — requires POSSVM to be run first">D/S nodes</button>
                <button class="ctrl-btn" id="btn-lengths" onclick="toggleLengths()" title="Draw branches proportional to evolutionary distance instead of a cladogram">Branch lengths</button>
                <span class="ctrl-field">
                  Tip:
                  <label style="display:flex;align-items:center;gap:2px;cursor:pointer" title="Show gene identifier in tip labels"><input type="checkbox" id="chk-geneid"> ID</label>
                  <label style="display:flex;align-items:center;gap:2px;cursor:pointer" title="Show orthogroup name in tip labels"><input type="checkbox" id="chk-og"> OG</label>
                  <label style="display:flex;align-items:center;gap:2px;cursor:pointer" title="Show reference orthologue in tip labels"><input type="checkbox" id="chk-ref"> ref</label>
                  <label style="display:flex;align-items:center;gap:2px;cursor:pointer" title="Hide tip labels for species not in the current highlight set"><input type="checkbox" id="chk-hide-nonhl"> hide non-hl</label>
                </span>
                <div class="ctrl-inline-group">
                  <span class="ctrl-subtle">Non-highlighted branches</span>
                  <button id="btn-focus-collapse-style" class="ctrl-btn active-btn" onclick="toggleFocusCollapseStyle()" title="Toggle how non-highlighted branches are collapsed: MRCA triangle or single circle node" style="padding:1px 6px;font-size:10px">&#9660; MRCA</button>
                </div>
              </div>
            </details>

            <details class="ctrl-panel">
              <summary>Style</summary>
              <div class="ctrl-panel-body">
                <label class="ctrl-field" title="Colour leaves by species, by orthogroup, or by a predefined clade grouping">Color:
                  <select id="color-by" style="font-size:11px;padding:2px 4px;border:1px solid #bbb;border-radius:3px">
                    <option value="species">by species</option>
                  </select>
                </label>
                <label class="ctrl-field" title="Font size of tip labels">
                  Label:
                  <input type="range" id="tip-font-slider" min="6" max="24" step="1" value="11" style="width:60px;cursor:pointer;accent-color:#4a90d9">
                  <span id="tip-font-val" style="width:20px;text-align:right">11</span>px
                </label>
                <label class="ctrl-field" title="Thickness of tree branches">
                  Lines:
                  <input type="range" id="line-width-slider" min="1" max="8" step="0.5" value="1.3" style="width:60px;cursor:pointer;accent-color:#4a90d9">
                  <span id="line-width-val" style="width:20px;text-align:right">1.3</span>px
                </label>
                <label class="ctrl-field" title="Vertical spacing multiplier — increase to spread leaves further apart">
                  Height:
                  <input type="range" id="tree-height-slider" min="0.5" max="5" step="0.1" value="1" style="width:60px;cursor:pointer;accent-color:#4a90d9">
                  <span id="tree-height-val" style="width:24px;text-align:right">1</span>&times;
                </label>
                <label class="ctrl-field" title="Horizontal branch-length multiplier — increase to extend branches">
                  Width:
                  <input type="range" id="tree-width-slider" min="0.3" max="4" step="0.1" value="1" style="width:60px;cursor:pointer;accent-color:#4a90d9">
                  <span id="tree-width-val" style="width:24px;text-align:right">1.0</span>&times;
                </label>
                <label class="ctrl-field" title="Vertical space (as a fraction of normal row height) reserved for collapsed OG triangles">
                  Collapsed:
                  <input type="range" id="collapsed-frac-slider" min="0.1" max="1" step="0.05" value="1" style="width:60px;cursor:pointer;accent-color:#4a90d9">
                  <span id="collapsed-frac-val" style="width:28px;text-align:right">1.0</span>
                </label>
                <label class="ctrl-field" title="Opacity of clade highlight rectangles (right-click a node → Highlight to add one)">
                  Hl&#945;:
                  <input type="range" id="clade-hl-alpha-slider" min="0.05" max="0.8" step="0.05" value="0.22" style="width:55px;cursor:pointer;accent-color:#e67e22">
                  <span id="clade-hl-alpha-val" style="width:28px;text-align:right">0.22</span>
                </label>
                <label class="ctrl-field" title="How far (px) the clade highlight rectangle extends past the rightmost leaf">
                  Hl&#8594;:
                  <input type="range" id="clade-hl-extend-slider" min="0" max="300" step="10" value="20" style="width:55px;cursor:pointer;accent-color:#e67e22">
                  <span id="clade-hl-extend-val" style="width:28px;text-align:right">20</span>
                </label>
              </div>
            </details>

            <details class="ctrl-panel">
              <summary>Highlights</summary>
              <div class="ctrl-panel-body vert">
                <div class="ctrl-inline-group">
                  <span class="ctrl-subtle">Species</span>
                  <input id="hl-search" list="hl-list" placeholder="Species… (Enter)" title="Type a species name and press Enter to highlight it">
                  <datalist id="hl-list"></datalist>
                  <button id="hl-clear" onclick="clearHighlight()" title="Remove all species highlights">&#10005;</button>
                  <button class="ctrl-btn" id="btn-mini-sp" onclick="toggleMiniSpPanel(event)" title="Open a mini species tree — click a named node to highlight that entire clade at once">&#x1F333; Species tree</button>
                  <button class="ctrl-btn" id="btn-focus-hl" onclick="focusHighlighted()" style="display:none" title="Collapse all branches not leading to highlighted tips, keeping only the relevant subtree visible">Focus</button>
                </div>
                <div id="hl-tags"></div>
                <div class="ctrl-inline-group">
                  <span class="ctrl-subtle">Orthogroups</span>
                  <div style="position:relative;display:inline-block">
                    <input id="og-hl-search" autocomplete="off" placeholder="OG name… (Enter)" title="Type an OG name and press Enter to highlight it"
                      style="width:180px;font-size:11px;padding:3px 6px;border:1px solid #bbb;border-radius:3px"
                      oninput="ogHlSearchInput(this.value)" onkeydown="ogHlSearchKey(event)">
                    <div id="og-hl-dd" style="display:none;position:absolute;top:100%;left:0;z-index:520;background:#fff;border:1px solid #ccc;border-radius:0 0 4px 4px;max-height:200px;overflow-y:auto;min-width:200px;box-shadow:0 3px 8px rgba(0,0,0,.12);font-size:11px"></div>
                  </div>
                  <button id="og-hl-clear" onclick="clearOgHighlight()" title="Remove all OG highlights">&#10005;</button>
                </div>
                <div id="og-hl-tags"></div>
              </div>
            </details>
          </div>
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
  <div style="font-size:10px;color:#888;margin-bottom:2px;font-weight:600" id="ccp-title">Node options</div>
  <div style="display:flex;gap:6px;flex-wrap:wrap">
    <button id="ccp-triangle" style="padding:4px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer" title="Collapse to a filled triangle (proportional size)">&#x25BD; Triangle</button>
    <button id="ccp-circle" style="padding:4px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer" title="Collapse to a circle with leaf count">&#x25EF; Circle</button>
    <button id="ccp-reroot" style="padding:4px 10px;font-size:11px;border:1px solid #27ae60;border-radius:4px;background:#f0fff4;color:#1a7a42;cursor:pointer" title="Reroot the tree at this node (the selected node becomes one child of the new root)">&#x21C5; Reroot here</button>
    <button id="ccp-focus" style="padding:4px 10px;font-size:11px;border:1px solid #8e44ad;border-radius:4px;background:#faf5ff;color:#6c3483;cursor:pointer" title="Limit the tree view to this subtree">⌕ Focus subtree</button>
    <button id="ccp-compare" style="padding:4px 10px;font-size:11px;border:1px solid #4a90d9;border-radius:4px;background:#f0f6ff;color:#2c6090;cursor:pointer" title="Compare species coverage with another node">&#x2316; Compare</button>
    <button id="ccp-highlight" style="padding:4px 10px;font-size:11px;border:1px solid #e67e22;border-radius:4px;background:#fff8f0;color:#c0622a;cursor:pointer" title="Highlight subtree background with a colour">&#x25A0; Highlight</button>
  </div>
</div>
<div id="tri-action-popup" style="position:fixed;display:none;background:#fff;border:1px solid #bbb;border-radius:6px;padding:6px 8px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.18);z-index:260;flex-direction:column;gap:5px">
  <div style="font-size:10px;color:#888;margin-bottom:2px;font-weight:600" id="tri-popup-title"></div>
  <div style="display:flex;gap:5px;flex-wrap:wrap">
    <button id="tap-expand" style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x25B7; Expand</button>
    <button id="tap-focus" style="padding:3px 10px;font-size:11px;border:1px solid #8e44ad;border-radius:4px;background:#faf5ff;color:#6c3483;cursor:pointer">⌕ Focus</button>
    <button id="tap-rename" style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x270F; Rename</button>
    <button id="tap-color"  style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x1F3A8; Color</button>
    <button id="tap-compare" style="padding:3px 10px;font-size:11px;border:1px solid #4a90d9;border-radius:4px;background:#f0f6ff;color:#2c6090;cursor:pointer">&#x2316; Compare</button>
    <input id="tap-color-input" type="color" style="display:none">
  </div>
</div>
<!-- POSSVM interactive orthogroup calling panel -->
<div id="possvm-panel" style="position:fixed;display:none;background:#fff;border:1px solid #27ae60;border-radius:6px;box-shadow:0 3px 14px rgba(0,0,0,.22);z-index:320;padding:10px 12px;min-width:270px;max-width:350px;font-size:11px">
  <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px">
    <span style="font-weight:700;font-size:12px;color:#1a6b4a">POSSVM Orthogroup Calling</span>
    <button onclick="document.getElementById('possvm-panel').style.display='none'" style="border:none;background:none;cursor:pointer;font-size:15px;color:#888;padding:0 2px;line-height:1">&#x2715;</button>
  </div>
  <!-- SOS threshold -->
  <div style="display:flex;align-items:center;gap:6px;margin-bottom:8px;padding-bottom:8px;border-bottom:1px solid #eee">
    <span style="color:#555;white-space:nowrap">SOS threshold:</span>
    <input type="range" id="possvm-sos" min="0" max="1" step="0.05" value="0" style="flex:1;cursor:pointer;accent-color:#27ae60" oninput="document.getElementById('possvm-sos-val').textContent=parseFloat(this.value).toFixed(2)">
    <span id="possvm-sos-val" style="width:30px;font-weight:600;text-align:right">0.00</span>
    <span style="color:#bbb;cursor:help;font-size:12px" title="Species Overlap Score threshold. 0 = strict: any species overlap → duplication node. Higher values allow overlap before calling a duplication.">&#x24D8;</span>
  </div>
  <!-- Ingroup species -->
  <div style="margin-bottom:6px">
    <div style="color:#555;font-weight:600;margin-bottom:5px">Ingroup species <span id="possvm-sel-count" style="font-weight:normal;color:#888"></span></div>
    <div style="display:flex;gap:4px;margin-bottom:5px;flex-wrap:wrap">
      <button onclick="pvmSelectAll(true)"  style="padding:2px 7px;font-size:10px;border:1px solid #aaa;border-radius:3px;cursor:pointer;background:#f5f5f5">All</button>
      <button onclick="pvmSelectAll(false)" style="padding:2px 7px;font-size:10px;border:1px solid #aaa;border-radius:3px;cursor:pointer;background:#f5f5f5">None</button>
      <button id="possvm-clade-btn" onclick="pvmToggleCladeTree()" style="padding:2px 7px;font-size:10px;border:1px solid #27ae60;border-radius:3px;cursor:pointer;background:#f0fff6;color:#1a6b4a" title="Click a clade in the species tree to use as ingroup">&#x1F333; Pick clade</button>
    </div>
    <!-- mini clade-picker tree (shown on demand) -->
    <div id="possvm-clade-wrap" style="display:none;border:1px solid #d5ead5;border-radius:4px;padding:4px;margin-bottom:5px;overflow:auto;background:#f6fbf6;max-height:170px"></div>
    <!-- species checkboxes -->
    <div id="possvm-sp-list" style="max-height:150px;overflow-y:auto;border:1px solid #e8e8e8;border-radius:3px;padding:3px 5px;background:#fafafa;line-height:1.8"></div>
  </div>
  <!-- Buttons row -->
  <div style="display:flex;gap:6px;align-items:center;flex-wrap:wrap;padding-top:6px;border-top:1px solid #eee">
    <button onclick="pvmMidpointRoot()" style="padding:3px 10px;font-size:11px;border:1px solid #2980b9;border-radius:4px;background:#eaf4fb;color:#1a5276;cursor:pointer" title="Reroot tree at midpoint of longest leaf-to-leaf path (recommended before running POSSVM)">&#x21BB; Midpoint root</button>
    <button onclick="runPossvm()" style="padding:4px 14px;font-size:11px;border:1px solid #27ae60;border-radius:4px;background:#e8f8f0;color:#1a6b4a;cursor:pointer;font-weight:600">&#x25BA; Run</button>
    <button id="possvm-reset-btn" onclick="pvmReset()" style="padding:3px 9px;font-size:11px;border:1px solid #bbb;border-radius:4px;background:#fff;cursor:pointer;display:none">&#x21BA; Reset OGs</button>
  </div>
  <div id="possvm-result" style="margin-top:5px;color:#1a6b4a;font-size:10px;min-height:14px"></div>
  <div style="margin-top:10px;padding-top:8px;border-top:1px solid #eee">
    <details id="hog-details">
      <summary style="color:#555;font-weight:600;cursor:pointer;user-select:none">Hierarchical orthogroups</summary>
      <div style="margin-top:8px">
        <div style="font-size:10px;color:#7b8a93;line-height:1.45;margin-bottom:6px">
          Pick named species-tree clades and run POSSVM on each nested ingroup. The resulting hOG map shows how OGs split or merge across clade levels.
        </div>
        <div style="display:flex;gap:5px;align-items:center;flex-wrap:wrap;margin-bottom:6px">
          <input id="hog-node-search" list="hog-node-list" placeholder="Species-tree node…" style="flex:1;min-width:140px;font-size:11px;padding:3px 6px;border:1px solid #bbb;border-radius:3px">
          <datalist id="hog-node-list"></datalist>
          <button onclick="addHogNodeFromInput()" style="padding:3px 8px;font-size:10px;border:1px solid #aaa;border-radius:3px;background:#f5f5f5;cursor:pointer">Add</button>
          <button id="hog-clade-btn" onclick="pvmToggleHogTree()" style="padding:3px 8px;font-size:10px;border:1px solid #2980b9;border-radius:3px;background:#edf6fd;color:#1a5276;cursor:pointer" title="Pick named clades from the species tree">&#x1F333; Pick from tree</button>
          <button onclick="clearHogNodes()" style="padding:3px 8px;font-size:10px;border:1px solid #aaa;border-radius:3px;background:#f5f5f5;cursor:pointer">Clear</button>
          <button onclick="runHierPossvm()" style="padding:3px 8px;font-size:10px;border:1px solid #2980b9;border-radius:3px;background:#edf6fd;color:#1a5276;cursor:pointer">Run hOGs</button>
        </div>
        <div id="hog-clade-wrap" style="display:none;border:1px solid #d7e7f5;border-radius:4px;padding:4px;margin-bottom:6px;overflow:auto;background:#f7fbff;height:170px;min-height:80px;min-width:200px;resize:both"></div>
        <div id="hog-node-tags" style="display:flex;flex-wrap:wrap;gap:4px;margin-bottom:6px"></div>
        <div id="hog-result" style="color:#1a5276;font-size:10px;min-height:14px"></div>
      </div>
    </details>
  </div>
</div>
<div id="hog-map-panel" style="position:fixed;display:none;top:72px;right:18px;width:min(760px,calc(100vw - 80px));height:min(460px,calc(100vh - 110px));background:#fff;border:1px solid #cfd8de;border-radius:10px;box-shadow:0 10px 30px rgba(0,0,0,.2);z-index:330;overflow:hidden">
  <div id="hog-map-header" style="display:flex;justify-content:space-between;align-items:center;padding:10px 12px;border-bottom:1px solid #e8edf1;background:#f7fafc;cursor:move;user-select:none">
    <div>
      <div style="font-weight:700;font-size:12px;color:#314553">Hierarchical orthogroups</div>
      <div id="hog-map-subtitle" style="font-size:10px;color:#7c8b95;margin-top:2px"></div>
    </div>
    <button onclick="closeHogMapPanel()" style="border:none;background:none;cursor:pointer;font-size:16px;color:#8896a0;line-height:1">&times;</button>
  </div>
  <div id="hog-map-wrap" style="width:100%;height:calc(100% - 52px);overflow:auto;background:#fff"></div>
  <div id="hog-map-resize" title="Drag to resize" style="position:absolute;bottom:0;right:0;width:16px;height:16px;cursor:nwse-resize;opacity:.5"
    onmouseenter="this.style.opacity=1" onmouseleave="this.style.opacity=.5">
    <svg width="16" height="16"><polyline points="6,14 14,6" stroke="#88a" stroke-width="1.5" stroke-linecap="round"/><polyline points="10,14 14,10" stroke="#88a" stroke-width="1.5" stroke-linecap="round"/></svg>
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
const NO_TREE_GENES = %%NO_TREE_GENES_JSON%%;  // {hg_id: {species:[gene_ids]}} for HGs without trees
const DOMAIN_DATA   = %%DOMAIN_DATA_JSON%%;    // {gene_id: [{name,start,end}, ...]}
const GENE_META     = %%GENE_META_JSON%%;      // {gene_id: {length, og_support, ref_ortholog, ref_support}}
const REFNAME_MAP   = %%REFNAME_MAP_JSON%%;    // {gene_id: reference_name}

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
let showOGLabels   = false;       // toggle OG-node labels on expanded internals
let tipFontSize    = null;        // null = auto; number = user override (px)
let treeLinkWidth  = 1.3;        // branch stroke-width in screen px
let treeHeightMult = 1.0;        // vertical stretch multiplier for the tree
let treeWidthMult  = 1.0;        // horizontal stretch multiplier for the tree
let showGeneId    = false;       // tip label parts
let showOGName    = false;
let showRefOrtho  = false;
let showSupport        = false;   // show internal node support values
let showDSNodes        = false;   // show D/S event type at every internal node (POSSVM)
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
let SP_TREE_PRUNED = null;       // pruned tree used for the current display (set by drawSpeciesTree)
const hlTagColors = ["#e74c3c","#3498db","#27ae60","#f39c12","#8e44ad","#16a085","#e67e22","#c0392b"];
function hlTagColor(gi){ const q=hlQueries[gi]; return (q&&hlQueryColors[q])||hlTagColors[gi%hlTagColors.length]; }
function ogHlTagColor(gi){ const q=ogHlQueries[gi]; return (q&&ogHlQueryColors[q])||hlTagColors[gi%hlTagColors.length]; }
const _textMeasureCanvas=document.createElement("canvas");
const _textMeasureCtx=_textMeasureCanvas.getContext("2d");

function measureTextPx(text, fontSize, fontFamily="sans-serif", fontWeight="400"){
  if(!_textMeasureCtx||!text) return 0;
  _textMeasureCtx.font=`${fontWeight} ${fontSize}px ${fontFamily}`;
  return _textMeasureCtx.measureText(text).width;
}

function isLeafLabelVisible(d){
  if(!d||!d.data||!d.data.leaf) return false;
  if(hideNonHl){
    const gid=d.data.gene_id||d.data.name||"";
    if(hmFocusGids!==null){
      if(!hmFocusGids.has(gid)) return false;
    } else {
      const og=leafOgName(d);
      if(hlSet!==null&&!hlSet.has(d.data.species||"")) return false;
      if(ogHlSet!==null&&!ogHlSet.has(og)) return false;
    }
  }
  return true;
}

function leafLabelText(d){
  if(!isLeafLabelVisible(d)) return "";
  const gid=d.data.gene_id||d.data.name||"";
  const og=leafOgName(d);
  const ref=d.data.ref||"";
  const parts=[];
  if(showGeneId&&gid) parts.push(gid);
  if(showOGName&&!showOGLabels&&og) parts.push(og);
  if(showRefOrtho&&ref) parts.push(ref);
  return parts.join(" \u00b7 ");
}

function leafLabelParts(d){
  if(!isLeafLabelVisible(d)) return {gid:"", og:"", ref:""};
  return {
    gid:(d.data.gene_id||d.data.name||""),
    og:leafOgName(d),
    ref:(d.data.ref||""),
  };
}

function treeLabelFontSVG(){ return (ogHlSet!==null ? 2 : 1) * tipFontSVG(); }

function computeLeafLabelLayout(leaves){
  const fs=treeLabelFontSVG();
  const layout={
    gidW:0, ogW:0, refW:0,
    sepW:measureTextPx(" \u00b7 ", fs, "monospace", "400"),
    hasGid:false, hasOg:false, hasRef:false,
  };
  for(const d of leaves){
    const {gid, og, ref}=leafLabelParts(d);
    if(showGeneId&&gid){
      layout.hasGid=true;
      layout.gidW=Math.max(layout.gidW, measureTextPx(gid, fs, "monospace", "400"));
    }
    if(showOGName&&!showOGLabels&&og){
      layout.hasOg=true;
      layout.ogW=Math.max(layout.ogW, measureTextPx(og, fs, "monospace", "400"));
    }
    if(showRefOrtho&&ref){
      layout.hasRef=true;
      layout.refW=Math.max(layout.refW, measureTextPx(ref, fs, "monospace", "400"));
    }
  }
  return layout;
}

function leafLabelRightPx(d, mg, layout){
  const base=nodeX(d,mg)+7;
  const parts=leafLabelParts(d);
  if(!parts.gid&&!parts.og&&!parts.ref) return base;
  if(!layout) layout=computeLeafLabelLayout([d]);
  let width=0;
  if(layout.hasGid) width+=layout.gidW;
  if(layout.hasOg){
    if(layout.hasGid) width+=layout.sepW;
    width+=layout.ogW;
  }
  if(layout.hasRef){
    if(layout.hasGid||layout.hasOg) width+=layout.sepW;
    width+=layout.refW;
  }
  return base + width;
}

function leafLabelColumnXs(layout){
  const gidX=7;
  const ogSepX=gidX + (layout.hasGid ? layout.gidW : 0);
  const ogX=ogSepX + (layout.hasGid ? layout.sepW : 0);
  const refSepX=ogX + (layout.hasOg ? layout.ogW : 0);
  const refX=refSepX + ((layout.hasGid || layout.hasOg) ? layout.sepW : 0);
  return {gidX, ogSepX, ogX, refSepX, refX};
}

function nodeSpeciesSet(node){
  const sps=new Set();
  (function walk(n){
    const ch=n.children||n._children;
    if(!ch||!ch.length){
      const sp=n.data&&n.data.species ? n.data.species : getSpeciesPfx((n.data&&n.data.gene_id)||n.data&&n.data.name||"");
      if(sp) sps.add(sp);
      return;
    }
    for(const c of ch) walk(c);
  })(node);
  return sps;
}

function ogHighlightSubtitle(node){
  if(!node) return "";
  return spMRCAName(nodeSpeciesSet(node)) || "";
}

function leafOgName(d){
  const gid=d&&d.data ? (d.data.gene_id||d.data.name||"") : "";
  const sourceMap=sourceOgGeneMapForCurrentTree();
  return sourceMap[gid]||((d&&d.data&&d.data.og)||"");
}

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
    if (gi !== undefined) return hlTagColor(gi);
    if (window._hogClickGenes && !window._hogClickGenes.has(geneId)) return "#ccc";
    return ogLeaf2Color[geneId] || "#ccc";
  }
  if (window._hogClickGenes && !window._hogClickGenes.has(geneId)) return "#ccc";
  return ogLeaf2Color[geneId] || "#ccc";
}

function ogBaseColor(ogName, fallbackIndex){
  if(ogName && ogName2Color[ogName]) return ogName2Color[ogName];
  if(Number.isFinite(fallbackIndex)) return palette[fallbackIndex % palette.length];
  return "#4a7aad";
}

function naturalSortStrings(arr){
  return [...arr].sort((a,b)=>a.localeCompare(b, undefined, {numeric:true, sensitivity:"base"}));
}

function refNameTokens(values){
  return naturalSortStrings([...new Set(
    values
      .flatMap(v=>String(v||"").split("/"))
      .map(v=>v.trim())
      .filter(Boolean)
  )]);
}

function hogLabelFromGenes(baseOg, genes){
  // Name the OG by the reference gene names present in it, following the same
  // logic as POSSVM's ref_tagcluster: collect ref names from direct REFNAME_MAP
  // entries first, then from per-gene ref_ortholog annotations (which are
  // genuine orthologs assigned by the pipeline run — not "like", just real ones).
  // Only fall back to the raw OG id when no reference information is available.
  const direct=refNameTokens(genes.map(gid=>REFNAME_MAP[gid]).filter(Boolean));
  if(direct.length) return direct.join("/");
  const inferred=refNameTokens(
    genes
      .map(gid=>((GENE_META[gid]||{}).ref_ortholog||"").trim())
      .filter(name=>name && name!=="NA")
  );
  if(inferred.length) return inferred.join("/");
  return baseOg;
}

function hogReferenceSummary(genes){
  const direct=refNameTokens(genes.map(gid=>REFNAME_MAP[gid]).filter(Boolean));
  const inferred=refNameTokens(
    genes
      .map(gid=>((GENE_META[gid]||{}).ref_ortholog||"").trim())
      .filter(name=>name && name!=="NA")
  );
  return {direct, inferred};
}

function relabelOgMapWithReferences(rawOgMap){
  const used=new Set();
  const relabeled={};
  Object.entries(rawOgMap)
    .sort((a,b)=>a[0].localeCompare(b[0], undefined, {numeric:true, sensitivity:"base"}))
    .forEach(([rawOg, genes])=>{
      const baseLabel=hogLabelFromGenes(rawOg, genes);
      let label=baseLabel;
      let suffix=2;
      while(used.has(label)){
        label=`${baseLabel}#${suffix++}`;
      }
      used.add(label);
      relabeled[label]=[...genes];
    });
  return relabeled;
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
let hmBatchSelection = new Set(); // search dropdown items staged for bulk add

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
      const gid = d.data.gene_id||d.data.name||"";
      const sp = d.data.species || "?";
      const meta = GENE_META[gid] || {};
      html += '<div class="tt-row"><span>Species</span><strong style="color:'+spColor(sp)+'">'+sp+'</strong></div>';
      if (REFNAME_MAP[gid]) {
        html += '<div class="tt-row"><span>Reference gene</span><strong style="color:#8e44ad">'+REFNAME_MAP[gid]+'</strong></div>';
      }
      if (currentDetail) {
        for (const [og,genes] of Object.entries(sourceOgsForCurrentTree())) {
          if (genes.includes(gid)) {
            const ogCol=ogBaseColor(og);
            html += '<div class="tt-row"><span>OG</span><strong style="color:'+ogCol+'">'+og+'</strong></div>'; break;
          }
        }
      }
      if (meta.length != null) {
        html += '<div class="tt-row"><span>Protein length</span><strong>'+meta.length+' aa</strong></div>';
      }
      if (meta.og_support) {
        html += '<div class="tt-row"><span>OG support</span><strong>'+meta.og_support+'</strong></div>';
      }
      if (meta.ref_ortholog) {
        html += '<div class="tt-row"><span>Reference ortholog</span><strong>'+meta.ref_ortholog+'</strong></div>';
      }
      if (meta.ref_support) {
        html += '<div class="tt-row"><span>Reference support</span><strong>'+meta.ref_support+'</strong></div>';
      }
      const domains = DOMAIN_DATA[gid] || [];
      if (domains.length) {
        html += '<div class="tt-row"><span>Domain hits</span><strong>'+domains.length+'</strong></div>';
        html += '<div style="margin-top:4px;font-size:10px;color:#666;line-height:1.45">'
          + domains.map(dom =>
              '<div><span style="color:#2c3e50;font-weight:600">'+dom.name
              +'</span> <span style="color:#7f8c8d">'+dom.start+'-'+dom.end+'</span></div>'
            ).join("")
          + '</div>';
      }
    } else {
      const nL = d._children ? countDescLeaves(d._children) : d.leaves().length;
      html += '<div class="tt-row"><span>Subtree leaves</span><strong>'+nL+'</strong></div>';
      html += '<div class="tt-row" style="color:#888"><span>'+(d._children?"collapsed":"expanded")+'</span><span>click to '+(d._children?"expand":"collapse")+'</span></div>';
      if (isOGNode(d) && currentDetail && sourceOgsForCurrentTree()[d.data.name])
        html += '<div class="tt-row"><span>OG members</span><strong>'+sourceOgsForCurrentTree()[d.data.name].length+'</strong></div>';
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
let _origTreeDictForFocus=null;
let _isSubtreeFocused=false;

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
  document.getElementById("ccp-reroot").addEventListener("click",()=>{
    const node=_ccpNode; if(!node) return; hide();
    _userRerooted=true;
    rerootAtNode(node);
  });
  document.getElementById("ccp-focus").addEventListener("click",()=>{
    const node=_ccpNode; if(!node) return; hide();
    focusOnNode(node);
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
    const titleEl=document.getElementById("ccp-title");
    if(titleEl) titleEl.textContent=(d.data.name||"internal node")+" ("+(d.leaves().length)+" leaves)";
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
let cladeHlAlpha  = 0.22;        // global opacity for all clade highlight rectangles
let cladeHlExtend = 20;          // how many px past the rightmost leaf the highlight extends
(function(){
  const pop=document.getElementById("tri-action-popup");
  const title=document.getElementById("tri-popup-title");
  let _nd=null;
  function hide(){ pop.style.display="none"; _nd=null; }
  document.getElementById("tap-expand").addEventListener("click",()=>{
    const node=_nd; if(!node) return; hide();
    node.children=node._children; node._children=null; node._isOgCol=false; renderTree(true);
  });
  document.getElementById("tap-focus").addEventListener("click",()=>{
    const node=_nd; if(!node) return; hide();
    focusOnNode(node);
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
  SP_TREE_PRUNED = tree;
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
    FAMILY_DATA.filter(d=>d.annotated&&(prefix==="all"||d.pref===prefix)).forEach(d=>{
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

  const batchBar=ogSection.length ? `
    <div style="display:flex;align-items:center;gap:6px;padding:5px 10px;background:#fafcff;border-bottom:1px solid #e7eef7;font-size:10px;position:sticky;top:0;z-index:1">
      <button type="button" onclick="hmSearchSelectVisible()" style="padding:1px 6px;font-size:10px;border:1px solid #bfd0e3;border-radius:3px;background:#fff;cursor:pointer">Select visible</button>
      <button type="button" onclick="hmSearchClearSelected()" style="padding:1px 6px;font-size:10px;border:1px solid #d6d6d6;border-radius:3px;background:#fff;cursor:pointer">Clear</button>
      <button type="button" onclick="hmSearchAddSelected()" style="margin-left:auto;padding:1px 8px;font-size:10px;border:1px solid #4a90d9;color:#4a90d9;border-radius:3px;background:#fff;cursor:pointer">Add selected</button>
    </div>` : "";
  dd.innerHTML=batchBar+_hmSearchHits.map((h,i)=>{
    if(h.kind==="sep") return `<div style="padding:3px 10px;font-size:10px;color:#aaa;background:#f5f5f5;border-top:1px solid #eee;border-bottom:1px solid #eee;font-weight:600;letter-spacing:.05em;text-transform:uppercase">${h.label}</div>`;
    if(h.kind==="nav") return `<div data-i="${i}" style="padding:4px 10px;cursor:pointer;border-bottom:1px solid #f0f0f0;color:${h.type==="hg"?"#333":h.type==="family"?"#555":"#888"}" onmousedown="hmTextSearchGo(${i})">&#8594; ${h.label}</div>`;
    if(h.type==="group"){
      const done=hmSearchItemDone(h);
      const icon=h.groupType==="class"?"\uD83D\uDCC2":"\uD83E\uDDF9";
      const checked=hmBatchSelection.has(hmSearchItemKey(h));
      return `<div data-i="${i}" style="padding:4px 10px;cursor:pointer;border-bottom:1px solid #f0f0f0;background:${done?"#e8f5e9":checked?"#eaf3ff":"#f0f6ff"};color:${done?"#2e7d32":"#1a4a7a"};display:flex;align-items:center;gap:7px" onmousedown="hmTextSearchGo(${i})">
        <input type="checkbox" ${checked?'checked':''} ${done?'disabled':''} onmousedown="event.stopPropagation()" onclick="event.stopPropagation()" onchange="hmSearchToggle(${i})">
        <span>${icon} <b>${h.key}</b><span style="color:#aaa;font-size:10px"> \u00b7 ${h.groupType} \u00b7 ${h.count} OGs</span></span>
        ${done?`<span style="margin-left:auto">&#10003;</span>`:""}
      </div>`;
    }
    const r=hmOGIndex[h.og]||{family:"",total:0}; const already=hmSearchItemDone(h);
    const checked=hmBatchSelection.has(hmSearchItemKey(h));
    const lbl=h.matchedGene?`<span style="color:#888;font-size:10px">gene</span> <b>${h.matchedGene}</b><span style="color:#aaa;font-size:10px"> \u2192 ${h.og} \u00b7 ${r.family}</span>`:`<b>${h.og}</b><span style="color:#aaa;font-size:10px"> \u00b7 ${r.family} \u00b7 ${r.total} genes</span>`;
    return `<div data-i="${i}" style="padding:4px 10px;cursor:pointer;border-bottom:1px solid #f0f0f0;background:${already?"#e8f5e9":checked?"#eaf3ff":"#fff"};color:${already?"#2e7d32":"#222"};display:flex;align-items:center;gap:7px" onmousedown="hmTextSearchGo(${i})">
      <input type="checkbox" ${checked?'checked':''} ${already?'disabled':''} onmousedown="event.stopPropagation()" onclick="event.stopPropagation()" onchange="hmSearchToggle(${i})">
      <span>${lbl}</span>
      ${already?`<span style="margin-left:auto;color:#27ae60">&#10003;</span>`:""}
    </div>`;
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
  if(h.kind==="nav"){
    document.getElementById("hm-search-dd").style.display="none";
    document.getElementById("hm-text-search").value=""; hmTextFilter="";
    if(h.type==="class"){ hmViewMode="class"; hmActiveClass=h.cls; hmActiveFamily=hmActiveHG=hmActiveHGRec=null; }
    else if(h.type==="family"){ hmViewMode="family"; hmActiveClass=h.cls; hmActiveFamily=h.fam; hmActiveHG=hmActiveHGRec=null; }
    else if(h.type==="hg"){
      const rec=TREE_INDEX.find(r=>r.id===h.hg_id); if(!rec) return;
      hmViewMode="og"; hmActiveClass=h.cls; hmActiveFamily=h.fam; hmActiveHG=h.hg_id; hmActiveHGRec=rec;
    }
  } else {
    hmSearchToggle(i);
    return;
  }
  drawHeatmap();
}

function hmSearchItemKey(h){
  if(!h) return "";
  if(h.type==="group") return `group:${h.groupType}:${h.key}`;
  if(h.type==="og") return `og:${h.og}`;
  return "";
}

function hmSearchItemDone(h){
  if(!h || h.kind==="nav" || h.kind==="sep") return false;
  const effOGs=new Set(getEffectiveCustomOGs());
  if(h.type==="group")
    return hmCustomGroups.some(g=>g.groupType===h.groupType&&g.key===h.key) || h.ogs.every(og=>effOGs.has(og));
  return effOGs.has(h.og);
}

function hmSearchToggle(i){
  const h=_hmSearchHits[i];
  if(!h || h.kind==="nav" || h.kind==="sep" || hmSearchItemDone(h)) return;
  const key=hmSearchItemKey(h);
  if(!key) return;
  if(hmBatchSelection.has(key)) hmBatchSelection.delete(key);
  else hmBatchSelection.add(key);
  hmTextSearchInput(document.getElementById("hm-text-search").value||"");
}

function hmAddSearchHit(h){
  if(!h || hmSearchItemDone(h)) return false;
  if(h.type==="group"){
    if(!hmCustomGroups.some(g=>g.groupType===h.groupType&&g.key===h.key))
      hmCustomGroups.push({groupType:h.groupType,key:h.key,label:h.key,ogs:h.ogs});
  } else if(h.type==="og"){
    if(!hmCustomOGs.includes(h.og)) hmCustomOGs.push(h.og);
  } else {
    return false;
  }
  return true;
}

function hmSearchAddSelected(){
  let changed=false;
  _hmSearchHits.forEach(h=>{
    const key=hmSearchItemKey(h);
    if(key && hmBatchSelection.has(key)) changed = hmAddSearchHit(h) || changed;
  });
  hmBatchSelection.clear();
  if(!changed){
    hmTextSearchInput(document.getElementById("hm-text-search").value||"");
    return;
  }
  hmViewMode="custom";
  const bar=document.getElementById("hm-custom-bar");
  bar.style.display="flex";
  _positionCustomBar();
  renderCustomChips();
  drawHeatmap();
  hmTextSearchInput(document.getElementById("hm-text-search").value||"");
}

function hmSearchSelectVisible(){
  _hmSearchHits.forEach(h=>{
    const key=hmSearchItemKey(h);
    if(key && !hmSearchItemDone(h)) hmBatchSelection.add(key);
  });
  hmTextSearchInput(document.getElementById("hm-text-search").value||"");
}

function hmSearchClearSelected(){
  hmBatchSelection.clear();
  hmTextSearchInput(document.getElementById("hm-text-search").value||"");
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
  // Add genes from HGs that have no gene tree at all
  const NO_ANNO="(no annotation)";
  for(const [hgId, spGenes] of Object.entries(NO_TREE_GENES)){
    const hgRec=HG_DATA.find(d=>d.id===hgId)||{};
    const sc={};
    for(const [sp,genes] of Object.entries(spGenes)){
      sc[sp]=genes.length;
      for(const g of genes){
        if(!hmGeneIndex[g]) hmGeneIndex[g]=NO_ANNO+"::"+hgId;
      }
    }
    const key=NO_ANNO+"::"+hgId;
    hmOGIndex[key]={hgId, family:hgRec.family||"", cls:hgRec.class||hgRec.pref||"",
                    total:Object.values(sc).reduce((a,b)=>a+b,0), species_counts:sc};
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
document.getElementById("hog-node-search").addEventListener("keydown",function(e){
  if(e.key==="Enter"&&this.value.trim()){ addHogNodeFromInput(); e.preventDefault(); }
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
document.getElementById("clade-hl-extend-slider").addEventListener("input",function(){
  cladeHlExtend=+this.value;
  document.getElementById("clade-hl-extend-val").textContent=this.value;
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
syncOgTextControl();

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
  // Use pruned tree when prune-to-data is active, otherwise full tree
  const root = (spPruneToData && SP_TREE_PRUNED) ? SP_TREE_PRUNED : SP_TREE_DATA;
  const nwk=treeToNewick(root)+";";
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
    // Strip internal "::hg_id" suffix used as a unique key for unannotated genes
    const ogDisplay=og.includes("::")? og.split("::")[0] : og;
    rows.push([gene,ogDisplay,r.hgId||"",r.family||"",r.cls||""].join("\t"));
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
  const classCounts={};   // cls → gene count for this species
  let total=0;
  for(const [gene,og] of Object.entries(hmGeneIndex)){
    if(getSpeciesPfx(gene)!==sp) continue;
    const cls=(hmOGIndex[og]||{}).cls||"other";
    classCounts[cls]=(classCounts[cls]||0)+1;
    total++;
  }
  const sorted=Object.keys(classCounts).sort();
  if(!sorted.length){ downloadSpeciesGenes(sp); return; }
  const pop=document.getElementById("sp-annot-popup");
  document.getElementById("sp-annot-popup-title").textContent=sp;
  const btns=document.getElementById("sp-annot-popup-btns"); btns.innerHTML="";
  const btnStyle="padding:3px 8px;font-size:11px;border:1px solid #aaa;border-radius:3px;background:#f8f8f8;cursor:pointer;text-align:left;width:100%;display:flex;justify-content:space-between;gap:8px";
  const mkBtn=(label,cls,count)=>{
    const b=document.createElement("button"); b.style.cssText=btnStyle;
    const nameSpan=document.createElement("span"); nameSpan.textContent=label;
    const countSpan=document.createElement("span");
    countSpan.textContent=count; countSpan.style.cssText="color:#888;font-variant-numeric:tabular-nums;flex-shrink:0";
    b.appendChild(nameSpan); b.appendChild(countSpan);
    b.onclick=()=>{ pop.style.display="none"; downloadSpeciesGenes(sp,cls); };
    btns.appendChild(b);
  };
  mkBtn("All classes", undefined, total);
  sorted.forEach(c=>mkBtn(c, c, classCounts[c]));
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
  _pvmActive=false; _pvmIngroupSps=null; _pvmOgs=null;
  window._hogClickGenes=null;
  cladeSp2Color={}; cladeSp2Group={}; cladeGrpColor={}; ogLeaf2Color={}; ogName2Color={}; ogGene2Name={};
  cladeHighlights.clear();
  _pvmHogModel=null; closeHogMapPanel();
  document.getElementById("possvm-reset-btn").style.display="none";
  document.getElementById("possvm-result").textContent="";
  document.getElementById("hog-result").textContent="";
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
  window._hogClickGenes=null;
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
  _pvmActive=false; _pvmIngroupSps=null; _pvmOgs=null;
  _pvmHogModel=null; closeHogMapPanel();
  document.getElementById("possvm-reset-btn").style.display="none";
  document.getElementById("possvm-result").textContent="";
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
  if(_pvmActive&&_pvmOgs) return _pvmOgs;
  if(treeSource==="original"&&currentDetail)
    return currentDetail.prev_ogs||currentDetail.ogs||{};
  return currentDetail?currentDetail.ogs:{};
}

function sourceOgsForCurrentTree(){
  if(treeSource==="original"&&currentDetail)
    return currentDetail.prev_ogs||currentDetail.ogs||{};
  return currentDetail?currentDetail.ogs:{};
}

function sourceOgGeneMapForCurrentTree(){
  const gene2og={};
  for(const [og, genes] of Object.entries(sourceOgsForCurrentTree()||{})){
    (genes||[]).forEach(gid=>{ if(!(gid in gene2og)) gene2og[gid]=og; });
  }
  return gene2og;
}

const treeSvg = d3.select("#tree-svg");
let rootNode=null, gMain=null, _uid=0, _zoom=null, _zoomScale=1;
let _compareNode1=null;  // first node selected for species comparison
let useBranchLen=false, _phyloScale=1;

function tipFontSVG(){ return (tipFontSize!==null?tipFontSize:11)/_zoomScale; }
function applyTipFontSize(){
  if(!gMain) return;
  const fs=tipFontSVG();
  const labelFs=treeLabelFontSVG();
  gMain.selectAll(".leaf-label").attr("font-size",d=>d&&d.data&&d.data.leaf?labelFs:0);
  gMain.selectAll(".og-label").attr("font-size",labelFs);
  // update MRCA sub-label tspan inside og-label
  gMain.selectAll(".og-label tspan:nth-child(2)").attr("font-size",Math.max(7,labelFs*0.82));
  const colR2=Math.max(10,fs*0.9), leafR2=fs*0.36, colCx2=-(colR2-leafR2);
  gMain.selectAll(".count-label").attr("font-size",fs).attr("x",d=>d&&d._children&&!d._isOgCol?colCx2:null);
  gMain.selectAll("circle").filter(d=>d&&d._children&&!d._isOgCol).attr("r",colR2).attr("cx",colCx2);
  gMain.selectAll("circle").filter(d=>d&&!d._children).attr("r",d=>d.data&&d.data.leaf?fs*0.36:(_pvmActive&&d._pvmEvent==="D")||isOGNode(d)||d.data._og_label?fs*0.55:fs*0.26);
  gMain.selectAll(".link").attr("stroke-width",treeLinkWidth/_zoomScale);
}

function isOGNode(d){ return !d.data.leaf && d.data.name && sourceOgsForCurrentTree()[d.data.name]!==undefined; }

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

function leafGeneId(d){
  return d&&d.data ? (d.data.gene_id||d.data.name||"") : "";
}

function collectLeafGeneIds(node, out){
  if(!node) return out;
  const ch=node.children||node._children;
  if(!ch||!ch.length){
    const gid=leafGeneId(node);
    if(gid) out.push(gid);
    return out;
  }
  for(const c of ch) collectLeafGeneIds(c, out);
  return out;
}

function findExactOgRoot(leaves){
  if(!leaves||!leaves.length) return null;
  const node=findMRCA(leaves);
  if(!node||node.data.leaf) return null;
  const want=new Set(leaves.map(leafGeneId).filter(Boolean));
  if(!want.size) return null;
  const have=collectLeafGeneIds(node, []);
  if(have.length!==want.size) return null;
  for(const gid of have){
    if(!want.has(gid)) return null;
  }
  return node;
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

/** Orthogonal elbow path for a cladogram/phylogram link.
 *  Phylogenetic trees are rendered with hard right-angle joints:
 *  vertical stem first, then horizontal to child. */
function elbowPath(s, t, mg, r){
  const sx=nodeX(s,mg), sy=s.x+mg.top;
  const tx=nodeX(t,mg), ty=t.x+mg.top;
  return `M${sx},${sy}V${ty}H${tx}`;
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
  const sourceMap=sourceOgGeneMapForCurrentTree();
  (function walkAll(n){
    if(n.data.leaf){
      const gid=n.data.gene_id||n.data.name||"";
      const og=sourceMap[gid]||(n.data.og||"");
      if(og)(ogGroups[og]=ogGroups[og]||[]).push(n);
      return;
    }
    for(const c of (n.children||[])) walkAll(c);
    for(const c of (n._children||[])) walkAll(c);
  })(rootNode);
  for(const [og,leaves] of Object.entries(ogGroups)){
    const node=findExactOgRoot(leaves);
    if(visibleNodes.has(node)&&!node.data.leaf&&!node._children){ node.data._og_label=og; }
  }
}

function toggleOGLabels(){
  showOGLabels=!showOGLabels;
  document.getElementById("btn-og-labels").classList.toggle("active-btn", showOGLabels);
  syncOgTextControl();
  if(rootNode) renderTree(false);
}

function syncOgTextControl(){
  const chk=document.getElementById("chk-og");
  if(!chk) return;
  chk.disabled=showOGLabels;
  chk.title=showOGLabels
    ? "Tip-level OG text is hidden while internal OG labels are enabled."
    : "";
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

function toggleDSNodes(){
  showDSNodes=!showDSNodes;
  document.getElementById("btn-ds-nodes").classList.toggle("active-btn", showDSNodes);
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
  // collect only species present in current tree
  const activeSp=currentIndex?new Set(currentIndex.species):new Set(ALL_SPECIES);
  function prune(n){
    if(!n.children) return activeSp.has(n.name)?n:null;
    const k=n.children.map(prune).filter(Boolean);
    if(!k.length) return null;
    if(k.length===1) return k[0];
    n.children=k; return n;
  }
  const tree=prune(clone(SP_TREE_DATA));
  if(!tree){
    wrap.innerHTML='<span style="color:#888;font-size:10px;padding:4px">No matching species.</span>';
    return;
  }
  const allN=flat(tree);
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
  leaves.forEach(l=>{ svg.append("circle").attr("cx",sx(maxD)).attr("cy",l._y).attr("r",3).attr("fill","#2c3e50"); svg.append("text").attr("x",sx(maxD)+6).attr("y",l._y).attr("dy","0.35em").attr("font-size",9).attr("fill","#333").attr("font-family","monospace").text(l.name); });
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

// ── POSSVM interactive orthogroup calling ─────────────────────────────────────
let _pvmActive=false;          // true while POSSVM OGs are displayed
let _pvmIngroupSps=null;       // Set of ingroup species used in last POSSVM run
let _pvmOgs=null;              // interactive POSSVM OGs for the current tree
let _pvmMidpointApplied=false; // midpoint root was applied for current tree
let _userRerooted=false;       // user explicitly rerooted via the clade popup
const _pvmOgHlPalette=["#a8d8ea","#ffcef3","#d4f1c4","#ffd6a5","#e2d0f8","#c8f0de","#ffe4a8","#ffd0d0","#b3e5d4","#fce4b4"];
let _pvmCladeOpen=false;       // whether the clade-picker tree is visible
let _pvmHogCladeOpen=false;    // whether the hOG clade-picker tree is visible
let _pvmHogNodes=[];           // selected named species-tree clades for hierarchical OGs
let _pvmHogModel=null;         // rendered hOG metro-map model

function togglePossvmPanel(ev){
  const panel=document.getElementById("possvm-panel");
  if(panel.style.display==="block"){ panel.style.display="none"; return; }
  const btn=document.getElementById("btn-possvm");
  const r=btn.getBoundingClientRect();
  panel.style.top=(r.bottom+6)+"px";
  panel.style.left=Math.max(4,r.left-60)+"px";
  panel.style.display="block";
  pvmBuildSpList();
  pvmBuildHogNodeList();
  renderHogNodeTags();
}
document.addEventListener("click",ev=>{
  const panel=document.getElementById("possvm-panel");
  if(panel.style.display==="block"&&!panel.contains(ev.target)&&ev.target.id!=="btn-possvm")
    panel.style.display="none";
});

// Populate the species checkbox list from the current gene tree
function pvmBuildSpList(){
  const list=document.getElementById("possvm-sp-list");
  list.innerHTML="";
  const treeSpecies=pvmGetTreeSpecies();
  treeSpecies.forEach(sp=>{
    const lbl=document.createElement("label");
    lbl.style.cssText="display:flex;align-items:center;gap:5px;padding:1px 2px;cursor:pointer";
    const cb=document.createElement("input");
    cb.type="checkbox"; cb.value=sp; cb.checked=true;
    cb.setAttribute("data-pvm-sp",sp);
    cb.addEventListener("change",pvmUpdateSelCount);
    const dot=document.createElement("span");
    dot.style.cssText=`display:inline-block;width:9px;height:9px;border-radius:50%;flex-shrink:0;background:${spColor(sp)}`;
    lbl.append(cb,dot,document.createTextNode("\u00a0"+sp));
    list.appendChild(lbl);
  });
  pvmUpdateSelCount();
}

function getCurrentSpeciesSet(){
  return new Set(pvmGetTreeSpecies());
}

function getNamedSpeciesTreeClades(){
  const currentSpecies=getCurrentSpeciesSet();
  const out=[];
  (function walk(node, depth){
    if(!node) return [];
    if(!node.children||!node.children.length){
      const keep=currentSpecies.has(node.name);
      return keep ? [node.name] : [];
    }
    let leaves=[];
    for(const child of node.children) leaves=leaves.concat(walk(child, depth+1));
    if(node.name && leaves.length){
      out.push({name:node.name, depth, species:[...new Set(leaves)]});
    }
    return leaves;
  })(SP_TREE_DATA, 0);
  out.sort((a,b)=>a.depth-b.depth || a.name.localeCompare(b.name));
  return out;
}

function pvmBuildHogNodeList(){
  const dl=document.getElementById("hog-node-list");
  if(!dl) return;
  dl.innerHTML=getNamedSpeciesTreeClades()
    .map(rec=>`<option value="${rec.name.replace(/"/g,"&quot;")}"></option>`)
    .join("");
}

function renderHogNodeTags(){
  const wrap=document.getElementById("hog-node-tags");
  if(!wrap) return;
  wrap.innerHTML=_pvmHogNodes.map(name=>
    `<span class="hl-tag" style="background:#4a90d9">${name}<span class="hl-tag-x" onclick="removeHogNode('${name.replace(/'/g,"\\'")}')">&times;</span></span>`
  ).join("");
}

function addHogNode(name){
  const q=(name||"").trim();
  if(!q) return;
  const avail=getNamedSpeciesTreeClades();
  const match=avail.find(rec=>rec.name===q);
  if(!match){
    document.getElementById("hog-result").textContent=`Unknown clade: ${q}`;
    return;
  }
  if(!_pvmHogNodes.includes(match.name)) _pvmHogNodes.push(match.name);
  renderHogNodeTags();
  document.getElementById("hog-result").textContent=`${_pvmHogNodes.length} clade${_pvmHogNodes.length!==1?"s":""} selected`;
  if(_pvmHogCladeOpen) pvmDrawHogCladeTree();
}

function addHogNodeFromInput(){
  const inp=document.getElementById("hog-node-search");
  addHogNode(inp.value);
  inp.value="";
}

function removeHogNode(name){
  _pvmHogNodes=_pvmHogNodes.filter(n=>n!==name);
  renderHogNodeTags();
  document.getElementById("hog-result").textContent=_pvmHogNodes.length
    ? `${_pvmHogNodes.length} clade${_pvmHogNodes.length!==1?"s":""} selected`
    : "";
  if(_pvmHogCladeOpen) pvmDrawHogCladeTree();
}

function clearHogNodes(){
  _pvmHogNodes=[];
  renderHogNodeTags();
  const res=document.getElementById("hog-result");
  if(res) res.textContent="";
  if(_pvmHogCladeOpen) pvmDrawHogCladeTree();
}

// Return species in the current gene tree, in speciesOrder order
function pvmGetTreeSpecies(){
  if(!rootNode) return [];
  const sps=new Set();
  rootNode.each(d=>{if(d.data.leaf&&d.data.species) sps.add(d.data.species);});
  const ordered=speciesOrder.filter(s=>sps.has(s));
  sps.forEach(s=>{if(!ordered.includes(s)) ordered.push(s);});
  return ordered;
}

function pvmSelectAll(val){
  document.querySelectorAll("[data-pvm-sp]").forEach(cb=>cb.checked=val);
  pvmUpdateSelCount();
}

function pvmUpdateSelCount(){
  const n=document.querySelectorAll("[data-pvm-sp]:checked").length;
  const tot=document.querySelectorAll("[data-pvm-sp]").length;
  document.getElementById("possvm-sel-count").textContent="("+n+"/"+tot+")";
}

function pvmToggleCladeTree(){
  _pvmCladeOpen=!_pvmCladeOpen;
  document.getElementById("possvm-clade-wrap").style.display=_pvmCladeOpen?"block":"none";
  document.getElementById("possvm-clade-btn").classList.toggle("active-btn",_pvmCladeOpen);
  if(_pvmCladeOpen) pvmDrawCladeTree();
}

function pvmToggleHogTree(){
  _pvmHogCladeOpen=!_pvmHogCladeOpen;
  document.getElementById("hog-clade-wrap").style.display=_pvmHogCladeOpen?"block":"none";
  document.getElementById("hog-clade-btn").classList.toggle("active-btn",_pvmHogCladeOpen);
  if(_pvmHogCladeOpen) pvmDrawHogCladeTree();
}

// Draw a mini species tree in the clade-picker area; clicking a named node
// selects that clade's species as the ingroup
function pvmDrawCladeTree(){
  const wrap=document.getElementById("possvm-clade-wrap");
  wrap.innerHTML="";
  if(!SP_TREE_DATA){ wrap.innerHTML='<span style="color:#888;font-size:10px;padding:4px">No species tree loaded.</span>'; return; }
  function flat(n){return[n].concat(n.children?n.children.flatMap(flat):[]);}
  function clone(n){return Object.assign({},n,{children:n.children?n.children.map(clone):null});}
  const treeSpSet=new Set(pvmGetTreeSpecies());
  function prune(n){
    if(!n.children) return treeSpSet.has(n.name)?n:null;
    const k=n.children.map(prune).filter(Boolean);
    if(!k.length) return null;
    if(k.length===1) return k[0];
    n.children=k; return n;
  }
  const tree=prune(clone(SP_TREE_DATA));
  if(!tree){ wrap.innerHTML='<span style="color:#888;font-size:10px;padding:4px">No matching species.</span>'; return; }
  const allN=flat(tree);
  const leaves=allN.filter(n=>!n.children);
  const leafH=14, lM=6, W=270;
  leaves.forEach((l,i)=>{l._y=i*leafH+leafH/2;});
  function assignY(n){if(n.children){n.children.forEach(assignY);n._y=(n.children[0]._y+n.children[n.children.length-1]._y)/2;}}
  assignY(tree);
  let maxD=0;
  function assignD(n,d){n._d=d;maxD=Math.max(maxD,d);if(n.children)n.children.forEach(c=>assignD(c,d+1));}
  assignD(tree,0);
  const sx=d=>lM+(d/Math.max(1,maxD))*(W-lM-80);
  const H=leaves.length*leafH+10;
  const svg=d3.select(wrap).append("svg").attr("width",W).attr("height",H);
  function drawB(n){
    if(!n.children) return;
    const ys=n.children.map(c=>c._y);
    svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(n._d)).attr("y1",d3.min(ys)).attr("y2",d3.max(ys)).attr("stroke","#ccc");
    n.children.forEach(c=>{svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(c._d)).attr("y1",c._y).attr("y2",c._y).attr("stroke","#ccc");drawB(c);});
  }
  drawB(tree);
  leaves.forEach(l=>{
    svg.append("circle").attr("cx",sx(maxD)).attr("cy",l._y).attr("r",3).attr("fill",spColor(l.name));
    svg.append("text").attr("x",sx(maxD)+6).attr("y",l._y).attr("dy","0.35em").attr("font-size",9).attr("fill","#333").attr("font-family","monospace").text(l.name);
  });
  function drawInternals(n){
    if(!n.children) return;
    if(n.name){
      function sp2(nd){return nd.children?nd.children.flatMap(sp2):[nd.name];}
      const cladeSps=sp2(n).filter(s=>treeSpSet.has(s));
      // Highlight if currently selected
      const checked=new Set([...document.querySelectorAll("[data-pvm-sp]:checked")].map(c=>c.value));
      const allSelected=cladeSps.length>0&&cladeSps.every(s=>checked.has(s));
      svg.append("text")
        .attr("x",sx(n._d)).attr("y",n._y-5)
        .attr("text-anchor","middle")
        .attr("font-size",9).attr("fill",allSelected?"#e67e22":"#27ae60").attr("font-style","italic")
        .style("cursor","pointer").style("font-weight",allSelected?"700":"400").text(n.name)
        .on("click",(ev)=>{
          ev.stopPropagation();
          pvmSelectAll(false);
          document.querySelectorAll("[data-pvm-sp]").forEach(cb=>{if(cladeSps.includes(cb.value)) cb.checked=true;});
          pvmUpdateSelCount();
          pvmDrawCladeTree(); // refresh colours
        });
    }
    n.children.forEach(drawInternals);
  }
  drawInternals(tree);
}

function pvmDrawHogCladeTree(){
  const wrap=document.getElementById("hog-clade-wrap");
  wrap.innerHTML="";
  if(!SP_TREE_DATA){ wrap.innerHTML='<span style="color:#888;font-size:10px;padding:4px">No species tree loaded.</span>'; return; }
  // Redraw when the user resizes the pane
  if(!wrap._hogResizeObs){
    wrap._hogResizeObs=new ResizeObserver(()=>pvmDrawHogCladeTree());
    wrap._hogResizeObs.observe(wrap);
  }
  function flat(n){return[n].concat(n.children?n.children.flatMap(flat):[]);}
  function clone(n){return Object.assign({},n,{children:n.children?n.children.map(clone):null});}
  const treeSpSet=new Set(pvmGetTreeSpecies());
  function prune(n){
    if(!n.children) return treeSpSet.has(n.name)?n:null;
    const k=n.children.map(prune).filter(Boolean);
    if(!k.length) return null;
    if(k.length===1) return k[0];
    n.children=k; return n;
  }
  const tree=prune(clone(SP_TREE_DATA));
  if(!tree){ wrap.innerHTML='<span style="color:#888;font-size:10px;padding:4px">No matching species.</span>'; return; }
  const allN=flat(tree);
  const leaves=allN.filter(n=>!n.children);
  const leafH=14, lM=6, W=Math.max(220, wrap.clientWidth-10);
  leaves.forEach((l,i)=>{l._y=i*leafH+leafH/2;});
  function assignY(n){if(n.children){n.children.forEach(assignY);n._y=(n.children[0]._y+n.children[n.children.length-1]._y)/2;}}
  assignY(tree);
  let maxD=0;
  function assignD(n,d){n._d=d;maxD=Math.max(maxD,d);if(n.children)n.children.forEach(c=>assignD(c,d+1));}
  assignD(tree,0);
  const sx=d=>lM+(d/Math.max(1,maxD))*(W-lM-80);
  const H=leaves.length*leafH+10;
  const svg=d3.select(wrap).append("svg").attr("width",W).attr("height",H);
  function drawB(n){
    if(!n.children) return;
    const ys=n.children.map(c=>c._y);
    svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(n._d)).attr("y1",d3.min(ys)).attr("y2",d3.max(ys)).attr("stroke","#ccc");
    n.children.forEach(c=>{svg.append("line").attr("x1",sx(n._d)).attr("x2",sx(c._d)).attr("y1",c._y).attr("y2",c._y).attr("stroke","#ccc");drawB(c);});
  }
  drawB(tree);
  leaves.forEach(l=>{
    svg.append("circle").attr("cx",sx(maxD)).attr("cy",l._y).attr("r",3).attr("fill",spColor(l.name));
    svg.append("text").attr("x",sx(maxD)+6).attr("y",l._y).attr("dy","0.35em").attr("font-size",9).attr("fill","#333").attr("font-family","monospace").text(l.name);
  });
  function drawInternals(n){
    if(!n.children) return;
    if(n.name){
      function sp2(nd){return nd.children?nd.children.flatMap(sp2):[nd.name];}
      const cladeSps=sp2(n).filter(s=>treeSpSet.has(s));
      const isSelected=_pvmHogNodes.includes(n.name);
      svg.append("text")
        .attr("x",sx(n._d)).attr("y",n._y-5)
        .attr("text-anchor","middle")
        .attr("font-size",9).attr("fill",isSelected?"#e67e22":"#2980b9").attr("font-style","italic")
        .style("cursor","pointer").style("font-weight",isSelected?"700":"400").text(n.name)
        .on("click",(ev)=>{
          ev.stopPropagation();
          if(cladeSps.length) addHogNode(n.name);
          pvmDrawHogCladeTree();
        });
    }
    n.children.forEach(drawInternals);
  }
  drawInternals(tree);
}

// ── Midpoint rooting ──────────────────────────────────────────────────────────
// Roots the current gene tree at the node closest to the midpoint of the
// longest leaf-to-leaf path.  Uses branch lengths when available.
function pvmMidpointRoot(){
  if(!rootNode) return;
  // Expand all so we see the full topology
  rootNode.each(n=>{if(n._children){n.children=n._children;n._children=null;}});
  // Assign cumulative root-distances
  rootNode._rDist=0;
  rootNode.each(n=>{
    if(n.children) n.children.forEach(c=>{
      c._rDist=(n._rDist||0)+(useBranchLen?(c.data.dist||0):1);
    });
  });
  const leaves=rootNode.leaves();
  if(leaves.length<2) return;
  // Two-sweep to find diameter endpoints
  function farthestFrom(src){
    let best=src, bestD=-1;
    // Build ancestor set once for src
    const ancSrc=new Set();
    let tmp=src; while(tmp){ancSrc.add(tmp);tmp=tmp.parent;}
    leaves.forEach(l=>{
      let cur=l;
      while(cur){if(ancSrc.has(cur)){const d=(src._rDist+l._rDist-2*(cur._rDist||0));if(d>bestD){bestD=d;best=l;}break;}cur=cur.parent;}
    });
    return best;
  }
  const L1=farthestFrom(leaves[0]);
  const L2=farthestFrom(L1);
  // Compute actual diameter
  const ancL1=[];
  let cur=L1; while(cur){ancL1.push(cur);cur=cur.parent;}
  const ancL1Set=new Set(ancL1);
  let lcaNode=L2; while(lcaNode&&!ancL1Set.has(lcaNode)) lcaNode=lcaNode.parent;
  if(!lcaNode) return;
  const diam=L1._rDist+L2._rDist-2*(lcaNode._rDist||0);
  const half=diam/2;
  // Build the full path L1→LCA→L2 and pick node nearest to half from L1
  const pathUp=[];
  cur=L1; while(cur&&cur!==lcaNode){pathUp.push(cur);cur=cur.parent;} pathUp.push(lcaNode);
  const pathDown=[];
  cur=L2; while(cur&&cur!==lcaNode){pathDown.unshift(cur);cur=cur.parent;}
  const fullPath=[...pathUp,...pathDown];
  // For each node on the path, compute d(L1, node)
  function distL1(node){
    // find lca(L1, node)
    const anc=new Set(); let c=L1; while(c){anc.add(c);c=c.parent;}
    c=node; while(c){if(anc.has(c)) return L1._rDist+node._rDist-2*(c._rDist||0); c=c.parent;}
    return 0;
  }
  const dists=fullPath.map(n=>distL1(n));
  let bestNode=null, bestSplit=0.5, bestDiff=Infinity;
  for(let i=1;i<fullPath.length;i++){
    const a=fullPath[i-1], b=fullPath[i];
    const da=dists[i-1], db=dists[i];
    const lo=Math.min(da,db), hi=Math.max(da,db);
    if(half<lo||half>hi) continue;
    if(b.parent===a){
      const edge=Math.max(db-da, 0);
      const split=edge>0 ? (half-da)/edge : 0.5;
      bestNode=b; bestSplit=split; bestDiff=0;
      break;
    }
    if(a.parent===b){
      const edge=Math.max(da-db, 0);
      const split=edge>0 ? (half-db)/edge : 0.5;
      bestNode=a; bestSplit=split; bestDiff=0;
      break;
    }
  }
  if(bestNode===null){
    fullPath.forEach((n, i)=>{
      if(n===rootNode) return; // skip existing root
      const diff=Math.abs(dists[i]-half);
      if(diff<bestDiff){bestDiff=diff;bestNode=n;bestSplit=0.5;}
    });
  }
  if(bestNode&&bestNode.parent) rerootAtNode(bestNode, bestSplit);
  _pvmMidpointApplied=true;
  document.getElementById("possvm-result").textContent="Midpoint root applied.";
}

// ── POSSVM species-overlap OG assignment ──────────────────────────────────────
function cloneCurrentTreeForPvm(){
  return d3.hierarchy(_h2d(rootNode), d=>d.children||null);
}

function getCurrentGeneSpeciesMap(){
  const m=new Map();
  if(!rootNode) return m;
  rootNode.leaves().forEach(l=>{
    const gid=l.data.gene_id||l.data.name||"";
    if(gid) m.set(gid, l.data.species||"");
  });
  return m;
}

function pruneTreeToGeneSet(node, geneSet){
  if(!node) return null;
  const ch=node.children||node._children||null;
  if(!ch||!ch.length){
    const gid=node.data ? (node.data.gene_id||node.data.name||"") : "";
    return geneSet.has(gid) ? node : null;
  }
  const kept=ch.map(c=>pruneTreeToGeneSet(c, geneSet)).filter(Boolean);
  if(!kept.length) return null;
  node.children=kept;
  delete node._children;
  kept.forEach(k=>{ k.parent=node; });
  if(kept.length===1){
    const child=kept[0];
    child.parent=node.parent||null;
    return child;
  }
  return node;
}

function cloneTreeDict(treeDict){
  return JSON.parse(JSON.stringify(treeDict));
}

function hierarchyFromTreeData(treeData){
  if(!treeData) return null;
  // Accept either a plain tree dict or an existing d3.hierarchy node.
  if(typeof treeData.each==="function") return d3.hierarchy(_h2d(treeData), d=>d.children||null);
  return d3.hierarchy(cloneTreeDict(treeData), d=>d.children||null);
}

function treeDictForGeneSet(treeRoot, genes){
  const geneSet=genes instanceof Set ? genes : new Set(genes||[]);
  if(!treeRoot || !geneSet.size) return null;
  const cloned=hierarchyFromTreeData(treeRoot);
  const pruned=pruneTreeToGeneSet(cloned, geneSet);
  return pruned ? _h2d(pruned) : null;
}

function computePossvmAssignments(treeRoot, ingroupSps, sos){
  function naturalCmp(a, b){
    return String(a).localeCompare(String(b), undefined, {numeric:true, sensitivity:"base"});
  }

  function collect(node){
    const ch=node.children||node._children;
    if(!ch||!ch.length){
      const gid=node.data.gene_id||node.data.name||"";
      const sp=node.data.species||"";
      node._pvmAllSps=sp ? new Set([sp]) : new Set();
      node._pvmAllLeaves=gid ? [gid] : [];
      node._pvmInLeaves=(gid && ingroupSps.has(sp)) ? [gid] : [];
      node._pvmEvent=null;
      return;
    }
    node._pvmAllSps=new Set();
    node._pvmAllLeaves=[];
    node._pvmInLeaves=[];
    for(const c of ch){
      collect(c);
      c._pvmAllSps.forEach(s=>node._pvmAllSps.add(s));
      if(c._pvmAllLeaves&&c._pvmAllLeaves.length) node._pvmAllLeaves.push(...c._pvmAllLeaves);
      if(c._pvmInLeaves&&c._pvmInLeaves.length) node._pvmInLeaves.push(...c._pvmInLeaves);
    }
  }

  function classify(node){
    const ch=node.children||node._children;
    if(!ch||ch.length<2){ node._pvmEvent=null; return; }
    for(const c of ch) classify(c);
    let isD=false;
    outer: for(let i=0;i<ch.length&&!isD;i++){
      for(let j=i+1;j<ch.length&&!isD;j++){
        const si=ch[i]._pvmAllSps, sj=ch[j]._pvmAllSps;
        if(!si.size||!sj.size) continue;
        let ovl=0;
        si.forEach(s=>{ if(sj.has(s)) ovl++; });
        if(ovl/Math.min(si.size, sj.size) > sos){ isD=true; break outer; }
      }
    }
    if(isD && node._pvmAllSps.size===1) isD=false;
    node._pvmEvent=isD?"D":"S";
  }

  function isOrthologyEvent(node){
    if(node._pvmEvent==="S") return true;
    if(node._pvmEvent!=="D") return false;
    const ch=node.children||node._children;
    if(!ch||ch.length!==2) return false;
    const left=[...ch[0]._pvmAllSps];
    const right=[...ch[1]._pvmAllSps];
    return left.length===1 && right.length===1 && left[0]===right[0];
  }

  const adjacency=new Map();
  function ensureNode(gid){
    if(gid && !adjacency.has(gid)) adjacency.set(gid, new Set());
  }
  function addEdge(a, b){
    if(!a||!b||a===b) return;
    ensureNode(a); ensureNode(b);
    adjacency.get(a).add(b);
    adjacency.get(b).add(a);
  }

  function buildOrthologyGraph(node){
    const ch=node.children||node._children;
    if(!ch||ch.length<2) return;
    for(const c of ch) buildOrthologyGraph(c);
    if(!isOrthologyEvent(node)) return;
    for(let i=0;i<ch.length;i++){
      for(let j=i+1;j<ch.length;j++){
        const left=ch[i]._pvmInLeaves||[];
        const right=ch[j]._pvmInLeaves||[];
        if(!left.length||!right.length) continue;
        for(const a of left){
          for(const b of right){
            addEdge(a,b);
          }
        }
      }
    }
  }

  function seededOrder(nodes){
    function key(s){
      let h=11;
      const str=String(s);
      for(let i=0;i<str.length;i++) h=((h*131)+str.charCodeAt(i))>>>0;
      return h;
    }
    return [...nodes].sort((a,b)=>key(a)-key(b) || naturalCmp(a,b));
  }

  function lpaClusters(nodes){
    if(!nodes.length) return {};
    const labels=new Map(nodes.map(gid=>[gid,gid]));
    const baseOrder=seededOrder(nodes);
    const maxIter=Math.max(25, nodes.length*3);
    for(let iter=0; iter<maxIter; iter++){
      let changed=false;
      const order=(iter%2===0) ? baseOrder : [...baseOrder].reverse();
      for(const gid of order){
        const neigh=[...(adjacency.get(gid)||[])];
        if(!neigh.length) continue;
        const counts=new Map();
        for(const nb of neigh){
          const lbl=labels.get(nb)||nb;
          counts.set(lbl, (counts.get(lbl)||0)+1);
        }
        let best=labels.get(gid)||gid;
        let bestCount=-1;
        counts.forEach((count, lbl)=>{
          if(count>bestCount || (count===bestCount && naturalCmp(lbl, best)<0)){
            best=lbl;
            bestCount=count;
          }
        });
        if(best!==labels.get(gid)){
          labels.set(gid, best);
          changed=true;
        }
      }
      if(!changed) break;
    }
    const groups=new Map();
    for(const gid of nodes){
      const lbl=labels.get(gid)||gid;
      if(!groups.has(lbl)) groups.set(lbl, []);
      groups.get(lbl).push(gid);
    }
    const ordered=[...groups.values()].map(genes=>genes.sort(naturalCmp))
      .sort((a,b)=>b.length-a.length || naturalCmp(a[0]||"", b[0]||""));
    const ogs={};
    ordered.forEach((genes, idx)=>{
      ogs["OG_"+String(idx+1).padStart(4,"0")]=genes;
    });
    return ogs;
  }

  collect(treeRoot);
  classify(treeRoot);
  const ingroupGenes=(treeRoot._pvmInLeaves||[]).slice().sort(naturalCmp);
  ingroupGenes.forEach(ensureNode);
  buildOrthologyGraph(treeRoot);
  return {ogs:lpaClusters(ingroupGenes), root:treeRoot};
}

function applyPossvmRun(newOgs, ingroupSps){
  const namedOgs=relabelOgMapWithReferences(newOgs);
  const nOgs=Object.keys(namedOgs).length;
  ogLeaf2Color={}; ogName2Color={}; ogGene2Name={};
  Object.keys(namedOgs).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true, sensitivity:"base"})).forEach((og,i)=>{
    const col=palette[i%palette.length];
    ogName2Color[og]=col;
    for(const gid of namedOgs[og]){ogLeaf2Color[gid]=col;ogGene2Name[gid]=og;}
  });
  _pvmActive=true;
  _pvmIngroupSps=ingroupSps;
  _pvmOgs=namedOgs;
  colorMode="og";
  // Sync the colour-by dropdown so the UI reflects the active mode
  const _cbs=document.getElementById("color-by");
  if(_cbs) _cbs.value="og";

  // ── Auto-highlight each POSSVM OG clade ──────────────────────────────────────
  // Remove any previous OG-shading highlights, then re-add one per POSSVM OG.
  cladeHighlights.forEach((rec,uid)=>{ if(rec._fromOgHl) cladeHighlights.delete(uid); });
  // Fully expand the tree first so MRCA lookup works on all leaves.
  rootNode.each(n=>{if(n._children){n.children=n._children;n._children=null;n._isOgCol=false;}});
  rootNode.each(n=>{n._uid=n._uid||0;}); // ensure _uid is set (set during renderTree; defensive)
  Object.keys(namedOgs).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true, sensitivity:"base"})).forEach((og,i)=>{
    const leaves=namedOgs[og].map(gid=>{
      let found=null;
      rootNode.each(n=>{if(n.data.leaf&&(n.data.gene_id||n.data.name)===gid) found=n;});
      return found;
    }).filter(Boolean);
    if(!leaves.length) return;
    const mrca=findExactOgRoot(leaves);
    if(!mrca) return;
    cladeHighlights.set(mrca._uid,{
      color:ogBaseColor(og, i),
      label:og,
      subtitle:ogHighlightSubtitle(mrca),
      _fromOgHl:true
    });
  });
  _ogHlActive=true;
  document.getElementById("btn-highlight-ogs").classList.add("active-btn");

  document.getElementById("n-ogs-label").textContent=nOgs+" orthogroups (POSSVM)";
  document.getElementById("possvm-result").textContent="\u2714 "+nOgs+" OG"+(nOgs!==1?"s":"")+" found";
  document.getElementById("possvm-reset-btn").style.display="inline";
  renderTree(false);
}

function runPossvm(){
  if(!rootNode){document.getElementById("possvm-result").textContent="No tree loaded."; return;}
  const ingroupSps=new Set([...document.querySelectorAll("[data-pvm-sp]:checked")].map(cb=>cb.value));
  if(!ingroupSps.size){document.getElementById("possvm-result").textContent="Select at least one species."; return;}
  // Auto-apply midpoint root if neither the user nor pvmMidpointRoot has rerooted this tree.
  // POSSVM D/S classification is root-dependent; midpoint root gives the most neutral starting point.
  if(!_pvmMidpointApplied&&!_userRerooted) pvmMidpointRoot();
  const sos=parseFloat(document.getElementById("possvm-sos").value);
  const result=computePossvmAssignments(rootNode, ingroupSps, sos);
  applyPossvmRun(result.ogs, ingroupSps);
}

function buildHierPossvmModel(levelNames, sos){
  const clades=getNamedSpeciesTreeClades();
  const byName=new Map(clades.map(rec=>[rec.name, rec]));
  const gene2sp=getCurrentGeneSpeciesMap();
  const currentOgs=sourceOgsForCurrentTree();
  const seen=new Set();
  const selected=levelNames.map(name=>byName.get(name)).filter(Boolean).filter(rec=>{
    const key=rec.species.slice().sort().join("|");
    if(seen.has(key)) return false;
    seen.add(key);
    return true;
  });
  selected.sort((a,b)=>b.species.length-a.species.length || a.depth-b.depth || a.name.localeCompare(b.name));
  for(let i=1;i<selected.length;i++){
    const prev=new Set(selected[i-1].species);
    const cur=selected[i].species;
    if(!cur.every(sp=>prev.has(sp))){
      throw new Error(`Selected hOG clades are not nested: ${selected[i-1].name} -> ${selected[i].name}`);
    }
  }

  const levels=[];
  const links=[];
  if(!selected.length) return {levels, links, sos};

  const rootRec=selected[0];
  const rootSpecies=new Set(rootRec.species);
  const rootOgs=Object.entries(currentOgs)
    .map(([og, genes])=>[og, genes.filter(gid=>rootSpecies.has(gene2sp.get(gid)))])
    .filter(([, genes])=>genes.length)
    .map(([og, genes])=>{
    const species=[...new Set(genes.map(gid=>gene2sp.get(gid)).filter(Boolean))].sort((a,b)=>speciesOrder.indexOf(a)-speciesOrder.indexOf(b));
    const rawId=og;
    return {
      id:`${rootRec.name}::${og}`,
      og:hogLabelFromGenes(og, genes),
      raw_id:rawId,
      emergent:true,
      clade:rootRec.name,
      depth:rootRec.depth,
      size:genes.length,
      genes:[...genes],
      species,
      subtitle:spMRCAName(new Set(species)) || rootRec.name,
      parent_id:null,
      tree_dict:treeDictForGeneSet(rootNode, genes),
    };
  }).sort((a,b)=>b.size-a.size || a.og.localeCompare(b.og));
  levels.push({
    name:rootRec.name,
    depth:rootRec.depth,
    species:[...rootRec.species],
    species_count:rootRec.species.length,
    ogs:rootOgs,
  });

  for(let li=1; li<selected.length; li++){
    const rec=selected[li];
    const recSpecies=new Set(rec.species);
    const prevLevel=levels[li-1];
    const nextOgs=[];
    prevLevel.ogs.forEach(parentOg=>{
      const keepGenes=parentOg.genes.filter(gid=>recSpecies.has(gene2sp.get(gid)));
      if(!keepGenes.length) return;
      if(keepGenes.length===1){
        const gid=keepGenes[0];
        const species=[gene2sp.get(gid)].filter(Boolean);
        const singletonTree=parentOg.tree_dict ? treeDictForGeneSet(parentOg.tree_dict, [gid]) : null;
        nextOgs.push({
          id:`${rec.name}::${parentOg.id}::1`,
          og:parentOg.og,
          emergent:(parentOg.genes||[]).length>1,
          clade:rec.name,
          depth:rec.depth,
          size:1,
          genes:[gid],
          species,
          subtitle:spMRCAName(new Set(species)) || rec.name,
          parent_id:parentOg.id,
          raw_id:parentOg.raw_id,
          tree_dict:singletonTree || parentOg.tree_dict,
        });
        links.push({
          source:parentOg.id,
          target:`${rec.name}::${parentOg.id}::1`,
          source_level:li-1,
          target_level:li,
          source_og:parentOg.og,
          target_og:parentOg.og,
          weight:1,
          genes:[gid],
        });
        return;
      }
      const parentTree=parentOg.tree_dict ? hierarchyFromTreeData(parentOg.tree_dict) : cloneCurrentTreeForPvm();
      const pruned=pruneTreeToGeneSet(parentTree, new Set(keepGenes));
      if(!pruned) return;
      const result=computePossvmAssignments(pruned, recSpecies, sos);
      const entries=Object.entries(result.ogs).map(([og, genes])=>[og, genes.filter(gid=>keepGenes.includes(gid))]).filter(([, genes])=>genes.length);
      entries.sort((a,b)=>b[1].length-a[1].length || a[0].localeCompare(b[0]));
      const usedChildLabels=new Set();
      entries.forEach(([og, genes], idx)=>{
        const species=[...new Set(genes.map(gid=>gene2sp.get(gid)).filter(Boolean))].sort((a,b)=>speciesOrder.indexOf(a)-speciesOrder.indexOf(b));
        // Always base child labels on the parent OG's raw_id so they inherit
        // the pipeline OG context.  For splits, hogLabelFromGenes appends the
        // reference-gene content of each sub-cluster (e.g. OG_0001:like:Pax5).
        // Deduplicate within the same parent in case two children share refs.
        const rawId=parentOg.raw_id;
        const baseLabel=entries.length===1 ? parentOg.og : hogLabelFromGenes(rawId, genes);
        let label=baseLabel; let sfx=2;
        while(usedChildLabels.has(label)){ label=`${baseLabel}#${sfx++}`; }
        usedChildLabels.add(label);
        const childId=`${rec.name}::${parentOg.id}::${idx+1}`;
        nextOgs.push({
          id:childId,
          og:label,
          raw_id:rawId,
          raw_og:og,
          emergent:entries.length>1,
          clade:rec.name,
          depth:rec.depth,
          size:genes.length,
          genes:[...genes],
          species,
          subtitle:spMRCAName(new Set(species)) || rec.name,
          parent_id:parentOg.id,
          tree_dict:treeDictForGeneSet(result.root, genes),
        });
        links.push({
          source:parentOg.id,
          target:childId,
          source_level:li-1,
          target_level:li,
          source_og:parentOg.og,
          target_og:label,
          weight:genes.length,
          genes:[...genes],
        });
      });
    });
    levels.push({
      name:rec.name,
      depth:rec.depth,
      species:[...rec.species],
      species_count:rec.species.length,
      ogs:nextOgs.sort((a,b)=>b.size-a.size || a.og.localeCompare(b.og)),
    });
  }
  levels.forEach((level, li)=>level.ogs.forEach(og=>{ og._levelIndex=li; }));

  // ── Keep only OG chains that span ALL selected levels ────────────────────
  // Back-propagate from the last level: an OG "reaches the end" if it is at
  // the last level, or if at least one of its children reaches the end.
  // This lets the user see only meaningful duplication/split patterns that
  // are traceable across every clade they selected.
  if(levels.length>1){
    const childrenOf=new Map();
    links.forEach(link=>{
      if(!childrenOf.has(link.source)) childrenOf.set(link.source,[]);
      childrenOf.get(link.source).push(link.target);
    });
    const reachesEnd=new Set();
    levels[levels.length-1].ogs.forEach(og=>reachesEnd.add(og.id));
    for(let li=levels.length-2;li>=0;li--){
      levels[li].ogs.forEach(og=>{
        if((childrenOf.get(og.id)||[]).some(cid=>reachesEnd.has(cid))) reachesEnd.add(og.id);
      });
    }
    levels.forEach(level=>{ level.ogs=level.ogs.filter(og=>reachesEnd.has(og.id)); });
    const kept=links.filter(l=>reachesEnd.has(l.source)&&reachesEnd.has(l.target));
    links.length=0; kept.forEach(l=>links.push(l));
  }
  // ── End span-all-levels filter ────────────────────────────────────────────

  const allById=new Map();
  levels.forEach(level=>level.ogs.forEach(og=>allById.set(og.id, og)));
  const bySource=new Map();
  links.forEach(link=>{
    if(!bySource.has(link.source)) bySource.set(link.source, []);
    bySource.get(link.source).push(link);
  });

  const visibleLevels=levels.map((level, li)=>{
    const showAllForLevel=(li===0 || li===levels.length-1);
    let ogs=level.ogs.filter(og=>showAllForLevel || og.emergent);
    // If a selected intermediate level only carries parent OGs forward and
    // does not introduce a new split, still render those carried-through OGs
    // so the level does not misleadingly appear as "0 OGs".
    if(li>0 && !ogs.length && level.ogs.length){
      ogs=level.ogs.map(og=>({ ...og, carried:true }));
    } else {
      ogs=ogs.map(og=>({ ...og, carried:!!og.carried }));
    }
    return {
      ...level,
      ogs,
    };
  });
  const visibleSet=new Set();
  visibleLevels.forEach(level=>level.ogs.forEach(og=>visibleSet.add(og.id)));

  function collectVisibleTargets(sourceId){
    const out=bySource.get(sourceId)||[];
    const acc=[];
    for(const link of out){
      const child=allById.get(link.target);
      if(!child) continue;
      if(visibleSet.has(child.id)){
        acc.push({
          target:child,
          weight:child.size,
          genes:[...(child.genes||[])],
        });
      } else {
        acc.push(...collectVisibleTargets(child.id));
      }
    }
    return acc;
  }

  const visibleLinks=[];
  visibleLevels.forEach((level, li)=>{
    level.ogs.forEach(sourceOg=>{
      const targets=collectVisibleTargets(sourceOg.id);
      targets.forEach(({target, weight, genes})=>{
        visibleLinks.push({
          source:sourceOg.id,
          target:target.id,
          source_level:li,
          target_level:target._levelIndex,
          source_og:sourceOg.og,
          target_og:target.og,
          weight,
          genes,
        });
      });
    });
  });

  return {levels:visibleLevels, links:visibleLinks, sos};
}

function hogChildSortForParent(parent, children){
  return children.slice().sort((a,b)=>{
    const aCarry=(a.og.raw_id&&parent.raw_id&&a.og.raw_id===parent.raw_id) || a.og.og===parent.og;
    const bCarry=(b.og.raw_id&&parent.raw_id&&b.og.raw_id===parent.raw_id) || b.og.og===parent.og;
    if(aCarry!==bCarry) return aCarry ? -1 : 1;
    return b.og.size-a.og.size || a.og.og.localeCompare(b.og.og) || a.idx-b.idx;
  });
}

function assignHogOrders(levels, links){
  const incomingByTarget=new Map();
  links.forEach(link=>{
    if(!incomingByTarget.has(link.target)) incomingByTarget.set(link.target, []);
    incomingByTarget.get(link.target).push(link);
  });
  levels.forEach((level, li)=>{
    if(li===0){
      level.ogs.forEach((og, idx)=>{ og._order=idx; });
      return;
    }
    const prevLevel=levels[li-1];
    const childBuckets=new Map();
    const orphans=[];
    level.ogs.forEach((og, idx)=>{
      const inbound=(incomingByTarget.get(og.id)||[]).filter(link=>link.target_level===li);
      if(!inbound.length){
        orphans.push({og, idx});
        return;
      }
      const primary=inbound.slice().sort((a,b)=>b.weight-a.weight || a.source.localeCompare(b.source))[0];
      if(!childBuckets.has(primary.source)) childBuckets.set(primary.source, []);
      childBuckets.get(primary.source).push({og, idx, parentId:primary.source});
    });

    const assignments=[];
    let nextFree=0;
    prevLevel.ogs
      .slice()
      .sort((a,b)=>(a._order??0)-(b._order??0) || a.og?.localeCompare?.(b.og||"") || 0)
      .forEach(parent=>{
        const children=hogChildSortForParent(parent, childBuckets.get(parent.id)||[]);
        if(!children.length) return;
        const blockLen=children.length;
        const anchor=Math.max(0, Math.round(parent._order ?? nextFree));
        const idealStart=Math.max(0, Math.round(anchor - (blockLen-1)/2));
        const start=Math.max(nextFree, idealStart);
        const positions=Array.from({length:blockLen}, (_,i)=>start+i);
        nextFree=start+blockLen;

        const carryIdx=children.findIndex(ch=>
          (ch.og.raw_id&&parent.raw_id&&ch.og.raw_id===parent.raw_id) || ch.og.og===parent.og
        );
        const centerOrder=positions
          .map((pos,i)=>({pos,i,diff:Math.abs(pos-anchor)}))
          .sort((a,b)=>a.diff-b.diff || a.pos-b.pos);
        const placed=new Array(blockLen);
        const usedPosIdx=new Set();
        const usedChildIdx=new Set();
        if(carryIdx>=0){
          const slot=centerOrder[0].i;
          placed[slot]=children[carryIdx];
          usedPosIdx.add(slot);
          usedChildIdx.add(carryIdx);
        }
        const remainingChildren=children
          .map((child,i)=>({child,i}))
          .filter(rec=>!usedChildIdx.has(rec.i))
          .map(rec=>rec.child);
        const remainingSlots=centerOrder
          .map(rec=>rec.i)
          .filter(i=>!usedPosIdx.has(i));
        remainingSlots.forEach((slot, idx)=>{
          placed[slot]=remainingChildren[idx];
        });
        positions.forEach((pos, idx)=>{
          const child=placed[idx];
          if(!child) return;
          child.og._order=pos;
          assignments.push(child);
        });
      });

    orphans.sort((a,b)=>a.idx-b.idx).forEach(orphan=>{
      orphan.og._order=nextFree++;
      assignments.push(orphan);
    });
    assignments.sort((a,b)=>a.og._order-b.og._order || a.idx-b.idx);
    level.ogs=assignments.map(rec=>rec.og);
  });
}

function closeHogMapPanel(){
  const panel=document.getElementById("hog-map-panel");
  if(panel) panel.style.display="none";
  if(window._hogClickGenes){
    window._hogClickGenes=null;
    cladeHighlights.forEach((rec,uid)=>{ if(rec._fromHogClick) cladeHighlights.delete(uid); });
    if(rootNode) renderTree(false);
  }
}

let _hogPanelDragInit=false;
function initHogMapPanelDrag(){
  if(_hogPanelDragInit) return;
  const panel=document.getElementById("hog-map-panel");
  const header=document.getElementById("hog-map-header");
  if(!panel || !header) return;
  _hogPanelDragInit=true;
  let dragging=false, startX=0, startY=0, baseLeft=0, baseTop=0;
  header.addEventListener("mousedown",(ev)=>{
    if(ev.target.closest("button")) return;
    const rect=panel.getBoundingClientRect();
    dragging=true;
    startX=ev.clientX;
    startY=ev.clientY;
    baseLeft=rect.left;
    baseTop=rect.top;
    panel.style.left=rect.left+"px";
    panel.style.top=rect.top+"px";
    panel.style.right="auto";
    ev.preventDefault();
  });
  window.addEventListener("mousemove",(ev)=>{
    if(!dragging) return;
    const panelRect=panel.getBoundingClientRect();
    const nextLeft=Math.min(
      Math.max(8, baseLeft + (ev.clientX-startX)),
      Math.max(8, window.innerWidth - panelRect.width - 8)
    );
    const nextTop=Math.min(
      Math.max(8, baseTop + (ev.clientY-startY)),
      Math.max(8, window.innerHeight - panelRect.height - 8)
    );
    panel.style.left=nextLeft+"px";
    panel.style.top=nextTop+"px";
  });
  window.addEventListener("mouseup",()=>{ dragging=false; });

  // ── Resize handle ──────────────────────────────────────────────────────────
  const resizeHandle=document.getElementById("hog-map-resize");
  if(resizeHandle){
    let resizing=false, rStartX=0, rStartY=0, rBaseW=0, rBaseH=0;
    resizeHandle.addEventListener("mousedown",(ev)=>{
      const rect=panel.getBoundingClientRect();
      // Pin position so the panel doesn't jump when we switch from right: to left:
      panel.style.left=rect.left+"px";
      panel.style.top=rect.top+"px";
      panel.style.right="auto";
      resizing=true;
      rStartX=ev.clientX; rStartY=ev.clientY;
      rBaseW=rect.width;   rBaseH=rect.height;
      ev.preventDefault(); ev.stopPropagation();
    });
    window.addEventListener("mousemove",(ev)=>{
      if(!resizing) return;
      panel.style.width =Math.max(420, rBaseW+(ev.clientX-rStartX))+"px";
      panel.style.height=Math.max(260, rBaseH+(ev.clientY-rStartY))+"px";
    });
    window.addEventListener("mouseup",()=>{ resizing=false; });
  }
}

function highlightHogClade(og, color){
  if(!rootNode || !og || !og.genes || !og.genes.length) return;
  // Clear previous hOG-click clade backgrounds
  cladeHighlights.forEach((rec,uid)=>{ if(rec._fromHogClick) cladeHighlights.delete(uid); });
  // Clear previous hOG-click leaf coloring
  if(window._hogClickGenes){
    window._hogClickGenes.forEach(gid=>{ delete ogLeaf2Color[gid]; delete ogGene2Name[gid]; });
  }
  delete ogName2Color[window._hogClickOgName];

  // Register this OG's genes in the leaf-color maps so they light up
  // in the gene tree while all other leaves dim to the default grey.
  const fill=color||"#4a90d9";
  window._hogClickGenes=new Set(og.genes);
  window._hogClickOgName=og.og;
  og.genes.forEach(gid=>{ ogLeaf2Color[gid]=fill; ogGene2Name[gid]=og.og; });
  ogName2Color[og.og]=fill;
  colorMode="og";
  const cbs=document.getElementById("color-by"); if(cbs) cbs.value="og";

  // Expand any collapsed nodes so MRCA lookup works
  rootNode.each(d=>{ if(d._children){ d.children=d._children; d._children=null; d._isOgCol=false; } });
  const want=new Set(og.genes);
  const leaves=[];
  rootNode.each(d=>{
    if(d.data && d.data.leaf){
      const gid=d.data.gene_id||d.data.name||"";
      if(want.has(gid)) leaves.push(d);
    }
  });
  if(!leaves.length){ renderTree(false); return; }
  const mrca=findExactOgRoot(leaves);
  if(mrca){
    cladeHighlights.set(mrca._uid,{
      color: fill,
      label: og.og,
      subtitle: ogHighlightSubtitle(mrca),
      _fromHogClick:true,
    });
  }
  renderTree(false);
}

function renderHogMap(model){
  initHogMapPanelDrag();
  const wrap=document.getElementById("hog-map-wrap");
  const subtitle=document.getElementById("hog-map-subtitle");
  wrap.innerHTML="";
  if(!model||!model.levels.length){
    wrap.innerHTML='<div style="padding:18px;color:#7c8b95;font-size:12px">No hOG data to display.</div>';
    subtitle.textContent="";
    return;
  }
  const colGap=210, rowGap=54, pad={top:48,left:42,right:240,bottom:30};
  const maxRows=Math.max(...model.levels.map(l=>Math.max(1,l.ogs.length)));
  const width=pad.left+pad.right+Math.max(1, model.levels.length-1)*colGap+120;
  const height=pad.top+pad.bottom+Math.max(1,maxRows-1)*rowGap+70;
  const linkMap=new Map();
  model.links.forEach(link=>{
    const key=`${link.source_level}:${link.target_level}:${link.target}`;
    if(!linkMap.has(key)) linkMap.set(key, []);
    linkMap.get(key).push(link);
  });
  assignHogOrders(model.levels, model.links);

  model.levels.forEach((level, li)=>{
    level.x=pad.left + li*colGap;
    level.ogs.forEach((og)=>{ og.x=level.x; og.y=pad.top + og._order*rowGap; });
  });

  const colorScale=d3.scaleOrdinal(palette.concat(_pvmOgHlPalette));
  const id2og=new Map();
  model.levels.forEach(level=>level.ogs.forEach(og=>id2og.set(og.id, og)));
  const svg=d3.select(wrap).append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("display","block");

  const axis=svg.append("g");
  model.levels.forEach(level=>{
    axis.append("text")
      .attr("x", level.x)
      .attr("y", 18)
      .attr("text-anchor","middle")
      .attr("font-size", 12)
      .attr("font-weight", 700)
      .attr("fill", "#314553")
      .text(level.name);
    axis.append("text")
      .attr("x", level.x)
      .attr("y", 33)
      .attr("text-anchor","middle")
      .attr("font-size", 10)
      .attr("fill", "#7c8b95")
      .text(`${level.ogs.length} OG${level.ogs.length!==1?"s":""} · ${level.species_count} spp`);
  });

  const relatedMap=new Map();
  function relAdd(a,b){
    if(!relatedMap.has(a)) relatedMap.set(a, new Set([a]));
    if(!relatedMap.has(b)) relatedMap.set(b, new Set([b]));
    relatedMap.get(a).add(b);
    relatedMap.get(b).add(a);
  }
  model.levels.forEach(level=>level.ogs.forEach(og=>{
    if(!relatedMap.has(og.id)) relatedMap.set(og.id, new Set([og.id]));
  }));
  model.links.forEach(link=>relAdd(link.source, link.target));

  const linkG=svg.append("g").attr("fill","none");
  const linkSel=linkG.selectAll(".hog-link")
    .data(model.links)
    .enter()
    .append("path")
      .attr("class","hog-link")
      .attr("data-source", link=>link.source)
      .attr("data-target", link=>link.target);
  linkSel.each(function(link){
    const s=id2og.get(link.source), t=id2og.get(link.target);
    if(!s||!t){ d3.select(this).remove(); return; }
    const x1=s.x+9, x2=t.x-9, y1=s.y, y2=t.y;
    const dx=Math.max(30, (x2-x1)*0.45);
    d3.select(this)
      .attr("d", `M${x1},${y1}C${x1+dx},${y1} ${x2-dx},${y2} ${x2},${y2}`)
      .attr("stroke", colorScale(link.source_og))
      .attr("stroke-width", Math.max(1.5, Math.sqrt(link.weight)))
      .attr("stroke-opacity", 0.55)
      .append("title")
      .text(`${link.source_og} → ${link.target_og} (${link.weight} genes)`);
  });

  const nodeG=svg.append("g");
  const nodeSel=nodeG.selectAll(".hog-node")
    .data(model.levels.flatMap(level=>level.ogs.map(og=>({level, og}))))
    .enter()
    .append("g")
      .attr("class","hog-node")
      .attr("data-id", d=>d.og.id)
      .attr("transform", d=>`translate(${d.og.x},${d.og.y})`);

  function hogTooltipHtml(level, og){
    const refSummary=hogReferenceSummary(og.genes);
    const speciesList=og.species.length ? og.species.join(", ") : "NA";
    const genePreview=og.genes.slice(0,6).join(", ");
    const moreGenes=og.genes.length>6 ? `, +${og.genes.length-6} more` : "";
    let html='<div class="tt-name">'+og.og+'</div>';
    html += '<div class="tt-row"><span>Level</span><strong>'+level.name+'</strong></div>';
    html += '<div class="tt-row"><span>Genes</span><strong>'+og.size+'</strong></div>';
    html += '<div class="tt-row"><span>Species</span><strong>'+og.species.length+'</strong></div>';
    if(og.subtitle) html += '<div class="tt-row"><span>MRCA</span><strong>'+og.subtitle+'</strong></div>';
    if(refSummary.direct.length){
      html += '<div class="tt-row"><span>Reference genes</span><strong style="color:#8e44ad">'+refSummary.direct.join("/")+'</strong></div>';
    } else if(refSummary.inferred.length){
      html += '<div class="tt-row"><span>Reference-like</span><strong>'+refSummary.inferred.join("/")+'</strong></div>';
    }
    html += '<div style="margin-top:4px;font-size:10px;color:#666;line-height:1.45"><div><strong>Species:</strong> '+speciesList+'</div><div><strong>Genes:</strong> '+genePreview+moreGenes+'</div></div>';
    return html;
  }

  function clearHogHover(){
    linkSel
      .attr("stroke-opacity", 0.55)
      .attr("stroke-width", link=>Math.max(1.5, Math.sqrt(link.weight)));
    nodeSel.style("opacity", 1);
    nodeSel.select("circle")
      .attr("stroke", "#fff")
      .attr("stroke-width", 1.5);
    nodeSel.selectAll("text")
      .attr("opacity", 1);
    hideTip();
  }

  function applyHogHover(level, og, event){
    const related=relatedMap.get(og.id) || new Set([og.id]);
    nodeSel.style("opacity", d=>related.has(d.og.id)?1:0.22);
    nodeSel.select("circle")
      .attr("stroke", d=>d.og.id===og.id?"#111":"#fff")
      .attr("stroke-width", d=>d.og.id===og.id?3:1.5);
    nodeSel.selectAll("text")
      .attr("opacity", d=>related.has(d.og.id)?1:0.28);
    linkSel
      .attr("stroke-opacity", link=>(link.source===og.id||link.target===og.id||(related.has(link.source)&&related.has(link.target)))?0.95:0.08)
      .attr("stroke-width", link=>{
        const base=Math.max(1.5, Math.sqrt(link.weight));
        return (link.source===og.id||link.target===og.id)?base+1.3:base;
      });
    showTip(event, hogTooltipHtml(level, og));
  }

  model.levels.forEach(level=>{
    nodeSel.filter(d=>d.level===level).each(function(d){
      const og=d.og;
      const g=d3.select(this);
      const fill=colorScale(og.og);
      g.append("circle")
        .attr("r", Math.max(4, Math.min(11, 3 + Math.sqrt(og.size))))
        .attr("fill", fill)
        .attr("stroke", "#fff")
        .attr("stroke-width", 1.5);
      g.append("text")
        .attr("x", 12)
        .attr("y", -2)
        .attr("font-size", 11)
        .attr("font-weight", 600)
        .attr("fill", "#314553")
        .text(`${og.og} [${og.size}]`);
      g.append("text")
        .attr("x", 12)
        .attr("y", 11)
        .attr("font-size", 9)
        .attr("fill", "#8a98a3")
        .text(og.subtitle || "");
      g.on("mouseover", function(event){ applyHogHover(level, og, event); })
       .on("mousemove", moveTip)
       .on("mouseout", clearHogHover)
       .on("click", function(event){
         event.stopPropagation();
         highlightHogClade(og, fill);
       });
    });
  });

  subtitle.textContent=`${model.levels.length} clade levels · broad to nested · SOS ${model.sos.toFixed(2)}`;
}

function runHierPossvm(){
  if(!rootNode){
    document.getElementById("hog-result").textContent="No tree loaded.";
    return;
  }
  if(_pvmHogNodes.length<2){
    document.getElementById("hog-result").textContent="Select at least two named clades.";
    return;
  }
  const sos=parseFloat(document.getElementById("possvm-sos").value);
  let model;
  try{
    model=buildHierPossvmModel(_pvmHogNodes, sos);
  }catch(err){
    document.getElementById("hog-result").textContent=String(err&&err.message?err.message:err);
    return;
  }
  _pvmHogModel=model;
  renderHogMap(model);
  document.getElementById("hog-map-panel").style.display="block";
  const totalOgs=model.levels.reduce((acc, level)=>acc+level.ogs.length, 0);
  document.getElementById("hog-result").textContent=`${model.levels.length} levels · ${totalOgs} hOG nodes`;
}

function pvmReset(){
  _pvmActive=false; _pvmIngroupSps=null; _pvmOgs=null;
  cladeHighlights.forEach((rec,uid)=>{ if(rec._fromOgHl) cladeHighlights.delete(uid); });
  _ogHlActive=false;
  document.getElementById("btn-highlight-ogs").classList.remove("active-btn");
  if(showDSNodes){ showDSNodes=false; document.getElementById("btn-ds-nodes").classList.remove("active-btn"); }
  colorMode="og";
  rebuildOgColors();
  const n=Object.keys(activeOgs()).length;
  document.getElementById("n-ogs-label").textContent=n+" orthogroups";
  document.getElementById("possvm-reset-btn").style.display="none";
  document.getElementById("possvm-result").textContent="";
  renderTree(false);
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

function drawGeneTree(treeData, opts){
  opts=opts||{};
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
  // reset OG collapse/highlight toggle states for the new tree
  _ogCollapseActive=false; document.getElementById("btn-collapse-ogs").classList.remove("active-btn");
  _ogHlActive=false;       document.getElementById("btn-highlight-ogs").classList.remove("active-btn");
  if(opts.preserveRerootState){
    document.getElementById("btn-reset-root").style.display="inline";
  } else {
    _isRerooted=false; _origTreeDictForReroot=null;
    _pvmMidpointApplied=false; _userRerooted=false;
    document.getElementById("btn-reset-root").style.display="none";
  }
  if(opts.preserveFocusState){
    document.getElementById("btn-reset-focus").style.display="inline";
  } else {
    _isSubtreeFocused=false; _origTreeDictForFocus=null;
    document.getElementById("btn-reset-focus").style.display="none";
  }
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

  const visibleLeafNodes=rootNode.descendants().filter(d=>d.data&&d.data.leaf);
  const leafLabelLayout=computeLeafLabelLayout(visibleLeafNodes);
  const visLeafRight=d3.max(
    visibleLeafNodes,
    d=>leafLabelRightPx(d, mg, leafLabelLayout)
  ) || (nodeX(rootNode,mg)+iWeff);
  const maxCladeLabelWidth=cladeHighlights.size
    ? d3.max(Array.from(cladeHighlights.values()), rec=>Math.max(
        measureTextPx((rec&&rec.label)||"", 16, "sans-serif", "600"),
        measureTextPx((rec&&rec.subtitle)||"", 12, "sans-serif", "400")
      ))
    : 0;
  treeSvg.attr("width", Math.max(W, visLeafRight + cladeHlExtend + 14 + (maxCladeLabelWidth||0) + 24));

  const dur=animate?240:0;

  // ── Clade highlight backgrounds ──
  {
    let hlLayer=gMain.select(".clade-hl-layer");
    if(hlLayer.empty()) hlLayer=gMain.insert("g",":first-child").attr("class","clade-hl-layer");
    hlLayer.lower();
    hlLayer.selectAll("*").remove();
    if(cladeHighlights.size){
      const allNodes=rootNode.descendants();
      // Extend the highlight box past the actual rendered leaf labels.
      const boxRight=visLeafRight+cladeHlExtend;
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
        const subtitle=rec.subtitle||"";
        hlLayer.append("rect")
          .attr("x",x1).attr("y",y1)
          .attr("width",Math.max(0,boxRight-x1)).attr("height",Math.max(0,y2-y1))
          .attr("fill",color).attr("opacity",cladeHlAlpha).attr("rx",4)
          .style("cursor","pointer")
          .on("click",(ev)=>{ ev.stopPropagation(); showCladeHlPopup(ev,uid); })
          .on("contextmenu",(ev)=>showCladeHlPopup(ev,uid));
        if(label){
          const labelFontSize=Math.max(11, Math.min(13, tipFontSVG()*1.05));
          const textY=(y1+y2)/2 - (subtitle ? labelFontSize*0.42 : 0);
          const txt=hlLayer.append("text")
            .attr("x", boxRight+8)
            .attr("y", textY)
            .attr("text-anchor","start")
            .attr("font-size", labelFontSize)
            .attr("font-weight","600")
            .attr("fill", "#111")
            .attr("opacity", 1)
            .style("pointer-events","none");
          txt.append("tspan").text(label);
          if(subtitle){
            txt.append("tspan")
              .attr("x", boxRight+8)
              .attr("dy","1.15em")
              .attr("font-size",Math.max(10, labelFontSize*0.86))
              .attr("font-weight","400")
              .attr("fill","#8a97a1")
              .text(subtitle);
          }
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
      if(showDSNodes&&_pvmActive&&d._pvmEvent) return tipFontSVG()*0.55;
      if(_pvmActive&&d._pvmEvent==="D") return tipFontSVG()*0.55;
      const isOGlbl=showOGLabels&&(isOGNode(d)||d.data._og_label);
      if(isOGlbl) return tipFontSVG()*0.55;
      if(isOGNode(d)) return tipFontSVG()*0.5;
      return tipFontSVG()*0.26;
    })
    .attr("fill",d=>{
      if(d.data.leaf){
        const gid2=d.data.gene_id||d.data.name;
        const og2=leafOgName(d);
        if(ogHlSet!==null){
          if(!ogHlSet.has(og2)) return "#e8e8e8";
          const gi=ogHlGroupIndex.get(og2)??0;
          return ogHlTagColor(gi);
        }
        if(colorMode==="og") return ogLeafColor(gid2, d.data.species);
        const c=leafColor(d.data.species||"");
        return hlSet&&!hlSet.has(d.data.species||"")?"#e8e8e8":c;
      }
      if(showDSNodes&&_pvmActive&&d._pvmEvent==="D") return "#e67e22";
      if(showDSNodes&&_pvmActive&&d._pvmEvent==="S") return "#27ae60";
      if(_pvmActive&&d._pvmEvent==="D") return "#111";
      if(showOGLabels&&(isOGNode(d)||d.data._og_label)) return "#222";
      if(isOGNode(d))return "#e74c3c";
      return "#ccc";
    })
    .attr("stroke",d=>{
      if(d.data.leaf) return "none";
      if(showDSNodes&&_pvmActive&&d._pvmEvent==="D") return "#c0392b";
      if(showDSNodes&&_pvmActive&&d._pvmEvent==="S") return "#1e8449";
      if(_pvmActive&&d._pvmEvent==="D") return "#000";
      if(showOGLabels&&(isOGNode(d)||d.data._og_label)) return "#000";
      if(isOGNode(d)) return "#b03a2e";
      return "#aaa";
    })
    .attr("stroke-width",d=>{
      if(d.data.leaf) return null;
      if(showDSNodes&&_pvmActive&&d._pvmEvent) return 1.5;
      if(_pvmActive&&d._pvmEvent==="D") return 2;
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
      if(!d.data.leaf&&d.children){
        const evLabel=_pvmActive&&d._pvmEvent?'<span style="font-size:9px;margin-left:4px;padding:1px 4px;border-radius:3px;background:'+(d._pvmEvent==='D'?'#111':'#2ecc71')+';color:#fff">'+d._pvmEvent+'</span>':'';
        showTip(event,'<b>'+(d.data.name||"internal")+'</b>'+evLabel+'<div style="font-size:9px;color:#aaa;margin-top:3px">click to choose collapse / reroot style</div>');
      }
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
    .attr("x",7).attr("y",0).attr("dy",null).attr("dominant-baseline","middle").attr("text-anchor","start")
    .attr("font-size",d=>d.data.leaf?treeLabelFontSVG():0)
    .attr("font-family","monospace")
    .style("cursor",()=>_compareNode1?"crosshair":"default")
    .on("click",(event,d)=>{
      // In compare mode, leaf labels also act as a click target for the second node
      if(_compareNode1!==null){ event.stopPropagation(); if(_compareNode1!==d) runSpeciesComparison(d); hideTip(); }
    })
    .attr("display",d=>(isLeafLabelVisible(d)?null:"none"))
    .text("")
    .each(function(d){
      if(!d.data.leaf) return;
      const el=d3.select(this);
      el.selectAll("tspan").remove();
      if(!isLeafLabelVisible(d)) return;
      const cols=leafLabelColumnXs(leafLabelLayout);
      const {gid, og, ref}=leafLabelParts(d);
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
      const hasOgText=showOGName&&!showOGLabels&&og;
      const hasRefText=showRefOrtho&&ref;
      if(leafLabelLayout.hasGid){
        if(showGeneId&&gid){
          el.append("tspan").attr("x",cols.gidX).attr("fill",baseCol).text(gid);
        }
      }
      if(leafLabelLayout.hasOg){
        if(leafLabelLayout.hasGid && (hasOgText || hasRefText)){
          el.append("tspan").attr("x",cols.ogSepX).attr("fill",sepCol).text(" \u00b7 ");
        }
        if(hasOgText){
          el.append("tspan").attr("x",cols.ogX).attr("fill",ogCol).text(og);
        }
      }
      if(leafLabelLayout.hasRef){
        if((leafLabelLayout.hasGid||leafLabelLayout.hasOg) && hasRefText){
          el.append("tspan").attr("x",cols.refSepX).attr("fill",sepCol).text(" \u00b7 ");
        }
        if(hasRefText){
          el.append("tspan").attr("x",cols.refX).attr("fill",refCol).text(ref);
        }
      }
    });

  // OG labels: beside OG-collapsed triangle, or beside expanded OG-named internal
  nodeSel.select(".og-label")
    .attr("x",d=>d._children?BADGE_W+6:-7)
    .attr("dy","0.35em")
    .attr("text-anchor",d=>d._children?"start":"end")
    .attr("font-size",treeLabelFontSVG())
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
            .attr("font-size",Math.max(7,treeLabelFontSVG()*0.82))
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

// ── reroot ──────────────────────────────────────────────────────────────────
let _origTreeDictForReroot=null;   // saved before the first reroot
let _isRerooted=false;
function resetFocus(){
  if(!_origTreeDictForFocus)return;
  const fullTree=_origTreeDictForFocus;
  _origTreeDictForFocus=null;
  _isSubtreeFocused=false;
  document.getElementById("btn-reset-focus").style.display="none";
  drawGeneTree(fullTree);
}
function resetRoot(){
  if(!_origTreeDictForReroot)return;
  _isRerooted=false;
  document.getElementById("btn-reset-root").style.display="none";
  drawGeneTree(_origTreeDictForReroot, {preserveFocusState:_isSubtreeFocused});
  _origTreeDictForReroot=null;
}

function focusOnTreeDict(treeDict){
  if(!treeDict||!rootNode) return;
  if(!_isSubtreeFocused) _origTreeDictForFocus=_h2d(rootNode);
  // Explicit subtree focus should show the focused clade itself, not inherit a
  // previous heatmap-driven gene filter that can hide all tip labels.
  hmFocusGids=null;
  _isSubtreeFocused=true;
  document.getElementById("btn-reset-focus").style.display="inline";
  drawGeneTree(cloneTreeDict(treeDict), {preserveFocusState:true});
}

function focusOnNode(d){
  if(!d||!rootNode) return;
  if(d.data&&d.data.leaf) return;
  focusOnTreeDict(_h2d(d));
}

// Convert a d3-hierarchy node to a plain dict.
// Explicitly strips d.data.children (original JSON) before setting from d3 hierarchy,
// avoiding stale branches being copied through Object.assign.
function _h2d(d){
  const ch=d.children||d._children;
  const obj=Object.assign({},d.data);
  delete obj.children;            // remove original JSON children – we'll set from d3
  if(ch&&ch.length) obj.children=ch.map(_h2d);
  return obj;
}

// Deep-copy a plain JSON-serialisable dict (avoids shared-reference mutations).
function _deepCopyDict(o){ return JSON.parse(JSON.stringify(o)); }

function _contractUnaryNodes(node, keepRoot){
  if(!node||!node.children||!node.children.length) return node;
  node.children=node.children.map(ch=>_contractUnaryNodes(ch,false));
  if(!keepRoot&&node.children.length===1){
    const child=node.children[0];
    child.dist=(Number(child.dist||0)+Number(node.dist||0));
    return child;
  }
  return node;
}

// Reroot a plain tree dict so that the node at `path` (array of child indices from root)
// becomes one of the two children of a new synthetic root.
// Algorithm: convert the tree into an undirected graph, split the selected edge,
// then rebuild a rooted tree away from the new synthetic root.
function _rerootDictAt(root, path, splitFrac){
  if(!path.length) return _deepCopyDict(root);
  if(splitFrac==null||!Number.isFinite(splitFrac)) splitFrac=0.5;
  splitFrac=Math.max(0, Math.min(1, splitFrac));

  const work=_deepCopyDict(root);
  let nextId=0;
  function assignIds(node){
    node._rid="n"+(++nextId);
    (node.children||[]).forEach(assignIds);
  }
  assignIds(work);

  let parent=null;
  let target=work;
  for(let i=0;i<path.length;i++){
    const idx=path[i];
    if(!target.children||idx>=target.children.length) return _deepCopyDict(root);
    parent=target;
    target=target.children[idx];
  }
  if(!parent) return _deepCopyDict(root);

  const nodes=new Map();
  const adj=new Map();
  function indexTree(node){
    const data=Object.assign({}, node);
    delete data.children;
    delete data._rid;
    nodes.set(node._rid, data);
    if(!adj.has(node._rid)) adj.set(node._rid, []);
    for(const ch of (node.children||[])){
      indexTree(ch);
      const len=Number(ch.dist||0);
      adj.get(node._rid).push({id:ch._rid, len});
      adj.get(ch._rid).push({id:node._rid, len});
    }
  }
  indexTree(work);

  function buildFrom(nodeId, fromId, distFromParent){
    const data=Object.assign({}, nodes.get(nodeId) || {});
    if(distFromParent!=null) data.dist=distFromParent;
    const children=(adj.get(nodeId)||[])
      .filter(edge=>edge.id!==fromId)
      .map(edge=>buildFrom(edge.id, nodeId, edge.len));
    if(children.length){
      data.children=children;
      delete data.leaf;
    } else {
      delete data.children;
    }
    return data;
  }

  const edgeLen=Number(target.dist||0);
  const aboveLen=edgeLen*splitFrac;
  const targetLen=edgeLen-aboveLen;
  const newRoot={
    name:"",
    leaf:false,
    children:[
      buildFrom(target._rid, parent._rid, targetLen),
      buildFrom(parent._rid, target._rid, aboveLen),
    ],
  };
  return _contractUnaryNodes(newRoot, true);
}

function rerootAtNode(d, splitFrac){
  if(!rootNode||!d.parent) return;   // can't reroot at existing root
  // Expand all collapsed nodes so the full topology is available
  rootNode.each(n=>{if(n._children){n.children=n._children;n._children=null;}});
  // Save original tree only on the first reroot (after expanding)
  if(!_isRerooted) _origTreeDictForReroot=_h2d(rootNode);
  // Build path (array of child indices) from root down to d
  const path=[];
  let cur=d;
  while(cur.parent){
    const idx=(cur.parent.children||[]).indexOf(cur);
    if(idx<0){ console.warn("reroot: node not found in parent.children"); return; }
    path.unshift(idx);
    cur=cur.parent;
  }
  const currentDict=_h2d(rootNode);
  const newRootDict=_rerootDictAt(currentDict,path, splitFrac);
  _isRerooted=true;
  // Clear POSSVM state — it's invalidated when topology changes (re-run after rerooting)
  if(_pvmActive){
    _pvmActive=false; _pvmIngroupSps=null; _pvmOgs=null;
    cladeHighlights.forEach((rec,uid)=>{ if(rec._fromOgHl) cladeHighlights.delete(uid); });
    if(_ogHlActive){ _ogHlActive=false; document.getElementById("btn-highlight-ogs").classList.remove("active-btn"); }
    document.getElementById("possvm-reset-btn").style.display="none";
    document.getElementById("possvm-result").textContent="";
  }
  drawGeneTree(newRootDict,{preserveRerootState:true, preserveFocusState:_isSubtreeFocused});
}
// ── tree controls ──
function expandAll(){
  if(!rootNode)return;
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;d._isOgCol=false;}});
  renderTree(true);
}
let _ogCollapseActive=false;
function toggleCollapseToOGs(){
  if(_ogCollapseActive){ _expandAllOgCol(); return; }
  _collapseToOGs();
}
function _expandAllOgCol(){
  if(!rootNode)return;
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;} d._isOgCol=false; delete d.data._og_label;});
  _ogCollapseActive=false;
  document.getElementById("btn-collapse-ogs").classList.remove("active-btn");
  renderTree(true); setTimeout(fitTree,260);
}
function _collapseToOGs(){
  if(!rootNode)return;
  // pass 1: expand all, clear flags
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;} d._isOgCol=false;});
  // pass 2a: collapse named OG internal nodes (POSSVM trees with annotated internals)
  let found=false;
  rootNode.each(d=>{
    if(!d.data.leaf&&isOGNode(d)&&d.children){d._children=d.children;d.children=null;d._isOgCol=true;found=true;}
  });
  // pass 2b: derive OGs from leaf og field or ogGene2Name map for any OGs not
  // already collapsed by pass 2a.  Leaves inside pass-2a triangles are invisible
  // to rootNode.leaves(), so there is no double-processing.
  {
    const ogGroups={};
    rootNode.leaves().forEach(l=>{
      const gid=l.data.gene_id||l.data.name||"";
      const og=leafOgName(l);
      if(og)(ogGroups[og]=ogGroups[og]||[]).push(l);
    });
    for(const [og,leaves] of Object.entries(ogGroups)){
      const mrca=findExactOgRoot(leaves);
      if(mrca&&!mrca.data.leaf&&mrca.children){
        mrca.data._og_label=og; mrca._isOgCol=true;
        mrca._children=mrca.children; mrca.children=null;
      }
    }
  }
  _ogCollapseActive=true;
  document.getElementById("btn-collapse-ogs").classList.add("active-btn");
  renderTree(true); setTimeout(fitTree,260);
}
// keep the bare name available for internal callers (e.g. re-collapse on tree load)
function collapseToOGs(){ _collapseToOGs(); }

let _ogHlActive=false;
function toggleHighlightOGs(){
  if(!rootNode)return;
  if(_ogHlActive){
    // remove OG highlights
    cladeHighlights.forEach((rec,uid)=>{ if(rec._fromOgHl) cladeHighlights.delete(uid); });
    _ogHlActive=false;
    document.getElementById("btn-highlight-ogs").classList.remove("active-btn");
    renderTree(false); return;
  }
  // Ensure tree is fully expanded so we can find all OG nodes
  rootNode.each(d=>{if(d._children){d.children=d._children;d._children=null;} d._isOgCol=false;});
  // Collect OG root nodes (same logic as _collapseToOGs pass 2a/2b)
  const ogNodes=[];   // [{node, label}]
  const ogs=activeOgs();
  const foundOgNames=new Set();
  rootNode.each(d=>{
    if(!d.data.leaf&&isOGNode(d)){ogNodes.push({node:d,label:d.data.name}); foundOgNames.add(d.data.name);}
  });
  // Also collect OGs present only in leaf annotations (not named on internal nodes)
  {
    const ogGroups={};
    rootNode.leaves().forEach(l=>{const og=leafOgName(l); if(og&&!foundOgNames.has(og))(ogGroups[og]=ogGroups[og]||[]).push(l);});
    for(const [og,leaves] of Object.entries(ogGroups)){
      const mrca=findExactOgRoot(leaves);
      if(mrca&&!mrca.data.leaf) ogNodes.push({node:mrca,label:og});
    }
  }
  if(!ogNodes.length){alert("No OG nodes found."); return;}
  // Remove any previous OG highlights, then add one per OG
  cladeHighlights.forEach((rec,uid)=>{ if(rec._fromOgHl) cladeHighlights.delete(uid); });
  ogNodes.forEach(({node,label},i)=>{
    const col=ogBaseColor(label, i);
    cladeHighlights.set(node._uid,{
      color:col,
      label,
      subtitle:ogHighlightSubtitle(node),
      _fromOgHl:true
    });
  });
  _ogCollapseActive=false;
  document.getElementById("btn-collapse-ogs").classList.remove("active-btn");
  _ogHlActive=true;
  document.getElementById("btn-highlight-ogs").classList.add("active-btn");
  renderTree(true); setTimeout(fitTree,260);
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
        const og=leafOgName(d);
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
    p.add_argument("--family_info", default="genefam.csv",
                   help="TSV with family metadata (genefam.csv or gene_families_searchinfo.csv)")
    p.add_argument("--species_tree", default=None,
                   help="Newick species tree (named internal nodes for clade colouring)")
    p.add_argument("--refnames", default=None,
                   help="POSSVM refnames TSV (gene_id<TAB>name) for OG naming and tooltip annotation")
    p.add_argument("--refsps", default=None,
                   help="Comma-separated reference species filter, matching POSSVM --refsps")
    p.add_argument("--output", required=True, help="Output HTML file path")
    return p.parse_args(argv)


def _load_species_tree_bundle(species_tree_path):
    species_order: list = []
    tree_dict: dict = {}
    newick_raw = ""
    clade_groupings: list = []
    if species_tree_path and Path(species_tree_path).exists():
        species_order, tree_dict = load_tree_data(species_tree_path)
        newick_raw = Path(species_tree_path).read_text().strip()
        clade_groupings = parse_clade_groupings(Path(species_tree_path))
        print(f"Extracted {len(clade_groupings)} clade groupings.", file=sys.stderr)
    return species_order, tree_dict, newick_raw, clade_groupings


def _load_tree_records(possvm_dir: Path, possvm_prev_dir):
    if not possvm_dir.exists():
        print(f"WARN: {possvm_dir} does not exist – no gene trees.", file=sys.stderr)

    records, all_species, gene_meta = (
        load_possvm_trees(possvm_dir, source="generax") if possvm_dir.exists() else ([], [], {})
    )
    print(f"Loaded {len(records)} gene trees, {len(all_species)} species.", file=sys.stderr)

    generax_ids = {r["id"] for r in records}
    prev_records: dict = {}
    if possvm_prev_dir:
        prev_dir = Path(possvm_prev_dir)
        if prev_dir.is_dir():
            prev_list, prev_sp, prev_gene_meta = load_possvm_trees(prev_dir, source="prev")
            prev_records = {r["id"]: r for r in prev_list}
            for gene_id, meta in prev_gene_meta.items():
                gene_meta.setdefault(gene_id, {}).update(meta)
            all_species = sorted(set(all_species) | set(prev_sp))
            records.extend(r for r in prev_list if r["id"] not in generax_ids)
            print(
                f"Loaded {len(prev_records)} prev gene trees (original IQ-TREE2).",
                file=sys.stderr,
            )

    return records, all_species, prev_records, gene_meta


def _keep_record_for_family_info(rec: dict, family_info: dict) -> bool:
    return rec["prefix"] in family_info or rec["family"] in family_info


def _filter_report_inputs(records, prev_records, family_records, hg_records, family_info):
    if not family_info:
        return records, prev_records, family_records, hg_records

    before = len(records), len(family_records), len(hg_records)
    records = [r for r in records if _keep_record_for_family_info(r, family_info)]
    prev_records = {
        rec_id: rec
        for rec_id, rec in prev_records.items()
        if _keep_record_for_family_info(rec, family_info)
    }
    family_records = [r for r in family_records if r["family"] in family_info]
    hg_records = [r for r in hg_records if r["family"] in family_info]
    print(
        f"Filtered to genefam families: "
        f"{before[0]}→{len(records)} trees, "
        f"{before[1]}→{len(family_records)} families, "
        f"{before[2]}→{len(hg_records)} HGs.",
        file=sys.stderr,
    )
    return records, prev_records, family_records, hg_records


def _build_index_records(records, prev_records, family_info):
    index_records = []
    for rec in records:
        idx = {k: v for k, v in rec.items() if k not in ("tree_dict", "ogs")}
        idx["has_prev"] = rec["id"] in prev_records
        idx["source"] = rec.get("source", "generax")
        fam_key = rec["family"] if rec["family"] in family_info else rec["prefix"]
        idx["class"] = family_info.get(fam_key, rec.get("prefix", ""))
        index_records.append(idx)
    return index_records


def _build_no_tree_genes(hg_records, tree_hg_ids, cluster_dir: Path) -> dict:
    no_tree_genes: dict = {}
    for rec in hg_records:
        if rec["id"] in tree_hg_ids:
            continue
        fasta = cluster_dir / (rec["id"] + ".fasta")
        genes_by_sp = parse_fasta_genes(fasta)
        if genes_by_sp:
            no_tree_genes[rec["id"]] = genes_by_sp
    return no_tree_genes


def _build_family_info_records(
    family_details: dict,
    family_records: list,
    hg_records: list,
    records: list,
    prev_records: dict,
):
    fam_hg_counts = defaultdict(int)
    fam_gene_counts = defaultdict(int)
    fam_species_sets = defaultdict(set)
    for rec in hg_records:
        fam_hg_counts[rec["family"]] += 1
    for rec in family_records:
        fam_gene_counts[rec["family"]] = rec.get("total", 0)
        fam_species_sets[rec["family"]] = set(rec.get("species_counts", {}).keys())

    fam_generax_counts = defaultdict(int)
    fam_tree_hg_ids = defaultdict(set)
    for rec in records:
        fam_tree_hg_ids[rec["family"]].add(rec["id"])
        if rec.get("source") == "generax":
            fam_generax_counts[rec["family"]] += 1
    for rec in prev_records.values():
        fam_tree_hg_ids[rec["family"]].add(rec["id"])

    family_info_records = []
    all_families = sorted(
        set(family_details.keys()) | set(fam_hg_counts.keys()) | set(fam_gene_counts.keys())
    )
    for fam in all_families:
        det = family_details.get(fam, {})
        family_info_records.append({
            "family": fam,
            "pfam": det.get("pfam", []),
            "category": det.get("category", ""),
            "cls": det.get("cls", ""),
            "n_hgs": fam_hg_counts.get(fam, 0),
            "total": fam_gene_counts.get(fam, 0),
            "n_species": len(fam_species_sets.get(fam, set())),
            "n_trees": len(fam_tree_hg_ids.get(fam, set())),
            "n_generax": fam_generax_counts.get(fam, 0),
        })

    have_generax = any(rec.get("source") == "generax" for rec in records)
    return family_info_records, have_generax


def _build_lazy_scripts(records, prev_records):
    lazy_parts = []
    for rec in records:
        detail = {"tree": rec["tree_dict"], "ogs": rec["ogs"]}
        prev = prev_records.get(rec["id"])
        if prev:
            detail["prev_tree"] = prev["tree_dict"]
            detail["prev_ogs"] = prev["ogs"]
        tag_id = _html.escape(rec["id"], quote=True)
        lazy_parts.append(
            f'<script type="application/json" id="treedata-{tag_id}">'
            + json.dumps(detail, separators=(",", ":"))
            + "</script>"
        )
    return "\n".join(lazy_parts)


def build_report_context(args) -> dict:
    possvm_dir = Path(args.possvm_dir)
    search_dir = Path(args.search_dir)
    cluster_dir = Path(args.cluster_dir)

    family_info = load_family_info(args.family_info)
    family_details = load_family_details(args.family_info)
    family_records = build_family_records(search_dir, family_info)
    hg_records = build_hg_records(cluster_dir, family_info)

    species_order, tree_dict, newick_raw, clade_groupings = _load_species_tree_bundle(
        args.species_tree
    )

    records, all_species, prev_records, gene_meta = _load_tree_records(possvm_dir, args.possvm_prev_dir)
    records, prev_records, family_records, hg_records = _filter_report_inputs(
        records, prev_records, family_records, hg_records, family_info
    )
    print(
        f"Loaded {len(family_records)} families, {len(hg_records)} HGs for heatmap.",
        file=sys.stderr,
    )

    index_records = _build_index_records(records, prev_records, family_info)
    tree_hg_ids = {rec["id"] for rec in records}
    no_tree_genes = _build_no_tree_genes(hg_records, tree_hg_ids, cluster_dir)
    family_info_records, have_generax = _build_family_info_records(
        family_details, family_records, hg_records, records, prev_records
    )
    domain_hits = load_domain_hits(search_dir)
    gene_lengths = load_gene_lengths(cluster_dir)
    refname_map = load_reference_names(args.refnames, args.refsps)
    for gene_id, length in gene_lengths.items():
        gene_meta.setdefault(gene_id, {})["length"] = length
    for gene_id, ref_name in refname_map.items():
        gene_meta.setdefault(gene_id, {})["is_reference_gene"] = True
        gene_meta[gene_id]["reference_gene_name"] = ref_name

    return {
        "species_order": species_order,
        "tree_dict": tree_dict,
        "family_records": family_records,
        "hg_records": hg_records,
        "index_records": index_records,
        "all_species": all_species,
        "clade_groupings": clade_groupings,
        "newick_raw": newick_raw,
        "family_info_records": family_info_records,
        "have_generax": have_generax,
        "no_tree_genes": no_tree_genes,
        "domain_hits": domain_hits,
        "gene_meta": gene_meta,
        "refname_map": refname_map,
        "records": records,
        "prev_records": prev_records,
    }


def render_report_html(context: dict) -> str:
    return (
        HTML_TEMPLATE
        .replace("%%LAZY_SCRIPTS%%", _build_lazy_scripts(context["records"], context["prev_records"]))
        .replace("%%SPECIES_ORDER%%", json.dumps(context["species_order"]))
        .replace("%%TREE_DATA%%", json.dumps(context["tree_dict"]))
        .replace("%%FAMILY_DATA%%", json.dumps(context["family_records"]))
        .replace("%%HG_DATA%%", json.dumps(context["hg_records"]))
        .replace("%%TREE_INDEX_JSON%%", json.dumps(context["index_records"]))
        .replace("%%SPECIES_JSON%%", json.dumps(context["all_species"]))
        .replace("%%CLADE_DATA_JSON%%", json.dumps(context["clade_groupings"]))
        .replace("%%NEWICK_RAW%%", json.dumps(context["newick_raw"]))
        .replace("%%FAMILY_INFO_JSON%%", json.dumps(context["family_info_records"]))
        .replace("%%HAVE_GENERAX_JSON%%", json.dumps(context["have_generax"]))
        .replace("%%NO_TREE_GENES_JSON%%", json.dumps(context["no_tree_genes"]))
        .replace("%%DOMAIN_DATA_JSON%%", json.dumps(context["domain_hits"]))
        .replace("%%GENE_META_JSON%%", json.dumps(context["gene_meta"]))
        .replace("%%REFNAME_MAP_JSON%%", json.dumps(context["refname_map"]))
    )


def main(argv=None):
    args = parse_args(argv)
    html = render_report_html(build_report_context(args))
    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
