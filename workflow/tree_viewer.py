#!/usr/bin/env python3
"""Standalone interactive gene-tree viewer.

Takes a single gene tree (Newick) and optionally a species tree, POSSVM OG CSV,
an alternative tree (e.g. pre-GeneRax), and a reference-name table.  Generates
a self-contained HTML with the same interactive features as the Gene Trees tab
in the full step-2 report (POSSVM, hOG map, collapse/expand, clade colouring,
highlights, species-tree mini panel, exports, …).

Requires report_step2.py to be in the same directory (the JS engine is shared).

Usage
-----
    python tree_viewer.py --tree my.newick --species_tree sptree.nwk \\
        --og_csv my.ortholog_groups.csv --output viewer.html
"""

import argparse
import html as _html
import importlib.util
import json
import sys
from pathlib import Path
from typing import Optional


# ── Import shared helpers from report_step2.py ───────────────────────────────

def _load_report_module():
    here = Path(__file__).parent
    path = here / "report_step2.py"
    if not path.exists():
        print(f"ERROR: report_step2.py not found at {path}", file=sys.stderr)
        sys.exit(1)
    spec = importlib.util.spec_from_file_location("_report_step2", path)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


_r2 = _load_report_module()
get_species_prefix    = _r2.get_species_prefix
gene_tree_to_dict     = _r2.gene_tree_to_dict
load_og_csv           = _r2.load_og_csv
load_reference_names  = _r2.load_reference_names
load_tree_data        = _r2.load_tree_data
parse_clade_groupings = _r2.parse_clade_groupings


# ── Gene-tree loading helper ──────────────────────────────────────────────────

def _load_gene_tree(path: str):
    """Return (tree_dict, leaves, species, stem) from a gene-tree Newick."""
    try:
        from ete3 import Tree  # type: ignore
    except ImportError:
        print("ERROR: ete3 not installed.  Run: pip install ete3", file=sys.stderr)
        sys.exit(1)
    t = Tree(str(path), format=1)
    tree_dict = gene_tree_to_dict(t)
    leaves = t.get_leaves()
    species = sorted({get_species_prefix(lf.name) for lf in leaves if lf.name})
    stem = Path(path).stem
    for suffix in (".ortholog_groups", ".treefile.ortholog_groups", ".generax.tree", ".generax"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return tree_dict, leaves, species, stem


# ── Template adaptation ───────────────────────────────────────────────────────

def _adapt_template(tmpl: str) -> str:
    """Strip heatmap / families / sidebar HTML and patch JS for standalone use."""

    # ── HTML structural changes ──────────────────────────────────────────────

    # Title
    tmpl = tmpl.replace("<title>Step 2 Report</title>", "<title>%%VIEWER_TITLE%%</title>")
    tmpl = tmpl.replace("  <h2>Step 2 Report</h2>\n", "  <h2>%%VIEWER_TITLE%%</h2>\n")

    # Remove heatmap top-bar buttons
    for frag in (
        '  <button class="btn" id="hm-back" style="display:none" onclick="hmBack()">&#8592; Back</button>\n',
        '  <span id="hm-breadcrumb" style="font-size:11px"></span>\n',
        '  <button id="hm-expand-og" style="display:none;padding:2px 8px;font-size:11px;border:1px solid #27ae60;color:#27ae60;border-radius:3px;background:#fff;cursor:pointer" title="Switch to Custom view with all OGs from all visible HGs" onclick="hmExpandToOGs()">&#43; Expand all to OGs</button>\n',
        '  <span id="tree-count" style="font-size:11px;color:#95a5a6;display:none"></span>\n',
    ):
        tmpl = tmpl.replace(frag, "")

    # Tab-strip: remove families + heatmap buttons; trees active by default
    tmpl = tmpl.replace(
        "    <button class=\"tab-btn\" data-tab=\"families\" onclick=\"switchTab('families')\">&#9783; Families</button>\n",
        "",
    )
    tmpl = tmpl.replace(
        "    <button class=\"tab-btn\" data-tab=\"heatmap\" onclick=\"switchTab('heatmap')\">&#9639; Counts</button>\n",
        "",
    )
    tmpl = tmpl.replace(
        "    <button class=\"tab-btn active\" data-tab=\"sptree\" onclick=\"switchTab('sptree')\">&#10022; Species Tree</button>\n",
        "    <button class=\"tab-btn\" data-tab=\"sptree\" onclick=\"switchTab('sptree')\">&#10022; Species Tree</button>\n",
    )
    tmpl = tmpl.replace(
        "    <button class=\"tab-btn\" data-tab=\"trees\" onclick=\"switchTab('trees')\">&#11044; Gene Trees</button>\n",
        "    <button class=\"tab-btn active\" data-tab=\"trees\" onclick=\"switchTab('trees')\">&#11044; Gene Trees</button>\n",
    )

    # Remove the entire heatmap pane
    hm_start = "  <!-- ── Heatmap pane ── -->\n"
    hm_end   = "\n  <!-- ── Gene tree pane ── -->"
    idx_s = tmpl.find(hm_start)
    idx_e = tmpl.find(hm_end)
    if idx_s != -1 and idx_e != -1:
        tmpl = tmpl[:idx_s] + tmpl[idx_e:]

    # Remove the families pane
    fam_start = "  <!-- ── Families pane ── -->\n"
    fam_end   = "\n  <!-- ── Species tree pane ── -->"
    idx_s = tmpl.find(fam_start)
    idx_e = tmpl.find(fam_end)
    if idx_s != -1 and idx_e != -1:
        tmpl = tmpl[:idx_s] + tmpl[idx_e:]

    # Remove sidebar within the gene-tree pane (keep #main only)
    sb_start = "      <div id=\"sidebar\">\n"
    sb_end   = "      <div id=\"main\">"
    idx_s = tmpl.find(sb_start)
    idx_e = tmpl.find(sb_end)
    if idx_s != -1 and idx_e != -1:
        tmpl = tmpl[:idx_s] + "      " + tmpl[idx_e:]

    # Remove "Annotations TSV" button from species-tree controls
    tmpl = tmpl.replace(
        '      <button class="ctrl-btn" id="btn-dl-anno" onclick="downloadAnnotations()">&#11015; Annotations TSV</button>\n',
        "",
    )

    # Single lazy-script tag (not plural)
    tmpl = tmpl.replace("%%LAZY_SCRIPTS%%", "%%LAZY_SCRIPT%%")

    # ── JS data constants ────────────────────────────────────────────────────

    tmpl = tmpl.replace(
        "const FAMILY_DATA   = %%FAMILY_DATA%%;\n",
        "const FAMILY_DATA   = [];\n",
    )
    tmpl = tmpl.replace(
        "const HG_DATA       = %%HG_DATA%%;\n",
        "const HG_DATA       = [];\n",
    )
    tmpl = tmpl.replace(
        "const FAMILY_INFO   = %%FAMILY_INFO_JSON%%;\n",
        "const FAMILY_INFO   = [];\n",
    )
    tmpl = tmpl.replace(
        "const HAVE_GENERAX  = %%HAVE_GENERAX_JSON%%;\n",
        "const HAVE_GENERAX  = false;\n",
    )
    tmpl = tmpl.replace(
        "const NO_TREE_GENES = %%NO_TREE_GENES_JSON%%;  // {hg_id: {species:[gene_ids]}} for HGs without trees\n",
        "const NO_TREE_GENES = {};\n",
    )

    # ── JS variable initialisations ──────────────────────────────────────────

    # spMeta: was computed from FAMILY_DATA/HG_DATA → empty object
    tmpl = tmpl.replace(
        "const spMeta = (()=>{\n"
        "  const m={};\n"
        "  FAMILY_DATA.forEach(d=>{ for(const [sp,n] of Object.entries(d.species_counts)){ if(!m[sp]) m[sp]={genes:0,families:0,hgs:0}; m[sp].genes+=n; m[sp].families++; } });\n"
        "  HG_DATA.forEach(d=>{ for(const sp of Object.keys(d.species_counts)){ if(!m[sp]) m[sp]={genes:0,families:0,hgs:0}; m[sp].hgs++; } });\n"
        "  return m;\n"
        "})();\n",
        "const spMeta = {};\n",
    )

    # dataSpecies + speciesOrder: derived from gene-tree ALL_SPECIES instead
    tmpl = tmpl.replace(
        "const dataSpecies = new Set([...FAMILY_DATA,...HG_DATA].flatMap(d=>Object.keys(d.species_counts)));\n"
        "let speciesOrder = SPECIES_ORDER.filter(s=>dataSpecies.has(s));\n"
        "for (const s of dataSpecies) { if (!speciesOrder.includes(s)) speciesOrder.push(s); }\n",
        "const dataSpecies = new Set(ALL_SPECIES);\n"
        "let speciesOrder = [...SPECIES_ORDER];\n",
    )

    # ── switchTab: strip heatmap / families branches ─────────────────────────

    tmpl = tmpl.replace(
        "function switchTab(name) {\n"
        "  document.querySelectorAll(\".tab-btn\").forEach(b => b.classList.toggle(\"active\", b.dataset.tab===name));\n"
        "  document.querySelectorAll(\".tab-pane\").forEach(p => p.classList.toggle(\"active\", p.id===\"pane-\"+name));\n"
        "  const tc  = document.getElementById(\"tree-count\");\n"
        "  const hb  = document.getElementById(\"hm-back\");\n"
        "  const cr  = document.getElementById(\"hm-breadcrumb\");\n"
        "  const pfx = document.getElementById(\"prefixSelect\").parentElement; // the <label>\n"
        "  if (name===\"trees\") {\n"
        "    tc.style.display = \"inline\"; pfx.style.display = \"none\";\n"
        "    hb.style.display = \"none\"; cr.textContent = \"\";\n"
        "    if (!currentIndex && TREE_INDEX.length) { renderSidebar(\"\"); selectTree(TREE_INDEX[0]); }\n"
        "  } else if (name===\"sptree\") {\n"
        "    tc.style.display = \"none\"; pfx.style.display = \"none\";\n"
        "    drawSpeciesTree();\n"
        "  } else if (name===\"families\") {\n"
        "    tc.style.display = \"none\"; pfx.style.display = \"none\";\n"
        "    hb.style.display = \"none\"; cr.textContent = \"\";\n"
        "    drawFamilyTable();\n"
        "  } else {\n"
        "    tc.style.display = \"none\"; pfx.style.display = \"\";\n"
        "    drawHeatmap(); // drawCladogram() is called at the end of drawHeatmap()\n"
        "  }\n"
        "}\n",
        "function switchTab(name) {\n"
        "  document.querySelectorAll(\".tab-btn\").forEach(b => b.classList.toggle(\"active\", b.dataset.tab===name));\n"
        "  document.querySelectorAll(\".tab-pane\").forEach(p => p.classList.toggle(\"active\", p.id===\"pane-\"+name));\n"
        "  if (name===\"trees\") {\n"
        "    if (!currentIndex && TREE_INDEX.length) selectTree(TREE_INDEX[0]);\n"
        "  } else if (name===\"sptree\") {\n"
        "    drawSpeciesTree();\n"
        "  }\n"
        "}\n",
    )

    # selectTree: guard renderSidebar call (hg-search element absent in standalone viewer)
    tmpl = tmpl.replace(
        "  renderSidebar(document.getElementById(\"hg-search\").value);\n",
        "  // renderSidebar: sidebar not present in standalone viewer\n",
    )

    # ── drawSpeciesTree patches ──────────────────────────────────────────────

    # Leaf text click: replace showSpAnnotPopup with color picker
    tmpl = tmpl.replace(
        "        .on(\"click\",(ev)=>{ hideTip(); showSpAnnotPopup(ev, n.name); });\n",
        "        .on(\"click\",(ev)=>{ hideTip();\n"
        "          openColorPicker(spColor(n.name),c=>{ SP_COLORS[n.name]=c; drawSpeciesTree(); if(rootNode) renderTree(false); });\n"
        "        });\n",
    )

    # Internal node dblclick: remove drawCladogram / drawHeatmap calls
    tmpl = tmpl.replace(
        "            orig.name=v.trim();\n"
        "            drawSpeciesTree(); drawCladogram();\n"
        "            if(document.getElementById(\"pane-heatmap\").classList.contains(\"active\")) drawHeatmap();\n",
        "            orig.name=v.trim();\n"
        "            drawSpeciesTree();\n",
    )

    # Internal node click: remove shift+click heatmap-split branch
    tmpl = tmpl.replace(
        "          if(ev.shiftKey){\n"
        "            ev.stopPropagation();\n"
        "            if(isSplit){ hmSplitSets.splice(splitIdx,1); hmSplitLabels.splice(splitIdx,1); }\n"
        "            else{ hmSplitSets.push(cleavesSet); hmSplitLabels.push(nodeLabel); }\n"
        "            updateHmSplitBar();\n"
        "            hideTip(); drawSpeciesTree(); drawHeatmap();\n"
        "          } else {\n"
        "            spCollapsed.add(n._id); hideTip(); drawSpeciesTree();\n"
        "          }\n",
        "          spCollapsed.add(n._id); hideTip(); drawSpeciesTree();\n",
    )

    # Tooltip hint: remove shift+click mention
    tmpl = tmpl.replace(
        "      const tip=nodeLabel+'<div style=\"font-size:9px;color:#aaa;margin-top:3px\">click to collapse &nbsp;·&nbsp; shift+click to split heatmap</div>';\n",
        "      const tip=nodeLabel+'<div style=\"font-size:9px;color:#aaa;margin-top:3px\">click to collapse</div>';\n",
    )

    # ── updateHmSplitBar: null-guard bar element ─────────────────────────────
    tmpl = tmpl.replace(
        "function updateHmSplitBar() {\n"
        "  const bar = document.getElementById(\"hm-split-bar\");\n"
        "  if (!hmSplitSets.length) { bar.style.display = \"none\"; return; }\n",
        "function updateHmSplitBar() {\n"
        "  const bar = document.getElementById(\"hm-split-bar\");\n"
        "  if (!bar || !hmSplitSets.length) { if(bar) bar.style.display = \"none\"; return; }\n",
    )

    # ── IIFEs that access absent heatmap / sidebar elements ─────────────────

    # Prefix selector (heatmap only)
    tmpl = tmpl.replace(
        "(function(){\n"
        "  const sel = document.getElementById(\"prefixSelect\");\n"
        "  const prefs = [\"all\",",
        "(function(){\n"
        "  const sel = document.getElementById(\"prefixSelect\");\n"
        "  if(!sel) return;\n"
        "  const prefs = [\"all\",",
    )
    # Sort-by-species selector (heatmap only)
    tmpl = tmpl.replace(
        "(function(){\n"
        "  const sel=document.getElementById(\"hm-col-sort-sp\");\n"
        "  ALL_SPECIES.forEach(",
        "(function(){\n"
        "  const sel=document.getElementById(\"hm-col-sort-sp\");\n"
        "  if(!sel) return;\n"
        "  ALL_SPECIES.forEach(",
    )
    # Class-filter selector (sidebar only)
    tmpl = tmpl.replace(
        "(function(){\n"
        "  const sel=document.getElementById(\"class-filter\");\n"
        "  const classes=[",
        "(function(){\n"
        "  const sel=document.getElementById(\"class-filter\");\n"
        "  if(!sel) return;\n"
        "  const classes=[",
    )

    # ── Event listeners that target absent heatmap elements ──────────────────

    for el_id in ("hg-search", "hm-col-font-slider", "hm-col-rot-slider", "hm-color-mode"):
        tmpl = tmpl.replace(
            f"document.getElementById(\"{el_id}\").addEventListener(",
            f"(document.getElementById(\"{el_id}\")||{{addEventListener:()=>{{}}}}).addEventListener(",
        )

    # ── Remove heatmap scroll-sync init block ────────────────────────────────

    tmpl = tmpl.replace(
        "// scroll sync between species-tree panel and heatmap panel\n"
        "(function(){\n"
        "  const tp=document.getElementById(\"tree-panel\");\n"
        "  const hp=document.getElementById(\"heatmap-panel\");\n"
        "  let syncing=false;\n"
        "  tp.addEventListener(\"scroll\",()=>{ if(!syncing){syncing=true;hp.scrollTop=tp.scrollTop;syncing=false;} });\n"
        "  hp.addEventListener(\"scroll\",()=>{ if(!syncing){syncing=true;tp.scrollTop=hp.scrollTop;syncing=false;} });\n"
        "})();\n"
        "\n",
        "",
    )

    # ── Replace init block ───────────────────────────────────────────────────

    tmpl = tmpl.replace(
        "document.getElementById(\"tree-count\").textContent =\n"
        "  TREE_INDEX.length+\" gene tree\"+(TREE_INDEX.length!==1?\"s\":\"\");\n"
        "\n"
        "document.getElementById(\"btn-og-labels\").classList.toggle(\"active-btn\", showOGLabels);\n"
        "\n"
        "if (hasHeatmapData || TREE_INDEX.length > 0) {\n"
        "  switchTab(\"sptree\");\n"
        "} else {\n"
        "  document.getElementById(\"pane-heatmap\").innerHTML =\n"
        "    '<div style=\"padding:40px;color:#999;text-align:center\">No data found.<br>'+\n"
        "    'Pass <code>--possvm_dir</code> and/or <code>--search_dir</code> / <code>--cluster_dir</code>.</div>';\n"
        "}\n",
        "document.getElementById(\"btn-og-labels\").classList.toggle(\"active-btn\", showOGLabels);\n"
        "if (TREE_INDEX.length > 0) {\n"
        "  selectTree(TREE_INDEX[0]);\n"
        "  setTimeout(fitTree, 260);\n"
        "}\n",
    )

    return tmpl


# ── Build data context ────────────────────────────────────────────────────────

def build_context(args) -> dict:
    # Gene tree
    try:
        tree_dict, leaves, species, stem = _load_gene_tree(args.tree)
    except Exception as exc:
        print(f"ERROR loading gene tree: {exc}", file=sys.stderr)
        sys.exit(1)

    # OG CSV (POSSVM output)
    og_members: dict = {}
    gene_meta:  dict = {}
    if args.og_csv and Path(args.og_csv).exists():
        og_members, gene_meta = load_og_csv(Path(args.og_csv))
        print(f"Loaded {len(og_members)} orthogroups from {args.og_csv}", file=sys.stderr)

    # Optional alternative tree (e.g. pre-GeneRax IQ-TREE2)
    prev_tree_dict: Optional[dict] = None
    prev_ogs:       dict = {}
    if args.prev_tree and Path(args.prev_tree).exists():
        prev_tree_dict, _, _, _ = _load_gene_tree(args.prev_tree)
        if args.prev_og_csv and Path(args.prev_og_csv).exists():
            prev_ogs, prev_gm = load_og_csv(Path(args.prev_og_csv))
            for gid, m in prev_gm.items():
                gene_meta.setdefault(gid, {}).update(m)

    # Reference names
    refname_map = load_reference_names(args.refnames, args.refsps)

    # Species tree
    species_order: list = []
    tree_data: dict = {}
    newick_raw = ""
    clade_groupings: list = []
    if args.species_tree and Path(args.species_tree).exists():
        species_order, tree_data = load_tree_data(args.species_tree)
        newick_raw = Path(args.species_tree).read_text().strip()
        clade_groupings = parse_clade_groupings(Path(args.species_tree))
        print(f"Species tree: {len(species_order)} tips, {len(clade_groupings)} clade groupings.", file=sys.stderr)

    # Use gene-tree species if no species order from species-tree
    if not species_order:
        species_order = list(species)

    # Derive all_species from gene tree (+ any in species_order)
    all_species = sorted(set(species) | set(species_order))

    # Build a single tree record
    rec = {
        "id":       stem,
        "hg":       stem,
        "family":   "",
        "prefix":   "",
        "source":   "generax" if not args.prev_tree else "generax",
        "n_leaves": len(leaves),
        "species":  species,
        "og_names": sorted(og_members.keys()),
        "n_ogs":    len(og_members),
        "has_prev": prev_tree_dict is not None,
        "class":    "",
    }

    # Lazy detail payload
    detail: dict = {"tree": tree_dict, "ogs": og_members}
    if prev_tree_dict is not None:
        detail["prev_tree"] = prev_tree_dict
        detail["prev_ogs"]  = prev_ogs

    title = args.title or stem

    return {
        "title":           title,
        "stem":            stem,
        "species_order":   species_order,
        "tree_data":       tree_data,
        "newick_raw":      newick_raw,
        "clade_groupings": clade_groupings,
        "all_species":     all_species,
        "index_records":   [rec],
        "detail":          detail,
        "gene_meta":       gene_meta,
        "refname_map":     refname_map,
    }


# ── Render HTML ───────────────────────────────────────────────────────────────

_ADAPTED_TEMPLATE: Optional[str] = None


def _get_template() -> str:
    global _ADAPTED_TEMPLATE
    if _ADAPTED_TEMPLATE is None:
        _ADAPTED_TEMPLATE = _adapt_template(_r2.HTML_TEMPLATE)
    return _ADAPTED_TEMPLATE


def render_html(ctx: dict) -> str:
    stem    = ctx["stem"]
    detail  = ctx["detail"]
    tag_id  = _html.escape(stem, quote=True)
    lazy    = (
        f'<script type="application/json" id="treedata-{tag_id}">'
        + json.dumps(detail, separators=(",", ":"))
        + "</script>"
    )

    return (
        _get_template()
        .replace("%%LAZY_SCRIPT%%",       lazy)
        .replace("%%VIEWER_TITLE%%",      _html.escape(ctx["title"]))
        .replace("%%SPECIES_ORDER%%",     json.dumps(ctx["species_order"]))
        .replace("%%TREE_DATA%%",         json.dumps(ctx["tree_data"]))
        .replace("%%TREE_INDEX_JSON%%",   json.dumps(ctx["index_records"]))
        .replace("%%SPECIES_JSON%%",      json.dumps(ctx["all_species"]))
        .replace("%%CLADE_DATA_JSON%%",   json.dumps(ctx["clade_groupings"]))
        .replace("%%NEWICK_RAW%%",        json.dumps(ctx["newick_raw"]))
        .replace("%%DOMAIN_DATA_JSON%%",  json.dumps({}))
        .replace("%%GENE_META_JSON%%",    json.dumps(ctx["gene_meta"]))
        .replace("%%REFNAME_MAP_JSON%%",  json.dumps(ctx["refname_map"]))
    )


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate a standalone interactive gene-tree viewer HTML.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tree", required=True,
                   help="Gene tree in Newick format (POSSVM-annotated or plain)")
    p.add_argument("--species_tree", default=None,
                   help="Species tree Newick (named internal nodes for clade colouring)")
    p.add_argument("--og_csv", default=None,
                   help="POSSVM ortholog_groups.csv (gene<TAB>OG<TAB>…)")
    p.add_argument("--prev_tree", default=None,
                   help="Alternative gene tree (e.g. pre-GeneRax IQ-TREE2); enables toggle")
    p.add_argument("--prev_og_csv", default=None,
                   help="OG CSV for --prev_tree")
    p.add_argument("--refnames", default=None,
                   help="POSSVM refnames TSV (gene_id<TAB>name)")
    p.add_argument("--refsps", default=None,
                   help="Comma-separated reference species (matches POSSVM --refsps)")
    p.add_argument("--title", default=None,
                   help="Report title shown in the header (default: tree file stem)")
    p.add_argument("--output", required=True,
                   help="Output HTML file")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    ctx  = build_context(args)
    html = render_html(ctx)
    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Viewer written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
