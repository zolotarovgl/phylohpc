#!/usr/bin/env python3
"""Generate a standalone alignment viewer based on the step2 Alignments tab.

This reuses the HTML/JS template from ../phylohpc/workflow/report_step2.py and
injects a single alignment plus a single gene tree. An optional species tree
can also be provided to enable the mini clade-filter panel used in the original
report.

Usage
-----
python alignment_viewer.py \
    --alignment data/example.aln.fasta \
    --gene-tree data/example.treefile \
    --species-tree data/species_tree.nwk \
    --output example.alignment_viewer.html
"""

from __future__ import annotations

import argparse
import base64
import gzip
import html as _html
import importlib.util
import json
import sys
from pathlib import Path
from typing import Optional


def _load_report_module():
    here = Path(__file__).resolve().parent
    path = here.parent / "phylohpc" / "workflow" / "report_step2.py"
    if not path.exists():
        print(f"ERROR: report_step2.py not found at {path}", file=sys.stderr)
        sys.exit(1)
    spec = importlib.util.spec_from_file_location("_report_step2", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


_r2 = _load_report_module()
get_species_prefix = _r2.get_species_prefix
gene_tree_to_dict = _r2.gene_tree_to_dict
load_tree_data = _r2.load_tree_data
parse_clade_groupings = _r2.parse_clade_groupings


def _parse_fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                token = line[1:].strip().split()[0]
                if token:
                    ids.append(token)
    if not ids:
        raise ValueError(f"no FASTA records found in {path}")
    return ids


def _load_gene_tree(path: Path, keep_ids: set[str]):
    try:
        from ete3 import Tree  # type: ignore
    except ImportError:
        print("ERROR: ete3 not installed. Run: pip install ete3", file=sys.stderr)
        sys.exit(1)

    tree = Tree(str(path), format=1)
    original_leaf_ids = [leaf.name for leaf in tree.get_leaves() if leaf.name]
    retained = [gid for gid in original_leaf_ids if gid in keep_ids]
    dropped = sorted(set(original_leaf_ids) - set(retained))
    if not retained:
        raise ValueError(
            "gene tree has no leaf IDs in common with the alignment headers"
        )
    if dropped:
        print(
            f"Pruned {len(dropped)} gene-tree leaves not present in the alignment.",
            file=sys.stderr,
        )
        tree.prune(retained, preserve_branch_length=True)

    tree_dict = gene_tree_to_dict(tree)
    species = sorted(
        {get_species_prefix(leaf.name) for leaf in tree.get_leaves() if leaf.name}
    )
    return tree_dict, species, original_leaf_ids


def _gzip_json(payload: dict) -> str:
    raw = json.dumps(payload, separators=(",", ":")).encode("utf-8")
    gz = gzip.compress(raw, compresslevel=9)
    return base64.b64encode(gz).decode("ascii")


def _gzip_bytes(raw: bytes) -> str:
    gz = gzip.compress(raw, compresslevel=9)
    return base64.b64encode(gz).decode("ascii")


def _adapt_template(tmpl: str) -> str:
    tmpl = tmpl.replace("<title>Step 2 Report</title>", "<title>%%VIEWER_TITLE%%</title>")
    tmpl = tmpl.replace("  <h2>Step 2 Report</h2>\n", "  <h2>%%VIEWER_TITLE%%</h2>\n")
    tmpl = tmpl.replace(
        "</style>\n</head>",
        "\n#tab-strip{display:none!important}\n"
        "#pane-align{display:flex!important}\n"
        "#body-wrap{overflow:hidden}\n"
        "</style>\n</head>",
    )

    # Keep the original pane implementation, but remove the alignment sidebar.
    sb_start = '    <div id="aln-sidebar">\n'
    sb_end = '    <div id="aln-main">\n'
    idx_s = tmpl.find(sb_start)
    idx_e = tmpl.find(sb_end)
    if idx_s != -1 and idx_e != -1:
        tmpl = tmpl[:idx_s] + tmpl[idx_e:]

    # Strip report-only alignment controls.
    for frag in (
        '        <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:5px">Structure Height:\n'
        '          <input type="range" id="aln-domain-chip-slider" min="0.6" max="2.2" step="0.1" value="2.0" style="width:78px;cursor:pointer;accent-color:#48a868" oninput="alnSetDomainChipScale(this.value)">\n'
        '          <span id="aln-domain-chip-val">2.0</span>x\n'
        '        </label>\n',
        '        <button class="ctrl-btn" onclick="alnToggleCollapseAllOGs()" id="btn-aln-collapse-ogs" title="Collapse all orthogroups to single rows">Collapse OGs</button>\n',
        '        <button class="ctrl-btn" onclick="alnToggleRefCol()" id="btn-aln-ref" title="Show/hide reference ortholog name column">Ref</button>\n',
        '        <button class="ctrl-btn active-btn" onclick="alnToggleSpeciesCol()" id="btn-aln-species" title="Show/hide species column">Species</button>\n',
        '        <button class="ctrl-btn" onclick="alnToggleSpeciesGroupCol()" id="btn-aln-species-group" title="Show/hide species-group column" style="display:none">Species Group</button>\n',
        '        <button class="ctrl-btn active-btn" onclick="alnToggleMrcaCol()" id="btn-aln-mrca" title="Show/hide MRCA clade column">MRCA</button>\n',
        '        <button class="ctrl-btn active-btn" onclick="alnToggleRangeCol()" id="btn-aln-range" title="Show/hide compact domain-architecture column">Domain Architecture</button>\n',
        '        <button class="ctrl-btn" onclick="alnToggleSource()" id="btn-aln-src" title="Toggle GeneRax / IQ-TREE tree source" style="display:none">IQ-TREE</button>\n',
    ):
        tmpl = tmpl.replace(frag, "")

    # Default to a minimal alignment-only metadata layout.
    tmpl = tmpl.replace("let alnShowSpeciesCol = true;\n", "let alnShowSpeciesCol = false;\n")
    tmpl = tmpl.replace("let alnShowMrcaCol  = true;\n", "let alnShowMrcaCol  = false;\n")
    tmpl = tmpl.replace(
        "let alnNameW        = 220;\n",
        "let alnNameW        = 320;\n"
        "let _alnAutoNameWidth = true;\n"
        "let _alnPercIdDialogPos = null;\n",
    )
    tmpl = tmpl.replace(
        "let alnShowRangeCol = true;\n",
        "let alnShowRangeCol = false;\n"
        "let alnShowGapCol   = false;\n"
        "let alnShowPercIdCol = false;\n"
        "let alnShowPercIdDock = false;\n",
    )
    tmpl = tmpl.replace(
        "let alnRangeColW    = 96;\n",
        "let alnRangeColW    = 96;\n"
        "let alnGapColW      = 96;\n"
        "let alnPercIdColW   = 104;\n",
    )
    tmpl = tmpl.replace(
        "let alnMetaColOrder = ['ref', 'species', 'group', 'mrca', 'range', 'id'];\n",
        "let alnMetaColOrder = ['id', 'gap', 'percid'];\n",
    )

    tmpl = tmpl.replace(
        "const HAVE_ARCHITECTURES = %%HAVE_ARCHITECTURES_JSON%%;\n",
        "const HAVE_ARCHITECTURES = %%HAVE_ARCHITECTURES_JSON%%;\n"
        "const HAVE_SPECIES_TREE = %%HAVE_SPECIES_TREE_JSON%%;\n",
    )

    tmpl = tmpl.replace(
        'placeholder="Enter to add…"',
        'placeholder="Enter to highlight…"'
    )
    tmpl = tmpl.replace(
        'title="Open a mini species tree — click a named node to filter the alignment to that clade"',
        'title="Open a mini species tree — click a named node to highlight matching alignment sequences for that clade"',
    )
    tmpl = tmpl.replace(
        'title:"Species tree — click a named node to filter the alignment to that clade",',
        'title:"Species tree — click a named node to highlight matching alignment sequences for that clade",',
    )
    tmpl = tmpl.replace(
        '        <button class="ctrl-btn" onclick="alnToggleConsensus()" id="btn-aln-cons" title="Show/hide consensus row">Consensus</button>\n',
        '        <button class="ctrl-btn" onclick="alnToggleConsensus()" id="btn-aln-cons" title="Show/hide consensus row">Consensus</button>\n'
        '        <button class="ctrl-btn" onclick="alnToggleGapCol()" id="btn-aln-gap" title="Show/hide gap-percentage column">Gap %</button>\n',
    )
    tmpl = tmpl.replace(
        '        <button class="ctrl-btn" onclick="alnToggleGapCol()" id="btn-aln-gap" title="Show/hide gap-percentage column">Gap %</button>\n',
        '        <button class="ctrl-btn" onclick="alnToggleGapCol()" id="btn-aln-gap" title="Show/hide gap-percentage column">Gap %</button>\n'
        '        <button class="ctrl-btn" onclick="alnTogglePercIdCol()" id="btn-aln-percid-col" title="Show/hide percentile of mean non-self percent identity">PercID</button>\n'
        '        <button class="ctrl-btn" onclick="alnToggleDSNodes()" id="btn-aln-ds" title="Show duplication and speciation nodes using the POSSVM species-overlap algorithm">D/S nodes</button>\n',
    )
    tmpl = tmpl.replace(
        '        <button class="ctrl-btn" onclick="alnToggleDSNodes()" id="btn-aln-ds" title="Show duplication and speciation nodes using the POSSVM species-overlap algorithm">D/S nodes</button>\n',
        '        <button class="ctrl-btn" onclick="alnToggleDSNodes()" id="btn-aln-ds" title="Show duplication and speciation nodes using the POSSVM species-overlap algorithm">D/S nodes</button>\n'
        '        <button class="ctrl-btn" onclick="alnOpenPercIdPopup()" id="btn-aln-percid" title="Open pairwise percent-identity heatmap popup">&#9635; PercID Heatmap</button>\n',
    )
    tmpl = tmpl.replace(
        "let alnSpFilters   = [];   // committed species filter tags\n"
        "let alnUseGeneRax  = false;\n",
        "let alnSpFilters   = [];   // committed species filter tags\n"
        "let alnShowDSNodes = false;\n"
        "const ALN_POSSVM_SOS = 0.0;\n"
        "let _alnTreeEventNodes = [];\n"
        "let _alnSelectedTreeNodeKey = null;\n"
        "let _alnPercIdState = null;\n"
        "let _alnPercIdDockState = null;\n"
        "let _alnClickedSeqGid = null;\n"
        "let alnUseGeneRax  = false;\n",
    )
    tmpl = tmpl.replace(
        "function parseFasta(text) {\n"
        "  const seqs = [];\n"
        "  let name = null, buf = [];\n"
        "  for (const line of text.split(/\\r?\\n/)) {\n"
        "    if (line.startsWith('>')) {\n"
        "      if (name !== null) seqs.push({name, seq: buf.join('')});\n"
        "      name = line.slice(1).trim();\n"
        "      buf = [];\n"
        "    } else { buf.push(line.trim()); }\n"
        "  }\n"
        "  if (name !== null) seqs.push({name, seq: buf.join('')});\n"
        "  return seqs;\n"
        "}\n",
        "function _alnHeaderGeneId(name) {\n"
        "  return String(name || '').trim().split(/\\s/)[0] || '';\n"
        "}\n"
        "\n"
        "function _alnSpeciesFromHeader(name) {\n"
        "  const gid = _alnHeaderGeneId(name);\n"
        "  const idx = gid.indexOf('_');\n"
        "  return idx > 0 ? gid.slice(0, idx) : gid;\n"
        "}\n"
        "\n"
        "function parseFasta(text) {\n"
        "  const seqs = [];\n"
        "  let name = null, buf = [];\n"
        "  for (const line of text.split(/\\r?\\n/)) {\n"
        "    if (line.startsWith('>')) {\n"
        "      if (name !== null) {\n"
        "        const gid = _alnHeaderGeneId(name);\n"
        "        seqs.push({name, seq: buf.join(''), gid, species: _alnSpeciesFromHeader(name)});\n"
        "      }\n"
        "      name = line.slice(1).trim();\n"
        "      buf = [];\n"
        "    } else { buf.push(line.trim()); }\n"
        "  }\n"
        "  if (name !== null) {\n"
        "    const gid = _alnHeaderGeneId(name);\n"
        "    seqs.push({name, seq: buf.join(''), gid, species: _alnSpeciesFromHeader(name)});\n"
        "  }\n"
        "  return seqs;\n"
        "}\n",
    )
    tmpl = tmpl.replace(
        "function _alnActiveSpeciesSet() {\n"
        "  const set = new Set();\n"
        "  (alnSeqs || []).forEach(seq => {\n"
        "    const gid = (seq.name || '').split(/\\s/)[0];\n"
        "    const sp = getSpeciesPfx(gid);\n"
        "    if (sp) set.add(sp);\n"
        "  });\n"
        "  return set.size ? set : new Set(ALL_SPECIES);\n"
        "}\n",
        "function _alnActiveSpeciesSet() {\n"
        "  const set = new Set();\n"
        "  (alnSeqs || []).forEach(seq => {\n"
        "    const sp = seq.species || _alnSpeciesFromHeader(seq.name || '');\n"
        "    if (sp) set.add(sp);\n"
        "  });\n"
        "  return set.size ? set : new Set(ALL_SPECIES);\n"
        "}\n",
    )

    tmpl = tmpl.replace(
        "function _alnFilteredSeqs() {\n"
        "  let seqs = alnSeqs;\n"
        "  const spSet = _alnResolveSpFilter();\n"
        "  if (spSet) {\n"
        "    seqs = seqs.filter(s => {\n"
        "      const gid = s.name.split(/\\s/)[0];\n"
        "      const sp = gid.split('_')[0];\n"
        "      return spSet.has(sp);\n"
        "    });\n"
        "  }\n"
        "  if (alnHiddenOGs.size) {\n"
        "    seqs = seqs.filter(s => {\n"
        "      const og = alnOgMap[s.name.split(/\\s/)[0]] || '';\n"
        "      return !alnHiddenOGs.has(og);\n"
        "    });\n"
        "  }\n"
        "  return seqs;\n"
        "}\n",
        "function _alnFilteredSeqs() {\n"
        "  let seqs = alnSeqs;\n"
        "  if (alnHiddenOGs.size) {\n"
        "    seqs = seqs.filter(s => {\n"
        "      const og = alnOgMap[s.name.split(/\\s/)[0]] || '';\n"
        "      return !alnHiddenOGs.has(og);\n"
        "    });\n"
        "  }\n"
        "  return seqs;\n"
        "}\n",
    )
    tmpl = tmpl.replace(
        "function _alnColumnDefs() {\n"
        "  return {\n",
        "function _alnGapFrac(seq) {\n"
        "  const text = seq || '';\n"
        "  if (!text.length) return 0;\n"
        "  let gaps = 0;\n"
        "  for (const aa of text) if (aa === '-' || aa === '.') gaps += 1;\n"
        "  return gaps / text.length;\n"
        "}\n"
        "\n"
        "function _alnDrawGapCell(ctx, x0, y0, w, h, frac, rowBg, isCollapsed=false) {\n"
        "  if (!Number.isFinite(frac) || w <= 10) return;\n"
        "  ctx.fillStyle = rowBg || '#edf1f4';\n"
        "  ctx.fillRect(x0, y0, w, h);\n"
        "  const padX = 8;\n"
        "  const barW = Math.max(8, w - padX * 2);\n"
        "  const barH = Math.max(5, Math.min(9, Math.floor(h * 0.42)));\n"
        "  const barX = x0 + padX;\n"
        "  const barY = Math.round(y0 + (h - barH) / 2) + 0.5;\n"
        "  ctx.fillStyle = '#d7dde5';\n"
        "  ctx.fillRect(barX, barY, barW, barH);\n"
        "  const fillW = Math.max(0, Math.min(barW, barW * frac));\n"
        "  if (fillW > 0) {\n"
        "    const colorFrac = 1 - Math.max(0, Math.min(1, frac));\n"
        "    ctx.fillStyle = isCollapsed ? 'rgba(199,151,0,0.78)' : _alnSpectralColor(colorFrac);\n"
        "    ctx.fillRect(barX, barY, fillW, barH);\n"
        "  }\n"
        "  ctx.strokeStyle = '#9aa5b1';\n"
        "  ctx.strokeRect(barX + 0.5, barY + 0.5, Math.max(0, barW - 1), Math.max(0, barH - 1));\n"
        "}\n"
        "\n"
        "function _alnDrawPercIdCell(ctx, x0, y0, w, h, pct, rowBg, isCollapsed=false) {\n"
        "  if (!Number.isFinite(pct) || w <= 10) return;\n"
        "  const frac = Math.max(0, Math.min(1, pct));\n"
        "  ctx.fillStyle = rowBg || '#edf1f4';\n"
        "  ctx.fillRect(x0, y0, w, h);\n"
        "  const padX = 8;\n"
        "  const barW = Math.max(8, w - padX * 2);\n"
        "  const barH = Math.max(5, Math.min(9, Math.floor(h * 0.42)));\n"
        "  const barX = x0 + padX;\n"
        "  const barY = Math.round(y0 + (h - barH) / 2) + 0.5;\n"
        "  ctx.fillStyle = '#d7dde5';\n"
        "  ctx.fillRect(barX, barY, barW, barH);\n"
        "  const fillW = Math.max(0, Math.min(barW, barW * frac));\n"
        "  if (fillW > 0) {\n"
        "    ctx.fillStyle = isCollapsed ? '#5787a8' : _alnSpectralColor(frac);\n"
        "    ctx.fillRect(barX, barY, fillW, barH);\n"
        "  }\n"
        "  ctx.strokeStyle = '#9aa5b1';\n"
        "  ctx.strokeRect(barX + 0.5, barY + 0.5, Math.max(0, barW - 1), Math.max(0, barH - 1));\n"
        "}\n"
        "\n"
        "function _alnEnsureNameColWidth(seqs, cons) {\n"
        "  if (!_alnAutoNameWidth) return;\n"
        "  const probe = document.createElement('canvas').getContext('2d');\n"
        "  if (!probe) return;\n"
        "  probe.font = '12px monospace';\n"
        "  let maxW = cons ? probe.measureText('Consensus').width : 0;\n"
        "  (seqs || []).forEach(seq => {\n"
        "    const label = String(seq && (seq.name || seq.gid) || '');\n"
        "    maxW = Math.max(maxW, probe.measureText(label).width);\n"
        "  });\n"
        "  alnNameW = Math.max(220, Math.min(900, Math.ceil(maxW + 28)));\n"
        "}\n"
        "\n"
        "function _alnColumnDefs() {\n"
        "  return {\n",
    )
    tmpl = tmpl.replace(
        "    range: {\n"
        "      key: 'range',\n"
        "      label: 'Domain Architecture',\n"
        "      width: alnShowRangeCol ? alnRangeColW : 0,\n"
        "      minWidth: 70,\n"
        "      resizerId: 'aln-range-resize-bar',\n"
        "    },\n"
        "    id: {\n",
        "    range: {\n"
        "      key: 'range',\n"
        "      label: 'Domain Architecture',\n"
        "      width: alnShowRangeCol ? alnRangeColW : 0,\n"
        "      minWidth: 70,\n"
        "      resizerId: 'aln-range-resize-bar',\n"
        "    },\n"
        "    gap: {\n"
        "      key: 'gap',\n"
        "      label: 'Gap %',\n"
        "      width: alnShowGapCol ? alnGapColW : 0,\n"
        "      minWidth: 64,\n"
        "    },\n"
        "    percid: {\n"
        "      key: 'percid',\n"
        "      label: 'PercID',\n"
        "      width: alnShowPercIdCol ? alnPercIdColW : 0,\n"
        "      minWidth: 72,\n"
        "    },\n"
        "    id: {\n",
    )
    tmpl = tmpl.replace(
        "    case 'range': return alnRangeColW;\n"
        "    case 'id': return alnNameW;\n",
        "    case 'range': return alnRangeColW;\n"
        "    case 'gap': return alnGapColW;\n"
        "    case 'percid': return alnPercIdColW;\n"
        "    case 'id': return alnNameW;\n",
    )
    tmpl = tmpl.replace(
        "    case 'range': alnRangeColW = width; break;\n"
        "    case 'id': alnNameW = width; break;\n",
        "    case 'range': alnRangeColW = width; break;\n"
        "    case 'gap': alnGapColW = width; break;\n"
        "    case 'percid': alnPercIdColW = width; break;\n"
        "    case 'id': alnNameW = width; break;\n",
    )
    tmpl = tmpl.replace(
        "function _alnColTooltip(row, key) {\n"
        "  if (!row || row.isCons || row.isCollapsed) return '';\n"
        "  if (key !== 'range') return '';\n"
        "  const meta = _alnRangeMeta(row.gid);\n"
        "  if (!meta) return '';\n"
        "  if (meta.hits && meta.hits.length && meta.start != null && meta.end != null) {\n"
        "    return `Domain architecture: ${meta.start}-${meta.end} / ${meta.length || '?'} aa`;\n"
        "  }\n"
        "  if (meta.length != null && meta.start != null && meta.end != null) {\n"
        "    return `Range: ${meta.start}-${meta.end} / ${meta.length} aa`;\n"
        "  }\n"
        "  if (meta.length != null) return `Protein length: ${meta.length} aa`;\n"
        "  if (meta.start != null && meta.end != null) return `Range: ${meta.start}-${meta.end}`;\n"
        "  return '';\n"
        "}\n",
        "function _alnColTooltip(row, key) {\n"
        "  if (!row || row.isCons) return '';\n"
        "  if (key === 'gap') {\n"
        "    if (!Number.isFinite(row.gapFrac)) return '';\n"
        "    const pct = (row.gapFrac * 100).toFixed(1);\n"
        "    return row.isCollapsed ? `Mean gap fraction: ${pct}%` : `Gap fraction: ${pct}%`;\n"
        "  }\n"
        "  if (key === 'percid') {\n"
        "    if (!Number.isFinite(row.percIdPct)) return '';\n"
        "    const pct = (row.percIdPct * 100).toFixed(1);\n"
        "    const mean = Number.isFinite(row.percIdMean) ? `${(row.percIdMean * 100).toFixed(1)}%` : 'n/a';\n"
        "    return row.isCollapsed\n"
        "      ? `Mean non-self identity percentile: ${pct}% · mean identity ${mean}`\n"
        "      : `Non-self identity percentile: ${pct}% · mean identity ${mean}`;\n"
        "  }\n"
        "  if (row.isCollapsed) return '';\n"
        "  if (key !== 'range') return '';\n"
        "  const meta = _alnRangeMeta(row.gid);\n"
        "  if (!meta) return '';\n"
        "  if (meta.hits && meta.hits.length && meta.start != null && meta.end != null) {\n"
        "    return `Domain architecture: ${meta.start}-${meta.end} / ${meta.length || '?'} aa`;\n"
        "  }\n"
        "  if (meta.length != null && meta.start != null && meta.end != null) {\n"
        "    return `Range: ${meta.start}-${meta.end} / ${meta.length} aa`;\n"
        "  }\n"
        "  if (meta.length != null) return `Protein length: ${meta.length} aa`;\n"
        "  if (meta.start != null && meta.end != null) return `Range: ${meta.start}-${meta.end}`;\n"
        "  return '';\n"
        "}\n",
    )
    tmpl = tmpl.replace(
        "function renderAlignment() {\n",
        "function _alnLeafKey(ids){\n"
        "  return [...new Set((ids || []).filter(Boolean))]\n"
        "    .sort((a,b)=>String(a).localeCompare(String(b), undefined, {numeric:true, sensitivity:\"base\"}))\n"
        "    .join(\"\\x1f\");\n"
        "}\n"
        "\n"
        "function _alnTreeFocusColor(groupIndex) {\n"
        "  return (Number(groupIndex) % 2 === 1) ? '#c96a10' : '#b42318';\n"
        "}\n"
        "\n"
        "function _alnTreeFocusFill(groupIndex, alpha=0.16) {\n"
        "  return (Number(groupIndex) % 2 === 1)\n"
        "    ? `rgba(201,106,16,${alpha})`\n"
        "    : `rgba(217,45,32,${alpha})`;\n"
        "}\n"
        "\n"
        "function _alnTreeSubsetVisible(leaves, selectedGenes) {\n"
        "  if (!selectedGenes) return true;\n"
        "  const ids = (leaves || []).filter(Boolean);\n"
        "  return !!ids.length && ids.every(gid => selectedGenes.has(gid));\n"
        "}\n"
        "\n"
        "function _alnSpeciesByGene(seqs){\n"
        "  const map=new Map();\n"
        "  (seqs||[]).forEach(seq=>{\n"
        "    const gid=seq.gid || _alnHeaderGeneId(seq.name || \"\");\n"
        "    const sp=seq.species || _alnSpeciesFromHeader(seq.name || \"\");\n"
        "    if(gid && sp) map.set(gid, sp);\n"
        "  });\n"
        "  return map;\n"
        "}\n"
        "\n"
        "function _alnPossvmEventMap(treeDict, seqs, sos){\n"
        "  if(!treeDict) return new Map();\n"
        "  const hier=hierarchyFromTreeData(treeDict);\n"
        "  if(!hier) return new Map();\n"
        "  const gene2sp=_alnSpeciesByGene(seqs);\n"
        "  const naturalCmp=(a,b)=>String(a).localeCompare(String(b), undefined, {numeric:true, sensitivity:\"base\"});\n"
        "  function overlapSpecies(node){\n"
        "    const ch=node.children||node._children||[];\n"
        "    const dup=new Set();\n"
        "    for(let i=0;i<ch.length;i++){\n"
        "      for(let j=i+1;j<ch.length;j++){\n"
        "        const si=ch[i]._pvmAllSps||new Set();\n"
        "        const sj=ch[j]._pvmAllSps||new Set();\n"
        "        if(!si.size||!sj.size) continue;\n"
        "        const shared=[];\n"
        "        si.forEach(sp=>{ if(sj.has(sp)) shared.push(sp); });\n"
        "        if(shared.length && (shared.length / Math.min(si.size, sj.size)) > (sos == null ? 0 : sos)) {\n"
        "          shared.forEach(sp=>dup.add(sp));\n"
        "        }\n"
        "      }\n"
        "    }\n"
        "    return [...dup].sort(naturalCmp);\n"
        "  }\n"
        "  function focusGroups(node, dupSpecies){\n"
        "    const ch=node.children||node._children||[];\n"
        "    if(!dupSpecies || !dupSpecies.length || !ch.length) return [];\n"
        "    const dupSet=new Set(dupSpecies);\n"
        "    const groups=[];\n"
        "    ch.forEach(child=>{\n"
        "      const genes=(child.leaves ? child.leaves() : [])\n"
        "        .map(leaf=>leaf.data.gene_id || leaf.data.name || \"\")\n"
        "        .filter(gid=>gid && dupSet.has(gene2sp.get(gid) || _alnSpeciesFromHeader(gid)));\n"
        "      if(genes.length) groups.push([...new Set(genes)].sort(naturalCmp));\n"
        "    });\n"
        "    return groups;\n"
        "  }\n"
        "  const ingroup=new Set();\n"
        "  hier.each(node=>{\n"
        "    if(node.data && node.data.leaf){\n"
        "      const gid=node.data.gene_id || node.data.name || \"\";\n"
        "      const sp=gene2sp.get(gid) || node.data.species || _alnSpeciesFromHeader(gid);\n"
        "      node.data.species=sp || \"\";\n"
        "      if(sp) ingroup.add(sp);\n"
        "    }\n"
        "  });\n"
        "  if(!ingroup.size) return new Map();\n"
        "  computePossvmAssignments(hier, ingroup, sos == null ? 0 : sos);\n"
        "  const out=new Map();\n"
        "  hier.descendants().forEach(node=>{\n"
        "    if(!node.data || node.data.leaf || !node._pvmEvent) return;\n"
        "    const genes=node.leaves().map(leaf=>leaf.data.gene_id || leaf.data.name || \"\").filter(Boolean);\n"
        "    const dupSpecies=node._pvmEvent === 'D' ? overlapSpecies(node) : [];\n"
        "    const focusGeneGroups=node._pvmEvent === 'D' ? focusGroups(node, dupSpecies) : [];\n"
        "    const focusGenes=focusGeneGroups.flat();\n"
        "    const key=_alnLeafKey(genes);\n"
        "    if(key) out.set(key, {event: node._pvmEvent, genes, focusGenes, focusGeneGroups, dupSpecies});\n"
        "  });\n"
        "  return out;\n"
        "}\n"
        "\n"
        "function _alnLayoutTreeWithNodes(treeDict, gidY, treeW, useBrLen, eventMap){\n"
        "  if (!treeDict || !Object.keys(gidY).length) return null;\n"
        "  eventMap = eventMap instanceof Map ? eventMap : new Map();\n"
        "\n"
        "  function _alnMergeTreeCoords(coords) {\n"
        "    const byY = new Map();\n"
        "    coords.forEach(cc => {\n"
        "      if (!cc) return;\n"
        "      const key = Number(cc.y);\n"
        "      const prev = byY.get(key);\n"
        "      if (!prev || cc.x > prev.x) byY.set(key, cc);\n"
        "    });\n"
        "    return [...byY.values()].sort((a, b) => a.y - b.y);\n"
        "  }\n"
        "\n"
        "  function pushEventNode(nodes, x, y, leaves){\n"
        "    const key=_alnLeafKey(leaves);\n"
        "    const rec=eventMap.get(key) || null;\n"
        "    if(rec) nodes.push({x, y, key, event: rec.event, genes: rec.genes || leaves, focusGenes: rec.focusGenes || [], focusGeneGroups: rec.focusGeneGroups || [], dupSpecies: rec.dupSpecies || []});\n"
        "  }\n"
        "\n"
        "  if (!useBrLen) {\n"
        "    function maxDepth(n, d) {\n"
        "      if (n.leaf) return d;\n"
        "      let mx = d;\n"
        "      for (const c of (n.children || [])) mx = Math.max(mx, maxDepth(c, d + 1));\n"
        "      return mx;\n"
        "    }\n"
        "    const md = maxDepth(treeDict, 0) || 1;\n"
        "    const xScale = (treeW - 4) / md;\n"
        "    const lines = [];\n"
        "    const nodes = [];\n"
        "    function lay(node, depth) {\n"
        "      if (node.leaf) {\n"
        "        const gid = node.gene_id || node.name || '';\n"
        "        const y = gidY[gid];\n"
        "        if (y === undefined) return null;\n"
        "        return {x: treeW - 2, y, leaves:[gid]};\n"
        "      }\n"
        "      const childCoords = [];\n"
        "      for (const c of (node.children || [])) {\n"
        "        const cc = lay(c, depth + 1);\n"
        "        if (cc) childCoords.push(cc);\n"
        "      }\n"
        "      if (!childCoords.length) return null;\n"
        "      const mergedCoords = _alnMergeTreeCoords(childCoords.map(cc => ({x:cc.x, y:cc.y})));\n"
        "      const leaves = childCoords.flatMap(cc => cc.leaves || []);\n"
        "      if (!mergedCoords.length) return null;\n"
        "      if (mergedCoords.length === 1) return {x: mergedCoords[0].x, y: mergedCoords[0].y, leaves};\n"
        "      const x = 2 + depth * xScale;\n"
        "      const y = (mergedCoords[0].y + mergedCoords[mergedCoords.length - 1].y) / 2;\n"
        "      lines.push({x1: x, y1: mergedCoords[0].y, x2: x, y2: mergedCoords[mergedCoords.length - 1].y, leaves});\n"
        "      for (const cc of mergedCoords) lines.push({x1: x, y1: cc.y, x2: cc.x, y2: cc.y, leaves: cc.leaves || []});\n"
        "      pushEventNode(nodes, x, y, leaves);\n"
        "      return {x, y, leaves};\n"
        "    }\n"
        "    lay(treeDict, 0);\n"
        "    return {lines, nodes};\n"
        "  }\n"
        "\n"
        "  function maxRootDist(n, d) {\n"
        "    const dd = d + (n.dist || 0);\n"
        "    if (n.leaf) return dd;\n"
        "    let mx = dd;\n"
        "    for (const c of (n.children || [])) mx = Math.max(mx, maxRootDist(c, dd));\n"
        "    return mx;\n"
        "  }\n"
        "  const maxD = maxRootDist(treeDict, 0) || 1;\n"
        "  const xScale = (treeW - 4) / maxD;\n"
        "  const lines = [];\n"
        "  const nodes = [];\n"
        "  function lay(node, rootDist) {\n"
        "    const dist = rootDist + (node.dist || 0);\n"
        "    const x = 2 + dist * xScale;\n"
        "    if (node.leaf) {\n"
        "      const gid = node.gene_id || node.name || '';\n"
        "      const y = gidY[gid];\n"
        "      if (y === undefined) return null;\n"
        "      return {x, y, leaves:[gid]};\n"
        "    }\n"
        "    const childCoords = [];\n"
        "    for (const c of (node.children || [])) {\n"
        "      const cc = lay(c, dist);\n"
        "      if (cc) childCoords.push(cc);\n"
        "    }\n"
        "    if (!childCoords.length) return null;\n"
        "    const mergedCoords = _alnMergeTreeCoords(childCoords.map(cc => ({x:cc.x, y:cc.y})));\n"
        "    const leaves = childCoords.flatMap(cc => cc.leaves || []);\n"
        "    if (!mergedCoords.length) return null;\n"
        "    if (mergedCoords.length === 1) return {x: mergedCoords[0].x, y: mergedCoords[0].y, leaves};\n"
        "    const y = (mergedCoords[0].y + mergedCoords[mergedCoords.length - 1].y) / 2;\n"
        "    lines.push({x1: x, y1: mergedCoords[0].y, x2: x, y2: mergedCoords[mergedCoords.length - 1].y, leaves});\n"
        "    for (const cc of mergedCoords) lines.push({x1: x, y1: cc.y, x2: cc.x, y2: cc.y, leaves: cc.leaves || []});\n"
        "    pushEventNode(nodes, x, y, leaves);\n"
        "    return {x, y, leaves};\n"
        "  }\n"
        "  lay(treeDict, 0);\n"
        "  return {lines, nodes};\n"
        "}\n"
        "\n"
        "function renderAlignment() {\n",
    )

    tmpl = tmpl.replace(
        "function renderAlignment() {\n"
        "  if (!alnSeqs.length) return;\n"
        "  const seqs = _alnFilteredSeqs();\n"
        "  const consObj = alnShowCons ? _alnConsensus(seqs) : null;\n",
        "function renderAlignment() {\n"
        "  if (!alnSeqs.length) return;\n"
        "  const seqs = _alnFilteredSeqs();\n"
        "  const spHlSet = _alnResolveSpFilter();\n"
        "  const eventMap = alnShowDSNodes ? _alnPossvmEventMap(alnTreeDict, seqs, ALN_POSSVM_SOS) : null;\n"
        "  const eventEntries = eventMap || new Map();\n"
        "  const selectedEvent = _alnSelectedTreeNodeKey && eventEntries.has(_alnSelectedTreeNodeKey)\n"
        "    ? eventEntries.get(_alnSelectedTreeNodeKey)\n"
        "    : null;\n"
        "  const normalizeGeneId = gid => _alnHeaderGeneId(gid || '');\n"
        "  const selectedGenes = selectedEvent\n"
        "    ? new Set((selectedEvent.genes || []).map(normalizeGeneId).filter(Boolean))\n"
        "    : null;\n"
        "  const selectedFocusGenes = selectedEvent && selectedEvent.event === 'D'\n"
        "    ? new Set((selectedEvent.focusGenes || []).map(normalizeGeneId).filter(Boolean))\n"
        "    : null;\n"
        "  const selectedFocusGroups = selectedEvent && selectedEvent.event === 'D'\n"
        "    ? (selectedEvent.focusGeneGroups || []).map(group => new Set((group || []).map(normalizeGeneId).filter(Boolean)))\n"
        "    : [];\n"
        "  if (_alnSelectedTreeNodeKey && !selectedGenes) _alnSelectedTreeNodeKey = null;\n"
        "  const consObj = alnShowCons ? _alnConsensus(seqs) : null;\n",
    )
    tmpl = tmpl.replace(
        "  const nCols = (seqs[0]?.seq.length) || 0;\n",
        "  const nCols = (seqs[0]?.seq.length) || 0;\n"
        "  _alnEnsureNameColWidth(seqs, cons);\n",
    )
    tmpl = tmpl.replace(
        "  const rows = _alnRows;\n"
        "  const totalH = rows.length * cH;\n",
        "  const rows = _alnRows;\n"
        "  const percIdRowState = _alnBuildPercIdRowState(rows, seqs);\n"
        "  rows.forEach((row, idx) => {\n"
        "    const stats = percIdRowState.rowStats[idx] || {};\n"
        "    row.percIdMean = Number.isFinite(stats.mean) ? stats.mean : null;\n"
        "    row.percIdPct = Number.isFinite(stats.percentile) ? stats.percentile : null;\n"
        "  });\n"
        "  rows.forEach(row => { row._alnClicked = !!(_alnClickedSeqGid && row.gid && _alnClickedSeqGid === row.gid); });\n"
        "  const clickedRowIndex = rows.findIndex(row => row._alnClicked);\n"
        "  const totalH = rows.length * cH;\n"
        "  const percIdDockGap = alnShowPercIdDock ? 14 : 0;\n"
        "  const percIdDockW = alnShowPercIdDock ? (rows.length * cH) : 0;\n"
        "  const seqContentW = totalW + percIdDockGap + percIdDockW;\n",
    )
    tmpl = tmpl.replace(
        "    rulerCanvas.width  = totalW * dpr;\n"
        "    rulerCanvas.height = rulerH * dpr;\n"
        "    rulerCanvas.style.width  = totalW + 'px';\n",
        "    rulerCanvas.width  = seqContentW * dpr;\n"
        "    rulerCanvas.height = rulerH * dpr;\n"
        "    rulerCanvas.style.width  = seqContentW + 'px';\n",
    )
    tmpl = tmpl.replace(
        "    rCtx.clearRect(0, 0, totalW, rulerH);\n",
        "    rCtx.clearRect(0, 0, seqContentW, rulerH);\n",
    )
    tmpl = tmpl.replace(
        "    if (conservation && consBarH > 0) {\n"
        "      const barTop = 24;\n"
        "      const barH = consBarH - 2;\n"
        "      for (let c = 0; c < nCols; c++) {\n"
        "        const v = conservation[c] || 0;\n"
        "        const h = v * barH;\n"
        "        if (h < 0.5) continue;\n"
        "        const x = c * cW;\n"
        "        const r = Math.round(180 - v * 130);\n"
        "        const g = Math.round(180 - v * 100);\n"
        "        const b = Math.round(200 + v * 55);\n"
        "        rCtx.fillStyle = `rgb(${r},${g},${b})`;\n"
        "        rCtx.fillRect(x, barTop + barH - h, cW - (cW > 3 ? 1 : 0), h);\n"
        "      }\n"
        "    }\n",
        "    if (conservation && consBarH > 0) {\n"
        "      const barTop = 24;\n"
        "      const barH = consBarH - 2;\n"
        "      for (let c = 0; c < nCols; c++) {\n"
        "        const v = conservation[c] || 0;\n"
        "        const h = v * barH;\n"
        "        if (h < 0.5) continue;\n"
        "        const x = c * cW;\n"
        "        const r = Math.round(180 - v * 130);\n"
        "        const g = Math.round(180 - v * 100);\n"
        "        const b = Math.round(200 + v * 55);\n"
        "        rCtx.fillStyle = `rgb(${r},${g},${b})`;\n"
        "        rCtx.fillRect(x, barTop + barH - h, cW - (cW > 3 ? 1 : 0), h);\n"
        "      }\n"
        "    }\n"
        "    if (alnShowPercIdDock) {\n"
        "      const dockX = totalW + percIdDockGap;\n"
        "      rCtx.fillStyle = '#f5f7fa';\n"
        "      rCtx.fillRect(dockX, 0, percIdDockW, rulerH);\n"
        "      rCtx.strokeStyle = '#c9d3dd';\n"
        "      rCtx.strokeRect(dockX + 0.5, 0.5, Math.max(0, percIdDockW - 1), Math.max(0, rulerH - 1));\n"
        "      rCtx.fillStyle = '#44515e';\n"
        "      rCtx.font = '600 11px sans-serif';\n"
        "      rCtx.textAlign = 'center';\n"
        "      rCtx.fillText('PercID Heatmap', dockX + percIdDockW / 2, 18);\n"
        "      rCtx.font = '10px sans-serif';\n"
        "      rCtx.fillStyle = '#6b7683';\n"
        "      rCtx.fillText('Current row order and row height', dockX + percIdDockW / 2, 34);\n"
        "    }\n",
    )
    tmpl = tmpl.replace(
        "    seqCanvas.width  = totalW * dpr;\n"
        "    seqCanvas.height = totalH * dpr;\n"
        "    seqCanvas.style.width  = totalW + 'px';\n",
        "    seqCanvas.width  = seqContentW * dpr;\n"
        "    seqCanvas.height = totalH * dpr;\n"
        "    seqCanvas.style.width  = seqContentW + 'px';\n",
    )
    tmpl = tmpl.replace(
        "    sCtx.clearRect(0, 0, totalW, totalH);\n",
        "    sCtx.clearRect(0, 0, seqContentW, totalH);\n",
    )
    tmpl = tmpl.replace(
        "    ogBreakYs.forEach(y => {\n"
        "      sCtx.beginPath(); sCtx.moveTo(0, y); sCtx.lineTo(totalW, y); sCtx.stroke();\n"
        "    });\n",
        "    ogBreakYs.forEach(y => {\n"
        "      sCtx.beginPath(); sCtx.moveTo(0, y); sCtx.lineTo(seqContentW, y); sCtx.stroke();\n"
        "    });\n"
        "    if (alnShowPercIdDock) {\n"
        "      _alnPercIdDockState = percIdRowState;\n"
        "      _alnDrawDockedPercIdHeatmap(sCtx, totalW + percIdDockGap, rows, percIdRowState, cH, totalH);\n"
        "    } else {\n"
        "      _alnPercIdDockState = null;\n"
        "    }\n",
    )
    tmpl = tmpl.replace(
        "  if (cons) _alnRows.push({name:'Consensus', seq: cons, isCons: true, og:'', ogColor:''});\n",
        "  if (cons) _alnRows.push({name:'Consensus', seq: cons, isCons: true, og:'', ogColor:'', spMatched:false, spMatchCount:0, gapFrac:null, percIdGenes:[], percIdMean:null, percIdPct:null});\n",
    )
    tmpl = tmpl.replace(
        "        const ogCount = seqs.filter(ss => (alnOgMap[ss.name.split(/\\s/)[0]] || '') === og).length;\n"
        "        _alnRows.push({name: og, seq: '', isCons: false, isCollapsed: true, og, ogColor, ogCount, species:'', speciesLabel:'', speciesGroup:'', mrca: ogMrca[og] || '', gid:''});\n",
        "        const ogSeqs = seqs.filter(ss => (alnOgMap[ss.name.split(/\\s/)[0]] || '') === og);\n"
        "        const ogCount = ogSeqs.length;\n"
        "        const ogGapFrac = ogCount ? ogSeqs.reduce((sum, ss) => sum + _alnGapFrac(ss.seq), 0) / ogCount : null;\n"
        "        const ogMatchCount = spHlSet\n"
        "          ? ogSeqs.filter(ss => {\n"
        "              const ssGid = ss.name.split(/\\s/)[0];\n"
        "              return spHlSet.has(ss.species || _alnSpeciesFromHeader(ss.name || ssGid));\n"
        "            }).length\n"
        "          : 0;\n"
        "        const ogTreeSelCount = selectedGenes\n"
        "          ? ogSeqs.filter(ss => {\n"
        "              const ssGid = normalizeGeneId(ss.gid || ss.name || '');\n"
        "              return ssGid && selectedGenes.has(ssGid);\n"
        "            }).length\n"
        "          : 0;\n"
        "        const ogTreeFocusCount = selectedFocusGenes\n"
        "          ? ogSeqs.filter(ss => {\n"
        "              const ssGid = normalizeGeneId(ss.gid || ss.name || '');\n"
        "              return ssGid && selectedFocusGenes.has(ssGid);\n"
        "            }).length\n"
        "          : 0;\n"
        "        let ogTreeFocusGroup = -1;\n"
        "        if (selectedFocusGroups.length) {\n"
        "          for (let gi = 0; gi < selectedFocusGroups.length; gi++) {\n"
        "            if (ogSeqs.some(ss => {\n"
        "              const ssGid = normalizeGeneId(ss.gid || ss.name || '');\n"
        "              return ssGid && selectedFocusGroups[gi].has(ssGid);\n"
        "            })) { ogTreeFocusGroup = gi; break; }\n"
        "          }\n"
        "        }\n"
        "        _alnRows.push({name: og, seq: '', isCons: false, isCollapsed: true, og, ogColor, ogCount, species:'', speciesLabel:'', speciesGroup:'', mrca: ogMrca[og] || '', gid:'', geneKey:null, spMatched: !!(spHlSet && ogMatchCount), spMatchCount: ogMatchCount, gapFrac: ogGapFrac, percIdGenes: ogSeqs.map(ss => ss.gid || _alnHeaderGeneId(ss.name || '')).filter(Boolean), percIdMean:null, percIdPct:null, treeSelected: !!ogTreeSelCount, treeSelCount: ogTreeSelCount, treeFocus: !!ogTreeFocusCount, treeFocusCount: ogTreeFocusCount, treeFocusGroup: ogTreeFocusGroup});\n",
    )
    tmpl = tmpl.replace(
        "    const sp = getSpeciesPfx(gid);\n",
        "    const sp = s.species || _alnSpeciesFromHeader(s.name || gid);\n",
    )
    tmpl = tmpl.replace(
        "    _alnRows.push({...s, gid, og, ogColor, species: sp || '', speciesLabel: (SPECIES_INFO[sp] || sp || ''), speciesGroup: (SPECIES_GROUPS[sp] || ''), mrca: og ? (ogMrca[og] || '') : ''});\n",
        "    const geneKey = normalizeGeneId(gid);\n"
        "    let treeFocusGroup = -1;\n"
        "    if (geneKey && selectedFocusGroups.length) {\n"
        "      for (let gi = 0; gi < selectedFocusGroups.length; gi++) {\n"
        "        if (selectedFocusGroups[gi].has(geneKey)) { treeFocusGroup = gi; break; }\n"
        "      }\n"
        "    }\n"
        "    const isSelectedGene = geneKey && selectedGenes && selectedGenes.has(geneKey);\n"
        "    const isFocusGene = geneKey && selectedFocusGenes && selectedFocusGenes.has(geneKey);\n"
        "    _alnRows.push({...s, gid, og, ogColor, species: sp || '', speciesLabel: (SPECIES_INFO[sp] || sp || ''), speciesGroup: (SPECIES_GROUPS[sp] || ''), mrca: og ? (ogMrca[og] || '') : '', spMatched: !!(spHlSet && sp && spHlSet.has(sp)), spMatchCount: 0, gapFrac: _alnGapFrac(s.seq), percIdGenes:[gid], percIdMean:null, percIdPct:null, geneKey, treeSelected: !!isSelectedGene, treeSelCount: 0, treeFocus: !!isFocusGene, treeFocusCount: 0, treeFocusGroup});\n",
    )
    tmpl = tmpl.replace(
        "      nCtx.fillStyle = row.ogColor || '#bbb';\n"
        "      nCtx.globalAlpha = 0.45;\n",
        "      nCtx.fillStyle = row.ogColor || '#bbb';\n"
        "      nCtx.globalAlpha = selectedGenes\n"
        "        ? (row.treeSelected ? (row.treeFocus ? 0.58 : 0.45) : 0.08)\n"
        "        : (spHlSet ? (row.spMatched ? 0.62 : 0.2) : 0.45);\n",
    )
    tmpl = tmpl.replace(
        "    const lines = _alnLayoutTree(alnTreeDict, gidY, treeW, alnUseBrLen);\n"
        "    if (lines) {\n"
        "      nCtx.strokeStyle = '#333';\n"
        "      nCtx.lineWidth = 1;\n"
        "      for (const l of lines) {\n"
        "        nCtx.beginPath();\n"
        "        nCtx.moveTo(l.x1, l.y1);\n"
        "        nCtx.lineTo(l.x2, l.y2);\n"
        "        nCtx.stroke();\n"
        "      }\n"
        "    }\n",
        "    const layout = _alnLayoutTreeWithNodes(alnTreeDict, gidY, treeW, alnUseBrLen, eventMap);\n"
        "    _alnTreeEventNodes = layout && layout.nodes ? layout.nodes : [];\n"
        "    if (layout) {\n"
        "      for (const l of (layout.lines || [])) {\n"
        "        const inSelected = _alnTreeSubsetVisible(l.leaves, selectedGenes);\n"
        "        nCtx.globalAlpha = inSelected ? 1 : 0.14;\n"
        "        nCtx.strokeStyle = inSelected ? '#333' : '#9aa5b1';\n"
        "        nCtx.lineWidth = inSelected ? 1 : 0.8;\n"
        "        nCtx.beginPath();\n"
        "        nCtx.moveTo(l.x1, l.y1);\n"
        "        nCtx.lineTo(l.x2, l.y2);\n"
        "        nCtx.stroke();\n"
        "      }\n"
        "      nCtx.globalAlpha = 1;\n"
        "      if (alnShowDSNodes && layout.nodes && layout.nodes.length) {\n"
        "        const dsFont = Math.max(7, Math.min(10, cH - 4));\n"
        "        nCtx.textAlign = 'center';\n"
        "        nCtx.textBaseline = 'middle';\n"
        "        nCtx.font = `700 ${dsFont}px sans-serif`;\n"
        "        for (const nd of layout.nodes) {\n"
        "          const inSelected = _alnTreeSubsetVisible(nd.genes, selectedGenes);\n"
        "          const isSelectedNode = _alnSelectedTreeNodeKey && nd.key === _alnSelectedTreeNodeKey;\n"
        "          const fill = nd.event === 'D' ? '#e67e22' : '#27ae60';\n"
        "          const radius = nd.event === 'D' ? 5.5 : 5.0;\n"
        "          nCtx.globalAlpha = inSelected ? 1 : 0.18;\n"
        "          nCtx.fillStyle = fill;\n"
        "          nCtx.beginPath();\n"
        "          nCtx.arc(nd.x, nd.y, radius, 0, Math.PI * 2);\n"
        "          nCtx.fill();\n"
        "          nCtx.strokeStyle = isSelectedNode ? '#111' : '#fff';\n"
        "          nCtx.lineWidth = isSelectedNode ? 1.6 : 1.2;\n"
        "          nCtx.stroke();\n"
        "          nCtx.fillStyle = '#fff';\n"
        "          nCtx.fillText(nd.event, nd.x, nd.y + 0.2);\n"
        "        }\n"
        "        nCtx.globalAlpha = 1;\n"
        "      }\n"
        "    }\n",
    )
    tmpl = tmpl.replace(
        "          nCtx.fillText(_alnFitText(nCtx, `\\u25B6 ${row.name} (${row.ogCount})`, col.width - 8), col.x1 - 4, y0 + cH / 2);\n",
        "          const collapsedLabel = spHlSet && row.spMatched\n"
        "            ? `\\u25B6 ${row.name} (${row.spMatchCount}/${row.ogCount})`\n"
        "            : `\\u25B6 ${row.name} (${row.ogCount})`;\n"
        "          const focusColor = _alnTreeFocusColor(row.treeFocusGroup);\n"
        "          const nameSelected = row.treeSelected || (row.geneKey && selectedGenes && selectedGenes.has(row.geneKey));\n"
        "          nCtx.fillStyle = row.treeFocus ? focusColor : (nameSelected ? '#b24500' : (selectedGenes && !row.treeSelected ? '#8b95a1' : '#222'));\n"
        "          nCtx.fillText(_alnFitText(nCtx, collapsedLabel, col.width - 8), col.x1 - 4, y0 + cH / 2);\n",
    )
    tmpl = tmpl.replace(
        "        } else if (col.key === 'group' && row.speciesGroup) {\n"
        "          nCtx.fillStyle = groupColors[row.speciesGroup] || '#555';\n"
        "          nCtx.textAlign = 'right';\n"
        "          nCtx.fillText(_alnFitText(nCtx, row.speciesGroup, col.width - 8), col.x1 - 4, y0 + cH / 2);\n"
        "        } else if (col.key === 'id') {\n",
        "        } else if (col.key === 'group' && row.speciesGroup) {\n"
        "          nCtx.fillStyle = groupColors[row.speciesGroup] || '#555';\n"
        "          nCtx.textAlign = 'right';\n"
        "          nCtx.fillText(_alnFitText(nCtx, row.speciesGroup, col.width - 8), col.x1 - 4, y0 + cH / 2);\n"
        "        } else if (col.key === 'gap') {\n"
        "          _alnDrawGapCell(nCtx, col.x0, y0, col.width, cH, row.gapFrac, spHlSet && row.spMatched ? '#fff3bf' : '#f2f3f5', true);\n"
        "        } else if (col.key === 'percid') {\n"
        "          _alnDrawPercIdCell(nCtx, col.x0, y0, col.width, cH, row.percIdPct, spHlSet && row.spMatched ? '#fff3bf' : '#f2f3f5', true);\n"
        "        } else if (col.key === 'id') {\n",
    )
    tmpl = tmpl.replace(
        "      if (r % 2 === (cons ? 1 : 0)) {\n"
        "        nCtx.fillStyle = '#f7f7f7';\n"
        "        nCtx.fillRect(metaX0, y0, totalNameW, cH);\n"
        "      }\n",
        "      const striped = r % 2 === (cons ? 1 : 0);\n"
        "      const rowGeneKey = row.geneKey;\n"
        "      const rowIsSelected = rowGeneKey && selectedGenes && selectedGenes.has(rowGeneKey);\n"
        "      const highlightRow = row.treeSelected || rowIsSelected;\n"
        "      const rowBg = highlightRow && !row.treeFocus\n"
        "        ? 'rgba(255,215,64,0.18)'\n"
        "        : row.treeFocus\n"
        "          ? _alnTreeFocusFill(row.treeFocusGroup, 0.12)\n"
        "          : (selectedGenes && !row.treeSelected)\n"
        "            ? (striped ? '#f4f6f8' : '#fafbfc')\n"
        "            : (spHlSet\n"
        "              ? (row.spMatched ? '#fff3bf' : (striped ? '#f2f3f5' : '#fafbfc'))\n"
        "              : (striped ? '#f7f7f7' : ''));\n"
        "      if (rowBg) {\n"
        "        nCtx.fillStyle = rowBg;\n"
        "        nCtx.fillRect(metaX0, y0, totalNameW, cH);\n"
        "      }\n"
        "      if (row.treeSelected) {\n"
        "        const prevTreeSel = r > 0 && !!(rows[r - 1] && rows[r - 1].treeSelected);\n"
        "        const nextTreeSel = r + 1 < rows.length && !!(rows[r + 1] && rows[r + 1].treeSelected);\n"
        "        nCtx.strokeStyle = '#111';\n"
        "        nCtx.lineWidth = 0.9;\n"
        "        if (!prevTreeSel) {\n"
        "          nCtx.beginPath();\n"
        "          nCtx.moveTo(metaX0, y0 + 0.5);\n"
        "          nCtx.lineTo(metaX0 + totalNameW, y0 + 0.5);\n"
        "          nCtx.stroke();\n"
        "        }\n"
        "        if (!nextTreeSel) {\n"
        "          nCtx.beginPath();\n"
        "          nCtx.moveTo(metaX0, y0 + cH - 0.5);\n"
        "          nCtx.lineTo(metaX0 + totalNameW, y0 + cH - 0.5);\n"
        "          nCtx.stroke();\n"
        "        }\n"
        "      }\n",
    )
    tmpl = tmpl.replace(
        "          _alnDrawRangeCell(nCtx, col.x0, y0, col.width, cH, gid, r % 2 === (cons ? 1 : 0) ? '#eceff2' : '#edf1f4');\n",
        "          _alnDrawRangeCell(nCtx, col.x0, y0, col.width, cH, gid, spHlSet\n"
        "            ? (row.spMatched ? '#fff3bf' : (striped ? '#eef0f2' : '#f7f8fa'))\n"
        "            : (striped ? '#eceff2' : '#edf1f4'));\n",
    )
    tmpl = tmpl.replace(
        "        } else if (col.key === 'mrca' && row.mrca) {\n"
        "          nCtx.fillStyle = '#555';\n"
        "          nCtx.fillText(_alnFitText(nCtx, row.mrca, col.width - 8), col.x1 - 4, y0 + cH / 2);\n"
        "        } else if (col.key === 'range') {\n",
        "        } else if (col.key === 'mrca' && row.mrca) {\n"
        "          nCtx.fillStyle = '#555';\n"
        "          nCtx.fillText(_alnFitText(nCtx, row.mrca, col.width - 8), col.x1 - 4, y0 + cH / 2);\n"
        "        } else if (col.key === 'gap') {\n"
        "          _alnDrawGapCell(nCtx, col.x0, y0, col.width, cH, row.gapFrac, rowBg || (striped ? '#eceff2' : '#edf1f4'));\n"
        "        } else if (col.key === 'percid') {\n"
        "          _alnDrawPercIdCell(nCtx, col.x0, y0, col.width, cH, row.percIdPct, rowBg || (striped ? '#eceff2' : '#edf1f4'));\n"
        "        } else if (col.key === 'range') {\n",
    )
    tmpl = tmpl.replace(
        "          nCtx.fillStyle = isRefSeq ? '#cc0000' : '#333';\n"
        "          nCtx.fillText(_alnFitText(nCtx, row.name, col.width - 8), col.x1 - 4, y0 + cH / 2);\n",
        "          nCtx.fillStyle = row.treeFocus\n"
        "            ? _alnTreeFocusColor(row.treeFocusGroup)\n"
        "            : (selectedGenes && !row.treeSelected ? '#8b95a1' : (isRefSeq ? '#cc0000' : '#333'));\n"
        "          nCtx.fillText(_alnFitText(nCtx, row.name, col.width - 8), col.x1 - 4, y0 + cH / 2);\n",
    )
    tmpl = tmpl.replace(
        "        sCtx.font = `bold ${Math.max(10, cH - 2)}px sans-serif`;\n"
        "        sCtx.fillStyle = '#444';\n"
        "        sCtx.textAlign = 'center';\n"
        "        sCtx.textBaseline = 'middle';\n"
        "        sCtx.fillText(`${row.name} \\u00b7 ${row.ogCount} sequences`, totalW / 2, y0 + cH / 2);\n",
        "        sCtx.font = `bold ${Math.max(10, cH - 2)}px sans-serif`;\n"
        "        sCtx.fillStyle = '#444';\n"
        "        sCtx.textAlign = 'center';\n"
        "        sCtx.textBaseline = 'middle';\n"
        "        const collapsedSeqLabel = spHlSet && row.spMatched\n"
        "          ? `${row.name} \\u00b7 ${row.spMatchCount}/${row.ogCount} highlighted`\n"
        "          : `${row.name} \\u00b7 ${row.ogCount} sequences`;\n"
        "        sCtx.fillStyle = row.treeFocus ? _alnTreeFocusColor(row.treeFocusGroup) : (selectedGenes && !row.treeSelected ? '#8b95a1' : '#444');\n"
        "        sCtx.fillText(collapsedSeqLabel, totalW / 2, y0 + cH / 2);\n",
    )
    tmpl = tmpl.replace(
        "      for (let c = 0; c < nCols; c++) {\n"
        "        const aa = (row.seq[c] || '-').toUpperCase();\n"
        "        const x = c * cW;\n"
        "        let bg = AA_COLORS[aa];\n"
        "        if (!bg && aa !== '-' && aa !== '.') bg = '#e0e0e0';\n"
        "        if (isConsRow) {\n"
        "          sCtx.fillStyle = bg || '#ddf0ff'; sCtx.fillRect(x, y0, cW, cH);\n"
        "        } else {\n"
        "          if (r % 2 === (cons ? 1 : 0)) { sCtx.fillStyle = '#f7f7f7'; sCtx.fillRect(x, y0, cW, cH); }\n"
        "          if (bg) { sCtx.fillStyle = bg; sCtx.fillRect(x, y0, cW, cH); }\n"
        "        }\n"
        "        if (showLetters && aa !== '-' && aa !== '.') {\n"
        "          sCtx.fillStyle = isConsRow ? '#000' : '#222';\n"
        "          sCtx.fillText(aa, x + cW / 2, y0 + cH / 2);\n"
        "        }\n"
        "      }\n",
        "      const striped = r % 2 === (cons ? 1 : 0);\n"
        "      for (let c = 0; c < nCols; c++) {\n"
        "        const aa = (row.seq[c] || '-').toUpperCase();\n"
        "        const x = c * cW;\n"
        "        let bg = AA_COLORS[aa];\n"
        "        if (!bg && aa !== '-' && aa !== '.') bg = '#e0e0e0';\n"
        "        if (isConsRow) {\n"
        "          sCtx.fillStyle = bg || '#ddf0ff'; sCtx.fillRect(x, y0, cW, cH);\n"
        "        } else {\n"
        "          if (striped) { sCtx.fillStyle = '#f7f7f7'; sCtx.fillRect(x, y0, cW, cH); }\n"
        "          if (bg) { sCtx.fillStyle = bg; sCtx.fillRect(x, y0, cW, cH); }\n"
        "        }\n"
        "        if (showLetters && aa !== '-' && aa !== '.') {\n"
        "          sCtx.fillStyle = isConsRow ? '#000' : '#222';\n"
        "          sCtx.fillText(aa, x + cW / 2, y0 + cH / 2);\n"
        "        }\n"
        "      }\n"
        "      if (!isConsRow && spHlSet) {\n"
        "        if (row.spMatched) {\n"
        "          sCtx.fillStyle = 'rgba(255,215,64,0.14)';\n"
        "          sCtx.fillRect(0, y0, totalW, cH);\n"
        "          sCtx.strokeStyle = '#c99700';\n"
        "          sCtx.lineWidth = 2;\n"
        "          sCtx.strokeRect(1, y0 + 1, Math.max(0, totalW - 2), Math.max(0, cH - 2));\n"
        "        } else {\n"
        "          sCtx.fillStyle = 'rgba(255,255,255,0.32)';\n"
        "          sCtx.fillRect(0, y0, totalW, cH);\n"
        "        }\n"
        "      }\n"
        "      if (!isConsRow && row.treeFocus) {\n"
        "        sCtx.fillStyle = _alnTreeFocusFill(row.treeFocusGroup, 0.18);\n"
        "        sCtx.fillRect(0, y0, totalW, cH);\n"
        "      }\n"

    )
    tmpl = tmpl.replace(
        "function alnToggleRangeCol() {\n"
        "  alnShowRangeCol = !alnShowRangeCol;\n"
        "  const btn = document.getElementById('btn-aln-range');\n"
        "  if (btn) { btn.style.background = alnShowRangeCol?'#e8f0fe':'#fff'; btn.style.color = alnShowRangeCol?'#1a56c4':'#555'; btn.style.borderColor = alnShowRangeCol?'#1a56c4':'#888'; btn.style.fontWeight = alnShowRangeCol?'600':''; }\n"
        "  renderAlignment();\n"
        "}\n",
        "function alnToggleRangeCol() {\n"
        "  alnShowRangeCol = !alnShowRangeCol;\n"
        "  const btn = document.getElementById('btn-aln-range');\n"
        "  if (btn) { btn.style.background = alnShowRangeCol?'#e8f0fe':'#fff'; btn.style.color = alnShowRangeCol?'#1a56c4':'#555'; btn.style.borderColor = alnShowRangeCol?'#1a56c4':'#888'; btn.style.fontWeight = alnShowRangeCol?'600':''; }\n"
        "  renderAlignment();\n"
        "}\n"
        "function alnToggleGapCol() {\n"
        "  alnShowGapCol = !alnShowGapCol;\n"
        "  const btn = document.getElementById('btn-aln-gap');\n"
        "  if (btn) { btn.style.background = alnShowGapCol?'#e8f0fe':'#fff'; btn.style.color = alnShowGapCol?'#1a56c4':'#555'; btn.style.borderColor = alnShowGapCol?'#1a56c4':'#888'; btn.style.fontWeight = alnShowGapCol?'600':''; }\n"
        "  renderAlignment();\n"
        "}\n"
        "function alnTogglePercIdCol() {\n"
        "  alnShowPercIdCol = !alnShowPercIdCol;\n"
        "  const btn = document.getElementById('btn-aln-percid-col');\n"
        "  if (btn) { btn.style.background = alnShowPercIdCol?'#e8f0fe':'#fff'; btn.style.color = alnShowPercIdCol?'#1a56c4':'#555'; btn.style.borderColor = alnShowPercIdCol?'#1a56c4':'#888'; btn.style.fontWeight = alnShowPercIdCol?'600':''; }\n"
        "  renderAlignment();\n"
        "}\n"
        "function alnToggleDSNodes() {\n"
        "  alnShowDSNodes = !alnShowDSNodes;\n"
        "  const btn = document.getElementById('btn-aln-ds');\n"
        "  if (btn) { btn.style.background = alnShowDSNodes?'#e8f0fe':'#fff'; btn.style.color = alnShowDSNodes?'#1a56c4':'#555'; btn.style.borderColor = alnShowDSNodes?'#1a56c4':'#888'; btn.style.fontWeight = alnShowDSNodes?'600':''; }\n"
        "  renderAlignment();\n"
        "}\n"
        "function _alnSetPercIdBtn(active) {\n"
        "  const btn = document.getElementById('btn-aln-percid');\n"
        "  if (!btn) return;\n"
        "  btn.style.background = active ? '#e8f0fe' : '#fff';\n"
        "  btn.style.color = active ? '#1a56c4' : '#555';\n"
        "  btn.style.borderColor = active ? '#1a56c4' : '#888';\n"
        "  btn.style.fontWeight = active ? '600' : '';\n"
        "}\n"
        "\n"
        "function _alnIsGapChar(ch) {\n"
        "  return ch === '-' || ch === '.';\n"
        "}\n"
        "\n"
        "function _alnPercId(seqA, seqB) {\n"
        "  const a = String(seqA || '');\n"
        "  const b = String(seqB || '');\n"
        "  const n = Math.max(a.length, b.length);\n"
        "  let matches = 0;\n"
        "  let compared = 0;\n"
        "  for (let i = 0; i < n; i++) {\n"
        "    const aa = (a[i] || '-').toUpperCase();\n"
        "    const bb = (b[i] || '-').toUpperCase();\n"
        "    const gapA = _alnIsGapChar(aa);\n"
        "    const gapB = _alnIsGapChar(bb);\n"
        "    if (gapA && gapB) continue;\n"
        "    compared += 1;\n"
        "    if (aa === bb) matches += 1;\n"
        "  }\n"
        "  return compared ? (matches / compared) : null;\n"
        "}\n"
        "\n"
        "function _alnPercIdLabel(seq) {\n"
        "  return seq.gid || _alnHeaderGeneId(seq.name || '') || seq.name || '';\n"
        "}\n"
        "\n"
        "function _alnQuantile(values, q) {\n"
        "  if (!values || !values.length) return null;\n"
        "  if (values.length === 1) return values[0];\n"
        "  const pos = Math.max(0, Math.min(values.length - 1, (values.length - 1) * q));\n"
        "  const lo = Math.floor(pos);\n"
        "  const hi = Math.ceil(pos);\n"
        "  if (lo === hi) return values[lo];\n"
        "  const t = pos - lo;\n"
        "  return values[lo] * (1 - t) + values[hi] * t;\n"
        "}\n"
        "\n"
        "function _alnPercentileRank(values, value) {\n"
        "  if (!values || !values.length || !Number.isFinite(value)) return 0;\n"
        "  let lo = 0, hi = values.length;\n"
        "  while (lo < hi) {\n"
        "    const mid = (lo + hi) >> 1;\n"
        "    if (values[mid] <= value) lo = mid + 1;\n"
        "    else hi = mid;\n"
        "  }\n"
        "  return values.length === 1 ? 1 : Math.max(0, Math.min(1, (lo - 1) / (values.length - 1)));\n"
        "}\n"
        "\n"
        "function _alnMixColor(a, b, t) {\n"
        "  const tt = Math.max(0, Math.min(1, t));\n"
        "  const out = [0, 1, 2].map(i => Math.round(a[i] + (b[i] - a[i]) * tt));\n"
        "  return `rgb(${out[0]}, ${out[1]}, ${out[2]})`;\n"
        "}\n"

        "function _alnSpectralColor(t) {\n"
        "  return _alnRampColor([\n"
        "    [0.00, [180, 35, 24]],\n"
        "    [0.28, [227, 138, 77]],\n"
        "    [0.52, [244, 213, 120]],\n"
        "    [0.74, [135, 205, 184]],\n"
        "    [1.00, [44, 92, 170]],\n"
        "  ], t);\n"
        "}\n"
        "\n"
        "function _alnRampColor(stops, t) {\n"
        "  const tt = Math.max(0, Math.min(1, t));\n"
        "  for (let i = 1; i < stops.length; i++) {\n"
        "    if (tt <= stops[i][0]) {\n"
        "      const [p0, c0] = stops[i - 1];\n"
        "      const [p1, c1] = stops[i];\n"
        "      const local = p1 === p0 ? 0 : (tt - p0) / (p1 - p0);\n"
        "      return _alnMixColor(c0, c1, local);\n"
        "    }\n"
        "  }\n"
        "  return _alnMixColor(stops[stops.length - 1][1], stops[stops.length - 1][1], 0);\n"
        "}\n"
        "\n"
        "function _alnBuildPercIdState(seqs) {\n"
        "  const n = (seqs || []).length;\n"
        "  const labels = seqs.map(_alnPercIdLabel);\n"
        "  const matrix = Array.from({length: n}, () => Array(n).fill(null));\n"
        "  const values = [];\n"
        "  for (let i = 0; i < n; i++) {\n"
        "    for (let j = i; j < n; j++) {\n"
        "      const frac = _alnPercId(seqs[i]?.seq, seqs[j]?.seq);\n"
        "      matrix[i][j] = frac;\n"
        "      matrix[j][i] = frac;\n"
        "      if (Number.isFinite(frac) && i !== j) values.push(frac);\n"
        "    }\n"
        "  }\n"
        "  values.sort((a, b) => a - b);\n"
        "  const percentileScale = {\n"
        "    values,\n"
        "    q05: _alnQuantile(values, 0.05),\n"
        "    q25: _alnQuantile(values, 0.25),\n"
        "    q50: _alnQuantile(values, 0.50),\n"
        "    q75: _alnQuantile(values, 0.75),\n"
        "    q95: _alnQuantile(values, 0.95),\n"
        "  };\n"
        "  return {seqs, labels, matrix, percentileScale};\n"
        "}\n"
        "\n"
        "function _alnPercIdRowLabel(row) {\n"
        "  if (!row) return '';\n"
        "  if (row.isCons) return 'Consensus';\n"
        "  if (row.isCollapsed) return row.name || 'Collapsed';\n"
        "  return row.gid || _alnPercIdLabel(row);\n"
        "}\n"
        "\n"
        "function _alnBuildPercIdRowState(rows, seqs) {\n"
        "  const seqState = _alnBuildPercIdState(seqs || []);\n"
        "  const gidToIdx = new Map();\n"
        "  (seqState.seqs || []).forEach((seq, idx) => {\n"
        "    const gid = seq.gid || _alnHeaderGeneId(seq.name || '');\n"
        "    if (gid) gidToIdx.set(gid, idx);\n"
        "  });\n"
        "  const geneMeans = [];\n"
        "  const geneMeanById = new Map();\n"
        "  for (let i = 0; i < seqState.seqs.length; i++) {\n"
        "    let sum = 0;\n"
        "    let count = 0;\n"
        "    for (let j = 0; j < seqState.seqs.length; j++) {\n"
        "      if (i === j) continue;\n"
        "      const frac = seqState.matrix[i][j];\n"
        "      if (!Number.isFinite(frac)) continue;\n"
        "      sum += frac;\n"
        "      count += 1;\n"
        "    }\n"
        "    const gid = seqState.seqs[i]?.gid || _alnHeaderGeneId(seqState.seqs[i]?.name || '');\n"
        "    const mean = count ? (sum / count) : null;\n"
        "    if (gid) geneMeanById.set(gid, mean);\n"
        "    if (Number.isFinite(mean)) geneMeans.push(mean);\n"
        "  }\n"
        "  geneMeans.sort((a, b) => a - b);\n"
        "  const entries = (rows || []).map(row => {\n"
        "    const genes = row && !row.isCons\n"
        "      ? [...new Set((row.percIdGenes || (row.gid ? [row.gid] : [])).filter(Boolean))]\n"
        "      : [];\n"
        "    return {\n"
        "      label: _alnPercIdRowLabel(row),\n"
        "      genes,\n"
        "      isCons: !!(row && row.isCons),\n"
        "      isCollapsed: !!(row && row.isCollapsed),\n"
        "    };\n"
        "  });\n"
        "  function rowMean(genes) {\n"
        "    const vals = genes.map(gid => geneMeanById.get(gid)).filter(v => Number.isFinite(v));\n"
        "    return vals.length ? vals.reduce((a, b) => a + b, 0) / vals.length : null;\n"
        "  }\n"
        "  function agg(aGenes, bGenes) {\n"
        "    if (!aGenes.length || !bGenes.length) return null;\n"
        "    let sum = 0;\n"
        "    let count = 0;\n"
        "    for (const ga of aGenes) {\n"
        "      const ia = gidToIdx.get(ga);\n"
        "      if (ia == null) continue;\n"
        "      for (const gb of bGenes) {\n"
        "        const ib = gidToIdx.get(gb);\n"
        "        if (ib == null) continue;\n"
        "        if (ga === gb) continue;\n"
        "        const frac = seqState.matrix[ia][ib];\n"
        "        if (!Number.isFinite(frac)) continue;\n"
        "        sum += frac;\n"
        "        count += 1;\n"
        "      }\n"
        "    }\n"
        "    if (count) return sum / count;\n"
        "    return aGenes.length === 1 && bGenes.length === 1 && aGenes[0] === bGenes[0] ? 1 : null;\n"
        "  }\n"
        "  const nRows = entries.length;\n"
        "  const matrix = Array.from({length: nRows}, () => Array(nRows).fill(null));\n"
        "  const values = [];\n"
        "  const rowStats = entries.map(entry => {\n"
        "    const mean = rowMean(entry.genes);\n"
        "    return {\n"
        "      mean,\n"
        "      percentile: Number.isFinite(mean) ? _alnPercentileRank(geneMeans, mean) : null,\n"
        "      genes: entry.genes,\n"
        "    };\n"
        "  });\n"
        "  for (let r = 0; r < nRows; r++) {\n"
        "    for (let c = r; c < nRows; c++) {\n"
        "      const frac = agg(entries[r].genes, entries[c].genes);\n"
        "      matrix[r][c] = frac;\n"
        "      matrix[c][r] = frac;\n"
        "      if (Number.isFinite(frac) && r !== c) values.push(frac);\n"
        "    }\n"
        "  }\n"
        "  values.sort((a, b) => a - b);\n"
        "  const percentileScale = {\n"
        "    values,\n"
        "    q05: _alnQuantile(values, 0.05),\n"
        "    q25: _alnQuantile(values, 0.25),\n"
        "    q50: _alnQuantile(values, 0.50),\n"
        "    q75: _alnQuantile(values, 0.75),\n"
        "    q95: _alnQuantile(values, 0.95),\n"
        "  };\n"
        "  return {entries, labels: entries.map(entry => entry.label), matrix, rowStats, percentileScale, n: nRows};\n"
        "}\n"
        "\n"
        "function _alnPercIdColor(frac, percentileScale) {\n"
        "  if (!Number.isFinite(frac)) return '#eef2f6';\n"
        "  const t = _alnPercentileRank((percentileScale && percentileScale.values) || [], frac);\n"
        "  return _alnSpectralColor(t);\n"
        "}\n"
        "\n"
        "function _alnShortPercIdLabel(label) {\n"
        "  const text = String(label || '');\n"
        "  return text.length > 28 ? text.slice(0, 27) + '\\u2026' : text;\n"
        "}\n"
        "\n"
        "function _alnDrawDockedPercIdHeatmap(ctx, x0, rows, state, cH, totalH) {\n"
        "  if (!ctx || !state || !rows || !rows.length) return;\n"
        "  const side = rows.length * cH;\n"
        "  ctx.fillStyle = '#fbfcfd';\n"
        "  ctx.fillRect(x0, 0, side, totalH);\n"
        "  for (let r = 0; r < rows.length; r++) {\n"
        "    for (let c = 0; c < rows.length; c++) {\n"
        "      const frac = state.matrix[r]?.[c];\n"
        "      const x = x0 + c * cH;\n"
        "      const y = r * cH;\n"
        "      if (Number.isFinite(frac)) {\n"
        "        ctx.fillStyle = _alnPercIdColor(frac, state.percentileScale);\n"
        "      } else {\n"
        "        ctx.fillStyle = rows[r]?.isCons || rows[c]?.isCons ? '#eef3f7' : '#f6f8fa';\n"
        "      }\n"
        "      ctx.fillRect(x, y, cH, cH);\n"
        "    }\n"
        "  }\n"
        "  ctx.strokeStyle = 'rgba(255,255,255,.65)';\n"
        "  ctx.lineWidth = 1;\n"
        "  for (let i = 1; i < rows.length; i++) {\n"
        "    const off = i * cH + 0.5;\n"
        "    ctx.beginPath(); ctx.moveTo(x0 + off, 0); ctx.lineTo(x0 + off, totalH); ctx.stroke();\n"
        "    ctx.beginPath(); ctx.moveTo(x0, off); ctx.lineTo(x0 + side, off); ctx.stroke();\n"
        "  }\n"
        "  ctx.strokeStyle = '#93a1af';\n"
        "  ctx.strokeRect(x0 + 0.5, 0.5, Math.max(0, side - 1), Math.max(0, totalH - 1));\n"
        "}\n"
        "\n"
        "function _alnEnsurePercIdModal() {\n"
        "  let modal = document.getElementById('aln-percid-modal');\n"
        "  if (modal) return modal;\n"
        "  modal = document.createElement('div');\n"
        "  modal.id = 'aln-percid-modal';\n"
        "  modal.style.cssText = 'position:fixed;inset:0;display:none;background:rgba(22,26,31,.45);z-index:10020;padding:18px;box-sizing:border-box;';\n"
        "  modal.innerHTML = `\n"
        "<div id='aln-percid-dialog' style='position:absolute;width:min(96vw,1220px);height:min(92vh,980px);display:flex;flex-direction:column;background:#fff;border:1px solid #b9c2cb;border-radius:10px;box-shadow:0 16px 40px rgba(0,0,0,.24);overflow:hidden'>\n"
        "  <div id='aln-percid-head' style='display:flex;align-items:center;justify-content:space-between;gap:16px;padding:10px 14px;border-bottom:1px solid #e3e8ed;background:#f7f9fb;cursor:move;user-select:none'>\n"
        "    <div>\n"
        "      <div style='font-size:13px;font-weight:700;color:#1b2733'>PercID heatmap</div>\n"
        "      <div id='aln-percid-subtitle' style='font-size:11px;color:#677381'></div>\n"
        "    </div>\n"
        "    <div style='display:flex;align-items:center;gap:12px'>\n"
        "      <div id='aln-percid-hover' style='font-size:11px;color:#44515e;min-height:14px;white-space:nowrap'></div>\n"
        "      <button id='aln-percid-dock-btn' type='button' onclick='alnTogglePercIdDock()' style='padding:4px 10px;font-size:11px;border:1px solid #98a3af;border-radius:5px;background:#fff;color:#44515e;cursor:pointer'>Append Right</button>\n"
        "      <button type='button' onclick='alnClosePercIdPopup()' style='padding:4px 10px;font-size:11px;border:1px solid #98a3af;border-radius:5px;background:#fff;color:#44515e;cursor:pointer'>Close</button>\n"
        "    </div>\n"
        "  </div>\n"
        "  <div id='aln-percid-body' style='flex:1;min-height:0;overflow:auto;background:#fff;padding:14px'>\n"
        "    <canvas id='aln-percid-canvas' style='display:block'></canvas>\n"
        "  </div>\n"
        "  <div style='padding:8px 14px;border-top:1px solid #e3e8ed;background:#fafbfd;font-size:11px;color:#677381'>Gap-gap columns are ignored. Gap-residue pairs count as mismatches.</div>\n"
        "</div>`;\n"
        "  document.body.appendChild(modal);\n"
        "  modal.addEventListener('mousedown', ev => {\n"
        "    if (ev.target === modal) alnClosePercIdPopup();\n"
        "  });\n"
        "  document.addEventListener('keydown', ev => {\n"
        "    if (ev.key === 'Escape') alnClosePercIdPopup();\n"
        "  });\n"
        "  const canvas = modal.querySelector('#aln-percid-canvas');\n"
        "  const dialog = modal.querySelector('#aln-percid-dialog');\n"
        "  const head = modal.querySelector('#aln-percid-head');\n"
        "  let dragState = null;\n"
        "  if (canvas) {\n"
        "    canvas.addEventListener('mousemove', _alnHandlePercIdHover);\n"
        "    canvas.addEventListener('mouseleave', () => {\n"
        "      const hover = document.getElementById('aln-percid-hover');\n"
        "      if (hover) hover.textContent = '';\n"
        "    });\n"
        "  }\n"
        "  function _alnClampPercIdDialog() {\n"
        "    if (!dialog) return;\n"
        "    const pad = 18;\n"
        "    const rect = dialog.getBoundingClientRect();\n"
        "    const width = rect.width || Math.min(window.innerWidth - pad * 2, 1220);\n"
        "    const height = rect.height || Math.min(window.innerHeight - pad * 2, 980);\n"
        "    const fallback = {\n"
        "      left: Math.round((window.innerWidth - width) / 2),\n"
        "      top: Math.round((window.innerHeight - height) / 2),\n"
        "    };\n"
        "    const pos = _alnPercIdDialogPos || fallback;\n"
        "    const maxLeft = Math.max(pad, window.innerWidth - width - pad);\n"
        "    const maxTop = Math.max(pad, window.innerHeight - height - pad);\n"
        "    const left = Math.max(pad, Math.min(maxLeft, pos.left));\n"
        "    const top = Math.max(pad, Math.min(maxTop, pos.top));\n"
        "    _alnPercIdDialogPos = {left, top};\n"
        "    dialog.style.left = `${left}px`;\n"
        "    dialog.style.top = `${top}px`;\n"
        "  }\n"
        "  modal._alnClampPercIdDialog = _alnClampPercIdDialog;\n"
        "  modal._alnSyncPercIdDockBtn = () => {\n"
        "    const btn = document.getElementById('aln-percid-dock-btn');\n"
        "    if (!btn) return;\n"
        "    btn.textContent = alnShowPercIdDock ? 'Remove Right' : 'Append Right';\n"
        "    btn.style.background = alnShowPercIdDock ? '#e8f0fe' : '#fff';\n"
        "    btn.style.color = alnShowPercIdDock ? '#1a56c4' : '#44515e';\n"
        "    btn.style.borderColor = alnShowPercIdDock ? '#1a56c4' : '#98a3af';\n"
        "    btn.style.fontWeight = alnShowPercIdDock ? '600' : '';\n"
        "  };\n"
        "  modal._alnSyncPercIdDockBtn();\n"
        "  if (head && dialog) {\n"
        "    head.addEventListener('mousedown', ev => {\n"
        "      if (ev.target && typeof ev.target.closest === 'function' && ev.target.closest('button')) return;\n"
        "      const rect = dialog.getBoundingClientRect();\n"
        "      dragState = {dx: ev.clientX - rect.left, dy: ev.clientY - rect.top};\n"
        "      ev.preventDefault();\n"
        "    });\n"
        "    document.addEventListener('mousemove', ev => {\n"
        "      if (!dragState || modal.style.display !== 'block') return;\n"
        "      _alnPercIdDialogPos = {left: ev.clientX - dragState.dx, top: ev.clientY - dragState.dy};\n"
        "      _alnClampPercIdDialog();\n"
        "    });\n"
        "    document.addEventListener('mouseup', () => { dragState = null; });\n"
        "  }\n"
        "  window.addEventListener('resize', () => {\n"
        "    if (modal.style.display === 'block') {\n"
        "      _alnRenderPercIdHeatmap();\n"
        "      _alnClampPercIdDialog();\n"
        "    }\n"
        "  });\n"
        "  return modal;\n"
        "}\n"
        "\n"
        "function _alnRenderPercIdHeatmap() {\n"
        "  _alnEnsurePercIdModal();\n"
        "  const canvas = document.getElementById('aln-percid-canvas');\n"
        "  const subtitle = document.getElementById('aln-percid-subtitle');\n"
        "  const hover = document.getElementById('aln-percid-hover');\n"
        "  if (!canvas || !subtitle) return;\n"
        "  const seqs = _alnFilteredSeqs();\n"
        "  const rows = (_alnRows || []).slice();\n"
        "  if (!seqs.length || !rows.length) {\n"
        "    _alnPercIdState = null;\n"
        "    subtitle.textContent = 'No sequences available';\n"
        "    if (hover) hover.textContent = '';\n"
        "    canvas.width = 1;\n"
        "    canvas.height = 1;\n"
        "    canvas.style.width = '1px';\n"
        "    canvas.style.height = '1px';\n"
        "    return;\n"
        "  }\n"
        "  const state = _alnBuildPercIdRowState(rows, seqs);\n"
        "  const n = state.labels.length;\n"
        "  const cell = Math.max(8, alnCellH);\n"
        "  const maxLabelChars = state.labels.reduce((mx, label) => Math.max(mx, String(label || '').length), 0);\n"
        "  const leftPad = Math.min(320, Math.max(110, Math.round(maxLabelChars * Math.max(6, cell * 0.72))));\n"
        "  const topPad = leftPad;\n"
        "  const width = leftPad + n * cell + 22;\n"
        "  const height = topPad + n * cell + 22;\n"
        "  const dpr = window.devicePixelRatio || 1;\n"
        "  canvas.width = Math.max(1, Math.ceil(width * dpr));\n"
        "  canvas.height = Math.max(1, Math.ceil(height * dpr));\n"
        "  canvas.style.width = `${width}px`;\n"
        "  canvas.style.height = `${height}px`;\n"
        "  const ctx = canvas.getContext('2d');\n"
        "  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);\n"
        "  ctx.clearRect(0, 0, width, height);\n"
        "  ctx.fillStyle = '#fff';\n"
        "  ctx.fillRect(0, 0, width, height);\n"
        "  ctx.fillStyle = '#f4f7fa';\n"
        "  ctx.fillRect(leftPad, topPad, n * cell, n * cell);\n"
        "  ctx.strokeStyle = '#d7dee6';\n"
        "  ctx.strokeRect(leftPad + 0.5, topPad + 0.5, Math.max(0, n * cell - 1), Math.max(0, n * cell - 1));\n"
        "  for (let r = 0; r < n; r++) {\n"
        "    for (let c = 0; c < n; c++) {\n"
        "      const frac = state.matrix[r][c];\n"
        "      const x = leftPad + c * cell;\n"
        "      const y = topPad + r * cell;\n"
        "      ctx.fillStyle = _alnPercIdColor(frac, state.percentileScale);\n"
        "      ctx.fillRect(x, y, cell, cell);\n"
        "    }\n"
        "  }\n"
        "  ctx.strokeStyle = 'rgba(255,255,255,.55)';\n"
        "  ctx.lineWidth = 1;\n"
        "  for (let i = 1; i < n; i++) {\n"
        "    const off = i * cell + 0.5;\n"
        "    ctx.beginPath();\n"
        "    ctx.moveTo(leftPad + off, topPad);\n"
        "    ctx.lineTo(leftPad + off, topPad + n * cell);\n"
        "    ctx.stroke();\n"
        "    ctx.beginPath();\n"
        "    ctx.moveTo(leftPad, topPad + off);\n"
        "    ctx.lineTo(leftPad + n * cell, topPad + off);\n"
        "    ctx.stroke();\n"
        "  }\n"
        "  if (cell >= 18 && n <= 42) {\n"
        "    const valFont = Math.max(8, Math.min(11, Math.floor(cell * 0.58)));\n"
        "    ctx.font = `600 ${valFont}px sans-serif`;\n"
        "    ctx.textAlign = 'center';\n"
        "    ctx.textBaseline = 'middle';\n"
        "    for (let r = 0; r < n; r++) {\n"
        "      for (let c = 0; c < n; c++) {\n"
        "        const frac = state.matrix[r][c];\n"
        "        if (!Number.isFinite(frac)) continue;\n"
        "        const pct = Math.round(frac * 100);\n"
        "        ctx.fillStyle = frac >= 0.66 ? '#fff' : '#163149';\n"
        "        ctx.fillText(String(pct), leftPad + c * cell + cell / 2, topPad + r * cell + cell / 2 + 0.2);\n"
        "      }\n"
        "    }\n"
        "  }\n"
        "  const labelFont = Math.max(10, Math.min(12, cell + 1));\n"
        "  ctx.font = `600 ${labelFont}px sans-serif`;\n"
        "  ctx.fillStyle = '#344251';\n"
        "  ctx.textBaseline = 'middle';\n"
        "  ctx.textAlign = 'right';\n"
        "  for (let r = 0; r < n; r++) {\n"
        "    ctx.fillText(_alnShortPercIdLabel(state.labels[r]), leftPad - 8, topPad + r * cell + cell / 2);\n"
        "  }\n"
        "  for (let c = 0; c < n; c++) {\n"
        "    ctx.save();\n"
        "    ctx.translate(leftPad + c * cell + cell / 2, topPad - 8);\n"
        "    ctx.rotate(-Math.PI / 3);\n"
        "    ctx.textAlign = 'left';\n"
        "    ctx.textBaseline = 'middle';\n"
        "    ctx.fillText(_alnShortPercIdLabel(state.labels[c]), 0, 0);\n"
        "    ctx.restore();\n"
        "  }\n"
        "  const scale = state.percentileScale || {};\n"
        "  const fmtQ = v => Number.isFinite(v) ? `${(v * 100).toFixed(1)}%` : 'n/a';\n"
        "  subtitle.textContent = `${n} visible rows \\u00b7 row order and row height match the alignment view \\u00b7 percentile-scaled colors (P5 ${fmtQ(scale.q05)}, P50 ${fmtQ(scale.q50)}, P95 ${fmtQ(scale.q95)})`;\n"
        "  if (hover) hover.textContent = 'Hover a cell for exact identity';\n"
        "  state.n = n;\n"
        "  state.cell = cell;\n"
        "  state.leftPad = leftPad;\n"
        "  state.topPad = topPad;\n"
        "  _alnPercIdState = state;\n"
        "  const modal = document.getElementById('aln-percid-modal');\n"
        "  if (modal && typeof modal._alnSyncPercIdDockBtn === 'function') modal._alnSyncPercIdDockBtn();\n"
        "}\n"
        "\n"
        "function _alnHandlePercIdHover(ev) {\n"
        "  const hover = document.getElementById('aln-percid-hover');\n"
        "  const canvas = document.getElementById('aln-percid-canvas');\n"
        "  const state = _alnPercIdState;\n"
        "  if (!hover || !canvas || !state) return;\n"
        "  const rect = canvas.getBoundingClientRect();\n"
        "  const x = ev.clientX - rect.left;\n"
        "  const y = ev.clientY - rect.top;\n"
        "  const col = Math.floor((x - state.leftPad) / state.cell);\n"
        "  const row = Math.floor((y - state.topPad) / state.cell);\n"
        "  if (row < 0 || col < 0 || row >= state.n || col >= state.n) {\n"
        "    hover.textContent = '';\n"
        "    return;\n"
        "  }\n"
        "  const frac = state.matrix[row][col];\n"
        "  const a = state.labels[row] || `seq ${row + 1}`;\n"
        "  const b = state.labels[col] || `seq ${col + 1}`;\n"
        "  if (!Number.isFinite(frac)) {\n"
        "    hover.textContent = `${a} vs ${b}: n/a`;\n"
        "    return;\n"
        "  }\n"
        "  hover.textContent = `${a} vs ${b}: ${(frac * 100).toFixed(1)}%`;\n"
        "}\n"
        "\n"
        "function alnTogglePercIdDock() {\n"
        "  alnShowPercIdDock = !alnShowPercIdDock;\n"
        "  const modal = document.getElementById('aln-percid-modal');\n"
        "  if (modal && typeof modal._alnSyncPercIdDockBtn === 'function') modal._alnSyncPercIdDockBtn();\n"
        "  renderAlignment();\n"
        "}\n"
        "\n"
        "function alnOpenPercIdPopup() {\n"
        "  if (!_alnFilteredSeqs().length) return;\n"
        "  const modal = _alnEnsurePercIdModal();\n"
        "  modal.style.display = 'block';\n"
        "  if (typeof modal._alnClampPercIdDialog === 'function') modal._alnClampPercIdDialog();\n"
        "  if (typeof modal._alnSyncPercIdDockBtn === 'function') modal._alnSyncPercIdDockBtn();\n"
        "  _alnSetPercIdBtn(true);\n"
        "  _alnRenderPercIdHeatmap();\n"
        "  if (typeof modal._alnClampPercIdDialog === 'function') modal._alnClampPercIdDialog();\n"
        "}\n"
        "\n"
        "function alnClosePercIdPopup() {\n"
        "  const modal = document.getElementById('aln-percid-modal');\n"
        "  if (modal) modal.style.display = 'none';\n"
        "  const hover = document.getElementById('aln-percid-hover');\n"
        "  if (hover) hover.textContent = '';\n"
        "  _alnSetPercIdBtn(false);\n"
        "}\n",
    )
    tmpl = tmpl.replace(
        "  const info = document.getElementById('aln-info');\n"
        "  const ogInfo = Object.keys(alnOgMap).length ? ` \\u00b7 ${[...new Set(Object.values(alnOgMap))].length} OGs` : '';\n"
        "  const sourceInfo = _alnHasGeneRaxSource() && _alnHasPrevSource()\n"
        "    ? ` \\u00b7 ${alnUseGeneRax ? 'GeneRax' : 'IQ-TREE'} source`\n"
        "    : (_alnHasGeneRaxSource() ? ' \\u00b7 GeneRax source' : (_alnHasPrevSource() ? ' \\u00b7 IQ-TREE source' : ''));\n"
        "  info.textContent = `${seqs.length} seqs \\u00d7 ${nCols} cols${ogInfo}${sourceInfo}${alnShowSeqPanel ? '' : ' \\u00b7 alignment hidden'}`;\n",
        "  const info = document.getElementById('aln-info');\n"
        "  const ogInfo = Object.keys(alnOgMap).length ? ` \\u00b7 ${[...new Set(Object.values(alnOgMap))].length} OGs` : '';\n"
        "  const sourceInfo = _alnHasGeneRaxSource() && _alnHasPrevSource()\n"
        "    ? ` \\u00b7 ${alnUseGeneRax ? 'GeneRax' : 'IQ-TREE'} source`\n"
        "    : (_alnHasGeneRaxSource() ? ' \\u00b7 GeneRax source' : (_alnHasPrevSource() ? ' \\u00b7 IQ-TREE source' : ''));\n"
        "  const hlInfo = spHlSet\n"
        "    ? ` \\u00b7 ${rows.filter(row => !row.isCons && row.spMatched).length} highlighted`\n"
        "    : '';\n"
        "  const dsInfo = alnShowDSNodes ? ' \\u00b7 POSSVM D/S' : '';\n"
        "  info.textContent = `${seqs.length} seqs \\u00d7 ${nCols} cols${ogInfo}${sourceInfo}${hlInfo}${dsInfo}${alnShowSeqPanel ? '' : ' \\u00b7 alignment hidden'}`;\n"
        "  const percIdModal = document.getElementById('aln-percid-modal');\n"
        "  if (percIdModal && percIdModal.style.display === 'block') _alnRenderPercIdHeatmap();\n",
    )

    tmpl = tmpl.replace(
        "    if (_drag === 'aln-resize-bar') {\n"
        "      alnNameW = Math.max(60, _startW + delta);\n",
        "    if (_drag === 'aln-resize-bar') {\n"
        "      _alnAutoNameWidth = false;\n"
        "      alnNameW = Math.max(60, _startW + delta);\n",
    )
    tmpl = tmpl.replace(
        "function _alnHandleNameHover(ev) {\n",
        "function _alnNearestTreeEventNode(x, y) {\n"
        "  if (!alnShowDSNodes || !_alnTreeEventNodes || !_alnTreeEventNodes.length) return null;\n"
        "  let best = null;\n"
        "  let bestD2 = 81;\n"
        "  for (const nd of _alnTreeEventNodes) {\n"
        "    const dx = x - nd.x;\n"
        "    const dy = y - nd.y;\n"
        "    const d2 = dx * dx + dy * dy;\n"
        "    if (d2 <= bestD2) { best = nd; bestD2 = d2; }\n"
        "  }\n"
        "  return best;\n"
        "}\n"
        "\n"
        "function _alnHandleNameHover(ev) {\n",
    )
    tmpl = tmpl.replace(
        "  const row = _alnRows[r];\n"
        "  if (row.isCollapsed) {\n",
        "  const row = _alnRows[r];\n"
        "  const treeNodeHit = x < dims.treeW ? _alnNearestTreeEventNode(x, y) : null;\n"
        "  if (treeNodeHit) {\n"
        "    if (treeNodeHit.event === 'D' && treeNodeHit.dupSpecies && treeNodeHit.dupSpecies.length) {\n"
        "      _alnShowTip(ev, `Duplication node · overlap species: ${treeNodeHit.dupSpecies.join(', ')}`);\n"
        "    } else {\n"
        "      _alnShowTip(ev, `${treeNodeHit.event === 'D' ? 'Duplication' : 'Speciation'} node · ${treeNodeHit.genes.length} downstream seqs`);\n"
        "    }\n"
        "  } else if (row.isCollapsed) {\n",
    )
    tmpl = tmpl.replace(
        "document.getElementById(\"tree-count\").textContent =\n"
        "  TREE_INDEX.length+\" gene tree\"+(TREE_INDEX.length!==1?\"s\":\"\");\n"
        "\n"
        "document.getElementById(\"btn-og-labels\").classList.toggle(\"active-btn\", showOGLabels);\n"
        "\n"
        "if (!HAVE_ALIGNMENTS) {\n"
        "  const alignBtn = document.getElementById(\"tab-btn-align\");\n"
        "  if (alignBtn) alignBtn.style.display = \"none\";\n"
        "}\n"
        "if (!HAVE_ARCHITECTURES) {\n"
        "  const archBtn = document.querySelector('.tab-btn[data-tab=\"architectures\"]');\n"
        "  if (archBtn) archBtn.style.display = \"none\";\n"
        "}\n"
        "\n"
        "if (hasHeatmapData || TREE_INDEX.length > 0) {\n"
        "  switchTab(\"sptree\");\n"
        "} else {\n"
        "  document.getElementById(\"pane-heatmap\").innerHTML =\n"
        "    '<div style=\"padding:40px;color:#999;text-align:center\">No data found.<br>'+\n"
        "    'Pass <code>--possvm_dir</code> and/or <code>--search_dir</code> / <code>--cluster_dir</code>.</div>';\n"
        "}\n",
        "document.getElementById(\"tree-count\").textContent =\n"
        "  TREE_INDEX.length+\" alignment\"+(TREE_INDEX.length!==1?\"s\":\"\");\n"
        "\n"
        "document.getElementById(\"btn-og-labels\").classList.toggle(\"active-btn\", showOGLabels);\n"
        "\n"
        "const alignBtn = document.getElementById(\"tab-btn-align\");\n"
        "if (alignBtn) alignBtn.style.display = HAVE_ALIGNMENTS ? \"\" : \"none\";\n"
        "const archBtn = document.querySelector('.tab-btn[data-tab=\"architectures\"]');\n"
        "if (archBtn) archBtn.style.display = \"none\";\n"
        "const miniBtn = document.getElementById(\"btn-aln-mini-sp\");\n"
        "if (miniBtn) miniBtn.style.display = HAVE_SPECIES_TREE ? \"\" : \"none\";\n"
        "\n"
        "if (TREE_INDEX.length > 0 && HAVE_ALIGNMENTS) {\n"
        "  switchTab(\"align\");\n"
        "} else {\n"
        "  document.getElementById(\"pane-align\").innerHTML =\n"
        "    '<div style=\"padding:40px;color:#999;text-align:center\">No alignment/tree data found.</div>';\n"
        "  document.getElementById(\"pane-align\").classList.add(\"active\");\n"
        "}\n",
    )

    # Final cleanup pass for downstream-tip dimming after D-node selection.
    tmpl = tmpl.replace(
        "          nCtx.fillStyle = row.treeFocus ? '#b42318' : '#222';\n",
        "          nCtx.fillStyle = row.treeFocus ? _alnTreeFocusColor(row.treeFocusGroup) : (selectedGenes && !row.treeSelected ? '#8b95a1' : '#222');\n",
    )
    tmpl = tmpl.replace(
        "      const rowBg = row.treeFocus\n"
        "        ? _alnTreeFocusFill(row.treeFocusGroup, 0.12)\n"
        "        : (spHlSet\n"
        "          ? (row.spMatched ? '#fff3bf' : (striped ? '#f2f3f5' : '#fafbfc'))\n"
        "          : (striped ? '#f7f7f7' : ''));\n",
        "      const rowBg = row.treeFocus\n"
        "        ? '#ffe3e0'\n"
        "        : (selectedGenes && !row.treeSelected)\n"
        "          ? (striped ? '#f4f6f8' : '#fafbfc')\n"
        "          : (spHlSet\n"
        "            ? (row.spMatched ? '#fff3bf' : (striped ? '#f2f3f5' : '#fafbfc'))\n"
        "            : (striped ? '#f7f7f7' : ''));\n",
    )
    tmpl = tmpl.replace(
        "          nCtx.fillStyle = row.treeFocus\n"
        "            ? _alnTreeFocusColor(row.treeFocusGroup)\n"
        "            : (isRefSeq ? '#cc0000' : '#333');\n",
        "          nCtx.fillStyle = row.treeFocus\n"
        "            ? _alnTreeFocusColor(row.treeFocusGroup)\n"
        "            : (selectedGenes && !row.treeSelected ? '#8b95a1' : (isRefSeq ? '#cc0000' : '#333'));\n",
    )
    tmpl = tmpl.replace(
        "        sCtx.fillStyle = row.treeFocus ? _alnTreeFocusColor(row.treeFocusGroup) : '#444';\n",
        "        sCtx.fillStyle = row.treeFocus ? _alnTreeFocusColor(row.treeFocusGroup) : (selectedGenes && !row.treeSelected ? '#8b95a1' : '#444');\n",
    )
    tmpl = tmpl.replace(
        "      if (!isConsRow && row.treeFocus) {\n"
        "        sCtx.fillStyle = _alnTreeFocusFill(row.treeFocusGroup, 0.18);\n"
        "        sCtx.fillRect(0, y0, totalW, cH);\n"
        "        const prevTreeSel = r > 0 && !!(rows[r - 1] && rows[r - 1].treeFocus);\n"
        "        const nextTreeSel = r + 1 < rows.length && !!(rows[r + 1] && rows[r + 1].treeFocus);\n"
        "        sCtx.strokeStyle = '#111';\n"
        "        sCtx.lineWidth = 1.0;\n"
        "        if (!prevTreeSel) {\n"
        "          sCtx.beginPath();\n"
        "          sCtx.moveTo(0, y0 + 0.5);\n"
        "          sCtx.lineTo(totalW, y0 + 0.5);\n"
        "          sCtx.stroke();\n"
        "        }\n"
        "        if (!nextTreeSel) {\n"
        "          sCtx.beginPath();\n"
        "          sCtx.moveTo(0, y0 + cH - 0.5);\n"
        "          sCtx.lineTo(totalW, y0 + cH - 0.5);\n"
        "          sCtx.stroke();\n"
        "        }\n"
        "      }\n",
        "      if (!isConsRow && row.treeFocus) {\n"
        "        sCtx.fillStyle = _alnTreeFocusFill(row.treeFocusGroup, 0.18);\n"
        "        sCtx.fillRect(0, y0, totalW, cH);\n"
        "        const prevTreeSel = r > 0 && !!(rows[r - 1] && rows[r - 1].treeFocus);\n"
        "        const nextTreeSel = r + 1 < rows.length && !!(rows[r + 1] && rows[r + 1].treeFocus);\n"
        "        sCtx.strokeStyle = '#111';\n"
        "        sCtx.lineWidth = 1.0;\n"
        "        if (!prevTreeSel) {\n"
        "          sCtx.beginPath();\n"
        "          sCtx.moveTo(0, y0 + 0.5);\n"
        "          sCtx.lineTo(totalW, y0 + 0.5);\n"
        "          sCtx.stroke();\n"
        "        }\n"
        "        if (!nextTreeSel) {\n"
        "          sCtx.beginPath();\n"
        "          sCtx.moveTo(0, y0 + cH - 0.5);\n"
        "          sCtx.lineTo(totalW, y0 + cH - 0.5);\n"
        "          sCtx.stroke();\n"
        "        }\n"
        "      } else if (!isConsRow && selectedGenes && !row.treeSelected) {\n"
        "        sCtx.fillStyle = 'rgba(255,255,255,0.62)';\n"
        "        sCtx.fillRect(0, y0, totalW, cH);\n"
        "      }\n",
    )

    old_click = (
        "function _alnHandleNameClick(ev) {\n"
        "  const canvas = document.getElementById('aln-names-canvas');\n"
        "  if (!canvas || !_alnRows.length) return;\n"
        "  const rect = canvas.getBoundingClientRect();\n"
        "  const y = ev.clientY - rect.top;\n"
        "  const cH = alnCellH;\n"
        "  const dims = _alnLayoutMetrics();\n"
        "  const x = ev.clientX - rect.left;\n"
        "  const r = Math.floor(y / cH);\n"
        "  if (r < 0 || r >= _alnRows.length) return;\n"
        "  const row = _alnRows[r];\n"
        "  const onOgBand = x >= dims.treeW && x < dims.treeW + dims.ogW;\n"
        "  const idCol = dims.byKey.id;\n"
        "  if (!row.isCons && !row.isCollapsed && idCol && x >= idCol.x0 && x < idCol.x1) {\n"
        "    const gid=row.gid || row.name.split(/\\s/)[0];\n"
        "    _copyTextToClipboard(gid).then(()=>_alnShowTip(ev, `Copied: ${gid}`));\n"
        "    setTimeout(_alnHideTip, 900);\n"
        "    return;\n"
        "  }\n"
        "  if ((onOgBand || row.isCollapsed) && row.og) _alnShowOgPopup(ev, row.og);\n"
        "}\n"
    )
    new_click = (
        "function _alnHandleNameClick(ev) {\n"
        "  const canvas = document.getElementById('aln-names-canvas');\n"
        "  if (!canvas || !_alnRows.length) return;\n"
        "  const rect = canvas.getBoundingClientRect();\n"
        "  const y = ev.clientY - rect.top;\n"
        "  const cH = alnCellH;\n"
        "  const dims = _alnLayoutMetrics();\n"
        "  const x = ev.clientX - rect.left;\n"
        "  const r = Math.floor(y / cH);\n"
        "  if (r < 0 || r >= _alnRows.length) return;\n"
        "  const row = _alnRows[r];\n"
        "  const treeNodeHit = x < dims.treeW ? _alnNearestTreeEventNode(x, y) : null;\n"
        "  if (treeNodeHit) {\n"
        "    _alnSelectedTreeNodeKey = _alnSelectedTreeNodeKey === treeNodeHit.key ? null : treeNodeHit.key;\n"
        "    renderAlignment();\n"
        "    if (treeNodeHit.event === 'D' && treeNodeHit.dupSpecies && treeNodeHit.dupSpecies.length) {\n"
        "      _alnShowTip(ev, `Duplication node · overlap species: ${treeNodeHit.dupSpecies.join(', ')}`);\n"
        "    } else {\n"
        "      _alnShowTip(ev, `${treeNodeHit.event === 'D' ? 'Duplication' : 'Speciation'} node · ${treeNodeHit.genes.length} downstream seqs`);\n"
        "    }\n"
        "    setTimeout(_alnHideTip, 900);\n"
        "    return;\n"
        "  }\n"
        "  if (x < dims.treeW && _alnSelectedTreeNodeKey) {\n"
        "    _alnSelectedTreeNodeKey = null;\n"
        "    renderAlignment();\n"
        "    return;\n"
        "  }\n"
        "  const onOgBand = x >= dims.treeW && x < dims.treeW + dims.ogW;\n"
        "  const idCol = dims.byKey.id;\n"
        "  if (!row.isCons && !row.isCollapsed && idCol && x >= idCol.x0 && x < idCol.x1) {\n"
        "    const gid=row.gid || row.name.split(/\\s/)[0];\n"
        "    _copyTextToClipboard(gid).then(()=>_alnShowTip(ev, `Copied: ${gid}`));\n"
        "    setTimeout(_alnHideTip, 900);\n"
        "    return;\n"
        "  }\n"
        "  if ((onOgBand || row.isCollapsed) && row.og) _alnShowOgPopup(ev, row.og);\n"
        "}\n"
    )
    if old_click in tmpl:
        tmpl = tmpl.replace(old_click, new_click, 1)

    return tmpl


def build_context(args) -> dict:
    alignment_path = Path(args.alignment)
    gene_tree_path = Path(args.gene_tree)
    species_tree_path = Path(args.species_tree) if args.species_tree else None

    alignment_ids = _parse_fasta_ids(alignment_path)
    alignment_id_set = set(alignment_ids)
    alignment_species = sorted({get_species_prefix(gid) for gid in alignment_ids if gid})

    tree_dict, tree_species, original_tree_ids = _load_gene_tree(
        gene_tree_path, alignment_id_set
    )
    tree_id_set = set(original_tree_ids)
    missing_in_tree = sorted(alignment_id_set - tree_id_set)
    if missing_in_tree:
        print(
            f"Alignment contains {len(missing_in_tree)} sequence IDs not present in the gene tree; "
            "they will be shown after the tree-ordered sequences.",
            file=sys.stderr,
        )

    species_order: list[str] = []
    species_tree_dict: dict = {}
    newick_raw = ""
    clade_groupings: list = []
    if species_tree_path and species_tree_path.exists():
        species_order, species_tree_dict = load_tree_data(str(species_tree_path))
        newick_raw = species_tree_path.read_text(encoding="utf-8").strip()
        clade_groupings = parse_clade_groupings(species_tree_path)

    if not species_order:
        species_order = list(dict.fromkeys(tree_species + alignment_species))

    all_species = sorted(set(species_order) | set(tree_species) | set(alignment_species))

    stem = args.stem or alignment_path.name
    for suffix in (".aln.fasta", ".fasta", ".fa", ".faa", ".fas"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break

    record = {
        "id": stem,
        "hg": stem,
        "family": args.family or stem,
        "prefix": "",
        "source": "standalone",
        "n_leaves": len(alignment_ids),
        "species": all_species,
        "og_names": [],
        "n_ogs": 0,
        "has_prev": False,
        "class": "",
    }

    return {
        "title": args.title or stem,
        "stem": stem,
        "record": record,
        "alignment_b64_gz": _gzip_bytes(alignment_path.read_bytes()),
        "detail_b64_gz": _gzip_json({"tree": tree_dict, "ogs": {}}),
        "species_order": species_order,
        "species_tree_dict": species_tree_dict,
        "all_species": all_species,
        "clade_groupings": clade_groupings,
        "newick_raw": newick_raw,
        "have_species_tree": bool(species_tree_dict),
    }


_ADAPTED_TEMPLATE: Optional[str] = None


def _get_template() -> str:
    global _ADAPTED_TEMPLATE
    if _ADAPTED_TEMPLATE is None:
        _ADAPTED_TEMPLATE = _adapt_template(_r2.HTML_TEMPLATE)
    return _ADAPTED_TEMPLATE


def render_html(ctx: dict) -> str:
    tag_id = _html.escape(ctx["stem"], quote=True)
    aln_script = (
        f'<script type="application/json" id="alndata-{tag_id}">'
        + json.dumps({"gz": ctx["alignment_b64_gz"]}, separators=(",", ":"))
        + "</script>"
    )
    lazy_script = (
        f'<script type="application/json" id="treedata-{tag_id}">'
        + json.dumps({"gz": ctx["detail_b64_gz"]}, separators=(",", ":"))
        + "</script>"
    )

    return (
        _get_template()
        .replace("%%VIEWER_TITLE%%", _html.escape(ctx["title"]))
        .replace("%%LAZY_SCRIPTS%%", lazy_script)
        .replace("%%ALN_SCRIPTS%%", aln_script)
        .replace("%%PROTEIN_DOMAIN_SCRIPT%%", '<script type="application/json" id="protein-domain-data">{"gz":null}</script>')
        .replace("%%ARCHITECTURE_SCRIPT%%", '<script type="application/json" id="architecture-data">{"gz":null}</script>')
        .replace("%%SPECIES_ORDER%%", json.dumps(ctx["species_order"]))
        .replace("%%TREE_DATA%%", json.dumps(ctx["species_tree_dict"]))
        .replace("%%FAMILY_DATA%%", "[]")
        .replace("%%HG_DATA%%", "[]")
        .replace("%%TREE_INDEX_JSON%%", json.dumps([ctx["record"]]))
        .replace("%%SPECIES_JSON%%", json.dumps(ctx["all_species"]))
        .replace("%%CLADE_DATA_JSON%%", json.dumps(ctx["clade_groupings"]))
        .replace("%%NEWICK_RAW%%", json.dumps(ctx["newick_raw"]))
        .replace("%%FAMILY_INFO_JSON%%", "[]")
        .replace("%%HAVE_GENERAX_JSON%%", "false")
        .replace("%%NO_TREE_GENES_JSON%%", "{}")
        .replace("%%DOMAIN_DATA_JSON%%", "{}")
        .replace("%%GENE_META_JSON%%", "{}")
        .replace("%%REFNAME_MAP_JSON%%", "{}")
        .replace("%%SPECIES_INFO_JSON%%", "{}")
        .replace("%%SPECIES_GROUPS_JSON%%", "{}")
        .replace("%%SPECIES_IMAGES_JSON%%", "{}")
        .replace("%%HAVE_ALIGNMENTS_JSON%%", "true")
        .replace("%%HAVE_PROTEIN_DOMAINS_JSON%%", "false")
        .replace("%%HAVE_ARCHITECTURES_JSON%%", "false")
        .replace("%%HAVE_SPECIES_TREE_JSON%%", json.dumps(ctx["have_species_tree"]))
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate a standalone HTML viewer for one alignment and one gene tree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-a",
        "--alignment",
        required=True,
        help="Aligned FASTA file to render",
    )
    parser.add_argument(
        "-g",
        "--gene-tree",
        required=True,
        help="Gene tree Newick whose leaf IDs match the alignment headers",
    )
    parser.add_argument(
        "--species-tree",
        default=None,
        help="Optional species tree Newick for the mini clade-filter panel",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Viewer title shown in the header",
    )
    parser.add_argument(
        "--family",
        default=None,
        help="Optional family label to inject into the synthetic alignment record",
    )
    parser.add_argument(
        "--stem",
        default=None,
        help="Synthetic HG/alignment ID used inside the viewer",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output HTML path (default: <alignment stem>.alignment_viewer.html)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.output is None:
        alignment_path = Path(args.alignment)
        stem = alignment_path.name
        for suffix in (".aln.fasta", ".fasta", ".fa", ".faa", ".fas"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        output = Path.cwd() / f"{stem}.alignment_viewer.html"
    else:
        output = Path(args.output)

    ctx = build_context(args)
    output.write_text(render_html(ctx), encoding="utf-8")
    print(f"Viewer written to {output}", file=sys.stderr)


if __name__ == "__main__":
    main()
