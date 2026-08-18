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

from __future__ import annotations

import argparse
import html as _html
import json
import math
import re
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
    ansi_re = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
    try:
        with open(csv_path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 7:
                    continue
                family = ansi_re.sub("", cols[0]).strip()
                if not family:
                    continue
                yield {
                    "family": family,
                    "pfam_raw": ansi_re.sub("", cols[1]).strip(),
                    "category": ansi_re.sub("", cols[5]).strip(),
                    "cls": ansi_re.sub("", cols[6]).strip(),
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


def _family_info_allows(prefix: str, family: str, family_info: dict | None) -> bool:
    if not family_info:
        return True
    return prefix in family_info or family in family_info


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


def load_domain_hits(search_dir: Path, family_info: dict | None = None) -> dict:
    """Return {gene_id: [{name,start,end}, ...]} from *.domains.csv files."""
    hits: dict = defaultdict(list)
    if not search_dir.is_dir():
        return {}
    for domains_file in sorted(search_dir.glob("*.domains.csv")):
        track = domains_file.name[: -len(".domains.csv")]
        parts = track.split(".", 1)
        if len(parts) < 2:
            continue
        pref, family = parts
        if not _family_info_allows(pref, family, family_info):
            continue
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


def build_domain_spans(domain_hits: dict) -> dict:
    """Return compact {gene_id: [start,end]} spans for alignment-range rendering."""
    spans: dict = {}
    for gene_id, hits in domain_hits.items():
        if not hits:
            continue
        starts = [int(h["start"]) for h in hits if "start" in h]
        ends = [int(h["end"]) for h in hits if "end" in h]
        if not starts or not ends:
            continue
        spans[gene_id] = [min(starts), max(ends)]
    return spans


def load_search_gene_lengths(search_dir: Path, family_info: dict | None = None) -> dict:
    """Return {gene_id: protein_length} from *.seqs.fasta.fai sidecar files."""
    lengths: dict = {}
    if not search_dir.is_dir():
        return lengths
    for fai_file in sorted(search_dir.glob("*.seqs.fasta.fai")):
        track = fai_file.name[: -len(".seqs.fasta.fai")]
        parts = track.split(".", 1)
        if len(parts) < 2:
            continue
        pref, family = parts
        if not _family_info_allows(pref, family, family_info):
            continue
        try:
            with open(fai_file) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 2:
                        continue
                    gene_id = cols[0].strip()
                    if not gene_id:
                        continue
                    try:
                        lengths[gene_id] = int(cols[1])
                    except ValueError:
                        continue
        except OSError:
            continue
    return lengths


def load_gene_to_hg_map(cluster_dir: Path, family_info: dict | None = None) -> tuple[dict, dict]:
    """Return ({gene_id: hg_id}, {hg_id: size}) from cluster FASTA files."""
    gene_to_hg: dict = {}
    hg_sizes: dict = {}
    if not cluster_dir.is_dir():
        return gene_to_hg, hg_sizes
    for fasta_file in sorted(cluster_dir.glob("*.fasta")):
        hg_id = fasta_file.stem
        parts = hg_id.split(".", 2)
        if len(parts) < 3:
            continue
        pref, family, _ = parts
        if not _family_info_allows(pref, family, family_info):
            continue
        size = 0
        try:
            with open(fasta_file) as fh:
                for line in fh:
                    if not line.startswith(">"):
                        continue
                    gene_id = line[1:].strip().split()[0]
                    if not gene_id:
                        continue
                    gene_to_hg[gene_id] = hg_id
                    size += 1
        except OSError:
            continue
        hg_sizes[hg_id] = size
    return gene_to_hg, hg_sizes


def build_exact_domain_catalog(search_dir: Path, gene_lengths: dict, family_info: dict | None = None) -> dict:
    """Return a compact per-protein catalog from *.domains_ummerged.csv and *.domains.csv."""
    if not search_dir.is_dir():
        return {
            "genes": [],
            "lengths": [],
            "tracks": [],
            "families": [],
            "names": [],
            "pfams": [],
        }

    merged_ranges: dict = defaultdict(dict)
    for merged_file in sorted(search_dir.glob("*.domains.csv")):
        track = merged_file.name[: -len(".domains.csv")]
        parts = track.split(".", 1)
        if len(parts) < 2:
            continue
        pref, family = parts
        if not _family_info_allows(pref, family, family_info):
            continue
        try:
            with open(merged_file) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 3:
                        continue
                    gene_id = cols[0].strip()
                    if not gene_id:
                        continue
                    try:
                        start = int(float(cols[1]))
                        end = int(float(cols[2]))
                    except ValueError:
                        continue
                    merged_ranges[track][gene_id] = [start, end]
        except OSError:
            continue

    family_dict: list[str] = []
    family_ix: dict = {}
    name_dict: list[str] = []
    name_ix: dict = {}
    pfam_dict: list[str] = []
    pfam_ix: dict = {}

    def _intern(value: str, labels: list[str], index: dict) -> int:
        value = (value or "").strip()
        idx = index.get(value)
        if idx is None:
            idx = len(labels)
            labels.append(value)
            index[value] = idx
        return idx

    gene_tracks: dict = defaultdict(list)
    for ummerged_file in sorted(search_dir.glob("*.domains_ummerged.csv")):
        track = ummerged_file.name[: -len(".domains_ummerged.csv")]
        parts = track.split(".", 1)
        if len(parts) < 2:
            continue
        pref, family = parts
        if not _family_info_allows(pref, family, family_info):
            continue
        per_gene_hits: dict = defaultdict(list)
        try:
            with open(ummerged_file) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 6:
                        continue
                    gene_id = cols[0].strip()
                    if not gene_id:
                        continue
                    try:
                        start = int(float(cols[1]))
                        end = int(float(cols[2]))
                    except ValueError:
                        continue
                    hit_name = cols[3].strip()
                    pfam_id = cols[4].strip()
                    per_gene_hits[gene_id].append([
                        start,
                        end,
                        _intern(hit_name, name_dict, name_ix),
                        _intern(pfam_id, pfam_dict, pfam_ix),
                    ])
        except OSError:
            continue

        fam_idx = _intern(track, family_dict, family_ix)
        for gene_id, hits in per_gene_hits.items():
            hits.sort(key=lambda rec: (rec[0], rec[1], rec[2], rec[3]))
            merged = merged_ranges.get(track, {}).get(gene_id)
            if merged:
                range_start, range_end = merged
            else:
                range_start = min(h[0] for h in hits)
                range_end = max(h[1] for h in hits)
            gene_tracks[gene_id].append([fam_idx, range_start, range_end, hits])

    genes = sorted(gene_tracks.keys())
    tracks = []
    lengths = []
    for gene_id in genes:
        track_rows = gene_tracks[gene_id]
        track_rows.sort(key=lambda rec: (family_dict[rec[0]], rec[1], rec[2]))
        tracks.append(track_rows)
        lengths.append(int(gene_lengths.get(gene_id, 0) or 0))

    return {
        "genes": genes,
        "lengths": lengths,
        "tracks": tracks,
        "families": family_dict,
        "names": name_dict,
        "pfams": pfam_dict,
    }


def _median(values: list[float]) -> float:
    if not values:
        return 0.0
    vals = sorted(values)
    mid = len(vals) // 2
    if len(vals) % 2:
        return float(vals[mid])
    return (float(vals[mid - 1]) + float(vals[mid])) / 2.0


def _interval_overlap_fraction(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    overlap = min(a_end, b_end) - max(a_start, b_start)
    if overlap <= 0:
        return 0.0
    a_len = max(1, a_end - a_start)
    b_len = max(1, b_end - b_start)
    return overlap / min(a_len, b_len)


def _normalize_architecture_hits(hits: list[dict]) -> list[dict]:
    """Collapse strongly overlapping same-name hits, keeping the best c-evalue."""
    if not hits:
        return []
    ordered = sorted(
        hits,
        key=lambda rec: (
            int(rec["start"]),
            int(rec["end"]),
            str(rec.get("name") or ""),
            float(rec.get("evalue", math.inf)),
        ),
    )
    kept: list[dict] = []
    for hit in ordered:
        if not kept:
            kept.append(dict(hit))
            continue
        prev = kept[-1]
        if (
            (prev.get("name") or "") == (hit.get("name") or "")
            and _interval_overlap_fraction(
                int(prev["start"]),
                int(prev["end"]),
                int(hit["start"]),
                int(hit["end"]),
            ) >= 0.5
        ):
            prev_ev = float(prev.get("evalue", math.inf))
            hit_ev = float(hit.get("evalue", math.inf))
            if hit_ev < prev_ev:
                kept[-1] = dict(hit)
        else:
            kept.append(dict(hit))
    return kept


def build_domain_architecture_catalog(search_dir: Path, gene_lengths: dict, gene_to_hg: dict, hg_sizes: dict, family_info: dict | None = None) -> dict:
    """Return a compact family-centric catalog of exact-hit domain architectures."""
    empty = {
        "families": [],
        "domains": [],
        "species": [],
        "hgs": [],
        "entries": [],
        "overview": {"nodes": [], "links": []},
    }
    if not search_dir.is_dir():
        return empty

    family_hits: dict = defaultdict(lambda: defaultdict(list))
    shared_domains: dict = defaultdict(set)
    for ummerged_file in sorted(search_dir.glob("*.domains_ummerged.csv")):
        family = ummerged_file.name[: -len(".domains_ummerged.csv")]
        parts = family.split(".", 1)
        if len(parts) < 2:
            continue
        pref, fam_name = parts
        if not _family_info_allows(pref, fam_name, family_info):
            continue
        try:
            with open(ummerged_file) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 6:
                        continue
                    gene_id = cols[0].strip()
                    hit_name = cols[3].strip()
                    if not gene_id or not hit_name:
                        continue
                    try:
                        start = int(float(cols[1]))
                        end = int(float(cols[2]))
                    except ValueError:
                        continue
                    try:
                        c_evalue = float(cols[5])
                    except ValueError:
                        c_evalue = math.inf
                    family_hits[family][gene_id].append(
                        {
                            "start": start,
                            "end": end,
                            "name": hit_name,
                            "evalue": c_evalue,
                        }
                    )
                    shared_domains[hit_name].add(family)
        except OSError:
            continue

    families = sorted(family_hits.keys())
    if not families:
        return empty

    domain_labels: list[str] = []
    domain_index: dict[str, int] = {}
    species_labels: list[str] = []
    species_index: dict[str, int] = {}
    hg_labels: list[str] = []
    hg_index: dict[str, int] = {}

    def _intern(value: str, labels: list[str], index: dict[str, int]) -> int:
        value = (value or "").strip()
        idx = index.get(value)
        if idx is None:
            idx = len(labels)
            labels.append(value)
            index[value] = idx
        return idx

    entries = []
    overview_node_counts: dict[str, int] = defaultdict(int)
    overview_family_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    overview_edge_counts: dict[tuple[str, str], int] = defaultdict(int)
    repeated_domain_families = sorted(
        d for d, fams in shared_domains.items() if len(fams) > 1
    )
    if repeated_domain_families:
        print(
            f"WARN: {len(repeated_domain_families)} PFAM names occur in multiple families; architecture calling remains family-scoped.",
            file=sys.stderr,
        )

    for family in families:
        arch_map: dict = {}
        for gene_id, hits in family_hits[family].items():
            norm_hits = _normalize_architecture_hits(hits)
            if not norm_hits:
                continue
            uniq_domains = []
            seen_domains = set()
            for hit in norm_hits:
                domain_name = hit["name"]
                if domain_name in seen_domains:
                    continue
                seen_domains.add(domain_name)
                uniq_domains.append(domain_name)
            for domain_name in uniq_domains:
                overview_node_counts[domain_name] += 1
                overview_family_counts[domain_name][family] += 1
            for i in range(len(uniq_domains)):
                for j in range(i + 1, len(uniq_domains)):
                    a, b = sorted((uniq_domains[i], uniq_domains[j]))
                    overview_edge_counts[(a, b)] += 1
            series = tuple(hit["name"] for hit in norm_hits)
            rec = arch_map.get(series)
            if rec is None:
                rec = {"members": [], "layouts": []}
                arch_map[series] = rec
            species = get_species_prefix(gene_id)
            protein_length = int(gene_lengths.get(gene_id, 0) or 0)
            max_end = max(int(hit["end"]) for hit in norm_hits)
            denom = protein_length if protein_length > 0 else max_end
            denom = max(denom, 1)
            rec["members"].append((species, gene_id))
            rec["layouts"].append(
                [
                    (
                        int(round((int(hit["start"]) / denom) * 1000)),
                        int(round((int(hit["end"]) / denom) * 1000)),
                        hit["name"],
                    )
                    for hit in norm_hits
                ]
            )

        arch_rows = []
        for series, rec in arch_map.items():
            members = sorted(rec["members"], key=lambda row: (row[0], row[1]))
            count = len(members)
            species_ids = sorted(
                {
                    _intern(species, species_labels, species_index)
                    for species, _gene in members
                }
            )
            hg_counts: dict[str, int] = defaultdict(int)
            for _species, gene_id in members:
                hg_id = gene_to_hg.get(gene_id)
                if hg_id:
                    hg_counts[hg_id] += 1
            hg_rows = []
            for hg_id, arch_count in sorted(hg_counts.items(), key=lambda row: (-row[1], row[0])):
                hg_rows.append([
                    _intern(hg_id, hg_labels, hg_index),
                    arch_count,
                    int(hg_sizes.get(hg_id, arch_count) or arch_count),
                ])
            consensus_hits = []
            if rec["layouts"]:
                for pos, domain_name in enumerate(series):
                    starts = [layout[pos][0] for layout in rec["layouts"] if len(layout) > pos]
                    ends = [layout[pos][1] for layout in rec["layouts"] if len(layout) > pos]
                    consensus_hits.append(
                        [
                            int(round(_median(starts))),
                            int(round(_median(ends))),
                            _intern(domain_name, domain_labels, domain_index),
                        ]
                    )
            arch_rows.append(
                [
                    [_intern(name, domain_labels, domain_index) for name in series],
                    count,
                    species_ids,
                    hg_rows,
                    consensus_hits,
                ]
            )
        arch_rows.sort(
            key=lambda rec: (
                -int(rec[1]),
                len(rec[0]),
                [domain_labels[idx] for idx in rec[0]],
            )
        )
        entries.append(arch_rows)

    family_lookup = {family: idx for idx, family in enumerate(families)}
    overview_nodes = []
    for domain_name, count in overview_node_counts.items():
        fam_rows = sorted(
            overview_family_counts[domain_name].items(),
            key=lambda row: (-row[1], row[0]),
        )
        dominant_family_idx = family_lookup.get(fam_rows[0][0], -1) if fam_rows else -1
        overview_nodes.append(
            [
                _intern(domain_name, domain_labels, domain_index),
                count,
                dominant_family_idx,
                [[family_lookup[fam], fam_count] for fam, fam_count in fam_rows if fam in family_lookup],
            ]
        )
    overview_nodes.sort(
        key=lambda rec: (
            -int(rec[1]),
            domain_labels[rec[0]],
        )
    )
    overview_links = []
    for (a, b), count in sorted(
        overview_edge_counts.items(),
        key=lambda row: (-row[1], row[0][0], row[0][1]),
    ):
        overview_links.append(
            [
                _intern(a, domain_labels, domain_index),
                _intern(b, domain_labels, domain_index),
                count,
            ]
        )

    return {
        "families": families,
        "domains": domain_labels,
        "species": species_labels,
        "hgs": hg_labels,
        "entries": entries,
        "overview": {"nodes": overview_nodes, "links": overview_links},
    }


def load_species_images(img_dir: Optional[str]) -> dict:
    """Return {species_prefix: data_uri} for PNG images found in img_dir."""
    import base64
    if not img_dir:
        return {}
    p = Path(img_dir)
    if not p.is_dir():
        return {}
    result = {}
    for png in sorted(p.glob("*.png")):
        try:
            data = base64.b64encode(png.read_bytes()).decode("ascii")
            result[png.stem] = f"data:image/png;base64,{data}"
        except OSError:
            continue
    return result


def load_species_info(path: Optional[str]) -> tuple:
    """Return ({prefix: full_name}, {prefix: group}) from a TSV with 2–3 columns."""
    if not path:
        return {}, {}
    try:
        names, groups = {}, {}
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0].strip():
                    prefix = parts[0].strip()
                    names[prefix] = parts[1].strip()
                    if len(parts) >= 3 and parts[2].strip():
                        groups[prefix] = parts[2].strip()
        return names, groups
    except (FileNotFoundError, OSError):
        return {}, {}


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
        if not _family_info_allows(pref, family, family_info):
            continue
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
        if not _family_info_allows(pref, family, family_info):
            continue
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
        # With ete3 format=1, numeric internal-node labels land in node.name rather
        # than node.support (support stays at the default 1.0).  Try node.name as a
        # fallback so that IQ-TREE bootstrap and GeneRax PP values are preserved.
        if "support" not in d and node.name:
            try:
                sv2 = float(node.name)
                if sv2 != 1.0:
                    d["support"] = round(sv2, 2)
                    d["name"] = ""  # was a numeric label, not an OG/clade name
            except (TypeError, ValueError):
                pass  # genuine string label (OG name, clade name) — leave as name
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


def load_gene_lengths(cluster_dir: Path, family_info: dict | None = None) -> dict:
    """Return {gene_id: protein_length} from cluster FASTAs."""
    lengths: dict = {}
    if not cluster_dir.is_dir():
        return lengths
    for fasta_file in sorted(cluster_dir.glob("*.fasta")):
        parts = fasta_file.stem.split(".", 2)
        if len(parts) < 3:
            continue
        pref, family, _ = parts
        if not _family_info_allows(pref, family, family_info):
            continue
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


def load_possvm_trees(possvm_dir: Path, source: str = "generax", family_info: dict | None = None, detect_method: bool = False) -> tuple[list, list, dict]:
    """Return (tree_records, all_species, gene_meta).  source='generax' or 'prev'.

    When detect_method is True, each tree's source is inferred from its POSSVM
    output filename (``.generax`` -> 'generax', ``.treefile`` -> 'iqtree') instead
    of the passed ``source``, so raw trees that land in results/possvm/ when the
    pipeline ran without --run_generax are not mislabelled as GeneRax."""
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
        if not _family_info_allows(prefix, family, family_info):
            continue

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

        if detect_method:
            nm = nwk.name
            rec_source = "generax" if ".generax" in nm else ("iqtree" if ".treefile" in nm else source)
        else:
            rec_source = source

        records.append({
            "id":               stem,
            "hg":               hg,
            "family":           family,
            "prefix":           prefix,
            "source":           rec_source,
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


_ASSET_DIR = Path(__file__).resolve().parent / "report_assets"


def _load_template() -> str:
    """Assemble the self-contained HTML template from report_assets/."""
    html = (_ASSET_DIR / "report.html").read_text()
    css  = (_ASSET_DIR / "report.css").read_text()
    js   = (_ASSET_DIR / "report.js").read_text()
    return html.replace("%%CSS%%", css).replace("%%JS%%", js)




# ── Argument parsing and main ─────────────────────────────────────────────────

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate interactive HTML report for step2 outputs.")
    p.add_argument("--results_dir", default=None,
                   help="Base results directory containing possvm/, possvm_prev/, search/, and clusters/")
    p.add_argument("--possvm_dir", default=None,
                   help="Directory with POSSVM *.ortholog_groups.newick files")
    p.add_argument("--possvm_prev_dir", default=None,
                   help="Directory with POSSVM results on original IQ-TREE2 trees (pre-GeneRax); "
                        "enables toggle between GeneRax and original trees in the report")
    p.add_argument("--search_dir", default=None,
                   help="Directory with *.genes.list files (for heatmap)")
    p.add_argument("--cluster_dir", default=None,
                   help="Directory with *.fasta HG files (for heatmap)")
    p.add_argument("--family_info", default="genefam.csv",
                   help="TSV with family metadata (genefam.csv or gene_families_searchinfo.csv)")
    p.add_argument("--species_tree", default=None,
                   help="Newick species tree (named internal nodes for clade colouring)")
    p.add_argument("--species_images", default="img/phylo",
                   help="Directory of per-species PNG images named {prefix}.png")
    p.add_argument("--species_info", default="data/species_info.tsv",
                   help="Two-column TSV mapping species prefix to full species name")
    p.add_argument("--refnames", default=None,
                   help="POSSVM refnames TSV (gene_id<TAB>name) for OG naming and tooltip annotation")
    p.add_argument("--refsps", default=None,
                   help="Comma-separated reference species filter, matching POSSVM --refsps")
    p.add_argument("--align_dir", default=None,
                   help="Directory with *.aln.fasta alignment files (enables the Alignments tab only when explicitly provided)")
    p.add_argument("--output", required=True, help="Output HTML file path")
    return p.parse_args(argv)


def _resolve_results_paths(args):
    results_dir = Path(args.results_dir) if args.results_dir else None

    def pick(explicit_value, default_suffix):
        if explicit_value:
            return Path(explicit_value)
        if results_dir:
            return results_dir / default_suffix
        return Path("results") / default_suffix

    return {
        "possvm_dir": pick(args.possvm_dir, "possvm"),
        "possvm_prev_dir": Path(args.possvm_prev_dir) if args.possvm_prev_dir else (
            (results_dir / "possvm_prev") if results_dir else None
        ),
        "search_dir": pick(args.search_dir, "search"),
        "cluster_dir": pick(args.cluster_dir, "clusters"),
        "align_dir": Path(args.align_dir) if args.align_dir else None,
    }


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


def _load_tree_records(possvm_dir: Path, possvm_prev_dir, family_info: dict | None = None):
    if not possvm_dir.exists():
        print(f"WARN: {possvm_dir} does not exist – no gene trees.", file=sys.stderr)

    raw_records, all_species, gene_meta = (
        load_possvm_trees(possvm_dir, source="generax", family_info=family_info, detect_method=True) if possvm_dir.exists() else ([], [], {})
    )
    generax_gene_meta = dict(gene_meta)
    # Deduplicate: same possvm dir may have both .generax.tree and .treefile outputs.
    # Keep first occurrence per id (generax files sort before treefile).
    seen_ids: set = set()
    records = []
    for r in raw_records:
        if r["id"] not in seen_ids:
            seen_ids.add(r["id"])
            records.append(r)
    print(f"Loaded {len(records)} gene trees, {len(all_species)} species.", file=sys.stderr)

    generax_ids = {r["id"] for r in records}
    prev_records: dict = {}
    if possvm_prev_dir:
        prev_dir = Path(possvm_prev_dir)
        if prev_dir.is_dir():
            prev_list, prev_sp, prev_gene_meta = load_possvm_trees(prev_dir, source="prev", family_info=family_info)
            prev_records = {r["id"]: r for r in prev_list}
            for gene_id, meta in prev_gene_meta.items():
                gene_meta.setdefault(gene_id, {}).update(meta)
            all_species = sorted(set(all_species) | set(prev_sp))
            records.extend(r for r in prev_list if r["id"] not in generax_ids)
            print(
                f"Loaded {len(prev_records)} prev gene trees (original IQ-TREE2).",
                file=sys.stderr,
            )
        else:
            prev_gene_meta = {}
    else:
        prev_gene_meta = {}

    return records, all_species, prev_records, gene_meta, generax_gene_meta, prev_gene_meta


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


def _build_aln_scripts(align_dir, records):
    import gzip as _gzip
    import base64 as _base64
    if not align_dir:
        return ""
    parts = []
    seen_ids = set()
    adir = Path(align_dir) if align_dir else None
    for rec in records:
        if rec["id"] in seen_ids:
            continue
        seen_ids.add(rec["id"])
        tag_id = _html.escape(rec["id"], quote=True)
        fasta_path = adir / f"{rec['id']}.aln.fasta" if adir else None
        if fasta_path and fasta_path.exists():
            raw = fasta_path.read_bytes()
            gz  = _gzip.compress(raw, compresslevel=9)
            b64 = _base64.b64encode(gz).decode()
            payload = f'{{"gz":"{b64}"}}'
        else:
            payload = '{"gz":null}'
        parts.append(
            f'<script type="application/json" id="alndata-{tag_id}">'
            + payload
            + "</script>"
        )
    return "\n".join(parts)


def _build_protein_domain_script(catalog):
    import gzip as _gzip
    import base64 as _base64

    if not catalog.get("genes"):
        payload = '{"gz":null}'
    else:
        raw = json.dumps(catalog, separators=(",", ":")).encode("utf-8")
        gz = _gzip.compress(raw, compresslevel=9)
        b64 = _base64.b64encode(gz).decode("ascii")
        payload = f'{{"gz":"{b64}"}}'
    return (
        '<script type="application/json" id="protein-domain-data">'
        + payload
        + "</script>"
    )


def _build_architecture_script(catalog):
    import gzip as _gzip
    import base64 as _base64

    if not catalog.get("families"):
        payload = '{"gz":null}'
    else:
        raw = json.dumps(catalog, separators=(",", ":")).encode("utf-8")
        gz = _gzip.compress(raw, compresslevel=9)
        b64 = _base64.b64encode(gz).decode("ascii")
        payload = f'{{"gz":"{b64}"}}'
    return (
        '<script type="application/json" id="architecture-data">'
        + payload
        + "</script>"
    )


def _build_lazy_scripts(records, prev_records, generax_gene_meta, prev_gene_meta, domain_spans):
    import gzip as _gzip
    import base64 as _base64
    lazy_parts = []
    seen_ids = set()
    for rec in records:
        if rec["id"] in seen_ids:
            continue
        seen_ids.add(rec["id"])
        detail = {"tree": rec["tree_dict"], "ogs": rec["ogs"]}
        rec_genes = {g for genes in rec["ogs"].values() for g in genes}
        if rec_genes:
            detail["meta"] = {g: generax_gene_meta[g] for g in rec_genes if g in generax_gene_meta}
        prev = prev_records.get(rec["id"])
        all_genes = set(rec_genes)
        if prev:
            detail["prev_tree"] = prev["tree_dict"]
            detail["prev_ogs"] = prev["ogs"]
            prev_genes = {g for genes in prev["ogs"].values() for g in genes}
            all_genes |= prev_genes
            if prev_genes:
                detail["prev_meta"] = {g: prev_gene_meta[g] for g in prev_genes if g in prev_gene_meta}
        if all_genes:
            detail["range"] = {g: domain_spans[g] for g in all_genes if g in domain_spans}
        tag_id = _html.escape(rec["id"], quote=True)
        raw = json.dumps(detail, separators=(",", ":")).encode("utf-8")
        gz = _gzip.compress(raw, compresslevel=9)
        b64 = _base64.b64encode(gz).decode("ascii")
        lazy_parts.append(
            f'<script type="application/json" id="treedata-{tag_id}">'
            + f'{{"gz":"{b64}"}}'
            + "</script>"
        )
    return "\n".join(lazy_parts)


def build_report_context(args) -> dict:
    resolved_paths = _resolve_results_paths(args)
    possvm_dir = resolved_paths["possvm_dir"]
    possvm_prev_dir = resolved_paths["possvm_prev_dir"]
    search_dir = resolved_paths["search_dir"]
    cluster_dir = resolved_paths["cluster_dir"]

    family_info = load_family_info(args.family_info)
    family_details = load_family_details(args.family_info)
    family_records = build_family_records(search_dir, family_info)
    hg_records = build_hg_records(cluster_dir, family_info)

    species_order, tree_dict, newick_raw, clade_groupings = _load_species_tree_bundle(
        args.species_tree
    )

    records, all_species, prev_records, gene_meta, generax_gene_meta, prev_gene_meta = _load_tree_records(
        possvm_dir, possvm_prev_dir, family_info
    )
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
    domain_hits = load_domain_hits(search_dir, family_info)
    domain_spans = build_domain_spans(domain_hits)
    gene_lengths = load_gene_lengths(cluster_dir, family_info)
    gene_to_hg, hg_sizes = load_gene_to_hg_map(cluster_dir, family_info)
    search_gene_lengths = load_search_gene_lengths(search_dir, family_info)
    protein_lengths = dict(gene_lengths)
    protein_lengths.update(search_gene_lengths)
    species_info, species_groups = load_species_info(args.species_info)
    species_images = load_species_images(args.species_images)
    refname_map = load_reference_names(args.refnames, args.refsps)
    for gene_id, length in gene_lengths.items():
        gene_meta.setdefault(gene_id, {})["length"] = length
    for gene_id, length in search_gene_lengths.items():
        gene_meta.setdefault(gene_id, {})["length"] = length
    for gene_id, ref_name in refname_map.items():
        gene_meta.setdefault(gene_id, {})["is_reference_gene"] = True
        gene_meta[gene_id]["reference_gene_name"] = ref_name
    protein_domain_catalog = build_exact_domain_catalog(search_dir, protein_lengths, family_info)
    architecture_catalog = build_domain_architecture_catalog(
        search_dir, protein_lengths, gene_to_hg, hg_sizes, family_info
    )

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
        "domain_spans": domain_spans,
        "gene_meta": gene_meta,
        "refname_map": refname_map,
        "species_info": species_info,
        "species_groups": species_groups,
        "species_images": species_images,
        "records": records,
        "prev_records": prev_records,
        "generax_gene_meta": generax_gene_meta,
        "prev_gene_meta": prev_gene_meta,
        "protein_domain_catalog": protein_domain_catalog,
        "architecture_catalog": architecture_catalog,
        "align_dir": resolved_paths["align_dir"],
        "have_alignments": bool(args.align_dir),
    }


def render_report_html(context: dict) -> str:
    return (
        _load_template()
        .replace("%%LAZY_SCRIPTS%%", _build_lazy_scripts(context["records"], context["prev_records"], context["generax_gene_meta"], context["prev_gene_meta"], context["domain_spans"]))
        .replace("%%ALN_SCRIPTS%%", _build_aln_scripts(context["align_dir"], context["records"]))
        .replace("%%PROTEIN_DOMAIN_SCRIPT%%", _build_protein_domain_script(context["protein_domain_catalog"]))
        .replace("%%ARCHITECTURE_SCRIPT%%", _build_architecture_script(context["architecture_catalog"]))
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
        .replace("%%SPECIES_INFO_JSON%%", json.dumps(context["species_info"]))
        .replace("%%SPECIES_GROUPS_JSON%%", json.dumps(context["species_groups"]))
        .replace("%%SPECIES_IMAGES_JSON%%", json.dumps(context["species_images"]))
        .replace("%%HAVE_ALIGNMENTS_JSON%%", json.dumps(context["have_alignments"]))
        .replace("%%HAVE_PROTEIN_DOMAINS_JSON%%", json.dumps(bool(context["protein_domain_catalog"].get("genes"))))
        .replace("%%HAVE_ARCHITECTURES_JSON%%", json.dumps(bool(context["architecture_catalog"].get("families"))))
    )


def main(argv=None):
    args = parse_args(argv)
    html = render_report_html(build_report_context(args))
    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
