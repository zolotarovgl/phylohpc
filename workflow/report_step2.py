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


def load_search_gene_lengths(search_dir: Path) -> dict:
    """Return {gene_id: protein_length} from *.seqs.fasta.fai sidecar files."""
    lengths: dict = {}
    if not search_dir.is_dir():
        return lengths
    for fai_file in sorted(search_dir.glob("*.seqs.fasta.fai")):
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


def load_gene_to_hg_map(cluster_dir: Path) -> tuple[dict, dict]:
    """Return ({gene_id: hg_id}, {hg_id: size}) from cluster FASTA files."""
    gene_to_hg: dict = {}
    hg_sizes: dict = {}
    if not cluster_dir.is_dir():
        return gene_to_hg, hg_sizes
    for fasta_file in sorted(cluster_dir.glob("*.fasta")):
        hg_id = fasta_file.stem
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


def build_exact_domain_catalog(search_dir: Path, gene_lengths: dict) -> dict:
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


def build_domain_architecture_catalog(search_dir: Path, gene_lengths: dict, gene_to_hg: dict, hg_sizes: dict) -> dict:
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
<script src="https://cdn.jsdelivr.net/npm/pako@2.1.0/dist/pako.min.js"></script>
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

/* ── Alignment pane ── */
#pane-align{flex-direction:row}
#aln-sidebar{flex:0 0 220px;display:flex;flex-direction:column;border-right:1px solid #ccc;background:#fff;overflow:hidden}
#aln-sidebar-top{padding:8px;border-bottom:1px solid #eee;flex-shrink:0}
#aln-search{width:100%;padding:5px 7px;font-size:11px;border:1px solid #ccc;border-radius:3px;box-sizing:border-box}
#aln-list{flex:1;overflow-y:auto}
#aln-main{flex:1;display:flex;flex-direction:column;overflow:hidden;min-width:0}
#pane-proteins{flex-direction:row}
#prot-sidebar{flex:0 0 250px;display:flex;flex-direction:column;border-right:1px solid #ccc;background:#fff;overflow:hidden}
#prot-sidebar-top{padding:8px;border-bottom:1px solid #eee;display:flex;flex-direction:column;gap:6px}
#prot-search{width:100%;padding:5px 7px;font-size:11px;border:1px solid #ccc;border-radius:3px;box-sizing:border-box}
#prot-count{font-size:11px;color:#7f8c8d}
#prot-list{flex:1;overflow-y:auto}
.prot-item{padding:7px 10px;border-bottom:1px solid #f1f3f5;cursor:pointer;font-size:11px;color:#2c3e50;display:flex;align-items:center;justify-content:space-between;gap:8px}
.prot-item:hover{background:#f8fbff}
.prot-item.active{background:#eaf3ff;color:#144a75;font-weight:600}
.prot-badge{font-size:10px;color:#7f8c8d;background:#f3f5f7;border:1px solid #d8dde3;border-radius:999px;padding:1px 6px;white-space:nowrap}
#prot-main{flex:1;display:flex;flex-direction:column;overflow:hidden;min-width:0;background:#fcfcfd}
#prot-empty{flex:1;color:#999;font-size:13px;padding:40px;text-align:center;display:flex;align-items:center;justify-content:center;flex-direction:column;gap:8px}
#prot-view{flex:1;overflow:auto;padding:18px 20px 22px}
#prot-summary{display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-bottom:14px}
#prot-legend{display:flex;flex-wrap:wrap;gap:10px;margin:-2px 0 14px}
.prot-legend-item{display:flex;align-items:center;gap:7px;font-size:11px;color:#455563;background:#f8fafc;border:1px solid #dce4eb;border-radius:999px;padding:4px 10px}
.prot-legend-swatch{display:inline-block;flex:0 0 auto}
.prot-legend-swatch.phylo{width:22px;height:12px;background:#cfcfcf;border:1px solid #111}
.prot-legend-swatch.hitspan{width:22px;height:10px;border:1px solid #111;background:#e2e2e2}
.prot-legend-swatch.hit{width:22px;height:12px;background:#72d25c;border:1px solid rgba(44,62,80,.6)}
.prot-legend-swatch.guide{width:22px;height:0;border-top:2px dashed #8fa4b9}
.prot-summary-chip{font-size:11px;color:#2c3e50;background:#f4f7fb;border:1px solid #dbe4ef;border-radius:999px;padding:4px 10px}
.prot-track{background:#fff;border:1px solid #e3e8ee;border-radius:12px;padding:10px 12px 12px;margin-bottom:12px;box-shadow:0 1px 2px rgba(0,0,0,.03)}
.prot-track-head{display:flex;align-items:center;justify-content:space-between;gap:12px;flex-wrap:wrap;margin-bottom:8px}
.prot-track-title{font-size:12px;font-weight:700;color:#234}
.prot-track-meta{font-size:11px;color:#708090}
.prot-svg-wrap{overflow-x:auto;padding-bottom:2px}
.prot-svg-wrap.compact{overflow:visible}
.prot-track-svg{display:block;min-width:920px}
.prot-track-svg.compact{min-width:0}
.prot-domain-legend{display:flex;flex-wrap:wrap;gap:6px;margin-top:8px}
.prot-domain-chip{display:inline-flex;align-items:center;gap:5px;font-size:10px;color:#304354;background:#f8fafc;border:1px solid #d7dde5;border-radius:999px;padding:2px 8px;white-space:nowrap}
.prot-domain-chip-swatch{display:inline-block;width:12px;height:12px;border-radius:3px;border:1px solid rgba(17,17,17,.55);flex:0 0 auto}
.prot-track-svg .prot-hit-group{transition:opacity .12s ease}
.prot-track-svg:hover .prot-hit-group{opacity:.18}
.prot-track-svg:hover .prot-hit-group:hover{opacity:1}
.prot-track-svg .prot-hit-guide{stroke:#8fa4b9;stroke-width:1;stroke-dasharray:4 3;opacity:.7}
.prot-track-svg .prot-span-guide{stroke:#5f6b76;stroke-width:1;stroke-dasharray:5 4;opacity:.8}
.prot-track-svg .prot-span-group{transition:opacity .12s ease}
.prot-track-svg .prot-hit-box{stroke:#111;stroke-width:1}
.prot-track-svg .prot-hit-span{fill:#e2e2e2;stroke:#111;stroke-width:1.2}
.prot-track-svg .prot-phylo-span{fill:#cfcfcf;stroke:#111;stroke-width:1.2}
.prot-track-svg .prot-hit-group:hover .prot-hit-box{stroke:#203040;stroke-width:1.2}
.prot-track-svg .prot-hit-coord{opacity:0;transition:opacity .12s ease;fill:#405162;font-size:10px;font-weight:700}
.prot-track-svg .prot-hit-group:hover .prot-hit-coord{opacity:1}
.prot-track-svg .prot-span-coord{opacity:0;transition:opacity .12s ease;fill:#2f3c46;font-size:10px;font-weight:700}
.prot-track-svg .prot-span-group:hover .prot-span-coord{opacity:1}
.prot-track-svg .prot-hit-hover-label{opacity:0;transition:opacity .12s ease;fill:#102030;font-size:10px;font-weight:700}
.prot-track-svg .prot-hit-group:hover .prot-hit-hover-label{opacity:1}
.prot-track-svg .prot-ruler-tick{stroke:#8a98a8;stroke-width:1}
.prot-track-svg .prot-ruler-tick.major{stroke:#516173;stroke-width:1.2}
.prot-track-svg .prot-ruler-label{fill:#667788;font-size:10px}
.prot-track-svg .prot-row-label{fill:#6f7f8d;font-size:10px;font-weight:700;letter-spacing:.04em;text-transform:uppercase}
.prot-hits{display:flex;flex-wrap:wrap;gap:6px;margin-top:8px}
.prot-hit-chip{font-size:10px;color:#34495e;background:#f8fafc;border:1px solid #d7dde5;border-radius:999px;padding:2px 8px;white-space:nowrap}
#pane-architectures{flex-direction:row}
#arch-fam-sidebar,#arch-list-sidebar{display:flex;flex-direction:column;background:#fff;overflow:hidden}
#arch-fam-sidebar{flex:0 0 240px;border-right:1px solid #ccc}
#arch-list-sidebar{flex:0 0 280px;border-right:1px solid #ccc}
#arch-fam-top,#arch-list-top{padding:8px;border-bottom:1px solid #eee;display:flex;flex-direction:column;gap:6px}
#arch-family-search{width:100%;padding:5px 7px;font-size:11px;border:1px solid #ccc;border-radius:3px;box-sizing:border-box}
#arch-family-count,#arch-arch-count{font-size:11px;color:#7f8c8d}
#arch-family-list,#arch-arch-list{flex:1;overflow-y:auto}
.arch-item{padding:7px 10px;border-bottom:1px solid #f1f3f5;cursor:pointer;font-size:11px;color:#2c3e50;display:flex;align-items:center;justify-content:space-between;gap:8px}
.arch-item:hover{background:#f8fbff}
.arch-item.active{background:#eaf3ff;color:#144a75;font-weight:600}
.arch-item-main{display:flex;flex-direction:column;gap:4px;min-width:0;flex:1 1 auto}
.arch-chip-list{display:flex;align-items:center;gap:0;flex-wrap:wrap;min-width:0}
.arch-domain-chip{display:inline-flex;align-items:center;max-width:100%;padding:2px 8px;border:1px solid rgba(17,17,17,.5);border-radius:999px;font-size:10px;font-weight:700;line-height:1.2;color:#102030;white-space:nowrap}
.arch-domain-chip.empty{background:#eef2f6;color:#5d6b78;border-color:#c5d0da}
.arch-domain-link{width:16px;height:2px;background:#9aa5b1;flex:0 0 auto;margin:0 2px}
.arch-domain-stack{position:relative;display:inline-block;vertical-align:middle;margin-right:2px}
.arch-domain-chip-copy{position:absolute;top:0;display:inline-flex;z-index:1;pointer-events:none}
.arch-domain-chip-copy .arch-domain-chip{color:transparent !important;text-shadow:none}
.arch-domain-chip-front{position:relative;z-index:3;display:inline-flex}
.arch-badge{font-size:10px;color:#7f8c8d;background:#f3f5f7;border:1px solid #d8dde3;border-radius:999px;padding:1px 6px;white-space:nowrap}
#arch-main{flex:1;display:flex;flex-direction:column;overflow:hidden;min-width:0;background:#fcfcfd}
#arch-empty{flex:1;color:#999;font-size:13px;padding:40px;text-align:center;display:flex;align-items:center;justify-content:center;flex-direction:column;gap:8px}
#arch-view{flex:1;overflow:auto;padding:18px 20px 22px}
#arch-summary{display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-bottom:14px}
.arch-summary-chip{font-size:11px;color:#2c3e50;background:#f4f7fb;border:1px solid #dbe4ef;border-radius:999px;padding:4px 10px}
#arch-network-wrap{background:#fff;border:1px solid #e3e8ee;border-radius:12px;padding:12px 14px;margin-bottom:14px;box-shadow:0 1px 2px rgba(0,0,0,.03)}
#arch-network{margin-top:10px}
.arch-net-svg{display:block;width:100%;min-height:360px}
.arch-net-link{stroke:#a9b3bc;stroke-opacity:.7}
.arch-net-node{stroke:#111;stroke-width:1.2;cursor:pointer}
.arch-net-node.active{stroke:#12385a;stroke-width:2.2}
.arch-net-label{fill:#30404f;font-size:11px;font-weight:700;pointer-events:none}
#arch-detail{display:grid;grid-template-columns:minmax(320px,1fr) minmax(360px,1fr);gap:14px;align-items:start}
#arch-svg-wrap{background:#fff;border:1px solid #e3e8ee;border-radius:12px;padding:12px 14px;box-shadow:0 1px 2px rgba(0,0,0,.03);overflow-x:auto;margin-bottom:0}
.arch-svg{display:block;min-width:920px}
.arch-domain-box{stroke:#111;stroke-width:1}
.arch-domain-label{fill:#102030;font-size:10px;font-weight:700}
.arch-domain-label.above{fill:#102030}
.arch-connector{stroke:#98a3ae;stroke-width:3;stroke-linecap:round}
.arch-row-label{fill:#6f7f8d;font-size:10px;font-weight:700;letter-spacing:.04em;text-transform:uppercase}
#arch-hgs-wrap{background:#fff;border:1px solid #e3e8ee;border-radius:12px;padding:12px 14px;box-shadow:0 1px 2px rgba(0,0,0,.03);display:flex;flex-direction:column;height:min(72vh,820px);min-height:420px}
#arch-hgs{display:flex;flex-direction:column;gap:6px;margin-top:10px;overflow:auto;flex:1}
.arch-hg-row{display:flex;align-items:center;justify-content:space-between;gap:10px;padding:7px 9px;border:1px solid #d7dde5;border-radius:8px;background:#f8fafc;cursor:pointer}
.arch-hg-row:hover{background:#edf5ff;border-color:#b8cde5}
.arch-hg-main{display:flex;flex-direction:column;gap:2px;min-width:0}
.arch-hg-title{font-size:11px;font-weight:700;color:#234}
.arch-hg-meta{font-size:10px;color:#6d7c8a}
.arch-hg-badges{display:flex;align-items:center;gap:6px;flex:0 0 auto}
.arch-share-badge{font-size:10px;color:#35506a;background:#eef3f8;border:1px solid #c9d4df;border-radius:999px;padding:1px 6px;white-space:nowrap;font-weight:700}
#arch-singletons-wrap{display:flex;align-items:center;gap:6px;font-size:11px;color:#506070}
#aln-controls{flex-shrink:0;padding:6px 12px;background:#f5f5f5;border-bottom:1px solid #ddd;display:flex;align-items:center;gap:10px;flex-wrap:wrap}
#aln-viewer{flex:1;display:none;grid-template-columns:220px 5px 1fr;grid-template-rows:24px 1fr;overflow:hidden;background:#fff;min-height:0}
#aln-corner{grid-column:1;grid-row:1;background:#f0f0f0;border-right:1px solid #ccc;border-bottom:1px solid #ccc;position:relative;overflow:hidden;z-index:4}
.aln-col-header{position:absolute;top:0;height:100%;display:flex;align-items:center;justify-content:flex-start;padding:0 6px;box-sizing:border-box;color:#44515e;font-size:10px;font-weight:700;letter-spacing:.02em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;border-right:1px solid #d9d9d9;background:rgba(255,255,255,.55);user-select:none}
.aln-col-header.drag{cursor:grab}
#aln-resize-bar{grid-column:2;grid-row:1/3;background:#e8e8e8;border-left:1px solid #ccc;border-right:1px solid #ccc;cursor:col-resize;z-index:2}
#aln-resize-bar:hover,#aln-resize-bar.dragging{background:#bbb}
#aln-tree-resize-bar,#aln-ref-resize-bar,#aln-species-resize-bar,#aln-mrca-resize-bar,#aln-range-resize-bar{position:absolute;top:0;bottom:0;width:5px;background:#e8e8e8;border-left:1px solid #ccc;border-right:1px solid #ccc;cursor:col-resize;z-index:3}
#aln-tree-resize-bar:hover,#aln-tree-resize-bar.dragging,#aln-ref-resize-bar:hover,#aln-ref-resize-bar.dragging,#aln-species-resize-bar:hover,#aln-species-resize-bar.dragging,#aln-mrca-resize-bar:hover,#aln-mrca-resize-bar.dragging,#aln-range-resize-bar:hover,#aln-range-resize-bar.dragging{background:#bbb}
#aln-ruler-wrap{grid-column:3;grid-row:1;overflow:hidden;border-bottom:1px solid #ccc;background:#f8f8f8;position:relative}
#aln-names-wrap{grid-column:1;grid-row:2;overflow-y:auto;overflow-x:hidden;border-right:1px solid #ccc;background:#fafafa;position:relative}
#aln-names-wrap::-webkit-scrollbar{width:10px;height:10px}
#aln-names-wrap::-webkit-scrollbar-thumb{background:#bbb;border-radius:5px}
#aln-seq-wrap{grid-column:3;grid-row:2;overflow:auto;position:relative}
#aln-seq-wrap::-webkit-scrollbar{width:10px;height:10px}
#aln-seq-wrap::-webkit-scrollbar-thumb{background:#bbb;border-radius:5px}
#aln-empty{flex:1;color:#999;font-size:13px;padding:40px;text-align:center;display:flex;align-items:center;justify-content:center;flex-direction:column;gap:8px}

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
#tooltip{position:fixed;display:none;pointer-events:none;background:rgba(255,255,255,.97);border:1px solid #bbb;border-radius:5px;padding:8px 10px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.15);z-index:100;max-width:min(920px,calc(100vw - 40px));max-height:min(85vh,900px);overflow:auto}
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

/* ── Status tab ── */
#pane-status{flex-direction:column;overflow:auto;padding:20px 24px;gap:20px;background:#f7f7f7}
.status-cards{display:flex;gap:14px;flex-wrap:wrap;flex-shrink:0}
.status-card{background:#fff;border:1px solid #dde3e8;border-radius:8px;padding:14px 20px;min-width:130px;display:flex;flex-direction:column;gap:4px;box-shadow:0 1px 3px rgba(0,0,0,.06)}
.status-card .sc-val{font-size:26px;font-weight:700;color:#1e3a5f;line-height:1.1}
.status-card .sc-lbl{font-size:11px;color:#7a8898}
.status-card .sc-sub{font-size:10px;color:#a0aab3}
.status-section{background:#fff;border:1px solid #dde3e8;border-radius:8px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.06)}
.status-section-hdr{padding:8px 14px;background:#f0f4f8;border-bottom:1px solid #dde3e8;font-size:12px;font-weight:700;color:#2c3e50}
.status-tbl{border-collapse:collapse;width:100%;font-size:11px}
.status-tbl th{padding:5px 10px;background:#f8fafc;border-bottom:2px solid #e4e9ef;text-align:left;color:#596674;font-weight:600;white-space:nowrap;cursor:pointer;user-select:none}
.status-tbl th:hover{background:#edf2f7}
.status-tbl td{padding:5px 10px;border-bottom:1px solid #f0f4f7;vertical-align:middle}
.status-tbl tr:last-child td{border-bottom:none}
.status-tbl tr:hover td{background:#f8fafc}
.prog-bar-wrap{display:flex;align-items:center;gap:7px}
.prog-bar-bg{flex:1;height:7px;background:#e8ecf0;border-radius:4px;overflow:hidden;min-width:60px}
.prog-bar-fill{height:7px;border-radius:4px;transition:width .2s}
.prog-lbl{font-size:10px;color:#596674;white-space:nowrap;min-width:40px;text-align:right}
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
    <button class="tab-btn" data-tab="architectures" onclick="switchTab('architectures')">&#9776; Domain Architectures</button>
    <button class="tab-btn active" data-tab="sptree" onclick="switchTab('sptree')">&#10022; Species Tree</button>
    <button class="tab-btn" data-tab="heatmap" onclick="switchTab('heatmap')">&#9639; Counts</button>
    <button class="tab-btn" data-tab="trees" onclick="switchTab('trees')">&#11044; Gene Trees</button>
    <button class="tab-btn" id="tab-btn-align" data-tab="align" onclick="switchTab('align')">&#9644; Alignments</button>
  </div>

  <!-- ── Heatmap pane ── -->
  <div class="tab-pane" id="pane-heatmap">
    <div id="hm-custom-bar" style="display:none;position:fixed;z-index:500;border:1px solid #c8b96e;border-top:none;background:#fffbf0;box-shadow:0 4px 12px rgba(0,0,0,.15);border-radius:0 0 6px 6px;flex-direction:column;gap:0">
      <div style="display:flex;align-items:center;gap:6px;padding:3px 10px;cursor:pointer" onclick="hmToggleCustomBar()">
        <span id="hm-custom-toggle" style="font-size:10px;color:#888">&#9654;</span>
        <span id="hm-custom-summary" style="font-size:11px;color:#888;font-weight:600">Custom selection:</span>
        <button onclick="event.stopPropagation();hmCustomClear()" style="margin-left:auto;padding:1px 8px;font-size:10px;border:1px solid #ccc;border-radius:3px;background:#fff;cursor:pointer">Clear all</button>
        <button onclick="event.stopPropagation();hmExitCustom()" style="padding:1px 8px;font-size:10px;border:1px solid #4a90d9;color:#4a90d9;border-radius:3px;background:#fff;cursor:pointer">&#8592; Browse</button>
      </div>
      <div id="hm-custom-chips-wrap" style="display:none;padding:2px 10px 4px;max-height:120px;overflow-y:auto">
        <span id="hm-custom-chips" style="display:flex;gap:4px;flex-wrap:wrap;align-items:center"></span>
      </div>
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
            <input type="range" id="hm-col-rot-slider" min="0" max="90" step="5" value="90" style="width:70px;cursor:pointer;accent-color:#4a90d9">
            <span id="hm-col-rot-val" style="width:24px;text-align:right">90</span>&deg;
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
          <button id="btn-hm-sp-logos" onclick="hmToggleSpLogos()" title="Toggle species images in cladogram" style="padding:2px 8px;font-size:11px;border:1px solid #888;color:#555;border-radius:3px;background:#fff;cursor:pointer">&#128444; Logos</button>
          <button id="btn-hm-group-hg" onclick="hmToggleGroupByHG()" title="Group columns by Homology Group (drag group headers to reorder)" style="padding:2px 8px;font-size:11px;border:1px solid #888;color:#555;border-radius:3px;background:#fff;cursor:pointer">&#8801; Group HG</button>
          <button id="btn-hm-flip" onclick="hmToggleFlip()" title="Flip heatmap: species as columns, HGs/OGs as rows with horizontal species tree" style="padding:2px 8px;font-size:11px;border:1px solid #888;color:#555;border-radius:3px;background:#fff;cursor:pointer">&#8646; Flip</button>
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
            <th class="fam-th" data-col="family"   style="text-align:left;padding:6px 8px;cursor:pointer;white-space:nowrap">Category / Family &#8597;</th>
            <th class="fam-th" data-col="pfam"     style="text-align:left;padding:6px 8px;cursor:pointer">PFAM Domain(s)</th>
            <th class="fam-th" data-col="category" style="text-align:left;padding:6px 8px;cursor:pointer;white-space:nowrap">Class &#8597;</th>
            <th class="fam-th" data-col="n_species"style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap" title="Number of species">Species &#8597;</th>
            <th style="text-align:right;padding:6px 8px;white-space:nowrap;color:#888;font-weight:600">Fam.</th>
            <th class="fam-th" data-col="n_hgs"    style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap">HGs &#8597;</th>
            <th class="fam-th" data-col="n_trees"  style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap" title="HGs with any gene tree">Trees &#8597;</th>
            <th class="fam-th fam-th-generax" data-col="n_generax" style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap" title="HGs with a GeneRax tree">GeneRax &#8597;</th>
            <th class="fam-th" data-col="total"    style="text-align:right;padding:6px 8px;cursor:pointer;white-space:nowrap">Genes &#8597;</th>
          </tr>
        </thead>
        <tbody id="fam-tbody"></tbody>
      </table>
    </div>
  </div>

  <!-- ── Alignment pane ── -->
  <div class="tab-pane" id="pane-align">
    <div id="aln-sidebar">
      <div id="aln-sidebar-top">
        <input id="aln-search" type="text" placeholder="Search HG / family…" oninput="renderAlnSidebar(this.value)">
      </div>
      <div id="aln-list"></div>
    </div>
    <div id="aln-main">
      <div id="aln-controls">
        <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:5px">Cell:
          <input type="range" id="aln-cell-slider" min="4" max="18" step="1" value="9" style="width:70px;cursor:pointer;accent-color:#4a90d9" oninput="alnCellW=+this.value;alnCellH=Math.round(alnCellW*1.6);document.getElementById('aln-cell-val').textContent=alnCellW;renderAlignment()">
          <span id="aln-cell-val">9</span>px
        </label>
        <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:5px">Domains:
          <input type="range" id="aln-domain-chip-slider" min="0.6" max="2.2" step="0.1" value="1.0" style="width:78px;cursor:pointer;accent-color:#48a868" oninput="alnSetDomainChipScale(this.value)">
          <span id="aln-domain-chip-val">1.0</span>x
        </label>
        <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:5px">Species:
          <input type="text" id="aln-seq-filter" placeholder="Enter to add…" style="font-size:11px;padding:3px 6px;border:1px solid #ccc;border-radius:3px;width:120px" onkeydown="_alnSpFilterKey(event)">
        </label>
        <span id="aln-sp-tags" style="display:flex;align-items:center;flex-wrap:wrap;gap:2px"></span>
        <button class="ctrl-btn" onclick="alnToggleBrLen()" id="btn-aln-brlen" title="Toggle branch lengths on/off">Br. len</button>
        <button class="ctrl-btn active-btn" onclick="alnToggleSeqPanel()" id="btn-aln-seqpanel" title="Show/hide sequence alignment pane">Alignment</button>
        <button class="ctrl-btn" onclick="alnToggleConsensus()" id="btn-aln-cons" title="Show/hide consensus row">Consensus</button>
        <button class="ctrl-btn" onclick="alnToggleCollapseAllOGs()" id="btn-aln-collapse-ogs" title="Collapse all orthogroups to single rows">Collapse OGs</button>
        <button class="ctrl-btn active-btn" onclick="alnToggleSeqId()" id="btn-aln-seqid" title="Show/hide sequence ID column">ID</button>
        <button class="ctrl-btn" onclick="alnToggleRefCol()" id="btn-aln-ref" title="Show/hide reference ortholog name column">Ref</button>
        <button class="ctrl-btn active-btn" onclick="alnToggleSpeciesCol()" id="btn-aln-species" title="Show/hide species column">Species</button>
        <button class="ctrl-btn active-btn" onclick="alnToggleMrcaCol()" id="btn-aln-mrca" title="Show/hide MRCA clade column">MRCA</button>
        <button class="ctrl-btn active-btn" onclick="alnToggleRangeCol()" id="btn-aln-range" title="Show/hide compact domain-architecture column">Domain Architecture</button>
        <button class="ctrl-btn" onclick="alnToggleSource()" id="btn-aln-src" title="Toggle GeneRax / IQ-TREE tree source" style="display:none">IQ-TREE</button>
        <button class="ctrl-btn" onclick="alnDownload()" title="Download alignment as FASTA">&#11015; FASTA</button>
        <button class="ctrl-btn" onclick="alnDownloadPng()" title="Download alignment as PNG">&#11015; PNG</button>
        <span id="aln-info" style="font-size:11px;color:#888;margin-left:8px"></span>
      </div>
      <div id="aln-viewer">
        <div id="aln-corner"></div>
        <div id="aln-resize-bar"></div>
        <div id="aln-ruler-wrap"><canvas id="aln-ruler"></canvas></div>
        <div id="aln-names-wrap" onscroll="_alnSyncFromNamesWrap(this)" onmousemove="_alnHandleNameHover(event)" onmouseleave="_alnHideTip()" onclick="_alnHandleNameClick(event)"><canvas id="aln-names-canvas"></canvas><div id="aln-tree-resize-bar" style="display:none"></div><div id="aln-ref-resize-bar" style="display:none"></div><div id="aln-species-resize-bar" style="display:none"></div><div id="aln-mrca-resize-bar" style="display:none"></div><div id="aln-range-resize-bar" style="display:none"></div></div>
        <div id="aln-seq-wrap" onscroll="syncAlnScroll(this)" onclick="_alnHandleSeqClick(event)"><canvas id="aln-seq-canvas"></canvas></div>
      </div>
      <div id="aln-empty">
        <div style="font-size:32px">&#9644;</div>
        <div>Select an HG to view its alignment</div>
      </div>
    </div>
  </div>

  <!-- ── Domain-architecture pane ── -->
  <div class="tab-pane" id="pane-architectures">
    <div id="arch-fam-sidebar">
      <div id="arch-fam-top">
        <input id="arch-family-search" type="text" placeholder="Search family…" oninput="renderArchitectureFamilySidebar(this.value)">
        <div id="arch-family-count"></div>
      </div>
      <div id="arch-family-list"></div>
    </div>
    <div id="arch-list-sidebar">
      <div id="arch-list-top">
        <div id="arch-singletons-wrap">
          <label style="display:flex;align-items:center;gap:6px;cursor:pointer"><input id="arch-show-singletons" type="checkbox" onchange="architectureShowSingletons=this.checked;renderArchitectureList();renderArchitectureDetail()"> Include singletons</label>
        </div>
        <div id="arch-arch-count"></div>
      </div>
      <div id="arch-arch-list"></div>
    </div>
    <div id="arch-main">
      <div id="arch-empty">
        <div style="font-size:34px">&#9776;</div>
        <div>Select a family to explore its PFAM co-occurrence network</div>
      </div>
      <div id="arch-view" style="display:none">
        <div id="arch-summary"></div>
        <div id="arch-network-wrap">
          <div style="font-size:12px;font-weight:700;color:#234">PFAM co-occurrence network</div>
          <div id="arch-network"></div>
        </div>
        <div id="arch-detail-empty" style="display:none;color:#71808f;font-size:12px;background:#fff;border:1px dashed #cfd8e3;border-radius:12px;padding:18px 16px;text-align:center">Click a PFAM node to list architectures, then select an architecture to inspect its consensus and HG distribution</div>
        <div id="arch-detail" style="display:none">
          <div id="arch-svg-wrap"></div>
          <div id="arch-hgs-wrap">
            <div style="font-size:12px;font-weight:700;color:#234">Homology groups carrying this architecture</div>
            <div id="arch-hgs"></div>
          </div>
        </div>
      </div>
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
      <label style="font-size:11px;color:#555;display:flex;align-items:center;gap:6px">
        Species colours:
        <select id="sp-palette-select" style="font-size:11px;padding:2px 4px;border:1px solid #bbb;border-radius:3px;background:#fff">
          <option value="turbo">Turbo</option>
          <option value="viridis">Viridis</option>
          <option value="warm">Warm</option>
          <option value="cool">Cool</option>
          <option value="tableau">Tableau</option>
          <option value="set3">Set3</option>
        </select>
      </label>
      <span id="sp-palette-preview" style="display:inline-flex;align-items:center;gap:2px;min-height:16px"></span>
      <button class="ctrl-btn" id="btn-sp-apply-palette" onclick="applySelectedSpeciesPalette()" title="Apply the selected multi-colour palette to all species">Apply palette</button>
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
<div id="hm-open-popup" style="position:fixed;display:none;background:#fff;border:1px solid #bbb;border-radius:6px;padding:8px 10px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.18);z-index:320;min-width:160px">
  <div style="font-size:10px;color:#888;margin-bottom:6px;font-weight:600" id="hm-open-popup-title">Open in</div>
  <div style="display:flex;flex-direction:column;gap:5px">
    <button id="hm-open-tree" style="padding:4px 10px;font-size:11px;border:1px solid #4a90d9;border-radius:4px;background:#f0f6ff;color:#2c6090;cursor:pointer;text-align:left">Show Tree</button>
    <button id="hm-open-align" style="padding:4px 10px;font-size:11px;border:1px solid #27ae60;border-radius:4px;background:#f0fff4;color:#1a7a42;cursor:pointer;text-align:left">Show Alignment</button>
  </div>
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
<div id="sptree-node-popup" style="position:fixed;display:none;background:#fff;border:1px solid #bbb;border-radius:6px;padding:6px 8px;font-size:11px;box-shadow:0 2px 8px rgba(0,0,0,.18);z-index:260;flex-direction:column;gap:5px">
  <div style="font-size:10px;color:#888;margin-bottom:2px;font-weight:600" id="snp-title"></div>
  <div style="display:flex;gap:5px;flex-wrap:wrap">
    <button id="snp-collapse" style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x25B7; Collapse</button>
    <button id="snp-flip"     style="padding:3px 10px;font-size:11px;border:1px solid #aaa;border-radius:4px;background:#f8f8f8;cursor:pointer">&#x21C5; Flip</button>
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

<!-- ── Per-HG alignment data ── -->
<div id="aln-data" style="display:none">
%%ALN_SCRIPTS%%
</div>

<!-- ── Protein-domain catalog ── -->
<div id="protein-domain-data-wrap" style="display:none">
%%PROTEIN_DOMAIN_SCRIPT%%
</div>

<!-- ── Domain-architecture catalog ── -->
<div id="architecture-data-wrap" style="display:none">
%%ARCHITECTURE_SCRIPT%%
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
const SPECIES_INFO   = %%SPECIES_INFO_JSON%%;   // {species_prefix: full_species_name}
const SPECIES_GROUPS = %%SPECIES_GROUPS_JSON%%; // {species_prefix: group_name} (optional 3rd col of species_info.tsv)
const SPECIES_IMAGES = %%SPECIES_IMAGES_JSON%%; // {species_prefix: data_uri}
const HAVE_ALIGNMENTS = %%HAVE_ALIGNMENTS_JSON%%;
const HAVE_PROTEIN_DOMAINS = %%HAVE_PROTEIN_DOMAINS_JSON%%;
const HAVE_ARCHITECTURES = %%HAVE_ARCHITECTURES_JSON%%;

// ═══════════════════════════════════════════════════════════════════════════════
// COLOUR SYSTEM
// ═══════════════════════════════════════════════════════════════════════════════
const palette = [
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
  "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
  "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5"
];
const SPECIES_COLOR_PRESETS = {
  turbo:   {label:"Turbo", type:"interp", fn:d3.interpolateTurbo},
  viridis: {label:"Viridis", type:"interp", fn:d3.interpolateViridis},
  warm:    {label:"Warm", type:"interp", fn:d3.interpolateWarm},
  cool:    {label:"Cool", type:"interp", fn:d3.interpolateCool},
  tableau: {label:"Tableau", type:"array", colors:d3.schemeTableau10},
  set3:    {label:"Set3", type:"array", colors:d3.schemeSet3}
};
const SP_COLORS = {};
function spColor(sp) { return SP_COLORS[sp] || "#aaa"; }

// ── Group colors (average RGB of member species colors) ───────────────────────
function _cssToRgb(css){
  const m=css.match(/^#([0-9a-f]{6})$/i);
  if(m) return [parseInt(m[1].slice(0,2),16),parseInt(m[1].slice(2,4),16),parseInt(m[1].slice(4,6),16)];
  const m2=css.match(/rgb\s*\(\s*(\d+)[,\s]+(\d+)[,\s]+(\d+)/i);
  if(m2) return [+m2[1],+m2[2],+m2[3]];
  return null;
}
function _rgbToHex(r,g,b){ return '#'+[r,g,b].map(v=>Math.round(v).toString(16).padStart(2,'0')).join(''); }
let groupColors={};
function recomputeGroupColors(){
  const members={};
  Object.entries(SPECIES_GROUPS).forEach(([sp,grp])=>{ (members[grp]||(members[grp]=[])).push(sp); });
  const out={};
  Object.entries(members).forEach(([grp,sps])=>{
    const rgbs=sps.map(sp=>_cssToRgb(spColor(sp))).filter(Boolean);
    if(!rgbs.length) return;
    out[grp]=_rgbToHex(d3.mean(rgbs,c=>c[0]),d3.mean(rgbs,c=>c[1]),d3.mean(rgbs,c=>c[2]));
  });
  groupColors=out;
}
function groupColor(sp){ const g=SPECIES_GROUPS[sp]; return (g&&groupColors[g])||"#aaa"; }

function orderedSpeciesPaletteList(){
  const seen=new Set(), out=[];
  SPECIES_ORDER.forEach(sp=>{ if(!seen.has(sp)){ seen.add(sp); out.push(sp); } });
  ALL_SPECIES.forEach(sp=>{ if(!seen.has(sp)){ seen.add(sp); out.push(sp); } });
  return out;
}

function _sampleSpeciesPreset(name, i, n){
  const preset=SPECIES_COLOR_PRESETS[name]||SPECIES_COLOR_PRESETS.turbo;
  if(preset.type==="interp") return preset.fn(n>1?i/(n-1):0.5);
  const cols=preset.colors||palette;
  if(!cols.length) return "#888";
  if(n<=cols.length) return cols[Math.round((cols.length-1)*(n>1?i/(n-1):0.5))];
  return d3.interpolateRgbBasis(cols)(n>1?i/(n-1):0.5);
}

function orderedSpeciesGroupBuckets(){
  const buckets=[], byKey=new Map();
  orderedSpeciesPaletteList().forEach(sp=>{
    const grp=(SPECIES_GROUPS[sp]||"").trim();
    const key=grp?`grp:${grp}`:`solo:${sp}`;
    if(!byKey.has(key)){
      const rec={key, group:grp, species:[]};
      byKey.set(key, rec);
      buckets.push(rec);
    }
    byKey.get(key).species.push(sp);
  });
  return buckets;
}

function _speciesVariantAroundBase(base, i, n){
  if(n<=1) return base;
  const col=d3.color(base);
  if(!col) return base;
  const lo=0.18, hi=0.82;
  const t=lo + (hi-lo)*(n>1?i/(n-1):0.5);
  const light=col.brighter(1.1).formatHex();
  const dark=col.darker(0.9).formatHex();
  return d3.interpolateLab(light, dark)(t);
}

function refreshSpeciesColorViews(){
  recomputeGroupColors();
  drawSpeciesTree();
  if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
  else drawCladogram();
  if(rootNode) renderTree(false);
  if(currentAlnId) renderAlignment();
  const mini=document.getElementById("mini-sp-panel");
  if(mini&&mini.style.display==="block") drawMiniSpTree();
}

function applySpeciesPalette(name, redraw=true){
  const buckets=orderedSpeciesGroupBuckets();
  buckets.forEach((bucket,i)=>{
    const base=_sampleSpeciesPreset(name, i, buckets.length);
    bucket.species.forEach((sp,j)=>{
      SP_COLORS[sp]=_speciesVariantAroundBase(base, j, bucket.species.length);
    });
  });
  ALL_SPECIES.forEach(sp=>{ if(!SP_COLORS[sp]) SP_COLORS[sp]="#aaa"; });
  if(redraw) refreshSpeciesColorViews();
}

function updateSpeciesPalettePreview(name){
  const wrap=document.getElementById("sp-palette-preview");
  if(!wrap) return;
  const n=6;
  wrap.innerHTML=Array.from({length:n},(_,i)=>{
    const c=_sampleSpeciesPreset(name,i,n);
    return `<span style="display:inline-block;width:12px;height:12px;border-radius:2px;background:${c};border:1px solid rgba(0,0,0,.15)"></span>`;
  }).join("");
}

function applySelectedSpeciesPalette(){
  const sel=document.getElementById("sp-palette-select");
  applySpeciesPalette(sel?sel.value:"turbo");
}

(function(){
  applySpeciesPalette("turbo", false);
  recomputeGroupColors();
})();

// ── Stable IDs for SP_TREE_DATA nodes (allows editing original by reference) ──
const spNodeById = new Map();   // _spId → SP_TREE_DATA node reference
(function tagSpNodes(n, i={v:0}){
  n._spId = String(i.v++);
  spNodeById.set(n._spId, n);
  (n.children||[]).forEach(c=>tagSpNodes(c,i));
})(SP_TREE_DATA || {});

function spTreeNodeKey(n){
  if(n&&n._spId!=null) return "sp:"+n._spId;
  return "leaves:"+spNodeLeaves(n).slice().sort().join(",");
}

function findSpMrcaForSpecies(root, speciesSet){
  let best = null;
  function walk(n){
    let hits = 0;
    if(n.children&&n.children.length){
      n.children.forEach(c=>{ hits += walk(c); });
    } else if(speciesSet.has(n.name)) {
      hits = 1;
    }
    if(!best&&hits===speciesSet.size) best = n;
    return hits;
  }
  walk(root);
  return best;
}

function treeToNewick(n){
  const name=(n.name||"").replace(/[(),:;]/g,"_");
  const dist=n.dist!=null?":"+n.dist:"";
  if(!n.children||!n.children.length) return name+dist;
  return "("+n.children.map(treeToNewick).join(",")+")"+(name||"")+dist;
}

let colorMode     = "species";   // "species" | "clade" | "og" | "group"
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
let hmColRotation = 90;      // column label rotation angle (degrees)
let hmColorMode   = "zscore"; // "zscore" | "absolute"
let hmTextFilter  = "";       // free-text filter for heatmap columns
const spCollapsed = new Set();   // stable species-tree node keys collapsed in species tree
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
  const baseColor = colorMode === "species" ? spColor(sp)
                  : colorMode === "group"   ? groupColor(sp)
                  : (cladeSp2Color[sp] || "#ccc");
  if (hlSet !== null) {
    if (!hlSet.has(sp)) return "#ccc";
    const gi = hlGroupIndex.get(sp);
    return gi !== undefined ? hlTagColor(gi) : baseColor;
  }
  return baseColor;
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
let hmGroupByHG    = false;   // group columns by parent HG / family / class
let hmGroupOrderOverride = null; // ordered array of group keys for drag-reorder of groups
let hmFlipped      = false;   // transpose: species become columns, HG/OGs become rows
let hmOGIndex      = null;   // lazy: {og_name → {hgId,family,cls,total,species_counts}}
let hmGeneIndex    = null;   // lazy: {gene_id → og_name}
let hmSearchIndex  = null;   // lazy: [{label,type,cls,fam,hg_id}] for global search
let hmBatchSelection = new Set(); // search dropdown items staged for bulk add
let hmCustomBarExpanded = false;   // whether the custom-selection chip row is visible
let hmShowSpLogos = false;        // show species images in cladogram

function hmToggleCustomBar(){
  hmCustomBarExpanded = !hmCustomBarExpanded;
  const wrap = document.getElementById("hm-custom-chips-wrap");
  const arrow = document.getElementById("hm-custom-toggle");
  if(wrap) wrap.style.display = hmCustomBarExpanded ? "block" : "none";
  if(arrow) arrow.innerHTML   = hmCustomBarExpanded ? "&#9660;" : "&#9654;";
}
function hmToggleSpLogos(){
  hmShowSpLogos = !hmShowSpLogos;
  const btn = document.getElementById("btn-hm-sp-logos");
  if(btn){ btn.style.background = hmShowSpLogos ? "#d0e8ff" : "#fff"; btn.style.fontWeight = hmShowSpLogos ? "600" : ""; }
  drawCladogram();
}

function hmToggleGroupByHG(){
  hmGroupByHG = !hmGroupByHG;
  hmGroupOrderOverride = null;
  const btn = document.getElementById("btn-hm-group-hg");
  if(btn){
    btn.style.background   = hmGroupByHG ? "#e8f0fe" : "#fff";
    btn.style.color        = hmGroupByHG ? "#1a56c4" : "#555";
    btn.style.borderColor  = hmGroupByHG ? "#1a56c4" : "#888";
    btn.style.fontWeight   = hmGroupByHG ? "600" : "";
  }
  drawHeatmap();
}

function hmToggleFlip(){
  hmFlipped = !hmFlipped;
  const tp = document.getElementById("tree-panel");
  if(tp) tp.style.display = hmFlipped ? "none" : "";
  const btn = document.getElementById("btn-hm-flip");
  if(btn){
    btn.style.background  = hmFlipped ? "#e8f0fe" : "#fff";
    btn.style.color       = hmFlipped ? "#1a56c4" : "#555";
    btn.style.borderColor = hmFlipped ? "#1a56c4" : "#888";
    btn.style.fontWeight  = hmFlipped ? "600" : "";
  }
  drawHeatmap();
}

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
  if (name==="align" && !HAVE_ALIGNMENTS) return;
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
  } else if (name==="architectures") {
    tc.style.display = "none"; pfx.style.display = "none";
    hb.style.display = "none"; cr.textContent = "";
    if (!_architectureSidebarBuilt) { renderArchitectureFamilySidebar(""); _architectureSidebarBuilt = true; }
    else renderArchitectureFamilySidebar(document.getElementById('arch-family-search')?.value || "");
    if (!architectureCurrentFamily && _architectureFamilyCount()) selectArchitectureFamily(_loadArchitectureData().families[0]);
    else { renderArchitectureList(); renderArchitectureDetail(); }
  } else if (name==="align") {
    tc.style.display = "none"; pfx.style.display = "none";
    hb.style.display = "none"; cr.textContent = "";
    if (!_alnSidebarBuilt) { renderAlnSidebar(""); _alnSidebarBuilt = true; }
    if (!currentAlnId && TREE_INDEX.length) selectAlignment(TREE_INDEX[0].id);
  } else {
    tc.style.display = "none"; pfx.style.display = "";
    drawHeatmap(); // drawCladogram() is called at the end of drawHeatmap()
  }
}

// ═══════════════════════════════════════════════════════════════════════════════
// DOMAIN ARCHITECTURE VIEWER
// ═══════════════════════════════════════════════════════════════════════════════

let _architectureSidebarBuilt = false;
let _architectureDataLoaded = false;
let _architectureData = null;
let _architectureFamilyIndex = new Map();
let architectureCurrentFamily = null;
let architectureCurrentIndex = null;
let architectureCurrentDomain = null;
let architectureShowSingletons = false;

function _archEmptyCatalog() {
  return {families:[], domains:[], species:[], hgs:[], entries:[]};
}

function _loadArchitectureData() {
  if (_architectureDataLoaded) return _architectureData;
  _architectureDataLoaded = true;
  const tel = document.getElementById('architecture-data');
  if (!tel) { _architectureData = _archEmptyCatalog(); return _architectureData; }
  let payload = {};
  try { payload = JSON.parse(tel.textContent || '{}'); } catch (_) { payload = {}; }
  if (!payload.gz) {
    _architectureData = _archEmptyCatalog();
    return _architectureData;
  }
  try {
    _architectureData = JSON.parse(_decompressGz(payload.gz));
  } catch (err) {
    console.error('Failed to load architecture catalog', err);
    _architectureData = _archEmptyCatalog();
  }
  _architectureFamilyIndex = new Map((_architectureData.families || []).map((f, i) => [f, i]));
  return _architectureData;
}

function _architectureFamilyCount() {
  return (_loadArchitectureData().families || []).length;
}

function _architectureFamilyEntry(family) {
  const data = _loadArchitectureData();
  const idx = _architectureFamilyIndex.get(family);
  if (idx == null) return null;
  const domains = data.domains || [];
  const species = data.species || [];
  const hgs = data.hgs || [];
  const archs = (data.entries?.[idx] || []).map(rec => ({
    key: rec[0] || [],
    count: Number(rec[1] || 0),
    species: (rec[2] || []).map(i => species[i] || '').filter(Boolean),
    hgs: (rec[3] || []).map(row => {
      const hg = hgs[row[0]] || '';
      const nArch = Number(row[1] || 0);
      const nTotal = Number(row[2] || 0);
      const pct = nTotal > 0 ? (100 * nArch / nTotal) : 0;
      return {hg, nArch, nTotal, pct};
    }).filter(rec => rec.hg),
    consensus: (rec[4] || []).map(hit => ({
      start: Number(hit[0] || 0),
      end: Number(hit[1] || 0),
      name: domains[hit[2]] || '',
    })),
    labels: (rec[0] || []).map(i => domains[i] || '').filter(Boolean),
  }));
  return {family, architectures: archs};
}

function _architectureVisibleArchitectures(entry) {
  const archs = entry?.architectures || [];
  return architectureShowSingletons ? archs : archs.filter(rec => rec.count > 1);
}

function _architectureFilteredArchitectures(entry) {
  const visible = _architectureVisibleArchitectures(entry);
  if (!architectureCurrentDomain) return [];
  return visible.filter(rec => (rec.labels || []).includes(architectureCurrentDomain));
}

function _architectureDomainGraph(entry) {
  const archs = _architectureVisibleArchitectures(entry);
  const nodeCounts = new Map();
  const edgeCounts = new Map();
  for (const arch of archs) {
    const uniq = [...new Set((arch.labels || []).filter(Boolean))];
    for (const label of uniq) {
      nodeCounts.set(label, (nodeCounts.get(label) || 0) + arch.count);
    }
    for (let i = 0; i < uniq.length; i++) {
      for (let j = i + 1; j < uniq.length; j++) {
        const a = uniq[i];
        const b = uniq[j];
        const key = a < b ? `${a}\t${b}` : `${b}\t${a}`;
        edgeCounts.set(key, (edgeCounts.get(key) || 0) + arch.count);
      }
    }
  }
  const nodes = [...nodeCounts.entries()]
    .map(([id, count]) => ({id, count}))
    .sort((a, b) => (b.count - a.count) || a.id.localeCompare(b.id));
  const nodeIndex = new Map(nodes.map((node, idx) => [node.id, idx]));
  const links = [...edgeCounts.entries()]
    .map(([key, count]) => {
      const [source, target] = key.split('\t');
      return {source, target, count};
    })
    .filter(link => nodeIndex.has(link.source) && nodeIndex.has(link.target))
    .sort((a, b) => (b.count - a.count) || a.source.localeCompare(b.source) || a.target.localeCompare(b.target));
  return {nodes, links, totalArchitectures: archs.length};
}

function _architectureDefaultDomain(entry) {
  const graph = _architectureDomainGraph(entry);
  return graph.nodes.length ? graph.nodes[0].id : null;
}

function _archSeriesText(labels) {
  return (labels && labels.length) ? labels.join(' \u2192 ') : 'No domains';
}

function _consecutiveRuns(items, keyFn) {
  const out = [];
  for (const item of (items || [])) {
    const key = keyFn(item);
    const prev = out[out.length - 1];
    if (prev && prev.key === key) prev.items.push(item);
    else out.push({key, items:[item]});
  }
  return out;
}

function _archStackChipHtml(label, count) {
  const bg = _proteinColor((label || '') + '');
  if ((count || 1) <= 1) {
    return `<span class="arch-domain-chip" style="background:${bg}">${_proteinEsc(label || 'domain')}</span>`;
  }
  const shift = 7;
  const depth = Math.max(0, (count || 1) - 1);
  const copies = Array.from({length: depth}, (_, idx) => {
    const left = (depth - idx) * shift;
    return `<span class="arch-domain-chip-copy" style="left:${left}px"><span class="arch-domain-chip" style="background:${bg}">${_proteinEsc(label || 'domain')}</span></span>`;
  }).join('');
  return `<span class="arch-domain-stack" style="padding-right:${depth * shift}px">${copies}<span class="arch-domain-chip arch-domain-chip-front" style="background:${bg}">${_proteinEsc(label || 'domain')}</span></span>`;
}

function _archSeriesChipHtml(labels) {
  if (!labels || !labels.length) {
    return '<span class="arch-domain-chip empty">No domains</span>';
  }
  return _consecutiveRuns(labels, label => label).map((run, idx) => {
    const chip = _archStackChipHtml(run.key, run.items.length);
    if (idx === 0) return chip;
    return `<span class="arch-domain-link"></span>${chip}`;
  }).join('');
}

function _archPctColor(pct) {
  const p = Math.max(0, Math.min(100, Number(pct) || 0));
  const hue = (p / 100) * 120;
  return `hsl(${hue}, 68%, 58%)`;
}

function _archSvg(consensus, labels) {
  const W = 980, H = 82, left = 90, right = 24, plotW = W - left - right;
  const rowY = 36, rowH = 18, lineY = rowY + rowH / 2;
  const scale = (v) => left + (Math.max(0, Math.min(1000, v || 0)) / 1000) * plotW;
  const parts = [];
  parts.push(`<svg class="arch-svg" viewBox="0 0 ${W} ${H}" width="100%" height="${H}" xmlns="http://www.w3.org/2000/svg">`);
  parts.push(`<text x="${left - 10}" y="${rowY + 6}" text-anchor="end" class="arch-row-label">Consensus</text>`);
  parts.push(`<line x1="${left}" y1="${lineY}" x2="${left + plotW}" y2="${lineY}" class="arch-connector"/>`);
  for (const hit of consensus) {
    const x = scale(hit.start);
    const w = Math.max(8, scale(hit.end) - x);
    const color = _proteinColor((hit.name || '') + '');
    const label = _proteinEsc(hit.name || 'domain');
    parts.push(`<g>`);
    parts.push(`<title>${label}: ${hit.start / 10}% - ${hit.end / 10}%</title>`);
    parts.push(`<rect class="arch-domain-box" x="${x.toFixed(1)}" y="${rowY}" width="${w.toFixed(1)}" height="${rowH}" rx="1" ry="1" fill="${color}"/>`);
    if (w > 72) {
      parts.push(`<text x="${(x + w / 2).toFixed(1)}" y="${rowY + 12}" text-anchor="middle" class="arch-domain-label">${label}</text>`);
    } else {
      const lx = Math.max(left + 4, Math.min(left + plotW - 4, x + w / 2));
      parts.push(`<text x="${lx.toFixed(1)}" y="${rowY - 8}" text-anchor="middle" class="arch-domain-label above">${label}</text>`);
    }
    parts.push(`</g>`);
  }
  parts.push(`</svg>`);
  return parts.join('');
}

function renderArchitectureFamilySidebar(query="") {
  const list = document.getElementById('arch-family-list');
  const count = document.getElementById('arch-family-count');
  if (!list || !count) return;
  if (!HAVE_ARCHITECTURES) {
    count.textContent = 'No exact-hit architectures available';
    list.innerHTML = '';
    return;
  }
  const data = _loadArchitectureData();
  const families = data.families || [];
  const q = (query || '').trim().toLowerCase();
  const matches = q ? families.filter(f => f.toLowerCase().includes(q)) : families;
  count.textContent = `${matches.length} families`;
  list.innerHTML = matches.map(fam => {
    const entry = _architectureFamilyEntry(fam);
    const nArch = (entry?.architectures || []).length;
    const famToken = encodeURIComponent(fam);
    return `<div class="arch-item ${fam===architectureCurrentFamily?'active':''}" data-family="${_proteinEsc(famToken)}" onclick="selectArchitectureFamily(decodeURIComponent(this.dataset.family || ''))"><span>${_proteinEsc(fam)}</span><span class="arch-badge">${nArch} arch${nArch===1?'':'s'}</span></div>`;
  }).join('');
}

function selectArchitectureFamily(family) {
  architectureCurrentFamily = family;
  architectureCurrentIndex = null;
  const entry = _architectureFamilyEntry(family);
  architectureCurrentDomain = _architectureDefaultDomain(entry);
  renderArchitectureFamilySidebar(document.getElementById('arch-family-search')?.value || '');
  renderArchitectureList();
  renderArchitectureDetail();
}

function renderArchitectureList() {
  const list = document.getElementById('arch-arch-list');
  const count = document.getElementById('arch-arch-count');
  if (!list || !count) return;
  if (!HAVE_ARCHITECTURES || !architectureCurrentFamily) {
    count.textContent = 'Select a family';
    list.innerHTML = '';
    return;
  }
  const entry = _architectureFamilyEntry(architectureCurrentFamily);
  if (!architectureCurrentDomain) {
    count.textContent = `Click a PFAM node to list architectures${!architectureShowSingletons ? ' (singletons hidden)' : ''}`;
    list.innerHTML = '';
    architectureCurrentIndex = null;
    return;
  }
  const archs = _architectureFilteredArchitectures(entry);
  const domainTotal = archs.reduce((sum, rec) => sum + (rec.count || 0), 0);
  count.textContent = `${archs.length} architectures containing ${architectureCurrentDomain} across ${domainTotal} proteins${!architectureShowSingletons ? ' (singletons hidden)' : ''}`;
  list.innerHTML = archs.map((rec, idx) => {
    const series = _archSeriesChipHtml(rec.labels);
    const active = architectureCurrentIndex === idx;
    const pct = domainTotal > 0 ? (100 * rec.count / domainTotal) : 0;
    const pctTxt = `${pct.toFixed(pct >= 10 ? 0 : 1)}%`;
    const badgeBg = _archPctColor(pct);
    return `
      <div class="arch-item ${active?'active':''}" onclick="selectArchitecture(${idx})">
        <div class="arch-item-main" title="${_proteinEsc(_archSeriesText(rec.labels))}">
          <div class="arch-chip-list">${series}</div>
        </div>
        <span class="arch-badge" style="background:${badgeBg};border-color:rgba(20,20,20,.28);color:#102030;font-weight:700" title="${rec.count} proteins, ${pctTxt} of proteins carrying ${_proteinEsc(architectureCurrentDomain)}">${rec.count} · ${pctTxt}</span>
      </div>`;
  }).join('');
  if (archs.length && (architectureCurrentIndex == null || architectureCurrentIndex >= archs.length)) {
    architectureCurrentIndex = 0;
  }
}

function selectArchitecture(idx) {
  architectureCurrentIndex = idx;
  renderArchitectureDetail();
  renderArchitectureList();
}

function selectArchitectureDomain(label) {
  const domain = label || null;
  architectureCurrentDomain = (architectureCurrentDomain === domain) ? null : domain;
  architectureCurrentIndex = null;
  renderArchitectureList();
  renderArchitectureDetail();
}

function _architectureTreeRecForHg(hgId) {
  return TREE_INDEX.find(rec => rec.id === hgId || rec.hg === hgId) || null;
}

function architectureOpenCountsForHg(hgId) {
  const rec = _architectureTreeRecForHg(hgId);
  const hgRec = HG_DATA.find(rec2 => rec2.id === hgId || rec2.hg === hgId) || null;
  const family = rec?.family || hgRec?.family || null;
  if (family) famGoToFamily(family);
}

function showArchitectureHgPopup(ev, hgId) {
  const rec = _architectureTreeRecForHg(hgId);
  showHmOpenPopup(
    ev,
    hgId,
    () => architectureOpenCountsForHg(hgId),
    () => hmOpenAlignmentForHG(rec?.id || hgId),
    'Show Counts',
    'Show Alignment'
  );
}

function renderArchitectureDetail() {
  const empty = document.getElementById('arch-empty');
  const view = document.getElementById('arch-view');
  const summary = document.getElementById('arch-summary');
  const networkWrap = document.getElementById('arch-network');
  const detailEmpty = document.getElementById('arch-detail-empty');
  const detailWrap = document.getElementById('arch-detail');
  const svgWrap = document.getElementById('arch-svg-wrap');
  const hgWrap = document.getElementById('arch-hgs');
  if (!empty || !view || !summary || !networkWrap || !detailEmpty || !detailWrap || !svgWrap || !hgWrap) return;
  if (!HAVE_ARCHITECTURES || !architectureCurrentFamily) {
    empty.style.display = 'flex';
    view.style.display = 'none';
    return;
  }
  const entry = _architectureFamilyEntry(architectureCurrentFamily);
  const visibleArchs = _architectureVisibleArchitectures(entry);
  const graph = _architectureDomainGraph(entry);
  const archs = _architectureFilteredArchitectures(entry);
  const arch = archs[architectureCurrentIndex ?? 0];
  if (!arch) {
    architectureCurrentIndex = null;
  }
  empty.style.display = 'none';
  view.style.display = 'block';
  summary.innerHTML = [
    `<div class="arch-summary-chip"><strong>${_proteinEsc(architectureCurrentFamily)}</strong></div>`,
    `<div class="arch-summary-chip">PFAM domains: <strong>${graph.nodes.length}</strong></div>`,
    `<div class="arch-summary-chip">Co-occurrence edges: <strong>${graph.links.length}</strong></div>`,
    `<div class="arch-summary-chip">Architectures: <strong>${visibleArchs.length}</strong></div>`,
    architectureCurrentDomain ? `<div class="arch-summary-chip">Selected PFAM: <strong>${_proteinEsc(architectureCurrentDomain)}</strong></div>` : '',
  ].join('');
  renderArchitectureNetwork(networkWrap, graph);
  if (!arch) {
    detailEmpty.style.display = 'block';
    detailWrap.style.display = 'none';
    return;
  }
  detailEmpty.style.display = 'none';
  detailWrap.style.display = 'block';
  summary.innerHTML = [
    `<div class="arch-summary-chip"><strong>${_proteinEsc(architectureCurrentFamily)}</strong></div>`,
    architectureCurrentDomain ? `<div class="arch-summary-chip">Selected PFAM: <strong>${_proteinEsc(architectureCurrentDomain)}</strong></div>` : '',
    `<div class="arch-summary-chip">Architecture: <strong>${_proteinEsc(_archSeriesText(arch.labels))}</strong></div>`,
    `<div class="arch-summary-chip">Proteins: <strong>${arch.count}</strong></div>`,
    `<div class="arch-summary-chip">Species: <strong>${arch.species.length}</strong></div>`,
    `<div class="arch-summary-chip">HGs: <strong>${arch.hgs.length}</strong></div>`,
  ].join('');
  svgWrap.innerHTML = _archSvg(arch.consensus, arch.labels);
  hgWrap.innerHTML = arch.hgs.map(rec => {
    const pct = `${rec.pct.toFixed(rec.pct >= 10 ? 0 : 1)}%`;
    const archShare = arch.count > 0 ? (100 * rec.nArch / arch.count) : 0;
    const archShareTxt = `${archShare.toFixed(archShare >= 10 ? 0 : 1)}%`;
    const pctBg = _archPctColor(rec.pct);
    const hgToken = encodeURIComponent(rec.hg);
    return `
      <div class="arch-hg-row" data-hg="${_proteinEsc(hgToken)}" onclick="showArchitectureHgPopup(event, decodeURIComponent(this.dataset.hg || ''))">
        <div class="arch-hg-main">
          <div class="arch-hg-title">${_proteinEsc(rec.hg)}</div>
          <div class="arch-hg-meta">Within HG: ${rec.nArch} / ${rec.nTotal} proteins · Of architecture: ${rec.nArch} / ${arch.count} proteins</div>
        </div>
        <div class="arch-hg-badges">
          <div class="arch-badge" style="background:${pctBg};border-color:rgba(20,20,20,.28);color:#102030;font-weight:700" title="Share within this HG">${pct}</div>
          <div class="arch-share-badge" title="Share of this architecture found in the HG">${archShareTxt}</div>
        </div>
      </div>`;
  }).join('');
}

function renderArchitectureNetwork(container, graph) {
  if (!container) return;
  container.innerHTML = '';
  const width = Math.max(480, container.clientWidth || 960);
  const height = Math.max(340, Math.min(620, 220 + Math.sqrt(Math.max(1, graph.nodes.length)) * 54));
  const svg = d3.select(container)
    .append('svg')
    .attr('class', 'arch-net-svg')
    .attr('viewBox', `0 0 ${width} ${height}`)
    .attr('width', '100%')
    .attr('height', height);
  if (!graph.nodes.length) {
    svg.append('text')
      .attr('x', width / 2)
      .attr('y', height / 2)
      .attr('text-anchor', 'middle')
      .attr('fill', '#7b8794')
      .attr('font-size', 13)
      .text('No PFAM co-occurrence data available for this family');
    return;
  }
  const nodeMax = d3.max(graph.nodes, d => d.count) || 1;
  const edgeMax = d3.max(graph.links, d => d.count) || 1;
  const rScale = d3.scaleSqrt().domain([1, nodeMax]).range([9, 30]);
  const wScale = d3.scaleSqrt().domain([1, edgeMax]).range([1.2, 6]);
  const nodes = graph.nodes.map(node => ({...node}));
  const links = graph.links.map(link => ({...link}));
  const sim = d3.forceSimulation(nodes)
    .force('link', d3.forceLink(links).id(d => d.id).distance(link => Math.max(42, 170 - wScale(link.count) * 18)).strength(0.34))
    .force('charge', d3.forceManyBody().strength(d => -180 - rScale(d.count) * 24))
    .force('center', d3.forceCenter(width / 2, height / 2))
    .force('x', d3.forceX(width / 2).strength(0.03))
    .force('y', d3.forceY(height / 2).strength(0.03))
    .force('collide', d3.forceCollide().radius(d => rScale(d.count) + 34));
  for (let i = 0; i < 280; i++) sim.tick();
  sim.stop();
  const pad = 28;
  const extents = nodes.map(d => ({
    minX: d.x - rScale(d.count),
    maxX: d.x + rScale(d.count) + 8 + Math.max(24, String(d.id || '').length * 7),
    minY: d.y - rScale(d.count) - 14,
    maxY: d.y + rScale(d.count) + 14,
  }));
  const minX = d3.min(extents, d => d.minX) ?? 0;
  const maxX = d3.max(extents, d => d.maxX) ?? width;
  const minY = d3.min(extents, d => d.minY) ?? 0;
  const maxY = d3.max(extents, d => d.maxY) ?? height;
  const spanW = Math.max(1, maxX - minX);
  const spanH = Math.max(1, maxY - minY);
  const fitScale = Math.min((width - pad * 2) / spanW, (height - pad * 2) / spanH, 1.25);
  const tx = (width - spanW * fitScale) / 2 - minX * fitScale;
  const ty = (height - spanH * fitScale) / 2 - minY * fitScale;
  const baseTransform = d3.zoomIdentity.translate(tx, ty).scale(fitScale);
  const viewport = svg.append('g').attr('transform', baseTransform.toString());
  const plot = viewport.append('g');
  plot.append('g')
    .selectAll('line')
    .data(links)
    .enter()
    .append('line')
    .attr('class', 'arch-net-link')
    .attr('stroke-width', d => wScale(d.count))
    .attr('x1', d => d.source.x)
    .attr('y1', d => d.source.y)
    .attr('x2', d => d.target.x)
    .attr('y2', d => d.target.y)
    .append('title')
    .text(d => `${d.source.id} + ${d.target.id}: ${d.count} proteins`);
  const nodeG = plot.append('g')
    .selectAll('g')
    .data(nodes)
    .enter()
    .append('g')
    .attr('transform', d => `translate(${d.x},${d.y})`)
    .style('cursor', 'pointer')
    .on('click', (event, d) => selectArchitectureDomain(d.id));
  nodeG.append('circle')
    .attr('class', d => `arch-net-node${architectureCurrentDomain === d.id ? ' active' : ''}`)
    .attr('r', d => rScale(d.count))
    .attr('fill', d => _proteinColor(d.id + ''))
    .append('title')
    .text(d => `${d.id}: ${d.count} proteins`);
  nodeG.append('text')
    .attr('class', 'arch-net-label')
    .attr('x', d => rScale(d.count) + 6)
    .attr('y', 4)
    .text(d => d.id);
  const zoom = d3.zoom()
    .scaleExtent([0.45, 4])
    .on('zoom', event => {
      viewport.attr('transform', event.transform.toString());
    });
  svg.call(zoom).call(zoom.transform, baseTransform);
}

// ═══════════════════════════════════════════════════════════════════════════════
// PROTEIN DOMAIN VIEWER
// ═══════════════════════════════════════════════════════════════════════════════

let _proteinSidebarBuilt = false;
let _proteinDomainDataLoaded = false;
let _proteinDomainData = null;
let _proteinDomainIndex = new Map();
const _treeDetailCache = new Map();
let proteinCurrentGene = null;

function _proteinEmptyCatalog() {
  return {genes:[], lengths:[], tracks:[], families:[], names:[], pfams:[]};
}

function _proteinEsc(text) {
  return String(text == null ? '' : text)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

function _loadProteinDomainData() {
  if (_proteinDomainDataLoaded) return _proteinDomainData;
  _proteinDomainDataLoaded = true;
  const tel = document.getElementById('protein-domain-data');
  if (!tel) { _proteinDomainData = _proteinEmptyCatalog(); return _proteinDomainData; }
  let payload = {};
  try { payload = JSON.parse(tel.textContent || '{}'); } catch (_) { payload = {}; }
  if (!payload.gz) {
    _proteinDomainData = _proteinEmptyCatalog();
    return _proteinDomainData;
  }
  try {
    _proteinDomainData = JSON.parse(_decompressGz(payload.gz));
  } catch (err) {
    console.error('Failed to load protein-domain catalog', err);
    _proteinDomainData = _proteinEmptyCatalog();
  }
  _proteinDomainIndex = new Map((_proteinDomainData.genes || []).map((g, i) => [g, i]));
  return _proteinDomainData;
}

function _proteinGeneCount() {
  return (_loadProteinDomainData().genes || []).length;
}

function _proteinEntry(gid) {
  const data = _loadProteinDomainData();
  const idx = _proteinDomainIndex.get(gid);
  if (idx == null) return null;
  const families = data.families || [];
  const names = data.names || [];
  const pfams = data.pfams || [];
  const tracks = (data.tracks?.[idx] || []).map(rec => ({
    family: families[rec[0]] || '',
    rangeStart: Number.isFinite(+rec[1]) ? +rec[1] : null,
    rangeEnd: Number.isFinite(+rec[2]) ? +rec[2] : null,
    hits: (rec[3] || []).map(hit => ({
      start: Number.isFinite(+hit[0]) ? +hit[0] : null,
      end: Number.isFinite(+hit[1]) ? +hit[1] : null,
      name: names[hit[2]] || '',
      pfam: pfams[hit[3]] || '',
    })),
  }));
  return {
    gene: gid,
    length: Number.isFinite(+data.lengths?.[idx]) && +data.lengths[idx] > 0 ? +data.lengths[idx] : null,
    tracks,
  };
}

function _proteinTracksForFamily(entry, family) {
  if (!entry || !entry.tracks || !entry.tracks.length) return [];
  const fam = (family || '').trim();
  if (!fam) return entry.tracks.slice();
  const matched = entry.tracks.filter(track => track.family === fam);
  if (matched.length) return matched;
  return entry.tracks.length === 1 ? entry.tracks.slice() : [];
}

function _proteinColor(key) {
  let hash = 0;
  for (let i = 0; i < key.length; i++) hash = ((hash * 131) + key.charCodeAt(i)) % 360;
  return `hsl(${hash},62%,60%)`;
}

function _proteinRulerTickStep(total) {
  if (total <= 300) return 25;
  if (total <= 800) return 50;
  if (total <= 2000) return 100;
  if (total <= 5000) return 200;
  return 500;
}

function _proteinDomainLegendHtml(track) {
  const seen = new Set();
  const items = [];
  for (const hit of (track?.hits || [])) {
    const key = `${hit.pfam || ''}\t${hit.name || ''}`;
    if (seen.has(key)) continue;
    seen.add(key);
    items.push({
      name: hit.name || 'domain',
      pfam: hit.pfam || '',
      color: _proteinColor((hit.pfam || hit.name || track.family || '') + ''),
    });
  }
  if (!items.length) return '';
  return `<div class="prot-domain-legend">${items.map(item => (
    `<span class="prot-domain-chip"><span class="prot-domain-chip-swatch" style="background:${item.color}"></span>${_proteinEsc(item.name)}${item.pfam ? ` <span style="color:#6f8090">(${_proteinEsc(item.pfam)})</span>` : ''}</span>`
  )).join('')}</div>`;
}

function _proteinTrackSvg(length, track, opts={}) {
  const sortedHits = [...(track.hits || [])]
    .filter(hit => hit.start != null && hit.end != null)
    .sort((a, b) => (a.start - b.start) || (a.end - b.end) || String(a.name || '').localeCompare(String(b.name || '')));
  const maxHitEnd = sortedHits.reduce((m, h) => Math.max(m, h.end || 0), 0);
  const total = Math.max(length || 0, track.rangeEnd || 0, maxHitEnd || 0, 1);
  const hitSpanStart = sortedHits.length ? Math.min(...sortedHits.map(h => h.start || total)) : null;
  const hitSpanEnd = sortedHits.length ? Math.max(...sortedHits.map(h => h.end || 0)) : null;
  const focusToPhylogeny = !!opts.focusToPhylogeny;
  const displayStart = focusToPhylogeny
    ? Math.max(1, Math.min(total, track.rangeStart || hitSpanStart || 1))
    : 1;
  const displayEnd = focusToPhylogeny
    ? Math.max(displayStart, Math.min(total, track.rangeEnd || hitSpanEnd || total))
    : total;
  const displaySpan = Math.max(1, displayEnd - displayStart + 1);
  const compact = !!opts.compact;
  const W = compact ? 720 : 980, left = compact ? 148 : 120, right = compact ? 22 : 24, plotW = W - left - right;
  const phyloY = 18, phyloH = 16;
  const hitSpanY = 48, hitSpanH = 12;
  const hitTop = 76, hitH = 18, laneGap = 7;
  const scale = (v) => {
    const clamped = Math.max(displayStart, Math.min(displayEnd, v || displayStart));
    if (displaySpan <= 1) return left;
    return left + ((clamped - displayStart) / (displaySpan - 1)) * plotW;
  };
  const rangeStart = track.rangeStart || displayStart;
  const rangeEnd = track.rangeEnd || displayEnd;
  const rangeX = scale(rangeStart);
  const rangeW = Math.max(8, scale(rangeEnd) - rangeX);
  const tickStep = _proteinRulerTickStep(displaySpan);
  const mediumStep = tickStep < 50 ? tickStep * 2 : 50;
  const majorStep = tickStep > 100 ? tickStep : (displaySpan >= 100 ? 100 : tickStep * 2);
  const laneEnds = [];
  const laidOutHits = sortedHits.map(hit => {
    let lane = 0;
    while (lane < laneEnds.length && hit.start <= laneEnds[lane]) lane += 1;
    laneEnds[lane] = hit.end;
    return {...hit, lane};
  });
  const laneCount = Math.max(1, laneEnds.length);
  const hitBandH = laneCount * hitH + Math.max(0, laneCount - 1) * laneGap;
  const rulerY = hitTop + hitBandH + 30;
  const H = rulerY + 38;
  const labelX = compact ? 18 : 14;
  const parts = [];
  parts.push(`<svg class="prot-track-svg${compact ? ' compact' : ''}" viewBox="0 0 ${W} ${H}" width="100%" height="${H}" xmlns="http://www.w3.org/2000/svg">`);
  parts.push(`<text x="${labelX}" y="${phyloY + 11}" class="prot-row-label">Phylogeny span</text>`);
  parts.push(`<text x="${labelX}" y="${hitSpanY + 9}" class="prot-row-label">Hit span</text>`);
  parts.push(`<text x="${labelX}" y="${hitTop + 12}" class="prot-row-label">Exact hits</text>`);
  parts.push(`<text x="${labelX}" y="${rulerY + 4}" class="prot-row-label">${focusToPhylogeny ? 'Phylogeny window' : 'Protein'}</text>`);
  parts.push(`<rect x="${left}" y="${phyloY + 4}" width="${plotW}" height="8" fill="#edf2f7"/>`);
  parts.push(`<rect x="${left}" y="${hitSpanY + 3}" width="${plotW}" height="6" fill="#f2f5f8"/>`);
  parts.push(`<rect x="${left}" y="${hitTop + 5}" width="${plotW}" height="${Math.max(8, hitBandH - 10)}" fill="#f3f6f9"/>`);
  parts.push(`<line x1="${left}" y1="${rulerY}" x2="${left + plotW}" y2="${rulerY}" stroke="#7f8f9f" stroke-width="2"/>`);
  parts.push(`<line x1="${scale(displayStart).toFixed(1)}" y1="${rulerY}" x2="${scale(displayStart).toFixed(1)}" y2="${rulerY + 12}" class="prot-ruler-tick major"/>`);
  parts.push(`<text x="${scale(displayStart).toFixed(1)}" y="${rulerY + 26}" text-anchor="middle" class="prot-ruler-label">${displayStart}</text>`);
  const tickBegin = Math.ceil(displayStart / tickStep) * tickStep;
  for (let pos = tickBegin; pos < displayEnd; pos += tickStep) {
    if (pos <= displayStart) continue;
    const isMajor = majorStep > 0 && pos % majorStep === 0;
    const isMedium = !isMajor && mediumStep > 0 && pos % mediumStep === 0;
    const tickH = isMajor ? 12 : (isMedium ? 8 : 5);
    parts.push(`<line x1="${scale(pos).toFixed(1)}" y1="${rulerY}" x2="${scale(pos).toFixed(1)}" y2="${(rulerY + tickH).toFixed(1)}" class="prot-ruler-tick${isMajor ? ' major' : ''}"/>`);
    if (isMajor) {
      parts.push(`<text x="${scale(pos).toFixed(1)}" y="${rulerY + 26}" text-anchor="middle" class="prot-ruler-label">${pos}</text>`);
    }
  }
  parts.push(`<line x1="${scale(displayEnd).toFixed(1)}" y1="${rulerY}" x2="${scale(displayEnd).toFixed(1)}" y2="${rulerY + 12}" class="prot-ruler-tick major"/>`);
  if (displayEnd > displayStart) {
    parts.push(`<text x="${scale(displayEnd).toFixed(1)}" y="${rulerY + 26}" text-anchor="middle" class="prot-ruler-label">${displayEnd}</text>`);
  }
  parts.push(`<g class="prot-span-group">`);
  parts.push(`<rect class="prot-phylo-span" x="${rangeX.toFixed(1)}" y="${phyloY}" width="${rangeW.toFixed(1)}" height="${phyloH}" rx="1" ry="1"/>`);
  if (hitSpanStart != null && hitSpanEnd != null) {
    const p1 = scale(rangeStart);
    const p2 = scale(rangeEnd);
    parts.push(`<line x1="${p1.toFixed(1)}" y1="${phyloY + phyloH}" x2="${p1.toFixed(1)}" y2="${rulerY}" class="prot-span-guide"/>`);
    parts.push(`<line x1="${p2.toFixed(1)}" y1="${phyloY + phyloH}" x2="${p2.toFixed(1)}" y2="${rulerY}" class="prot-span-guide"/>`);
    parts.push(`<text x="${p1.toFixed(1)}" y="${rulerY - 20}" text-anchor="middle" class="prot-span-coord">${rangeStart}</text>`);
    parts.push(`<text x="${p2.toFixed(1)}" y="${rulerY - 20}" text-anchor="middle" class="prot-span-coord">${rangeEnd}</text>`);
  }
  parts.push(`</g>`);
  if (hitSpanStart != null && hitSpanEnd != null) {
    const hitSpanX = scale(hitSpanStart);
    const hitSpanW = Math.max(8, scale(hitSpanEnd) - hitSpanX);
    const h1 = scale(hitSpanStart);
    const h2 = scale(hitSpanEnd);
    parts.push(`<g class="prot-span-group">`);
    parts.push(`<rect class="prot-hit-span" x="${hitSpanX.toFixed(1)}" y="${hitSpanY}" width="${hitSpanW.toFixed(1)}" height="${hitSpanH}" rx="1" ry="1"/>`);
    parts.push(`<line x1="${h1.toFixed(1)}" y1="${hitSpanY + hitSpanH}" x2="${h1.toFixed(1)}" y2="${rulerY}" class="prot-span-guide"/>`);
    parts.push(`<line x1="${h2.toFixed(1)}" y1="${hitSpanY + hitSpanH}" x2="${h2.toFixed(1)}" y2="${rulerY}" class="prot-span-guide"/>`);
    parts.push(`<text x="${h1.toFixed(1)}" y="${rulerY - 8}" text-anchor="middle" class="prot-span-coord">${hitSpanStart}</text>`);
    parts.push(`<text x="${h2.toFixed(1)}" y="${rulerY - 8}" text-anchor="middle" class="prot-span-coord">${hitSpanEnd}</text>`);
    parts.push(`</g>`);
  }
  for (const hit of laidOutHits) {
    const x = scale(hit.start);
    const w = Math.max(8, scale(hit.end) - x);
    const color = _proteinColor((hit.pfam || hit.name || track.family) + '');
    const label = _proteinEsc(hit.name || hit.pfam || 'domain');
    const pfam = _proteinEsc(hit.pfam || '');
    const x2 = x + w;
    const y = hitTop + hit.lane * (hitH + laneGap);
    parts.push(`<g class="prot-hit-group">`);
    parts.push(`<title>${label}${pfam ? ` (${pfam})` : ''}: ${hit.start}-${hit.end}</title>`);
    parts.push(`<line x1="${x.toFixed(1)}" y1="${(y + hitH).toFixed(1)}" x2="${x.toFixed(1)}" y2="${rulerY}" class="prot-hit-guide"/>`);
    parts.push(`<line x1="${x2.toFixed(1)}" y1="${(y + hitH).toFixed(1)}" x2="${x2.toFixed(1)}" y2="${rulerY}" class="prot-hit-guide"/>`);
    parts.push(`<rect class="prot-hit-box" x="${x.toFixed(1)}" y="${y}" width="${w.toFixed(1)}" height="${hitH}" rx="1" ry="1" fill="${color}"/>`);
    if (w > 72) {
      parts.push(`<text x="${(x + w / 2).toFixed(1)}" y="${y + 12}" text-anchor="middle" font-size="10" fill="#102030" font-weight="700">${label}</text>`);
    } else {
      parts.push(`<text x="${(x2 + 8).toFixed(1)}" y="${y + 12}" text-anchor="start" class="prot-hit-hover-label">${label}</text>`);
    }
    parts.push(`<text x="${x.toFixed(1)}" y="${rulerY - 8}" text-anchor="middle" class="prot-hit-coord">${hit.start}</text>`);
    parts.push(`<text x="${x2.toFixed(1)}" y="${rulerY - 8}" text-anchor="middle" class="prot-hit-coord">${hit.end}</text>`);
    parts.push(`</g>`);
  }
  parts.push(`<text x="${left + plotW + 10}" y="${phyloY + 12}" class="prot-ruler-label">${rangeStart}-${rangeEnd}</text>`);
  if (hitSpanStart != null && hitSpanEnd != null) {
    parts.push(`<text x="${left + plotW + 10}" y="${hitSpanY + 10}" class="prot-ruler-label">${hitSpanStart}-${hitSpanEnd}</text>`);
  }
  parts.push(`<text x="${left + plotW + 10}" y="${rulerY + 4}" class="prot-ruler-label">${focusToPhylogeny ? `${displayStart}-${displayEnd} aa window` : `${total} aa`}</text>`);
  parts.push(`</svg>`);
  return parts.join('');
}

function _proteinPopupHtml(gid, opts={}) {
  const entry = _proteinEntry(gid);
  const tracksForPopup = _proteinTracksForFamily(entry, opts.family);
  if (!entry || !tracksForPopup.length) return '';
  const hitCount = tracksForPopup.reduce((n, t) => n + t.hits.length, 0);
  const tracks = tracksForPopup.slice(0, 3).map(track => {
    const phyloRangeTxt = (track.rangeStart != null && track.rangeEnd != null)
      ? `${track.rangeStart}-${track.rangeEnd}`
      : 'NA';
    const hitStarts = track.hits.map(hit => hit.start).filter(v => v != null);
    const hitEnds = track.hits.map(hit => hit.end).filter(v => v != null);
    const hitRangeTxt = (hitStarts.length && hitEnds.length)
      ? `${Math.min(...hitStarts)}-${Math.max(...hitEnds)}`
      : 'NA';
    return `
      <div class="prot-track" style="margin-bottom:8px;padding:8px 10px 10px">
        <div class="prot-track-head" style="margin-bottom:6px">
          <div class="prot-track-title">${_proteinEsc(track.family)}</div>
          <div class="prot-track-meta">Hit span: ${hitRangeTxt} · phylogeny span: ${phyloRangeTxt}</div>
        </div>
        <div class="prot-svg-wrap compact">${_proteinTrackSvg(entry.length, track, {focusToPhylogeny:true, compact:true})}</div>
        ${_proteinDomainLegendHtml(track)}
      </div>`;
  }).join('');
  const hiddenCount = tracksForPopup.length > 3 ? tracksForPopup.length - 3 : 0;
  return `
    <div style="width:min(760px,calc(100vw - 56px));max-width:min(760px,calc(100vw - 56px))">
      <div class="tt-name">${_proteinEsc(entry.gene)}</div>
      <div style="display:flex;gap:8px;flex-wrap:wrap;margin-bottom:8px">
        <span class="prot-summary-chip">Protein length: <strong>${entry.length || 'NA'} aa</strong></span>
        <span class="prot-summary-chip">Family tracks: <strong>${tracksForPopup.length}</strong></span>
        <span class="prot-summary-chip">Exact hits: <strong>${hitCount}</strong></span>
      </div>
      <div class="prot-legend-popup" style="display:flex;flex-wrap:wrap;gap:10px;margin:0 0 8px">
        <div class="prot-legend-item"><span class="prot-legend-swatch phylo"></span><span>Phylogeny span</span></div>
        <div class="prot-legend-item"><span class="prot-legend-swatch hitspan"></span><span>Hit span</span></div>
        <div class="prot-legend-item"><span class="prot-legend-swatch hit"></span><span>Exact hits</span></div>
      </div>
      ${tracks}
      ${hiddenCount ? `<div style="font-size:10px;color:#7f8c8d">+ ${hiddenCount} additional family track${hiddenCount===1?'':'s'}</div>` : ''}
    </div>`;
}

function renderProteinSidebar(query="") {
  const list = document.getElementById('prot-list');
  const count = document.getElementById('prot-count');
  if (!list || !count) return;
  if (!HAVE_PROTEIN_DOMAINS) {
    count.textContent = 'No exact domain-hit data available';
    list.innerHTML = '';
    return;
  }
  const genes = _loadProteinDomainData().genes || [];
  const q = (query || '').trim().toLowerCase();
  const matches = q ? genes.filter(g => g.toLowerCase().includes(q)) : genes;
  const shown = matches.slice(0, 250);
  count.textContent = `${matches.length} proteins${matches.length > shown.length ? ` (showing first ${shown.length})` : ''}`;
  list.innerHTML = shown.map(gid => {
    const idx = _proteinDomainIndex.get(gid);
    const nTracks = ((_proteinDomainData.tracks || [])[idx] || []).length;
    const gidToken = encodeURIComponent(gid);
    return `<div class="prot-item ${gid===proteinCurrentGene?'active':''}" data-gene="${_proteinEsc(gidToken)}" onclick="selectProteinDomainGene(decodeURIComponent(this.dataset.gene || ''))"><span>${_proteinEsc(gid)}</span><span class="prot-badge">${nTracks} track${nTracks===1?'':'s'}</span></div>`;
  }).join('');
}

function selectProteinDomainGene(gid) {
  proteinCurrentGene = gid;
  renderProteinSidebar(document.getElementById('prot-search')?.value || '');
  renderProteinDomainView();
}

function renderProteinDomainView() {
  const empty = document.getElementById('prot-empty');
  const view = document.getElementById('prot-view');
  const summary = document.getElementById('prot-summary');
  const legend = document.getElementById('prot-legend');
  const tracksWrap = document.getElementById('prot-tracks');
  if (!empty || !view || !summary || !legend || !tracksWrap) return;
  if (!HAVE_PROTEIN_DOMAINS || !proteinCurrentGene) {
    empty.style.display = 'flex';
    view.style.display = 'none';
    return;
  }
  const entry = _proteinEntry(proteinCurrentGene);
  if (!entry) {
    empty.style.display = 'flex';
    view.style.display = 'none';
    return;
  }
  empty.style.display = 'none';
  view.style.display = 'block';
  const hitCount = entry.tracks.reduce((n, t) => n + t.hits.length, 0);
  summary.innerHTML = [
    `<div class="prot-summary-chip"><strong>${_proteinEsc(entry.gene)}</strong></div>`,
    `<div class="prot-summary-chip">Protein length: <strong>${entry.length || 'NA'} aa</strong></div>`,
    `<div class="prot-summary-chip">Family tracks: <strong>${entry.tracks.length}</strong></div>`,
    `<div class="prot-summary-chip">Exact hits: <strong>${hitCount}</strong></div>`,
  ].join('');
  legend.innerHTML = [
    `<div class="prot-legend-item"><span class="prot-legend-swatch phylo"></span><span>Phylogeny span: expanded range used downstream in alignments and trees</span></div>`,
    `<div class="prot-legend-item"><span class="prot-legend-swatch hitspan"></span><span>Hit span: min-to-max range covered by the exact hits</span></div>`,
    `<div class="prot-legend-item"><span class="prot-legend-swatch hit"></span><span>Exact hits: individual PFAM domain matches</span></div>`,
    `<div class="prot-legend-item"><span class="prot-legend-swatch guide"></span><span>Dashed guides connect hit boundaries to the protein ruler; hover a hit to show its coordinates</span></div>`,
  ].join('');
  tracksWrap.innerHTML = entry.tracks.map(track => {
    const phyloRangeTxt = (track.rangeStart != null && track.rangeEnd != null)
      ? `${track.rangeStart}-${track.rangeEnd}`
      : 'NA';
    const hitStarts = track.hits.map(hit => hit.start).filter(v => v != null);
    const hitEnds = track.hits.map(hit => hit.end).filter(v => v != null);
    const hitRangeTxt = (hitStarts.length && hitEnds.length)
      ? `${Math.min(...hitStarts)}-${Math.max(...hitEnds)}`
      : 'NA';
    const chips = track.hits.map(hit => {
      const label = _proteinEsc(hit.name || 'domain');
      const pfam = _proteinEsc(hit.pfam || '');
      const loc = `${hit.start || '?'}-${hit.end || '?'}`;
      return `<span class="prot-hit-chip">${label}${pfam ? ` (${pfam})` : ''} <span style="color:#7f8c8d">${loc}</span></span>`;
    }).join('');
    return `
      <div class="prot-track">
        <div class="prot-track-head">
          <div class="prot-track-title">${_proteinEsc(track.family)}</div>
          <div class="prot-track-meta">Hit span: ${hitRangeTxt} · phylogeny span: ${phyloRangeTxt} · exact hits: ${track.hits.length}</div>
        </div>
        <div class="prot-svg-wrap">${_proteinTrackSvg(entry.length, track)}</div>
        <div class="prot-hits">${chips || '<span class="prot-hit-chip">No exact hits recorded</span>'}</div>
      </div>`;
  }).join('');
}

// ═══════════════════════════════════════════════════════════════════════════════
// ALIGNMENT VIEWER
// ═══════════════════════════════════════════════════════════════════════════════

let currentAlnId   = null;
let _alnSidebarBuilt = false;
let alnSeqs        = [];
let alnShowCons    = false;
let alnCellW       = 9;
let alnCellH       = 14;
let _alnRawFasta   = null;
let alnOgMap       = {};
let alnTreeOrder   = [];
let alnTreeDict    = null;
let alnTreeW       = 80;
let alnUseBrLen    = false;
let _alnRows       = [];
let alnSpFilters   = [];   // committed species filter tags
let alnUseGeneRax  = false;
let _alnTreeDataCache = null;  // full treedata object for current HG
let alnGeneMeta     = {};
let alnRangeData    = {};
let alnHiddenOGs    = new Set();
let alnCollapsedOGs = new Set();
let alnShowSeqId    = true;
let alnShowRefCol   = false;
let alnShowSpeciesCol = true;
let alnShowMrcaCol  = true;
let alnShowRangeCol = true;
let alnShowSeqPanel = true;
let alnNameW        = 220;
let alnRefColW      = 100;
let alnSpeciesColW  = 110;
let alnMrcaColW     = 120;
let alnRangeColW    = 96;
let alnDomainChipScale = 1.0;
let alnMetaColOrder = ['ref', 'species', 'mrca', 'range', 'id'];
let _alnSyncingScroll = false;

const AA_COLORS = {
  A:'#C8C8C8',V:'#C8C8C8',I:'#C8C8C8',L:'#C8C8C8',M:'#C8C8C8',
  F:'#FFAAAA',W:'#FFAAAA',Y:'#FFAAAA',
  K:'#6699FF',R:'#6699FF',H:'#88CCEE',
  D:'#FF4444',E:'#FF4444',
  S:'#88CC88',T:'#88CC88',N:'#88CC88',Q:'#88CC88',
  C:'#FFEE44', P:'#EEAAEE', G:'#EEEEEE', '-':'','.':''
};

const _alnOgPalette = [
  '#4e79a7','#f28e2b','#e15759','#76b7b2','#59a14f',
  '#edc948','#b07aa1','#ff9da7','#9c755f','#bab0ac',
  '#af7aa1','#86bcb6','#d4a6c8','#8cd17d','#499894'
];

function parseFasta(text) {
  const seqs = [];
  let name = null, buf = [];
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith('>')) {
      if (name !== null) seqs.push({name, seq: buf.join('')});
      name = line.slice(1).trim();
      buf = [];
    } else { buf.push(line.trim()); }
  }
  if (name !== null) seqs.push({name, seq: buf.join('')});
  return seqs;
}

function _decompressGz(b64) {
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  const inflated = pako.inflate(bytes);
  return new TextDecoder().decode(inflated);
}

function _loadEmbeddedJson(id) {
  if (_treeDetailCache.has(id)) return _treeDetailCache.get(id);
  const el = document.getElementById(id);
  if (!el) return null;
  try {
    let parsed = JSON.parse(el.textContent);
    if (parsed && typeof parsed === 'object' && Object.prototype.hasOwnProperty.call(parsed, 'gz')) {
      parsed = parsed.gz ? JSON.parse(_decompressGz(parsed.gz)) : null;
    }
    _treeDetailCache.set(id, parsed);
    return parsed;
  } catch (e) {
    _treeDetailCache.set(id, null);
    return null;
  }
}

function _alnDfsLeaves(node) {
  if (node.leaf) return [node.gene_id || node.name || ''];
  const out = [];
  for (const c of (node.children || [])) out.push(..._alnDfsLeaves(c));
  return out;
}

function _alnLoadTreeData(id) {
  alnOgMap = {}; alnTreeOrder = []; alnTreeDict = null; alnGeneMeta = {}; alnRangeData = {};
  _alnTreeDataCache = loadDetail(id);
  if (!_alnTreeDataCache) return;
  _alnApplyTreeSource();
}
function _alnApplyTreeSource() {
  alnOgMap = {}; alnTreeOrder = []; alnTreeDict = null; alnGeneMeta = {}; alnRangeData = {};
  if (!_alnTreeDataCache) return;
  const td = _alnTreeDataCache;
  const tree = alnUseGeneRax ? (td.tree || td.prev_tree) : (td.prev_tree || td.tree);
  const ogs  = alnUseGeneRax ? (td.ogs  || td.prev_ogs)  : (td.prev_ogs  || td.ogs);
  const meta = alnUseGeneRax ? (td.meta || td.prev_meta || {}) : (td.prev_meta || td.meta || {});
  if (tree) { alnTreeDict = tree; alnTreeOrder = _alnDfsLeaves(tree); }
  if (meta) alnGeneMeta = meta;
  if (td.range) alnRangeData = td.range;
  if (ogs) {
    for (const [og, genes] of Object.entries(ogs))
      for (const g of genes) alnOgMap[g] = og;
  }
}

function _alnReorderByTree(seqs) {
  if (!alnTreeOrder.length) return seqs;
  const byGene = {};
  for (const s of seqs) byGene[s.name.split(/\s/)[0]] = s;
  const ordered = [];
  for (const gid of alnTreeOrder) if (byGene[gid]) ordered.push(byGene[gid]);
  const inTree = new Set(alnTreeOrder);
  for (const s of seqs) if (!inTree.has(s.name.split(/\s/)[0])) ordered.push(s);
  return ordered;
}

function _alnFitText(ctx, text, maxPx) {
  if (maxPx <= 0) return '';
  if (ctx.measureText(text).width <= maxPx) return text;
  let lo = 0, hi = text.length;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (ctx.measureText(text.slice(0, mid) + '\u2026').width <= maxPx) lo = mid; else hi = mid - 1;
  }
  return lo > 0 ? text.slice(0, lo) + '\u2026' : '';
}

function _alnGetOgColor(og) {
  const ogs = Object.keys(alnOgMap).length ? [...new Set(Object.values(alnOgMap))].sort() : [];
  const idx = ogs.indexOf(og);
  return idx >= 0 ? _alnOgPalette[idx % _alnOgPalette.length] : '#ccc';
}

function _alnColumnDefs() {
  return {
    ref: {
      key: 'ref',
      label: 'Ref',
      width: alnShowRefCol ? alnRefColW : 0,
      minWidth: 40,
      resizerId: 'aln-ref-resize-bar',
    },
    species: {
      key: 'species',
      label: 'Species',
      width: alnShowSpeciesCol ? alnSpeciesColW : 0,
      minWidth: 50,
      resizerId: 'aln-species-resize-bar',
    },
    mrca: {
      key: 'mrca',
      label: 'MRCA',
      width: alnShowMrcaCol ? alnMrcaColW : 0,
      minWidth: 50,
      resizerId: 'aln-mrca-resize-bar',
    },
    range: {
      key: 'range',
      label: 'Domain Architecture',
      width: alnShowRangeCol ? alnRangeColW : 0,
      minWidth: 70,
      resizerId: 'aln-range-resize-bar',
    },
    id: {
      key: 'id',
      label: 'ID',
      width: alnShowSeqId ? alnNameW : 0,
      minWidth: 60,
      resizerId: 'aln-resize-bar',
    },
  };
}

function _alnLayoutMetrics() {
  const treeW = alnTreeDict ? alnTreeW : 0;
  const ogW = Object.keys(alnOgMap).length ? 14 : 0;
  const defs = _alnColumnDefs();
  const metaCols = [];
  const byKey = {};
  let relX = 0;
  for (const key of alnMetaColOrder) {
    const def = defs[key];
    if (!def || def.width <= 0) continue;
    const col = {
      ...def,
      relX0: relX,
      relX1: relX + def.width,
      x0: treeW + ogW + relX,
      x1: treeW + ogW + relX + def.width,
    };
    metaCols.push(col);
    byKey[key] = col;
    relX += def.width;
  }
  return {
    treeW,
    ogW,
    metaCols,
    byKey,
    metaW: relX,
    leftW: treeW + ogW + relX,
  };
}

function _alnMoveMetaColumn(srcKey, dstKey, placeAfter) {
  if (!srcKey || !dstKey || srcKey === dstKey) return;
  const order = alnMetaColOrder.filter(k => k !== srcKey);
  let dstIdx = order.indexOf(dstKey);
  if (dstIdx < 0) return;
  if (placeAfter) dstIdx += 1;
  order.splice(dstIdx, 0, srcKey);
  alnMetaColOrder = order;
  _alnApplyViewerWidths();
  if (alnSeqs.length) renderAlignment();
}

function _alnRangeMeta(gid) {
  if (!gid) return null;
  const proteinEntry = _proteinEntry(gid);
  const currentFamily = (TREE_INDEX.find(rec => rec.id === currentAlnId)?.family) || null;
  const familyTrack = proteinEntry?.tracks?.find(track => track.family === currentFamily)
    || (proteinEntry?.tracks?.length === 1 ? proteinEntry.tracks[0] : null);
  const meta = (GENE_META[gid] || {});
  const lenRaw = proteinEntry?.length ?? meta.length;
  const len = Number.isFinite(+lenRaw) ? +lenRaw : null;
  const span = alnRangeData[gid];
  let start = null, end = null;
  let hits = [];
  if (familyTrack) {
    start = Number.isFinite(+familyTrack.rangeStart) ? +familyTrack.rangeStart : null;
    end = Number.isFinite(+familyTrack.rangeEnd) ? +familyTrack.rangeEnd : null;
    hits = (familyTrack.hits || [])
      .filter(hit => hit.start != null && hit.end != null)
      .map(hit => ({
        start: +hit.start,
        end: +hit.end,
        name: hit.name || '',
        pfam: hit.pfam || '',
      }));
  }
  if (Array.isArray(span) && span.length >= 2) {
    if (start == null) start = Number.isFinite(+span[0]) ? +span[0] : null;
    if (end == null) end = Number.isFinite(+span[1]) ? +span[1] : null;
  }
  if (len != null) {
    if (start != null) start = Math.max(1, Math.min(len, start));
    if (end != null) end = Math.max(1, Math.min(len, end));
  }
  if (start != null && end != null && start > end) {
    const tmp = start;
    start = end;
    end = tmp;
  }
  if (len == null && start == null && end == null) return null;
  return {length: len, start, end, hits, family: familyTrack?.family || currentFamily || ''};
}

function _alnDrawRangeCell(ctx, x0, y0, w, h, gid, rowBg) {
  const meta = _alnRangeMeta(gid);
  if (!meta || w <= 10) return;
  const padX = 8;
  const barW = Math.max(6, w - padX * 2);
  const baseChipH = Math.max(4, Math.min(7, Math.floor(h * 0.30)));
  const chipH = Math.max(3, Math.min(h - 4, Math.round(baseChipH * alnDomainChipScale)));
  const midY = Math.round(y0 + h / 2) + 0.5;
  const chipY = Math.round(midY - chipH / 2) + 0.5;
  const barX = x0 + padX;
  const start = meta.start != null ? meta.start : (meta.hits.length ? Math.min(...meta.hits.map(hit => hit.start)) : 1);
  const end = meta.end != null ? meta.end : (meta.hits.length ? Math.max(...meta.hits.map(hit => hit.end)) : start);
  const span = Math.max(1, end - start + 1);
  const scale = (v) => barX + ((Math.max(start, Math.min(end, v)) - start) / Math.max(1, span - 1)) * barW;
  ctx.fillStyle = rowBg || '#edf1f4';
  ctx.fillRect(x0, y0, w, h);
  ctx.strokeStyle = '#9aa5b1';
  ctx.lineWidth = 1.2;
  ctx.beginPath();
  ctx.moveTo(barX, midY);
  ctx.lineTo(barX + barW, midY);
  ctx.stroke();
  const hits = (meta.hits || []).slice().sort((a, b) => (a.start - b.start) || (a.end - b.end));
  for (const hit of hits) {
    const x1 = scale(hit.start);
    const x2 = scale(hit.end);
    const chipW = Math.max(6, x2 - x1);
    const rr = Math.min(3, chipH / 2);
    ctx.fillStyle = _proteinColor((hit.pfam || hit.name || meta.family || '') + '');
    ctx.beginPath();
    if (typeof ctx.roundRect === 'function') {
      ctx.roundRect(x1, chipY, chipW, chipH, rr);
    } else {
      ctx.rect(x1, chipY, chipW, chipH);
    }
    ctx.fill();
    ctx.strokeStyle = 'rgba(17,17,17,.72)';
    ctx.beginPath();
    if (typeof ctx.roundRect === 'function') {
      ctx.roundRect(x1 + 0.5, chipY + 0.5, Math.max(1, chipW - 1), Math.max(1, chipH - 1), rr);
    } else {
      ctx.rect(x1 + 0.5, chipY + 0.5, Math.max(1, chipW - 1), Math.max(1, chipH - 1));
    }
    ctx.stroke();
  }
}

function _alnColTooltip(row, key) {
  if (!row || row.isCons || row.isCollapsed) return '';
  if (key !== 'range') return '';
  const meta = _alnRangeMeta(row.gid);
  if (!meta) return '';
  if (meta.hits && meta.hits.length && meta.start != null && meta.end != null) {
    return `Domain architecture: ${meta.start}-${meta.end} / ${meta.length || '?'} aa`;
  }
  if (meta.length != null && meta.start != null && meta.end != null) {
    return `Range: ${meta.start}-${meta.end} / ${meta.length} aa`;
  }
  if (meta.length != null) return `Protein length: ${meta.length} aa`;
  if (meta.start != null && meta.end != null) return `Range: ${meta.start}-${meta.end}`;
  return '';
}

function _alnLayoutTree(treeDict, gidY, treeW, useBrLen) {
  if (!treeDict || !Object.keys(gidY).length) return null;

  if (!useBrLen) {
    // cladogram: leaves aligned to right edge, internal nodes by depth
    function maxDepth(n, d) {
      if (n.leaf) return d;
      let mx = d;
      for (const c of (n.children || [])) mx = Math.max(mx, maxDepth(c, d + 1));
      return mx;
    }
    const md = maxDepth(treeDict, 0) || 1;
    const xScale = (treeW - 4) / md;
    const lines = [];
    function lay(node, depth) {
      if (node.leaf) {
        const gid = node.gene_id || node.name || '';
        const y = gidY[gid];
        if (y === undefined) return null;
        return {x: treeW - 2, y};
      }
      const childCoords = [];
      for (const c of (node.children || [])) {
        const cc = lay(c, depth + 1);
        if (cc) childCoords.push(cc);
      }
      if (!childCoords.length) return null;
      const x = 2 + depth * xScale;
      const y = (childCoords[0].y + childCoords[childCoords.length - 1].y) / 2;
      lines.push({x1: x, y1: childCoords[0].y, x2: x, y2: childCoords[childCoords.length - 1].y});
      for (const cc of childCoords) lines.push({x1: x, y1: cc.y, x2: cc.x, y2: cc.y});
      return {x, y};
    }
    lay(treeDict, 0);
    return lines;
  } else {
    // phylogram: x proportional to root-to-tip distance
    function maxRootDist(n, d) {
      const dd = d + (n.dist || 0);
      if (n.leaf) return dd;
      let mx = dd;
      for (const c of (n.children || [])) mx = Math.max(mx, maxRootDist(c, dd));
      return mx;
    }
    const maxD = maxRootDist(treeDict, 0) || 1;
    const xScale = (treeW - 4) / maxD;
    const lines = [];
    function lay(node, rootDist) {
      const dist = rootDist + (node.dist || 0);
      const x = 2 + dist * xScale;
      if (node.leaf) {
        const gid = node.gene_id || node.name || '';
        const y = gidY[gid];
        if (y === undefined) return null;
        return {x, y};
      }
      const childCoords = [];
      for (const c of (node.children || [])) {
        const cc = lay(c, dist);
        if (cc) childCoords.push(cc);
      }
      if (!childCoords.length) return null;
      const y = (childCoords[0].y + childCoords[childCoords.length - 1].y) / 2;
      lines.push({x1: x, y1: childCoords[0].y, x2: x, y2: childCoords[childCoords.length - 1].y});
      for (const cc of childCoords) lines.push({x1: x, y1: cc.y, x2: cc.x, y2: cc.y});
      return {x, y};
    }
    lay(treeDict, 0);
    return lines;
  }
}

// --- species filter ---
function _alnResolveSpFilter() {
  if (!alnSpFilters.length) return null;
  const union = new Set();
  for (const q of alnSpFilters) {
    const lq = q.toLowerCase();
    ALL_SPECIES.forEach(s => { if (s.toLowerCase().includes(lq)) union.add(s); });
  }
  return union.size ? union : null;
}

function _alnFilteredSeqs() {
  let seqs = alnSeqs;
  const spSet = _alnResolveSpFilter();
  if (spSet) {
    seqs = seqs.filter(s => {
      const gid = s.name.split(/\s/)[0];
      const sp = gid.split('_')[0];
      return spSet.has(sp);
    });
  }
  if (alnHiddenOGs.size) {
    seqs = seqs.filter(s => {
      const og = alnOgMap[s.name.split(/\s/)[0]] || '';
      return !alnHiddenOGs.has(og);
    });
  }
  return seqs;
}

function _alnAddSpFilter(query) {
  query = (query || '').trim();
  if (!query || alnSpFilters.includes(query)) return;
  alnSpFilters.push(query);
  const inp = document.getElementById('aln-seq-filter');
  if (inp) inp.value = '';
  _alnRenderSpTags();
  renderAlignment();
}

function _alnRemoveSpFilter(i) {
  alnSpFilters.splice(i, 1);
  _alnRenderSpTags();
  renderAlignment();
}

function _alnClearSpFilters() {
  alnSpFilters = [];
  const inp = document.getElementById('aln-seq-filter');
  if (inp) inp.value = '';
  _alnRenderSpTags();
  renderAlignment();
}

function _alnRenderSpTags() {
  const el = document.getElementById('aln-sp-tags');
  if (!el) return;
  el.innerHTML = '';
  alnSpFilters.forEach((q, i) => {
    const chip = document.createElement('span');
    chip.style.cssText = 'display:inline-flex;align-items:center;gap:2px;padding:1px 6px;background:#e0ecf8;color:#1a56c4;border-radius:3px;font-size:10px;margin-right:3px;';
    chip.textContent = q + ' ';
    const x = document.createElement('span');
    x.textContent = '\u00d7';
    x.style.cssText = 'cursor:pointer;font-weight:bold;margin-left:2px;';
    x.onclick = () => _alnRemoveSpFilter(i);
    chip.appendChild(x);
    el.appendChild(chip);
  });
  if (alnSpFilters.length) {
    const clr = document.createElement('span');
    clr.textContent = '\u2715';
    clr.title = 'Clear all';
    clr.style.cssText = 'cursor:pointer;color:#888;font-size:10px;margin-left:2px;';
    clr.onclick = _alnClearSpFilters;
    el.appendChild(clr);
  }
}

function _alnSpFilterKey(ev) {
  if (ev.key === 'Enter') { _alnAddSpFilter(ev.target.value); ev.preventDefault(); }
}

function selectAlignment(id) {
  if (!HAVE_ALIGNMENTS) return;
  currentAlnId = id;
  alnHiddenOGs.clear(); alnCollapsedOGs.clear();
  document.querySelectorAll('#aln-list .hg-item').forEach(el =>
    el.classList.toggle('selected', el.dataset.id === id));
  const viewer = document.getElementById('aln-viewer');
  const empty  = document.getElementById('aln-empty');
  const info   = document.getElementById('aln-info');
  viewer.style.display  = 'grid';
  empty.style.display   = 'none';
  info.textContent = 'Loading\u2026';
  const el = document.getElementById('alndata-' + id);
  if (!el) { info.textContent = 'No alignment data for ' + id; return; }
  const obj = JSON.parse(el.textContent);
  if (!obj.gz) { info.textContent = 'Alignment not available'; return; }
  let fasta;
  try { fasta = _decompressGz(obj.gz); }
  catch(e) { info.textContent = 'Decompression error: ' + e.message; return; }
  _alnRawFasta = fasta;
  alnSeqs = parseFasta(fasta);
  _alnLoadTreeData(id);
  // show source toggle only when both GeneRax and IQ-TREE data exist
  const srcBtn = document.getElementById('btn-aln-src');
  if (srcBtn) {
    const hasBoth = _alnTreeDataCache && _alnTreeDataCache.tree && _alnTreeDataCache.prev_tree;
    srcBtn.style.display = hasBoth ? '' : 'none';
  }
  alnSeqs = _alnReorderByTree(alnSeqs);
  renderAlignment();
}

function _alnConsensus(seqs) {
  if (!seqs.length) return {seq: '', conservation: []};
  const len = seqs[0].seq.length;
  let cons = '';
  const conservation = [];
  for (let c = 0; c < len; c++) {
    const freq = {};
    let total = 0;
    for (const s of seqs) {
      const aa = (s.seq[c]||'-').toUpperCase();
      if (aa !== '-' && aa !== '.') { freq[aa] = (freq[aa]||0)+1; total++; }
    }
    const entries = Object.entries(freq).sort((a,b)=>b[1]-a[1]);
    if (!entries.length) { cons += '-'; conservation.push(0); continue; }
    const top = entries[0];
    const frac = total > 0 ? top[1] / total : 0;
    conservation.push(frac);
    if (frac >= 0.8) cons += top[0];
    else if (frac >= 0.5) cons += top[0].toLowerCase();
    else cons += '.';
  }
  return {seq: cons, conservation};
}

function renderAlignment() {
  if (!alnSeqs.length) return;
  const seqs = _alnFilteredSeqs();
  const consObj = alnShowCons ? _alnConsensus(seqs) : null;
  const cons = consObj ? consObj.seq : null;
  const conservation = consObj ? consObj.conservation : null;
  const nCols = (seqs[0]?.seq.length) || 0;
  alnCellH = Math.round(alnCellW * 1.6);
  const cW = alnCellW, cH = alnCellH;
  const dims = _alnLayoutMetrics();
  const ogBandW = dims.ogW;
  const treeW = dims.treeW;
  const consBarH = alnShowCons ? 20 : 0;
  const rulerH = 24 + consBarH;
  const totalW = nCols * cW;
  const ogMrca = {};
  for (const s of seqs) {
    const gid = s.name.split(/\s/)[0];
    const og = alnOgMap[gid] || '';
    if (!og) continue;
    if (!ogMrca[og]) ogMrca[og] = new Set();
    const sp = getSpeciesPfx(gid);
    if (sp) ogMrca[og].add(sp);
  }
  Object.keys(ogMrca).forEach(og => { ogMrca[og] = spMRCAName(ogMrca[og]) || ''; });

  _alnRows = [];
  if (cons) _alnRows.push({name:'Consensus', seq: cons, isCons: true, og:'', ogColor:''});
  const _seenColl = new Set();
  for (const s of seqs) {
    const gid = s.name.split(/\s/)[0];
    const og = alnOgMap[gid] || '';
    const ogColor = og ? _alnGetOgColor(og) : '';
    if (og && alnCollapsedOGs.has(og)) {
      if (!_seenColl.has(og)) {
        _seenColl.add(og);
        const ogCount = seqs.filter(ss => (alnOgMap[ss.name.split(/\s/)[0]] || '') === og).length;
        _alnRows.push({name: og, seq: '', isCons: false, isCollapsed: true, og, ogColor, ogCount, species:'', speciesLabel:'', mrca: ogMrca[og] || '', gid:''});
      }
      continue;
    }
    const sp = getSpeciesPfx(gid);
    _alnRows.push({...s, gid, og, ogColor, species: sp || '', speciesLabel: (SPECIES_INFO[sp] || sp || ''), mrca: og ? (ogMrca[og] || '') : ''});
  }
  const rows = _alnRows;
  const totalH = rows.length * cH;

  const dpr = window.devicePixelRatio || 1;

  // --- Ruler + conservation bar ---
  const rulerCanvas = document.getElementById('aln-ruler');
  if (alnShowSeqPanel) {
    rulerCanvas.width  = totalW * dpr;
    rulerCanvas.height = rulerH * dpr;
    rulerCanvas.style.width  = totalW + 'px';
    rulerCanvas.style.height = rulerH + 'px';
    rulerCanvas.style.transform = '';
    const rCtx = rulerCanvas.getContext('2d');
    rCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
    rCtx.clearRect(0, 0, totalW, rulerH);
    rCtx.font = '9px monospace';
    rCtx.textAlign = 'right';
    const step = cW < 6 ? 50 : cW < 9 ? 20 : 10;
    for (let c = 0; c < nCols; c++) {
      if ((c+1) % step === 0) {
        const x = c * cW + cW;
        rCtx.fillStyle = '#aaa'; rCtx.fillRect(x - 0.5, 18, 1, 6);
        rCtx.fillStyle = '#666'; rCtx.fillText(String(c+1), x + 1, 16);
      }
    }
    if (conservation && consBarH > 0) {
      const barTop = 24;
      const barH = consBarH - 2;
      for (let c = 0; c < nCols; c++) {
        const v = conservation[c] || 0;
        const h = v * barH;
        if (h < 0.5) continue;
        const x = c * cW;
        const r = Math.round(180 - v * 130);
        const g = Math.round(180 - v * 100);
        const b = Math.round(200 + v * 55);
        rCtx.fillStyle = `rgb(${r},${g},${b})`;
        rCtx.fillRect(x, barTop + barH - h, cW - (cW > 3 ? 1 : 0), h);
      }
    }
  } else {
    rulerCanvas.width = 1;
    rulerCanvas.height = 1;
    rulerCanvas.style.width = '1px';
    rulerCanvas.style.height = '1px';
  }

  // --- Left panel: tree + OG band + names ---
  const namesCanvas = document.getElementById('aln-names-canvas');
  const leftW = dims.leftW;
  namesCanvas.width  = leftW * dpr;
  namesCanvas.height = totalH * dpr;
  namesCanvas.style.width  = dims.leftW + 'px';
  namesCanvas.style.height = totalH + 'px';
  namesCanvas.style.transform = '';
  const nCtx = namesCanvas.getContext('2d');
  nCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
  nCtx.clearRect(0, 0, dims.leftW, totalH);
  const fontSize = Math.max(8, Math.min(cH - 2, 12));
  nCtx.textBaseline = 'middle';

  // draw cladogram/phylogram
  if (treeW > 0 && alnTreeDict) {
    const gidY = {};
    const _collOgY = {};
    for (let r = 0; r < rows.length; r++) {
      const row = rows[r];
      if (row.isCollapsed) {
        _collOgY[row.og] = r * cH + cH / 2;
      } else if (!row.isCons) {
        const gid = row.name.split(/\s/)[0];
        gidY[gid] = r * cH + cH / 2;
      }
    }
    // assign collapsed-OG genes to their collapsed row Y so branches still render
    for (const s of seqs) {
      const gid = s.name.split(/\s/)[0];
      const og = alnOgMap[gid] || '';
      if (og && alnCollapsedOGs.has(og) && _collOgY[og] !== undefined) gidY[gid] = _collOgY[og];
    }
    const lines = _alnLayoutTree(alnTreeDict, gidY, treeW, alnUseBrLen);
    if (lines) {
      nCtx.strokeStyle = '#333';
      nCtx.lineWidth = 1;
      for (const l of lines) {
        nCtx.beginPath();
        nCtx.moveTo(l.x1, l.y1);
        nCtx.lineTo(l.x2, l.y2);
        nCtx.stroke();
      }
    }
  }

  // draw OG band + names per row
  for (let r = 0; r < rows.length; r++) {
    const row = rows[r];
    const y0 = r * cH;

    if (row.isCollapsed) {
      // collapsed OG: full-width colored bar with label
      nCtx.fillStyle = row.ogColor || '#bbb';
      nCtx.globalAlpha = 0.45;
      nCtx.fillRect(treeW, y0, ogBandW + dims.metaW, cH);
      nCtx.globalAlpha = 1;
      nCtx.strokeStyle = '#555';
      nCtx.lineWidth = 1;
      nCtx.strokeRect(treeW + 0.5, y0 + 0.5, ogBandW + dims.metaW - 1, cH - 1);
      nCtx.font = `bold ${fontSize}px monospace`;
      for (const col of dims.metaCols) {
        if (col.key === 'mrca' && row.mrca) {
          nCtx.fillStyle = '#333';
          nCtx.textAlign = 'right';
          nCtx.fillText(_alnFitText(nCtx, row.mrca, col.width - 8), col.x1 - 4, y0 + cH / 2);
        } else if (col.key === 'id') {
          nCtx.fillStyle = '#222';
          nCtx.textAlign = 'right';
          nCtx.fillText(_alnFitText(nCtx, `\u25B6 ${row.name} (${row.ogCount})`, col.width - 8), col.x1 - 4, y0 + cH / 2);
        }
      }
      continue;
    }

    // OG color band with black border
    if (ogBandW) {
      if (row.ogColor) {
        nCtx.fillStyle = row.ogColor;
        nCtx.globalAlpha = 0.35;
        nCtx.fillRect(treeW, y0, ogBandW, cH);
        nCtx.globalAlpha = 1;
      }
    }
    // name background + text
    const metaX0 = treeW + ogBandW;
    const totalNameW = dims.metaW;
    nCtx.font = `${fontSize}px monospace`;
    nCtx.textAlign = 'right';
    if (row.isCons) {
      nCtx.fillStyle = '#1a6b8a';
      nCtx.fillRect(metaX0, y0, totalNameW, cH);
      nCtx.fillStyle = '#fff';
      const idCol = dims.byKey.id;
      if (idCol) nCtx.fillText(_alnFitText(nCtx, 'Consensus', idCol.width - 8), idCol.x1 - 4, y0 + cH / 2);
    } else {
      if (r % 2 === (cons ? 1 : 0)) {
        nCtx.fillStyle = '#f7f7f7';
        nCtx.fillRect(metaX0, y0, totalNameW, cH);
      }
      const gid = row.gid || row.name.split(/\s/)[0];
      const activeMeta = alnGeneMeta[gid] || {};
      const globalMeta = GENE_META[gid] || {};
      const isRefSeq = !!globalMeta.is_reference_gene || !!REFNAME_MAP[gid];
      const refName = activeMeta.ref_ortholog || globalMeta.ref_ortholog || REFNAME_MAP[gid] || '';
      for (const col of dims.metaCols) {
        if (col.key === 'ref' && refName) {
          nCtx.fillStyle = isRefSeq ? '#cc0000' : '#777';
          nCtx.fillText(_alnFitText(nCtx, refName, col.width - 8), col.x1 - 4, y0 + cH / 2);
        } else if (col.key === 'species' && row.speciesLabel) {
          nCtx.fillStyle = row.species ? spColor(row.species) : '#555';
          nCtx.fillText(_alnFitText(nCtx, row.speciesLabel, col.width - 8), col.x1 - 4, y0 + cH / 2);
        } else if (col.key === 'mrca' && row.mrca) {
          nCtx.fillStyle = '#555';
          nCtx.fillText(_alnFitText(nCtx, row.mrca, col.width - 8), col.x1 - 4, y0 + cH / 2);
        } else if (col.key === 'range') {
          _alnDrawRangeCell(nCtx, col.x0, y0, col.width, cH, gid, r % 2 === (cons ? 1 : 0) ? '#eceff2' : '#edf1f4');
        } else if (col.key === 'id') {
          nCtx.fillStyle = isRefSeq ? '#cc0000' : '#333';
          nCtx.fillText(_alnFitText(nCtx, row.name, col.width - 8), col.x1 - 4, y0 + cH / 2);
        }
      }
    }
  }

  for (let i = 0; i < dims.metaCols.length - 1; i++) {
    const sepX = dims.metaCols[i].x1 + 0.5;
    nCtx.strokeStyle = '#ddd'; nCtx.lineWidth = 1;
    nCtx.beginPath(); nCtx.moveTo(sepX, 0); nCtx.lineTo(sepX, totalH); nCtx.stroke();
  }

  // OG group borders: draw merged boxes (left/right/top/bottom per OG run)
  if (ogBandW) {
    nCtx.strokeStyle = '#000';
    nCtx.lineWidth = 1;
    let runStart = -1, runOg = '';
    for (let r = 0; r <= rows.length; r++) {
      const og = r < rows.length ? rows[r].og : '';
      if (og !== runOg || r === rows.length) {
        if (runOg && runStart >= 0) {
          const y1 = runStart * cH;
          const y2 = r * cH;
          nCtx.strokeRect(treeW + 0.5, y1 + 0.5, ogBandW - 1, y2 - y1 - 1);
        }
        runStart = r;
        runOg = og;
      }
    }
  }

  // OG separator lines on names canvas
  nCtx.strokeStyle = '#000';
  nCtx.lineWidth = 1;
  for (let r = 1; r < rows.length; r++) {
    if ((rows[r].og || rows[r-1].og) && rows[r].og !== rows[r-1].og) {
      const y = r * cH + 0.5;
      nCtx.beginPath(); nCtx.moveTo(0, y); nCtx.lineTo(dims.leftW, y); nCtx.stroke();
    }
  }

  // --- Sequences ---
  const seqCanvas = document.getElementById('aln-seq-canvas');
  const seqWrap   = document.getElementById('aln-seq-wrap');
  if (alnShowSeqPanel) {
    seqCanvas.width  = totalW * dpr;
    seqCanvas.height = totalH * dpr;
    seqCanvas.style.width  = totalW + 'px';
    seqCanvas.style.height = totalH + 'px';
    const sCtx = seqCanvas.getContext('2d');
    sCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
    sCtx.clearRect(0, 0, totalW, totalH);
    const showLetters = cW >= 9;
    if (showLetters) {
      sCtx.font = `bold ${Math.max(8, cW - 1)}px monospace`;
      sCtx.textAlign = 'center';
      sCtx.textBaseline = 'middle';
    }
    for (let r = 0; r < rows.length; r++) {
      const row = rows[r];
      const y0 = r * cH;
      const isConsRow = row.isCons;

      if (row.isCollapsed) {
        sCtx.fillStyle = row.ogColor || '#bbb';
        sCtx.globalAlpha = 0.18;
        sCtx.fillRect(0, y0, totalW, cH);
        sCtx.globalAlpha = 1;
        sCtx.strokeStyle = '#aaa';
        sCtx.lineWidth = 0.5;
        sCtx.beginPath(); sCtx.moveTo(0, y0); sCtx.lineTo(totalW, y0); sCtx.stroke();
        sCtx.beginPath(); sCtx.moveTo(0, y0 + cH - 0.5); sCtx.lineTo(totalW, y0 + cH - 0.5); sCtx.stroke();
        sCtx.font = `bold ${Math.max(10, cH - 2)}px sans-serif`;
        sCtx.fillStyle = '#444';
        sCtx.textAlign = 'center';
        sCtx.textBaseline = 'middle';
        sCtx.fillText(`${row.name} \u00b7 ${row.ogCount} sequences`, totalW / 2, y0 + cH / 2);
        continue;
      }

      for (let c = 0; c < nCols; c++) {
        const aa = (row.seq[c] || '-').toUpperCase();
        const x = c * cW;
        let bg = AA_COLORS[aa];
        if (!bg && aa !== '-' && aa !== '.') bg = '#e0e0e0';
        if (isConsRow) {
          sCtx.fillStyle = bg || '#ddf0ff'; sCtx.fillRect(x, y0, cW, cH);
        } else {
          if (r % 2 === (cons ? 1 : 0)) { sCtx.fillStyle = '#f7f7f7'; sCtx.fillRect(x, y0, cW, cH); }
          if (bg) { sCtx.fillStyle = bg; sCtx.fillRect(x, y0, cW, cH); }
        }
        if (showLetters && aa !== '-' && aa !== '.') {
          sCtx.fillStyle = isConsRow ? '#000' : '#222';
          sCtx.fillText(aa, x + cW / 2, y0 + cH / 2);
        }
      }
    }
    sCtx.strokeStyle = '#000';
    sCtx.lineWidth = 1;
    for (let r = 1; r < rows.length; r++) {
      if ((rows[r].og || rows[r-1].og) && rows[r].og !== rows[r-1].og) {
        const y = r * cH + 0.5;
        sCtx.beginPath(); sCtx.moveTo(0, y); sCtx.lineTo(totalW, y); sCtx.stroke();
      }
    }
  } else {
    seqCanvas.width = 1;
    seqCanvas.height = 1;
    seqCanvas.style.width = '1px';
    seqCanvas.style.height = '1px';
  }

  const info = document.getElementById('aln-info');
  const ogInfo = Object.keys(alnOgMap).length ? ` \u00b7 ${[...new Set(Object.values(alnOgMap))].length} OGs` : '';
  info.textContent = `${seqs.length} seqs \u00d7 ${nCols} cols${ogInfo}${alnShowSeqPanel ? '' : ' \u00b7 alignment hidden'}`;

  const viewer = document.getElementById('aln-viewer');
  _alnApplyViewerWidths();
  viewer.style.gridTemplateRows = `${alnShowSeqPanel ? rulerH : 0}px 1fr`;

  if (alnShowSeqPanel) {
    seqWrap.scrollLeft = 0; seqWrap.scrollTop = 0;
    syncAlnScroll(seqWrap);
  }
  _alnUpdateCollapseAllBtn();
}

function syncAlnScroll(seqWrap) {
  if (!alnShowSeqPanel) return;
  const ruler = document.getElementById('aln-ruler');
  const rulerWrap = document.getElementById('aln-ruler-wrap');
  const namesWrap = document.getElementById('aln-names-wrap');
  if (ruler) ruler.style.transform = `translateX(-${seqWrap.scrollLeft}px)`;
  if (_alnSyncingScroll) return;
  _alnSyncingScroll = true;
  if (rulerWrap) rulerWrap.scrollLeft = seqWrap.scrollLeft;
  if (namesWrap) namesWrap.scrollTop = seqWrap.scrollTop;
  _alnSyncingScroll = false;
}

function _alnSyncFromNamesWrap(namesWrap) {
  if (!alnShowSeqPanel || _alnSyncingScroll) return;
  const seqWrap = document.getElementById('aln-seq-wrap');
  if (!seqWrap || !namesWrap) return;
  _alnSyncingScroll = true;
  seqWrap.scrollTop = namesWrap.scrollTop;
  const rulerWrap = document.getElementById('aln-ruler-wrap');
  if (rulerWrap) rulerWrap.scrollLeft = seqWrap.scrollLeft;
  _alnSyncingScroll = false;
}

function _alnHandleNameHover(ev) {
  const canvas = document.getElementById('aln-names-canvas');
  if (!canvas || !_alnRows.length) { _alnHideTip(); return; }
  const rect = canvas.getBoundingClientRect();
  const y = ev.clientY - rect.top;
  const cH = alnCellH;
  const dims = _alnLayoutMetrics();
  const x = ev.clientX - rect.left;
  const r = Math.floor(y / cH);
  if (r < 0 || r >= _alnRows.length) { _alnHideTip(); return; }
  const row = _alnRows[r];
  if (row.isCollapsed) {
    _alnShowTip(ev, `${row.og} \u00b7 ${row.ogCount} seqs (collapsed)`);
  } else if (x >= dims.treeW && x < dims.treeW + dims.ogW && row.og) {
    _alnShowTip(ev, row.og);
  } else {
    let shown = false;
    for (const col of dims.metaCols) {
      if (x < col.x0 || x >= col.x1) continue;
      if (!row.isCons && col.key === 'id') {
        _alnShowTip(ev, `Click to copy: ${row.name.split(/\s/)[0]}`);
        shown = true;
      } else if (!row.isCons && !row.isCollapsed && col.key === 'range') {
        const popup = _proteinPopupHtml(row.gid, {
          family: (TREE_INDEX.find(rec => rec.id === currentAlnId)?.family) || ''
        });
        if (popup) {
          _alnShowTip(ev, popup);
          shown = true;
        } else {
          const tip = _alnColTooltip(row, col.key);
          if (tip) {
            _alnShowTip(ev, tip);
            shown = true;
          }
        }
      } else {
        const tip = _alnColTooltip(row, col.key);
        if (tip) {
          _alnShowTip(ev, tip);
          shown = true;
        }
      }
      break;
    }
    if (!shown) _alnHideTip();
  }
}

function _alnShowTip(ev, text) {
  let tip = document.getElementById('aln-tooltip');
  if (!tip) {
    tip = document.createElement('div');
    tip.id = 'aln-tooltip';
    tip.style.cssText = 'position:fixed;padding:6px 8px;background:rgba(255,255,255,.97);color:#222;font-size:11px;border-radius:5px;border:1px solid #bbb;box-shadow:0 2px 8px rgba(0,0,0,.15);pointer-events:none;z-index:9999;white-space:normal;max-width:min(820px,calc(100vw - 32px));';
    document.body.appendChild(tip);
  }
  if (typeof text === 'string' && text.indexOf('<') >= 0) tip.innerHTML = text;
  else tip.textContent = text;
  tip.style.display = 'block';
  const left = Math.min(ev.clientX + 12, window.innerWidth - tip.offsetWidth - 8);
  const top = Math.min(Math.max(8, ev.clientY - 8), window.innerHeight - tip.offsetHeight - 8);
  tip.style.left = Math.max(8, left) + 'px';
  tip.style.top  = Math.max(8, top) + 'px';
}

function _alnHideTip() {
  const tip = document.getElementById('aln-tooltip');
  if (tip) tip.style.display = 'none';
}

function _copyTextToClipboard(txt) {
  if (navigator.clipboard) return navigator.clipboard.writeText(txt).catch(()=>{});
  return new Promise(resolve=>{
    const ta=document.createElement('textarea');
    ta.value=txt;
    ta.style.position='fixed';
    ta.style.opacity='0';
    document.body.appendChild(ta);
    ta.select();
    try { document.execCommand('copy'); } catch(_) {}
    document.body.removeChild(ta);
    resolve();
  });
}

function _alnHandleNameClick(ev) {
  const canvas = document.getElementById('aln-names-canvas');
  if (!canvas || !_alnRows.length) return;
  const rect = canvas.getBoundingClientRect();
  const y = ev.clientY - rect.top;
  const cH = alnCellH;
  const dims = _alnLayoutMetrics();
  const x = ev.clientX - rect.left;
  const r = Math.floor(y / cH);
  if (r < 0 || r >= _alnRows.length) return;
  const row = _alnRows[r];
  const onOgBand = x >= dims.treeW && x < dims.treeW + dims.ogW;
  const idCol = dims.byKey.id;
  if (!row.isCons && !row.isCollapsed && idCol && x >= idCol.x0 && x < idCol.x1) {
    const gid=row.gid || row.name.split(/\s/)[0];
    _copyTextToClipboard(gid).then(()=>_alnShowTip(ev, `Copied: ${gid}`));
    setTimeout(_alnHideTip, 900);
    return;
  }
  if ((onOgBand || row.isCollapsed) && row.og) _alnShowOgPopup(ev, row.og);
}

function _alnHandleSeqClick(ev) {
  if (!_alnRows.length) return;
  const seqWrap = document.getElementById('aln-seq-wrap');
  const canvas = document.getElementById('aln-seq-canvas');
  if (!seqWrap || !canvas) return;
  const rect = canvas.getBoundingClientRect();
  const y = ev.clientY - rect.top;
  const cH = alnCellH;
  const r = Math.floor(y / cH);
  if (r < 0 || r >= _alnRows.length) return;
  const row = _alnRows[r];
  if (row.isCollapsed && row.og) _alnShowOgPopup(ev, row.og);
}

function _alnShowOgPopup(ev, og) {
  let popup = document.getElementById('aln-og-popup');
  if (!popup) {
    popup = document.createElement('div');
    popup.id = 'aln-og-popup';
    popup.style.cssText = 'position:fixed;background:#fff;border:1px solid #bbb;border-radius:6px;padding:8px 10px;font-size:11px;box-shadow:0 3px 10px rgba(0,0,0,.2);z-index:9999;min-width:170px;';
    document.body.appendChild(popup);
    document.addEventListener('mousedown', e => {
      if (!document.getElementById('aln-og-popup')?.contains(e.target))
        _alnCloseOgPopup();
    }, true);
  }
  const isCollapsed = alnCollapsedOGs.has(og);
  const isHidden    = alnHiddenOGs.has(og);
  const hasState    = alnHiddenOGs.size > 0 || alnCollapsedOGs.size > 0;
  const esc = og.replace(/\\/g,'\\\\').replace(/'/g,"\\'");
  const btnStyle = 'display:block;width:100%;margin-top:5px;padding:4px 8px;font-size:11px;border:1px solid #ccc;border-radius:4px;background:#f8f8f8;cursor:pointer;text-align:left;';
  let html = `<div style="font-weight:700;color:#333;margin-bottom:4px;font-size:12px">${og}</div>`;
  html += `<button style="${btnStyle}" onclick="_alnOgToggleCollapse('${esc}')">` +
          (isCollapsed ? '\u25BC Expand' : '\u25B6 Collapse') + '</button>';
  html += `<button style="${btnStyle}" onclick="_alnOgIsolate('${esc}')">&#9654; Show only this OG</button>`;
  if (isHidden)
    html += `<button style="${btnStyle}" onclick="_alnOgShow('${esc}')">&#128065; Show this OG</button>`;
  else
    html += `<button style="${btnStyle}" onclick="_alnOgHide('${esc}')">&#10006; Hide this OG</button>`;
  if (hasState)
    html += `<button style="${btnStyle}color:#2980b9;border-color:#2980b9;" onclick="_alnOgShowAll()">&#8635; Show all OGs</button>`;
  popup.innerHTML = html;
  popup.style.display = 'block';
  popup.style.left = Math.min(ev.clientX + 8, window.innerWidth - 185) + 'px';
  popup.style.top  = (ev.clientY - 8) + 'px';
}

function _alnCloseOgPopup() {
  const p = document.getElementById('aln-og-popup');
  if (p) p.style.display = 'none';
}

function _alnOgToggleCollapse(og) {
  _alnCloseOgPopup();
  if (alnCollapsedOGs.has(og)) alnCollapsedOGs.delete(og);
  else { alnCollapsedOGs.add(og); alnHiddenOGs.delete(og); }
  renderAlignment();
}

function _alnOgIsolate(og) {
  _alnCloseOgPopup();
  alnHiddenOGs.clear(); alnCollapsedOGs.clear();
  const allOGs = [...new Set(Object.values(alnOgMap))];
  allOGs.forEach(o => { if (o !== og) alnHiddenOGs.add(o); });
  renderAlignment();
}

function _alnOgHide(og) {
  _alnCloseOgPopup();
  alnHiddenOGs.add(og); alnCollapsedOGs.delete(og);
  renderAlignment();
}

function _alnOgShow(og) {
  _alnCloseOgPopup();
  alnHiddenOGs.delete(og);
  renderAlignment();
}

function _alnOgShowAll() {
  _alnCloseOgPopup();
  alnHiddenOGs.clear(); alnCollapsedOGs.clear();
  renderAlignment();
}

function _alnRenderCorner(){
  const corner=document.getElementById('aln-corner');
  if(!corner) return;
  const dims=_alnLayoutMetrics();
  corner.style.width=dims.leftW+'px';
  corner.style.minWidth=dims.leftW+'px';
  const labels=[];
  if(dims.treeW>0) labels.push({left:0,width:dims.treeW,label:'Gene Tree'});
  if(dims.ogW>0) labels.push({left:dims.treeW,width:dims.ogW,label:'OG'});
  dims.metaCols.forEach(col=>labels.push({left:col.x0,width:col.width,label:col.label,key:col.key}));
  corner.innerHTML = '';
  for (const rec of labels) {
    const el = document.createElement('div');
    el.className = `aln-col-header${rec.key ? ' drag' : ''}`;
    el.style.left = `${rec.left}px`;
    el.style.width = `${Math.max(0, rec.width)}px`;
    el.textContent = rec.label;
    el.title = rec.label;
    if (rec.key) {
      el.dataset.alnCol = rec.key;
      el.draggable = true;
    }
    corner.appendChild(el);
  }
}

function _alnProxyWheel(ev) {
  const seqWrap = document.getElementById('aln-seq-wrap');
  if (!alnShowSeqPanel || !seqWrap) return;
  const dx = Number(ev.deltaX) || 0;
  const dy = Number(ev.deltaY) || 0;
  if (!dx && !dy) return;
  if (Math.abs(dy) >= Math.abs(dx)) seqWrap.scrollTop += dy;
  else seqWrap.scrollLeft += dx;
  syncAlnScroll(seqWrap);
  ev.preventDefault();
}

(function(){
  const bind = id => {
    const el = document.getElementById(id);
    if (!el) return;
    el.addEventListener('wheel', _alnProxyWheel, {passive:false});
  };
  bind('aln-corner');
  bind('aln-ruler-wrap');
  bind('aln-names-wrap');
})();

function _alnApplyViewerWidths(){
  const viewer = document.getElementById('aln-viewer');
  if(!viewer) return;
  const dims = _alnLayoutMetrics();
  const seqVisible = !!alnShowSeqPanel;
  viewer.style.gridTemplateColumns = `${dims.leftW}px ${seqVisible ? '5px' : '0px'} ${seqVisible ? '1fr' : '0px'}`;
  _alnRenderCorner();
  const rulerWrap = document.getElementById('aln-ruler-wrap');
  const seqWrap = document.getElementById('aln-seq-wrap');
  const mainBar = document.getElementById('aln-resize-bar');
  if (rulerWrap) rulerWrap.style.display = seqVisible ? '' : 'none';
  if (seqWrap) seqWrap.style.display = seqVisible ? '' : 'none';
  if (mainBar) mainBar.style.display = seqVisible ? '' : 'none';
  const treeBar = document.getElementById('aln-tree-resize-bar');
  const refBar = document.getElementById('aln-ref-resize-bar');
  const speciesBar = document.getElementById('aln-species-resize-bar');
  const mrcaBar = document.getElementById('aln-mrca-resize-bar');
  const rangeBar = document.getElementById('aln-range-resize-bar');
  if (treeBar) {
    treeBar.style.display = alnTreeDict ? 'block' : 'none';
    treeBar.style.left = dims.treeW + 'px';
  }
  if (refBar) {
    refBar.style.display = dims.byKey.ref ? 'block' : 'none';
    refBar.style.left = dims.byKey.ref ? dims.byKey.ref.x1 + 'px' : '0px';
  }
  if (speciesBar) {
    speciesBar.style.display = dims.byKey.species ? 'block' : 'none';
    speciesBar.style.left = dims.byKey.species ? dims.byKey.species.x1 + 'px' : '0px';
  }
  if (mrcaBar) {
    mrcaBar.style.display = dims.byKey.mrca ? 'block' : 'none';
    mrcaBar.style.left = dims.byKey.mrca ? dims.byKey.mrca.x1 + 'px' : '0px';
  }
  if (rangeBar) {
    rangeBar.style.display = dims.byKey.range ? 'block' : 'none';
    rangeBar.style.left = dims.byKey.range ? dims.byKey.range.x1 + 'px' : '0px';
  }
}

(function(){
  let _drag = null, _startX = 0, _startW = 0;
  document.addEventListener('mousedown', ev => {
    const id = ev.target.id;
    if (id !== 'aln-resize-bar' && id !== 'aln-tree-resize-bar' && id !== 'aln-ref-resize-bar' && id !== 'aln-species-resize-bar' && id !== 'aln-mrca-resize-bar' && id !== 'aln-range-resize-bar') return;
    _drag = id; _startX = ev.clientX;
    _startW =
      id === 'aln-resize-bar' ? alnNameW :
      id === 'aln-tree-resize-bar' ? alnTreeW :
      id === 'aln-ref-resize-bar' ? alnRefColW :
      id === 'aln-species-resize-bar' ? alnSpeciesColW :
      id === 'aln-mrca-resize-bar' ? alnMrcaColW :
      alnRangeColW;
    ev.target.classList.add('dragging');
    document.body.style.cursor = 'col-resize';
    ev.preventDefault();
  });
  document.addEventListener('mousemove', ev => {
    if (!_drag) return;
    const delta = ev.clientX - _startX;
    if (_drag === 'aln-resize-bar') {
      alnNameW = Math.max(60, _startW + delta);
    } else if (_drag === 'aln-tree-resize-bar') {
      alnTreeW = Math.max(0, _startW + delta);
    } else if (_drag === 'aln-ref-resize-bar') {
      alnRefColW = Math.max(40, _startW + delta);
    } else if (_drag === 'aln-species-resize-bar') {
      alnSpeciesColW = Math.max(50, _startW + delta);
    } else if (_drag === 'aln-mrca-resize-bar') {
      alnMrcaColW = Math.max(50, _startW + delta);
    } else {
      alnRangeColW = Math.max(70, _startW + delta);
    }
    _alnApplyViewerWidths();
  });
  document.addEventListener('mouseup', ev => {
    if (!_drag) return;
    const bar = document.getElementById(_drag);
    if (bar) bar.classList.remove('dragging');
    _drag = null;
    document.body.style.cursor = '';
    if (alnSeqs.length) renderAlignment();
  });
})();

(function(){
  const corner = document.getElementById('aln-corner');
  if (!corner) return;
  let dragKey = null;
  function headerTarget(node){
    return node && typeof node.closest === 'function' ? node.closest('[data-aln-col]') : null;
  }
  corner.addEventListener('dragstart', ev => {
    const target = headerTarget(ev.target);
    if (!target) return;
    dragKey = target.dataset.alnCol || null;
    if (ev.dataTransfer) {
      ev.dataTransfer.effectAllowed = 'move';
      ev.dataTransfer.setData('text/plain', dragKey || '');
    }
    target.style.opacity = '0.55';
  });
  corner.addEventListener('dragover', ev => {
    const target = headerTarget(ev.target);
    if (!dragKey || !target || target.dataset.alnCol === dragKey) return;
    ev.preventDefault();
    if (ev.dataTransfer) ev.dataTransfer.dropEffect = 'move';
  });
  corner.addEventListener('drop', ev => {
    const target = headerTarget(ev.target);
    if (!dragKey || !target) return;
    const dstKey = target.dataset.alnCol || '';
    if (!dstKey || dstKey === dragKey) return;
    ev.preventDefault();
    const rect = target.getBoundingClientRect();
    const placeAfter = ev.clientX > rect.left + rect.width / 2;
    _alnMoveMetaColumn(dragKey, dstKey, placeAfter);
  });
  corner.addEventListener('dragend', ev => {
    const target = headerTarget(ev.target);
    if (target) target.style.opacity = '';
    dragKey = null;
  });
})();

function alnToggleSeqId() {
  alnShowSeqId = !alnShowSeqId;
  const btn = document.getElementById('btn-aln-seqid');
  if (btn) { btn.style.background = alnShowSeqId?'#e8f0fe':'#fff'; btn.style.color = alnShowSeqId?'#1a56c4':'#555'; btn.style.borderColor = alnShowSeqId?'#1a56c4':'#888'; btn.style.fontWeight = alnShowSeqId?'600':''; }
  renderAlignment();
}
function alnToggleRefCol() {
  alnShowRefCol = !alnShowRefCol;
  const btn = document.getElementById('btn-aln-ref');
  if (btn) { btn.style.background = alnShowRefCol?'#e8f0fe':'#fff'; btn.style.color = alnShowRefCol?'#1a56c4':'#555'; btn.style.borderColor = alnShowRefCol?'#1a56c4':'#888'; btn.style.fontWeight = alnShowRefCol?'600':''; }
  renderAlignment();
}
function alnToggleSpeciesCol() {
  alnShowSpeciesCol = !alnShowSpeciesCol;
  const btn = document.getElementById('btn-aln-species');
  if (btn) { btn.style.background = alnShowSpeciesCol?'#e8f0fe':'#fff'; btn.style.color = alnShowSpeciesCol?'#1a56c4':'#555'; btn.style.borderColor = alnShowSpeciesCol?'#1a56c4':'#888'; btn.style.fontWeight = alnShowSpeciesCol?'600':''; }
  renderAlignment();
}
function alnToggleMrcaCol() {
  alnShowMrcaCol = !alnShowMrcaCol;
  const btn = document.getElementById('btn-aln-mrca');
  if (btn) { btn.style.background = alnShowMrcaCol?'#e8f0fe':'#fff'; btn.style.color = alnShowMrcaCol?'#1a56c4':'#555'; btn.style.borderColor = alnShowMrcaCol?'#1a56c4':'#888'; btn.style.fontWeight = alnShowMrcaCol?'600':''; }
  renderAlignment();
}
function alnToggleRangeCol() {
  alnShowRangeCol = !alnShowRangeCol;
  const btn = document.getElementById('btn-aln-range');
  if (btn) { btn.style.background = alnShowRangeCol?'#e8f0fe':'#fff'; btn.style.color = alnShowRangeCol?'#1a56c4':'#555'; btn.style.borderColor = alnShowRangeCol?'#1a56c4':'#888'; btn.style.fontWeight = alnShowRangeCol?'600':''; }
  renderAlignment();
}
function alnSetDomainChipScale(value) {
  const v = Math.max(0.6, Math.min(2.2, Number(value) || 1));
  alnDomainChipScale = v;
  const out = document.getElementById('aln-domain-chip-val');
  if (out) out.textContent = v.toFixed(1);
  if (alnSeqs.length) renderAlignment();
}
function alnToggleSeqPanel() {
  alnShowSeqPanel = !alnShowSeqPanel;
  const btn = document.getElementById('btn-aln-seqpanel');
  if (btn) {
    btn.style.background = alnShowSeqPanel ? '#e8f0fe' : '#fff';
    btn.style.color = alnShowSeqPanel ? '#1a56c4' : '#555';
    btn.style.borderColor = alnShowSeqPanel ? '#1a56c4' : '#888';
    btn.style.fontWeight = alnShowSeqPanel ? '600' : '';
  }
  renderAlignment();
}

function _alnUpdateCollapseAllBtn() {
  const btn = document.getElementById('btn-aln-collapse-ogs');
  if (!btn) return;
  const allOGs = [...new Set(Object.values(alnOgMap))];
  const active = allOGs.length > 0 && allOGs.every(og => alnCollapsedOGs.has(og));
  btn.style.background  = active ? '#e8f0fe' : '#fff';
  btn.style.color       = active ? '#1a56c4' : '#555';
  btn.style.borderColor = active ? '#1a56c4' : '#888';
  btn.style.fontWeight  = active ? '600' : '';
}

function alnToggleCollapseAllOGs() {
  const allOGs = [...new Set(Object.values(alnOgMap))];
  if (!allOGs.length) return;
  const allCollapsed = allOGs.every(og => alnCollapsedOGs.has(og));
  alnCollapsedOGs.clear();
  if (!allCollapsed) { allOGs.forEach(og => alnCollapsedOGs.add(og)); alnHiddenOGs.clear(); }
  _alnUpdateCollapseAllBtn();
  renderAlignment();
}

function alnToggleConsensus() {
  alnShowCons = !alnShowCons;
  const btn = document.getElementById('btn-aln-cons');
  if (btn) {
    btn.style.background  = alnShowCons ? '#e8f0fe' : '#fff';
    btn.style.color       = alnShowCons ? '#1a56c4' : '#555';
    btn.style.borderColor = alnShowCons ? '#1a56c4' : '#888';
    btn.style.fontWeight  = alnShowCons ? '600' : '';
  }
  renderAlignment();
}

function alnToggleBrLen() {
  alnUseBrLen = !alnUseBrLen;
  const btn = document.getElementById('btn-aln-brlen');
  if (btn) {
    btn.style.background  = alnUseBrLen ? '#e8f0fe' : '#fff';
    btn.style.color       = alnUseBrLen ? '#1a56c4' : '#555';
    btn.style.borderColor = alnUseBrLen ? '#1a56c4' : '#888';
    btn.style.fontWeight  = alnUseBrLen ? '600' : '';
  }
  renderAlignment();
}

function alnDownload() {
  if (!_alnRawFasta || !currentAlnId) return;
  const blob = new Blob([_alnRawFasta], {type:'text/plain'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = currentAlnId + '.aln.fasta';
  a.click();
  URL.revokeObjectURL(a.href);
}

function alnToggleSource() {
  alnUseGeneRax = !alnUseGeneRax;
  const btn = document.getElementById('btn-aln-src');
  if (btn) {
    btn.textContent = alnUseGeneRax ? 'GeneRax' : 'IQ-TREE';
    btn.style.background  = alnUseGeneRax ? '#e8f0fe' : '#fff';
    btn.style.color       = alnUseGeneRax ? '#1a56c4' : '#555';
    btn.style.borderColor = alnUseGeneRax ? '#1a56c4' : '#888';
    btn.style.fontWeight  = alnUseGeneRax ? '600' : '';
  }
  _alnApplyTreeSource();
  if (alnSeqs.length) {
    alnSeqs = _alnReorderByTree(alnSeqs);
    renderAlignment();
  }
}

function alnDownloadPng() {
  if (!currentAlnId) return;
  const ruler = document.getElementById('aln-ruler');
  const names = document.getElementById('aln-names-canvas');
  const seqs  = document.getElementById('aln-seq-canvas');
  if (!ruler || !names || !seqs) return;
  const leftW = names.width, topH = ruler.height;
  const totalW = leftW + seqs.width, totalH = topH + seqs.height;
  const out = document.createElement('canvas');
  out.width = totalW; out.height = totalH;
  const ctx = out.getContext('2d');
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, totalW, totalH);
  ctx.drawImage(ruler, leftW, 0);
  ctx.drawImage(names, 0, topH);
  ctx.drawImage(seqs, leftW, topH);
  out.toBlob(blob => {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = currentAlnId + '.alignment.png';
    a.click();
    URL.revokeObjectURL(a.href);
  });
}

function renderAlnSidebar(query) {
  const q = (query || '').toLowerCase().trim();
  const list = document.getElementById('aln-list');
  if (!list) return;
  // group by family (reuse TREE_INDEX which has family/id fields)
  const byFam = {};
  for (const rec of TREE_INDEX) {
    const fam = rec.family || '?';
    if (!byFam[fam]) byFam[fam] = [];
    byFam[fam].push(rec);
  }
  let html = '';
  for (const [fam, recs] of Object.entries(byFam).sort((a,b)=>a[0].localeCompare(b[0]))) {
    const filtered = q ? recs.filter(r => r.id.toLowerCase().includes(q) || fam.toLowerCase().includes(q) || (r.og_names||[]).some(o => o.toLowerCase().includes(q))) : recs;
    if (!filtered.length) continue;
    filtered.sort((a,b)=>{ const na=parseInt((a.id).match(/\d+/)?.[0]||'0',10), nb=parseInt((b.id).match(/\d+/)?.[0]||'0',10); return na-nb; });
    const open = !q && Object.keys(byFam).length > 1 ? '' : ' open';
    html += `<div class="fam-header${open ? ' open' : ''}" onclick="this.classList.toggle('open');this.nextElementSibling.classList.toggle('open')">` +
            `<span class="fam-arrow">&#9658;</span>${fam} <span style="color:#999;font-weight:400">(${filtered.length})</span></div>` +
            `<div class="fam-body${open ? ' open' : ''}">`;
    for (const rec of filtered) {
      const sel = rec.id === currentAlnId ? ' selected' : '';
      const hasAln = !!document.getElementById('alndata-' + rec.id);
      const badge = hasAln ? '' : '<span style="font-size:9px;color:#aaa;margin-left:4px">no aln</span>';
      html += `<div class="hg-item${sel}" data-id="${rec.id}" onclick="selectAlignment('${rec.id}')">` +
              `<div class="hg-name">${rec.id}${badge}</div>` +
              `<div class="hg-meta">${rec.n_ogs||0} OGs · ${rec.species?.length||0} sps</div></div>`;
    }
    html += '</div>';
  }
  list.innerHTML = html || '<div style="padding:12px;color:#aaa;font-size:11px">No results</div>';
  // open first family if only one and query active
  if (q) {
    list.querySelectorAll('.fam-header').forEach(h => { h.classList.add('open'); h.nextElementSibling.classList.add('open'); });
  }
}

// ═══════════════════════════════════════════════════════════════════════════════
// FAMILIES TABLE
// ═══════════════════════════════════════════════════════════════════════════════

let _famSortCol = "family", _famSortAsc = true;
let _famFilter  = "";
let _famRendered = false;
let _famExpandedCats = new Set();  // category keys currently expanded

// pre-index FAMILY_DATA by family name for quick species-union lookup
const _famDataByName = {};
FAMILY_DATA.forEach(d=>{ _famDataByName[d.family]=d; });

function famToggleCat(cat){
  if(_famExpandedCats.has(cat)) _famExpandedCats.delete(cat);
  else _famExpandedCats.add(cat);
  drawFamilyTable();
}

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

  countEl.textContent = rows.length + " / " + FAMILY_INFO.length + " families";

  const _totalSp = ALL_SPECIES.length;
  function speciesFrac(n){
    if(!_totalSp) return `<span style="color:#555">${n}</span>`;
    const frac = n/_totalSp;
    const col = frac>=0.8?"#27ae60":frac>=0.4?"#e67e22":"#c0392b";
    return `<span style="color:${col};font-weight:600">${n}</span><span style="color:#aaa">/${_totalSp}</span>`;
  }
  function treeFrac(n, tot){
    if(!tot) return `<span style="color:#bbb">—</span>`;
    const frac = n/tot;
    const col = frac>=1 ? "#27ae60" : frac>0 ? "#e67e22" : "#c0392b";
    return `<span style="color:${col};font-weight:${frac<1?"600":"normal"}">${n}</span><span style="color:#aaa">/${tot}</span>`;
  }
  function pfamLinks(d){
    return (d.pfam||[]).map(p=>
      `<a href="https://www.ebi.ac.uk/interpro/search/text/${encodeURIComponent(p)}/"
          target="_blank" rel="noopener"
          style="color:#2980b9;text-decoration:none;margin-right:6px;white-space:nowrap"
          title="Search InterPro for ${p}">${p}</a>`
    ).join("") || '<span style="color:#bbb">—</span>';
  }

  // ── group by category (cls) ──────────────────────────────────────────────
  const catMap = {};
  rows.forEach(d=>{
    const cat = d.cls || "(other)";
    (catMap[cat]||(catMap[cat]=[])).push(d);
  });
  const catKeys = Object.keys(catMap).sort();

  // auto-expand all cats when searching
  const autoExpand = q.length > 0;

  let html = "";
  catKeys.forEach(cat=>{
    const fams = catMap[cat].slice().sort((a,b)=>{
      let va=a[_famSortCol], vb=b[_famSortCol];
      if(_famSortCol==="pfam"){ va=(va||[]).join(","); vb=(vb||[]).join(","); }
      if(typeof va==="number") return _famSortAsc ? va-vb : vb-va;
      va=String(va||"").toLowerCase(); vb=String(vb||"").toLowerCase();
      return _famSortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
    });
    const expanded = autoExpand || _famExpandedCats.has(cat);
    const arrow = expanded ? "&#9660;" : "&#9654;";

    // aggregate stats
    const nHgs     = fams.reduce((s,d)=>s+(d.n_hgs||0), 0);
    const nTrees   = fams.reduce((s,d)=>s+(d.n_trees||0), 0);
    const nGenerax = fams.reduce((s,d)=>s+(d.n_generax||0), 0);
    const nGenes   = fams.reduce((s,d)=>s+(d.total||0), 0);
    // union of species across all families in category
    const spUnion = new Set();
    fams.forEach(d=>{ const fd=_famDataByName[d.family]; if(fd) Object.keys(fd.species_counts).forEach(s=>spUnion.add(s)); });
    const nSp = spUnion.size;
    // representative cls badge (all fams in cat share same cls)
    const badge = _clsBadge(fams[0].cls);

    const catKey = cat.replace(/'/g,"\\'");
    html += `<tr class="fam-cat-row" style="background:#f0f4f8;border-bottom:2px solid #d4dde8;cursor:pointer"
        onclick="famToggleCat('${catKey}')">
      <td style="padding:6px 8px;font-weight:700;white-space:nowrap;color:#1a3a5c">
        <span style="display:inline-block;width:14px;text-align:center;font-size:10px;color:#555">${arrow}</span>
        ${cat}
      </td>
      <td style="padding:6px 8px;color:#aaa;font-size:11px">—</td>
      <td style="padding:6px 8px">${badge}</td>
      <td style="padding:6px 8px;text-align:right">${speciesFrac(nSp)}</td>
      <td style="padding:6px 8px;text-align:right;color:#555;font-weight:600">${fams.length}</td>
      <td style="padding:6px 8px;text-align:right;color:#555;font-weight:600">${nHgs}</td>
      <td style="padding:6px 8px;text-align:right">${treeFrac(nTrees, nHgs)}</td>
      <td class="fam-td-generax" style="padding:6px 8px;text-align:right">${treeFrac(nGenerax, nHgs)}</td>
      <td style="padding:6px 8px;text-align:right;color:#555">${nGenes}</td>
    </tr>`;

    if(expanded){
      fams.forEach((d,i)=>{
        const bg = i%2===0?"#fff":"#fafcff";
        html += `<tr class="fam-fam-row" style="background:${bg};border-bottom:1px solid #eee">
          <td style="padding:4px 8px 4px 28px;font-weight:600;cursor:pointer;color:#2980b9;white-space:nowrap;text-decoration:underline"
              onclick="famGoToFamily('${d.family.replace(/'/g,"\\'")}')"
              title="Open HG-level counts in Counts tab">${d.family}</td>
          <td style="padding:4px 8px;font-size:11px">${pfamLinks(d)}</td>
          <td style="padding:4px 8px"><span style="color:#888;font-size:11px">${d.category||""}</span></td>
          <td style="padding:4px 8px;text-align:right">${speciesFrac(d.n_species||0)}</td>
          <td style="padding:4px 8px;text-align:right;color:#aaa;font-size:11px">—</td>
          <td style="padding:4px 8px;text-align:right;color:#555">${d.n_hgs||0}</td>
          <td style="padding:4px 8px;text-align:right">${treeFrac(d.n_trees||0, d.n_hgs||0)}</td>
          <td class="fam-td-generax" style="padding:4px 8px;text-align:right">${treeFrac(d.n_generax||0, d.n_hgs||0)}</td>
          <td style="padding:4px 8px;text-align:right;color:#555">${d.total||0}</td>
        </tr>`;
      });
    }
  });
  tbody.innerHTML = html;

  // wire sort headers (only once)
  if(!_famRendered){
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
      const proteinPopup = _proteinPopupHtml(gid, {family: currentIndex?.family || ''});
      if (proteinPopup) html += '<div style="margin-top:8px">'+proteinPopup+'</div>';
      else {
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
    cladoCollapsed.delete(k);
    spCollapsed.delete(k);
    fn(); drawSpeciesTree();
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

// ── Species-tree node action popup (collapse / flip) ────────────────────────
(function(){
  const pop=document.getElementById("sptree-node-popup");
  let _n=null, _nodeKey=null, _spId=null;
  function hide(){ pop.style.display="none"; _n=null; _nodeKey=null; _spId=null; }
  document.getElementById("snp-collapse").addEventListener("click",()=>{
    if(_nodeKey===null||!_n) return;
    const nodeKey=_nodeKey;
    const nodeName=_n.name||"";
    hide();
    spCollapsed.add(nodeKey);
    cladoCollapsed.add(nodeKey);
    if(nodeName&&!cladoNames.has(nodeKey)) cladoNames.set(nodeKey, nodeName);
    drawSpeciesTree(); drawCladogram();
    if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
  });
  document.getElementById("snp-flip").addEventListener("click",()=>{
    if(_spId===null) return;
    const spId=_spId;
    hide();
    const orig=spNodeById.get(spId);
    if(orig&&orig.children) orig.children.reverse();
    recomputeSpeciesOrder();
    drawSpeciesTree(); drawCladogram();
    if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
  });
  document.addEventListener("click",(ev)=>{ if(!pop.contains(ev.target)) hide(); });
  window.showSpNodeActionPopup=function(ev, n){
    ev.stopPropagation(); hideTip();
    _n=n; _nodeKey=spTreeNodeKey(n); _spId=n._spId;
    document.getElementById("snp-title").textContent=n.name||(_nodeKey);
    pop.style.display="flex";
    const x=Math.min(ev.clientX+4,window.innerWidth-pop.offsetWidth-8);
    const y=Math.min(ev.clientY+4,window.innerHeight-pop.offsetHeight-8);
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

function recomputeSpeciesOrder(){
  function spLeaves(n){ return n.children?n.children.flatMap(spLeaves):[n.name]; }
  const treeOrder=spLeaves(SP_TREE_DATA);
  speciesOrder=treeOrder.filter(s=>dataSpecies.has(s));
  for(const s of dataSpecies){ if(!speciesOrder.includes(s)) speciesOrder.push(s); }
}

function _avgSpeciesColor(speciesList){
  const rgbs=(speciesList||[]).map(sp=>_cssToRgb(spColor(sp))).filter(Boolean);
  if(!rgbs.length) return colTriFill;
  return _rgbToHex(d3.mean(rgbs,c=>c[0]), d3.mean(rgbs,c=>c[1]), d3.mean(rgbs,c=>c[2]));
}

function _collapsedCladeFill(node, leafGetter){
  const key=spTreeNodeKey(node);
  const custom=cladoColors.get(key);
  if(custom) return custom;
  const leaves=(typeof leafGetter==="function" ? leafGetter(node) : []).filter(Boolean);
  return _avgSpeciesColor(leaves);
}

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

const cladoCollapsed  = new Set();   // stable species-tree node keys collapsed in cladograms
const cladoNames      = new Map();   // key → custom display name (user-editable)
const cladoColors     = new Map();   // key → custom fill color override
// map: clade label → sorted array of species (updated by drawCladogram, read by drawHeatmap)
const hmCollapsedBands = [];  // [{species:[...], label}] — reset on each drawCladogram call

(function seedSpeciesGroupCollapses(){
  if(!SP_TREE_DATA||!Object.keys(SPECIES_GROUPS).length) return;
  const treeLeaves = new Set(spNodeLeaves(SP_TREE_DATA));
  const groups = {};
  Object.entries(SPECIES_GROUPS).forEach(([sp, grp])=>{
    if(!grp||!treeLeaves.has(sp)) return;
    (groups[grp]||(groups[grp]=[])).push(sp);
  });
  Object.entries(groups).forEach(([grp, sps])=>{
    const uniq = [...new Set(sps)];
    if(uniq.length < 2) return;
    const node = findSpMrcaForSpecies(SP_TREE_DATA, new Set(uniq));
    if(!node||node===SP_TREE_DATA||!node.children||node.children.length<2) return;
    const key = spTreeNodeKey(node);
    spCollapsed.add(key);
    cladoCollapsed.add(key);
    if(!cladoNames.has(key)) cladoNames.set(key, grp);
  });
})();

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
  // y-range of all leaves in a subtree
  function yRange(n){ const ys=flat(n).filter(x=>!x.children).map(x=>leafY[x.name]||0); return [d3.min(ys),d3.max(ys)]; }

  assignY(tree); assignX(tree);
  const maxD=d3.max(flat(tree),d=>d._d)||1;
  flat(tree).forEach(n=>{ if(!n.children) n._d=maxD; }); // align all tips
  const sx=d=>10+(d/maxD)*150;

  // ── branches (skip into collapsed subtrees) ──────────────────────────────
  function drawB(n){
    if(!n.children) return;
    if(cladoCollapsed.has(spTreeNodeKey(n))) return;
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
    const key=spTreeNodeKey(n);
    if(cladoCollapsed.has(key)){
      const [y0,y1]=yRange(n);
      const displayName=cladoNames.get(key)||(n.name||"clade");
      // collect leaf species for heatmap shading
      function getLeafSp(nd){ return nd.children?nd.children.flatMap(getLeafSp):[nd.name]; }
      hmCollapsedBands.push({species:getLeafSp(n),label:displayName});
      const triColor=_collapsedCladeFill(n, getLeafSp);
      const redraw=()=>{ drawCladogram(); drawSpeciesTree(); };
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
    const key=spTreeNodeKey(n);
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
          cladoCollapsed.add(key); spCollapsed.add(key); drawCladogram(); drawSpeciesTree();
          if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
        }
      });
    n.children.forEach(drawCladoNodes);
  }
  drawCladoNodes(tree);

  // ── species logos (optional) ─────────────────────────────────────────────
  if(hmShowSpLogos){
    const imgSz=12, imgX=sx(maxD)+4;
    speciesOrder.forEach(sp=>{
      const src=SPECIES_IMAGES[sp];
      if(!src) return;
      const y=leafY[sp];
      if(y===undefined) return;
      svg.append("image")
        .attr("x",imgX).attr("y",y-imgSz/2)
        .attr("width",imgSz).attr("height",imgSz)
        .attr("href",src)
        .attr("preserveAspectRatio","xMidYMid meet");
    });
  }
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

(function(){
  const pop=document.getElementById("hm-open-popup");
  const title=document.getElementById("hm-open-popup-title");
  const treeBtn=document.getElementById("hm-open-tree");
  const alignBtnEl=document.getElementById("hm-open-align");
  let _onTree=null, _onAlign=null;
  function hide(){ pop.style.display="none"; _onTree=null; _onAlign=null; }
  treeBtn.addEventListener("click",()=>{ const fn=_onTree; hide(); if(fn) fn(); });
  alignBtnEl.addEventListener("click",()=>{ const fn=_onAlign; hide(); if(fn) fn(); });
  document.addEventListener("click",(ev)=>{ if(pop.style.display==="block"&&!pop.contains(ev.target)) hide(); });
  window.showHmOpenPopup=function(ev,label,onTree,onAlign,treeLabel,alignLabel){
    ev.stopPropagation();
    hideTip();
    _onTree=onTree; _onAlign=onAlign;
    title.textContent=label||"Open in";
    treeBtn.textContent = treeLabel || 'Show Tree';
    alignBtnEl.textContent = alignLabel || 'Show Alignment';
    treeBtn.style.display = typeof onTree === "function" ? "" : "none";
    if(alignBtnEl) alignBtnEl.style.display=(HAVE_ALIGNMENTS && typeof onAlign==="function") ? "" : "none";
    pop.style.display="block";
    const x=Math.min(ev.clientX+8,window.innerWidth-pop.offsetWidth-8);
    const y=Math.min(ev.clientY+8,window.innerHeight-pop.offsetHeight-8);
    pop.style.left=Math.max(4,x)+"px";
    pop.style.top=Math.max(4,y)+"px";
  };
})();

function hmOpenTreeForHG(tRec, sp){
  if(!tRec) return;
  hmFocusGids=null;
  switchTab("trees");
  selectTree(tRec);
  renderSidebar("");
  if(sp){
    hlQueries=[sp];
    rebuildHlSet();
    setTimeout(focusHighlighted,60);
  }
}

function hmOpenTreeForOG(treeRec, ogId, sp){
  if(!treeRec) return;
  const det=loadDetail(treeRec.id);
  const ogGids=(det&&det.ogs&&det.ogs[ogId])||[];
  const targetGids=sp?ogGids.filter(g=>getSpeciesPfx(g)===sp):ogGids;
  switchTab("trees");
  selectTree(treeRec);
  renderSidebar("");
  if(targetGids.length){
    hmFocusGids=new Set(targetGids);
    hideNonHl=true;
    const chk=document.getElementById("chk-hide-nonhl"); if(chk) chk.checked=true;
  }
  setTimeout(focusHighlighted,80);
}

function hmOpenAlignmentForHG(alnId){
  if(!HAVE_ALIGNMENTS) return;
  if(!alnId) return;
  switchTab("align");
  selectAlignment(alnId);
}

function hmOpenAlignmentForOG(treeRec, ogId){
  if(!HAVE_ALIGNMENTS) return;
  if(!treeRec) return;
  switchTab("align");
  selectAlignment(treeRec.id);
  const allOGs=[...new Set(Object.values(alnOgMap))];
  alnHiddenOGs.clear();
  alnCollapsedOGs.clear();
  allOGs.forEach(og=>{ if(og!==ogId) alnCollapsedOGs.add(og); });
  renderAlignment();
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

  // groupColors is recomputed whenever species colors change

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
  function avgTipColor(leaves){
    const rgbs=leaves.map(l=>_cssToRgb(tipColor(l))).filter(Boolean);
    if(!rgbs.length) return colTriFill;
    return _rgbToHex(d3.mean(rgbs,c=>c[0]), d3.mean(rgbs,c=>c[1]), d3.mean(rgbs,c=>c[2]));
  }

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
    if(spCollapsed.has(spTreeNodeKey(n))){
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
      const x0=sx(n._d), xTip=sx(maxD), x1=xTip-8;
      const y0=n._colY0, y1=n._colY0+n._colH, yM=n._y;
      const leaves=getLeaves(n);
      // horizontal branch from parent edge to apex already drawn; draw triangle
      const cladoKey = spTreeNodeKey(n);
      const displayName = cladoNames.get(cladoKey)||(n.name||"clade");
      const triFill = _collapsedCladeFill(n, node => getLeaves(node).map(l => l.name));
      const expandNode = ()=>{ spCollapsed.delete(cladoKey); cladoCollapsed.delete(cladoKey); drawSpeciesTree(); drawCladogram(); if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap(); };
      svg.append("polygon")
        .attr("points",`${x0},${yM} ${x1},${y0} ${x1},${y1}`)
        .attr("fill",triFill).attr("stroke","#222").attr("stroke-width",1)
        .style("cursor","pointer")
        .on("click",expandNode)
        .on("mouseover",ev=>showTip(ev,'<div class="tt-name">'+displayName+'</div><div style="font-size:9px;color:#aaa">click to expand</div>'))
        .on("mousemove",moveTip).on("mouseout",hideTip);
      // small circle at apex for easy uncollapse
      svg.append("circle")
        .attr("cx",x0).attr("cy",yM).attr("r",5)
        .attr("fill","#f5f5f5").attr("stroke","#555").attr("stroke-width",1.2)
        .style("cursor","pointer")
        .on("click",expandNode)
        .on("mouseover",ev=>showTip(ev,'<div class="tt-name">'+displayName+'</div><div style="font-size:9px;color:#aaa">click to expand</div>'))
        .on("mousemove",moveTip).on("mouseout",hideTip);
      // node name badge at apex
      if(displayName){
        svg.append("text").attr("x",x0+4).attr("y",yM-3)
          .attr("font-size",8).attr("fill","#5a7090").attr("font-style","italic")
          .text(displayName+" ["+leaves.length+"]");
      }
      // species names at triangle base — same layout as drawLeaves
      const nameH=n._colH/leaves.length;
      const imgW=24, imgH=18, imgGap=4;
      leaves.forEach((l,i)=>{
        const ly=n._colY0+i*nameH+nameH/2;
        const tc=tipColor(l);
        const imgSrc=SPECIES_IMAGES[l.name]||"";
        const showImg=imgSrc&&nameH>=14;
        const textX=xTip+10+(showImg?imgW+imgGap:0);
        svg.append("circle").attr("cx",xTip).attr("cy",ly).attr("r",5)
          .attr("fill",tc).attr("stroke","#fff").attr("stroke-width",0.8);
        if(showImg){
          svg.append("image")
            .attr("x",xTip+10).attr("y",ly-imgH/2)
            .attr("width",imgW).attr("height",imgH)
            .attr("href",imgSrc)
            .attr("preserveAspectRatio","xMidYMid meet");
        }
        svg.append("text").attr("x",textX).attr("y",ly).attr("dy","0.35em")
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
    // root node: draw a flip-only dot (no collapse)
    if(n._id===tree._id){
      svg.append("circle").attr("cx",sx(n._d)).attr("cy",n._y).attr("r",5)
        .attr("fill","#fff").attr("stroke","#ccc").attr("stroke-width",1.2)
        .style("cursor","pointer")
        .on("mouseover",ev=>showTip(ev,(n.name||"root")+'<div style="font-size:9px;color:#aaa;margin-top:3px">click to flip</div>'))
        .on("mousemove",moveTip).on("mouseout",hideTip)
        .on("click",(ev)=>{
          if(ev.shiftKey) return;
          ev.stopPropagation(); hideTip();
          const orig=spNodeById.get(n._spId);
          if(orig&&orig.children) orig.children.reverse();
          recomputeSpeciesOrder();
          drawSpeciesTree(); drawCladogram();
          if(document.getElementById("pane-heatmap").classList.contains("active")) drawHeatmap();
        });
      n.children.forEach(drawInternals);
      return;
    }
    if(n._id!==tree._id){
      const nodeLabel=n.name||"(unnamed)";
      // check if this clade is already a split group
      const cleavesSet=new Set(spNodeLeaves(n));
      const splitIdx=hmSplitSets.findIndex(s=>s.size===cleavesSet.size&&[...cleavesSet].every(v=>s.has(v)));
      const isSplit=splitIdx>=0;
      const nodeCol=isSplit?hmSplitLineColors[splitIdx%hmSplitLineColors.length]:"#999";
      const tip=nodeLabel+'<div style="font-size:9px;color:#aaa;margin-top:3px">click for options &nbsp;·&nbsp; shift+click to split heatmap</div>';
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
            showSpNodeActionPopup(ev, n);
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
            refreshSpeciesColorViews();
          });
        });
      const imgSrc=SPECIES_IMAGES[n.name]||"";
      const imgW=24, imgH=18, imgGap=4;
      const textX=sx(n._d)+10+(imgSrc?imgW+imgGap:0);
      if(imgSrc){
        svg.append("image")
          .attr("x",sx(n._d)+10).attr("y",n._y-imgH/2)
          .attr("width",imgW).attr("height",imgH)
          .attr("href",imgSrc)
          .attr("preserveAspectRatio","xMidYMid meet")
          .style("cursor","pointer")
          .on("mouseover",(ev)=>showTip(ev,makeTip())).on("mousemove",moveTip).on("mouseout",hideTip)
          .on("click",(ev)=>{ hideTip(); showSpAnnotPopup(ev, n.name); });
      }
      const grp=SPECIES_GROUPS[n.name]||"";
      const grpCol=grp?groupColors[grp]:"";
      function makeTip(){
        const m=spMeta[n.name]||{genes:0,families:0,hgs:0};
        const fullName=SPECIES_INFO[n.name]||"";
        return (imgSrc?'<img src="'+imgSrc+'" style="max-width:220px;max-height:150px;display:block;border-radius:3px;margin-bottom:6px">':"")
          +'<div class="tt-name" style="color:'+tc+'">'+n.name+'</div>'
          +(fullName?'<div style="font-style:italic;color:#555;font-size:11px;margin-bottom:4px">'+fullName+'</div>':"")
          +(grp?'<div style="font-size:10px;margin-bottom:4px"><span style="display:inline-block;width:9px;height:9px;border-radius:2px;background:'+grpCol+';margin-right:4px;vertical-align:middle"></span>'+grp+'</div>':"")
          +'<div class="tt-row"><span>Annotated genes</span><strong>'+m.genes+'</strong></div>'
          +'<div class="tt-row"><span>Families</span><strong>'+m.families+'</strong></div>'
          +'<div class="tt-row"><span>Homology groups</span><strong>'+m.hgs+'</strong></div>'
          +'<div style="font-size:9px;color:#aaa;margin-top:3px">click to export annotations</div>';
      }
      svg.append("text").attr("x",textX).attr("y",n._y).attr("dy","0.35em")
        .attr("font-size",12).attr("fill",inData?"#222":"#bbb").attr("font-family","monospace").text(n.name)
        .style("cursor","pointer")
        .on("mouseover",(ev)=>showTip(ev,makeTip())).on("mousemove",moveTip).on("mouseout",hideTip)
        .on("click",(ev)=>{ hideTip(); showSpAnnotPopup(ev, n.name); });
      if(grp && grpCol){
        const pillH=13, pillRx=3, badgeX=W-8;
        const estW=grp.length*6.5+10;
        svg.append("rect")
          .attr("x",badgeX-estW).attr("y",n._y-pillH/2)
          .attr("width",estW).attr("height",pillH).attr("rx",pillRx)
          .attr("fill",grpCol).attr("opacity",0.85);
        svg.append("text")
          .attr("x",badgeX-5).attr("y",n._y).attr("dy","0.35em")
          .attr("text-anchor","end").attr("font-size",9).attr("fill","#fff")
          .attr("font-weight","600").text(grp)
          .on("mouseover",(ev)=>showTip(ev,makeTip())).on("mousemove",moveTip).on("mouseout",hideTip);
      }
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

function hmIsNoAnnotationRec(d){
  return !!(d&&typeof d.id==="string"&&d.id.startsWith("(no annotation)"));
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
      const tRec=TREE_INDEX.find(t=>t.id===d.id)||TREE_INDEX.find(t=>t.family===d.family&&t.hg===d.hg);
      if(!tRec) return;
      showHmOpenPopup(
        ev,
        (d.hg||d.id)+" · "+sp,
        ()=>hmOpenTreeForHG(tRec, d.species_counts[sp]?sp:null),
        ()=>hmOpenAlignmentForHG(tRec.id)
      );
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
    clickHandler=(ev,_d,sp)=>{
      if(!treeRec) return;
      showHmOpenPopup(
        ev,
        _d.id+(sp?" · "+sp:""),
        ()=>hmOpenTreeForOG(treeRec, _d.id, sp),
        ()=>hmOpenAlignmentForOG(treeRec, _d.id)
      );
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
    clickHandler=(ev,d,sp)=>{
      const tRec=hmOGIndex[d.id]?TREE_INDEX.find(t=>t.id===hmOGIndex[d.id].hgId):null;
      if(!tRec) return;
      showHmOpenPopup(
        ev,
        d.id+(sp?" · "+sp:""),
        ()=>hmOpenTreeForOG(tRec, d.id, sp),
        ()=>hmOpenAlignmentForOG(tRec, d.id)
      );
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

  // ── column grouping by HG / family / class ─────────────────────────────
  let _colGrpKey = null; // fn: d → group key string
  if(hmGroupByHG){
    if(hmViewMode==="custom"||hmViewMode==="og"){ buildOGIndex(); _colGrpKey=d=>hmOGIndex[d.id]?.hgId||"?"; }
    else if(hmViewMode==="hg")     _colGrpKey=d=>d.family||"?";
    else if(hmViewMode==="family") _colGrpKey=d=>d.class||"?";
  }
  let hmGrpBoundaries=new Set(); // col indices where a new group starts (excl. 0)
  let hmGrpMeta=[];              // [{key, startIdx, count}]
  if(_colGrpKey && data.length){
    // collect ordered unique group keys
    const seenGrps=[], grpSeen=new Set();
    data.forEach(d=>{ const k=_colGrpKey(d); if(!grpSeen.has(k)){grpSeen.add(k);seenGrps.push(k);} });
    let orderedGrps=seenGrps;
    if(hmGroupOrderOverride&&hmGroupOrderOverride.length){
      const inOv=new Set(hmGroupOrderOverride);
      orderedGrps=[...hmGroupOrderOverride.filter(k=>grpSeen.has(k)),...seenGrps.filter(k=>!inOv.has(k))];
    }
    const grpRank=new Map(orderedGrps.map((k,i)=>[k,i]));
    // stable sort by group rank, preserving intra-group column order
    const origPos=new Map(data.map((d,i)=>[d.id,i]));
    data=[...data].sort((a,b)=>{
      const ga=grpRank.get(_colGrpKey(a))??999, gb=grpRank.get(_colGrpKey(b))??999;
      if(ga!==gb) return ga-gb;
      return (origPos.get(a.id)??0)-(origPos.get(b.id)??0);
    });
    // compute boundaries and group metadata
    let prevKey=_colGrpKey(data[0]), grpStart=0;
    hmGrpMeta.push({key:prevKey,startIdx:0,count:0});
    for(let i=1;i<data.length;i++){
      const k=_colGrpKey(data[i]);
      if(k!==prevKey){
        hmGrpBoundaries.add(i);
        hmGrpMeta[hmGrpMeta.length-1].count=i-grpStart;
        grpStart=i; prevKey=k;
        hmGrpMeta.push({key:k,startIdx:i,count:0});
      }
    }
    hmGrpMeta[hmGrpMeta.length-1].count=data.length-grpStart;
  }
  const HM_GROUP_H = hmGrpMeta.length>1 ? 20 : 0;

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

  // ── FLIPPED MODE: species as columns, HG/OGs as rows ─────────────────
  if(hmFlipped){
    const tp=document.getElementById("tree-panel"); if(tp) tp.style.display="none";
    const cW=18, cH=14;
    const colLabelMaxCharsF=Math.round(hmColFontSize*3.5);
    const colLabelTruncF=d=>{ const s=colLabel(d)||""; return s.length>colLabelMaxCharsF?s.slice(0,colLabelMaxCharsF-1)+"\u2026":s; };
    const maxHGLen=data.reduce((m,d)=>Math.max(m,(colLabelTruncF(d)||"").length),0);
    const ROW_LABEL_W_F=Math.max(110,Math.min(240,maxHGLen*hmColFontSize*0.65+16))+(hmGrpMeta.length>1?HM_GROUP_H:0);
    const SP_TREE_H=64;
    const maxSpLen=hmOrder.reduce((m,s)=>Math.max(m,s.length),0);
    const SP_LABEL_H=Math.ceil(maxSpLen*hmColFontSize*0.62)+8;
    const hmTM_F=SP_TREE_H+SP_LABEL_H+HM_GROUP_H;
    hmTopActual=hmTM_F;
    const nRows_F=data.length, nCols_F=hmOrder.length;
    const CELLS_BOTTOM_F=nRows_F*cH+hmTM_F;
    const BAR_LEFT=nCols_F*cW+ROW_LABEL_W_F+8;
    const BAR_W_MAX=50;
    const svgW_F=BAR_LEFT+BAR_W_MAX+32;
    const svgH_F=CELLS_BOTTOM_F+40;
    const panel=document.getElementById("heatmap-panel");
    const svg=d3.select(panel).html("").append("svg").attr("width",svgW_F).attr("height",svgH_F);
    if(!data.length){ svg.append("text").attr("x",40).attr("y",hmTM_F+30).attr("fill","#999").text("No data."); return; }

    // ── z-score per row ──
    const zMat_F=data.map(rec=>{ const v=hmOrder.map(s=>rec.species_counts[s]||0); const m=d3.mean(v),sd=d3.deviation(v)||1; return v.map(x=>(x-m)/sd); });
    const zAll_F=zMat_F.flat();
    const zMax_F=Math.max(Math.abs(d3.min(zAll_F)||0),Math.abs(d3.max(zAll_F)||0),0.01);
    const zColor_F=d3.scaleDiverging().domain([zMax_F,0,-zMax_F]).interpolator(d3.interpolateRdBu);
    const maxAbs_F=d3.max(data.flatMap(rec=>hmOrder.map(s=>rec.species_counts[s]||0)))||1;
    const absColor_F=d3.scaleSequential().domain([0,maxAbs_F]).interpolator(d3.interpolateBlues);
    const color_F=hmColorMode==="absolute"?absColor_F:zColor_F;

    // ── horizontal species tree ──
    (function(){
      function clF(n){ return JSON.parse(JSON.stringify(n)); }
      function prF(n){ if(!n.children) return hmOrder.includes(n.name)?n:null; const k=n.children.map(prF).filter(Boolean); if(!k.length) return null; if(k.length===1) return k[0]; n.children=k; return n; }
      const htree=prF(clF(SP_TREE_DATA)); if(!htree) return;
      const leafX={}; hmOrder.forEach((sp,ci)=>{ leafX[sp]=ROW_LABEL_W_F+ci*cW+cW/2; });
      function aX(n){ if(!n.children){n._lx=leafX[n.name]||0;return n._lx;} const xs=n.children.map(aX); n._lx=d3.mean(xs);return n._lx; }
      function aD(n,d=0){ n._d=d; if(n.children) n.children.forEach(c=>aD(c,d+1)); }
      function flatF(n){ return [n].concat(n.children?n.children.flatMap(flatF):[]); }
      function nkF(n){ return spTreeNodeKey(n); }
      aX(htree); aD(htree);
      const maxDF=d3.max(flatF(htree),n=>n._d)||1;
      flatF(htree).forEach(n=>{ if(!n.children) n._d=maxDF; });
      const syF=d=>4+(d/maxDF)*(SP_TREE_H-10);
      // branches
      function drawHB(n){ if(!n.children||cladoCollapsed.has(nkF(n))) return; const xs=n.children.map(c=>c._lx); svg.append("line").attr("x1",d3.min(xs)).attr("x2",d3.max(xs)).attr("y1",syF(n._d)).attr("y2",syF(n._d)).attr("stroke","#bbb").attr("stroke-width",1); n.children.forEach(c=>{ svg.append("line").attr("x1",c._lx).attr("x2",c._lx).attr("y1",syF(n._d)).attr("y2",syF(c._d)).attr("stroke","#bbb").attr("stroke-width",1); drawHB(c); }); }
      drawHB(htree);
      // interactive nodes
      function drawHN(n){ if(!n.children) return; const key=nkF(n); if(cladoCollapsed.has(key)) return; const clvs=new Set(flatF(n).filter(x=>!x.children).map(x=>x.name)); const si=hmSplitSets.findIndex(s=>s.size===clvs.size&&[...clvs].every(v=>s.has(v))); const isSpl=si>=0; const lbl=n.name||(clvs.size+" spp."); svg.append("circle").attr("cx",n._lx).attr("cy",syF(n._d)).attr("r",4).attr("fill",isSpl?hmSplitColors[si%hmSplitColors.length]:"#fff").attr("stroke",isSpl?hmSplitLineColors[si%hmSplitLineColors.length]:"#aaa").attr("stroke-width",isSpl?2:1).style("cursor","pointer").on("mouseover",ev=>showTip(ev,lbl+'<div style="font-size:9px;color:#aaa">click=collapse &nbsp;·&nbsp; shift+click=split</div>')).on("mousemove",moveTip).on("mouseout",hideTip).on("click",(ev)=>{ hideTip(); if(ev.shiftKey){ ev.stopPropagation(); if(isSpl){hmSplitSets.splice(si,1);hmSplitLabels.splice(si,1);}else{hmSplitSets.push(clvs);hmSplitLabels.push(lbl);} updateHmSplitBar(); drawHeatmap(); } else { cladoCollapsed.add(key); spCollapsed.add(key); drawHeatmap(); } }); n.children.forEach(drawHN); }
      drawHN(htree);
      // logos
      if(hmShowSpLogos){ hmOrder.forEach(sp=>{ const src=SPECIES_IMAGES[sp]; if(!src) return; svg.append("image").attr("x",leafX[sp]-6).attr("y",SP_TREE_H-13).attr("width",12).attr("height",12).attr("href",src).attr("preserveAspectRatio","xMidYMid meet"); }); }
      // species column labels (below tree, rotated)
      hmOrder.forEach((sp,ci)=>{ svg.append("text").attr("transform",`translate(${ROW_LABEL_W_F+ci*cW+cW/2},${SP_TREE_H+SP_LABEL_H-4}) rotate(-90)`).attr("text-anchor","start").attr("font-size",hmColFontSize).attr("fill","#333").text(sp); });
    })();

    // ── column split bands (species split groups become vertical bands) ──
    if(nSplitGroups>0){
      hmOrder.forEach((sp,ci)=>{ const gi=spSplitGroup.get(sp); if(gi===undefined) return; svg.append("rect").attr("x",ci*cW+ROW_LABEL_W_F).attr("y",hmTM_F).attr("width",cW).attr("height",nRows_F*cH).attr("fill",bandFillColors[gi%bandFillColors.length]).attr("pointer-events","none"); });
      for(let ci=1;ci<nCols_F;ci++){ const pg=spSplitGroup.get(hmOrder[ci-1]),cg=spSplitGroup.get(hmOrder[ci]); if(pg!==cg) svg.append("line").attr("x1",ci*cW+ROW_LABEL_W_F).attr("x2",ci*cW+ROW_LABEL_W_F).attr("y1",hmTM_F).attr("y2",CELLS_BOTTOM_F).attr("stroke","#999").attr("stroke-width",1.5).attr("stroke-dasharray","4,2"); }
    }

    // ── cells (transposed: rows=HGs, cols=species) ──
    data.forEach((rec,ri)=>{
      hmOrder.forEach((sp,ci)=>{
        const count=rec.species_counts[sp]||0, z=zMat_F[ri][ci];
        const cellX=ci*cW+ROW_LABEL_W_F, cellY=ri*cH+hmTM_F;
        const fill=hmColorMode==="absolute"?absColor_F(count):zColor_F(z);
        svg.append("rect").attr("class","hm-cell").attr("x",cellX).attr("y",cellY).attr("width",cW-2).attr("height",cH-2).attr("fill",count===0?"#fff":fill).style("cursor","pointer")
          .on("mouseover",ev=>showTip(ev,'<b>'+sp+'</b><br>'+rec.id+'<br>count: <b>'+count+'</b>'+(hmColorMode==="zscore"?'<br>z: <b>'+z.toFixed(2)+'</b>':"")))
          .on("mousemove",moveTip).on("mouseout",hideTip).on("click",ev=>clickHandler(ev,rec,sp));
        if(count>0){ const fc=d3.color(fill); const lum=fc?(0.2126*fc.r+0.7152*fc.g+0.0722*fc.b)/255:1; svg.append("text").attr("x",cellX+(cW-2)/2).attr("y",cellY+cH/2).attr("dy","0.35em").attr("text-anchor","middle").attr("font-size",6.5).attr("fill",lum>0.45?"#222":"#eee").attr("pointer-events","none").text(count>999?"\u226b":count); }
      });
    });

    // ── group HG horizontal separators + band labels ──
    if(hmGrpMeta.length>1){
      const GBW=HM_GROUP_H-2;
      let _gdk=null,_gdt=null,_gdm=false,_gdy=0;
      hmGrpMeta.forEach(g=>{
        const y0=g.startIdx*cH+hmTM_F, h=g.count*cH, midY=y0+h/2;
        const bx=ROW_LABEL_W_F-GBW-2;
        svg.append("rect").attr("x",bx).attr("y",y0+1).attr("width",GBW).attr("height",h-2).attr("rx",3).attr("fill","#e8eef6").attr("stroke","#b0bcd0").attr("stroke-width",1).style("cursor","grab")
          .on("mouseover",ev=>showTip(ev,"<b>"+g.key+"</b><br><span style='color:#aaa;font-size:10px'>"+g.count+" row"+(g.count!==1?"s":"")+" — drag to reorder</span>")).on("mousemove",moveTip).on("mouseout",hideTip)
          .call(d3.drag()
            .on("start",(ev)=>{ _gdk=g.key; _gdm=false; _gdy=ev.y; })
            .on("drag",(ev)=>{ if(Math.abs(ev.y-_gdy)>6) _gdm=true; if(!_gdm) return; let cy=hmTM_F,tgt=hmGrpMeta.length; for(let i=0;i<hmGrpMeta.length;i++){if(ev.y<cy+hmGrpMeta[i].count*cH/2){tgt=i;break;} cy+=hmGrpMeta[i].count*cH;} _gdt=tgt; })
            .on("end",()=>{ hideTip(); if(_gdm&&_gdk!==null&&_gdt!==null){ const ord=hmGrpMeta.map(x=>x.key); const si=ord.indexOf(_gdk); if(si>=0){ ord.splice(si,1); ord.splice(Math.max(0,Math.min(ord.length,_gdt>si?_gdt-1:_gdt)),0,_gdk); hmGroupOrderOverride=[...ord]; drawHeatmap(); } } _gdk=null;_gdt=null;_gdm=false; }));
        // vertical label text
        const maxC=Math.max(2,Math.floor(h/7));
        const lbl=g.key.length>maxC?g.key.slice(0,maxC-1)+"\u2026":g.key;
        svg.append("text").attr("x",bx+GBW/2).attr("y",midY).attr("text-anchor","middle").attr("dominant-baseline","middle").attr("transform",`rotate(-90,${bx+GBW/2},${midY})`).attr("font-size",8).attr("font-weight","600").attr("fill","#2c4a7c").attr("pointer-events","none").text(lbl);
      });
      hmGrpBoundaries.forEach(ri=>{ const y=ri*cH+hmTM_F; svg.append("line").attr("x1",ROW_LABEL_W_F-GBW-2).attr("x2",nCols_F*cW+ROW_LABEL_W_F).attr("y1",y).attr("y2",y).attr("stroke","#1a1a1a").attr("stroke-width",1.5).attr("pointer-events","none"); });
    }

    // ── group header band above cells (species grouping from hmSplitSets) ──
    if(HM_GROUP_H>0 && hmGrpMeta.length<=1){
      // no group HG active but HM_GROUP_H is 0 anyway
    }

    // ── row labels with drag (HG/OG names on left) ──
    const _dg=svg.append("rect").attr("x",ROW_LABEL_W_F).attr("width",nCols_F*cW).attr("fill","rgba(74,144,217,0.13)").attr("stroke","#4a90d9").attr("stroke-width",1).attr("opacity",0).attr("pointer-events","none");
    const _dm=svg.append("line").attr("x1",ROW_LABEL_W_F).attr("x2",nCols_F*cW+ROW_LABEL_W_F).attr("stroke","#4a90d9").attr("stroke-width",2.5).attr("opacity",0).attr("pointer-events","none");
    let _di=null,_dt=null,_dmv=false,_dsy=0;
    svg.selectAll("text.hm-row-lbl").data(data).enter().append("text").attr("class","hm-row-lbl")
      .attr("x",ROW_LABEL_W_F-(hmGrpMeta.length>1?HM_GROUP_H:0)-4)
      .attr("y",(d,ri)=>ri*cH+hmTM_F+cH/2)
      .attr("dominant-baseline","middle").attr("text-anchor","end").attr("font-size",hmColFontSize).style("cursor","grab")
      .attr("fill",d=>hmIsNoAnnotationRec(d)?"#9aa0a6":"#333")
      .text(colLabelTruncF)
      .call(d3.drag()
        .on("start",(ev,d)=>{ _di=data.findIndex(r=>r.id===d.id); _dmv=false; _dsy=ev.y; _dg.attr("y",_di*cH+hmTM_F).attr("height",cH).attr("opacity",1); svg.style("cursor","grabbing"); })
        .on("drag",(ev)=>{ if(Math.abs(ev.y-_dsy)>5) _dmv=true; if(!_dmv) return; const t=Math.max(0,Math.min(data.length,Math.round((ev.y-hmTM_F)/cH))); _dt=t; _dm.attr("y1",t*cH+hmTM_F).attr("y2",t*cH+hmTM_F).attr("opacity",1); })
        .on("end",(ev,d)=>{ _dg.attr("opacity",0); _dm.attr("opacity",0); svg.style("cursor",""); if(!_dmv){ clickHandler(ev,d,null); } else if(_di!==null&&_dt!==null&&_dt!==_di&&_dt!==_di+1){ const nd=[...data]; const [mv]=nd.splice(_di,1); nd.splice(_dt>_di?_dt-1:_dt,0,mv); hmColOrderOverride=nd.map(r=>r.id); drawHeatmap(); } _di=null;_dt=null; }));

    // ── bar chart right (per-row totals) ──
    const rowTots=data.map(rec=>d3.sum(hmOrder.map(s=>rec.species_counts[s]||0)));
    const maxRT=d3.max(rowTots)||1;
    svg.append("line").attr("x1",BAR_LEFT).attr("x2",BAR_LEFT).attr("y1",hmTM_F).attr("y2",CELLS_BOTTOM_F).attr("stroke","#ccc").attr("stroke-width",0.8);
    svg.append("text").attr("x",BAR_LEFT+2).attr("y",hmTM_F-4).attr("font-size",8).attr("fill","#aaa").text("total \u2192");
    data.forEach((rec,ri)=>{ const t=rowTots[ri]; const bW=Math.max(t>0?1:0,(t/maxRT)*BAR_W_MAX); svg.append("rect").attr("x",BAR_LEFT).attr("y",ri*cH+hmTM_F+1).attr("width",bW).attr("height",cH-2).attr("fill","#7fb3d3").attr("opacity",0.8); if(t>0) svg.append("text").attr("x",BAR_LEFT+bW+2).attr("y",ri*cH+hmTM_F+cH/2).attr("dominant-baseline","middle").attr("font-size",6.5).attr("fill","#555").text(t>=10000?(t/1000).toFixed(0)+"k":t); });

    // ── group header row above cells (for species column groups from hmGrpMeta when grouped by species-level) ──
    if(HM_GROUP_H>0 && hmGrpMeta.length>1){
      // already drawn group bands on left above; also draw column group header if species are grouped
    }

    // ── colorbar ──
    const cbX_F=ROW_LABEL_W_F, cbY_F=svgH_F-20, cbW_F=80, cbH_F=8;
    const gid_F="hmgf"+Date.now(); const df=svg.append("defs"); const gf=df.append("linearGradient").attr("id",gid_F).attr("x1","0%").attr("x2","100%");
    for(let i=0;i<=10;i++) gf.append("stop").attr("offset",(i*10)+"%").attr("stop-color",color_F(zMax_F-i*2*zMax_F/10));
    svg.append("rect").attr("x",cbX_F).attr("y",cbY_F).attr("width",cbW_F).attr("height",cbH_F).attr("fill","url(#"+gid_F+")").attr("rx",2);
    ["high","mean","low"].forEach((l,i)=>svg.append("text").attr("x",cbX_F+i*cbW_F/2).attr("y",cbY_F+cbH_F+9).attr("font-size",8).attr("fill","#888").attr("text-anchor","middle").text(l));

    return; // skip normal rendering; drawCladogram NOT called
  }
  // ── restore tree-panel visibility (normal mode) ──
  { const tp=document.getElementById("tree-panel"); if(tp) tp.style.display=""; }

  const cW=18, cH=12;
  const maxNameLen=hmOrder.reduce((m,s)=>Math.max(m,s.length),0);
  const ROW_LABEL_W=Math.max(110,Math.min(200,maxNameLen*7+14));
  // truncation limit scales with font size: bigger font = more characters shown
  const colLabelMaxChars = Math.round(hmColFontSize * 3.5);
  const colLabelTrunc = d => { const s=colLabel(d)||""; return s.length>colLabelMaxChars ? s.slice(0,colLabelMaxChars-1)+"\u2026" : s; };
  // top margin: enough room for rotated column labels (using truncated length)
  const maxColLen=data.reduce((m,d)=>Math.max(m,(colLabelTrunc(d)||"").length),0);
  const hmColLabelH=Math.ceil(Math.sin(hmColRotation*Math.PI/180)*maxColLen*hmColFontSize*0.6)+8;
  const hmTM=Math.max(HM_TOP, hmColLabelH+14)+HM_GROUP_H;
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
    .attr("transform",(d,i)=>`translate(${i*cW+ROW_LABEL_W+cW/2},${hmTM-6-HM_GROUP_H}) rotate(-${hmColRotation})`)
    .attr("font-size",hmColFontSize).style("cursor","grab")
    .attr("text-anchor","start")
    .attr("fill",d=>hmIsNoAnnotationRec(d)?"#9aa0a6":"#333")
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

  // ── group headers + vertical separators ───────────────────────────────
  if(hmGrpMeta.length>1){
    const grpHdrY=hmTM-HM_GROUP_H;
    let _hmGrpDragKey=null,_hmGrpDragTgt=null,_hmGrpDragMoved=false,_hmGrpDragSX=0;
    hmGrpMeta.forEach(g=>{
      const x0=g.startIdx*cW+ROW_LABEL_W;
      const w=g.count*cW;
      const midX=x0+w/2;
      // drag rect
      svg.append("rect")
        .attr("x",x0+1).attr("y",grpHdrY)
        .attr("width",w-2).attr("height",HM_GROUP_H-4).attr("rx",3)
        .attr("fill","#e8eef6").attr("stroke","#b0bcd0").attr("stroke-width",1)
        .style("cursor","grab")
        .on("mouseover",ev=>showTip(ev,"<b>"+g.key+"</b><br><span style='color:#aaa;font-size:10px'>"+g.count+" column"+(g.count!==1?"s":"")+" — drag to reorder group</span>"))
        .on("mousemove",moveTip).on("mouseout",hideTip)
        .call(d3.drag()
          .on("start",(event)=>{ _hmGrpDragKey=g.key; _hmGrpDragMoved=false; _hmGrpDragSX=event.x; })
          .on("drag",(event)=>{
            if(Math.abs(event.x-_hmGrpDragSX)>6) _hmGrpDragMoved=true;
            if(!_hmGrpDragMoved) return;
            // estimate target group index from mouse x
            let cumX=ROW_LABEL_W, tgt=hmGrpMeta.length;
            for(let gi=0;gi<hmGrpMeta.length;gi++){
              if(event.x<cumX+hmGrpMeta[gi].count*cW/2){ tgt=gi; break; }
              cumX+=hmGrpMeta[gi].count*cW;
            }
            _hmGrpDragTgt=tgt;
          })
          .on("end",()=>{
            hideTip();
            if(_hmGrpDragMoved&&_hmGrpDragKey!==null&&_hmGrpDragTgt!==null){
              const ordered=hmGrpMeta.map(x=>x.key);
              const srcIdx=ordered.indexOf(_hmGrpDragKey);
              if(srcIdx>=0){
                ordered.splice(srcIdx,1);
                const ins=Math.max(0,Math.min(ordered.length,_hmGrpDragTgt>srcIdx?_hmGrpDragTgt-1:_hmGrpDragTgt));
                ordered.splice(ins,0,_hmGrpDragKey);
                hmGroupOrderOverride=[...ordered];
                drawHeatmap();
              }
            }
            _hmGrpDragKey=null; _hmGrpDragTgt=null; _hmGrpDragMoved=false;
          }));
      // group key label (top line)
      const maxChars=Math.max(2,Math.floor(w/5.5));
      svg.append("text")
        .attr("x",midX).attr("y",grpHdrY+7)
        .attr("text-anchor","middle").attr("dominant-baseline","middle")
        .attr("font-size",9).attr("font-weight","600").attr("fill","#2c4a7c")
        .attr("pointer-events","none")
        .text(g.key.length>maxChars?g.key.slice(0,maxChars-1)+"\u2026":g.key);
      // column count badge (bottom line)
      svg.append("text")
        .attr("x",midX).attr("y",grpHdrY+HM_GROUP_H-6)
        .attr("text-anchor","middle").attr("dominant-baseline","middle")
        .attr("font-size",7).attr("fill","#6a7f99").attr("pointer-events","none")
        .text(g.count);
    });
    // black vertical separator lines at group boundaries
    hmGrpBoundaries.forEach(ci=>{
      const x=ci*cW+ROW_LABEL_W;
      svg.append("line")
        .attr("x1",x).attr("x2",x)
        .attr("y1",grpHdrY).attr("y2",CELLS_BOTTOM+HM_BAR_H+10)
        .attr("stroke","#1a1a1a").attr("stroke-width",1.5)
        .attr("pointer-events","none");
    });
  }

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
  const n=getEffectiveCustomOGs().length;
  const sumEl=document.getElementById("hm-custom-summary");
  if(sumEl) sumEl.textContent="Custom selection ("+n+" OG"+(n!==1?"s":"")+"):";
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
  return _loadEmbeddedJson("treedata-"+id);
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
    matching.sort((a,b)=>{ const na=parseInt((a.hg||a.id).match(/\d+/)?.[0]||'0',10), nb=parseInt((b.hg||b.id).match(/\d+/)?.[0]||'0',10); return na-nb; });
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
  // group option — only when SPECIES_GROUPS has data
  if(Object.keys(SPECIES_GROUPS).length>0){
    const oGrp=document.createElement("option"); oGrp.value="group"; oGrp.textContent="by group"; sel.appendChild(oGrp);
  }
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
  } else if(v==="group"){
    colorMode="group"; cladeSp2Color={}; cladeSp2Group={}; cladeGrpColor={}; ogLeaf2Color={}; ogName2Color={};
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
document.getElementById("sp-palette-select").addEventListener("change",function(){
  updateSpeciesPalettePreview(this.value);
});
document.getElementById("col-tri-fill").addEventListener("input",function(){
  colTriFill=this.value;
  // only gene-tree .col-tri and the species-tree view use this global colour
  document.querySelectorAll("#tree-svg .col-tri").forEach(el=>{ el.style.fill=colTriFill; });
  drawSpeciesTree();
});
updateSpeciesPalettePreview(document.getElementById("sp-palette-select").value);

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
  if(speciesSet.size===1) return leaves[0].data.name||null;
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

if (!HAVE_ALIGNMENTS) {
  const alignBtn = document.getElementById("tab-btn-align");
  if (alignBtn) alignBtn.style.display = "none";
}
if (!HAVE_ARCHITECTURES) {
  const archBtn = document.querySelector('.tab-btn[data-tab="architectures"]');
  if (archBtn) archBtn.style.display = "none";
}

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


def _load_tree_records(possvm_dir: Path, possvm_prev_dir):
    if not possvm_dir.exists():
        print(f"WARN: {possvm_dir} does not exist – no gene trees.", file=sys.stderr)

    raw_records, all_species, gene_meta = (
        load_possvm_trees(possvm_dir, source="generax") if possvm_dir.exists() else ([], [], {})
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

    records, all_species, prev_records, gene_meta, generax_gene_meta, prev_gene_meta = _load_tree_records(possvm_dir, possvm_prev_dir)
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
    domain_spans = build_domain_spans(domain_hits)
    gene_lengths = load_gene_lengths(cluster_dir)
    gene_to_hg, hg_sizes = load_gene_to_hg_map(cluster_dir)
    search_gene_lengths = load_search_gene_lengths(search_dir)
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
    protein_domain_catalog = build_exact_domain_catalog(search_dir, protein_lengths)
    architecture_catalog = build_domain_architecture_catalog(search_dir, protein_lengths, gene_to_hg, hg_sizes)

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
        HTML_TEMPLATE
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
