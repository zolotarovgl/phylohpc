#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Set


def iter_fasta_files(cluster_dir: Path) -> Iterable[Path]:
    exts = {".fa", ".fna", ".faa", ".fasta"}
    for p in cluster_dir.iterdir():
        if p.is_file() and p.suffix.lower() in exts:
            yield p


def analyze_fasta(
    fasta_path: Path, sois: List[str], soi_mode: str = "any"
) -> Tuple[int, bool, int, Dict[str, int], Set[str]]:
    """sois: species prefixes (empty = no soi filter). soi_mode 'any' passes an HG holding
    at least one soi, 'all' only one holding every soi. Species = text before the first '_'."""
    nseq = 0
    n_soi: Dict[str, int] = {s: 0 for s in sois}
    species: Set[str] = set()

    with fasta_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                nseq += 1
                header = line[1:].strip()

                if "_" in header:
                    sp = header.split("_", 1)[0]
                    species.add(sp)
                    if sp in n_soi:
                        n_soi[sp] += 1

    if not sois:
        has_soi = True
    elif soi_mode == "all":
        has_soi = all(n > 0 for n in n_soi.values())
    else:
        has_soi = any(n > 0 for n in n_soi.values())
    return nseq, has_soi, len(species), n_soi, species


def split_hg_id(stem: str) -> Tuple[str, str, str]:
    """'neu.Syntaxin.HG3' -> ('neu', 'Syntaxin', 'HG3'); anything else -> ('', stem, '')."""
    parts = stem.split(".")
    if len(parts) >= 3 and parts[-1].startswith("HG"):
        return parts[0], ".".join(parts[1:-1]), parts[-1]
    return "", stem, ""


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "List FASTA basenames (without extension) from cluster_dir. "
            "Optional filters: --soi, --min_seqs, --min_sps."
        )
    )
    parser.add_argument(
        "-d",
        "--cluster_dir",
        default="results/clusters",
        help="Directory containing FASTA files (default: results/clusters)",
    )
    parser.add_argument(
        "--soi",
        default=None,
        help="Species of interest, comma-separated species prefixes (e.g. Mlei or Mlei,Plebac).",
    )
    parser.add_argument(
        "--soi_mode",
        choices=["any", "all"],
        default="any",
        help="With several --soi: keep an HG holding ANY of them (default) or ALL of them.",
    )
    parser.add_argument("--min_seqs", type=int, default=1)
    parser.add_argument("--min_sps", type=int, default=1)
    parser.add_argument(
        "--out",
        default=None,
        help="Optional output TXT file (one ID per line).",
    )
    parser.add_argument(
        "--meta",
        default=None,
        help=(
            "Optional output TSV with one row per HG (ALL HGs, passed or not): "
            "hg, pref, family, hg_num, n_seqs, n_species, n_soi, passed, fail_reason, species."
        ),
    )

    args = parser.parse_args()
    sois = [s.strip() for s in args.soi.split(",") if s.strip()] if args.soi else []

    cluster_dir = Path(args.cluster_dir)
    if not cluster_dir.exists() or not cluster_dir.is_dir():
        raise SystemExit(f"ERROR: cluster_dir does not exist: {cluster_dir}")

    ids = []
    meta_rows = []
    total_files = 0

    for fasta in sorted(iter_fasta_files(cluster_dir)):
        total_files += 1
        nseq, has_soi, n_species, n_soi, species = analyze_fasta(fasta, sois, args.soi_mode)

        # every failed filter is recorded, not just the first
        reasons = []
        if nseq < args.min_seqs:
            reasons.append(f"n_seqs<{args.min_seqs}")
        if sois and not has_soi:
            missing = [s for s in sois if n_soi[s] == 0]
            reasons.append(f"no_{'+'.join(missing)}" if args.soi_mode == "all" else f"no_{'|'.join(sois)}")
        if n_species < args.min_sps:
            reasons.append(f"n_species<{args.min_sps}")

        pref, family, hg_num = split_hg_id(fasta.stem)
        meta_rows.append([
            fasta.stem, pref, family, hg_num, str(nseq), str(n_species),
            ",".join(f"{s}:{n_soi[s]}" for s in sois) if sois else "NA",
            "TRUE" if not reasons else "FALSE",
            ",".join(reasons), ",".join(sorted(species)),
        ])

        if reasons:
            continue

        ids.append(fasta.stem)

    if args.meta:
        meta_path = Path(args.meta)
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        header = ["hg", "pref", "family", "hg_num", "n_seqs", "n_species",
                  "n_soi", "passed", "fail_reason", "species"]
        with meta_path.open("w", encoding="utf-8") as fh:
            fh.write("\t".join(header) + "\n")
            for row in meta_rows:
                fh.write("\t".join(row) + "\n")
        # stdout carries the id list when --out is unset, so only report in --out mode
        if args.out:
            print(f"HG metadata: {len(meta_rows)} rows -> {meta_path}", flush=True)

    output_text = "\n".join(ids)

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(output_text + "\n", encoding="utf-8")

        # Reporting
        print(
            f"Filtering complete: {len(ids)} / {total_files} HGs passed.",
            flush=True,
        )
    else:
        print(output_text)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())