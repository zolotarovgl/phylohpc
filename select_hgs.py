#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Iterable, Optional, Tuple, Set


def iter_fasta_files(cluster_dir: Path) -> Iterable[Path]:
    exts = {".fa", ".fna", ".faa", ".fasta"}
    for p in cluster_dir.iterdir():
        if p.is_file() and p.suffix.lower() in exts:
            yield p


def analyze_fasta(
    fasta_path: Path, soi: Optional[str]
) -> Tuple[int, bool, int]:
    nseq = 0
    has_soi = False
    species: Set[str] = set()

    soi_prefix = f"{soi}_" if soi else None

    with fasta_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                nseq += 1
                header = line[1:].strip()

                if "_" in header:
                    sp = header.split("_", 1)[0]
                    species.add(sp)

                    if soi_prefix and header.startswith(soi_prefix):
                        has_soi = True

    return nseq, has_soi if soi else True, len(species)


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
    parser.add_argument("--soi", default=None)
    parser.add_argument("--min_seqs", type=int, default=1)
    parser.add_argument("--min_sps", type=int, default=1)
    parser.add_argument(
        "--out",
        default=None,
        help="Optional output TXT file (one ID per line).",
    )

    args = parser.parse_args()

    cluster_dir = Path(args.cluster_dir)
    if not cluster_dir.exists() or not cluster_dir.is_dir():
        raise SystemExit(f"ERROR: cluster_dir does not exist: {cluster_dir}")

    ids = []
    total_files = 0

    for fasta in sorted(iter_fasta_files(cluster_dir)):
        total_files += 1
        nseq, has_soi, n_species = analyze_fasta(fasta, args.soi)

        if nseq < args.min_seqs:
            continue
        if args.soi and not has_soi:
            continue
        if n_species < args.min_sps:
            continue

        ids.append(fasta.stem)

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