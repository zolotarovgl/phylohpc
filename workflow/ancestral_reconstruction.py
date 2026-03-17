#!/usr/bin/env python3
"""
Ancestral state reconstruction using PastML (MPPA + F81 model).

PastML implements the MPPA (Marginal Posterior Probabilities Approximation)
method, which computes P(state | all tip data) at every internal node using
a Bayesian integration over the model parameters.  The F81 model allows
asymmetric gain/loss rates via estimated stationary frequencies, making it
equivalent to a two-rate binary Mk model.

Reference:
    Ishikawa SA et al. (2019) A fast likelihood method to reconstruct and
    visualize ancestral scenarios. Mol Biol Evol 36(9):2069–2085.
    doi:10.1093/molbev/msz131

Input:
    --pam        Binary presence/absence matrix TSV (rows=species, cols=OGs)
    --tree       Pruned species tree (newick, with named internal nodes)

Output:
    --output     Per-OG summary: P(present at root/LCA), support label (TSV)
    --node_probs Per-OG × per-internal-node marginal probabilities (TSV)
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ete3 import Tree


# ── Support classification ────────────────────────────────────────────────────

def classify_support(p: float) -> str:
    if np.isnan(p):
        return "uncertain"
    if p >= 0.9:
        return "present"
    if p >= 0.5:
        return "likely_present"
    if p >= 0.1:
        return "likely_absent"
    return "absent"


# ── Tree helpers ──────────────────────────────────────────────────────────────

def get_tree_info(tree_path: str) -> tuple:
    """Return (root_name, internal_node_names_set, leaf_names_set)."""
    t = Tree(tree_path, format=1)
    root_name = t.name if t.name else None
    internal  = {n.name for n in t.traverse() if not n.is_leaf() and n.name}
    leaves    = set(t.get_leaf_names())
    return root_name, internal, leaves


# ── PastML runner ─────────────────────────────────────────────────────────────

def run_pastml(tree_path: str, data_path: str, work_dir: Path) -> None:
    """Run PastML MPPA+F81.  Raises RuntimeError on non-zero exit."""
    cmd = [
        "pastml",
        "--tree",              tree_path,
        "--data",              str(data_path),
        "--model",             "F81",
        "--prediction_method", "MPPA",
        "--work_dir",          str(work_dir),
    ]
    print(f"Running PastML: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"PastML exited with code {result.returncode}.\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )


# ── Output parser ─────────────────────────────────────────────────────────────

def parse_marginal_file(prob_path: Path) -> pd.Series:
    """
    Parse a PastML marginal_probabilities.{column}.tab file.

    The file has rows = nodes, columns = state labels ('0', '1').
    Returns a float Series indexed by node name with P(state='1') values.
    """
    try:
        df = pd.read_csv(prob_path, sep="\t", index_col=0)
    except Exception as exc:
        print(f"  WARNING: Cannot parse {prob_path}: {exc}", file=sys.stderr)
        return pd.Series(dtype=float)

    # Locate the column for state "1" (present)
    col1 = next((c for c in df.columns if str(c).strip() == "1"), None)
    if col1 is None:
        print(f"  WARNING: State '1' column not found in {prob_path}", file=sys.stderr)
        return pd.Series(dtype=float)

    return df[col1].rename(index=str)


from typing import Optional

def get_root_prob(
    probs: pd.Series,
    root_name: Optional[str],
    internal_nodes: set,
    leaf_names: set,
    og: str,
) -> float:
    """
    Extract P(present) at the root from a marginal probability Series.
    Falls back to searching for non-leaf nodes when root_name is absent.
    """
    if probs.empty:
        return float("nan")

    # Preferred: direct name lookup
    if root_name and root_name in probs.index:
        return float(probs[root_name])

    # Fallback: any internal-node entry not in the leaf set
    candidates = [n for n in probs.index if n not in leaf_names]
    if not candidates:
        print(f"  WARNING: No internal nodes found in PastML output for {og}", file=sys.stderr)
        return float("nan")

    if len(candidates) > 1:
        print(
            f"  WARNING: Root node '{root_name}' not found for {og}; "
            f"multiple candidates: {candidates[:4]}. "
            "Name the root in the pruned tree for unambiguous results.",
            file=sys.stderr,
        )

    # Return the first candidate (order in PastML output is typically top-down)
    return float(probs[candidates[0]])


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ancestral state reconstruction using PastML (MPPA + F81)."
    )
    parser.add_argument("--pam",        required=True, help="PAM TSV (rows=species, cols=OGs)")
    parser.add_argument("--tree",       required=True, help="Pruned species tree (newick)")
    parser.add_argument("--output",     required=True, help="Per-OG ancestral states TSV")
    parser.add_argument("--node_probs", required=True, help="Per-OG × per-node probabilities TSV")
    parser.add_argument("--node",       default="root", help="Clade name for labelling")
    args = parser.parse_args()

    # ── Load PAM ──────────────────────────────────────────────────────────────
    pam = pd.read_csv(args.pam, sep="\t", index_col=0)
    print(f"PAM: {pam.shape[0]} species × {pam.shape[1]} OGs", file=sys.stderr)

    empty_cols = ["og", "n_present", "n_total", "P_at_root", "support"]

    if pam.empty or pam.shape[1] == 0:
        pd.DataFrame(columns=empty_cols).to_csv(args.output, sep="\t", index=False)
        pd.DataFrame(columns=["og", "node", "P_present"]).to_csv(
            args.node_probs, sep="\t", index=False
        )
        return

    root_name, internal_nodes, leaf_names = get_tree_info(args.tree)
    print(
        f"Tree root: '{root_name}' | internal nodes: {len(internal_nodes)} "
        f"| leaves: {len(leaf_names)}",
        file=sys.stderr,
    )

    # ── Prepare PastML input ──────────────────────────────────────────────────
    # Rename OG columns to safe aliases (og0, og1, …) to avoid issues with
    # dots/special characters in OG names being mangled in output filenames.
    og_list   = pam.columns.tolist()
    col_alias = {og: f"og{i}" for i, og in enumerate(og_list)}
    alias_og  = {v: k for k, v in col_alias.items()}

    work_dir  = Path("pastml_work")
    work_dir.mkdir(exist_ok=True)

    pastml_data = pam.rename(columns=col_alias).astype(str)
    pastml_data.index.name = "node"
    data_path = work_dir / "input.tsv"
    pastml_data.to_csv(data_path, sep="\t")

    # ── Run PastML ────────────────────────────────────────────────────────────
    try:
        run_pastml(args.tree, str(data_path), work_dir)
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)

    # ── Parse outputs ─────────────────────────────────────────────────────────
    results   = []
    node_rows = []

    for alias, og in alias_og.items():
        n_present = int(pam[og].sum())
        n_total   = int(pam[og].notna().sum())

        prob_path = work_dir / f"marginal_probabilities.{alias}.tab"
        probs     = parse_marginal_file(prob_path) if prob_path.exists() else pd.Series(dtype=float)

        p_root = get_root_prob(probs, root_name, internal_nodes, leaf_names, og)

        results.append(dict(
            og=og,
            n_present=n_present,
            n_total=n_total,
            P_at_root=round(p_root, 4) if not np.isnan(p_root) else float("nan"),
            support=classify_support(p_root),
        ))

        # Store internal-node probabilities only (skip leaf rows)
        for node_id, p1 in probs.items():
            if node_id not in leaf_names:
                node_rows.append(dict(og=og, node=node_id, P_present=round(float(p1), 4)))

    # ── Write outputs ─────────────────────────────────────────────────────────
    pd.DataFrame(results).to_csv(args.output, sep="\t", index=False)
    pd.DataFrame(node_rows).to_csv(args.node_probs, sep="\t", index=False)

    print(f"Written: {args.output}", file=sys.stderr)
    print(f"Written: {args.node_probs}", file=sys.stderr)

    # Summary
    results_df = pd.DataFrame(results)
    if not results_df.empty:
        counts = results_df["support"].value_counts()
        print("\nSupport summary:", file=sys.stderr)
        for label, n in counts.items():
            print(f"  {label:<20} {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
