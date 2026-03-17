#!/usr/bin/env python3
"""
Ancestral state reconstruction for binary (presence/absence) orthogroup data
using the Mk model (continuous-time Markov chain on a binary character).

The model:
  States : 0 = absent,  1 = present
  Rate matrix:
        Q = [[-q01,  q01],
             [ q10, -q10]]
  where q01 = gain rate and q10 = loss rate.

Both rates are estimated jointly by maximum likelihood using the Felsenstein
pruning algorithm, then marginal ancestral probabilities P(state=1) are
computed at every internal node via the two-pass (pruning + peeling) algorithm.

Outputs:
  --output      : per-OG summary (root / LCA probability, rate estimates)
  --node_probs  : per-OG × per-internal-node probability table

Usage example:
  python ancestral_reconstruction.py \\
      --pam  Bilateria.pam.tsv \\
      --tree Bilateria.pruned.tree \\
      --output      Bilateria.ancestral_states.tsv \\
      --node_probs  Bilateria.node_probs.tsv \\
      --node        Bilateria
"""

import argparse
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from ete3 import Tree

warnings.filterwarnings("ignore")

# ── Mk model: transition probability matrix ───────────────────────────────────

def mk_transition(q01: float, q10: float, t: float) -> np.ndarray:
    """
    Return 2×2 transition probability matrix P(t) = exp(Q*t) for the
    two-rate binary Mk model.

    P[s_from, s_to]
    """
    r = q01 + q10
    if r < 1e-15 or t < 1e-15:
        return np.eye(2)

    e   = np.exp(-r * t)
    p01 = q01 / r * (1.0 - e)   # P(0 → 1)
    p10 = q10 / r * (1.0 - e)   # P(1 → 0)

    return np.array([[1.0 - p01, p01],
                     [p10,       1.0 - p10]])


# ── Felsenstein pruning (bottom-up conditional likelihoods) ───────────────────

def felsenstein_pruning(
    tree: Tree,
    tip_states: dict[str, int],
    q01: float,
    q10: float,
) -> dict[int, np.ndarray]:
    """
    Compute conditional likelihood vectors L[node] = [P(data | state=0),
    P(data | state=1)] for every node using the postorder (Felsenstein)
    pruning algorithm.

    tip_states : {leaf_name: 0 or 1}   (missing leaves → ambiguous)
    Returns    : {id(node): array([L0, L1])}
    """
    condL: dict[int, np.ndarray] = {}

    for node in tree.traverse("postorder"):
        nid = id(node)

        if node.is_leaf():
            state = tip_states.get(node.name, -1)
            if state == 0:
                condL[nid] = np.array([1.0, 0.0])
            elif state == 1:
                condL[nid] = np.array([0.0, 1.0])
            else:
                condL[nid] = np.array([1.0, 1.0])   # ambiguous / missing
        else:
            L = np.ones(2)
            for child in node.children:
                t        = max(child.dist, 1e-6)
                P        = mk_transition(q01, q10, t)   # shape (2,2)
                child_L  = condL[id(child)]
                # P[s_parent, s_child] summed over s_child
                L       *= P @ child_L
            condL[nid] = L

    return condL


# ── Negative log-likelihood ───────────────────────────────────────────────────

def neg_log_likelihood(
    log_params: np.ndarray,
    tree: Tree,
    tip_states: dict[str, int],
) -> float:
    q01, q10 = np.exp(log_params)
    condL    = felsenstein_pruning(tree, tip_states, q01, q10)

    r = q01 + q10
    pi = np.array([q10 / r, q01 / r]) if r > 1e-15 else np.array([0.5, 0.5])

    root_L = condL[id(tree)]
    lnL    = np.log(max(np.dot(pi, root_L), 1e-300))
    return -lnL


# ── ML optimisation ───────────────────────────────────────────────────────────

def fit_mk_model(
    tree: Tree,
    tip_states: dict[str, int],
    n_restarts: int = 5,
) -> tuple[float, float, float]:
    """
    Estimate q01 and q10 by maximum likelihood.

    Returns (q01, q10, lnL).
    """
    rng = np.random.default_rng(42)

    # Starting points: one near typical values + random restarts
    starts = [np.log([0.1, 0.1])] + [
        rng.uniform(-4.0, 1.0, size=2) for _ in range(n_restarts - 1)
    ]

    best_nll = np.inf
    best_res = None

    for p0 in starts:
        try:
            res = minimize(
                neg_log_likelihood,
                p0,
                args=(tree, tip_states),
                method="Nelder-Mead",
                options={"maxiter": 10_000, "xatol": 1e-7, "fatol": 1e-7},
            )
            if res.fun < best_nll:
                best_nll = res.fun
                best_res = res
        except Exception:
            pass

    if best_res is None:
        return 0.1, 0.1, float("nan")

    q01, q10 = np.exp(best_res.x)
    return float(q01), float(q10), float(-best_nll)


# ── Marginal ancestral reconstruction (two-pass) ─────────────────────────────

def marginal_ancestral_states(
    tree: Tree,
    tip_states: dict[str, int],
    q01: float,
    q10: float,
) -> dict[str, float]:
    """
    Compute marginal P(state=1) at each internal node using the
    two-pass (pruning + peeling) algorithm.

    Returns {node_label: P(present)}.
    Node labels are the node's .name attribute or a synthetic
    'node_<postorder_index>' for unnamed nodes.
    """
    condL = felsenstein_pruning(tree, tip_states, q01, q10)

    r  = q01 + q10
    pi = np.array([q10 / r, q01 / r]) if r > 1e-15 else np.array([0.5, 0.5])

    # margL[node] ≈ joint P(data, state) — not fully normalised but
    # proportional within each node.
    margL: dict[int, np.ndarray] = {}

    for node in tree.traverse("preorder"):
        nid = id(node)

        if node.is_root():
            margL[nid] = pi * condL[nid]
        else:
            parent   = node.up
            pid      = id(parent)
            t        = max(node.dist, 1e-6)
            P        = mk_transition(q01, q10, t)

            # Partial parent marginal: parent's margL divided by this child's
            # contribution to avoid double-counting.
            child_contrib = P @ condL[nid]   # shape (2,)
            parent_partial = np.where(
                child_contrib > 1e-300,
                margL[pid] / child_contrib,
                margL[pid],
            )

            # node_margL[s_child] ∝ condL[child][s_child]
            #                        × Σ_{s_parent} P[s_parent,s_child] × parent_partial[s_parent]
            node_margL = condL[nid] * (P.T @ parent_partial)
            margL[nid] = node_margL

    # Collect P(present) for each internal node
    result: dict[str, float] = {}
    counter = 0

    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        nid   = id(node)
        total = margL[nid].sum()

        if total > 1e-300:
            p1 = float(margL[nid][1] / total)
        else:
            p1 = float("nan")

        # Assign a stable label
        if node.name:
            label = node.name
        elif node.is_root():
            label = "root"
        else:
            label = f"node_{counter}"
        counter += 1

        result[label] = p1

    return result


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


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mk-model ancestral state reconstruction for binary OG presence/absence data."
    )
    parser.add_argument("--pam",        required=True, help="PAM TSV (rows=species, cols=OGs)")
    parser.add_argument("--tree",       required=True, help="Pruned species tree (newick)")
    parser.add_argument("--output",     required=True, help="Per-OG ancestral states at root (TSV)")
    parser.add_argument("--node_probs", required=True, help="Per-OG × per-node probabilities (TSV)")
    parser.add_argument("--node",       default="root", help="Clade name used for labelling")
    parser.add_argument(
        "--n_restarts", type=int, default=5,
        help="Number of ML optimisation restarts (default: 5)"
    )
    args = parser.parse_args()

    # ── Load data ─────────────────────────────────────────────────────────────
    pam = pd.read_csv(args.pam, sep="\t", index_col=0)
    print(f"PAM: {pam.shape[0]} species × {pam.shape[1]} OGs", file=sys.stderr)

    if pam.empty or pam.shape[1] == 0:
        print("WARNING: Empty PAM — nothing to reconstruct.", file=sys.stderr)
        pd.DataFrame(
            columns=["og", "n_present", "n_total", "q01", "q10", "lnL", "P_at_root", "support"]
        ).to_csv(args.output, sep="\t", index=False)
        pd.DataFrame(
            columns=["og", "node", "P_present"]
        ).to_csv(args.node_probs, sep="\t", index=False)
        return

    tree = Tree(args.tree, format=1)

    # Ensure non-zero branch lengths
    for node in tree.traverse():
        if node.dist == 0.0 and not node.is_root():
            node.dist = 1.0

    tree_species = set(tree.get_leaf_names())
    pam_species  = set(pam.index)

    missing_in_pam  = tree_species - pam_species
    missing_in_tree = pam_species  - tree_species

    if missing_in_pam:
        print(
            f"  {len(missing_in_pam)} species in tree but absent from PAM "
            "(treated as ambiguous): " + ", ".join(sorted(missing_in_pam)[:6]),
            file=sys.stderr,
        )
    if missing_in_tree:
        print(
            f"  {len(missing_in_tree)} species in PAM but absent from tree "
            "(ignored): " + ", ".join(sorted(missing_in_tree)[:6]),
            file=sys.stderr,
        )

    # ── Process each OG ───────────────────────────────────────────────────────
    og_list       = pam.columns.tolist()
    results       = []
    node_rows     = []

    for idx, og in enumerate(og_list, start=1):
        if idx % 100 == 0 or idx == 1:
            print(f"  OG {idx}/{len(og_list)}: {og}", file=sys.stderr)

        # Build tip_states
        tip_states: dict[str, int] = {}
        for sps in pam.index:
            if sps in tree_species:
                tip_states[sps] = int(pam.loc[sps, og])

        n_present = sum(v == 1 for v in tip_states.values())
        n_total   = len(tip_states)

        # ── Shortcut: trivial cases ───────────────────────────────────────────
        if n_present == 0:
            results.append(dict(
                og=og, n_present=0, n_total=n_total,
                q01=float("nan"), q10=float("nan"), lnL=float("nan"),
                P_at_root=0.0, support="absent",
            ))
            continue

        if n_present == n_total:
            results.append(dict(
                og=og, n_present=n_present, n_total=n_total,
                q01=float("nan"), q10=float("nan"), lnL=float("nan"),
                P_at_root=1.0, support="present",
            ))
            continue

        # ── ML fit ────────────────────────────────────────────────────────────
        try:
            q01, q10, lnL = fit_mk_model(tree, tip_states, n_restarts=args.n_restarts)

            node_p1 = marginal_ancestral_states(tree, tip_states, q01, q10)

            # Root probability: first entry in dict is the root (preorder)
            root_key = next(iter(node_p1))
            p_root   = node_p1[root_key]

            results.append(dict(
                og=og,
                n_present=n_present,
                n_total=n_total,
                q01=round(q01, 6),
                q10=round(q10, 6),
                lnL=round(lnL, 4),
                P_at_root=round(p_root, 4),
                support=classify_support(p_root),
            ))

            for node_label, p1 in node_p1.items():
                node_rows.append(dict(og=og, node=node_label, P_present=round(p1, 4)))

        except Exception as exc:
            print(f"  WARNING: OG {og} failed — {exc}", file=sys.stderr)
            results.append(dict(
                og=og, n_present=n_present, n_total=n_total,
                q01=float("nan"), q10=float("nan"), lnL=float("nan"),
                P_at_root=float("nan"), support="failed",
            ))

    # ── Write outputs ──────────────────────────────────────────────────────────
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output, sep="\t", index=False)
    print(f"Written: {args.output}", file=sys.stderr)

    node_df = pd.DataFrame(node_rows)
    node_df.to_csv(args.node_probs, sep="\t", index=False)
    print(f"Written: {args.node_probs}", file=sys.stderr)

    # ── Summary ───────────────────────────────────────────────────────────────
    if not results_df.empty:
        counts = results_df["support"].value_counts()
        print("\nSupport summary:", file=sys.stderr)
        for label, n in counts.items():
            print(f"  {label:<20} {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
