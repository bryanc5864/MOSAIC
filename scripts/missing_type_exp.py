#!/usr/bin/env python3
# MIT License
# Part of MOSAIC - Exp 3: missing cell type detection
"""Exp 3 - can MOSAIC's cluster-resolved entropy detect cells whose true
cluster is absent from the target modality?

Protocol:
  1. Load the trained IB-VAE embeddings from exp001_pbmc_final.
  2. Pick a target cluster k*.
  3. Remove all ATAC cells in cluster k* (simulate absence in target modality).
  4. Re-run Sinkhorn on the reduced set (RNA full, ATAC minus k*).
  5. Compute cluster-resolved entropy.
  6. Measure AUROC for the task:
         "is this RNA cell of the removed type" vs. "cluster entropy"
     Expectation: RNA cells of the removed type have higher cluster entropy
     (they have no good partner in ATAC, the OT plan spreads mass across
     wrong clusters). AUROC should be well above 0.5.

For a fair test, we repeat for each candidate cluster (one-at-a-time leave-out)
and report the AUROC distribution.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
from sklearn.metrics import roc_auc_score

from src.models.ot_align import sinkhorn_align
from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR


def cluster_marginal(plan: np.ndarray, cluster_ids_target: np.ndarray
                      ) -> tuple[np.ndarray, np.ndarray]:
    """Returns cluster entropy per source cell, and the per-cluster marginal
    probability (n_src x n_clusters)."""
    rs = plan.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1e-30
    P = plan / rs
    unique = np.unique(cluster_ids_target)
    idx_by_c = [np.where(cluster_ids_target == c)[0] for c in unique]
    cluster_marg = np.stack([P[:, idx].sum(axis=1) for idx in idx_by_c], axis=1)
    log_K = float(np.log(max(len(unique), 2)))
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(cluster_marg > 0, np.log(cluster_marg), 0.0)
    H = -(cluster_marg * logp).sum(axis=1) / log_K
    return H, cluster_marg, unique


def run_leave_out_cluster(
        Z_rna: np.ndarray, Z_atac: np.ndarray,
        labels: np.ndarray, target_cluster: str,
        epsilon: float = 0.05, subsample_n: int = 3000,
        seed: int = 0,
) -> dict:
    """Remove `target_cluster` from ATAC, run Sinkhorn, compute cluster
    entropy. Return AUROC for detecting the removed type via entropy.

    We keep the RNA side complete; the ATAC side has the target cluster's
    cells removed. For the AUROC, positive class = "this RNA cell is of
    the removed type" (no partner in ATAC), negative = any other cell.

    The subsample restricts the expensive Sinkhorn to `subsample_n` cells
    from the RNA side. We pick the subsample to include EVERY cell of the
    removed type (else there's nothing to detect) plus a random sample
    from the other clusters to reach `subsample_n`.
    """
    rng = np.random.default_rng(seed)

    # Cells in target cluster (to detect)
    target_mask = labels == target_cluster
    n_target = int(target_mask.sum())
    if n_target < 5:
        return {"target_cluster": target_cluster, "n_target": n_target,
                "auroc": float("nan"), "note": "too few target cells"}

    # RNA subsample: all target cells + random sample of others
    other_indices = np.where(~target_mask)[0]
    n_others = min(subsample_n - n_target, len(other_indices))
    other_sample = rng.choice(other_indices, size=n_others, replace=False)
    target_indices = np.where(target_mask)[0]
    rna_sub = np.concatenate([target_indices, other_sample])
    rng.shuffle(rna_sub)
    is_target_rna = labels[rna_sub] == target_cluster

    Z_rna_sub = Z_rna[rna_sub]

    # ATAC: keep everything EXCEPT the target cluster
    atac_keep_mask = labels != target_cluster
    # Also subsample ATAC to roughly the same size
    atac_kept_indices = np.where(atac_keep_mask)[0]
    atac_sample_n = min(subsample_n, len(atac_kept_indices))
    atac_sub_indices = rng.choice(atac_kept_indices, size=atac_sample_n, replace=False)
    Z_atac_sub = Z_atac[atac_sub_indices]
    labels_atac_sub = labels[atac_sub_indices]

    print(f"  [leave-out {target_cluster}] RNA {len(rna_sub)} (target={n_target}) "
          f"vs ATAC {len(atac_sub_indices)} (target removed)")

    # Run Sinkhorn
    result = sinkhorn_align(Z_rna_sub, Z_atac_sub, epsilon=epsilon,
                            n_iter=200, normalize="median")

    # Cluster-resolved entropy
    H_cluster, _, _ = cluster_marginal(result.plan, labels_atac_sub)

    # AUROC for detecting the removed type via entropy
    try:
        auroc = float(roc_auc_score(is_target_rna.astype(int), H_cluster))
    except ValueError:
        auroc = float("nan")

    return {
        "target_cluster": str(target_cluster),
        "n_target_rna": int(n_target),
        "n_rna_sub": int(len(rna_sub)),
        "n_atac_sub_after_removal": int(len(atac_sub_indices)),
        "auroc_cluster_entropy": auroc,
        "mean_entropy_target": float(H_cluster[is_target_rna].mean()) if is_target_rna.any() else float("nan"),
        "mean_entropy_other": float(H_cluster[~is_target_rna].mean()) if (~is_target_rna).any() else float("nan"),
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--exp", required=True)
    p.add_argument("--dataset", default="pbmc10k_multiome")
    p.add_argument("--subsample-n", type=int, default=2500)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    exp_dir = EXPERIMENTS_DIR / args.exp
    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac_path = exp_dir / "z_atac_aligned.npy"
    Z_atac = np.load(Z_atac_path if Z_atac_path.exists() else exp_dir / "z_atac.npy")

    rna = ad.read_h5ad(PROCESSED_DIR / f"{args.dataset}_rna.h5ad")
    labels = rna.obs["cell_type"].astype(str).values

    # Per-cluster stats to pick candidates
    unique_clusters = sorted(np.unique(labels), key=lambda c: int(c))
    sizes = {c: int((labels == c).sum()) for c in unique_clusters}
    # Only test clusters with enough cells to get a meaningful AUROC
    candidates = [c for c in unique_clusters if 30 <= sizes[c] <= 1500]
    print(f"dataset {args.dataset}: {len(unique_clusters)} clusters, "
          f"testing {len(candidates)} leave-out experiments")

    results = []
    for c in candidates:
        print(f"\n-- leave-out cluster {c} (n={sizes[c]}) --")
        r = run_leave_out_cluster(Z_rna, Z_atac, labels, c,
                                    subsample_n=args.subsample_n, seed=args.seed)
        print(f"   AUROC: {r.get('auroc_cluster_entropy', 'nan')}")
        print(f"   mean H target: {r.get('mean_entropy_target', 'nan'):.3f}, "
              f"other: {r.get('mean_entropy_other', 'nan'):.3f}")
        results.append(r)

    # Aggregate
    valid = [r["auroc_cluster_entropy"] for r in results
             if r.get("auroc_cluster_entropy") is not None
             and not np.isnan(r["auroc_cluster_entropy"])]
    aggregate = {
        "dataset": args.dataset,
        "exp": args.exp,
        "n_clusters_tested": len(valid),
        "mean_auroc": float(np.mean(valid)) if valid else float("nan"),
        "median_auroc": float(np.median(valid)) if valid else float("nan"),
        "min_auroc": float(min(valid)) if valid else float("nan"),
        "max_auroc": float(max(valid)) if valid else float("nan"),
        "per_cluster": results,
    }

    out_path = exp_dir / "exp003_missing_type.json"
    with out_path.open("w") as f:
        json.dump(aggregate, f, indent=2)

    print(f"\n=== Exp 3 summary ({args.dataset}) ===")
    print(f"clusters tested: {aggregate['n_clusters_tested']}")
    print(f"mean AUROC: {aggregate['mean_auroc']:.4f}")
    print(f"median AUROC: {aggregate['median_auroc']:.4f}")
    print(f"range: [{aggregate['min_auroc']:.4f}, {aggregate['max_auroc']:.4f}]")
    print(f"saved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
