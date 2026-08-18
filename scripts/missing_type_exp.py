#!/usr/bin/env python3
# MIT License
"""Leave-one-cluster-out: does cluster entropy spot a missing cell type?

For each candidate cluster, drop all its ATAC cells, rerun Sinkhorn against the
full RNA side, and score AUROC for "this RNA cell has no partner left" against
cluster entropy. Cells of the removed type should spread their transport mass
over the wrong clusters and come out high-entropy.

Repeated one cluster at a time; we report the AUROC distribution.
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
    """cluster entropy per source cell plus the (n_src, n_clusters) marginal"""
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
    """Drop target_cluster from ATAC, realign, score detection AUROC.

    RNA stays complete; positive class is "this RNA cell is of the removed
    type". The subsample keeps every cell of the removed type (otherwise
    there's nothing to detect) and fills up to subsample_n at random.
    """
    rng = np.random.default_rng(seed)

    target_mask = labels == target_cluster
    n_target = int(target_mask.sum())
    if n_target < 5:
        return {"target_cluster": target_cluster, "n_target": n_target,
                "auroc": float("nan"), "note": "too few target cells"}

    other_indices = np.where(~target_mask)[0]
    n_others = min(subsample_n - n_target, len(other_indices))
    other_sample = rng.choice(other_indices, size=n_others, replace=False)
    target_indices = np.where(target_mask)[0]
    rna_sub = np.concatenate([target_indices, other_sample])
    rng.shuffle(rna_sub)
    is_target_rna = labels[rna_sub] == target_cluster

    Z_rna_sub = Z_rna[rna_sub]

    # atac loses the target cluster entirely
    atac_keep_mask = labels != target_cluster
    atac_kept_indices = np.where(atac_keep_mask)[0]
    atac_sample_n = min(subsample_n, len(atac_kept_indices))
    atac_sub_indices = rng.choice(atac_kept_indices, size=atac_sample_n, replace=False)
    Z_atac_sub = Z_atac[atac_sub_indices]
    labels_atac_sub = labels[atac_sub_indices]

    print(f"  [leave-out {target_cluster}] RNA {len(rna_sub)} (target={n_target}) "
          f"vs ATAC {len(atac_sub_indices)} (target removed)")

    result = sinkhorn_align(Z_rna_sub, Z_atac_sub, epsilon=epsilon,
                            n_iter=200, normalize="median")

    H_cluster, _, _ = cluster_marginal(result.plan, labels_atac_sub)

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

    unique_clusters = sorted(np.unique(labels), key=lambda c: int(c))
    sizes = {c: int((labels == c).sum()) for c in unique_clusters}
    # too few cells and the auroc is noise
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

    print(f"\n[summary] missing-type, {args.dataset}")
    print(f"clusters tested: {aggregate['n_clusters_tested']}")
    print(f"mean AUROC: {aggregate['mean_auroc']:.4f}")
    print(f"median AUROC: {aggregate['median_auroc']:.4f}")
    print(f"range: [{aggregate['min_auroc']:.4f}, {aggregate['max_auroc']:.4f}]")
    print(f"saved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
