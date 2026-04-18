#!/usr/bin/env python3
# MIT License
# Part of MOSAIC - Exp 2 analysis
"""Cluster-resolved entropy analysis for MOSAIC's per-cell uncertainty.

Rationale: the per-cell row entropy on the entropic-OT transport plan
marginally captures "how many cells in the target modality are similar to
this source cell". That's confounded with cluster density: dense clusters
have high row entropy AND small distance-to-true-partner, giving a WRONG-
signed calibration.

The fix: for each source cell i, compute the CLUSTER-LEVEL marginal of
the transport plan:
    p(cluster k | i) = sum_{j in cluster k} P_ij   (after row normalization)
and the cluster-level entropy:
    H_cluster(i) = -sum_k p(cluster k | i) log p(cluster k | i)
                   / log(n_clusters)

Low H_cluster = confident about which cluster the partner is in (good
alignment). High H_cluster = spread over clusters (uncertain — should
correlate with genuine alignment failures).

This script:
  1. Loads an experiment's OT plan and subsample indices.
  2. Computes cluster-level entropy using the paired dataset's leiden labels.
  3. Compares cell-level and cluster-level entropy correlations with the
     true alignment distance.
  4. Reports Spearman rho for both measures, saves a JSON summary.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
from scipy.stats import spearmanr

from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR


def cluster_resolved_entropy(plan: np.ndarray, cluster_ids_target: np.ndarray
                              ) -> np.ndarray:
    """Given a transport plan (n_src, n_tgt) and a cluster label per target
    cell, return the per-source-cell cluster entropy in [0, 1].
    """
    # Row-normalize
    row_sums = plan.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1e-30
    P = plan / row_sums

    unique_clusters = np.unique(cluster_ids_target)
    n_clusters = len(unique_clusters)
    cluster_col = np.zeros((P.shape[1], n_clusters), dtype=np.float32)
    for k, c in enumerate(unique_clusters):
        cluster_col[cluster_ids_target == c, k] = 1.0
    # p(cluster k | i) = sum_{j in cluster k} P[i, j]
    P_cluster = P @ cluster_col               # (n_src, n_clusters)
    log_K = float(np.log(max(n_clusters, 2)))
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(P_cluster > 0, np.log(P_cluster), 0.0)
    H = -(P_cluster * logp).sum(axis=1) / log_K
    return H, P_cluster, unique_clusters


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--exp", required=True)
    p.add_argument("--dataset", default="pbmc10k_multiome")
    args = p.parse_args()

    exp_dir = EXPERIMENTS_DIR / args.exp
    plan = np.load(exp_dir / "alignment_plan_subsample.npy")
    sub_idx = np.load(exp_dir / "ot_subsample_indices.npy")
    entropy_cell = np.load(exp_dir / "alignment_entropy_subsample.npy")
    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac_path = exp_dir / "z_atac_aligned.npy"
    Z_atac = np.load(Z_atac_path if Z_atac_path.exists() else exp_dir / "z_atac.npy")

    rna = ad.read_h5ad(PROCESSED_DIR / f"{args.dataset}_rna.h5ad")
    atac = ad.read_h5ad(PROCESSED_DIR / f"{args.dataset}_atac.h5ad")
    labels = rna.obs["cell_type"].astype(str).values  # same for both (paired)

    # Subsample cluster labels (both modalities share indices in sub_idx)
    labels_sub = labels[sub_idx]
    Z_rna_sub = Z_rna[sub_idx]
    Z_atac_sub = Z_atac[sub_idx]

    # Cluster-resolved entropy
    H_cluster, P_cluster, uniq = cluster_resolved_entropy(plan, labels_sub)

    # True alignment error (subsample)
    errs = np.sqrt(((Z_rna_sub - Z_atac_sub) ** 2).sum(axis=1))

    rho_cell, p_cell = spearmanr(entropy_cell, errs)
    rho_cluster, p_cluster_val = spearmanr(H_cluster, errs)

    # Also: is the argmax cluster correct?
    argmax_cluster_idx = P_cluster.argmax(axis=1)
    argmax_cluster = uniq[argmax_cluster_idx]
    cluster_correct = (argmax_cluster == labels_sub)
    argmax_accuracy = float(cluster_correct.mean())

    # Key UQ test: does cluster entropy predict cluster-assignment errors?
    # AUROC for (is_wrong_cluster, cluster_entropy) — a well-calibrated
    # uncertainty measure should give AUROC >> 0.5 on this task.
    from sklearn.metrics import roc_auc_score
    is_wrong = (~cluster_correct).astype(int)
    if is_wrong.sum() >= 2 and is_wrong.sum() < len(is_wrong) - 1:
        auroc_entropy_wrong = float(roc_auc_score(is_wrong, H_cluster))
    else:
        auroc_entropy_wrong = float("nan")

    # Also: is the TOP-1 cell (not cluster) in the correct cluster?
    top_cell_idx = plan.argmax(axis=1)  # n_src per-cell argmax
    top_cell_cluster = labels_sub[top_cell_idx]
    top_cell_correct_cluster = (top_cell_cluster == labels_sub)
    top_cell_cluster_accuracy = float(top_cell_correct_cluster.mean())

    # Per-cluster stats — how many cells of each cluster were wrong?
    per_cluster_wrong_rate = {}
    for c in uniq:
        mask = labels_sub == c
        if mask.sum() > 0:
            per_cluster_wrong_rate[str(c)] = {
                "n": int(mask.sum()),
                "argmax_wrong_rate": float(is_wrong[mask].mean()),
                "mean_cluster_entropy": float(H_cluster[mask].mean()),
            }

    summary = {
        "exp": args.exp,
        "dataset": args.dataset,
        "n_subsample": int(len(sub_idx)),
        "n_clusters": int(len(uniq)),
        "cell_level_entropy": {
            "mean": float(entropy_cell.mean()),
            "std": float(entropy_cell.std()),
            "spearman_rho_vs_error": float(rho_cell),
            "spearman_p": float(p_cell),
        },
        "cluster_level_entropy": {
            "mean": float(H_cluster.mean()),
            "std": float(H_cluster.std()),
            "spearman_rho_vs_error": float(rho_cluster),
            "spearman_p": float(p_cluster_val),
            "argmax_cluster_accuracy": argmax_accuracy,
            "auroc_entropy_vs_wrong_cluster": auroc_entropy_wrong,
            "n_wrong_cluster_cells": int(is_wrong.sum()),
        },
        "top_cell_argmax": {
            "top1_cell_in_correct_cluster_rate": top_cell_cluster_accuracy,
        },
        "per_cluster_wrong_rate": per_cluster_wrong_rate,
    }

    out_path = exp_dir / "cluster_entropy_analysis.json"
    with out_path.open("w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))

    # Save per-cell cluster entropy and argmax cluster for downstream plots
    np.save(exp_dir / "alignment_entropy_cluster.npy", H_cluster.astype(np.float32))
    np.save(exp_dir / "alignment_plan_cluster_marginal.npy", P_cluster.astype(np.float32))

    print(f"\nsaved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
