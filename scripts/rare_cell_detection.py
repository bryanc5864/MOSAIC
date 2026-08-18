#!/usr/bin/env python3
# MIT License
# Part of MOSAIC
"""Does cluster entropy penalize rare cell types just for being rare?

Rare populations (pDCs, Tregs, plasma cells) matter diagnostically, and a low
cluster density means fewer candidates in the OT plan, so the entropy could
flag them as uncertain even when they are present. This checks whether that
happens, by correlating cluster size against mean H_cluster across all three
datasets.
"""
from __future__ import annotations
import json, sys, os

import numpy as np
import anndata as ad
from scipy import stats

from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR

CONFIGS = [
    ("pbmc10k_multiome", "exp001_pbmc_beta0001", "PBMC 10k (beta=0.001)"),
    ("brain3k_multiome", "exp001_brain_beta0001", "Brain 5k (beta=0.001)"),
    ("citeseq_pbmc", "exp001_citeseq", "CITE-seq PBMC (beta=0.001)"),
]


def analyze_dataset(dataset_id, exp_id, label):
    exp_dir = EXPERIMENTS_DIR / exp_id
    H_path = exp_dir / "alignment_entropy_cluster.npy"
    sub_idx_path = exp_dir / "ot_subsample_indices.npy"

    if not H_path.exists():
        return {"label": label, "error": "H_cluster not found"}

    H = np.load(H_path)
    sub_idx = np.load(sub_idx_path) if sub_idx_path.exists() else np.arange(len(H))

    rna = ad.read_h5ad(PROCESSED_DIR / f"{dataset_id}_rna.h5ad")
    labels = rna.obs["leiden"].astype(str).values[sub_idx]

    clusters = np.unique(labels)
    per_cluster = []
    for c in clusters:
        mask = labels == c
        n = int(mask.sum())
        full_n = int((rna.obs["leiden"].astype(str).values == c).sum())
        if n < 3:
            continue
        per_cluster.append({
            "cluster": c,
            "n_in_subsample": n,
            "n_full": full_n,
            "mean_H": float(H[mask].mean()),
            "std_H": float(H[mask].std()),
        })

    per_cluster.sort(key=lambda x: x["n_full"])

    sizes = np.array([r["n_full"] for r in per_cluster])
    entropies = np.array([r["mean_H"] for r in per_cluster])
    rho, pval = stats.spearmanr(sizes, entropies)

    median_size = np.median(sizes)
    rare_mask = sizes < median_size
    rare_H = entropies[rare_mask]
    abundant_H = entropies[~rare_mask]
    stat, mw_pval = stats.mannwhitneyu(rare_H, abundant_H, alternative="greater") if (len(rare_H) > 1 and len(abundant_H) > 1) else (0, 1.0)

    return {
        "label": label,
        "dataset": dataset_id,
        "exp": exp_id,
        "n_clusters": len(per_cluster),
        "size_entropy_spearman_rho": float(rho),
        "size_entropy_spearman_p": float(pval),
        "rare_vs_abundant_H": {
            "rare_mean_H": float(rare_H.mean()) if len(rare_H) > 0 else float("nan"),
            "abundant_mean_H": float(abundant_H.mean()) if len(abundant_H) > 0 else float("nan"),
            "mannwhitney_p_rare_greater": float(mw_pval),
        },
        "per_cluster": per_cluster,
        "rare_clusters": [r for r in per_cluster if r["n_full"] < median_size],
        "abundant_clusters": [r for r in per_cluster if r["n_full"] >= median_size],
    }


def main():
    results = []
    for dataset_id, exp_id, label in CONFIGS:
        print(f"\n[dataset] {label}")
        r = analyze_dataset(dataset_id, exp_id, label)
        results.append(r)
        if "error" in r:
            print(f"  failed: {r['error']}")
            continue
        print(f"  Size-entropy Spearman rho={r['size_entropy_spearman_rho']:.3f}  p={r['size_entropy_spearman_p']:.2e}")
        print(f"  Rare clusters: mean H={r['rare_vs_abundant_H']['rare_mean_H']:.4f}")
        print(f"  Abundant clusters: mean H={r['rare_vs_abundant_H']['abundant_mean_H']:.4f}")
        print(f"  Mann-Whitney p(rare > abundant): {r['rare_vs_abundant_H']['mannwhitney_p_rare_greater']:.3e}")
        print("  Per-cluster (smallest first):")
        for c in r['per_cluster'][:5]:
            print(f"    cluster {c['cluster']:3s}  n_full={c['n_full']:5d}  mean_H={c['mean_H']:.4f}")

    output = {
        "experiment": "rare_cell_detection",
        "datasets": results,
        "summary": "Spearman correlation between cluster size and mean H_cluster (negative = smaller clusters have higher entropy)",
        "correlations": {r["label"]: r.get("size_entropy_spearman_rho", float("nan")) for r in results if "error" not in r},
    }

    os.makedirs("experiments/rare_cell_detection", exist_ok=True)
    with open("experiments/rare_cell_detection/results.json", "w") as f:
        json.dump(output, f, indent=2)
    print("\nsaved experiments/rare_cell_detection/results.json")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
