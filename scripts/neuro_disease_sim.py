#!/usr/bin/env python3
# MIT License
# Part of MOSAIC
"""Neurological disease simulation on the E18 brain benchmark.

Same leave-out machinery as the immunodeficiency sim, with brain leiden
clusters mapped onto neurodegeneration scenarios. Runs on the best brain
config (beta=0.001, seed 0).
"""
from __future__ import annotations
import json, sys

import numpy as np
import anndata as ad
from sklearn.metrics import roc_auc_score

from src.models.ot_align import sinkhorn_align
from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR

EXP_ID = "exp001_brain_beta0001"
DATASET = "brain3k_multiome"
SUBSAMPLE_N = 2500
EPSILON = 0.05
SEED = 0

def cluster_marginal(plan, cluster_ids_target):
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


def identify_cell_types(rna, Z_rna):
    """cluster sizes, biggest first"""
    labels = rna.obs["leiden"].astype(str).values
    sizes = {c: int((labels==c).sum()) for c in np.unique(labels)}
    return labels, sizes


def run_scenario(Z_rna, Z_atac, labels, target_clusters, name, description, rng):
    target_mask = np.array([l in set(target_clusters) for l in labels])
    n_target = int(target_mask.sum())
    if n_target < 10:
        return {"name": name, "n_target": n_target, "auroc": float("nan"), "note": "too few"}

    other_idx = np.where(~target_mask)[0]
    n_others = min(SUBSAMPLE_N - n_target, len(other_idx))
    other_sample = rng.choice(other_idx, size=n_others, replace=False)
    target_idx = np.where(target_mask)[0]
    rna_sub = np.concatenate([target_idx, other_sample])
    rng.shuffle(rna_sub)
    is_target = target_mask[rna_sub]
    Z_rna_sub = Z_rna[rna_sub]

    atac_keep = np.where(~target_mask)[0]
    atac_n = min(SUBSAMPLE_N, len(atac_keep))
    atac_sub = rng.choice(atac_keep, size=atac_n, replace=False)
    Z_atac_sub = Z_atac[atac_sub]
    labels_atac_sub = labels[atac_sub]

    print(f"  [{name}] RNA {len(rna_sub)} (target={n_target}) vs ATAC {atac_n}", flush=True)
    result = sinkhorn_align(Z_rna_sub, Z_atac_sub, epsilon=EPSILON, n_iter=200, normalize="median")
    H, _, _ = cluster_marginal(result.plan, labels_atac_sub)
    try:
        auroc = float(roc_auc_score(is_target.astype(int), H))
    except ValueError:
        auroc = float("nan")

    h_t = float(H[is_target].mean()) if is_target.any() else float("nan")
    h_o = float(H[~is_target].mean()) if (~is_target).any() else float("nan")
    return {
        "name": name,
        "description": description,
        "target_clusters": target_clusters,
        "n_target": n_target,
        "auroc": auroc,
        "mean_h_target": h_t,
        "mean_h_other": h_o,
        "entropy_ratio": h_t / h_o if h_o > 0 else float("nan"),
    }


def main():
    rng = np.random.default_rng(SEED)
    exp_dir = EXPERIMENTS_DIR / EXP_ID
    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac_path = exp_dir / "z_atac_aligned.npy"
    Z_atac = np.load(Z_atac_path if Z_atac_path.exists() else exp_dir / "z_atac.npy")
    rna = ad.read_h5ad(PROCESSED_DIR / f"{DATASET}_rna.h5ad")
    labels = rna.obs["leiden"].astype(str).values

    sizes = {c: int((labels==c).sum()) for c in np.unique(labels)}
    print("Brain clusters:", sorted(sizes.items(), key=lambda x: -x[1])[:10])

    # 20 clusters here; cell types are assigned by rank in cluster size, which
    # is crude but matches the expected E18 composition (neurons dominate,
    # glia are the small tail)
    sorted_clusters = sorted(sizes.items(), key=lambda x: -x[1])
    cluster_ids = [c for c, _ in sorted_clusters]

    scenarios = [
        {
            "name": "Excitatory_neuron_loss",
            "description": "Excitatory neuron loss (Alzheimer's disease / neurodegeneration)",
            "target_clusters": [cluster_ids[0], cluster_ids[1]],  # two largest, excitatory
        },
        {
            "name": "Inhibitory_neuron_loss",
            "description": "Inhibitory neuron loss (epilepsy / interneuron dysfunction)",
            "target_clusters": [cluster_ids[2], cluster_ids[3]],
        },
        {
            "name": "Oligodendrocyte_loss",
            "description": "Oligodendrocyte loss (multiple sclerosis / leukodystrophy)",
            "target_clusters": [cluster_ids[5], cluster_ids[6]],
        },
        {
            "name": "Astrocyte_loss",
            "description": "Astrocyte depletion (reactive astrogliosis / ALS context)",
            "target_clusters": [cluster_ids[4]],
        },
        {
            "name": "Microglia_depletion",
            "description": "Microglia depletion (PLX5622 CSF1R inhibition, preclinical model)",
            "target_clusters": [cluster_ids[-2], cluster_ids[-1]],  # smallest, microglia/endothelial
        },
    ]

    results = []
    for s in scenarios:
        print(f"\n[scenario] {s['name']}")
        r = run_scenario(Z_rna, Z_atac, labels,
                         s["target_clusters"], s["name"], s["description"], rng)
        print(f"  AUROC={r.get('auroc','nan'):.4f}  "
              f"H_ratio={r.get('entropy_ratio','nan'):.2f}x  n={r.get('n_target','?')}")
        results.append(r)

    valid = [r["auroc"] for r in results if not np.isnan(r.get("auroc", float("nan")))]
    output = {
        "experiment": "neuro_disease_sim",
        "dataset": DATASET,
        "source_exp": EXP_ID,
        "mean_auroc": float(np.mean(valid)),
        "min_auroc": float(min(valid)),
        "max_auroc": float(max(valid)),
        "scenarios": results,
    }

    import os; os.makedirs(f"experiments/neuro_disease_sim", exist_ok=True)
    with open("experiments/neuro_disease_sim/results.json", "w") as f:
        json.dump(output, f, indent=2)

    print("\n[summary] neuro disease scenarios")
    print(f"Mean AUROC: {output['mean_auroc']:.4f}")
    for r in results:
        print(f"  {r['name']}: AUROC={r.get('auroc','nan'):.4f}  "
              f"ratio={r.get('entropy_ratio','nan'):.2f}x")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
