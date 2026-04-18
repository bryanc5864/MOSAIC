#!/usr/bin/env python3
# MIT License - Bryan Cheng, 2026
# Part of MOSAIC
"""Clinical disease simulation experiment.

Maps PBMC leiden clusters to clinical immune cell types and simulates
immunodeficiency / disease states by removing entire clinical populations
from the target (ATAC) modality. Tests whether MOSAIC's cluster-resolved
entropy flags RNA cells whose immune population is absent from the reference.

Disease scenarios modelled:
  CD8 T lymphopenia    — clusters 0, 14 (CD8A/B high, NKG7, GNLY)
  NK cell deficiency   — cluster 4 (NCAM1/CD56 high, GNLY high)
  B cell aplasia       — clusters 11, 15 (MS4A1/CD20 high)
  Monocytopenia        — cluster 6 (CD14 high, LYZ high)
  Treg deficiency      — cluster 5 (FOXP3 high, IL2RA high)

Each scenario: remove the target population from ATAC, rerun Sinkhorn on
full RNA vs. reduced ATAC, compute cluster-resolved entropy, report AUROC
for discriminating removed cells by entropy score.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
from sklearn.metrics import roc_auc_score

from src.models.ot_align import sinkhorn_align
from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR

# ── Clinical population definitions ──────────────────────────────────────────
CLINICAL_SCENARIOS = [
    {
        "name": "CD8_T_lymphopenia",
        "description": "CD8 T cell depletion (post-transplant immunosuppression / viral-induced)",
        "clusters": ["0", "14"],
        "markers": "CD8A+, CD8B+, NKG7+, GNLY+",
    },
    {
        "name": "NK_cell_deficiency",
        "description": "NK cell deficiency (chronic viral infection immunosuppression)",
        "clusters": ["4"],
        "markers": "NCAM1/CD56+, NKG7+, GNLY+, FCGR3A+",
    },
    {
        "name": "B_cell_aplasia",
        "description": "B cell aplasia (post-rituximab therapy / X-linked agammaglobulinemia)",
        "clusters": ["11", "15"],
        "markers": "MS4A1/CD20+",
    },
    {
        "name": "Monocytopenia",
        "description": "Monocytopenia (hairy cell leukemia / bone marrow failure)",
        "clusters": ["6"],
        "markers": "CD14+, LYZ+, CST3+",
    },
    {
        "name": "Treg_deficiency",
        "description": "Regulatory T cell deficiency (autoimmune disease context)",
        "clusters": ["5"],
        "markers": "FOXP3+, IL2RA/CD25+",
    },
]

EXP_ID = "exp001_pbmc_beta0001"   # use 10-seed best: seed 0
DATASET = "pbmc10k_multiome"
SUBSAMPLE_N = 3000
EPSILON = 0.05
SEED = 0


def cluster_marginal(plan: np.ndarray, cluster_ids_target: np.ndarray
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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


def run_clinical_scenario(
        Z_rna: np.ndarray,
        Z_atac: np.ndarray,
        labels: np.ndarray,
        scenario: dict,
        rng: np.random.Generator,
) -> dict:
    target_clusters = set(scenario["clusters"])
    target_mask = np.array([l in target_clusters for l in labels])
    n_target = int(target_mask.sum())

    if n_target < 10:
        return {**scenario, "n_target": n_target, "auroc": float("nan"),
                "note": "too few target cells"}

    # RNA subsample: ALL target cells + random others to reach SUBSAMPLE_N
    other_idx = np.where(~target_mask)[0]
    n_others = min(SUBSAMPLE_N - n_target, len(other_idx))
    other_sample = rng.choice(other_idx, size=n_others, replace=False)
    target_idx = np.where(target_mask)[0]
    rna_sub = np.concatenate([target_idx, other_sample])
    rng.shuffle(rna_sub)

    is_target_rna = target_mask[rna_sub]
    Z_rna_sub = Z_rna[rna_sub]

    # ATAC: remove the entire target clinical population
    atac_keep = np.where(~target_mask)[0]
    atac_n = min(SUBSAMPLE_N, len(atac_keep))
    atac_sub = rng.choice(atac_keep, size=atac_n, replace=False)
    Z_atac_sub = Z_atac[atac_sub]
    labels_atac_sub = labels[atac_sub]

    n_remaining_clusters = len(np.unique(labels_atac_sub))
    print(f"  [{scenario['name']}] RNA {len(rna_sub)} (target={n_target}) "
          f"vs ATAC {atac_n} ({n_remaining_clusters} clusters remain)")

    result = sinkhorn_align(Z_rna_sub, Z_atac_sub, epsilon=EPSILON,
                            n_iter=200, normalize="median")

    H, _, _ = cluster_marginal(result.plan, labels_atac_sub)

    try:
        auroc = float(roc_auc_score(is_target_rna.astype(int), H))
    except ValueError:
        auroc = float("nan")

    h_target = float(H[is_target_rna].mean()) if is_target_rna.any() else float("nan")
    h_other = float(H[~is_target_rna].mean()) if (~is_target_rna).any() else float("nan")

    return {
        **scenario,
        "n_target": n_target,
        "n_rna_sub": int(len(rna_sub)),
        "n_atac_sub": int(atac_n),
        "n_remaining_clusters_atac": n_remaining_clusters,
        "auroc": auroc,
        "mean_h_target": h_target,
        "mean_h_other": h_other,
        "entropy_ratio": h_target / h_other if h_other > 0 else float("nan"),
    }


def main() -> int:
    rng = np.random.default_rng(SEED)

    exp_dir = EXPERIMENTS_DIR / EXP_ID
    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac_path = exp_dir / "z_atac_aligned.npy"
    Z_atac = np.load(Z_atac_path if Z_atac_path.exists() else exp_dir / "z_atac.npy")

    rna = ad.read_h5ad(PROCESSED_DIR / f"{DATASET}_rna.h5ad")
    labels = rna.obs["leiden"].astype(str).values

    print(f"Dataset: {DATASET}, N={len(labels)}, clusters={len(np.unique(labels))}")
    print(f"Using embeddings from: {exp_dir}\n")

    results = []
    for scenario in CLINICAL_SCENARIOS:
        print(f"\n=== {scenario['name']} ===")
        print(f"  Clinical context: {scenario['description']}")
        print(f"  Removing ATAC clusters: {scenario['clusters']}")
        r = run_clinical_scenario(Z_rna, Z_atac, labels, scenario, rng)
        print(f"  AUROC: {r.get('auroc', 'nan'):.4f}")
        print(f"  Mean H (target): {r.get('mean_h_target', 'nan'):.4f}  "
              f"Mean H (other): {r.get('mean_h_other', 'nan'):.4f}  "
              f"ratio: {r.get('entropy_ratio', 'nan'):.2f}x")
        results.append(r)

    valid_aurocs = [r["auroc"] for r in results
                    if r.get("auroc") is not None and not np.isnan(r["auroc"])]

    output = {
        "experiment": "clinical_disease_sim",
        "dataset": DATASET,
        "source_exp": EXP_ID,
        "n_scenarios": len(CLINICAL_SCENARIOS),
        "mean_auroc": float(np.mean(valid_aurocs)) if valid_aurocs else float("nan"),
        "min_auroc": float(min(valid_aurocs)) if valid_aurocs else float("nan"),
        "max_auroc": float(max(valid_aurocs)) if valid_aurocs else float("nan"),
        "scenarios": results,
    }

    out_dir = EXPERIMENTS_DIR / "clinical_disease_sim"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "results.json"
    with out_path.open("w") as f:
        json.dump(output, f, indent=2)

    print(f"\n=== Clinical Disease Simulation Summary ===")
    print(f"Scenarios: {len(valid_aurocs)} valid")
    print(f"Mean AUROC: {output['mean_auroc']:.4f}")
    print(f"Range: [{output['min_auroc']:.4f}, {output['max_auroc']:.4f}]")
    for r in results:
        print(f"  {r['name']}: AUROC={r.get('auroc', 'nan'):.4f} "
              f"H_ratio={r.get('entropy_ratio', 'nan'):.2f}x "
              f"(n={r.get('n_target', '?')})")
    print(f"\nSaved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
