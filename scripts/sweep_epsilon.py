#!/usr/bin/env python3
# MIT License
# Part of MOSAIC
"""Sweep Sinkhorn epsilon over the saved embeddings of an existing experiment.

Useful because epsilon does not affect the IB-VAE training or the embeddings;
only the OT plan and per-cell entropy depend on it. We can thus retroactively
sweep epsilon on a single trained run without retraining.

Usage:
    python -m scripts.sweep_epsilon --exp run004_procrustes
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np

from src.evaluation.metrics import entropy_error_corr, foscttm
from src.models.ot_align import sinkhorn_align
from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--exp", required=True)
    p.add_argument("--dataset", default="pbmc10k_multiome")
    p.add_argument("--epsilons", nargs="*", type=float,
                   default=[0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1])
    args = p.parse_args()

    exp_dir = EXPERIMENTS_DIR / args.exp
    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac_aligned_path = exp_dir / "z_atac_aligned.npy"
    if Z_atac_aligned_path.exists():
        Z_atac = np.load(Z_atac_aligned_path)
    else:
        Z_atac = np.load(exp_dir / "z_atac.npy")
    print(f"loaded Z_rna {Z_rna.shape}, Z_atac {Z_atac.shape}")

    # Pair indices for entropy calibration
    rna = ad.read_h5ad(PROCESSED_DIR / f"{args.dataset}_rna.h5ad")
    atac = ad.read_h5ad(PROCESSED_DIR / f"{args.dataset}_atac.h5ad")
    pair_a = rna.obs["pair_idx"].values.astype(np.int64)
    pair_b = atac.obs["pair_idx"].values.astype(np.int64)

    # Diagnostic: distribution of cost matrix entries
    from src.models.ot_align import normalize_cost, pairwise_sqeuclidean
    C = pairwise_sqeuclidean(Z_rna, Z_atac)
    print(f"\ncost matrix raw: shape {C.shape}, "
          f"min {C.min():.3f}, mean {C.mean():.3f}, max {C.max():.3f}, "
          f"std {C.std():.3f}")
    Cn = normalize_cost(C)
    print(f"cost matrix normalized to [0,1]: "
          f"min {Cn.min():.3f}, mean {Cn.mean():.3f}, std {Cn.std():.3f}")

    # Per-row diagnostics on unnormalized cost
    row_min = C.min(axis=1)
    row_max = C.max(axis=1)
    row_range = row_max - row_min
    print(f"row range (max-min) on raw cost: "
          f"mean {row_range.mean():.3f}, std {row_range.std():.3f}")

    # FOSCTTM is constant across epsilon -- compute once.
    f = foscttm(Z_rna, Z_atac, pair_a, pair_b)
    print(f"\nFOSCTTM (constant across epsilon): {f['foscttm_mean']:.4f}")

    # Sweep
    print("\nepsilon  | mean H | std H  | top1 prob | corr  rho | n_iter ok")
    print("---------+--------+--------+-----------+-----------+----------")
    rows = []
    for eps in args.epsilons:
        try:
            res = sinkhorn_align(Z_rna, Z_atac, epsilon=eps)
            ent_corr = entropy_error_corr(res.entropy, Z_rna, Z_atac, pair_a, pair_b)
            row = {
                "epsilon": eps,
                "mean_entropy": float(res.entropy.mean()),
                "std_entropy": float(res.entropy.std()),
                "mean_top1_prob": float(res.top_match_prob.mean()),
                "spearman_rho": ent_corr["spearman_rho"],
                "ok": True,
            }
            print(f" {eps:.4f}  | {row['mean_entropy']:.4f} | "
                  f"{row['std_entropy']:.4f} | {row['mean_top1_prob']:.4f}    | "
                  f"{row['spearman_rho']:+.4f}   | yes")
        except Exception as e:
            row = {"epsilon": eps, "ok": False, "error": str(e)}
            print(f" {eps:.4f}  | FAIL: {e}")
        rows.append(row)

    out_path = exp_dir / "epsilon_sweep.json"
    with out_path.open("w") as f:
        json.dump({"sweep": rows, "foscttm": float(foscttm(Z_rna, Z_atac, pair_a, pair_b)["foscttm_mean"])},
                  f, indent=2)
    print(f"\nsaved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
