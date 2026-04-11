#!/usr/bin/env python3
# MIT License - Bryan Cheng, 2026
# Part of MOSAIC — Exp 2: entropy calibration figure
"""Generate the entropy-calibration figure from a completed experiment.

For each source cell, plots:
  - Scatter of per-cell alignment entropy vs. distance to true partner.
  - Binned calibration curve (mean error per entropy decile).

Inputs: one experiment directory containing alignment_entropy.npy,
        z_rna.npy, z_atac_aligned.npy (or z_atac.npy), and the processed
        dataset to load pair_idx.

Output: figures/exp002_calibration_<dataset>.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from src.utils.paths import EXPERIMENTS_DIR, FIGURES_DIR, PROCESSED_DIR


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--exp", required=True, help="experiment name")
    p.add_argument("--dataset", default="pbmc10k_multiome")
    p.add_argument("--out", default=None)
    args = p.parse_args()

    exp_dir = EXPERIMENTS_DIR / args.exp
    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac_aligned_path = exp_dir / "z_atac_aligned.npy"
    if Z_atac_aligned_path.exists():
        Z_atac = np.load(Z_atac_aligned_path)
    else:
        Z_atac = np.load(exp_dir / "z_atac.npy")
    entropy = np.load(exp_dir / "alignment_entropy.npy")

    rna = ad.read_h5ad(PROCESSED_DIR / f"{args.dataset}_rna.h5ad")
    atac = ad.read_h5ad(PROCESSED_DIR / f"{args.dataset}_atac.h5ad")
    pair_a = rna.obs["pair_idx"].values.astype(np.int64)
    pair_b = atac.obs["pair_idx"].values.astype(np.int64)
    labels = rna.obs["cell_type"].astype(str).values
    b_lookup = {int(k): i for i, k in enumerate(pair_b)}

    # For each cell, compute distance to its true partner in aligned space.
    errs = np.full(len(Z_rna), np.nan)
    for i, k in enumerate(pair_a):
        j = b_lookup.get(int(k))
        if j is None:
            continue
        errs[i] = float(np.sqrt(((Z_rna[i] - Z_atac[j]) ** 2).sum()))

    keep = ~np.isnan(errs)
    H = entropy[keep]
    E = errs[keep]
    L = labels[keep]

    rho, pval = spearmanr(H, E)

    # Binned calibration: 10 equal-count entropy deciles
    order = np.argsort(H)
    H_sorted = H[order]
    E_sorted = E[order]
    n_bins = 10
    bin_size = len(H) // n_bins
    bin_H = np.zeros(n_bins)
    bin_E = np.zeros(n_bins)
    bin_Estd = np.zeros(n_bins)
    for bi in range(n_bins):
        lo = bi * bin_size
        hi = (bi + 1) * bin_size if bi < n_bins - 1 else len(H)
        bin_H[bi] = H_sorted[lo:hi].mean()
        bin_E[bi] = E_sorted[lo:hi].mean()
        bin_Estd[bi] = E_sorted[lo:hi].std()

    # Figure
    plt.rcParams.update({
        "font.size": 12, "axes.labelsize": 14, "axes.titlesize": 14,
        "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10,
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8))

    # LEFT: scatter colored by cluster (gray dots with alpha)
    ax1.scatter(H, E, s=3, alpha=0.25, c="#3770B0", edgecolors="none")
    ax1.set_xlabel("normalized alignment entropy H")
    ax1.set_ylabel("distance to true partner (aligned latent)")
    ax1.set_title(f"per-cell calibration\nSpearman rho = {rho:.3f} (p={pval:.1e}), n={len(H)}")
    ax1.grid(alpha=0.3)

    # RIGHT: binned calibration curve with error bars
    ax2.errorbar(bin_H, bin_E, yerr=bin_Estd, marker="o", lw=1.5, capsize=3,
                 color="#C03030")
    ax2.set_xlabel("mean entropy per decile")
    ax2.set_ylabel("mean distance to true partner")
    ax2.set_title("binned calibration (10 entropy deciles)")
    ax2.grid(alpha=0.3)

    plt.suptitle(f"MOSAIC entropy calibration - {args.dataset} / {args.exp}",
                 y=1.02)
    plt.tight_layout()

    out_path = (FIGURES_DIR / f"exp002_calibration_{args.dataset}_{args.exp}.png") \
        if args.out is None else Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)
    plt.close(fig)
    print(f"saved {out_path}")
    print(f"Spearman rho = {rho:.4f}, p = {pval:.2e}, n = {len(H)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
