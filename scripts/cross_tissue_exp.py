#!/usr/bin/env python3
# MIT License
# Part of MOSAIC - Exp 6 cross-tissue negative-control experiment
"""Cross-tissue alignment as a negative control for the UQ signal.

Question: if we align MOSAIC embeddings from two totally disjoint tissues
(PBMC and brain — no shared cell types), does the cluster-resolved
entropy correctly flag every cell as high-uncertainty?

A well-calibrated UQ signal should say "I'm not confident" on every cell
when there's no real underlying structure to align. A miscalibrated
signal (like the naive per-row entropy we showed is wrong-sign) would
either pretend to be confident or have no informative pattern.

Protocol:
  1. Take PBMC 10k Multiome's trained IB-VAE RNA embedding (from
     exp001_pbmc_final).
  2. Take Brain 5k Multiome's trained IB-VAE ATAC embedding (from
     exp001_brain_beta0001).
  3. Fit a Procrustes rotation on the UNION cluster centroids, treating
     PBMC leiden clusters and Brain leiden clusters as the label space.
     (The cluster labels between datasets don't overlap; we use cluster
     IDs as opaque tokens.)
  4. Run Sinkhorn alignment.
  5. Compute cluster-resolved entropy.
  6. Compare the distribution of cluster entropies to the within-dataset
     PBMC and Brain distributions (which have low mean ~0.08-0.15).
     Prediction: the cross-tissue mean entropy should be substantially
     higher (cells are uncertain about cluster assignment because there's
     nothing meaningful to match to).

This is a negative control, not a positive result. The demonstration is
that the UQ signal correctly reports high uncertainty in a setting where
there's no ground truth to match.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.evaluation.metrics import foscttm, label_transfer_accuracy
from src.models.align_post import apply_alignment, fit_orthogonal_procrustes
from src.models.ot_align import sinkhorn_align
from src.utils.paths import EXPERIMENTS_DIR, FIGURES_DIR, PROCESSED_DIR


def cluster_marginal(plan, cluster_ids_target):
    rs = plan.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1e-30
    P = plan / rs
    unique = np.unique(cluster_ids_target)
    cluster_col = np.zeros((P.shape[1], len(unique)), dtype=np.float32)
    for k, c in enumerate(unique):
        cluster_col[cluster_ids_target == c, k] = 1.0
    cluster_marg = P @ cluster_col
    log_K = float(np.log(max(len(unique), 2)))
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(cluster_marg > 0, np.log(cluster_marg), 0.0)
    H = -(cluster_marg * logp).sum(axis=1) / log_K
    return H, cluster_marg, unique


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--rna-exp", default="exp001_pbmc_final",
                   help="Experiment whose z_rna to use")
    p.add_argument("--rna-dataset", default="pbmc10k_multiome")
    p.add_argument("--atac-exp", default="exp001_brain_beta0001",
                   help="Experiment whose z_atac_aligned to use")
    p.add_argument("--atac-dataset", default="brain3k_multiome")
    p.add_argument("--subsample", type=int, default=3000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    # Load RNA embeddings and labels
    z_rna_full = np.load(EXPERIMENTS_DIR / args.rna_exp / "z_rna.npy")
    rna_ad = ad.read_h5ad(PROCESSED_DIR / f"{args.rna_dataset}_rna.h5ad")
    labels_rna = rna_ad.obs["cell_type"].astype(str).values
    # Prepend "rna_" so PBMC labels are distinct from brain labels
    labels_rna_tagged = np.array([f"rna_{c}" for c in labels_rna])

    # Load ATAC embeddings and labels
    atac_z_path = EXPERIMENTS_DIR / args.atac_exp / "z_atac_aligned.npy"
    if not atac_z_path.exists():
        atac_z_path = EXPERIMENTS_DIR / args.atac_exp / "z_atac.npy"
    z_atac_full = np.load(atac_z_path)
    atac_ad = ad.read_h5ad(PROCESSED_DIR / f"{args.atac_dataset}_atac.h5ad")
    labels_atac = atac_ad.obs["cell_type"].astype(str).values
    labels_atac_tagged = np.array([f"atac_{c}" for c in labels_atac])

    print(f"RNA from {args.rna_dataset}: {z_rna_full.shape}, "
          f"{len(np.unique(labels_rna))} clusters")
    print(f"ATAC from {args.atac_dataset}: {z_atac_full.shape}, "
          f"{len(np.unique(labels_atac))} clusters")

    # Subsample each side
    rng = np.random.default_rng(args.seed)
    n_rna = min(args.subsample, z_rna_full.shape[0])
    n_atac = min(args.subsample, z_atac_full.shape[0])
    idx_rna = rng.choice(z_rna_full.shape[0], n_rna, replace=False)
    idx_atac = rng.choice(z_atac_full.shape[0], n_atac, replace=False)
    z_rna = z_rna_full[idx_rna]
    z_atac = z_atac_full[idx_atac]
    lbl_rna = labels_rna_tagged[idx_rna]
    lbl_atac = labels_atac_tagged[idx_atac]

    # NOTE: we don't try to fit Procrustes here because there are NO shared
    # cluster labels between PBMC and brain (they're tagged with different
    # prefixes and represent completely different cell types). The whole
    # point of this experiment is to see what happens when the two modalities
    # have NO meaningful alignment target. Procrustes is skipped; Sinkhorn
    # runs directly on the two latents.

    print(f"\nRunning Sinkhorn on {n_rna} RNA x {n_atac} ATAC cells...")
    align = sinkhorn_align(z_rna, z_atac, epsilon=0.05)
    print(f"  raw cell-level entropy: mean={align.entropy.mean():.4f}, "
          f"std={align.entropy.std():.4f}")

    # Cluster-resolved entropy — clusters come from the ATAC side since
    # that's what the plan marginalizes over.
    H_cluster, P_cluster, unique = cluster_marginal(align.plan, lbl_atac)
    print(f"  cluster-resolved entropy (over ATAC clusters): "
          f"mean={H_cluster.mean():.4f}, std={H_cluster.std():.4f}")

    # Compare to within-dataset baselines (loaded from the per-seed
    # cluster_entropy_analysis.json files).
    pbmc_cea_path = EXPERIMENTS_DIR / args.rna_exp / "cluster_entropy_analysis.json"
    atac_cea_path = EXPERIMENTS_DIR / args.atac_exp / "cluster_entropy_analysis.json"
    pbmc_within = None
    atac_within = None
    if pbmc_cea_path.exists():
        with pbmc_cea_path.open() as f:
            d = json.load(f)
        pbmc_within = d["cluster_level_entropy"]["mean"]
    if atac_cea_path.exists():
        with atac_cea_path.open() as f:
            d = json.load(f)
        atac_within = d["cluster_level_entropy"]["mean"]

    summary = {
        "experiment": "exp006_cross_tissue",
        "rna_source": args.rna_exp,
        "atac_source": args.atac_exp,
        "n_rna": n_rna,
        "n_atac": n_atac,
        "cross_tissue_cluster_entropy_mean": float(H_cluster.mean()),
        "cross_tissue_cluster_entropy_std": float(H_cluster.std()),
        "cross_tissue_cell_entropy_mean": float(align.entropy.mean()),
        "cross_tissue_cell_entropy_std": float(align.entropy.std()),
        "within_pbmc_cluster_entropy_mean": pbmc_within,
        "within_brain_cluster_entropy_mean": atac_within,
        "n_atac_clusters": int(len(unique)),
    }
    out_dir = EXPERIMENTS_DIR / "exp006_cross_tissue"
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "results.json").open("w") as f:
        json.dump(summary, f, indent=2)

    # Figure: histogram comparing within-PBMC, within-Brain, and
    # cross-tissue cluster entropy distributions.
    plt.rcParams.update({
        "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13,
        "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    fig, ax = plt.subplots(figsize=(9, 5))
    bins = np.linspace(0, 1, 40)
    if pbmc_within is not None:
        # We don't have full distribution for within-PBMC, just mean. Draw a
        # vertical line for the within means instead.
        ax.axvline(pbmc_within, color="#3770B0", linestyle="--", linewidth=2,
                   label=f"within PBMC mean H = {pbmc_within:.3f}")
    if atac_within is not None:
        ax.axvline(atac_within, color="#5AA06A", linestyle="--", linewidth=2,
                   label=f"within Brain mean H = {atac_within:.3f}")
    ax.hist(H_cluster, bins=bins, alpha=0.7, color="#C03030",
            label=f"cross-tissue H_cluster (n={n_rna}, mean {H_cluster.mean():.3f})",
            edgecolor="black", linewidth=0.3)
    ax.set_xlabel("cluster-resolved alignment entropy $H_{\\mathrm{cluster}}$")
    ax.set_ylabel("count")
    ax.set_title("Exp 6 — cross-tissue negative control\n"
                 "PBMC 10k RNA vs Brain 5k ATAC (no shared cell types)")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(axis="y", alpha=0.3)
    fig_path = FIGURES_DIR / "fig6_cross_tissue_negative_control.png"
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig_path)
    plt.savefig(fig_path.with_suffix(".pdf"))
    plt.close(fig)
    print(f"\n[fig] saved {fig_path}")

    print("\n=== Cross-tissue negative control summary ===")
    print(f"Cross-tissue H_cluster  : {H_cluster.mean():.4f} ± {H_cluster.std():.4f}")
    print(f"Within-PBMC H_cluster   : {pbmc_within:.4f}")
    print(f"Within-Brain H_cluster  : {atac_within:.4f}")
    print(f"Cross / within ratio    : {H_cluster.mean() / max(pbmc_within, atac_within, 1e-9):.2f}×")
    print(f"\nExpected: cross-tissue >> within (MOSAIC correctly signals no shared structure)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
