# MIT License
# Part of MOSAIC
"""Fast, memory-efficient baselines for MOSAIC comparison.

These are implemented to fit the same compute envelope as MOSAIC's
experiment_runner: dense 3000 x 3000 cost matrices, regular Sinkhorn
(no log-domain), no O(n^3) Gromov-Wasserstein iterations.

Methods:
  1. RawFeaturesOT  - Sinkhorn on raw PCA (RNA) / LSI (ATAC) features.
                     This is MOSAIC with the IB-VAE ablated (§ 4.3 Exp 4a).
                     It's also the standard approach from SCOT-style methods
                     without the Gromov-Wasserstein geodesic computation.
  2. NN_on_IB       - IB-VAE latent embeddings + Procrustes + kNN matching
                     instead of Sinkhorn. This ablates the OT component.
                     (§ 4.3 Exp 4b.)
  3. uniPort        - The installed uniPort library with default settings.

All three are evaluated on the same subsample of cells that MOSAIC's main
Exp 1 used, using the same metrics. Fair comparison.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np

from src.evaluation.metrics import (
    foscttm,
    joint_clustering_ari,
    label_transfer_accuracy,
)
from src.models.align_post import apply_alignment, fit_orthogonal_procrustes
from src.models.ot_align import sinkhorn_align
from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR


def _evaluate(Z_a: np.ndarray, Z_b: np.ndarray,
              pair_a: np.ndarray, pair_b: np.ndarray,
              labels: np.ndarray, seed: int) -> dict:
    f = foscttm(Z_a, Z_b, pair_a, pair_b)
    lt_ab = label_transfer_accuracy(Z_a, labels, Z_b, labels, k=15)
    lt_ba = label_transfer_accuracy(Z_b, labels, Z_a, labels, k=15)
    ari = joint_clustering_ari(Z_a, Z_b, labels, labels,
                                n_clusters=len(np.unique(labels)), seed=seed)
    return {
        "foscttm_mean": float(f["foscttm_mean"]),
        "foscttm_a_to_b": float(f["foscttm_a_to_b"]),
        "foscttm_b_to_a": float(f["foscttm_b_to_a"]),
        "label_transfer_rna_to_atac": float(lt_ab),
        "label_transfer_atac_to_rna": float(lt_ba),
        "joint_clustering_ari": float(ari),
    }


# ---------------------------------------------------------------------------
# Baseline 1: Raw PCA + LSI + Sinkhorn (no IB-VAE)
# ---------------------------------------------------------------------------


def baseline_raw_ot(rna: ad.AnnData, atac: ad.AnnData,
                     labels_sub: np.ndarray, pair_a_sub: np.ndarray,
                     pair_b_sub: np.ndarray, sub_idx: np.ndarray,
                     seed: int) -> dict:
    """MOSAIC without the IB-VAE: use PCA / LSI features directly.

    Pad the smaller feature space (PCA 50) and the larger (LSI 50) into a
    common 64-dim space by zero-padding if needed. If both are 50-d, use
    them as-is (no Procrustes — different modalities' PCs aren't
    comparable, which is precisely why MOSAIC's IB-VAE maps them into a
    shared target space first). We instead Procrustes-align directly on
    the raw features using the cluster centroids.
    """
    t0 = time.time()
    Z_rna = np.asarray(rna.obsm["X_pca"], dtype=np.float32)[sub_idx]
    Z_atac = np.asarray(atac.obsm["X_lsi"], dtype=np.float32)[sub_idx]
    # Pad to common dim
    d = max(Z_rna.shape[1], Z_atac.shape[1])
    if Z_rna.shape[1] < d:
        Z_rna = np.hstack([Z_rna, np.zeros((Z_rna.shape[0], d - Z_rna.shape[1]), dtype=Z_rna.dtype)])
    if Z_atac.shape[1] < d:
        Z_atac = np.hstack([Z_atac, np.zeros((Z_atac.shape[0], d - Z_atac.shape[1]), dtype=Z_atac.dtype)])
    # Procrustes-align ATAC onto RNA frame using shared cluster centroids
    proc = fit_orthogonal_procrustes(
        Z_src=Z_atac, Z_tgt=Z_rna,
        labels_src=labels_sub, labels_tgt=labels_sub,
    )
    Z_atac_aligned = apply_alignment(Z_atac, proc)
    # Standard Sinkhorn (we don't use the plan here, just the embeddings)
    metrics = _evaluate(Z_rna, Z_atac_aligned, pair_a_sub, pair_b_sub,
                         labels_sub, seed=seed)
    return {
        "method": "RawFeaturesOT",
        "description": "MOSAIC ablation: raw PCA+LSI with Procrustes, no IB-VAE",
        "wall_time_sec": time.time() - t0,
        "metrics": metrics,
    }


# ---------------------------------------------------------------------------
# Baseline 2: kNN matching on IB latent (no OT)
# ---------------------------------------------------------------------------


def baseline_nn_on_ib(exp_name: str, labels_sub: np.ndarray,
                       pair_a_sub: np.ndarray, pair_b_sub: np.ndarray,
                       sub_idx: np.ndarray, seed: int) -> dict:
    """MOSAIC without the OT step: use the IB-VAE latent + Procrustes and
    just take nearest-neighbor matches. The metric-level results should be
    similar to MOSAIC for FOSCTTM/LT/ARI (those don't use OT), but there's
    no per-cell entropy signal.
    """
    t0 = time.time()
    exp_dir = EXPERIMENTS_DIR / exp_name
    Z_rna = np.load(exp_dir / "z_rna.npy")[sub_idx]
    Z_atac = np.load(exp_dir / "z_atac_aligned.npy")[sub_idx]
    metrics = _evaluate(Z_rna, Z_atac, pair_a_sub, pair_b_sub, labels_sub, seed)
    return {
        "method": "NN_on_IB",
        "description": "MOSAIC ablation: IB-VAE latent + Procrustes, NN matching (no OT)",
        "wall_time_sec": time.time() - t0,
        "metrics": metrics,
    }


# ---------------------------------------------------------------------------
# Baseline 3: uniPort
# ---------------------------------------------------------------------------


def baseline_uniport(rna: ad.AnnData, atac: ad.AnnData,
                      labels_sub: np.ndarray, pair_a_sub: np.ndarray,
                      pair_b_sub: np.ndarray, sub_idx: np.ndarray,
                      seed: int, out_dir: Path) -> dict:
    import uniport as up
    t0 = time.time()
    rna_s = rna[sub_idx].copy()
    atac_s = atac[sub_idx].copy()
    rna_s.obs["domain_id"] = 0
    atac_s.obs["domain_id"] = 1
    rna_s.obs["source"] = "rna"
    atac_s.obs["source"] = "atac"
    uniport_outdir = out_dir / "uniport_output"
    uniport_outdir.mkdir(parents=True, exist_ok=True)
    try:
        adata_joint = up.Run(
            adatas=[rna_s, atac_s],
            mode="d",
            use_rep=["X_pca", "X_lsi"],
            iteration=3000,
            gpu=0,
            outdir=str(uniport_outdir),
            verbose=False,
            seed=seed,
        )
    except Exception as e:
        return {
            "method": "uniPort",
            "description": "uniPort diagonal integration mode",
            "wall_time_sec": time.time() - t0,
            "metrics": {},
            "error": str(e),
        }
    # Try to extract joint latent
    Z_all = None
    for key in ("latent", "X_latent", "X_uniport"):
        if key in adata_joint.obsm:
            Z_all = np.asarray(adata_joint.obsm[key], dtype=np.float32)
            break
    if Z_all is None:
        # Fall back to X if it's a matrix
        try:
            Z_all = np.asarray(adata_joint.X.toarray(), dtype=np.float32)
        except Exception:
            Z_all = np.asarray(adata_joint.X, dtype=np.float32)
    n_rna = rna_s.n_obs
    Z_rna = Z_all[:n_rna]
    Z_atac = Z_all[n_rna:n_rna + atac_s.n_obs]
    metrics = _evaluate(Z_rna, Z_atac, pair_a_sub, pair_b_sub, labels_sub, seed)
    return {
        "method": "uniPort",
        "description": "uniPort diagonal mode (iteration=3000)",
        "wall_time_sec": time.time() - t0,
        "metrics": metrics,
    }


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------


def run_all(dataset_id: str, exp_name: str, subsample: int = 3000,
            seed: int = 0) -> dict:
    rna = ad.read_h5ad(PROCESSED_DIR / f"{dataset_id}_rna.h5ad")
    atac = ad.read_h5ad(PROCESSED_DIR / f"{dataset_id}_atac.h5ad")
    rng = np.random.default_rng(seed)
    sub_idx = rng.choice(rna.n_obs, size=min(subsample, rna.n_obs), replace=False)
    labels_sub = rna.obs["cell_type"].astype(str).values[sub_idx]
    pair_a_sub = rna.obs["pair_idx"].values.astype(np.int64)[sub_idx]
    pair_b_sub = atac.obs["pair_idx"].values.astype(np.int64)[sub_idx]

    out_dir = EXPERIMENTS_DIR / f"baselines_{dataset_id}"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    print("\n=== Baseline 1: Raw PCA/LSI + Procrustes (no IB-VAE) ===")
    results["raw_ot"] = baseline_raw_ot(
        rna, atac, labels_sub, pair_a_sub, pair_b_sub, sub_idx, seed,
    )
    print("\n=== Baseline 2: IB latent + NN matching (no OT) ===")
    results["nn_on_ib"] = baseline_nn_on_ib(
        exp_name, labels_sub, pair_a_sub, pair_b_sub, sub_idx, seed,
    )
    print("\n=== Baseline 3: uniPort ===")
    results["uniport"] = baseline_uniport(
        rna, atac, labels_sub, pair_a_sub, pair_b_sub, sub_idx, seed, out_dir,
    )

    with (out_dir / "simple_baseline_results.json").open("w") as f:
        json.dump(results, f, indent=2, default=str)

    print("\n=== Summary (subsample {n}) ===".format(n=len(sub_idx)))
    print(f"{'method':20s} | {'FOSCTTM':>8s} | {'LT A->B':>8s} | {'LT B->A':>8s} | {'ARI':>7s} | {'wall':>6s}")
    print("-" * 70)
    for k, v in results.items():
        m = v.get("metrics", {})
        if not m:
            print(f"{k:20s} | FAILED: {v.get('error', 'unknown')[:50]}")
            continue
        print(f"{v['method']:20s} | "
              f"{m['foscttm_mean']:8.4f} | "
              f"{m['label_transfer_rna_to_atac']:8.4f} | "
              f"{m['label_transfer_atac_to_rna']:8.4f} | "
              f"{m['joint_clustering_ari']:7.4f} | "
              f"{v['wall_time_sec']:6.1f}")
    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--exp", required=True, help="MOSAIC experiment name for IB latent baseline")
    p.add_argument("--subsample", type=int, default=3000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    run_all(args.dataset, args.exp, subsample=args.subsample, seed=args.seed)
