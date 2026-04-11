# MIT License - Bryan Cheng, 2026
# Part of MOSAIC
"""Fair-comparison baseline runner for MOSAIC benchmarks.

Runs SCOT (GW-OT on intra-modality geodesic distances) and uniPort (VAE + OT)
on the SAME preprocessed AnnData files that MOSAIC uses, and evaluates them
with the SAME metrics code (src/evaluation/metrics.py).

This enforces the fair-comparison protocol from RESEARCH_PLAN.md section 4.2:
    - Same data
    - Same preprocessing
    - Same splits
    - Same evaluation code
    - Same random seed

Usage:
    python -m src.baselines.run_baselines --dataset pbmc10k_multiome --methods scot uniport
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np

from src.evaluation.metrics import (
    entropy_error_corr,
    foscttm,
    joint_clustering_ari,
    label_transfer_accuracy,
)
from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR


def compute_metrics_on_aligned(Z_a: np.ndarray, Z_b: np.ndarray,
                                 pair_a: np.ndarray, pair_b: np.ndarray,
                                 labels_a: np.ndarray, labels_b: np.ndarray,
                                 seed: int = 0,
                                 entropy: np.ndarray | None = None) -> dict:
    """Compute the standard MOSAIC evaluation metrics on a pair of aligned
    embeddings. Used by both SCOT and uniPort baselines so results are
    directly comparable to run_experiment.py.
    """
    foscttm_res = foscttm(Z_a, Z_b, pair_a, pair_b)
    lt_a_to_b = label_transfer_accuracy(Z_a, labels_a, Z_b, labels_b, k=15)
    lt_b_to_a = label_transfer_accuracy(Z_b, labels_b, Z_a, labels_a, k=15)
    ari = joint_clustering_ari(Z_a, Z_b, labels_a, labels_b,
                                n_clusters=len(np.unique(labels_a)), seed=seed)
    metrics = {
        "foscttm": foscttm_res,
        "label_transfer_rna_to_atac": lt_a_to_b,
        "label_transfer_atac_to_rna": lt_b_to_a,
        "joint_clustering_ari": ari,
    }
    if entropy is not None:
        metrics["entropy_error_corr"] = entropy_error_corr(entropy, Z_a, Z_b, pair_a, pair_b)
    return metrics


# ----------------------------------------------------------------------------
# SCOT baseline
# ----------------------------------------------------------------------------


def run_scot_baseline(rna: ad.AnnData, atac: ad.AnnData, out_dir: Path,
                       seed: int = 0, subsample: int = 3000) -> dict:
    """Run SCOT v1 on the processed AnnDatas and save results.

    SCOT solves entropic Gromov-Wasserstein OT between the two modalities'
    intra-modality distance matrices. The inner iterations allocate multiple
    N x N matrices and, at float64, 11k x 11k GW reliably OOM-kills Python
    at 16 GB RAM. We subsample to `subsample` cells per modality with the
    same seeded index as the rest of the benchmark.
    """
    from src.baselines.scot_gw import run_scot

    # SCOT uses PCA of RNA and LSI of ATAC as input features.
    Z_rna = np.asarray(rna.obsm["X_pca"], dtype=np.float32)
    Z_atac = np.asarray(atac.obsm["X_lsi"], dtype=np.float32)
    print(f"[scot] full shapes: RNA {Z_rna.shape}, ATAC {Z_atac.shape}")

    rng = np.random.default_rng(seed)
    n_total = Z_rna.shape[0]
    n_sub = min(subsample, n_total)
    idx = rng.choice(n_total, size=n_sub, replace=False)
    Z_rna_sub = Z_rna[idx]
    Z_atac_sub = Z_atac[idx]
    print(f"[scot] subsampled to {n_sub}/{n_total}")

    t0 = time.time()
    result = run_scot(Z_rna_sub, Z_atac_sub, k=10, epsilon=5e-3, n_iter=300, verbose=True)
    wall_time = time.time() - t0

    # Barycentric projection of RNA into ATAC's LSI space (used for metrics).
    Z_rna_aligned = result.barycentric_embedding
    Z_atac_target = Z_atac_sub

    pair_a = rna.obs["pair_idx"].values.astype(np.int64)[idx]
    pair_b = atac.obs["pair_idx"].values.astype(np.int64)[idx]
    labels_a = rna.obs["cell_type"].astype(str).values[idx]
    labels_b = atac.obs["cell_type"].astype(str).values[idx]
    metrics = compute_metrics_on_aligned(
        Z_rna_aligned, Z_atac_target, pair_a, pair_b, labels_a, labels_b,
        seed=seed,
    )
    metrics["_evaluation_n"] = int(n_sub)
    # SCOT doesn't produce per-cell UQ, so no entropy calibration.
    np.save(out_dir / "scot_plan.npy", result.plan.astype(np.float32))
    np.save(out_dir / "scot_z_rna_aligned.npy", Z_rna_aligned)

    return {
        "method": "SCOT",
        "wall_time_sec": wall_time,
        "metrics": metrics,
    }


# ----------------------------------------------------------------------------
# uniPort baseline
# ----------------------------------------------------------------------------


def run_uniport_baseline(rna: ad.AnnData, atac: ad.AnnData, out_dir: Path,
                          seed: int = 0) -> dict:
    """Run uniPort (diagonal integration mode) on the processed AnnDatas.

    uniPort's 'd' mode is the unpaired integration setting. It takes two
    AnnDatas with DIFFERENT feature spaces and produces a shared latent.
    """
    import uniport as up

    # uniPort expects raw-ish features with a shared low-dim representation
    # via the 'use_rep' argument. We use precomputed PCA / LSI.
    rna_copy = rna.copy()
    atac_copy = atac.copy()
    # uniPort requires domain_id metadata
    rna_copy.obs["domain_id"] = 0
    atac_copy.obs["domain_id"] = 1
    rna_copy.obs["source"] = "rna"
    atac_copy.obs["source"] = "atac"

    uniport_outdir = out_dir / "uniport_output"
    uniport_outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    try:
        adata_joint = up.Run(
            adatas=[rna_copy, atac_copy],
            mode="d",
            use_rep=["X_pca", "X_lsi"],
            iteration=3000,           # much smaller than the default 30000 to keep runtime in budget
            gpu=0,
            outdir=str(uniport_outdir),
            verbose=False,
            seed=seed,
        )
    except Exception as e:
        wall_time = time.time() - t0
        print(f"[uniport] ERROR: {e}")
        return {
            "method": "uniPort",
            "wall_time_sec": wall_time,
            "metrics": {},
            "error": str(e),
        }
    wall_time = time.time() - t0

    # uniPort produces a joint AnnData with the aligned latent in .obsm['latent']
    if "latent" in adata_joint.obsm:
        Z = adata_joint.obsm["latent"]
    elif "X_latent" in adata_joint.obsm:
        Z = adata_joint.obsm["X_latent"]
    else:
        # Fall back: use .X
        Z = np.asarray(adata_joint.X, dtype=np.float32)
    n_rna = rna.n_obs
    Z_rna = Z[:n_rna]
    Z_atac = Z[n_rna:]

    pair_a = rna.obs["pair_idx"].values.astype(np.int64)
    pair_b = atac.obs["pair_idx"].values.astype(np.int64)
    labels_a = rna.obs["cell_type"].astype(str).values
    labels_b = atac.obs["cell_type"].astype(str).values
    metrics = compute_metrics_on_aligned(
        Z_rna, Z_atac, pair_a, pair_b, labels_a, labels_b, seed=seed,
    )
    np.save(out_dir / "uniport_z_rna.npy", Z_rna.astype(np.float32))
    np.save(out_dir / "uniport_z_atac.npy", Z_atac.astype(np.float32))

    return {
        "method": "uniPort",
        "wall_time_sec": wall_time,
        "metrics": metrics,
    }


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--methods", nargs="+", choices=["scot", "uniport"],
                        default=["scot", "uniport"])
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    rna_path = PROCESSED_DIR / f"{args.dataset}_rna.h5ad"
    atac_path = PROCESSED_DIR / f"{args.dataset}_atac.h5ad"
    rna = ad.read_h5ad(rna_path)
    atac = ad.read_h5ad(atac_path)

    out_dir = EXPERIMENTS_DIR / f"baselines_{args.dataset}"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    for method in args.methods:
        print(f"\n=== Running {method} on {args.dataset} ===")
        if method == "scot":
            results["scot"] = run_scot_baseline(rna, atac, out_dir, seed=args.seed)
        elif method == "uniport":
            results["uniport"] = run_uniport_baseline(rna, atac, out_dir, seed=args.seed)

    with (out_dir / "baseline_results.json").open("w") as f:
        json.dump(results, f, indent=2, default=str)

    # Print summary
    print("\n=== Baseline results summary ===")
    for name, res in results.items():
        m = res.get("metrics", {})
        if not m:
            print(f"  {name:10s}  FAILED: {res.get('error', 'unknown')}")
            continue
        f_ = m.get("foscttm", {}).get("foscttm_mean", float("nan"))
        lt_ab = m.get("label_transfer_rna_to_atac", float("nan"))
        lt_ba = m.get("label_transfer_atac_to_rna", float("nan"))
        ari = m.get("joint_clustering_ari", float("nan"))
        wall = res.get("wall_time_sec", 0.0)
        print(f"  {name:10s}  FOSCTTM {f_:.4f}  LT {lt_ab:.3f}/{lt_ba:.3f}  "
              f"ARI {ari:.4f}  wall {wall:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
