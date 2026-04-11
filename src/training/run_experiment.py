# MIT License - Bryan Cheng, 2026
# Part of MOSAIC
"""Top-level experiment runner for MOSAIC.

An "experiment" here is:
  1. Load the processed RNA and ATAC AnnDatas for a dataset.
  2. Train IBVAE_RNA (via train_ibvae).
  3. Train IBVAE_ATAC (via train_ibvae).
  4. Load both final embeddings (z_rna, z_atac).
  5. Run entropic OT alignment on the two embeddings.
  6. Compute metrics:
       - FOSCTTM
       - label transfer accuracy (RNA -> ATAC via top-k in latent)
       - joint clustering ARI
       - entropy-error Spearman correlation
  7. Save results.json and append a row to TRAINING_LOG.md / RESULTS.md.

This is the workhorse of Exp 1 (and forms the basis for Exp 2-4).
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict
from pathlib import Path

import anndata as ad
import numpy as np

from src.evaluation.metrics import (
    entropy_error_corr,
    foscttm,
    joint_clustering_ari,
    label_transfer_accuracy,
)
from src.models.align_post import apply_alignment, fit_orthogonal_procrustes
from src.models.ot_align import sinkhorn_align
from src.training.train_ibvae import TrainConfig, train
from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR


def run_experiment(dataset_id: str, exp_name: str, *,
                    epochs: int = 200, batch_size: int = 512, beta: float = 0.01,
                    lambda_pred: float = 1.0, patience: int = 15,
                    seed: int = 0, device: str = "cuda",
                    epsilon: float = 0.05,
                    procrustes: bool = True,
                    ot_subsample: int = 5000) -> dict:
    exp_dir = EXPERIMENTS_DIR / exp_name
    exp_dir.mkdir(parents=True, exist_ok=True)

    rna_path = PROCESSED_DIR / f"{dataset_id}_rna.h5ad"
    atac_path = PROCESSED_DIR / f"{dataset_id}_atac.h5ad"
    assert rna_path.exists() and atac_path.exists(), \
        f"processed data not found for {dataset_id}"

    t0 = time.time()

    # ---- 1. Train RNA IB-VAE ----
    print("\n--- Training RNA IB-VAE ---")
    rna_cfg = TrainConfig(
        modality="rna",
        processed_path=str(rna_path),
        out_dir=str(exp_dir),
        n_epochs=epochs, batch_size=batch_size, beta=beta,
        lambda_pred=lambda_pred, patience=patience,
        seed=seed, device=device,
    )
    rna_result = train(rna_cfg)

    # ---- 2. Train ATAC IB-VAE ----
    print("\n--- Training ATAC IB-VAE ---")
    atac_cfg = TrainConfig(
        modality="atac",
        processed_path=str(atac_path),
        out_dir=str(exp_dir),
        n_epochs=epochs, batch_size=batch_size, beta=beta,
        lambda_pred=lambda_pred, patience=patience,
        seed=seed, device=device,
    )
    atac_result = train(atac_cfg)

    # ---- 3. Load embeddings and metadata ----
    Z_rna = np.load(rna_result["embeddings"])
    Z_atac = np.load(atac_result["embeddings"])
    rna = ad.read_h5ad(rna_path)
    atac = ad.read_h5ad(atac_path)
    pair_a = rna.obs["pair_idx"].values.astype(np.int64)
    pair_b = atac.obs["pair_idx"].values.astype(np.int64)
    labels_a = rna.obs["cell_type"].astype(str).values
    labels_b = atac.obs["cell_type"].astype(str).values

    # ---- 4a. Optional Procrustes alignment of ATAC latent onto RNA frame ----
    procrustes_info: dict = {}
    if procrustes:
        print("\n--- Procrustes alignment (ATAC -> RNA frame) ---")
        proc = fit_orthogonal_procrustes(
            Z_src=Z_atac, Z_tgt=Z_rna,
            labels_src=labels_b, labels_tgt=labels_a,
        )
        Z_atac_aligned = apply_alignment(Z_atac, proc)
        procrustes_info = {
            "n_clusters": proc.n_clusters,
            "centroid_residual": proc.residual,
        }
        print(f"  fit on {proc.n_clusters} cluster centroids, "
              f"residual {proc.residual:.4f}")
    else:
        Z_atac_aligned = Z_atac
        procrustes_info = {"applied": False}

    # ---- 4b. Sinkhorn alignment on a subsample ----
    # Memory note: for N ~ 11000 the full NxN cost matrix is ~1 GB at float64
    # and POT's Sinkhorn iterations allocate several times that in intermediate
    # log/exp buffers. On a 16 GB workstation (with ~5 GB already held by the
    # training process) this reliably OOM-kills the interpreter silently. We
    # therefore run OT on a reproducible subsample of `ot_subsample` cells
    # per modality. FOSCTTM / label transfer / ARI are computed on the FULL
    # embeddings since they don't need the OT plan. The per-cell entropy
    # calibration metric (Exp 2) is computed on the subsample.
    n_cells_total = Z_rna.shape[0]
    n_sub = min(ot_subsample, n_cells_total)
    sub_rng = np.random.default_rng(seed)
    sub_idx = sub_rng.choice(n_cells_total, size=n_sub, replace=False)
    print(f"\n--- Running Sinkhorn alignment on {n_sub}/{n_cells_total} subsample ---")
    align = sinkhorn_align(
        Z_rna[sub_idx], Z_atac_aligned[sub_idx], epsilon=epsilon,
    )
    np.save(exp_dir / "alignment_plan_subsample.npy", align.plan.astype(np.float32))
    np.save(exp_dir / "alignment_entropy_subsample.npy", align.entropy.astype(np.float32))
    np.save(exp_dir / "alignment_top_match_subsample.npy", align.top_match.astype(np.int64))
    np.save(exp_dir / "ot_subsample_indices.npy", sub_idx.astype(np.int64))
    np.save(exp_dir / "z_atac_aligned.npy", Z_atac_aligned.astype(np.float32))

    # ---- 5. Compute metrics on the post-alignment embeddings ----
    # FOSCTTM / label transfer / ARI use the FULL aligned embeddings.
    # Entropy calibration uses only the OT subsample (entropy is only
    # defined for the cells on which we ran Sinkhorn).
    print("\n--- Computing metrics ---")
    foscttm_res = foscttm(Z_rna, Z_atac_aligned, pair_a, pair_b)
    lt_rna_to_atac = label_transfer_accuracy(Z_rna, labels_a, Z_atac_aligned, labels_b, k=15)
    lt_atac_to_rna = label_transfer_accuracy(Z_atac_aligned, labels_b, Z_rna, labels_a, k=15)
    ari = joint_clustering_ari(Z_rna, Z_atac_aligned, labels_a, labels_b,
                                n_clusters=len(np.unique(labels_a)), seed=seed)

    # For entropy calibration, restrict the full-set evaluation to the
    # subsampled cells. pair_a/pair_b for the RNA subsample point into the
    # ATAC subsample's pair_idx (the subsample is the same indices in both).
    pair_a_sub = pair_a[sub_idx]
    pair_b_sub = pair_b[sub_idx]
    ent_corr = entropy_error_corr(
        align.entropy,
        Z_rna[sub_idx], Z_atac_aligned[sub_idx],
        pair_a_sub, pair_b_sub,
    )

    results = {
        "dataset_id": dataset_id,
        "experiment": exp_name,
        "wall_time_sec": time.time() - t0,
        "rna_train": rna_result,
        "atac_train": atac_result,
        "alignment": {
            "epsilon": align.epsilon,
            "cost_scale": align.cost_scale,
            "entropy_mean": float(align.entropy.mean()),
            "entropy_std": float(align.entropy.std()),
            "top_match_prob_mean": float(align.top_match_prob.mean()),
            "procrustes": procrustes_info,
            "ot_subsample_n": int(n_sub),
            "ot_subsample_of": int(n_cells_total),
        },
        "metrics": {
            "foscttm": foscttm_res,
            "label_transfer_rna_to_atac": lt_rna_to_atac,
            "label_transfer_atac_to_rna": lt_atac_to_rna,
            "joint_clustering_ari": ari,
            "entropy_error_corr": ent_corr,
        },
    }

    with (exp_dir / "results.json").open("w") as f:
        json.dump(results, f, indent=2, default=str)

    print("\n=== Results ===")
    print(f"  FOSCTTM (mean)              : {foscttm_res['foscttm_mean']:.4f}")
    print(f"  Label transfer RNA->ATAC    : {lt_rna_to_atac:.4f}")
    print(f"  Label transfer ATAC->RNA    : {lt_atac_to_rna:.4f}")
    print(f"  Joint clustering ARI        : {ari:.4f}")
    print(f"  Entropy-error Spearman rho  : {ent_corr['spearman_rho']:.4f}")
    print(f"  Mean alignment entropy      : {align.entropy.mean():.4f}")
    print(f"  Wall time                   : {time.time() - t0:.1f} sec")
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Run a MOSAIC experiment end to end")
    p.add_argument("--dataset", required=True)
    p.add_argument("--name", required=True, help="experiment directory name")
    p.add_argument("--epochs", type=int, default=200)
    p.add_argument("--batch-size", type=int, default=512)
    p.add_argument("--beta", type=float, default=0.01)
    p.add_argument("--lambda-pred", type=float, default=1.0,
                   dest="lambda_pred")
    p.add_argument("--patience", type=int, default=15)
    p.add_argument("--epsilon", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--device", default="cuda")
    p.add_argument("--no-procrustes", dest="procrustes", action="store_false")
    p.add_argument("--ot-subsample", type=int, default=5000,
                   dest="ot_subsample")
    p.set_defaults(procrustes=True)
    args = p.parse_args()

    run_experiment(
        dataset_id=args.dataset,
        exp_name=args.name,
        epochs=args.epochs,
        batch_size=args.batch_size,
        beta=args.beta,
        lambda_pred=args.lambda_pred,
        patience=args.patience,
        seed=args.seed,
        device=args.device,
        epsilon=args.epsilon,
        procrustes=args.procrustes,
        ot_subsample=args.ot_subsample,
    )
