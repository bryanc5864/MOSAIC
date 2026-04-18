#!/usr/bin/env python3
"""Run OT alignment + metrics for CITE-seq seed 1 (training already done)."""
import sys, json, time
sys.path.insert(0, 'C:/Users/Maozer/projects/MOSAIC')

import numpy as np
import anndata as ad

from src.evaluation.metrics import (
    entropy_error_corr, foscttm, joint_clustering_ari, label_transfer_accuracy,
)
from src.models.align_post import apply_alignment, fit_orthogonal_procrustes
from src.models.ot_align import sinkhorn_align
from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR

EXP = "exp001_citeseq_seed1"
DATASET = "citeseq_pbmc"
SEED = 1
EPSILON = 0.05
OT_SUBSAMPLE = 3000

exp_dir = EXPERIMENTS_DIR / EXP
t0 = time.time()

Z_rna = np.load(exp_dir / "z_rna.npy")
Z_atac = np.load(exp_dir / "z_atac.npy")

rna = ad.read_h5ad(PROCESSED_DIR / f"{DATASET}_rna.h5ad")
atac = ad.read_h5ad(PROCESSED_DIR / f"{DATASET}_atac.h5ad")
pair_a = rna.obs["pair_idx"].values.astype(np.int64)
pair_b = atac.obs["pair_idx"].values.astype(np.int64)
labels_a = rna.obs["cell_type"].astype(str).values
labels_b = atac.obs["cell_type"].astype(str).values

print("Procrustes alignment...", flush=True)
proc = fit_orthogonal_procrustes(Z_atac, Z_rna, labels_src=labels_b, labels_tgt=labels_a)
Z_atac_aligned = apply_alignment(Z_atac, proc)
print(f"  clusters={proc.n_clusters}, residual={proc.residual:.4f}")

n_cells = Z_rna.shape[0]
n_sub = min(OT_SUBSAMPLE, n_cells)
sub_rng = np.random.default_rng(SEED)
sub_idx = sub_rng.choice(n_cells, size=n_sub, replace=False)

print(f"Sinkhorn on {n_sub}/{n_cells}...", flush=True)
align = sinkhorn_align(Z_rna[sub_idx], Z_atac_aligned[sub_idx], epsilon=EPSILON)

np.save(exp_dir / "alignment_plan_subsample.npy", align.plan.astype(np.float32))
np.save(exp_dir / "alignment_entropy_subsample.npy", align.entropy.astype(np.float32))
np.save(exp_dir / "alignment_top_match_subsample.npy", align.top_match.astype(np.int64))
np.save(exp_dir / "ot_subsample_indices.npy", sub_idx.astype(np.int64))
np.save(exp_dir / "z_atac_aligned.npy", Z_atac_aligned.astype(np.float32))

print("Computing metrics...", flush=True)
foscttm_res = foscttm(Z_rna, Z_atac_aligned, pair_a, pair_b)
lt_rna_to_atac = label_transfer_accuracy(Z_rna, labels_a, Z_atac_aligned, labels_b, k=15)
lt_atac_to_rna = label_transfer_accuracy(Z_atac_aligned, labels_b, Z_rna, labels_a, k=15)
ari = joint_clustering_ari(Z_rna, Z_atac_aligned, labels_a, labels_b,
                            n_clusters=len(np.unique(labels_a)), seed=SEED)
pair_a_sub = pair_a[sub_idx]
pair_b_sub = pair_b[sub_idx]
ent_corr = entropy_error_corr(align.entropy, Z_rna[sub_idx], Z_atac_aligned[sub_idx],
                               pair_a_sub, pair_b_sub)

results = {
    "dataset_id": DATASET,
    "experiment": EXP,
    "wall_time_sec": time.time() - t0,
    "alignment": {
        "epsilon": EPSILON,
        "ot_subsample_n": n_sub,
        "ot_subsample_of": n_cells,
        "procrustes": {"n_clusters": proc.n_clusters, "centroid_residual": proc.residual},
    },
    "metrics": {
        "foscttm": {
            "foscttm_a_to_b": float(foscttm_res["foscttm_a_to_b"]),
            "foscttm_b_to_a": float(foscttm_res["foscttm_b_to_a"]),
            "foscttm_mean": float(foscttm_res["foscttm_mean"]),
            "n_paired": int(foscttm_res["n_paired"]),
        },
        "label_transfer_rna_to_atac": float(lt_rna_to_atac),
        "label_transfer_atac_to_rna": float(lt_atac_to_rna),
        "joint_clustering_ari": float(ari),
        "entropy_error_corr": {
            "spearman_rho": float(ent_corr["spearman_rho"]),
            "spearman_p": float(ent_corr["spearman_p"]),
            "n": int(ent_corr.get("n", len(sub_idx))),
        },
    },
}

out_path = exp_dir / "results.json"
with out_path.open("w") as f:
    json.dump(results, f, indent=2)

m = results["metrics"]
print(f"\nDone! FOSCTTM={m['foscttm']['foscttm_mean']:.4f}")
print(f"LT RNA->ADT: {m['label_transfer_rna_to_atac']:.4f}")
print(f"LT ADT->RNA: {m['label_transfer_atac_to_rna']:.4f}")
print(f"ARI: {m['joint_clustering_ari']:.4f}")
print(f"Saved: {out_path}")
