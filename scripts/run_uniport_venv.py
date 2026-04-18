#!/usr/bin/env python3
# MIT License
# Part of MOSAIC - uniPort baseline in isolated venv
"""Run uniPort as a baseline in a venv with numpy 1.x.

uniPort 1.3 uses `np.Inf` which was removed in numpy 2.0. This script is
designed to be run from `venv_uniport/Scripts/python.exe` rather than the
main Python, so it has its own numpy 1.x and scanpy installation.

Usage:
    venv_uniport/Scripts/python.exe scripts/run_uniport_venv.py \\
        --dataset pbmc10k_multiome --subsample 3000 --seed 0

Output written to experiments/baselines_<dataset>/uniport_venv_results.json
and also saves per-metric values so the main environment's aggregator can
pick them up.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np

# Path fix for imports from the main repo's src/
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


def _sqdist(A, B):
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    a_sq = (A ** 2).sum(1, keepdims=True)
    b_sq = (B ** 2).sum(1, keepdims=True).T
    D = a_sq + b_sq - 2.0 * A @ B.T
    np.maximum(D, 0.0, out=D)
    return D


def foscttm(Z_a, Z_b, pair_a, pair_b):
    b_lookup = {int(k): i for i, k in enumerate(pair_b)}
    paired_a, paired_b = [], []
    for i, k in enumerate(pair_a):
        j = b_lookup.get(int(k))
        if j is not None:
            paired_a.append(i); paired_b.append(j)
    paired_a = np.asarray(paired_a); paired_b = np.asarray(paired_b)
    n = len(paired_a)
    if n == 0: return {"mean": float("nan")}
    Z_a_p = Z_a[paired_a]; Z_b_p = Z_b
    D = _sqdist(Z_a_p, Z_b_p)
    true_d = D[np.arange(n), paired_b]
    frac_a = (D < true_d[:, None]).sum(1) / max(Z_b.shape[0] - 1, 1)
    Z_b_q = Z_b[paired_b]
    D2 = _sqdist(Z_b_q, Z_a)
    true_d_ba = D2[np.arange(n), paired_a]
    frac_b = (D2 < true_d_ba[:, None]).sum(1) / max(Z_a.shape[0] - 1, 1)
    return {
        "foscttm_mean": float(0.5 * (frac_a.mean() + frac_b.mean())),
        "foscttm_a_to_b": float(frac_a.mean()),
        "foscttm_b_to_a": float(frac_b.mean()),
    }


def label_transfer(Z_src, labels_src, Z_tgt, labels_tgt, k=15):
    clf = KNeighborsClassifier(n_neighbors=k, metric="euclidean", n_jobs=-1)
    clf.fit(Z_src, labels_src)
    pred = clf.predict(Z_tgt)
    return float((pred == labels_tgt).mean())


def joint_ari(Z_a, Z_b, labels_a, labels_b, k, seed=0):
    Z = np.concatenate([Z_a, Z_b], axis=0)
    lbl = np.concatenate([labels_a, labels_b])
    pred = KMeans(n_clusters=k, random_state=seed, n_init=10).fit_predict(Z)
    return float(adjusted_rand_score(lbl, pred))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--subsample", type=int, default=3000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    import uniport as up

    processed = ROOT / "data" / "processed"
    # Prefer the _uniport variant (uns-stripped for venv compatibility) if present
    rna_path = processed / f"{args.dataset}_rna_uniport.h5ad"
    atac_path = processed / f"{args.dataset}_atac_uniport.h5ad"
    if not rna_path.exists():
        rna_path = processed / f"{args.dataset}_rna.h5ad"
        atac_path = processed / f"{args.dataset}_atac.h5ad"
    rna = ad.read_h5ad(rna_path)
    atac = ad.read_h5ad(atac_path)

    rng = np.random.default_rng(args.seed)
    n = min(args.subsample, rna.n_obs)
    idx = rng.choice(rna.n_obs, size=n, replace=False)

    rna_s = rna[idx].copy()
    atac_s = atac[idx].copy()
    rna_s.obs["domain_id"] = 0
    atac_s.obs["domain_id"] = 1
    rna_s.obs["source"] = "rna"
    atac_s.obs["source"] = "atac"

    out_dir = ROOT / "experiments" / f"baselines_{args.dataset}"
    out_dir.mkdir(parents=True, exist_ok=True)
    uniport_outdir = out_dir / "uniport_output"
    uniport_outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    try:
        adata_joint = up.Run(
            adatas=[rna_s, atac_s],
            mode="d",
            use_rep=["X_pca", "X_lsi"],
            iteration=3000,
            gpu=-1,  # CPU only
            outdir=str(uniport_outdir),
            verbose=False,
            seed=args.seed,
        )
    except Exception as e:
        out = {"method": "uniPort", "error": str(e), "wall_time_sec": time.time() - t0}
        with (out_dir / "uniport_venv_results.json").open("w") as f:
            json.dump(out, f, indent=2)
        print(f"uniPort FAILED: {e}")
        return 1
    wall = time.time() - t0

    # Extract latent
    Z_all = None
    for key in ("latent", "X_latent", "X_uniport"):
        if key in adata_joint.obsm:
            Z_all = np.asarray(adata_joint.obsm[key], dtype=np.float32)
            break
    if Z_all is None:
        X = adata_joint.X
        Z_all = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=np.float32)

    n_rna = rna_s.n_obs
    Z_rna = Z_all[:n_rna]
    Z_atac = Z_all[n_rna:n_rna + atac_s.n_obs]

    pair_a = rna_s.obs["pair_idx"].values.astype(np.int64)
    pair_b = atac_s.obs["pair_idx"].values.astype(np.int64)
    labels_a = rna_s.obs["cell_type"].astype(str).values
    labels_b = atac_s.obs["cell_type"].astype(str).values

    f = foscttm(Z_rna, Z_atac, pair_a, pair_b)
    lt_ab = label_transfer(Z_rna, labels_a, Z_atac, labels_b, k=15)
    lt_ba = label_transfer(Z_atac, labels_b, Z_rna, labels_a, k=15)
    ari = joint_ari(Z_rna, Z_atac, labels_a, labels_b, k=len(np.unique(labels_a)), seed=args.seed)

    out = {
        "method": "uniPort",
        "dataset": args.dataset,
        "subsample_n": n,
        "seed": args.seed,
        "wall_time_sec": wall,
        "metrics": {
            "foscttm": f,
            "label_transfer_rna_to_atac": lt_ab,
            "label_transfer_atac_to_rna": lt_ba,
            "joint_clustering_ari": ari,
        },
    }
    with (out_dir / "uniport_venv_results.json").open("w") as f_out:
        json.dump(out, f_out, indent=2)

    print(f"\n=== uniPort on {args.dataset} ({n} cells, seed {args.seed}) ===")
    print(f"  FOSCTTM     : {f['foscttm_mean']:.4f}")
    print(f"  LT RNA->ATAC: {lt_ab:.4f}")
    print(f"  LT ATAC->RNA: {lt_ba:.4f}")
    print(f"  Joint ARI   : {ari:.4f}")
    print(f"  Wall time   : {wall:.1f} sec")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
