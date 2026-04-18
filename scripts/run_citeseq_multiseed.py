#!/usr/bin/env python3
# MIT License - Bryan Cheng, 2026
# Part of MOSAIC
"""Run CITE-seq experiment at multiple seeds for 3-seed statistical reporting."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.training.run_experiment import run_experiment
from src.utils.paths import EXPERIMENTS_DIR

DATASET = "citeseq_pbmc"
SEEDS = [1, 2]
BETA = 0.001
LAMBDA_PRED = 5.0


def main() -> int:
    for seed in SEEDS:
        exp_name = f"exp001_citeseq_seed{seed}"
        out_dir = EXPERIMENTS_DIR / exp_name

        if (out_dir / "results.json").exists():
            print(f"Seed {seed} already done, skipping.")
            continue

        print(f"\n{'='*60}")
        print(f"Running CITE-seq, seed={seed}")
        print(f"{'='*60}")

        result = run_experiment(
            dataset_id=DATASET,
            exp_name=exp_name,
            epochs=200,
            batch_size=512,
            beta=BETA,
            lambda_pred=LAMBDA_PRED,
            patience=999,
            seed=seed,
            device="cuda",
            epsilon=0.05,
            ot_subsample=3000,
        )
        print(f"Seed {seed} done. FOSCTTM={result['metrics']['foscttm']['foscttm_mean']:.4f}")

    # Aggregate 3 seeds (0, 1, 2)
    print("\n--- Aggregating 3 seeds ---")
    exp_names = ["exp001_citeseq", "exp001_citeseq_seed1", "exp001_citeseq_seed2"]
    results = []
    for name in exp_names:
        p = EXPERIMENTS_DIR / name / "results.json"
        if p.exists():
            results.append(json.loads(p.read_text()))
        else:
            print(f"WARNING: {p} not found, skipping from aggregate")

    if not results:
        print("No results to aggregate.")
        return 1

    def collect(key_path: list[str]) -> list[float]:
        out = []
        for r in results:
            val = r
            for k in key_path:
                val = val[k]
            out.append(float(val))
        return out

    foscttm = collect(["metrics", "foscttm", "foscttm_mean"])
    lt_rna = collect(["metrics", "label_transfer_rna_to_atac"])
    lt_atac = collect(["metrics", "label_transfer_atac_to_rna"])
    ari = collect(["metrics", "joint_clustering_ari"])

    aggregate = {
        "dataset": DATASET,
        "n_seeds": len(results),
        "seeds": [0, 1, 2][:len(results)],
        "foscttm_mean": float(np.mean(foscttm)),
        "foscttm_std": float(np.std(foscttm)),
        "lt_rna_to_atac_mean": float(np.mean(lt_rna)),
        "lt_rna_to_atac_std": float(np.std(lt_rna)),
        "lt_atac_to_rna_mean": float(np.mean(lt_atac)),
        "lt_atac_to_rna_std": float(np.std(lt_atac)),
        "joint_ari_mean": float(np.mean(ari)),
        "joint_ari_std": float(np.std(ari)),
        "per_seed": [
            {
                "seed": i, "exp": exp_names[i],
                "foscttm": foscttm[i], "lt_rna": lt_rna[i],
                "lt_atac": lt_atac[i], "ari": ari[i],
            }
            for i in range(len(results))
        ],
    }

    out_path = EXPERIMENTS_DIR / "aggregate_citeseq_3seed.json"
    with out_path.open("w") as f:
        json.dump(aggregate, f, indent=2)

    print(f"\n=== CITE-seq 3-seed Aggregate ===")
    print(f"FOSCTTM:  {aggregate['foscttm_mean']:.4f} ± {aggregate['foscttm_std']:.4f}")
    print(f"LT RNA→P: {aggregate['lt_rna_to_atac_mean']:.4f} ± {aggregate['lt_rna_to_atac_std']:.4f}")
    print(f"LT P→RNA: {aggregate['lt_atac_to_rna_mean']:.4f} ± {aggregate['lt_atac_to_rna_std']:.4f}")
    print(f"ARI:      {aggregate['joint_ari_mean']:.4f} ± {aggregate['joint_ari_std']:.4f}")
    print(f"Saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
