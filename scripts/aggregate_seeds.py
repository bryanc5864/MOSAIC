#!/usr/bin/env python3
# MIT License - Bryan Cheng, 2026
# Part of MOSAIC - multi-seed aggregation for W5 fix
"""Aggregate multi-seed results into mean +/- std summaries.

Given a list of (seed, exp_name) tuples for a single dataset, loads each
results.json and computes mean +/- std for every primary metric. Writes
a single aggregate JSON and a markdown table to stdout.

Usage:
    python -m scripts.aggregate_seeds --dataset pbmc10k_multiome \\
        --runs exp001_pbmc_final:0 exp001_pbmc_seed1:1 exp001_pbmc_seed2:2
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from src.utils.paths import EXPERIMENTS_DIR


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--runs", nargs="+", required=True,
                   help="List of 'exp_name:seed' pairs")
    p.add_argument("--label", default=None,
                   help="Optional pretty label for the dataset (defaults to dataset id)")
    p.add_argument("--out-suffix", default="",
                   help="Suffix appended to the aggregate JSON filename, e.g. '_beta0001'")
    args = p.parse_args()

    label = args.label or args.dataset

    seeds = []
    metrics_per_seed = []
    cluster_metrics_per_seed = []
    for run_spec in args.runs:
        exp_name, seed_str = run_spec.split(":")
        seed = int(seed_str)
        seeds.append(seed)
        results_path = EXPERIMENTS_DIR / exp_name / "results.json"
        with results_path.open() as f:
            res = json.load(f)
        m = res["metrics"]
        metrics_per_seed.append({
            "seed": seed,
            "exp": exp_name,
            "foscttm_mean": m["foscttm"]["foscttm_mean"],
            "foscttm_a_to_b": m["foscttm"]["foscttm_a_to_b"],
            "foscttm_b_to_a": m["foscttm"]["foscttm_b_to_a"],
            "lt_rna_to_atac": m["label_transfer_rna_to_atac"],
            "lt_atac_to_rna": m["label_transfer_atac_to_rna"],
            "joint_ari": m["joint_clustering_ari"],
            "entropy_mean": res["alignment"]["entropy_mean"],
            "entropy_rho": m["entropy_error_corr"]["spearman_rho"],
            "wall_time_sec": res["wall_time_sec"],
        })
        # Also try to load cluster entropy analysis
        cea = EXPERIMENTS_DIR / exp_name / "cluster_entropy_analysis.json"
        if cea.exists():
            with cea.open() as f:
                cea_res = json.load(f)
            cluster_metrics_per_seed.append({
                "seed": seed,
                "argmax_cluster_accuracy": cea_res["cluster_level_entropy"]["argmax_cluster_accuracy"],
                "mean_h_cluster": cea_res["cluster_level_entropy"]["mean"],
                "auroc_wrong_cluster": cea_res["cluster_level_entropy"]["auroc_entropy_vs_wrong_cluster"],
                "n_wrong_cluster": cea_res["cluster_level_entropy"]["n_wrong_cluster_cells"],
            })

    def mean_std(key: str) -> tuple[float, float]:
        vals = np.array([m[key] for m in metrics_per_seed])
        return float(vals.mean()), float(vals.std(ddof=1) if len(vals) > 1 else 0.0)

    aggregated = {
        "dataset": args.dataset,
        "seeds": seeds,
        "n_seeds": len(seeds),
        "per_seed": metrics_per_seed,
    }
    for key in ["foscttm_mean", "foscttm_a_to_b", "foscttm_b_to_a",
                "lt_rna_to_atac", "lt_atac_to_rna", "joint_ari",
                "entropy_mean", "entropy_rho", "wall_time_sec"]:
        mu, sd = mean_std(key)
        aggregated[f"{key}_mean"] = mu
        aggregated[f"{key}_std"] = sd

    if cluster_metrics_per_seed:
        aggregated["cluster_per_seed"] = cluster_metrics_per_seed
        for key in ["argmax_cluster_accuracy", "mean_h_cluster", "auroc_wrong_cluster"]:
            vals = np.array([m[key] for m in cluster_metrics_per_seed])
            aggregated[f"{key}_mean"] = float(vals.mean())
            aggregated[f"{key}_std"] = float(vals.std(ddof=1) if len(vals) > 1 else 0.0)

    out_path = EXPERIMENTS_DIR / f"aggregate_{args.dataset}{args.out_suffix}.json"
    with out_path.open("w") as f:
        json.dump(aggregated, f, indent=2)

    print(f"\n=== {label} (n_seeds = {len(seeds)}) ===")
    print(f"{'metric':28s} | {'mean':>9s} +/- {'std':>7s}")
    print("-" * 60)
    for key, pretty in [
        ("foscttm_mean",     "FOSCTTM mean"),
        ("foscttm_a_to_b",   "  A->B"),
        ("foscttm_b_to_a",   "  B->A"),
        ("lt_rna_to_atac",   "LT RNA->ATAC"),
        ("lt_atac_to_rna",   "LT ATAC->RNA"),
        ("joint_ari",        "Joint ARI"),
        ("entropy_mean",     "Mean H_cell"),
        ("entropy_rho",      "Spearman rho"),
    ]:
        print(f"{pretty:28s} | {aggregated[f'{key}_mean']:9.4f} +/- {aggregated[f'{key}_std']:7.4f}")
    if cluster_metrics_per_seed:
        print(f"{'Argmax cluster acc':28s} | "
              f"{aggregated['argmax_cluster_accuracy_mean']:9.4f} +/- "
              f"{aggregated['argmax_cluster_accuracy_std']:7.4f}")
        print(f"{'Mean H_cluster':28s} | "
              f"{aggregated['mean_h_cluster_mean']:9.4f} +/- "
              f"{aggregated['mean_h_cluster_std']:7.4f}")
        print(f"{'AUROC wrong cluster':28s} | "
              f"{aggregated['auroc_wrong_cluster_mean']:9.4f} +/- "
              f"{aggregated['auroc_wrong_cluster_std']:7.4f}")
    print(f"\nsaved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
