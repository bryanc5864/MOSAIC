#!/usr/bin/env python3
# MIT License
# Part of MOSAIC - Exp 4: ablation sweeps
"""Run a set of MOSAIC variants to quantify each design choice's contribution.

Sweep:
  (a) beta in {0.001, 0.01, 0.1}                         -- IB bottleneck strength
  (b) lambda_pred in {1.0, 5.0, 20.0}                    -- cross-modal prediction weight
  (c) no cross-modal head (lambda_pred=0)                -- ablates cross-modal pred entirely
  (d) epsilon in {0.02, 0.05, 0.2}                       -- OT sharpness

Each variant is a single call to run_experiment with different kwargs.
Runs are small (200 epochs each, ~2 min per variant) so the full sweep
takes ~15-20 minutes.

Results written to experiments/ablation_<dataset>/ and summarized at end.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from src.training.run_experiment import run_experiment


SWEEP = [
    # Base config (should match exp001_pbmc_final)
    {"name": "base", "beta": 0.01, "lambda_pred": 5.0, "epsilon": 0.05},
    # beta sweep
    {"name": "beta_0.001", "beta": 0.001, "lambda_pred": 5.0, "epsilon": 0.05},
    {"name": "beta_0.1", "beta": 0.1, "lambda_pred": 5.0, "epsilon": 0.05},
    # lambda_pred sweep
    {"name": "lambda_1", "beta": 0.01, "lambda_pred": 1.0, "epsilon": 0.05},
    {"name": "lambda_20", "beta": 0.01, "lambda_pred": 20.0, "epsilon": 0.05},
    # no cross-modal head
    {"name": "no_cross_head", "beta": 0.01, "lambda_pred": 0.0, "epsilon": 0.05},
    # epsilon sweep
    {"name": "eps_0.02", "beta": 0.01, "lambda_pred": 5.0, "epsilon": 0.02},
    {"name": "eps_0.2", "beta": 0.01, "lambda_pred": 5.0, "epsilon": 0.2},
]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", default="pbmc10k_multiome")
    p.add_argument("--epochs", type=int, default=200)
    p.add_argument("--subsample", type=int, default=4000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--only", nargs="*", default=None,
                   help="run only these variant names")
    args = p.parse_args()

    variants = SWEEP if args.only is None else \
        [v for v in SWEEP if v["name"] in args.only]

    all_results = {}
    for v in variants:
        exp_name = f"ablation_{args.dataset}_{v['name']}"
        print(f"\n{'='*60}")
        print(f"  variant: {v['name']}")
        print(f"  beta={v['beta']}, lambda_pred={v['lambda_pred']}, epsilon={v['epsilon']}")
        print(f"{'='*60}", flush=True)
        try:
            res = run_experiment(
                dataset_id=args.dataset,
                exp_name=exp_name,
                epochs=args.epochs,
                batch_size=512,
                beta=v["beta"],
                lambda_pred=v["lambda_pred"],
                patience=999,
                epsilon=v["epsilon"],
                seed=args.seed,
                device="cuda",
                procrustes=True,
                ot_subsample=args.subsample,
            )
        except Exception as e:
            print(f"VARIANT FAILED: {e}", flush=True)
            res = {"error": str(e)}
        # Keep only the key metrics for the summary
        m = res.get("metrics", {})
        all_results[v["name"]] = {
            "config": {"beta": v["beta"], "lambda_pred": v["lambda_pred"],
                       "epsilon": v["epsilon"]},
            "foscttm": m.get("foscttm", {}).get("foscttm_mean") if m else None,
            "lt_a_to_b": m.get("label_transfer_rna_to_atac") if m else None,
            "lt_b_to_a": m.get("label_transfer_atac_to_rna") if m else None,
            "ari": m.get("joint_clustering_ari") if m else None,
            "entropy_mean": res.get("alignment", {}).get("entropy_mean") if m else None,
            "entropy_rho": m.get("entropy_error_corr", {}).get("spearman_rho") if m else None,
            "wall_time_sec": res.get("wall_time_sec"),
            "error": res.get("error"),
        }

    out = Path("experiments") / f"ablation_{args.dataset}_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        json.dump(all_results, f, indent=2, default=str)

    print("\n" + "="*70)
    print(f"{'variant':20s} | FOSCTTM | LT A->B | LT B->A |   ARI  |  H mean | H rho ")
    print("-" * 80)
    for k, v in all_results.items():
        if v["foscttm"] is None:
            print(f"{k:20s} | FAILED: {v.get('error', 'unknown')[:40]}")
            continue
        print(f"{k:20s} | "
              f"{v['foscttm']:7.4f} | "
              f"{v['lt_a_to_b']:7.4f} | "
              f"{v['lt_b_to_a']:7.4f} | "
              f"{v['ari']:6.4f} | "
              f"{v['entropy_mean']:7.4f} | "
              f"{v['entropy_rho']:+7.4f}")
    print(f"\nsaved {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
