#!/usr/bin/env python3
# MIT License
# Part of MOSAIC
"""Which surface markers mark the cells the alignment is unsure about.

High cluster entropy on the RNA -> protein alignment means the transcriptome
doesn't cleanly predict the surface phenotype: transitional B cells, exhausted
T cells with odd markers, myeloid boundary cells.

Takes the top and bottom entropy deciles and compares each protein marker's
mean expression between them, then summarizes entropy per leiden cluster.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy import stats

from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR

EXP_ID = "exp001_citeseq"
DATASET = "citeseq_pbmc"

MARKER_ANNOTATIONS = {
    "CD3_TotalSeqB":   "T cell lineage (pan-T)",
    "CD4_TotalSeqB":   "Helper T / Treg",
    "CD8a_TotalSeqB":  "Cytotoxic T",
    "CD14_TotalSeqB":  "Classical monocyte",
    "CD15_TotalSeqB":  "Granulocyte",
    "CD16_TotalSeqB":  "NK / non-classical monocyte",
    "CD56_TotalSeqB":  "NK cell",
    "CD19_TotalSeqB":  "B cell",
    "CD25_TotalSeqB":  "Treg / activated T",
    "CD45RA_TotalSeqB":"Naive T / terminally differentiated",
    "CD45RO_TotalSeqB":"Memory T",
    "PD-1_TotalSeqB":  "T cell exhaustion checkpoint",
    "TIGIT_TotalSeqB": "T cell exhaustion checkpoint",
    "CD127_TotalSeqB": "Memory T / T cell survival",
}


def main() -> int:
    exp_dir = EXPERIMENTS_DIR / EXP_ID

    H_full_path = exp_dir / "alignment_entropy_cluster.npy"
    if not H_full_path.exists():
        print(f"{H_full_path} not found; run exp001_citeseq first")
        return 1

    H = np.load(H_full_path)
    print(f"Loaded H_cluster: shape={H.shape}, mean={H.mean():.4f}, std={H.std():.4f}")

    # protein sits in the atac slot for cite-seq
    prot = ad.read_h5ad(PROCESSED_DIR / f"{DATASET}_atac.h5ad")
    rna = ad.read_h5ad(PROCESSED_DIR / f"{DATASET}_rna.h5ad")

    # H is over the ot subsample, so we need its indices
    sub_idx_path = exp_dir / "ot_subsample_indices.npy"
    if sub_idx_path.exists():
        sub_idx = np.load(sub_idx_path)
        print(f"OT subsample: {len(sub_idx)} cells")
        H_sub = H
        prot_sub = prot[sub_idx]
        rna_sub = rna[sub_idx]
    else:
        # no subsample file, assume H covers everything
        H_sub = H
        prot_sub = prot
        rna_sub = rna

    n = len(H_sub)
    pct10 = int(0.10 * n)
    pct10 = max(pct10, 10)

    order = np.argsort(H_sub)
    low_idx = order[:pct10]
    high_idx = order[-pct10:]

    X_prot = prot_sub.X
    if sp.issparse(X_prot):
        X_prot = X_prot.toarray()

    marker_results = []
    for i, marker in enumerate(prot_sub.var_names):
        vals = X_prot[:, i]
        low_vals = vals[low_idx]
        high_vals = vals[high_idx]
        stat, pval = stats.mannwhitneyu(high_vals, low_vals, alternative="two-sided")
        fc = (high_vals.mean() - low_vals.mean())  # clr units, not a fold change
        annotation = MARKER_ANNOTATIONS.get(marker, "")
        marker_results.append({
            "marker": marker,
            "clinical_interpretation": annotation,
            "mean_high_entropy": float(high_vals.mean()),
            "mean_low_entropy": float(low_vals.mean()),
            "expression_diff_high_minus_low": float(fc),
            "mannwhitney_p": float(pval),
            "is_enriched_in_uncertain": bool(fc > 0),
        })

    marker_results.sort(key=lambda x: abs(x["expression_diff_high_minus_low"]), reverse=True)

    leiden = rna_sub.obs["leiden"].astype(str).values
    clusters = sorted(np.unique(leiden), key=lambda c: int(c))
    per_cluster = []
    for c in clusters:
        mask = leiden == c
        if mask.sum() < 5:
            continue
        h_c = H_sub[mask]
        per_cluster.append({
            "cluster": c,
            "n": int(mask.sum()),
            "mean_h": float(h_c.mean()),
            "std_h": float(h_c.std()),
            "high_entropy_frac": float((h_c > np.percentile(H_sub, 90)).mean()),
        })

    per_cluster.sort(key=lambda x: -x["mean_h"])

    output = {
        "experiment": "protein_uq_analysis",
        "dataset": DATASET,
        "source_exp": EXP_ID,
        "n_cells_analyzed": n,
        "top10pct_threshold": float(H_sub[order[-pct10]]),
        "bottom10pct_threshold": float(H_sub[order[pct10]]),
        "marker_analysis": marker_results,
        "per_cluster_entropy": per_cluster,
        "summary": {
            "markers_enriched_in_uncertain": [
                r["marker"] for r in marker_results if r["is_enriched_in_uncertain"]
                and r["mannwhitney_p"] < 0.05
            ],
            "markers_enriched_in_certain": [
                r["marker"] for r in marker_results if not r["is_enriched_in_uncertain"]
                and r["mannwhitney_p"] < 0.05
            ],
        }
    }

    out_dir = EXPERIMENTS_DIR / "protein_uq_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "results.json"
    with out_path.open("w") as f:
        json.dump(output, f, indent=2)

    print("\n[summary] protein marker uq")
    print(f"Cells analyzed: {n} (top/bottom 10% = {pct10} each)")
    print(f"\nTop markers by expression difference (uncertain vs. certain cells):")
    for r in marker_results[:7]:
        direction = "^ uncertain" if r["is_enriched_in_uncertain"] else "^ certain"
        direction_str = "^ uncertain" if r["is_enriched_in_uncertain"] else "^ certain"
        print(f"  {r['marker']:25s} {direction_str:15s} "
              f"diff={r['expression_diff_high_minus_low']:+.3f}  "
              f"p={r['mannwhitney_p']:.2e}  [{r['clinical_interpretation']}]")

    print(f"\nHighest-entropy clusters:")
    for r in per_cluster[:5]:
        print(f"  Cluster {r['cluster']:3s}  n={r['n']:4d}  mean_H={r['mean_h']:.4f}  "
              f"top10%_frac={r['high_entropy_frac']:.2f}")

    print(f"\nsaved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
