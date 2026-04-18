#!/usr/bin/env python3
# MIT License - Bryan Cheng, 2026
# Part of MOSAIC
"""Immune checkpoint biomarker analysis (CITE-seq).

PD-1 and TIGIT are immune checkpoint receptors expressed on exhausted T cells.
Their expression predicts response to checkpoint inhibitor immunotherapy
(anti-PD-1/anti-TIGIT antibodies). MOSAIC's cluster-resolved entropy for
these cells reveals whether transcriptome-to-protein alignment is reliable
for checkpoint immunotherapy biomarker detection.

Analysis:
1. Identify PD-1+ and TIGIT+ cells from CITE-seq protein expression
2. Compare their cluster-resolved entropy to PD-1-/TIGIT- cells
3. Test whether high-entropy cells are enriched among checkpoint-high cells
4. Analyze CD4+PD-1+ (Treg-exhausted) vs CD8a+PD-1+ (effector-exhausted) subsets
5. Report AUROC for "is this cell a checkpoint-high cell" using H_cluster

Clinical relevance: cells where transcriptome and protein disagree on checkpoint
status are the most uncertain for immunotherapy response prediction.
"""
from __future__ import annotations
import json, sys, os
sys.path.insert(0, 'C:/Users/Maozer/projects/MOSAIC')

import numpy as np
import scipy.sparse as sp
import anndata as ad
from scipy import stats
from sklearn.metrics import roc_auc_score

from src.utils.paths import EXPERIMENTS_DIR, PROCESSED_DIR

EXP_ID = "exp001_citeseq"
DATASET = "citeseq_pbmc"

# Protein markers relevant to immunotherapy
CHECKPOINT_MARKERS = {
    "PD-1_TotalSeqB": "T cell exhaustion / checkpoint target",
    "TIGIT_TotalSeqB": "T cell exhaustion / checkpoint target",
    "CD25_TotalSeqB": "Treg / activated T (IL-2Ra)",
    "CD127_TotalSeqB": "Memory T survival / IL-7Ra",
}
LINEAGE_MARKERS = {
    "CD3_TotalSeqB": "T lineage",
    "CD4_TotalSeqB": "Helper/Treg",
    "CD8a_TotalSeqB": "Cytotoxic T",
}


def main():
    exp_dir = EXPERIMENTS_DIR / EXP_ID
    H_full = np.load(exp_dir / "alignment_entropy_cluster.npy")
    sub_idx = np.load(exp_dir / "ot_subsample_indices.npy")

    prot = ad.read_h5ad(PROCESSED_DIR / f"{DATASET}_atac.h5ad")
    rna = ad.read_h5ad(PROCESSED_DIR / f"{DATASET}_rna.h5ad")

    H = H_full  # 3000 cells (OT subsample)
    prot_sub = prot[sub_idx]
    rna_sub = rna[sub_idx]

    X_prot = prot_sub.X
    if sp.issparse(X_prot):
        X_prot = X_prot.toarray()

    marker_names = list(prot_sub.var_names)

    def get_marker(name):
        if name in marker_names:
            return X_prot[:, marker_names.index(name)]
        return None

    pd1 = get_marker("PD-1_TotalSeqB")
    tigit = get_marker("TIGIT_TotalSeqB")
    cd3 = get_marker("CD3_TotalSeqB")
    cd4 = get_marker("CD4_TotalSeqB")
    cd8a = get_marker("CD8a_TotalSeqB")
    cd25 = get_marker("CD25_TotalSeqB")
    cd127 = get_marker("CD127_TotalSeqB")

    n = len(H)
    results = {"experiment": "checkpoint_immunotherapy", "dataset": DATASET,
               "source_exp": EXP_ID, "n_cells": n, "analyses": []}

    # 1. PD-1 high vs low entropy
    for marker_name, vals, label in [
        ("PD-1_TotalSeqB", pd1, "PD-1+ exhausted T"),
        ("TIGIT_TotalSeqB", tigit, "TIGIT+ exhausted T"),
    ]:
        if vals is None:
            continue
        thresh = np.percentile(vals, 75)
        high_mask = vals >= thresh
        low_mask = vals < np.percentile(vals, 25)
        stat, pval = stats.mannwhitneyu(H[high_mask], H[low_mask], alternative="greater")
        try:
            auroc = float(roc_auc_score(high_mask.astype(int), H))
        except: auroc = float("nan")
        results["analyses"].append({
            "analysis": f"{marker_name}_entropy",
            "label": label,
            "n_high": int(high_mask.sum()),
            "n_low": int(low_mask.sum()),
            "mean_H_high": float(H[high_mask].mean()),
            "mean_H_low": float(H[low_mask].mean()),
            "H_ratio": float(H[high_mask].mean() / H[low_mask].mean()) if H[low_mask].mean() > 0 else float("nan"),
            "mannwhitney_p_greater": float(pval),
            "auroc_H_predicts_checkpoint_high": auroc,
            "clinical_interpretation": f"Higher entropy in checkpoint-high cells suggests transcriptome-proteome discordance for {label} identification",
        })

    # 2. PD-1+TIGIT+ double positive (exhausted) vs PD-1-TIGIT- (not exhausted)
    if pd1 is not None and tigit is not None:
        pd1_75 = np.percentile(pd1, 75)
        tigit_75 = np.percentile(tigit, 75)
        exhausted = (pd1 >= pd1_75) & (tigit >= tigit_75)
        fresh = (pd1 < np.percentile(pd1, 25)) & (tigit < np.percentile(tigit, 25))
        stat, pval = stats.mannwhitneyu(H[exhausted], H[fresh], alternative="greater")
        results["analyses"].append({
            "analysis": "PD1_TIGIT_double_positive_exhausted",
            "label": "Exhausted T (PD-1+TIGIT+ double-positive vs PD-1-TIGIT- fresh)",
            "n_exhausted": int(exhausted.sum()),
            "n_fresh": int(fresh.sum()),
            "mean_H_exhausted": float(H[exhausted].mean()),
            "mean_H_fresh": float(H[fresh].mean()),
            "H_ratio": float(H[exhausted].mean() / H[fresh].mean()) if H[fresh].mean() > 0 else float("nan"),
            "mannwhitney_p_greater": float(pval),
            "clinical_interpretation": "Exhausted T cells (checkpoint immunotherapy targets) have higher alignment uncertainty",
        })

    # 3. CD4+PD-1+ (Treg/exhausted CD4) vs CD8a+PD-1+ (effector-exhausted CD8)
    if pd1 is not None and cd4 is not None and cd8a is not None:
        pd1_med = np.percentile(pd1, 50)
        cd4_med = np.percentile(cd4, 50)
        cd8a_med = np.percentile(cd8a, 50)
        cd4_pd1 = (cd4 >= cd4_med) & (pd1 >= pd1_med) & (cd8a < cd8a_med)
        cd8_pd1 = (cd8a >= cd8a_med) & (pd1 >= pd1_med) & (cd4 < cd4_med)
        if cd4_pd1.sum() > 10 and cd8_pd1.sum() > 10:
            stat, pval = stats.mannwhitneyu(H[cd4_pd1], H[cd8_pd1], alternative="two-sided")
            results["analyses"].append({
                "analysis": "CD4_vs_CD8_PD1_entropy",
                "label": "CD4+PD-1+ vs CD8a+PD-1+ (Treg-exhausted vs effector-exhausted)",
                "n_cd4_pd1": int(cd4_pd1.sum()),
                "n_cd8_pd1": int(cd8_pd1.sum()),
                "mean_H_cd4_pd1": float(H[cd4_pd1].mean()),
                "mean_H_cd8_pd1": float(H[cd8_pd1].mean()),
                "mannwhitney_p_twosided": float(pval),
                "clinical_interpretation": "Differential uncertainty between CD4 and CD8 exhausted T cell subsets",
            })

    # 4. CD127+ memory T (favorable prognosis) vs CD127- effector/exhausted
    if cd127 is not None and cd3 is not None:
        cd3_med = np.percentile(cd3, 50)
        cd127_75 = np.percentile(cd127, 75)
        cd127_25 = np.percentile(cd127, 25)
        memory_t = (cd3 >= cd3_med) & (cd127 >= cd127_75)
        effector_t = (cd3 >= cd3_med) & (cd127 < cd127_25)
        if memory_t.sum() > 10 and effector_t.sum() > 10:
            stat, pval = stats.mannwhitneyu(H[memory_t], H[effector_t], alternative="greater")
            results["analyses"].append({
                "analysis": "CD127_memory_vs_effector_entropy",
                "label": "CD3+CD127+ memory T vs CD3+CD127- effector/exhausted T",
                "n_memory": int(memory_t.sum()),
                "n_effector": int(effector_t.sum()),
                "mean_H_memory": float(H[memory_t].mean()),
                "mean_H_effector": float(H[effector_t].mean()),
                "H_ratio": float(H[memory_t].mean() / H[effector_t].mean()) if H[effector_t].mean() > 0 else float("nan"),
                "mannwhitney_p_greater": float(pval),
                "clinical_interpretation": "Memory T (CD127+) may have higher proteome-transcriptome discordance than terminally differentiated effectors",
            })

    os.makedirs("experiments/checkpoint_immunotherapy", exist_ok=True)
    out_path = "experiments/checkpoint_immunotherapy/results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print("\n=== Checkpoint Immunotherapy Analysis ===")
    for a in results["analyses"]:
        print(f"\n{a['analysis']}:")
        print(f"  {a['label']}")
        if "mean_H_high" in a:
            print(f"  H checkpoint-high={a['mean_H_high']:.4f}  H checkpoint-low={a['mean_H_low']:.4f}  ratio={a.get('H_ratio','nan'):.2f}x")
            print(f"  Mann-Whitney p (high > low): {a['mannwhitney_p_greater']:.2e}")
            print(f"  AUROC (H predicts checkpoint-high): {a.get('auroc_H_predicts_checkpoint_high','nan'):.4f}")
        elif "mean_H_exhausted" in a:
            print(f"  H exhausted={a['mean_H_exhausted']:.4f}  H fresh={a['mean_H_fresh']:.4f}  ratio={a.get('H_ratio','nan'):.2f}x")
            print(f"  Mann-Whitney p: {a['mannwhitney_p_greater']:.2e}")
    print(f"\nSaved: {out_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
