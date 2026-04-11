#!/usr/bin/env python3
# MIT License - Bryan Cheng, 2026
# Part of MOSAIC - paper figure generation
"""Generate the main figures for the MOSAIC paper.

Figures:
  Fig1 — overview UMAP of aligned latent (RNA ∪ ATAC) colored by leiden cluster,
         plus a paired-matching inset showing per-cell distance histogram.
  Fig2 — cell-level vs cluster-level entropy calibration. Two-panel figure
         showing the cell-level scatter (with negative Spearman) and the
         cluster-level histogram of H values per cluster-correct vs wrong.
  Fig3 — missing-type detection bar chart. Per-cluster AUROC for each leave-out
         experiment on PBMC and brain, highlighting the consistency.
  Fig4 — baseline comparison grouped bar chart. FOSCTTM / LT / ARI across
         MOSAIC, NN-on-IB, SCOT, raw-features for PBMC 10k.

All figures use a consistent publication style (matplotlib 14pt labels,
300 dpi savefig, tight bounding box). Output to figures/*.png and *.pdf.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from src.utils.paths import EXPERIMENTS_DIR, FIGURES_DIR, PROCESSED_DIR

# Publication-quality rcParams
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def _save(fig, name: str):
    """Save figure as PNG and PDF in the FIGURES_DIR."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        path = FIGURES_DIR / f"{name}.{ext}"
        fig.savefig(path)
    plt.close(fig)
    print(f"[fig] saved figures/{name}.png & .pdf")


# ---------------------------------------------------------------------------
# Figure 1 — aligned UMAP (RNA ∪ ATAC)
# ---------------------------------------------------------------------------


def figure_1_aligned_umap(exp_name: str, dataset_id: str):
    """Joint UMAP of the aligned latent, colored by leiden cluster."""
    exp_dir = EXPERIMENTS_DIR / exp_name
    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac = np.load(exp_dir / "z_atac_aligned.npy")

    rna = ad.read_h5ad(PROCESSED_DIR / f"{dataset_id}_rna.h5ad")
    labels = rna.obs["cell_type"].astype(str).values
    # 2D projection of the joint 64-d latent via PCA (fast, deterministic)
    # Use one common PCA fitted on the concatenation.
    joint = np.concatenate([Z_rna, Z_atac], axis=0)
    pca = PCA(n_components=2, random_state=0)
    joint_2d = pca.fit_transform(joint)
    rna_2d = joint_2d[:len(Z_rna)]
    atac_2d = joint_2d[len(Z_rna):]

    unique = sorted(np.unique(labels), key=lambda c: int(c))
    cmap = plt.get_cmap("tab20", len(unique))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    titles = ["RNA latent (PC-projected)",
              "ATAC latent aligned to RNA",
              "Joint (both modalities)"]
    coords = [rna_2d, atac_2d, joint_2d]
    labels_list = [labels, labels, np.concatenate([labels, labels])]
    for ax, title, c, lab in zip(axes, titles, coords, labels_list):
        for i, cl in enumerate(unique):
            m = lab == cl
            ax.scatter(c[m, 0], c[m, 1], s=3, alpha=0.7,
                       c=[cmap(i)], label=cl, edgecolors="none")
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("PC 1")
        ax.set_ylabel("PC 2" if ax is axes[0] else "")

    # Legend on the right, compact
    handles, hlabels = axes[-1].get_legend_handles_labels()
    fig.legend(handles[:min(20, len(unique))], hlabels[:min(20, len(unique))],
               loc="center right", bbox_to_anchor=(1.03, 0.5),
               frameon=False, ncol=1, markerscale=2)
    plt.suptitle(f"Aligned latent — {dataset_id}", y=1.02)
    plt.tight_layout()
    _save(fig, f"fig1_aligned_latent_{dataset_id}")


# ---------------------------------------------------------------------------
# Figure 2 — cell-level vs cluster-level entropy
# ---------------------------------------------------------------------------


def figure_2_entropy_comparison(exp_name: str, dataset_id: str):
    """Side-by-side comparison of cell-level vs cluster-level entropy."""
    exp_dir = EXPERIMENTS_DIR / exp_name
    plan = np.load(exp_dir / "alignment_plan_subsample.npy")
    sub_idx = np.load(exp_dir / "ot_subsample_indices.npy")
    H_cell = np.load(exp_dir / "alignment_entropy_subsample.npy")
    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac_path = exp_dir / "z_atac_aligned.npy"
    Z_atac = np.load(Z_atac_path if Z_atac_path.exists() else exp_dir / "z_atac.npy")

    rna = ad.read_h5ad(PROCESSED_DIR / f"{dataset_id}_rna.h5ad")
    labels = rna.obs["cell_type"].astype(str).values[sub_idx]

    # Recompute cluster-level H
    rs = plan.sum(axis=1, keepdims=True); rs[rs == 0] = 1e-30
    P = plan / rs
    unique = np.unique(labels)
    P_cluster = np.zeros((P.shape[0], len(unique)))
    for i, c in enumerate(unique):
        P_cluster[:, i] = P[:, labels == c].sum(axis=1)
    log_K = float(np.log(max(len(unique), 2)))
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(P_cluster > 0, np.log(P_cluster), 0.0)
    H_cluster = -(P_cluster * logp).sum(axis=1) / log_K

    # Argmax cluster correctness
    argmax_cluster_idx = P_cluster.argmax(axis=1)
    argmax_cluster = unique[argmax_cluster_idx]
    correct_mask = argmax_cluster == labels
    correct_accuracy = correct_mask.mean()

    # True alignment distance (per-cell error)
    errs = np.sqrt(((Z_rna[sub_idx] - Z_atac[sub_idx]) ** 2).sum(axis=1))

    from scipy.stats import spearmanr
    rho_cell, _ = spearmanr(H_cell, errs)
    rho_cluster, _ = spearmanr(H_cluster, errs)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    ax = axes[0]
    ax.scatter(H_cell, errs, s=3, alpha=0.2, c="#3770B0", edgecolors="none")
    ax.set_xlabel("cell-level row entropy $H_{\\mathrm{cell}}$")
    ax.set_ylabel("distance to true partner")
    ax.set_title(f"Cell-level entropy (wrong sign)\n"
                 f"Spearman $\\rho$ = {rho_cell:.3f}")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.scatter(H_cluster, errs, s=3, alpha=0.2, c="#C03030", edgecolors="none")
    ax.set_xlabel("cluster-resolved entropy $H_{\\mathrm{cluster}}$")
    ax.set_ylabel("distance to true partner")
    ax.set_title(f"Cluster-resolved entropy\n"
                 f"Spearman $\\rho$ = {rho_cluster:.3f}")
    ax.grid(alpha=0.3)

    ax = axes[2]
    # Histogram of H_cluster for correct vs wrong cluster assignments
    bins = np.linspace(0, 1, 40)
    ax.hist(H_cluster[correct_mask], bins=bins, color="#5AA06A",
            alpha=0.7, label=f"argmax correct ({correct_mask.sum()})")
    ax.hist(H_cluster[~correct_mask], bins=bins, color="#C03030",
            alpha=0.7, label=f"argmax wrong ({(~correct_mask).sum()})")
    ax.set_xlabel("$H_{\\mathrm{cluster}}$")
    ax.set_ylabel("count (log scale)")
    ax.set_yscale("log")
    ax.legend(frameon=False)
    ax.set_title(f"cluster-argmax correctness vs $H_{{\\mathrm{{cluster}}}}$\n"
                 f"{correct_accuracy*100:.1f}% correct overall")

    plt.suptitle(f"MOSAIC alignment uncertainty — {dataset_id}", y=1.02)
    plt.tight_layout()
    _save(fig, f"fig2_entropy_comparison_{dataset_id}")


# ---------------------------------------------------------------------------
# Figure 3 — missing cell type detection AUROC bar chart
# ---------------------------------------------------------------------------


def figure_3_missing_type():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, (name, exp, dataset) in zip(axes, [
        ("PBMC 10k (β=0.01)", "exp001_pbmc_final", "pbmc10k_multiome"),
        ("Brain 5k (β=0.001)", "exp001_brain_beta0001", "brain3k_multiome"),
    ]):
        path = EXPERIMENTS_DIR / exp / "exp003_missing_type.json"
        if not path.exists():
            ax.text(0.5, 0.5, "not run", ha="center", va="center")
            continue
        with path.open() as f:
            data = json.load(f)
        per = data["per_cluster"]
        clusters = [p["target_cluster"] for p in per
                    if p.get("auroc_cluster_entropy") is not None]
        aurocs = [p["auroc_cluster_entropy"] for p in per
                  if p.get("auroc_cluster_entropy") is not None]
        order = np.argsort(aurocs)[::-1]
        clusters_s = [clusters[i] for i in order]
        aurocs_s = [aurocs[i] for i in order]
        bars = ax.bar(range(len(aurocs_s)), aurocs_s, color="#3770B0", edgecolor="black", linewidth=0.5)
        ax.axhline(0.5, color="gray", linestyle="--", label="random")
        ax.set_xticks(range(len(clusters_s)))
        ax.set_xticklabels([f"c{c}" for c in clusters_s], rotation=60, ha="right", fontsize=9)
        ax.set_ylim(0.5, 1.02)
        ax.set_title(f"{name} (mean AUROC = {data['mean_auroc']:.3f})")
        ax.set_ylabel("AUROC")
        ax.grid(axis="y", alpha=0.3)
    plt.suptitle("Exp 3 — detecting missing cell types via cluster entropy", y=1.02)
    plt.tight_layout()
    _save(fig, "fig3_missing_type_auroc")


# ---------------------------------------------------------------------------
# Figure 4 — baseline comparison
# ---------------------------------------------------------------------------


def figure_5_beta_tradeoff():
    """Figure showing the multi-seed beta comparison on both datasets."""
    import json as _json

    with (EXPERIMENTS_DIR / "aggregate_brain3k_multiome_beta0.01.json").open() as f:
        brain_b01 = _json.load(f)
    with (EXPERIMENTS_DIR / "aggregate_brain3k_multiome_beta0.001.json").open() as f:
        brain_b0001 = _json.load(f)
    with (EXPERIMENTS_DIR / "aggregate_pbmc10k_multiome.json").open() as f:
        pbmc_b01 = _json.load(f)
    with (EXPERIMENTS_DIR / "aggregate_pbmc10k_multiome_beta0.001.json").open() as f:
        pbmc_b0001 = _json.load(f)

    # Map (display title, aggregate key)
    metric_keys = [
        ("FOSCTTM\n(lower is better)", "foscttm_mean"),
        ("LT RNA→ATAC", "lt_rna_to_atac"),
        ("LT ATAC→RNA", "lt_atac_to_rna"),
        ("Joint ARI", "joint_ari"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(14, 4))

    for ax, (title, key) in zip(axes, metric_keys):
        pb01_mu = pbmc_b01[f"{key}_mean"]; pb01_err = pbmc_b01[f"{key}_std"]
        pb0001_mu = pbmc_b0001[f"{key}_mean"]; pb0001_err = pbmc_b0001[f"{key}_std"]
        br01_mu = brain_b01[f"{key}_mean"]; br01_err = brain_b01[f"{key}_std"]
        br0001_mu = brain_b0001[f"{key}_mean"]; br0001_err = brain_b0001[f"{key}_std"]

        x = np.arange(2)
        width = 0.35

        ax.bar(x - width/2, [pb01_mu, pb0001_mu], width,
               yerr=[pb01_err, pb0001_err], capsize=4,
               label="PBMC 10k (3 seeds)",
               color="#3770B0", edgecolor="black", linewidth=0.5)
        ax.bar(x + width/2, [br01_mu, br0001_mu], width,
               yerr=[br01_err, br0001_err], capsize=4,
               label="Brain 5k (3 seeds)",
               color="#5AA06A", edgecolor="black", linewidth=0.5)

        ax.set_xticks(x)
        ax.set_xticklabels(["β = 0.01", "β = 0.001"])
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
        if ax is axes[0]:
            ax.legend(frameon=False, loc="upper right", fontsize=9)

    plt.suptitle("MOSAIC β comparison (multi-seed) — β=0.001 wins or ties on every metric",
                 y=1.02, fontsize=14)
    plt.tight_layout()
    _save(fig, "fig5_beta_tradeoff")


def figure_4_baselines():
    """Baselines comparison with multi-seed error bars on MOSAIC.

    MOSAIC bar heights and error bars are computed from the aggregate JSONs
    (3-seed mean ± std). SCOT, NN-on-IB (single-seed ablation), and Raw
    baselines are single-seed and drawn without error bars — which is the
    honest rendering, since we only have one run of each baseline method.
    """
    import json as _json

    def _load_baselines(dataset):
        d = EXPERIMENTS_DIR / f"baselines_{dataset}"
        with (d / "baseline_results.json").open() as f:
            scot = _json.load(f)["scot"]["metrics"]
        with (d / "simple_baseline_results.json").open() as f:
            simple = _json.load(f)
        nn_m = simple["nn_on_ib"]["metrics"]
        raw_m = simple["raw_ot"]["metrics"]
        return scot, nn_m, raw_m

    def _load_aggregate(path):
        """Load mean ± std from an aggregate seeds JSON."""
        with path.open() as f:
            d = _json.load(f)
        return {
            "foscttm_mean": d["foscttm_mean_mean"],
            "foscttm_std": d["foscttm_mean_std"],
            "lt_a_b_mean": d["lt_rna_to_atac_mean"],
            "lt_a_b_std": d["lt_rna_to_atac_std"],
            "lt_b_a_mean": d["lt_atac_to_rna_mean"],
            "lt_b_a_std": d["lt_atac_to_rna_std"],
            "ari_mean": d["joint_ari_mean"],
            "ari_std": d["joint_ari_std"],
            "n_seeds": d["n_seeds"],
        }

    pbmc_scot, pbmc_nn, pbmc_raw = _load_baselines("pbmc10k_multiome")
    brain_scot, brain_nn, brain_raw = _load_baselines("brain3k_multiome")

    # MOSAIC multi-seed: both at β=0.001 (paper default after 10-seed resolution).
    # PBMC uses the 10-seed aggregate; Brain uses the 3-seed aggregate (already tight).
    pbmc_mosaic = _load_aggregate(EXPERIMENTS_DIR / "aggregate_pbmc10k_multiome_beta0001_10seed.json")
    brain_mosaic = _load_aggregate(EXPERIMENTS_DIR / "aggregate_brain3k_multiome_beta0.001.json")

    methods = ["MOSAIC\n(ours)", "NN on IB\n(no OT)", "SCOT\n(GW)", "Raw PCA/LSI\n(no IB)"]

    def _series(mosaic_agg, scot, nn_m, raw_m):
        # Heights and stds; std=0 for single-seed baselines so no error bar shows.
        foscttm_h = [1 - mosaic_agg["foscttm_mean"],
                     1 - nn_m["foscttm_mean"],
                     1 - scot["foscttm"]["foscttm_mean"],
                     1 - raw_m["foscttm_mean"]]
        foscttm_err = [mosaic_agg["foscttm_std"], 0.0, 0.0, 0.0]
        lt_a_b_h = [mosaic_agg["lt_a_b_mean"],
                    nn_m["label_transfer_rna_to_atac"],
                    scot["label_transfer_rna_to_atac"],
                    raw_m["label_transfer_rna_to_atac"]]
        lt_a_b_err = [mosaic_agg["lt_a_b_std"], 0.0, 0.0, 0.0]
        lt_b_a_h = [mosaic_agg["lt_b_a_mean"],
                    nn_m["label_transfer_atac_to_rna"],
                    scot["label_transfer_atac_to_rna"],
                    raw_m["label_transfer_atac_to_rna"]]
        lt_b_a_err = [mosaic_agg["lt_b_a_std"], 0.0, 0.0, 0.0]
        ari_h = [mosaic_agg["ari_mean"],
                 nn_m["joint_clustering_ari"],
                 scot["joint_clustering_ari"],
                 raw_m["joint_clustering_ari"]]
        ari_err = [mosaic_agg["ari_std"], 0.0, 0.0, 0.0]
        return {
            "foscttm_h": foscttm_h, "foscttm_err": foscttm_err,
            "lt_a_b_h": lt_a_b_h, "lt_a_b_err": lt_a_b_err,
            "lt_b_a_h": lt_b_a_h, "lt_b_a_err": lt_b_a_err,
            "ari_h": ari_h, "ari_err": ari_err,
        }

    pbmc = _series(pbmc_mosaic, pbmc_scot, pbmc_nn, pbmc_raw)
    brain = _series(brain_mosaic, brain_scot, brain_nn, brain_raw)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    x = np.arange(len(methods))
    width = 0.2

    titles = [
        f"PBMC 10k (MOSAIC β=0.001, n={pbmc_mosaic['n_seeds']} seeds)",
        f"Brain 5k (MOSAIC β=0.001, n={brain_mosaic['n_seeds']} seeds)",
    ]
    err_kw = dict(capsize=3, ecolor="#303030", elinewidth=1)
    for ax, data, title in zip(axes, [pbmc, brain], titles):
        ax.bar(x - 1.5*width, data["foscttm_h"], width,
               yerr=data["foscttm_err"], label="1 − FOSCTTM",
               color="#3770B0", error_kw=err_kw)
        ax.bar(x - 0.5*width, data["lt_a_b_h"], width,
               yerr=data["lt_a_b_err"], label="LT RNA→ATAC",
               color="#5AA06A", error_kw=err_kw)
        ax.bar(x + 0.5*width, data["lt_b_a_h"], width,
               yerr=data["lt_b_a_err"], label="LT ATAC→RNA",
               color="#C07030", error_kw=err_kw)
        ax.bar(x + 1.5*width, data["ari_h"], width,
               yerr=data["ari_err"], label="joint ARI",
               color="#8860B0", error_kw=err_kw)
        ax.set_xticks(x)
        ax.set_xticklabels(methods)
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
        ax.set_ylim(0, 1.05)
        if ax is axes[0]:
            ax.set_ylabel("metric (higher is better)")
            ax.legend(frameon=False, loc="upper right", fontsize=10)

    plt.suptitle("MOSAIC vs baselines (MOSAIC bars: mean ± std across seeds; "
                 "baselines: single-seed)",
                 y=1.02, fontsize=13)
    plt.tight_layout()
    _save(fig, "fig4_baselines_both")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    figure_1_aligned_umap("exp001_pbmc_final", "pbmc10k_multiome")
    figure_1_aligned_umap("exp001_brain_final", "brain3k_multiome")
    figure_2_entropy_comparison("exp001_pbmc_final", "pbmc10k_multiome")
    figure_2_entropy_comparison("exp001_brain_final", "brain3k_multiome")
    figure_3_missing_type()
    figure_4_baselines()
    figure_5_beta_tradeoff()
    print("\nAll figures generated.")
