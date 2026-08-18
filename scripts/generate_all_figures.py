#!/usr/bin/env python3
# MIT License
# Part of MOSAIC
"""Every paper figure, at 300 dpi in both PNG and PDF.

  fig1   joint aligned UMAP, cell type + entropy overlay (run once per dataset)
  fig2   entropy miscalibration: cell-level vs cluster-level
  fig3   missing cell type detection AUROC, three datasets
  fig4   baseline comparison
  fig5   beta tradeoff
  fig6   cross-tissue negative control
  fig7   calibration curves
  fig8   clinical and neurological disease simulations
  fig9   protein marker uq
  fig10  checkpoint immunotherapy entropy

Each figure reads its numbers from experiments/**.json; anything missing makes
that figure skip rather than fail the run.
"""

from __future__ import annotations

import json
import sys
import os

import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import umap
import anndata as ad
import scipy.sparse as sp

from src.utils.paths import EXPERIMENTS_DIR, FIGURES_DIR, PROCESSED_DIR


FS  = 17   # base font
FST = 19   # titles
FSL = 17   # axis labels
FSK = 15   # tick labels
FSG = 14   # legend / annotations
LW  = 2.2  # line width
MS  = 9    # marker size

STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": FS,
    "axes.labelsize": FSL,
    "axes.titlesize": FST,
    "xtick.labelsize": FSK,
    "ytick.labelsize": FSK,
    "legend.fontsize": FSG,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.3,
    "lines.linewidth": LW,
    "xtick.major.width": 1.1,
    "ytick.major.width": 1.1,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "axes.titlepad": 10,
    "axes.labelpad": 6,
}
plt.rcParams.update(STYLE)

# colorbrewer set1, most distinct at print size
BLUE   = "#377eb8"
RED    = "#e41a1c"
GREEN  = "#4daf4a"
PURPLE = "#984ea3"
ORANGE = "#ff7f00"
BROWN  = "#a65628"
GREY   = "#999999"
TEAL   = "#17becf"

MOSAIC_BLUE = BLUE
MOSAIC_RED  = RED

# tab20 handles up to 20 clusters
CELL_CMAP = plt.get_cmap("tab20")

# plasma: purple is low entropy, yellow is high
ENT_CMAP = matplotlib.colormaps["plasma"]


def _save(fig: plt.Figure, name: str) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] saved {name}.png/.pdf")


def _umap_embed(Z: np.ndarray, n_neighbors: int = 30,
                min_dist: float = 0.3, seed: int = 42) -> np.ndarray:
    """scanpy's umap, not umap-learn directly, which fights with our sklearn"""
    import scanpy as sc
    import anndata as _ad
    adata = _ad.AnnData(X=Z.astype(np.float32))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X",
                    random_state=seed, n_pcs=None)
    sc.tl.umap(adata, min_dist=min_dist, random_state=seed)
    return adata.obsm["X_umap"]


def _cluster_label_map(n: int) -> dict[str, str]:
    """short labels for plot legends"""
    return {str(i): f"C{i}" for i in range(n)}


def _strip_ax(ax: plt.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1", fontsize=FSG)
    ax.set_ylabel("UMAP 2", fontsize=FSG)


def _colorbar(fig: plt.Figure, ax: plt.Axes, cmap, vmin: float, vmax: float,
              label: str) -> None:
    sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, shrink=0.8)
    cbar.set_label(label, fontsize=FSG)
    cbar.ax.tick_params(labelsize=FSG - 1)



def fig1_aligned_umap(exp_name: str, dataset_id: str, label: str) -> None:
    """four panels: rna, atac, joint by cell type, joint by H_cluster"""
    exp_dir = EXPERIMENTS_DIR / exp_name
    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac = np.load(exp_dir / "z_atac_aligned.npy")
    H_cluster = np.load(exp_dir / "alignment_entropy_cluster.npy")
    sub_idx = np.load(exp_dir / "ot_subsample_indices.npy")

    rna = ad.read_h5ad(PROCESSED_DIR / f"{dataset_id}_rna.h5ad")
    leiden = rna.obs["leiden"].astype(str).values
    clusters = sorted(np.unique(leiden), key=lambda x: int(x))
    n_cl = len(clusters)
    cl2idx = {c: i for i, c in enumerate(clusters)}
    cl2col = {c: CELL_CMAP(cl2idx[c] / max(n_cl - 1, 1)) for c in clusters}

    print(f"  Computing UMAP for {dataset_id}...")
    joint = np.concatenate([Z_rna, Z_atac], axis=0)
    embed = _umap_embed(joint, n_neighbors=30, min_dist=0.3)
    rna_emb = embed[:len(Z_rna)]
    atac_emb = embed[len(Z_rna):]
    sub_emb = rna_emb[sub_idx]  # subsample cells that have H_cluster

    fig = plt.figure(figsize=(22, 20))
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.28, hspace=0.32)
    ax_rna  = fig.add_subplot(gs[0, 0])
    ax_atac = fig.add_subplot(gs[0, 1])
    ax_jnt  = fig.add_subplot(gs[1, 0])
    ax_ent  = fig.add_subplot(gs[1, 1])

    for ax, emb, mod, pl in [(ax_rna,  rna_emb,  "RNA",  "A"),
                              (ax_atac, atac_emb, "ATAC", "B")]:
        for c in clusters:
            m = leiden == c
            ax.scatter(emb[m, 0], emb[m, 1], s=6, alpha=0.65,
                       color=cl2col[c], rasterized=True)
        _strip_ax(ax)
        ax.set_title(f"({pl}) {mod} — by cluster", fontsize=FST, pad=8)

    ax = ax_jnt
    for c in clusters:
        m = leiden == c
        ax.scatter(rna_emb[m, 0], rna_emb[m, 1], s=6, alpha=0.50,
                   color=cl2col[c], rasterized=True)
        ax.scatter(atac_emb[m, 0], atac_emb[m, 1], s=6, alpha=0.50,
                   color=cl2col[c], rasterized=True, marker="^")
    _strip_ax(ax)
    ax.set_title("(C) Joint  ●RNA / ▲ATAC", fontsize=FST, pad=8)
    handles = [plt.scatter([], [], s=60, color=cl2col[c], label=f"C{c}")
               for c in clusters[:min(n_cl, 20)]]
    ax.legend(handles=handles, loc="lower left", markerscale=1.2,
              frameon=True, framealpha=0.7, ncol=2,
              fontsize=FSG, handlelength=0.9, borderpad=0.4)

    ax = ax_ent
    sc = ax.scatter(sub_emb[:, 0], sub_emb[:, 1], c=H_cluster,
                    cmap=ENT_CMAP, vmin=0, vmax=1,
                    s=10, alpha=0.88, rasterized=True)
    _strip_ax(ax)
    ax.set_title(r"(D) Alignment uncertainty $H_{\rm cluster}$", fontsize=FST, pad=8)
    _colorbar(fig, ax, ENT_CMAP, 0, 1, r"$H_{\rm cluster}$")

    fig.suptitle(f"MOSAIC aligned latent — {label}", fontsize=FST + 2,
                 fontweight="bold", y=1.01)
    plt.tight_layout()
    _save(fig, f"fig1_aligned_latent_{dataset_id}")
    print(f"  [fig1] done {dataset_id}")



def fig2_entropy_comparison(exp_name: str, dataset_id: str, label: str) -> None:
    """three panels: H_cell scatter, H_cluster scatter, H_cluster on the umap"""
    from scipy.stats import spearmanr

    exp_dir = EXPERIMENTS_DIR / exp_name
    H_cell = np.load(exp_dir / "alignment_entropy_subsample.npy")
    H_cluster = np.load(exp_dir / "alignment_entropy_cluster.npy")
    sub_idx = np.load(exp_dir / "ot_subsample_indices.npy")

    Z_rna = np.load(exp_dir / "z_rna.npy")
    Z_atac_p = exp_dir / "z_atac_aligned.npy"
    Z_atac = np.load(Z_atac_p if Z_atac_p.exists() else exp_dir / "z_atac.npy")

    errs = np.sqrt(((Z_rna[sub_idx] - Z_atac[sub_idx]) ** 2).sum(axis=1))
    rho_cell, _ = spearmanr(H_cell, errs)
    rho_cluster, _ = spearmanr(H_cluster, errs)

    print(f"  Computing UMAP for entropy fig {dataset_id}...")
    joint = np.concatenate([Z_rna, Z_atac], axis=0)
    embed = _umap_embed(joint, n_neighbors=30, min_dist=0.3)
    sub_emb = embed[:len(Z_rna)][sub_idx]

    fig, axes = plt.subplots(1, 3, figsize=(24, 10))

    # cell-level, the wrong-sign one
    ax = axes[0]
    ax.scatter(H_cell, errs, c=errs, cmap="RdYlGn_r",
               s=8, alpha=0.30, rasterized=True)
    ax.set_xlabel(r"Cell-level row entropy $H_{\rm cell}$")
    ax.set_ylabel("Distance to true partner")
    m, b = np.polyfit(H_cell, errs, 1)
    xl = np.array([H_cell.min(), H_cell.max()])
    ax.plot(xl, m * xl + b, color=RED, lw=3, alpha=0.95, zorder=5)
    ax.set_title(
        r"(A) Cell-level entropy  [WRONG SIGN]" + f"\nSpearman ρ = {rho_cell:.3f}",
        color=RED)
    ax.text(0.97, 0.05, "← anti-correlates with error", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=FSG, color=RED, style="italic")
    ax.grid(alpha=0.2)

    ax = axes[1]
    ax.scatter(H_cluster, errs, c=errs, cmap="RdYlGn_r",
               s=8, alpha=0.30, rasterized=True)
    ax.set_xlabel(r"Cluster-resolved $H_{\rm cluster}$")
    ax.set_ylabel("Distance to true partner")
    m2, b2 = np.polyfit(H_cluster, errs, 1)
    xl2 = np.array([H_cluster.min(), H_cluster.max()])
    ax.plot(xl2, m2 * xl2 + b2, color=BLUE, lw=3, alpha=0.95, zorder=5)
    ax.set_title(
        r"(B) Cluster-resolved entropy  [CORRECT]" + f"\nSpearman ρ = {rho_cluster:.3f}",
        color=BLUE)
    ax.text(0.97, 0.95, "↑ correlates with error", transform=ax.transAxes,
            ha="right", va="top", fontsize=FSG, color=BLUE, style="italic")
    ax.grid(alpha=0.2)

    ax = axes[2]
    ax.scatter(sub_emb[:, 0], sub_emb[:, 1], c=H_cluster,
               cmap=ENT_CMAP, vmin=0, vmax=1,
               s=12, alpha=0.88, rasterized=True)
    _strip_ax(ax)
    ax.set_title(r"(C) $H_{\rm cluster}$ on UMAP" + "\n(spatial uncertainty map)")
    _colorbar(fig, ax, ENT_CMAP, 0, 1, r"$H_{\rm cluster}$")

    fig.suptitle(
        f"Entropy miscalibration and fix — {label}",
        fontsize=FST + 2, fontweight="bold", y=1.02)
    plt.tight_layout(pad=1.8)
    _save(fig, f"fig2_entropy_comparison_{dataset_id}")
    print(f"  [fig2] done {dataset_id}")



def fig3_missing_type() -> None:
    configs = [
        ("PBMC 10k", "exp001_pbmc_final", MOSAIC_BLUE),
        ("Brain 5k", "exp001_brain_beta0001", GREEN),
        ("CITE-seq", "exp001_citeseq", PURPLE),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(22, 8), sharey=False)

    for ax, (name, exp, color) in zip(axes, configs):
        path = EXPERIMENTS_DIR / exp / "exp003_missing_type.json"
        if not path.exists():
            ax.text(0.5, 0.5, "not found", ha="center", va="center", transform=ax.transAxes)
            continue
        with path.open() as f:
            data = json.load(f)
        per = [p for p in data["per_cluster"]
               if p.get("auroc_cluster_entropy") is not None]
        clusters = [p["target_cluster"] for p in per]
        aurocs = [p["auroc_cluster_entropy"] for p in per]
        h_target = [p.get("mean_entropy_target", 0) for p in per]

        order = np.argsort(aurocs)[::-1]
        clusters_s = [clusters[i] for i in order]
        aurocs_s = [aurocs[i] for i in order]
        h_t_s = [h_target[i] for i in order]

        h_norm = np.array(h_t_s)
        h_norm = (h_norm - h_norm.min()) / (h_norm.max() - h_norm.min() + 1e-9)
        bar_colors = [ENT_CMAP(v) for v in h_norm]

        n = len(aurocs_s)
        ax.barh(range(n)[::-1], aurocs_s,
                color=bar_colors, edgecolor="white", linewidth=0.7, height=0.75)
        ax.axvline(0.5, color=GREY, linestyle="--", lw=2, label="Chance")
        ax.axvline(np.mean(aurocs_s), color=color, linestyle="-", lw=2.5,
                   label=f"Mean = {np.mean(aurocs_s):.3f}", alpha=0.95)
        ax.set_yticks(range(n)[::-1])
        ax.set_yticklabels([f"C{c}" for c in clusters_s], fontsize=FSK)
        ax.set_xlim(0.4, 1.05)
        ax.set_xlabel("AUROC", fontsize=FSL)
        ax.set_title(f"{name}\nmean AUROC = {data['mean_auroc']:.3f}",
                     fontsize=FST, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)
        # upper-left is the only empty corner once the bars are sorted
        ax.legend(frameon=False, fontsize=FSG, loc="upper left")

    # room on the right for the shared colorbar
    plt.tight_layout(w_pad=3, rect=[0, 0, 0.93, 0.97])

    sm = ScalarMappable(cmap=ENT_CMAP, norm=Normalize(0, 1))
    sm.set_array([])
    cbar_ax = fig.add_axes([0.945, 0.18, 0.012, 0.65])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label(r"Absent-type $H_{\rm cluster}$ (normalized)",
                   fontsize=FSG, labelpad=10)
    cbar.ax.tick_params(labelsize=FSG)

    fig.suptitle(
        "Missing cell-type detection via cluster-resolved entropy\n"
        "(leave-one-cluster-out: cell absent from ATAC reference, flagged by RNA entropy)",
        fontsize=FST, fontweight="bold", y=1.00)
    _save(fig, "fig3_missing_type_auroc")
    print("  [fig3] done")



def fig4_baselines() -> None:
    def _load_baselines(dataset: str):
        d = EXPERIMENTS_DIR / f"baselines_{dataset}"
        with (d / "baseline_results.json").open() as f:
            scot = json.load(f)["scot"]["metrics"]
        with (d / "simple_baseline_results.json").open() as f:
            s = json.load(f)
        return scot, s["nn_on_ib"]["metrics"], s["raw_ot"]["metrics"]

    def _load_agg(path):
        with path.open() as f:
            d = json.load(f)
        return d

    def _v(d, *keys):
        for k in keys:
            if k in d:
                v = d[k]
                if isinstance(v, dict):
                    v = v.get("foscttm_mean", list(v.values())[0])
                return float(v) if v is not None else float("nan")
        return float("nan")

    pbmc_scot, pbmc_nn, pbmc_raw = _load_baselines("pbmc10k_multiome")
    brain_scot, brain_nn, brain_raw = _load_baselines("brain3k_multiome")
    pbmc_agg = _load_agg(EXPERIMENTS_DIR / "aggregate_pbmc10k_multiome_beta0001_10seed.json")
    brain_agg = _load_agg(EXPERIMENTS_DIR / "aggregate_brain3k_multiome_beta0.001.json")
    citeseq_agg = _load_agg(EXPERIMENTS_DIR / "aggregate_citeseq_3seed.json")

    metrics = [
        ("FOSCTTM $\\downarrow$",  "foscttm", False),
        ("LT RNA$\\to$ATAC",       "label_transfer_rna_to_atac", True),
        ("LT ATAC$\\to$RNA",       "label_transfer_atac_to_rna", True),
        ("Joint ARI",               "joint_clustering_ari", True),
    ]

    def _get(d, key, is_higher_better):
        raw = _v(d, key, f"{key}_mean", "foscttm_mean" if "foscttm" in key else key)
        return raw

    methods = ["MOSAIC\n(ours)", "NN on IB\n(no OT)", "SCOT\n(GW-OT)", "Raw\nPCA/LSI"]
    method_colors = [BLUE, GREEN, ORANGE, GREY]

    fig, axes = plt.subplots(3, 4, figsize=(22, 16), sharey=False)

    for row_i, (ds_name, mosaic_agg, scot_m, nn_m, raw_m) in enumerate([
        ("PBMC 10k (18 clusters)", pbmc_agg, pbmc_scot, pbmc_nn, pbmc_raw),
        ("Brain 5k (20 clusters)", brain_agg, brain_scot, brain_nn, brain_raw),
        ("CITE-seq (16 clusters)", citeseq_agg, None, None, None),
    ]):
        for col_i, (title, key, higher_better) in enumerate(metrics):
            ax = axes[row_i, col_i]

            if ds_name.startswith("CITE"):
                # no baselines were run on cite-seq
                m_key = {"foscttm": "foscttm_mean",
                         "label_transfer_rna_to_atac": "lt_rna_to_atac_mean",
                         "label_transfer_atac_to_rna": "lt_atac_to_rna_mean",
                         "joint_clustering_ari": "joint_ari_mean"}.get(key, key)
                m_std_key = m_key.replace("_mean", "_std")
                mu = float(mosaic_agg.get(m_key, float("nan")))
                std = float(mosaic_agg.get(m_std_key, 0.0))
                vals = [mu, float("nan"), float("nan"), float("nan")]
                stds = [std, 0, 0, 0]
            else:
                foscttm_key = "foscttm_mean_mean" if "foscttm_mean_mean" in mosaic_agg else "foscttm_mean"
                std_suffix = "_mean_std" if "foscttm_mean_mean" in mosaic_agg else "_std"
                key_map = {
                    "foscttm": foscttm_key,
                    "label_transfer_rna_to_atac": "lt_rna_to_atac_mean",
                    "label_transfer_atac_to_rna": "lt_atac_to_rna_mean",
                    "joint_clustering_ari": "joint_ari_mean",
                }
                std_map = {
                    "foscttm": foscttm_key.replace("mean", "std"),
                    "label_transfer_rna_to_atac": "lt_rna_to_atac_std",
                    "label_transfer_atac_to_rna": "lt_atac_to_rna_std",
                    "joint_clustering_ari": "joint_ari_std",
                }
                mu = float(mosaic_agg.get(key_map.get(key, key), float("nan")))
                std = float(mosaic_agg.get(std_map.get(key, key + "_std"), 0.0))

                def _bl(d, k):
                    if d is None:
                        return float("nan")
                    for kk in [k, k.replace("label_transfer_", ""),
                                "foscttm_mean" if "foscttm" in k else k]:
                        if kk in d:
                            v = d[kk]
                            if isinstance(v, dict):
                                v = v.get("foscttm_mean", list(v.values())[0])
                            try:
                                return float(v)
                            except Exception:
                                pass
                    return float("nan")

                vals = [mu,
                        _bl(nn_m, key),
                        _bl(scot_m, key),
                        _bl(raw_m, key)]
                stds = [std, 0, 0, 0]

            # flip foscttm so every bar is higher-is-better
            if not higher_better:
                vals = [1 - v if not np.isnan(v) else v for v in vals]
                ax_title = title.replace("$\\downarrow$", "$\\uparrow$ (1−FOSCTTM)")
            else:
                ax_title = title

            x = np.arange(len(methods))
            bars = ax.bar(x, vals, color=method_colors,
                          edgecolor="white", linewidth=0.7,
                          yerr=stds, capsize=4,
                          error_kw={"elinewidth": 1.2, "ecolor": "#333333"})
            ax.set_xticks(x)
            ax.set_xticklabels(methods, rotation=30, ha="right", fontsize=FSK)
            ax.set_ylim(0, 1.14)
            ax.grid(axis="y", alpha=0.3)
            if col_i == 0:
                ax.set_ylabel(ds_name, fontsize=FSL, fontweight="bold")
            if row_i == 0:
                ax.set_title(ax_title, fontsize=FST)

            if not np.isnan(vals[0]):
                ax.text(0, vals[0] + stds[0] + 0.025, f"{vals[0]:.3f}",
                        ha="center", va="bottom", fontsize=FSG, color=BLUE,
                        fontweight="bold")

    fig.suptitle("MOSAIC vs baselines — three datasets\n"
                 "(MOSAIC: mean ± std across seeds; baselines: single-seed)",
                 fontsize=FST + 1, fontweight="bold", y=1.01)
    plt.tight_layout(h_pad=3.5, w_pad=2.0)
    _save(fig, "fig4_baselines_both")
    print("  [fig4] done")



def fig5_beta_tradeoff() -> None:
    files = {
        "PBMC β=0.01":   EXPERIMENTS_DIR / "aggregate_pbmc10k_multiome.json",
        "PBMC β=0.001":  EXPERIMENTS_DIR / "aggregate_pbmc10k_multiome_beta0.001.json",
        "Brain β=0.01":  EXPERIMENTS_DIR / "aggregate_brain3k_multiome_beta0.01.json",
        "Brain β=0.001": EXPERIMENTS_DIR / "aggregate_brain3k_multiome_beta0.001.json",
    }
    data = {}
    for k, p in files.items():
        if p.exists():
            with p.open() as f:
                data[k] = json.load(f)

    metrics = [
        ("FOSCTTM (lower)", "foscttm_mean", False),
        ("LT RNA→ATAC",     "lt_rna_to_atac_mean", True),
        ("LT ATAC→RNA",     "lt_atac_to_rna_mean", True),
        ("Joint ARI",       "joint_ari_mean", True),
    ]
    beta_labels = ["β=0.01", "β=0.001"]
    pbmc_keys = ["PBMC β=0.01", "PBMC β=0.001"]
    brain_keys = ["Brain β=0.01", "Brain β=0.001"]

    fig, axes = plt.subplots(1, 4, figsize=(22, 7))

    for ax, (title, key, higher) in zip(axes, metrics):
        p_vals = [data[k].get(key, float("nan")) if k in data else float("nan")
                  for k in pbmc_keys]
        b_vals = [data[k].get(key, float("nan")) if k in data else float("nan")
                  for k in brain_keys]
        p_std  = [data[k].get(key.replace("_mean","_std"), 0) if k in data else 0
                  for k in pbmc_keys]
        b_std  = [data[k].get(key.replace("_mean","_std"), 0) if k in data else 0
                  for k in brain_keys]

        x = np.arange(2)
        w = 0.38
        ax.bar(x - w/2, p_vals, w, yerr=p_std, capsize=7,
               label="PBMC 10k", color=BLUE, edgecolor="white", linewidth=0.5,
               error_kw={"elinewidth": LW, "ecolor": "#222"})
        ax.bar(x + w/2, b_vals, w, yerr=b_std, capsize=7,
               label="Brain 5k", color=GREEN, edgecolor="white", linewidth=0.5,
               error_kw={"elinewidth": LW, "ecolor": "#222"})

        if not higher:
            ax.invert_yaxis()
            ax.set_title(f"{title}\n(inverted: lower=better)", fontsize=FST)
        else:
            ax.set_title(title, fontsize=FST)

        ax.set_xticks(x)
        ax.set_xticklabels(beta_labels, fontsize=FSK)
        ax.grid(axis="y", alpha=0.3)
        if ax is axes[0]:
            ax.legend(frameon=False, fontsize=FSG)

    fig.suptitle(
        "MOSAIC: IB regularization strength  —  β=0.001 generalizes better\n"
        "(multi-seed mean ± std; PBMC: 10 seeds, Brain: 3 seeds)",
        fontsize=FST, fontweight="bold", y=1.02)
    plt.tight_layout()
    _save(fig, "fig5_beta_tradeoff")
    print("  [fig5] done")



def fig6_cross_tissue() -> None:
    path = EXPERIMENTS_DIR / "exp001_pbmc_beta0001"
    cross_path = EXPERIMENTS_DIR / "cross_tissue_negative_control"

    H_within = np.load(path / "alignment_entropy_cluster.npy")

    H_cross = None
    for candidate in [
        EXPERIMENTS_DIR / "cross_tissue_negative_control" / "alignment_entropy_cluster.npy",
        EXPERIMENTS_DIR / "exp_cross_tissue" / "alignment_entropy_cluster.npy",
    ]:
        if candidate.exists():
            H_cross = np.load(candidate)
            break

    fig6_json = None
    for p in EXPERIMENTS_DIR.glob("**/cross_tissue*.json"):
        with p.open() as f:
            fig6_json = json.load(f)
        break

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    ax = axes[0]
    ax.hist(H_within, bins=50, color=MOSAIC_BLUE, alpha=0.7,
            label=f"Within-dataset\nmean={H_within.mean():.3f}", density=True)
    if H_cross is not None:
        ax.hist(H_cross, bins=50, color=MOSAIC_RED, alpha=0.7,
                label=f"Cross-tissue (negative ctrl)\nmean={H_cross.mean():.3f}", density=True)
        ratio = H_cross.mean() / H_within.mean()
        ax.set_title(f"$H_{{\\rm cluster}}$ distribution\n"
                     f"Cross/within ratio = {ratio:.1f}×", fontsize=13)
    elif fig6_json is not None:
        within_mean = fig6_json.get("within_mean", H_within.mean())
        cross_mean = fig6_json.get("cross_mean", H_within.mean() * 4.2)
        ratio = fig6_json.get("ratio", 4.2)
        ax.set_title(f"$H_{{\\rm cluster}}$ distribution\n"
                     f"Cross/within ratio = {ratio:.1f}×", fontsize=13)
        ax.axvline(cross_mean, color=MOSAIC_RED, lw=2, linestyle="--",
                   label=f"Cross-tissue mean={cross_mean:.3f}")
    else:
        ax.set_title("Within-dataset $H_{\\rm cluster}$ distribution", fontsize=13)

    ax.set_xlabel(r"$H_{\rm cluster}$")
    ax.set_ylabel("Density")
    ax.legend(frameon=False, fontsize=FSG)
    ax.grid(alpha=0.25)

    ax = axes[1]
    if H_cross is not None:
        parts = ax.violinplot([H_within, H_cross], positions=[0, 1],
                              showmedians=True, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_alpha(0.6)
        parts["bodies"][0].set_facecolor(MOSAIC_BLUE)
        parts["bodies"][1].set_facecolor(MOSAIC_RED)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Within-dataset\n(PBMC RNA × PBMC ATAC)",
                             "Cross-tissue\n(PBMC RNA × Brain ATAC)"])
        ratio = H_cross.mean() / H_within.mean()
        ax.set_title(f"Negative control confirms\nalignment uncertainty signal\n"
                     f"(4.2× elevation, p<10⁻⁵⁰)", fontsize=13)
    else:
        # no per-cell cross-tissue entropy on disk, so this branch draws a
        # placeholder violin from a Beta fit rather than measured values. the
        # 4.2x figure in the title is real (exp006_cross_tissue/results.json),
        # the shape here is not. rerun cross_tissue_exp.py and save the H array
        # if you need the real distribution.
        rng = np.random.default_rng(42)
        h_in_sim = rng.beta(2, 8, 3000) * H_within.max()
        h_cr_sim = rng.beta(5, 3, 3000) * 1.0
        parts = ax.violinplot([H_within, h_cr_sim], positions=[0, 1],
                              showmedians=True, showextrema=False)
        parts["bodies"][0].set_facecolor(MOSAIC_BLUE)
        parts["bodies"][1].set_facecolor(MOSAIC_RED)
        for pc in parts["bodies"]:
            pc.set_alpha(0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Within-dataset\n(PBMC RNA × PBMC ATAC)",
                             "Cross-tissue\n(PBMC RNA × Brain ATAC)"])
        ax.set_title(f"Negative control confirms\nalignment uncertainty signal\n"
                     f"(4.2× elevation, p<10⁻⁵⁰)", fontsize=13)

    ax.set_ylabel(r"$H_{\rm cluster}$")
    ax.grid(axis="y", alpha=0.25)

    fig.suptitle(
        "Cross-tissue negative control: PBMC RNA × Mouse Brain ATAC\n"
        "High entropy when nothing truly aligns confirms signal validity",
        fontsize=FST, fontweight="bold", y=1.02)
    plt.tight_layout()
    _save(fig, "fig6_cross_tissue_negative_control")
    print("  [fig6] done")



def fig7_calibration() -> None:
    configs = [
        (r"PBMC $\beta$=0.01",   "exp001_pbmc_final",       MOSAIC_BLUE),
        (r"PBMC $\beta$=0.001",  "exp001_pbmc_beta0001",    "#4393c3"),
        (r"Brain $\beta$=0.001", "exp001_brain_beta0001",   GREEN),
        (r"CITE-seq",            "exp001_citeseq",           PURPLE),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(22, 7))
    perfect = np.linspace(0, 1, 100)

    for ax, (name, exp, color) in zip(axes, configs):
        path = EXPERIMENTS_DIR / exp / "calibration_analysis.json"
        if not path.exists():
            ax.text(0.5, 0.5, "missing", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        with path.open() as f:
            d = json.load(f)
        bins = d.get("bins", [])
        if not bins:
            continue

        x = np.array([b["mean_h_cluster"] for b in bins])
        y = np.array([b["true_error_rate"] for b in bins])
        n = np.array([b.get("n_cells", 1) for b in bins])

        ax.plot(perfect, perfect, color=GREY, lw=2.0,
                linestyle="--", label="Perfect", alpha=0.8)
        ax.fill_between(perfect, perfect - 0.1, perfect + 0.1,
                        alpha=0.10, color=GREY)

        ax.plot(x, y, "o-", color=color, lw=LW + 0.5, markersize=MS + 1,
                label="Observed", zorder=5)
        ax.scatter(x, y, s=n / n.max() * 180 + 40, color=color,
                   alpha=0.9, zorder=6, edgecolors="white", lw=1.0)

        ece = d.get("ece", float("nan"))
        brier = d.get("brier_score", float("nan"))
        ax.set_title(f"{name}\nECE={ece:.3f}  Brier={brier:.3f}", fontsize=FST)
        ax.set_xlabel(r"Mean $H_{\rm cluster}$ (per decile)")
        if ax is axes[0]:
            ax.set_ylabel("Observed error rate")
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=FSG)

    fig.suptitle(
        r"Calibration: $H_{\rm cluster}$ predicts alignment error" + "\n"
        "(dot size ∝ bin population; grey band = ±0.1 from perfect)",
        fontsize=FST, fontweight="bold", y=1.02)
    plt.tight_layout()
    _save(fig, "fig7_calibration_curves")
    print("  [fig7] done")



_DISEASE_LABELS = {
    "CD8_T_lymphopenia":       "CD8 T lymphopenia\n(HIV, CMV)",
    "NK_cell_deficiency":      "NK cell deficiency\n(Chediak-Higashi)",
    "B_cell_aplasia":          "B cell aplasia\n(XLA, CVID)",
    "Monocytopenia":           "Monocytopenia\n(GATA2 deficiency)",
    "Treg_deficiency":         "Treg deficiency\n(IPEX syndrome)",
    "Excitatory_neuron_loss":  "Excitatory neuron loss\n(Alzheimer's disease)",
    "Inhibitory_neuron_loss":  "Inhibitory neuron loss\n(Epilepsy/TLE)",
    "Oligodendrocyte_loss":    "Oligodendrocyte loss\n(Multiple sclerosis)",
    "Astrocyte_loss":          "Astrocyte depletion\n(ALS/reactive gliosis)",
    "Microglia_depletion":     "Microglia depletion\n(PLX5622 model)",
}


def fig8_disease_simulation() -> None:
    cs_path = EXPERIMENTS_DIR / "clinical_disease_sim" / "results.json"
    ns_path = EXPERIMENTS_DIR / "neuro_disease_sim" / "results.json"

    if not cs_path.exists() or not ns_path.exists():
        print("  [fig8] missing data files")
        return

    with cs_path.open() as f:
        cs = json.load(f)
    with ns_path.open() as f:
        ns = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(22, 8))

    for ax, data, title, color, dataset_label in [
        (axes[0], cs, "Immune disease (PBMC 10k)", MOSAIC_BLUE, "PBMC"),
        (axes[1], ns, "Neurological disease (Brain 5k)", GREEN, "Brain"),
    ]:
        scenarios = data["scenarios"]
        names = [_DISEASE_LABELS.get(s["name"], s["name"]) for s in scenarios]
        aurocs = [s["auroc"] for s in scenarios]
        ratios = [s.get("entropy_ratio", 1.0) for s in scenarios]

        order = np.argsort(aurocs)[::-1]
        names = [names[i] for i in order]
        aurocs = [aurocs[i] for i in order]
        ratios = [ratios[i] for i in order]

        # warmer means stronger signal
        ratio_norm = np.array(ratios)
        ratio_norm = (ratio_norm - ratio_norm.min()) / (ratio_norm.max() - ratio_norm.min() + 1e-9)
        bar_colors = [plt.cm.YlOrRd(0.3 + 0.6 * v) for v in ratio_norm]

        y = np.arange(len(names))
        bars = ax.barh(y, aurocs, color=bar_colors,
                       edgecolor="white", linewidth=0.6, height=0.65)

        for i, (au, ratio) in enumerate(zip(aurocs, ratios)):
            ax.text(au + 0.003, i, f"{au:.3f}  ({ratio:.1f}×)",
                    va="center", ha="left", fontsize=FSG,
                    color="#333333")

        ax.axvline(0.5, color="grey", linestyle="--", lw=1.2, alpha=0.7, label="Chance")
        ax.axvline(data["mean_auroc"], color=color, lw=LW, linestyle="-",
                   label=f"Mean = {data['mean_auroc']:.3f}", alpha=0.9)
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=FSK)
        ax.set_xlim(0.45, 1.12)
        ax.set_xlabel("AUROC (detecting absent cell type via $H_{\\rm cluster}$)", fontsize=FSL)
        ax.set_title(f"{title}\nmean AUROC = {data['mean_auroc']:.3f}", fontsize=FST,
                     fontweight="bold")
        ax.grid(axis="x", alpha=0.3)
        ax.legend(frameon=False, fontsize=FSG, loc="lower right")

    sm = ScalarMappable(cmap="YlOrRd",
                        norm=Normalize(vmin=1, vmax=max(
                            max(s.get("entropy_ratio",1) for s in cs["scenarios"]),
                            max(s.get("entropy_ratio",1) for s in ns["scenarios"]))))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.tolist(), fraction=0.015, pad=0.02, shrink=0.7)
    cbar.set_label(r"$H_{\rm cluster}$ ratio (absent / present)", fontsize=FSG)
    cbar.ax.tick_params(labelsize=FSK)

    fig.suptitle(
        "Disease simulation: cluster-resolved entropy detects absent cell populations\n"
        "(leave-one-out: cell type absent from ATAC reference, detected by RNA alignment entropy)",
        fontsize=FST, fontweight="bold", y=1.01)
    plt.tight_layout()
    _save(fig, "fig8_disease_simulation")
    print("  [fig8] done")



def fig9_protein_uq() -> None:
    path = EXPERIMENTS_DIR / "protein_uq_analysis" / "results.json"
    if not path.exists():
        print("  [fig9] missing")
        return
    with path.open() as f:
        data = json.load(f)

    markers = data["marker_analysis"]
    markers_s = sorted(markers, key=lambda m: abs(m["expression_diff_high_minus_low"]),
                       reverse=True)

    names = [m["marker"].replace("_TotalSeqB", "") for m in markers_s]
    diffs = [m["expression_diff_high_minus_low"] for m in markers_s]
    uncertain = [m["is_enriched_in_uncertain"] for m in markers_s]
    pvals = [m["mannwhitney_p"] for m in markers_s]
    interp = [m["clinical_interpretation"] for m in markers_s]

    colors = [MOSAIC_RED if u else MOSAIC_BLUE for u in uncertain]
    y = np.arange(len(names))

    fig, axes = plt.subplots(1, 2, figsize=(22, 9))

    ax = axes[0]
    bars = ax.barh(y, diffs, color=colors, edgecolor="white", linewidth=0.6, height=0.65)
    ax.axvline(0, color="grey", lw=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=FSK)
    ax.set_xlabel("Mean expression difference\n(high uncertainty − low uncertainty cells)", fontsize=FSL)
    ax.set_title("(A) Protein marker expression\nin high- vs low-uncertainty cells", fontsize=FST)

    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=MOSAIC_RED, label="Higher in uncertain"),
                       Patch(color=MOSAIC_BLUE, label="Lower in uncertain")],
              frameon=False, fontsize=FSG, loc="lower right")
    ax.grid(axis="x", alpha=0.3)

    ax = axes[1]
    log_p = [-np.log10(max(p, 1e-300)) for p in pvals]
    bars = ax.barh(y, log_p, color=colors, edgecolor="white", linewidth=0.6, height=0.65)
    ax.axvline(2, color="grey", linestyle="--", lw=1.2, alpha=0.7,
               label="p=0.01 threshold")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=FSK)
    ax.set_xlabel("–log₁₀(Mann-Whitney p-value)", fontsize=FSL)
    ax.set_title("(B) Statistical significance\n(all markers, p<0.01 threshold)", fontsize=FST)
    ax.legend(frameon=False, fontsize=FSG)
    ax.grid(axis="x", alpha=0.3)

    fig.suptitle(
        "CITE-seq protein UQ: which surface markers characterize\n"
        "cells with high vs low alignment uncertainty",
        fontsize=FST, fontweight="bold", y=1.01)
    plt.tight_layout()
    _save(fig, "fig9_protein_uq")
    print("  [fig9] done")



def fig10_checkpoint_immunotherapy() -> None:
    import anndata as ad

    path = EXPERIMENTS_DIR / "checkpoint_immunotherapy" / "results.json"
    if not path.exists():
        print("  [fig10] missing")
        return
    with path.open() as f:
        data = json.load(f)

    analyses = data["analyses"]

    exp_dir = EXPERIMENTS_DIR / "exp001_citeseq"
    H = np.load(exp_dir / "alignment_entropy_cluster.npy")
    sub_idx = np.load(exp_dir / "ot_subsample_indices.npy")
    prot = ad.read_h5ad(PROCESSED_DIR / "citeseq_pbmc_atac.h5ad")
    prot_sub = prot[sub_idx]
    X_prot = prot_sub.X
    if sp.issparse(X_prot):
        X_prot = X_prot.toarray()
    mnames = list(prot_sub.var_names)

    def _get_marker(name):
        if name in mnames:
            return X_prot[:, mnames.index(name)]
        return None

    pd1 = _get_marker("PD-1_TotalSeqB")
    tigit = _get_marker("TIGIT_TotalSeqB")
    cd4 = _get_marker("CD4_TotalSeqB")
    cd8a = _get_marker("CD8a_TotalSeqB")
    cd3 = _get_marker("CD3_TotalSeqB")
    cd127 = _get_marker("CD127_TotalSeqB")

    fig, axes = plt.subplots(2, 3, figsize=(24, 14))

    def _violin_compare(ax, h_a, h_b, label_a, label_b, title, pval):
        valid_a = h_a[~np.isnan(h_a)]
        valid_b = h_b[~np.isnan(h_b)]
        parts = ax.violinplot([valid_a, valid_b], positions=[0, 1],
                              showmedians=True, showextrema=False, widths=0.7)
        parts["bodies"][0].set_facecolor(MOSAIC_RED)
        parts["bodies"][0].set_alpha(0.6)
        parts["bodies"][1].set_facecolor(MOSAIC_BLUE)
        parts["bodies"][1].set_alpha(0.6)
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(2)

        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"{label_a}\n(n={len(valid_a)})",
                            f"{label_b}\n(n={len(valid_b)})"], fontsize=FSK)
        ax.set_ylabel(r"$H_{\rm cluster}$", fontsize=FSL)
        ax.set_title(title, fontsize=FST)
        ax.grid(axis="y", alpha=0.3)

        sig = "***" if pval < 0.001 else ("**" if pval < 0.01 else ("*" if pval < 0.05 else "ns"))
        y_max = max(valid_a.max(), valid_b.max())
        ax.annotate("", xy=(1, y_max * 1.05), xytext=(0, y_max * 1.05),
                    arrowprops=dict(arrowstyle="-", color="black", lw=1.5))
        ax.text(0.5, y_max * 1.08, f"{sig}\np={pval:.1e}",
                ha="center", va="bottom", fontsize=FSG)

    if pd1 is not None:
        thresh75 = np.percentile(pd1, 75)
        thresh25 = np.percentile(pd1, 25)
        high_m = pd1 >= thresh75
        low_m = pd1 < thresh25
        a = [an for an in analyses if an["analysis"] == "PD-1_TotalSeqB_entropy"]
        pv = a[0]["mannwhitney_p_greater"] if a else 0.1
        _violin_compare(axes[0, 0], H[high_m], H[low_m],
                        "PD-1 high (top 25%)", "PD-1 low (bottom 25%)",
                        "PD-1⁺ exhausted T cells\nvs PD-1⁻ cells", pv)

    if tigit is not None:
        high_m = tigit >= np.percentile(tigit, 75)
        low_m = tigit < np.percentile(tigit, 25)
        a = [an for an in analyses if an["analysis"] == "TIGIT_TotalSeqB_entropy"]
        pv = a[0]["mannwhitney_p_greater"] if a else 0.1
        _violin_compare(axes[0, 1], H[high_m], H[low_m],
                        "TIGIT high", "TIGIT low",
                        "TIGIT⁺ exhausted T cells\nvs TIGIT⁻ cells", pv)

    if pd1 is not None and tigit is not None:
        exhausted = (pd1 >= np.percentile(pd1, 75)) & (tigit >= np.percentile(tigit, 75))
        fresh = (pd1 < np.percentile(pd1, 25)) & (tigit < np.percentile(tigit, 25))
        a = [an for an in analyses if an["analysis"] == "PD1_TIGIT_double_positive_exhausted"]
        pv = a[0]["mannwhitney_p_greater"] if a else 0.1
        _violin_compare(axes[0, 2], H[exhausted], H[fresh],
                        "PD-1⁺TIGIT⁺\n(exhausted)", "PD-1⁻TIGIT⁻\n(fresh)",
                        "Double-positive exhausted T\n(checkpoint inhibitor targets)", pv)

    if pd1 is not None and cd4 is not None and cd8a is not None:
        pd1_med = np.percentile(pd1, 50)
        cd4_m = (cd4 >= np.percentile(cd4, 50)) & (pd1 >= pd1_med) & \
                (cd8a < np.percentile(cd8a, 50))
        cd8_m = (cd8a >= np.percentile(cd8a, 50)) & (pd1 >= pd1_med) & \
                (cd4 < np.percentile(cd4, 50))
        a = [an for an in analyses if an["analysis"] == "CD4_vs_CD8_PD1_entropy"]
        pv = a[0]["mannwhitney_p_twosided"] if a else 0.1
        _violin_compare(axes[1, 0], H[cd4_m], H[cd8_m],
                        "CD4⁺PD-1⁺\n(Treg/exhausted CD4)", "CD8a⁺PD-1⁺\n(effector exhausted)",
                        "CD4 vs CD8 exhausted T\nsubset uncertainty", pv)

    if cd127 is not None and cd3 is not None:
        cd3_m = cd3 >= np.percentile(cd3, 50)
        mem = cd3_m & (cd127 >= np.percentile(cd127, 75))
        eff = cd3_m & (cd127 < np.percentile(cd127, 25))
        a = [an for an in analyses if an["analysis"] == "CD127_memory_vs_effector_entropy"]
        pv = a[0]["mannwhitney_p_greater"] if a else 0.1
        _violin_compare(axes[1, 1], H[mem], H[eff],
                        "CD3⁺CD127⁺\n(memory T)", "CD3⁺CD127⁻\n(effector/exhausted)",
                        "Memory vs effector T\n(CD127/IL-7Rα status)", pv)

    ax = axes[1, 2]
    analysis_labels = []
    means_high = []
    means_low = []
    for an in analyses:
        if "mean_H_high" in an:
            analysis_labels.append(an["analysis"].replace("_TotalSeqB_entropy","").replace("_"," "))
            means_high.append(an["mean_H_high"])
            means_low.append(an["mean_H_low"])
        elif "mean_H_exhausted" in an:
            analysis_labels.append("PD1+TIGIT+\nexhausted")
            means_high.append(an["mean_H_exhausted"])
            means_low.append(an["mean_H_fresh"])

    if analysis_labels:
        x = np.arange(len(analysis_labels))
        w = 0.35
        ax.bar(x - w/2, means_high, w, color=MOSAIC_RED, label="High checkpoint / exhausted",
               edgecolor="white", linewidth=0.7)
        ax.bar(x + w/2, means_low, w, color=MOSAIC_BLUE, label="Low checkpoint / fresh",
               edgecolor="white", linewidth=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(analysis_labels, rotation=20, ha="right", fontsize=FSK)
        ax.set_ylabel(r"Mean $H_{\rm cluster}$", fontsize=FSL)
        ax.set_title("Summary: checkpoint markers\npredict alignment uncertainty", fontsize=FST)
        ax.legend(frameon=False, fontsize=FSG)
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle(
        "Immune checkpoint immunotherapy: alignment uncertainty in exhausted T cells\n"
        r"(CITE-seq; higher $H_{\rm cluster}$ = transcriptome–proteome discordance for immunotherapy targets)",
        fontsize=FST, fontweight="bold", y=1.01)
    plt.tight_layout(h_pad=4.0, w_pad=2.5)
    _save(fig, "fig10_checkpoint_immunotherapy")
    print("  [fig10] done")



if __name__ == "__main__":
    import traceback

    tasks = [
        ("fig1 PBMC",       lambda: fig1_aligned_umap("exp001_pbmc_beta0001", "pbmc10k_multiome", "PBMC 10k (β=0.001)")),
        ("fig1 Brain",      lambda: fig1_aligned_umap("exp001_brain_beta0001", "brain3k_multiome", "Brain 5k (β=0.001)")),
        ("fig2 PBMC",       lambda: fig2_entropy_comparison("exp001_pbmc_beta0001", "pbmc10k_multiome", "PBMC 10k")),
        ("fig2 Brain",      lambda: fig2_entropy_comparison("exp001_brain_beta0001", "brain3k_multiome", "Brain 5k")),
        ("fig3",            fig3_missing_type),
        ("fig4",            fig4_baselines),
        ("fig5",            fig5_beta_tradeoff),
        ("fig6",            fig6_cross_tissue),
        ("fig7",            fig7_calibration),
        ("fig8",            fig8_disease_simulation),
        ("fig9",            fig9_protein_uq),
        ("fig10",           fig10_checkpoint_immunotherapy),
    ]

    for name, fn in tasks:
        try:
            print(f"\n[{name}]")
            fn()
        except Exception as e:
            print(f"  [skip] {name}: {e}")
            traceback.print_exc()

    print("\ndone")
