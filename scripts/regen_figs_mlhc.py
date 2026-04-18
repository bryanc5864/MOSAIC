"""Regenerate fig 3 (missing type), fig 4 (baselines), fig 5 (beta), fig 7 (calibration) with larger fonts and blue/green palette for MLHC submission."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json
import numpy as np
from pathlib import Path

plt.rcParams.update({
    "font.size": 20,
    "axes.labelsize": 22,
    "axes.titlesize": 22,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.2,
    "lines.linewidth": 2.0,
})

BLUE = "#2c7fb8"
GREEN = "#41ae76"
BLUE_DARK = "#08519c"
GREEN_DARK = "#006d2c"

REPO = Path(__file__).resolve().parent.parent
E = REPO / "experiments"
F = REPO / "figures"


def fig3():
    fig, axes = plt.subplots(1, 2, figsize=(20, 8), sharey=True)
    for ax, (name, exp) in zip(axes, [
        ("PBMC 10k ($\\beta$=0.01)", "exp001_pbmc_final"),
        ("Brain 5k ($\\beta$=0.001)", "exp001_brain_beta0001"),
    ]):
        with (E / exp / "exp003_missing_type.json").open() as f:
            data = json.load(f)
        per = data["per_cluster"]
        cl = [p["target_cluster"] for p in per if p.get("auroc_cluster_entropy") is not None]
        au = [p["auroc_cluster_entropy"] for p in per if p.get("auroc_cluster_entropy") is not None]
        o = np.argsort(au)[::-1]
        cl_s = [cl[i] for i in o]
        au_s = [au[i] for i in o]
        bar_color = BLUE if "PBMC" in name else GREEN
        ax.bar(range(len(au_s)), au_s, color=bar_color, edgecolor="black", linewidth=0.8)
        ax.axhline(0.5, color="gray", linestyle="--", linewidth=1.5, label="random")
        ax.set_xticks(range(len(cl_s)))
        ax.set_xticklabels([f"c{c}" for c in cl_s], rotation=60, ha="right", fontsize=16)
        ax.set_ylim(0.5, 1.02)
        ax.set_title(f"{name}\nmean AUROC = {data['mean_auroc']:.3f}")
        ax.set_ylabel("AUROC")
        ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    fig.savefig(F / "fig3_missing_type_auroc.png")
    plt.close(fig)
    print("[fig3] saved")


def fig5():
    with (E / "aggregate_brain3k_multiome_beta0.01.json").open() as f:
        brain_b01 = json.load(f)
    with (E / "aggregate_brain3k_multiome_beta0.001.json").open() as f:
        brain_b0001 = json.load(f)
    with (E / "aggregate_pbmc10k_multiome.json").open() as f:
        pbmc_b01 = json.load(f)
    with (E / "aggregate_pbmc10k_multiome_beta0.001.json").open() as f:
        pbmc_b0001 = json.load(f)

    metric_keys = [
        ("FOSCTTM\n(lower better)", "foscttm_mean"),
        ("LT RNA$\\to$ATAC", "lt_rna_to_atac"),
        ("LT ATAC$\\to$RNA", "lt_atac_to_rna"),
        ("Joint ARI", "joint_ari"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(22, 7))
    for ax, (title, key) in zip(axes, metric_keys):
        p01 = (pbmc_b01[f"{key}_mean"], pbmc_b01[f"{key}_std"])
        p0001 = (pbmc_b0001[f"{key}_mean"], pbmc_b0001[f"{key}_std"])
        b01 = (brain_b01[f"{key}_mean"], brain_b01[f"{key}_std"])
        b0001 = (brain_b0001[f"{key}_mean"], brain_b0001[f"{key}_std"])
        x = np.arange(2)
        w = 0.35
        ax.bar(x - w / 2, [p01[0], p0001[0]], w, yerr=[p01[1], p0001[1]],
               capsize=6, label="PBMC 10k", color=BLUE, edgecolor="black", linewidth=0.8,
               error_kw={"ecolor": "black", "elinewidth": 1.2})
        ax.bar(x + w / 2, [b01[0], b0001[0]], w, yerr=[b01[1], b0001[1]],
               capsize=6, label="Brain 5k", color=GREEN, edgecolor="black", linewidth=0.8,
               error_kw={"ecolor": "black", "elinewidth": 1.2})
        ax.set_xticks(x)
        ax.set_xticklabels(["$\\beta$=0.01", "$\\beta$=0.001"])
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
        if ax is axes[0]:
            ax.legend(frameon=False, loc="best")
    plt.tight_layout()
    fig.savefig(F / "fig5_beta_tradeoff.png")
    plt.close(fig)
    print("[fig5] saved")


def fig4():
    def load_baselines(ds):
        base = E / f"baselines_{ds}"
        with (base / "baseline_results.json").open() as f:
            scot_raw = json.load(f)["scot"]["metrics"]
        with (base / "simple_baseline_results.json").open() as f:
            simple = json.load(f)
        def flat(m):
            out = {}
            if "foscttm" in m and isinstance(m["foscttm"], dict):
                out["foscttm_mean"] = m["foscttm"].get("foscttm_mean")
            else:
                out["foscttm_mean"] = m.get("foscttm_mean")
            out["lt_rna_to_atac"] = m.get("label_transfer_rna_to_atac", m.get("lt_rna_to_atac"))
            out["lt_atac_to_rna"] = m.get("label_transfer_atac_to_rna", m.get("lt_atac_to_rna"))
            out["joint_ari"] = m.get("joint_clustering_ari", m.get("joint_ari"))
            return out
        return flat(scot_raw), flat(simple["nn_on_ib"].get("metrics", {})), flat(simple["raw_ot"].get("metrics", {}))

    def load_aggr(path):
        with path.open() as f:
            return json.load(f)

    pbmc_scot, pbmc_nn, pbmc_raw = load_baselines("pbmc10k_multiome")
    brain_scot, brain_nn, brain_raw = load_baselines("brain3k_multiome")
    pbmc_m = load_aggr(E / "aggregate_pbmc10k_multiome_beta0001_10seed.json")
    brain_m = load_aggr(E / "aggregate_brain3k_multiome_beta0.001.json")

    metrics = [("FOSCTTM$\\downarrow$", "foscttm_mean"),
               ("LT RNA$\\to$ATAC", "lt_rna_to_atac"),
               ("LT ATAC$\\to$RNA", "lt_atac_to_rna"),
               ("ARI", "joint_ari")]

    def _flt(v):
        if v is None or isinstance(v, dict):
            return float("nan")
        try:
            return float(v)
        except Exception:
            return float("nan")

    labels = ["MOSAIC", "NN-on-IB", "SCOT", "Raw"]
    colors = [BLUE, GREEN, BLUE_DARK, GREEN_DARK]

    fig, axes = plt.subplots(2, 4, figsize=(22, 11))
    for row_i, (ds_name, mosaic, scot, nn, raw) in enumerate([
        ("PBMC 10k", pbmc_m, pbmc_scot, pbmc_nn, pbmc_raw),
        ("Brain 5k", brain_m, brain_scot, brain_nn, brain_raw),
    ]):
        for col_i, (title, key) in enumerate(metrics):
            ax = axes[row_i, col_i]
            vals = [
                _flt(mosaic.get(f"{key}_mean")),
                _flt(nn.get(key)),
                _flt(scot.get(key)),
                _flt(raw.get(key)),
            ]
            errs = [_flt(mosaic.get(f"{key}_std", 0)) or 0.0, 0.0, 0.0, 0.0]
            ax.bar(range(4), vals, yerr=errs, color=colors,
                   edgecolor="black", linewidth=0.8, capsize=6)
            ax.set_xticks(range(4))
            ax.set_xticklabels(labels, rotation=30, ha="right")
            if row_i == 0:
                ax.set_title(title)
            ax.grid(axis="y", alpha=0.3)
            if col_i == 0:
                ax.set_ylabel(f"{ds_name}\n{title.split(chr(36))[0].strip()}" if row_i == 1 else ds_name,
                              fontsize=20)
    plt.tight_layout()
    fig.savefig(F / "fig4_baselines_both.png")
    plt.close(fig)
    print("[fig4] saved")


def fig7():
    configs = [
        ("PBMC $\\beta$=0.01", "exp001_pbmc_final"),
        ("PBMC $\\beta$=0.001", "exp001_pbmc_beta0001"),
        ("Brain $\\beta$=0.001", "exp001_brain_beta0001"),
        ("CITE-seq $\\beta$=0.001", "exp001_citeseq"),
    ]
    fig, axes = plt.subplots(1, len(configs), figsize=(24, 7))
    for i, (ax, (name, exp)) in enumerate(zip(axes, configs)):
        path = E / exp / "calibration_analysis.json"
        if not path.exists():
            ax.text(0.5, 0.5, "missing", ha="center", va="center")
            continue
        with path.open() as f:
            d = json.load(f)
        bins = d.get("bins", [])
        if not bins:
            ax.text(0.5, 0.5, "no bins", ha="center", va="center")
            continue
        x = [b["mean_h_cluster"] for b in bins]
        y = [b["true_error_rate"] for b in bins]
        perfect_line, = ax.plot([0, 1], [0, 1], color=GREEN, linestyle="--",
                                linewidth=2.0, label="perfect")
        obs_line, = ax.plot(x, y, marker="o", color=BLUE, linewidth=2.5,
                            markersize=10, label="observed")
        ax.set_xlabel("mean $H_{\\mathrm{cluster}}$ per decile")
        if i == 0:
            ax.set_ylabel("observed error rate")
        ax.set_title(f"{name}\nECE={d.get('ece', np.nan):.3f}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
    # One shared legend anchored OUTSIDE the axes at the bottom
    fig.legend(handles=[perfect_line, obs_line], labels=["perfect", "observed"],
               loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=2, frameon=False, fontsize=20)
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(F / "fig7_calibration_curves.png")
    plt.close(fig)
    print("[fig7] saved")


if __name__ == "__main__":
    fig3()
    fig5()
    try:
        fig4()
    except Exception as e:
        print(f"[fig4] failed: {e}")
    try:
        fig7()
    except Exception as e:
        print(f"[fig7] failed: {e}")
