"""Generate fig_disease_simulation.png: AUROC per scenario for
clinical immunodeficiency simulation (PBMC, Exp 13) and neurological
disease simulation (Brain 5k, Exp 15). Blue bars for PBMC immune,
green bars for Brain neuro."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

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
})

BLUE = "#2c7fb8"
GREEN = "#41ae76"

REPO = Path(__file__).resolve().parent.parent
E = REPO / "experiments"
F = REPO / "figures"


def load(name):
    with (E / name / "results.json").open() as f:
        return json.load(f)


def short(scenario_name):
    m = {
        "CD8_T_lymphopenia": "CD8 T lymphopenia",
        "NK_cell_deficiency": "NK deficiency",
        "B_cell_aplasia": "B cell aplasia",
        "Monocytopenia": "Monocytopenia",
        "Treg_deficiency": "Treg deficiency",
        "Excitatory_neuron_loss": "Excitatory neuron (AD)",
        "Inhibitory_neuron_loss": "Inhibitory neuron (epilepsy)",
        "Oligodendrocyte_loss": "Oligodendrocyte (MS)",
        "Astrocyte_loss": "Astrocyte (ALS)",
        "Microglia_depletion": "Microglia (PLX5622)",
    }
    return m.get(scenario_name, scenario_name.replace("_", " "))


def main():
    immune = load("clinical_disease_sim")["scenarios"]
    neuro = load("neuro_disease_sim")["scenarios"]
    immune_sorted = sorted(immune, key=lambda s: -s["auroc"])
    neuro_sorted = sorted(neuro, key=lambda s: -s["auroc"])

    fig, axes = plt.subplots(1, 2, figsize=(20, 7.5), sharey=True)

    for ax, scenarios, color, label, mean in [
        (axes[0], immune_sorted, BLUE, "PBMC 10k — immunodeficiency",
         np.mean([s["auroc"] for s in immune_sorted])),
        (axes[1], neuro_sorted, GREEN, "Brain 5k — neurological disease",
         np.mean([s["auroc"] for s in neuro_sorted])),
    ]:
        names = [short(s["name"]) for s in scenarios]
        vals = [s["auroc"] for s in scenarios]
        ax.bar(range(len(vals)), vals, color=color, edgecolor="black", linewidth=0.8)
        ax.axhline(0.5, color="gray", linestyle="--", linewidth=1.5, label="random")
        ax.axhline(mean, color="black", linestyle=":", linewidth=1.8,
                   label=f"mean = {mean:.3f}")
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=40, ha="right", fontsize=16)
        ax.set_ylim(0.5, 1.02)
        ax.set_title(label)
        ax.set_ylabel("AUROC")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(frameon=False, loc="lower left")

    plt.tight_layout()
    fig.savefig(F / "fig_disease_simulation.png")
    plt.close(fig)
    print("saved fig_disease_simulation.png")


if __name__ == "__main__":
    main()
