# Calibrated Per-Cell Uncertainty for Unpaired Single-Cell Multi-Omics Integration

**Authors**: Bryan Cheng, Austin Jin, Joshua Chang, Brendan Lo  
**License**: MIT  
**Submitted to**: Machine Learning for Healthcare (MLHC) 2026

---

## What this is

Single-cell disease atlases are built by measuring different molecular properties of cells — gene expression (RNA), chromatin accessibility (ATAC), and surface proteins — but almost never on the same individual cell. Existing methods that computationally match cells across these modalities return a single best-guess alignment with no indication of confidence. A wrong match propagates silently into every downstream analysis: cell-type discovery, regulatory inference, and biomarker identification.

This project introduces a method that treats unpaired multi-omics alignment as a probabilistic problem and outputs a **per-cell uncertainty score** alongside each match. The key technical contribution is *cluster-resolved alignment entropy* — a calibrated signal obtained by marginalizing the optimal-transport plan over partner cluster labels rather than individual cells, which corrects a fundamental directional miscalibration in naive row-entropy scores.

---

## Key results

| Dataset | FOSCTTM ↓ | Label Transfer ↑ | ARI ↑ | Missing-type AUROC ↑ |
|---|---|---|---|---|
| PBMC immune (11,303 cells) | 0.194 | 0.664 / 0.694 | 0.651 | 0.960 |
| Mouse brain E18 (4,531 cells) | 0.052 | 0.960 / 0.967 | 0.932 | 0.995 |
| PBMC CITE-seq (10,000 cells) | 0.094 | 0.872 / 0.960 | 0.791 | 0.946 |

- Wrong-cluster detection AUROC: **0.88–0.95** across all datasets
- Cross-tissue negative control: **4.2× higher** entropy than valid within-dataset alignments
- PD-1⁺TIGIT⁺ exhausted T cells (checkpoint inhibitor targets) have significantly higher uncertainty than fresh T cells (*p* = 9.2 × 10⁻⁵)

---

## Method overview

The pipeline has three stages:

1. **Per-modality encoding** — Each modality (RNA or ATAC) is independently compressed to a 64-dimensional latent representation using an information-bottleneck variational autoencoder. A cluster-centroid cross-modal prediction head forces both representations into a compatible space without per-cell memorization.

2. **Procrustes alignment** — An orthogonal rotation is fit in closed form (SVD) to align the two latent coordinate systems using cluster-centroid landmarks.

3. **Entropic optimal transport + uncertainty** — A Sinkhorn transport plan is computed between the aligned latents. For each cell, the transport mass is marginalized over partner *cluster labels* (not individual cells) and the entropy of that cluster distribution is the uncertainty score $H_\text{cluster} \in [0, 1]$.

---

## Repository structure

```
.
├── LICENSE
├── README.md
├── RESEARCH_PLAN.md          # Original research plan
├── RESULTS.md                # Comprehensive experimental results
├── TRAINING_LOG.md           # Full log of all training runs
├── REVIEW_REPORT_PRE_TRAINING.md
├── REVIEW_REPORT_POST_RESULTS.md
├── src/
│   ├── models/               # IB-VAE encoder/decoder definitions
│   ├── data/                 # Preprocessing and data loading
│   ├── training/             # Training loops and logging
│   ├── evaluation/           # Alignment metrics (FOSCTTM, LT, ARI)
│   └── utils/                # Paths, helpers, Sinkhorn solver
├── scripts/                  # Experiment runners and figure generation
│   ├── run_experiment.py
│   ├── generate_all_figures.py
│   ├── rare_cell_detection.py
│   ├── checkpoint_immunotherapy_analysis.py
│   └── ...
├── experiments/              # JSON results for every experiment run
├── figures/                  # Publication-quality figures (PNG + PDF)
├── mlhc_submission/          # MLHC 2026 paper (LaTeX + compiled PDF)
│   ├── main.tex
│   ├── main.pdf
│   └── figures/
└── paper/                    # Draft paper sections (Markdown)
```

---

## Reproducing the results

All experiments were run on a single NVIDIA RTX 3080 (16 GB). Wall time: ~120 s (PBMC), ~50 s (brain), ~95 s (CITE-seq) per run.

```bash
# Install dependencies
pip install -r requirements.txt

# Run primary alignment experiment
python scripts/run_experiment.py --dataset pbmc10k_multiome --beta 0.001

# Regenerate all figures
python scripts/generate_all_figures.py
```

Aggregate results across seeds are in `experiments/aggregate_*.json`. All random seeds are set deterministically.

---

## License

MIT License — Copyright (c) 2026 Bryan Cheng, Austin Jin, Joshua Chang, Brendan Lo
