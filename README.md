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
├── CITATION.cff
├── requirements.txt
├── src/
│   ├── models/               # IB-VAE encoder/decoder definitions
│   ├── data/                 # Preprocessing and data loading
│   ├── training/             # Training loops and logging
│   ├── evaluation/           # FOSCTTM, label-transfer, ARI, cluster entropy
│   └── utils/                # Paths, Sinkhorn solver, helpers
├── scripts/                  # Experiment runners and figure generation
├── experiments/              # JSON results for every experiment
├── figures/                  # Publication-quality figures (PNG + PDF)
└── mlhc_submission/          # MLHC 2026 paper (LaTeX + compiled PDF)
    ├── main.tex
    ├── main.pdf
    ├── references.bib
    └── figures/
```

---

## Reproducing the results

All experiments were run on a single NVIDIA RTX 3080 (16 GB). Wall time: ~120 s (PBMC), ~50 s (brain), ~95 s (CITE-seq) per alignment run.

```bash
# 1. Install dependencies (Python 3.11+ recommended)
pip install -r requirements.txt

# 2. Run the primary alignment experiment for a dataset
python scripts/run_experiment.py --dataset pbmc10k_multiome --beta 0.001

# 3. Multi-seed replication (seeds 0-9 for PBMC; 0-2 for Brain/CITE-seq)
python scripts/run_experiment.py --dataset brain3k_multiome --beta 0.001 --seed 0
python scripts/aggregate_seeds.py --dataset brain3k_multiome --beta 0.001

# 4. Downstream experiments
python scripts/missing_type_exp.py           # leave-one-cluster-out
python scripts/cross_tissue_exp.py           # PBMC RNA x Brain ATAC control
python scripts/clinical_disease_sim.py       # 5 immunodeficiency scenarios
python scripts/neuro_disease_sim.py          # 5 CNS disease scenarios
python scripts/checkpoint_immunotherapy_analysis.py
python scripts/protein_uq_analysis.py
python scripts/rare_cell_detection.py

# 5. Regenerate every figure in the paper
python scripts/generate_all_figures.py
```

### Mapping paper claims to artifacts

| Paper claim | Aggregate JSON |
|---|---|
| PBMC 10-seed alignment metrics | `experiments/aggregate_pbmc10k_multiome_beta0001_10seed.json` |
| Brain 3-seed alignment metrics | `experiments/aggregate_brain3k_multiome_beta0.001.json` |
| CITE-seq 3-seed alignment metrics | `experiments/aggregate_citeseq_3seed.json` |
| Missing cell-type detection AUROC | `experiments/*/exp003_missing_type.json` |
| Cross-tissue negative control | `experiments/exp006_cross_tissue/results.json` |
| Calibration curves (ECE, Brier) | `experiments/*/calibration_analysis.json` |
| Immunodeficiency simulation | `experiments/clinical_disease_sim/results.json` |
| Neurological disease simulation | `experiments/neuro_disease_sim/results.json` |
| Checkpoint-immunotherapy entropy | `experiments/checkpoint_immunotherapy/results.json` |
| Protein marker UQ analysis | `experiments/protein_uq_analysis/results.json` |
| Rare-cell sanity check | `experiments/rare_cell_detection/results.json` |

All random seeds are set deterministically (NumPy, PyTorch, cuDNN). Baselines (SCOT, uniPort) use a separate virtual environment pinned to `numpy<2`; see `scripts/run_uniport_venv.py`.

---

## Data

All three benchmarks are publicly released demonstration datasets from 10x Genomics with no identifying patient information. Data paths are configured in `src/utils/paths.py`; raw data is expected under `data/raw/`.

| Dataset | Modalities | Cells (paired) | Clusters | Source |
|---|---|---|---|---|
| PBMC 10k | RNA + ATAC | 11,303 | 18 | 10x Genomics Multiome demo |
| E18 mouse brain 5k | RNA + ATAC | 4,531 | 20 | 10x Genomics Multiome demo |
| PBMC 10k CITE-seq | RNA + 14 ADT proteins | ~10,000 | 16 | 10x Genomics CITE-seq demo |

---

## License

MIT License — Copyright (c) 2026 Bryan Cheng, Austin Jin, Joshua Chang, Brendan Lo
