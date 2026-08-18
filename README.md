# Calibrated per-cell uncertainty for unpaired single-cell multi-omics integration

Single-cell atlases measure RNA, chromatin accessibility, and surface proteins on
overlapping but non-identical cells. Methods that computationally match cells across
those modalities return one best-guess alignment and say nothing about how much to
trust any individual match, so a wrong match propagates silently into cell-type calls,
regulatory inference, and biomarker discovery.

This repository treats the alignment as a probabilistic problem and attaches a per-cell
uncertainty score to every match. The score is the entropy of the entropic-OT transport
plan after marginalizing each row over *partner cluster labels* rather than over
individual partner cells. That marginalization matters: raw row entropy is
anti-correlated with true alignment error (Spearman -0.65 on PBMC), because cells in
dense clusters have both many plausible partners and a very close true partner. Taking
the cluster marginal first removes the confound and flips the sign.

The pipeline is three stages. Each modality is compressed to 64 dimensions by its own
information-bottleneck VAE with a cluster-centroid cross-modal prediction head. An
orthogonal Procrustes rotation, fit in closed form on cluster centroids, brings the two
latent spaces into a common frame. Sinkhorn OT then produces the transport plan, and the
cluster-marginal entropy of each row is the uncertainty score.

## Install

Python 3.11 and a single GPU. Runs take roughly two minutes per alignment.

```bash
pip install -r requirements.txt
```

SCOT and uniPort need a separate environment pinned to `numpy<2`, because uniPort 1.3
still calls `np.Inf`. See `scripts/run_uniport_venv.py`.

## Run

Everything imports from `src/`, so run from the repository root with the root on the
path.

```bash
# download and preprocess one benchmark
python -m src.data.datasets --which pbmc10k_multiome
python -m src.data.preprocess pbmc10k_multiome
python -m src.data.validate pbmc10k_multiome

# train both encoders, align, and score
python -m src.training.run_experiment \
    --dataset pbmc10k_multiome --name exp001_pbmc_beta0001 --beta 0.001

# per-cell cluster entropy for that run
PYTHONPATH=. python scripts/cluster_resolved_entropy.py \
    --exp exp001_pbmc_beta0001 --dataset pbmc10k_multiome

# multi-seed summary
PYTHONPATH=. python scripts/aggregate_seeds.py --dataset pbmc10k_multiome \
    --runs exp001_pbmc_beta0001:0 exp001_pbmc_beta0001_seed1:1 exp001_pbmc_beta0001_seed2:2
```

Downstream experiments follow the same pattern: `missing_type_exp.py` (leave one cluster
out), `cross_tissue_exp.py` (PBMC RNA against brain ATAC), `clinical_disease_sim.py` and
`neuro_disease_sim.py` (population knockouts), `protein_uq_analysis.py`,
`checkpoint_immunotherapy_analysis.py`, `rare_cell_detection.py`, `ablation_sweep.py`.
`generate_all_figures.py` rebuilds every paper figure from the JSON in `experiments/`.

## Results

Multi-seed alignment quality, mean plus or minus standard deviation. These come from
`experiments/aggregate_*.json`.

| dataset | config | FOSCTTM (lower better) | label transfer A→B | B→A | ARI |
|---|---|---|---|---|---|
| PBMC 10k multiome | β=0.001, 10 seeds | 0.118 ± 0.008 | 0.913 ± 0.074 | 0.909 ± 0.035 | 0.724 ± 0.103 |
| PBMC 10k multiome | β=0.01, 3 seeds | 0.188 ± 0.006 | 0.689 ± 0.005 | 0.771 ± 0.029 | 0.687 ± 0.009 |
| E18 brain 5k multiome | β=0.001, 3 seeds | 0.049 ± 0.002 | 0.962 ± 0.002 | 0.966 ± 0.010 | 0.935 ± 0.003 |
| PBMC CITE-seq | β=0.001, 3 seeds | 0.091 ± 0.003 | 0.874 ± 0.006 | 0.951 ± 0.009 | 0.766 ± 0.018 |

SCOT, the closest baseline, gets FOSCTTM 0.248 on PBMC and 0.475 on brain; its
Gromov-Wasserstein solver does not converge on the brain or protein data. uniPort comes
out worse than random (FOSCTTM above 0.5) on all three, which we attribute to its
diagonal-integration assumption of a shared feature space.

For the uncertainty signal itself: cluster entropy separates wrong from right cluster
assignments with AUROC 0.883 ± 0.010 (PBMC β=0.01), 0.896 ± 0.026 (PBMC β=0.001), and
0.946 ± 0.017 (brain). In leave-one-cluster-out, where a whole population is removed
from the target modality, mean detection AUROC is 0.960 on PBMC, 0.995 on brain, and
0.946 on CITE-seq; the one bad case on CITE-seq (AUROC 0.32) is a cluster nearly
redundant with another one still present. Aligning PBMC RNA against brain ATAC, where
nothing should match, raises mean cluster entropy to 0.635 against 0.153 within PBMC,
about a fourfold increase. Population-knockout simulations average AUROC 0.952 across
five immunodeficiency scenarios and 0.993 across five neurological ones.

The score discriminates well but is not calibrated as a probability: expected
calibration error runs 0.050 to 0.160 depending on configuration. Platt scaling or
isotonic regression on held-out data is the right fix if you need real probabilities.

## Caveats worth reading before you trust a number

The cross-modal prediction target is a cluster centroid, and the cluster labels on the
ATAC side are copied from the paired RNA side (`src/data/preprocess.py`,
`add_cross_modal_targets`). Cell-level pairing never reaches the training loop, and the
evaluation metrics are unaffected, but cluster-level paired information does enter
training. On a genuinely unpaired dataset you would have to cluster one modality alone
or supply external labels, and the numbers above should not be read as a fully unpaired
result.

The right-hand violin panel of `figures/fig6_cross_tissue_negative_control.png` is drawn
from a Beta distribution, not from measured data: `generate_all_figures.py` falls back to
a placeholder shape when the per-cell cross-tissue entropy array is absent, which it
always is because the `.npy` outputs are gitignored. The fourfold ratio in that figure's
title is real and traces to `experiments/exp006_cross_tissue/results.json`; the
distribution shape and the `p<10^-50` annotation are not. Rerun `cross_tissue_exp.py`
and save the entropy array if you need the real distribution.

In the checkpoint-immunotherapy analysis, only the PD-1+TIGIT+ double-positive versus
double-negative contrast is significant (p = 9.2e-5, entropy ratio 1.19). PD-1 alone
(p = 0.30) and TIGIT alone (p = 1.0) are not.

PBMC ARI at β=0.001 has a standard deviation of 0.10 across ten seeds, driven by KMeans
initialization on a more continuous latent. Treat that column as indicative.

## Data

Three public 10x Genomics demonstration datasets, no patient-identifying information.
URLs live in `src/data/datasets.py`; downloads land in `data/raw/` and processed
AnnData in `data/processed/`, both gitignored.

| dataset | modalities | cells | clusters |
|---|---|---|---|
| PBMC 10k multiome | RNA + ATAC | 11,303 | 18 |
| E18 mouse brain 5k multiome | RNA + ATAC | 4,531 | 20 |
| PBMC 10k CITE-seq | RNA + 14 ADT proteins | ~7,900 | 16 |

## Paper

`mlhc_submission/` holds the MLHC 2026 submission (LaTeX plus compiled PDF);
`icml_workshop/` holds a four-page short version. Both are anonymized for review.

MIT licensed.
