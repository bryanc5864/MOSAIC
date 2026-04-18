# Methods

*Calibrated Per-Cell Uncertainty for Unpaired Single-Cell Multi-Omics Integration*

## Overview

the method takes two single-cell modalities $A$ (RNA) and $B$ (ATAC) measured on disjoint cell populations (or paired cells whose correspondence is treated as unknown) and produces (i) a joint latent space in which paired cells map close together and (ii) a per-cell **cluster-resolved alignment uncertainty** score that is calibrated against true alignment error. The pipeline has three stages:

1. Per-modality **Information Bottleneck VAE** (IB-VAE) with a cross-modal prediction head, trained independently on each modality but with a shared cluster-centroid target that encourages cross-modal latent agreement.
2. **Orthogonal Procrustes** post-alignment of the ATAC latent onto the RNA latent frame using shared cluster centroids.
3. **Entropic Sinkhorn** optimal transport between the two aligned latent clouds, producing a transport plan whose row marginals over the partner modality's cluster labels yield the cluster-resolved entropy.

## Data preprocessing

For each benchmark dataset (10x Multiome PBMC 10k; 10x Multiome mouse brain 5k) we apply a standard scanpy pipeline. RNA counts are filtered to cells with ≥200 genes and ≥500 UMIs, normalized to 1e4 total per cell, log1p-transformed, scaled per-gene, and reduced to the top-2000 highly variable genes (scanpy `flavor="seurat_v3"`) on the raw counts. ATAC peak counts are TF-IDF normalized, log1p-transformed, and reduced via 50-component LSI (with the first LSI component dropped, standard practice). Each modality receives an independent Leiden clustering at resolution 1.0 on its own k=15 neighborhood graph. For the paired benchmark datasets, the ground-truth 1:1 pairing is stored in `obs["pair_idx"]` but not exposed to the training loop (see "Transparency note on cluster propagation" below).

## IB-VAE architecture and loss

Each modality's encoder is a 3-layer MLP $512 \to 256 \to 64$ (2 × 64 output for $\mu, \log \sigma$ of the latent posterior). The decoder mirrors the encoder: $64 \to 256 \to 512 \to d_x$. For RNA the decoder produces the three parameters of a zero-inflated negative binomial (ZINB) distribution over raw counts with cell-library-size scaling; the reconstruction loss is the ZINB negative log-likelihood. For ATAC the decoder produces per-peak Bernoulli logits against binarized peak accessibility; the reconstruction loss is binary cross-entropy.

A **cross-modal prediction head** $g_\phi: z \to \mathbb{R}^{d_y}$ is attached to the latent. Its target $y_{\text{cross}}(i)$ for cell $i$ is the centroid, in the partner modality's 50-dim LSI/PCA space, of all cells that share $i$'s cluster assignment. The full IB-VAE loss for modality $m$ is

$$\mathcal{L}_m = \mathcal{L}_{\text{recon},m} + \beta \cdot \mathrm{KL}(q(z|x) \,\|\, \mathcal{N}(0, I)) + \lambda_{\text{pred}} \cdot \| g_\phi(z) - y_{\text{cross}} \|_2^2$$

where $\beta$ is the information-bottleneck weight annealed over 30 epochs from 0 to its target value and $\lambda_{\text{pred}}$ is the cross-modal prediction weight. Training uses AdamW at lr $10^{-3}$ with a 10-epoch cosine warmup and 200 total epochs. We disable early stopping (patience = 999) because the cross-modal prediction loss behaves non-monotonically during the first half of training (see §Training log, Run 006). We use $\beta = 0.01$ and $\lambda_{\text{pred}} = 5$ as the default; the Brain 5k dataset is also reported at $\beta = 0.001$, which is strictly better on every Brain metric.

## Procrustes post-alignment

After independent IB-VAE training on each modality, the two latent clouds live in different regions of the same 64-dim space and must be brought into a common frame. We fit an **orthogonal Procrustes** rotation $R \in \mathrm{O}(64)$ that minimizes $\sum_k \| R \, \bar z^{\text{ATAC}}_k - \bar z^{\text{RNA}}_k \|_2^2$ where $\bar z^{\text{RNA}}_k, \bar z^{\text{ATAC}}_k$ are the RNA and ATAC latent centroids of cluster $k$. The rotation is applied to the ATAC latent.

We experimented with a similarity Procrustes variant that additionally fits an isotropic scale factor. It reduced the centroid-fit residual from 1.75 to 0.49 but collapsed label-transfer accuracy from 0.68 to 0.19 — the isotropic scale over-compresses the ATAC cloud, destroying within-cluster geometry. Rotation-only is the configuration we report.

## Entropic Sinkhorn alignment

On the post-Procrustes latents we run entropic OT. The cost matrix is the squared Euclidean distance between every RNA cell latent and every ATAC cell latent, divided by the **median** of the off-diagonal cost (we found median-normalization substantially more robust than max-normalization, whose scale is dominated by the tail). We then solve
$$\min_{P \in \Pi(\mu, \nu)} \langle C, P \rangle - \epsilon H(P)$$
via POT's sinkhorn (not the log-domain variant — the latter OOMs at 11k × 11k on 16 GB RAM because `scipy.logsumexp` allocates multi-GB intermediate buffers). We use $\epsilon = 0.05$ and, for the per-cell uncertainty analysis, restrict to a seeded 5000-cell subsample of each modality to keep the full Sinkhorn plan in memory. The downstream metrics FOSCTTM, label transfer, and ARI are computed on the full (not subsampled) aligned embeddings because they don't require the transport plan.

## Cluster-resolved uncertainty

For each source cell $i$ the transport plan row $P_{i \cdot}$ is a distribution over partner cells. Naively taking its entropy $H(P_{i \cdot}) = -\sum_j P_{ij} \log P_{ij}$ gives a signal that anti-correlates with true alignment error (Spearman $\rho = -0.65$ on PBMC) because cells in dense clusters have many valid within-cluster candidates and a closer true partner. We replace this with **cluster-resolved entropy**: marginalize the plan over the partner's cluster labels,
$$p(c | i) = \frac{\sum_{j : \text{cluster}(j) = c} P_{ij}}{\sum_j P_{ij}}, \qquad H_{\text{cluster}}(i) = -\sum_c p(c|i) \log p(c|i) \;/\; \log K$$
normalized to $[0, 1]$ by the number of partner clusters $K$. This marginalizes away within-cluster ambiguity and tracks only the cross-cluster alignment confidence. The argmax over $c$ is the method's predicted partner cluster for cell $i$, and its accuracy against the true partner cluster (recovered from `pair_idx` at evaluation time) is the primary validation of the UQ signal.

## Evaluation metrics

- **FOSCTTM** (Fraction Of Samples Closer Than the True Match): for each cell in modality $A$, the fraction of cells in modality $B$ that are closer to it (Euclidean on the aligned latents) than its true partner. We report the symmetric mean $\frac{1}{2}(\text{A→B} + \text{B→A})$. Random $\approx 0.5$, perfect = $0$.
- **Label transfer accuracy**: a $k$=15 nearest-neighbor classifier trained on one modality and evaluated on the other, in the aligned latent space.
- **Joint clustering ARI**: KMeans ($k$ = number of Leiden clusters) on the concatenated $[Z^{\text{RNA}}; Z^{\text{ATAC}}]$ latent, vs. the ground-truth Leiden labels.
- **Argmax cluster accuracy**: fraction of source cells whose cluster-resolved argmax partner cluster matches the ground-truth partner's cluster.
- **AUROC(H_cluster → wrong-cluster cell)**: for the subset of cells where the argmax cluster is wrong, can cluster entropy rank-order them above the correct-argmax cells? Higher AUROC = better-calibrated UQ.

## Transparency note on cluster propagation

We flag one methodological detail that affects how these results should be interpreted. The cross-modal prediction target $y_{\text{cross}}$ defined above is the centroid, in the partner modality's LSI/PCA space, over cells sharing a cluster assignment. In our paired benchmarks this cluster assignment was obtained by clustering one modality (RNA) with Leiden and **propagating the labels to the other modality via the paired ground truth** — that is, cell $i$ in RNA and the corresponding cell $i$ in ATAC (matched via `pair_idx`) are assigned the same cluster label for the purpose of computing $y_{\text{cross}}$. This is implemented as a simple row-copy in `src/data/preprocess.py:263`:

```python
atac.obs["leiden"] = rna.obs["leiden"].values
```

**The training loop itself never reads `pair_idx`.** We verified by inspection that `pair_idx` appears in `src/training/run_experiment.py` only in the post-training evaluation block (FOSCTTM, label transfer, ARI, cluster entropy evaluation), and not in `src/training/dataloader.py`, `src/training/train_ibvae.py`, or the IB-VAE loss computation. The dataloader exposes only `x`, the reconstruction target, and the precomputed `y_cross` to the training loop.

However, because `y_cross` was precomputed using paired cluster labels, the training signal does carry **cluster-level** paired information (not cell-level). For a fully unpaired application (two modalities profiled on different cell populations), a practitioner would need to replace the paired-GT cluster propagation with one of: (a) clustering each modality independently and then finding cluster correspondences via an external label source or a separate matching step, (b) using known biological cell-type labels from either modality as the cluster label, or (c) the same joint latent matching run iteratively in an EM-like outer loop. We leave these extensions to future work.

This caveat is important because it bounds the claims we can make: the FOSCTTM and label-transfer numbers reported in Table 1 reflect the method's ability to recover **cell-level pairing given cluster-level correspondences**, not to recover pairing from nothing. The strongest positive claim we make — and the one supported by the cross-tissue negative control in §"Cross-tissue negative control" — is that the cluster-resolved entropy is a well-calibrated per-cell uncertainty signal *given* the paired benchmark setup.

## Implementation and reproducibility

All code is available at the project repository, released under the MIT license. All seeds are set deterministically (NumPy, PyTorch, and cuDNN). Training uses a single NVIDIA RTX 3080 GPU. Wall-clock times are ≈ 120 seconds per dataset per seed for the full the method pipeline. Multi-seed results (reported as mean ± std over 3 seeds) are computed by running the pipeline with seeds 0, 1, 2 and aggregating via `scripts/aggregate_seeds.py`. The uniPort baseline is run from a separate Python virtual environment (`venv_uniport/`) with `numpy<2` and `anndata 0.11.4` because uniPort 1.3 uses `np.Inf`, which was removed in NumPy 2.0.
