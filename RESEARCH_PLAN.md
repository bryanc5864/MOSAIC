# Calibrated Per-Cell Uncertainty for Unpaired Single-Cell Multi-Omics Integration

**Field**: Computational genomics — single-cell multi-omics integration
**Constraints**: Data — public datasets only, <10 GB total on disk | Compute — 1× NVIDIA RTX 3080 (10 GB VRAM), local workstation
**Date**: 2026-04-10
**License**: MIT — Bryan Cheng, Austin Jin, Joshua Chang, Brendan Lo, 2026
**Target venue**: RECOMB 2026 (Thessaloniki); backups Cell Systems / Genome Research / JCB

> **Note on compute deviation from the author-provided plan.** The original plan assumed 10× RTX 2080 Ti. The executing environment has one RTX 3080 (10 GB). Implications, baked in below:
> - IB-VAEs are trained sequentially (per modality, per dataset) rather than in parallel.
> - OT for the cross-atlas experiment uses minibatch Sinkhorn + subsampling instead of the full 100K×80K cost matrix (which would need >60 GB RAM).
> - Wall-clock budget grows from ~2 days to ~1 week. Timeline in §9 reflects this.
> - 500K-cell scalability row in the original Exp 5 is replaced by a scaling *curve* up to the largest size the hardware allows (likely ~50–100K cells), with extrapolation reported honestly.

---

## 1. Abstract

Unpaired single-cell multi-omics integration — aligning scRNA-seq and scATAC-seq cells measured in different experiments — is routinely performed by methods (SCOT, uniPort, CMOT, MultiVI) that return point-estimate matchings with no per-cell confidence. Downstream analyses therefore cannot tell a trustworthy alignment from an artifact. the method closes this gap by combining **Information Bottleneck (IB) autoencoders** — which pre-filter each modality to retain only cross-modally predictive features — with **entropic optimal transport** (Sinkhorn), producing a transport plan whose per-row Shannon entropy serves as a calibrated per-cell alignment uncertainty score. We prove, under sub-Gaussian embedding assumptions, that alignment entropy converges to zero for shared cell types and remains bounded away from zero for cell types absent from one modality. Empirically, we evaluate on four paired ground-truth multi-omics datasets (10x Multiome PBMC, 10x Multiome Brain, SHARE-seq skin, CITE-seq PBMC), benchmark against seven baselines, demonstrate that alignment entropy correctly detects missing cell types (AUROC), and apply the method to an unpaired cross-atlas setting (Tabula Sapiens subset × ENCODE scATAC) to surface putative enhancer-priming regulatory states.

## 2. Background & Motivation

**The unpaired multi-omics alignment problem.** Modern single-cell assays profile distinct molecular layers (transcriptome, chromatin accessibility, surface protein, methylation) but rarely jointly — most atlases are single-modality. Integrating these into a common reference requires matching cells across modalities without ground-truth correspondence. Optimal-transport–based methods (SCOT [Demetci et al. 2022], uniPort [Cao et al. 2022], CMOT [Dai et al. 2023]) and deep generative models (MultiVI [Ashuach et al. 2023], LIGER [Welch et al. 2019], Harmony [Korsunsky et al. 2019], Seurat WNN [Hao et al. 2021]) are now standard. None quantify the confidence of individual matchings.

**Why per-cell uncertainty matters.** Silent misalignment corrupts joint clustering, regulatory inference, and spatial deconvolution downstream. A calibrated per-cell score would enable users to (a) filter low-confidence matches, (b) detect cell populations present in one modality but absent in the other, and (c) distinguish genuine biology from batch effects that resemble alignment uncertainty.

**The gap the method fills.** Entropic OT naturally produces a probability distribution over candidate matches. Prior work treats this as a computational convenience to be rounded into a hard assignment. We treat it as a first-class output: the Shannon entropy of each row *is* the per-cell alignment uncertainty, and under suitable conditions it has a formal interpretation. To make this work in practice, the latent space fed into OT must be clean — IB regularization (Alemi et al. 2017; Tishby & Zaslavsky 2015) provides the mechanism to discard modality-specific noise.

**Why no existing method combines IB + OT with per-cell UQ.** IB-VAEs are used for disentanglement and representation compression, not as a preprocessor for transport. OT methods in multi-omics use raw or VAE embeddings without a bottleneck, so their cost matrices are noisy and the resulting transport plans are high-entropy for reasons unrelated to true biological ambiguity. Coupling the two cleanly, and giving the row-entropy a formal meaning, is novel.

## 3. Technical Approach

### 3.1 Overview

The the method pipeline has four stages:

```
Raw counts ──► Preprocess ──► IB-VAE (per modality) ──► Entropic OT ──► H_i + matches
  (genes,     (HVG, LSI,      (64-dim latent z,          (N×M plan P    (per-cell
   peaks)      ZINB/BCE        cross-modal pred head)    via Sinkhorn)   entropy UQ)
```

**Stage 1 — Preprocessing.** Per-modality standard pipelines: scanpy for RNA (HVG → log-normalize), TF-IDF + LSI for ATAC. Converts raw count matrices to dense feature matrices of fixed dimensionality.

**Stage 2 — IB-VAE encoding.** Two modality-specific variational autoencoders with an information bottleneck regularizer. Each encoder compresses inputs to a 64-dim latent conditioned on a cross-modal prediction objective, so that latent coordinates retain only features that are predictive of the *other* modality.

**Stage 3 — Entropic OT alignment.** A Sinkhorn solver computes a coupling P ∈ ℝ^{N×M} between latent embeddings of the two modalities, minimizing expected squared-Euclidean cost plus entropic regularization ε·KL(P ‖ a⊗b).

**Stage 4 — Entropy scoring and downstream.** For each source cell i, H_i = −Σ_j P̃_ij log P̃_ij (where P̃_i is the row normalized to a probability distribution), divided by log M for scale invariance. Output: aligned embeddings + per-cell H_i ∈ [0, 1].

### 3.2 Architecture / Algorithm

#### 3.2.1 IB-VAE (per modality)

- **Encoder**: `input_dim → 512 → 256 → (μ, log σ²)`, where the latent head produces a 64-dim diagonal Gaussian posterior `q(z | x) = 𝒩(μ(x), diag σ²(x))`. MLP, GELU activations, LayerNorm between hidden layers.
- **Decoder**: `64 → 256 → 512 → input_dim`.
  - **RNA decoder loss**: ZINB (zero-inflated negative binomial) — parameters (μ, θ, π) produced by three heads sharing the trunk. Implementation via `scvi.distributions.ZeroInflatedNegativeBinomial`.
  - **ATAC decoder loss**: Bernoulli / BCE on binarized peak accessibility.
- **Cross-modal prediction head**: a linear layer `z_A → ŷ_B` where `y_B` is a low-dim summary of the other modality (for RNA, gene-activity scores computed from ATAC peaks; for ATAC, the top HVG expression vector). Loss: MSE on standardized `y_B`. This is the term that forces the bottleneck to keep cross-modally predictive axes.
- **Full IB-VAE loss** (per modality A, with target y_B):
  ```
  L_A = L_recon(x, x̂)                           (ZINB or BCE)
       + λ_pred · ‖ŷ_B − y_B‖²                    (cross-modal prediction)
       + β · KL(q(z|x) ‖ 𝒩(0, I))                 (IB bottleneck)
  ```
  Default β = 0.01, λ_pred = 1.0. β is tuned in Exp 4.

#### 3.2.2 Entropic OT via Sinkhorn

- Inputs: source embeddings Z_A ∈ ℝ^{N×64}, target embeddings Z_B ∈ ℝ^{M×64}, uniform marginals a = 1_N/N, b = 1_M/M.
- Cost matrix C_ij = ‖z^A_i − z^B_j‖². Normalized to unit max for numerical stability.
- Sinkhorn problem: P* = argmin ⟨C, P⟩ − ε·H(P), s.t. P1 = a, Pᵀ1 = b.
- Solver: `ot.sinkhorn` (POT library). Default ε = 0.05 (on normalized cost). Max 1000 iterations, tol 1e-7. Log-domain Sinkhorn for numerical stability.
- GPU path: if N·M ≤ 2×10⁸ entries (fits in ~1.5 GB VRAM at fp32), run on device. Otherwise, chunked CPU Sinkhorn or minibatch OT with averaging across minibatches (as in uniPort).

#### 3.2.3 Alignment entropy

- For each row i of the coupling, normalize to a probability: p_ij = P_ij / Σ_j P_ij (since uniform marginals, p_ij = M·P_ij).
- Shannon entropy: H_i = −Σ_j p_ij log p_ij.
- Normalize: H̃_i = H_i / log M ∈ [0, 1].
- Outputs: per-cell `entropy.npy`, matching embeddings, top-k matches per source cell.

### 3.3 Key Design Choices

1. **Why IB before OT?** Two alternatives were considered: (a) concatenate modalities and train a single joint VAE (MultiVI approach), (b) run OT directly on scVI embeddings. Both mix in modality-specific noise axes that inflate OT entropy for reasons unrelated to true biological ambiguity, breaking the theorem's assumptions. IB explicitly discards such axes.
2. **Why entropic OT (not exact OT or neural OT)?** Entropic OT yields a differentiable, stable coupling and — crucially — a distribution per row, which is exactly the per-cell uncertainty we need. Exact OT gives a 0/1 permutation, killing the UQ signal.
3. **Why uniform marginals?** Avoids fitting a prior over cell-type abundances, which would confound "absent cell type" signals. Non-uniform marginals are a natural extension but would introduce another tunable.
4. **Why 64-dim latent?** Matches common single-cell VAE practice (scVI default 10–30; we go slightly larger because cross-modal prediction needs headroom). Small enough that Sinkhorn cost-matrix distances are well-behaved.
5. **Why normalize entropy by log M?** Makes H̃ ∈ [0,1] comparable across datasets with different cell counts.

### 3.4 Training / Optimization

- **Optimizer**: AdamW, lr = 1e-3, weight_decay = 1e-4.
- **Schedule**: Cosine annealing over 200 epochs, 10-epoch linear warmup.
- **Batch size**: 512 cells.
- **Precision**: fp32 for stability with ZINB (fp16 loss components can blow up).
- **KL annealing**: β ramps linearly from 0 to target β over first 30 epochs to prevent posterior collapse.
- **Early stopping**: on held-out cross-modal prediction MSE (patience 15 epochs).
- **Seeds**: Every run seeded (numpy, torch, cuda). Recorded in config yaml.

## 4. Experimental Design

### 4.1 Experiments

Every experiment logs to `TRAINING_LOG.md` and updates `RESULTS.md`. Every metric is computed by code whose commit hash is pinned in the experiment config.

**Exp 1 — Alignment accuracy on paired data (primary).**
- *Objective*: Is the method competitive with or better than existing unpaired-integration methods when judged against ground-truth pairings?
- *Setup*: All four paired datasets, run sequentially: 10x Multiome PBMC → 10x Multiome Brain → CITE-seq PBMC → SHARE-seq skin. PBMC is the debugging + hyperparameter-tuning target; the other three are held for reporting (no hyperparameter tuning on them). Split into the two modalities, *drop pairing*, run the method, compare to true pairing.
- *Metrics*: FOSCTTM (lower = better), label-transfer accuracy (cell type classifier trained on RNA, evaluated on ATAC via aligned embedding, higher = better), ARI of joint clustering vs. known cell-type labels.
- *Success*: the method achieves FOSCTTM ≤ the best of (SCOT, uniPort) on at least 2 of 4 datasets, or a principled explanation if not.

**Exp 2 — Entropy calibration.**
- *Objective*: Is per-cell alignment entropy correlated with actual alignment error?
- *Setup*: Same paired dataset. For each cell, compute alignment error = latent distance to its true partner's latent. Correlate with H̃_i.
- *Metrics*: Spearman ρ between H̃_i and alignment error (positive means higher entropy = higher error = well-calibrated).
- *Success*: Spearman ρ ≥ 0.4 on PBMC. Also produce a calibration curve binning cells by H̃_i and plotting mean error per bin.

**Exp 3 — Missing cell type detection.**
- *Objective*: Does alignment entropy flag cells whose type is absent from the target modality?
- *Setup*: Remove one cell type from the ATAC modality (e.g., platelets from PBMC, which are small in count). Run the method. Measure whether RNA cells of the removed type have elevated H̃_i.
- *Metrics*: AUROC for (removed-type vs. kept-type) classification via H̃_i threshold.
- *Success*: AUROC ≥ 0.8 on at least one held-out type; report across all candidate holdouts.

**Exp 4 — Ablation study.**
- *Ablations*: (a) OT directly on raw features (no IB), (b) IB embeddings + nearest-neighbor matching (no OT), (c) IB-VAE without cross-modal prediction head, (d) β sweep ∈ {0.001, 0.01, 0.1, 1.0}, (e) ε sweep ∈ {0.01, 0.05, 0.1, 0.5}.
- *Metric*: FOSCTTM on PBMC.
- *Success criterion*: each component (IB, OT, cross-modal head) must be necessary — removing it should degrade FOSCTTM by a measurable, consistent amount.

**Exp 5 — Scalability curve.**
- *Setup*: Subsample PBMC to {2K, 5K, 10K}, measure wall-clock time for each pipeline stage. Report the largest N at which full-matrix Sinkhorn fits on the RTX 3080; document the minibatch-OT crossover.
- *Metrics*: wall-clock time, peak RAM, peak VRAM. (Note: scales down from original plan due to single-GPU constraint.)

**Exp 6 — Cross-atlas regulatory discovery (stretch).**
- *Setup*: Tabula Sapiens RNA subset (lung, ~10–20K cells) × ENCODE scATAC cortex subset. Apply the method, find high-entropy cell populations, GO-enrichment and enhancer-priming analysis.
- *Success*: at least one qualitatively interpretable high-entropy population surfaces, with biological interpretation supported by marker genes / known enhancer maps.

### 4.2 Baselines

Minimum viable comparison: **SCOT** (GW-OT) and **uniPort** (VAE + unbalanced OT), both of which directly target unpaired alignment. Stretch: MultiVI, Seurat v5 WNN.

Fair-comparison protocol:
- All methods receive the same preprocessed AnnData objects (same HVGs, same LSI dim for ATAC).
- Default hyperparameters from each method's paper.
- Same random seed for any stochastic step.
- Same evaluation code for all methods (single source of truth for FOSCTTM, etc.).

### 4.3 Ablation Studies

See Exp 4 above. Ablations are first-class experiments.

## 5. Dataset Strategy

### 5.1 Data Sources

| Dataset | Modalities | Cells | Source | Approx size |
|---|---|---|---|---|
| 10x Multiome PBMC (10k) | scRNA + scATAC (paired) | ~10K | 10x Genomics public demo | ~1.5 GB |
| 10x Multiome Brain | scRNA + scATAC (paired) | ~12K | 10x Genomics public demo | ~2 GB |
| SHARE-seq mouse skin | scRNA + scATAC (paired) | ~34K | GSE140203 (Ma et al. 2020) | ~4 GB |
| CITE-seq PBMC | scRNA + ADT protein (paired) | ~30K | GSE164378 (Hao et al. 2021) | ~1 GB |
| Tabula Sapiens (lung subset) | scRNA | ~10–20K | CZI cellxgene | ~1 GB |
| ENCODE scATAC (cortex subset) | scATAC | ~10–20K | ENCODE | ~1 GB |

**Scoping to the <10 GB budget**: all four paired datasets are in scope and will be run sequentially (user confirmation 2026-04-10). Order: (1) 10x Multiome PBMC (primary, smallest — debugging target), (2) 10x Multiome Brain, (3) CITE-seq PBMC, (4) SHARE-seq skin (largest — run last). SHARE-seq is subsampled to ~15K cells if raw download exceeds the budget. Tabula Sapiens is restricted to one tissue (lung).

### 5.2 Preprocessing Pipeline

- **scRNA-seq** (scanpy): filter cells with <200 genes; filter genes in <3 cells; normalize total counts to 1e4; log1p; select top 2000 HVGs via `scanpy.pp.highly_variable_genes(flavor='seurat_v3')`.
- **scATAC-seq** (episcanpy / scanpy with custom TF-IDF): filter cells with <1000 fragments; TF-IDF normalize peak matrix; select top 50K variable peaks; reduce to 50-dim LSI via TruncatedSVD (drop first component, standard ATAC practice).
- **Gene activity (for cross-modal head, ATAC → RNA)**: sum ATAC fragments within ±2 kb of gene TSS to produce a cell × gene activity matrix. Standardized per gene.
- **Paired vs. unpaired**: for benchmark experiments, keep pairings as metadata only; never pass them through the training loop.
- **Splits**: for each dataset, define a held-out test split (20% of cells, stratified by cell type) used only for final reporting. Validation split (10%) used for early stopping. All splits are seeded and saved to disk.

### 5.3 Data Validation

- After each preprocessing step: assert non-NaN, non-empty matrices, expected shapes, expected dtype.
- Visual spot-check: UMAP of each processed modality colored by cell type — sanity-check that clusters are present.
- Pair-leakage check (automated): a test that constructs a DataLoader from the training AnnData and verifies that no index maps a cell of modality A to its paired cell in modality B within any batch. Runs in CI for every experiment.
- Integrity hashes: every processed AnnData file is saved with an MD5 of its `.X`, logged in the experiment config.

## 6. Implementation Roadmap

### Milestone 1: Environment & data
- **Objective**: Reproducible Python env; PBMC dataset preprocessed and saved; validation pipeline working.
- **Deliverable**: `requirements.txt` (pinned), `data/processed/pbmc10k.h5ad`, passing data-validation tests.
- **Verification**: Running `python -m src.data.validate pbmc10k` prints all green, RNA and ATAC UMAPs saved to `figures/data_sanity/`.

### Milestone 2: IB-VAE implementation
- **Objective**: Working IB-VAE per modality with all three loss terms, testable end-to-end.
- **Deliverable**: `src/models/ib_vae.py` with unit tests on synthetic data.
- **Verification**: Unit test trains for 5 epochs on a toy dataset and shows strictly decreasing loss; cross-modal prediction head beats random baseline.

### Milestone 3: OT + entropy
- **Objective**: Sinkhorn wrapper, entropy scoring, top-k match extraction.
- **Deliverable**: `src/models/ot_align.py` with unit tests.
- **Verification**: On two Gaussian blobs (known matching), Sinkhorn transport plan recovers identity permutation and per-row entropy is near zero; when blobs are made identical, entropy is near uniform.

### Milestone 4: Evaluation metrics
- **Objective**: FOSCTTM, label transfer, ARI, entropy-error correlation, missing-type AUROC, all as pure functions.
- **Deliverable**: `src/evaluation/metrics.py` with unit tests against toy ground truth.
- **Verification**: Metric functions give expected values on hand-constructed toy data (identity matching → FOSCTTM 0, random matching → FOSCTTM ≈ 0.5).

### **Milestone 5: 🔍 Pre-Training Review Gate (mandatory)**
- **Objective**: Run `/review` in Pre-Training mode. Must pass before any training.
- **Deliverable**: `REVIEW_REPORT_PRE_TRAINING.md` with verdict 🟢 PASS (or 🟡 PASS WITH CONCERNS + written acknowledgment).
- **Verification**: No critical data-leakage or reproducibility issues open.

### Milestone 6: Exp 1 (primary alignment accuracy)
- **Objective**: Train IB-VAEs on PBMC, run OT, compute FOSCTTM / label transfer / ARI.
- **Deliverable**: `experiments/exp_001_pbmc_primary/` with config, logs, outputs; RESULTS.md updated.
- **Verification**: All three metrics computed from code and logged to TRAINING_LOG.md.

### Milestone 7: Baselines (SCOT, uniPort)
- **Objective**: Run baselines under fair-comparison protocol.
- **Deliverable**: Baseline entries in RESULTS.md with same metrics.
- **Verification**: Same evaluation code used for all methods.

### Milestone 8: Exp 2–4 (calibration, missing-type, ablations)
- **Deliverable**: Experiment dirs, figures, RESULTS.md updates.

### **Milestone 9: 🔍 Post-Results Review Gate (mandatory)**
- **Objective**: Run `/review` in Post-Results mode. Must pass before figures and paper drafting.
- **Deliverable**: `REVIEW_REPORT_POST_RESULTS.md` with no fabrication findings.

### Milestone 10: Figures and paper drafts
- **Deliverable**: `figures/` populated with publication-quality PNG+PDF; `paper/{introduction,methods,results,discussion}.md` drafted.

## 7. Evaluation Criteria

### 7.1 Primary Metrics

1. **FOSCTTM on 10x Multiome PBMC (Exp 1)** — the single number the paper hinges on.
2. **Spearman ρ between H̃_i and true alignment error (Exp 2)** — validates the central novelty claim.
3. **AUROC for missing cell type detection (Exp 3)** — validates the usefulness claim.

### 7.2 Secondary Metrics

- Label transfer accuracy, ARI, wall-clock, VRAM usage, sensitivity to β and ε.
- Qualitative: UMAP overlays, entropy-binned error curves.

### 7.3 Definition of Success

- **Full success**: the method beats or matches SCOT/uniPort on FOSCTTM AND entropy calibration ρ ≥ 0.4 AND missing-type AUROC ≥ 0.8 AND at least one ablation shows a clear component contribution AND cross-atlas application surfaces an interpretable finding.
- **Partial success**: The UQ component (ρ and AUROC) is calibrated, even if raw alignment accuracy only ties existing methods. This is still a clear paper because no baseline provides UQ at all.
- **Failure**: Alignment entropy is uncorrelated with error, and the method is strictly worse than baselines on FOSCTTM. In that case, report honestly, diagnose (is the IB collapsing? is ε poorly chosen?), and pivot scope — most likely to an analytical / theoretical paper about when entropic OT couplings do and do not carry calibrated uncertainty.

## 8. Risk Mitigation

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Single-GPU OOM on full cost matrix | H | M | Chunked Sinkhorn + minibatch OT fallback; cap max N per experiment; document crossover empirically. |
| IB bottleneck collapses posterior | M | H | KL warmup (30 epochs) + β sweep in Exp 4; monitor `KL / cells` during training. |
| Alignment entropy uncalibrated in practice | M | H | Exp 2 is exactly the test. If it fails, diagnose with ablations, consider temperature scaling of Sinkhorn or learned calibration layer. |
| Dataset download fails / license issue | L | M | Primary dataset (10x Multiome PBMC) is a long-standing public demo; backup is GEO downloads. |
| SHARE-seq / Tabula Sapiens >10 GB | M | L | Tissue / cell-count subsampling documented in §5.1; budget applies to *processed* footprint. |
| Baseline install conflicts (POT vs scvi-tools versions) | M | M | Separate conda environments for the method vs. SCOT / uniPort if needed; shared evaluation env reads outputs from disk. |
| Enhancer-priming discovery (Exp 6) is unsupported | M | L | Exp 6 is a stretch; negative result is still reported. Primary paper claims do not depend on it. |

## 9. Timeline

Compute deviation (single GPU, local) reflected throughout. Estimated calendar time rather than wall-clock.

| Phase | Estimated Duration | Dependencies |
|---|---|---|
| Milestone 1 — env + data | Days 1–3 | — |
| Milestone 2 — IB-VAE impl | Days 3–5 | M1 |
| Milestone 3 — OT + entropy | Days 5–6 | — |
| Milestone 4 — metrics | Days 5–6 | — |
| Milestone 5 — **Pre-Training Review** | Day 7 | M1–M4 |
| Milestone 6 — Exp 1 on PBMC (debugging target) | Days 8–10 | M5 |
| Milestone 6b — Exp 1 on Brain / CITE-seq / SHARE-seq | Days 10–14 | M6 |
| Milestone 7 — baselines on all four datasets | Days 13–16 | M1 |
| Milestone 8 — Exp 2–4 on PBMC (+ Exp 2 on others) | Days 15–19 | M6 |
| Milestone 9 — **Post-Results Review** | Day 20 | M8 |
| Milestone 10 — figures + paper | Days 20–24 | M9 |
| Exp 6 (stretch) | Days 22–26 | M6 (parallel to M10) |

Total: ~4 weeks calendar (extended from 3 to accommodate running all four paired datasets sequentially; compare to ~2 days wall-clock in the original 10-GPU plan).

## 10. Expected Deliverables

- [ ] Working codebase under MIT license, `src/` tree with tests.
- [ ] `requirements.txt` (pinned).
- [ ] Processed AnnData files under `experiments/<exp>/data/` (or `data/processed/` for shared).
- [ ] All trained IB-VAE checkpoints and seed configs.
- [ ] `RESULTS.md` with every experiment traced to code.
- [ ] `TRAINING_LOG.md` with every run logged in canonical format.
- [ ] `REVIEW_REPORT_PRE_TRAINING.md` and `REVIEW_REPORT_POST_RESULTS.md`.
- [ ] Publication-quality figures (PNG + PDF) in `figures/`.
- [ ] Paper section drafts in `paper/` (introduction, methods, results, discussion).
- [ ] Informal writeup of the convergence theorem and proof sketch in `paper/theorem.md`.
