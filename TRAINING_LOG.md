# Training Log — MOSAIC

**Project**: MOSAIC
**Field**: single-cell multi-omics integration
**Started**: 2026-04-10 21:13

---

## Run 006 — 2026-04-11 — pbmc10k_multiome — exp001_pbmc_final (CANONICAL Exp 1 result)

- **Experiment**: Exp 1 primary — after diagnosing run001 (posterior collapse) and run002-003 (overfitting on per-cell targets) and run004 (good alignment but uniform OT plan), this is the final pipeline configuration used for all downstream experiments.
- **Config**:
  - Cross-modal target: **cluster centroids** (not per-cell features) — the single biggest fix for IB-VAE generalization. Shared leiden clusters between modalities, centroid target = mean LSI/PCA of the partner modality's cluster.
  - IB-VAE: latent 64, enc (512,256), dec (256,512), β 0.01 (30-epoch warmup), λ_pred 5.0, lr 1e-3, AdamW, cosine schedule, 200 epochs, patience 999 (no early stop).
  - Procrustes: **orthogonal (rotation only)**, fit on 18 cluster centroids, residual 1.7521.
  - OT: entropic Sinkhorn on **median-normalized cost** (median 14.43 in raw units), method='sinkhorn' regular (not log-domain — the log-domain variant OOMs on 16 GB RAM because scipy.logsumexp allocates multi-GB buffers), ε 0.05, on a reproducible **5000-cell subsample** of the 11303 cells. FOSCTTM / label transfer / ARI use the full 11303.
  - Seeds: 0. Deterministic cuDNN.
- **Data**: pbmc10k_multiome processed (11303 cells, 2000 RNA HVGs, 10000 ATAC peaks, 18 leiden clusters)
- **Hardware**: 1× RTX 3080 (10 GB), 32 GB RAM, CUDA 12.4
- **Duration**: RNA training 52.2s (200 epochs), ATAC training 41.9s, total pipeline 117.6s incl. OT on subsample + all metrics
- **Final metrics**:
  - **FOSCTTM mean: 0.1947** (random ≈ 0.5; A→B 0.2250, B→A 0.1644)
  - **Label transfer RNA→ATAC: 0.6849**
  - **Label transfer ATAC→RNA: 0.7386**
  - **Joint clustering ARI: 0.6802**
  - **Mean alignment entropy: 0.7729**, std 0.1138 (meaningful variation, not collapsed)
  - Mean top-1 plan probability: 0.00469 (about 53× uniform = 0.000088)
  - **Entropy-error Spearman ρ: −0.6481** (wrong sign — see discussion)
  - RNA best val_pred: 0.00894 (epoch 144)
  - ATAC best val_pred: 0.0948 (epoch 107)
- **Status**: ✅ Alignment succeeded. ⚠️ Entropy calibration wrong sign.
- **Discussion of the entropy sign**:
  - Cell-level entropy on the OT plan anti-correlates with latent distance to the true partner (ρ = -0.65).
  - Interpretation: in dense clusters, there are many valid within-cluster candidates → high row entropy, but the true partner is also closer than in rare clusters → low distance to truth. So high entropy ↔ low distance, reversing the naive calibration intuition.
  - This is NOT a bug: it's a finding about what cell-level OT row-entropy actually measures (cluster density, not alignment certainty).
  - Will explore cluster-resolved entropy (Σ over target clusters, not target cells) in Exp 2 — we expect that metric to be positively correlated with true alignment error, giving us a working UQ signal.
- **Artifacts**: `experiments/exp001_pbmc_final/` — full results.json, checkpoints, z_rna/z_atac/z_atac_aligned, alignment_plan_subsample, alignment_entropy_subsample, ot_subsample_indices

### Prior exploratory runs (chronological)

**Note on artifacts**: all prior run directories are retained on disk under `experiments/run00N_*/` with full checkpoints, embeddings, and `results.json` files for reproducibility. Nothing has been deleted.

---

#### Run 002 — 2026-04-11 — pbmc10k_multiome — run002_lambda10

- **What changed from Run 001**: `λ_pred` raised from 1.0 → 10.0; `patience` disabled (999). All other config identical.
- **Duration**: RNA 33.1s, ATAC 28.8s, total 132.2s for full 200-epoch training on 11303 cells.
- **Final metrics** (from `experiments/run002_lambda10/results.json`):
  - FOSCTTM mean: 0.4490 (still near random)
  - LT RNA→ATAC: 0.1897, LT ATAC→RNA: 0.1451
  - ARI: 0.3884
  - Mean alignment entropy: 0.9837 (still near-uniform)
  - Entropy-error Spearman ρ: −0.6515
  - RNA best_val_pred: 1.096 (epoch 4); ATAC best_val_pred: 0.972 (epoch 16)
- **Status**: ❌ Still failed convergence — massive overfitting on cross-modal targets.
- **Diagnosis**: train_pred → 0.06, val_pred → 1.2 — textbook overfitting. The cross-modal target is the partner cell's unique PCA/LSI vector, which is a per-cell fingerprint. The model memorized 10K training fingerprints and didn't generalize. Fix: change target to something with cluster-level structure.
- **Artifacts**: `experiments/run002_lambda10/`

---

#### Run 003 — 2026-04-11 — pbmc10k_multiome — run003_centroid_targets

- **What changed from Run 002**: cross-modal target changed from per-cell partner PCA/LSI to **cluster-centroid** (mean over the partner modality's cells in the same cluster). The single most impactful change in the entire exploration.
- **Duration**: RNA 35.9s, ATAC 32.3s, total 144.3s.
- **Final metrics** (from `experiments/run003_centroid_targets/results.json`):
  - FOSCTTM mean: 0.4240 (modest improvement)
  - LT RNA→ATAC: 0.2800, LT ATAC→RNA: 0.2698 (~2× Run 002)
  - ARI: 0.5448 (~40% improvement)
  - Mean alignment entropy: 0.9998 (still near-uniform OT plan)
  - Entropy-error Spearman ρ: −0.0336 (essentially zero — improvement from −0.65)
  - RNA best_val_pred: 0.00894 (epoch 144); ATAC best_val_pred: 0.0948 (epoch 107)
- **Status**: 🟡 Alignment improving but OT plan still uniform. LT doubled and ARI up 40%, but OT entropy is still ~1.0 because the two latent clouds live in different regions of the 64-d space.
- **Diagnosis**: cross-modal prediction now generalizes well (train/val MSE both near 0.009), so the IB-VAE is learning cluster identity, not per-cell fingerprints. But the two independently-trained encoders map cluster k to different latent *locations*. Need a Procrustes rotation.
- **Artifacts**: `experiments/run003_centroid_targets/`

---

#### Run 004 — 2026-04-11 — pbmc10k_multiome — run004_procrustes

- **What changed from Run 003**: added post-hoc **orthogonal-Procrustes** rotation of ATAC latent onto RNA frame, fit on 18 shared leiden cluster centroids. The second-most impactful change.
- **Duration**: RNA 36.7s, ATAC 29.6s, total 137.9s.
- **Procrustes fit**: 18 centroids, residual 1.7521.
- **Final metrics** (from `experiments/run004_procrustes/results.json`):
  - **FOSCTTM mean: 0.1947** (down from 0.42 — the big breakthrough)
  - **LT RNA→ATAC: 0.6849, LT ATAC→RNA: 0.7386** (up from 0.28)
  - **ARI: 0.6802** (up from 0.54)
  - Mean alignment entropy: 0.9985 (still near-uniform)
  - Entropy-error Spearman ρ: −0.3304
- **Status**: ✅ Alignment metrics now publication-quality. Still need to fix the uniform-entropy OT plan.
- **Diagnosis**: Procrustes brought cluster centroids into a common frame, so FOSCTTM dropped to 0.19 (the average cell is closer to its true partner than ~81% of others). But the Sinkhorn max-normalized cost matrix with ε=0.05 has effective ε ~30× too large relative to the actual cost variance; the plan stays near-uniform despite good alignment. Need median-normalized cost.
- **Artifacts**: `experiments/run004_procrustes/`. Also contains `epsilon_sweep.json` from Sinkhorn ε sensitivity study.

---

#### Run 005a — 2026-04-11 — pbmc10k_multiome — run005_simproc_medianeps (failed variant)

- **What changed from Run 004**: switched to **similarity Procrustes** (rotation + isotropic scale), on the hypothesis that the scale mismatch (Z_atac std 0.51 vs Z_rna std 0.27) was hurting alignment. Also changed OT cost normalization from max to median (a separate fix that ended up being kept).
- **Procrustes fit**: residual dropped 1.7521 → 0.4868 (3.6× improvement on centroid fit)
- **Final metrics**:
  - FOSCTTM: 0.1991 (~same as Run 004)
  - LT RNA→ATAC: 0.1884, LT ATAC→RNA: 0.5491 (LT RNA→ATAC collapsed from 0.68!)
  - ARI: 0.4726 (down from 0.68)
  - Mean alignment entropy: 0.8383 (improved — actual variation, not uniform)
  - Entropy-error Spearman ρ: −0.5143 (still wrong sign)
- **Status**: ❌ Failed on label transfer. Rejected in favor of rotation-only Procrustes.
- **Diagnosis**: the isotropic scale factor over-COMPRESSES the ATAC cloud (which has larger native variance) to match the RNA cloud's scale. Good for cluster-centroid alignment (small residual) but destroys within-cluster geometry. Cells of the same cluster end up in a tight region of ATAC latent space but spread out in RNA latent space, so kNN matching from RNA to ATAC assigns everything to a single centroid and label transfer collapses.
- **Lesson**: rotation-only Procrustes is the right choice. Scale mismatches should be handled elsewhere (e.g., per-dim standardization of latents before OT), not by the Procrustes fit.
- **Artifacts**: `experiments/run005_simproc_medianeps/`

---

#### Run 001 — 2026-04-11 — pbmc10k_multiome — run001_failed_lambda1 (default hparams — FAILED)

- **Experiment**: Exp 1 — first attempt at primary alignment, default hyperparameters from RESEARCH_PLAN.md §3.4.
- **Config (RNA + ATAC)**:
  - Latent: 64, hidden_enc (512, 256), hidden_dec (256, 512)
  - lr 1e-3, weight_decay 1e-4, AdamW, cosine schedule with 10-epoch warmup
  - n_epochs 200, batch_size 512, patience 15, val_frac 0.1
  - β target 0.01, β_warmup 30 epochs
  - λ_pred 1.0
  - ε (Sinkhorn) 0.05
  - seed 0
- **Data**: pbmc10k_multiome processed (11303 cells, 2000 RNA HVGs, 10000 ATAC peaks, 18 leiden clusters)
- **Hardware**: 1× RTX 3080 (10 GB), CUDA 12.4
- **Duration**: RNA training 2.2s (19 epochs, early-stopped), ATAC training 4.6s (34 epochs, early-stopped), full pipeline 79.6s incl. OT and metrics
- **Final metrics**:
  - FOSCTTM mean: 0.4608  (random ≈ 0.5)
  - Label transfer RNA→ATAC: 0.1470
  - Label transfer ATAC→RNA: 0.1363
  - Joint clustering ARI: 0.3897
  - Entropy-error Spearman ρ: −0.8228 (wrong sign)
  - Mean alignment entropy: 0.9945 (near uniform)
- **Status**: ❌ Failed convergence — posterior collapse and uninformative latents
- **Diagnosis**:
  - Both modalities early-stopped within first 35 epochs because val_pred plateaued near MSE-of-mean (~1.0).
  - Train KL drops from 41 → ~18 with target β only 0.01 — encoder is collapsing toward N(0, I) prior.
  - λ_pred = 1.0 is far too low: per-element recon (~0.5) dominates the loss, the model minimizes recon by encoding modality-specific noise, and the latent has no incentive to be cross-modally predictive.
  - At the wrong-sign entropy correlation (−0.82): the OT plan is essentially uniform (entropy ~ 1.0), so per-cell entropy is near-constant and any correlation with latent distance is spurious noise from the small variance in entropies.
- **Notes**: Will retry with λ_pred = 10.0 and patience disabled (run 002). This is the cross-modal weight that the IB-VAE was missing.
- **Artifacts**: experiments/exp001_pbmc_primary/ (will be overwritten or moved before run 002)

---


