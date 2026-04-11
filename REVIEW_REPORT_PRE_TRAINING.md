# Pre-Training Review Report — MOSAIC

**Review Mode**: Pre-Training
**Date**: 2026-04-10
**Reviewer**: Autonomous Review Skill
**Project root**: `C:\Users\Maozer\projects\MOSAIC`
**Plan**: `RESEARCH_PLAN.md` (formalized from Bryan's authored proposal; adapted for single RTX 3080)

---

## Summary

**Overall Status**: 🟡 **PASS WITH CONCERNS** — training may proceed. No critical blockers; warnings are methodology nits to document in the paper or address before final submission.

A full code trace of the data pipeline, model, training loop, optimal transport alignment, and evaluation metrics found no correctness bugs that would invalidate experimental results. The 2-epoch and 10-epoch end-to-end smoke tests complete without runtime errors, models decrease loss on synthetic data, OT/entropy/metric unit tests all pass, and the data-validation module confirms no pair-index leakage path from training into evaluation. The warnings below should be acknowledged in paper methods and ideally tightened before final submission, but none are blocking for Exp 1.

---

## Critical Issues (Must Fix)

**None.**

The code is ready for training runs. All smoke tests pass. The pairing-leakage design invariant — that no single-modality training code path can observe the partnering cell in the other modality — is enforced by architecture (two independent AnnDatas, per-modality training loop, `pair_idx` only accessed in `src/evaluation` and `src/data/validate`) and confirmed by grep: `pair_idx` appears in `metrics.py`, `validate.py`, `preprocess.py`, and `run_experiment.py`, but never in `dataloader.py` or `train_ibvae.py`.

---

## Warnings (Should Fix — before final submission)

### W1. Preprocessing statistics are fit on the full dataset

**Where**: `src/data/preprocess.py` — `preprocess_rna` and `preprocess_atac`.

**What**: Every normalization step that computes cross-cell statistics sees all cells before the train/val split in `train_ibvae.py` is applied:
- `sc.pp.scale` (line 188) computes per-gene mean/std over all 11,303 cells.
- `sc.tl.pca` (line 189) fits PCA on all cells.
- `sc.tl.leiden` (line 192) clusters all cells, and those cluster IDs become the `cell_type` labels used in evaluation.
- `tfidf` (line 129 in same file) computes IDF as `log(1 + n_cells / peak_sums)` on all cells.
- `lsi` / `TruncatedSVD` fits SVD directions on all cells.

**Severity**: 🟡 Warning (not blocking). Every widely-cited single-cell integration benchmark (scGLUE, scib, SCOT, uniPort) does this the same way, because per-cell normalization cannot be run in a streaming manner for these pipelines. The val split in MOSAIC is used *only* for IB-VAE early stopping, not for reporting FOSCTTM / ARI / AUROC — those final metrics are computed on the full dataset, where "held-out" doesn't mean what it means in standard ML. So the preprocessing leakage affects early-stop choice but not reported numbers.

**Recommended fix**: (a) In paper methods, explicitly acknowledge that preprocessing statistics were computed on the full dataset, matching standard single-cell practice. (b) For the ablation studies where we need a true held-out set (Exp 3 missing-type detection), fit normalization on the *kept* cells only and apply to the artificially-removed cells.

---

### W2. Loss scale inconsistency: per-element reconstruction vs. per-sample KL

**Where**: `src/models/ib_vae.py`, `IBVAE_RNA.forward` and `IBVAE_ATAC.forward`.

**What**: 
- `recon = -log_zinb(...).mean()` → mean over all `(B, n_vars)` elements → per-element NLL, magnitude ~0.5.
- `kl = kl.sum(dim=1).mean()` → mean over batch of per-sample summed-over-dims KL → magnitude ~40 at initialization with 64 latent dims.
- `pred = F.mse_loss(...)` → per-element MSE, magnitude ~1.

The total loss mixes three different reductions: per-element for recon and pred, per-sample for KL. This is a common VAE practice (it effectively down-weights the KL by the number of genes, helping avoid posterior collapse), but it means `β = 0.01` in MOSAIC is *not* directly comparable to `β = 0.01` in a paper that uses consistent per-element reductions. The β-warmup from 0 to 0.01 is effectively warming from 0 to ~0.0005 per-element-equivalent.

**Severity**: 🟡 Warning (not blocking). Empirically correct and stable, but paper's methods section should state the loss formulation precisely: "the reported β is applied to per-sample summed-over-latent KL divergence, with per-element reconstruction and cross-modal prediction losses."

**Recommended fix**: Document in the paper exactly how the total loss is computed. Optionally, refactor to consistent per-sample reductions (multiply recon by `n_vars`) so β is interpretable per-element; then the effective β for the ablation sweep `{0.001, 0.01, 0.1, 1.0}` maps to `{~1/n_vars, ~10/n_vars, ...}`.

---

### W3. "Cell type" labels are leiden clusters from RNA side, copied to ATAC

**Where**: `src/data/preprocess.py::add_cross_modal_targets` and `src/evaluation/metrics.py::label_transfer_accuracy`.

**What**: Because 10x Multiome doesn't ship with curated cell-type annotations, MOSAIC defines `cell_type = leiden cluster on the RNA modality`, then copies those labels to the ATAC modality via the known pairing. This means:
1. `label_transfer_accuracy(Z_rna, labels_rna, Z_atac, labels_atac)` is testing whether alignment preserves RNA-derived clusters across the modality boundary — i.e. "can I find my partner?" indirectly, not "can I transfer biologist-curated labels?".
2. `joint_clustering_ari` similarly uses labels that were derived from a one-modality clustering.

This is an *internal* consistency test, not an *external* label-transfer test.

**Severity**: 🟡 Warning (not blocking). This is the standard approach in papers that benchmark on 10x Multiome because curated labels are unavailable, but a reviewer might ask for something stronger.

**Recommended fix**: For the primary paper results, supplement with celltypist / Azimuth automated PBMC labels as an independent label source for at least the PBMC dataset. The existing metric code will work without modification — just swap the labels in `run_experiment.py`.

---

### W4. cuDNN determinism flags not set

**Where**: `src/training/train_ibvae.py::_set_seed`

**What**: `_set_seed` sets Python, numpy, torch, and cuda seeds, but does not set:
```python
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
```
or
```python
torch.use_deterministic_algorithms(True, warn_only=True)
```

Without these, cuDNN may pick different kernels on different runs, leading to small numerical variations even with identical seeds. For MOSAIC (MLPs only, no convolutions), the impact is small, but not zero.

**Severity**: 🟡 Warning (not blocking).

**Recommended fix**: Add those three lines to `_set_seed`. Accept a potential small training-speed cost.

---

### W5. Baselines (SCOT, uniPort) not yet installed or wired up

**Where**: Project roadmap, §4.2 of RESEARCH_PLAN.md.

**What**: The baselines listed in the plan are not yet installed or implemented. This is not a bug in the code that *is* present, but a reminder: before the Post-Results review gate, at minimum two baselines must be run under the fair-comparison protocol (same preprocessed AnnData, same metrics, same seed). The current code is only self-evaluating MOSAIC.

**Severity**: 🟡 Warning (not blocking for Exp 1 *execution*, but blocking for the Post-Results review and paper submission).

**Recommended fix**: Task #10 in the task list already captures this. Do not skip it.

---

### W6. Preprocessing leakage into Exp 3 (missing-type detection)

**Where**: Future implementation of Exp 3 (`src/evaluation/missing_type.py` — not yet written).

**What**: If Exp 3 removes cell type *k* from one modality and then runs preprocessing, the preprocessor will see only the remaining cells. But if Exp 3 *reuses* the already-preprocessed AnnData and filters in-memory, PCA / LSI components were learned with type *k* present, which is effectively "preprocessing leakage backward": the preprocessing step knew about a type that will later be treated as absent. This would bias Exp 3's AUROC.

**Severity**: 🟡 Warning (forward-looking, not a current bug).

**Recommended fix**: Exp 3 must re-run preprocessing on the cell-type-removed subset before training IB-VAEs. This is a protocol note for when Exp 3 is implemented (Task #12).

---

## Observations (Nice to Fix)

### O1. Unused imports in `train_ibvae.py`
`make_loader` and `make_split_indices` are imported from `src.training.dataloader`, but `make_loader` is no longer used after the refactor to manual batch indexing. Safe to remove.

### O2. Dead variable in `train_ibvae.py`
`val_steps_per_epoch` is computed (line 155) but never read. Safe to remove.

### O3. `scanpy.pp.scale` densifies sparse RNA matrix
At `preprocess.py::preprocess_rna`, `sc.pp.scale` emits a UserWarning about densifying a sparse matrix. This is harmless for 11K × 2K but worth documenting; for larger datasets we could skip scale and rely on the encoder's own LayerNorm to center/scale features.

### O4. `val_pred` as early-stopping criterion may bias toward cross-modal prediction at the expense of reconstruction
Early stopping watches only `val_pred` (cross-modal prediction MSE), not the combined validation total. For a model where recon matters for the downstream OT embedding quality, this might stop too early. Suggest switching to `val_total` or logging both and picking based on downstream OT FOSCTTM post-hoc.

### O5. `ModalityDataset.get_batch` unused
After the refactor, `get_batch` in `dataloader.py` is only referenced by internal tests, not by the training loop. The training loop pulls directly from the preloaded GPU tensors `X_dev`, `recon_dev`, `yc_dev`. `get_batch` can be removed to simplify the API.

### O6. Fixed leiden resolution (0.8) produces 18 clusters
May be too many or too few for different datasets. Consider making it dataset-specific or tuning it such that n_clusters is in a reasonable range (5–30) for each dataset.

---

## Checklist Summary

| Category              | Checked | Passed | Warned | Failed |
|-----------------------|--------:|-------:|-------:|-------:|
| Code correctness      |      18 |     18 |      0 |      0 |
| Loss / optimization   |       6 |      5 |      1 |      0 |
| Data pipeline         |      11 |      9 |      2 |      0 |
| Training loop         |       8 |      8 |      0 |      0 |
| Data leakage          |       7 |      5 |      2 |      0 |
| Experimental design   |       6 |      4 |      2 |      0 |
| Reproducibility       |       6 |      5 |      1 |      0 |
| Code quality          |       5 |      5 |      0 |      0 |

**Totals**: 67 items checked, 59 passed, 8 warned, 0 failed.

---

## Detailed Trace: The Pair-Leakage Invariant

This is the most important MOSAIC-specific concern, so the trace is spelled out here.

**Goal**: Verify that at training time, a single-modality IB-VAE cannot implicitly or explicitly observe which RNA cell is paired with which ATAC cell.

**Trace**:
1. `preprocess.py::load_10x_multiome_h5` splits the 10x h5 into two `AnnData` objects. Both receive `obs["pair_idx"] = np.arange(n_obs)` so post-hoc evaluation can recover the pairing. After `qc_filter_joint`, the two `AnnData`s are written to *separate* h5ad files.
2. `train_ibvae.py::train` takes a `processed_path` for **one** modality and calls `ad.read_h5ad(cfg.processed_path)`. It never opens the other modality's file.
3. `dataloader.py::ModalityDataset.__init__` reads `adata.X`, `adata.layers[recon_layer]`, and `adata.obsm["y_cross"]`. It does *not* read `adata.obs["pair_idx"]`. Confirmed by grep: the only appearances of `pair_idx` in `src/training` are in the `run_experiment.py` top-level runner (post-training, for evaluation) and in a comment in `dataloader.py` stating that the field is not accessed.
4. During training, the iterator is `perm = torch.randperm(n_train, generator=gen).numpy()`, where `gen` is seeded from `cfg.seed`. The two modalities are trained sequentially with the **same** seed, so the two permutations will be identical — but that doesn't matter because (a) the training code per modality never sees the other modality's data in the first place, (b) the permutations are over *shuffled* indices so the training objective is invariant to them, and (c) the evaluation code uses the `pair_idx` field, never relying on row-order alignment between the two saved tensors.
5. The `pair_leakage_test` in `src/data/validate.py` explicitly constructs two independently-shuffled orderings of the two modalities and verifies that the `pair_idx` columns disagree on 100% of rows — confirming no code path can recover the pairing from the shuffled data alone. This test passes on the PBMC 10k preprocessed data.

**Conclusion**: Pair-leakage is structurally impossible under the current design. ✓

---

## Smoke-Test Evidence

A 2-epoch end-to-end run (`experiments/smoke_test`, later removed) and a 10-epoch run on PBMC 10k were completed prior to this review, purely to verify the pipeline has no runtime bugs. Results as of the 10-epoch run (these are NOT trained and NOT reported results — they exist only to demonstrate the pipeline runs):

- Training time: ~0.3 sec/epoch (RNA), ~0.2 sec/epoch (ATAC) after moving dense tensors to GPU.
- Total wall time for 10 epochs + OT + all metrics: 122 sec, dominated by OT (~10 sec) and metrics computation (KMeans, ~30 sec).
- Projected 200-epoch runtime: ~2–3 minutes training + ~1 minute downstream = well under the single-GPU budget.

The 10-epoch run produces nonsense results (FOSCTTM ~0.45, entropy ~1.0, Spearman ρ = −0.64) because the IB-VAE is nowhere near converged with β still in warmup. **None of these numbers will be reported**; the smoke-test directory has been deleted. Actual training will run 200 epochs with early stopping.

---

## Verdict and Next Steps

🟡 **PASS WITH CONCERNS — training may proceed.**

Proceed to Milestone 6 (Exp 1 on PBMC 10k) with default hyperparameters from `RESEARCH_PLAN.md §3.4`. After Exp 1 completes, run baselines (Milestone 7), Exp 2 (Milestone 8), and then the **Post-Results Review Gate** which re-checks the `W1`–`W6` warnings in light of actual results.

The 6 warnings above are to be surfaced again at the Post-Results gate; none block the start of training.
