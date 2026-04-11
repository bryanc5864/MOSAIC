# Post-Results Review Report — MOSAIC

**Review Mode**: Post-Results
**Date**: 2026-04-11 04:40
**Reviewer**: Autonomous Review Skill
**Project root**: `C:\Users\Maozer\projects\MOSAIC`

---

## Summary

**Overall Status**: 🟡 **PASS WITH CONCERNS** — paper drafting may proceed. No fabrication, no data leakage in the results, and every numerical claim in `RESULTS.md` traces to a specific JSON file on disk. Findings below are cosmetic rounding, documentation inconsistencies, and pre-existing limitations noted in the Pre-Training review; none call the scientific conclusions into question.

The core scientific story is sound:
- FOSCTTM, label transfer, ARI, and cluster-resolved entropy metrics on both PBMC 10k and Brain 5k are all reproduced exactly by their JSON source files.
- Cluster-resolved entropy argmax accuracy (98.76% PBMC, 96.15% Brain) is computed from `alignment_plan_subsample.npy` × `cell_type` labels and re-derived by `scripts/cluster_resolved_entropy.py` — the script is clean, no ground-truth pairing leaks into the computation.
- Missing-type AUROC mean 0.96 on PBMC and 0.96 on Brain is computed by `scripts/missing_type_exp.py`, which runs a legitimate leave-one-cluster-out loop on the already-trained embeddings.
- Baseline comparison against SCOT is fair (same preprocessed AnnData, same seeded 3000-cell subsample, same metric functions).

---

## Integrity Verdict

- **Fabrication detected**: No
- **Data leakage detected**: No — the pair-leakage invariant from the Pre-Training review still holds. `grep -rn pair_idx src/training/` confirms `pair_idx` is only read in `run_experiment.py` (post-training evaluation), never in `dataloader.py` or `train_ibvae.py`.
- **All results traceable to code**: Yes — see §6.1 below for the line-by-line audit.
- **Statistical claims valid**: Partially — the magnitudes are real, but the paper lacks confidence intervals or seed variance (single seed run only).

---

## Critical Issues (Must Fix)

**None.** No blockers for paper drafting.

---

## Warnings (Should Fix)

### W1. SCOT baseline numbers in `RESULTS.md` and `fig4_baselines_pbmc.png` are rounded incorrectly

**Where**: `RESULTS.md` baseline comparison table; `scripts/make_main_figures.py::figure_4_baselines` (hardcoded values).

**What**:
- SCOT `label_transfer_rna_to_atac` actual value in `experiments/baselines_pbmc10k_multiome/baseline_results.json` = **0.36066666…** (n=3000, 1082 correct hits ÷ 3000). `RESULTS.md` and the hardcoded figure value both say **0.3610**, which is 0.0003 higher than the true value.
- SCOT `label_transfer_atac_to_rna` actual value = **0.38633333…**. `RESULTS.md` says **0.3860**, which is 0.0003 lower than the true value.
- SCOT `foscttm_mean` actual = **0.24489268…**. `RESULTS.md` says **0.2449** — correctly rounded ✓.
- SCOT `joint_clustering_ari` actual = **0.34520127…**. `RESULTS.md` says **0.3452** ✓.

**Severity**: 🟡 Warning (not material to any conclusion). The SCOT vs MOSAIC gap is large — ~0.3 in label transfer, ~0.3 in ARI — and a 0.0003 correction in either direction leaves the comparison unchanged. But the rounding is inconsistent: the first value (0.3607 → 0.3610) rounds UP, the second (0.3863 → 0.3860) rounds DOWN. An attentive reviewer will flag this.

**Fix**: Replace with the exact 4-decimal values: 0.3607 and 0.3863 (or 3-decimal 0.361 and 0.386). Regenerate Fig 4.

---

### W2. "87% better ARI" claim is off by 1 percentage point

**Where**: `paper/results.md` line 29.

**What**: The paper claims "MOSAIC outperforms SCOT by … 87% better ARI." The actual improvement is `(0.6515 − 0.3452) / 0.3452 = 0.8874 = 88.7%`. The paper should say 89%, not 87%.

**Severity**: 🟡 Minor.

**Fix**: Change to "89% better ARI".

---

### W3. `TRAINING_LOG.md` is incomplete

**Where**: `TRAINING_LOG.md`.

**What**: The training log has detailed entries only for Run 001 and Run 005 (= `exp001_pbmc_final`). Runs 002, 003, 004, and the failed similarity-Procrustes experiment (`run005_simproc_medianeps`) are only described inline in the Run 005 entry's preamble, not as separate log entries with their own hyperparameters and metrics. The raw `experiments/run00N_*/` directories exist and contain their own `results.json`, so the data *is* on disk, but the log is out of sync with what's actually in `experiments/`.

**Additionally**: the log entry for Run 005 says "Artifacts: `experiments/exp001_pbmc_final/`" and the "Prior exploratory runs" section says these artifacts were "deleted to save disk". This is **inaccurate** — `run001_failed_lambda1/`, `run002_lambda10/`, `run003_centroid_targets/`, `run004_procrustes/`, and `run005_simproc_medianeps/` all still exist with full checkpoints, embeddings, and results.json files. Nothing has been deleted.

**Severity**: 🟡 Warning — documentation accuracy issue, not a fabrication issue. A reviewer who wanted to audit intermediate runs would be able to find the data; they'd just be confused by the log's claim that it was deleted.

**Fix**: Either (a) add formal TRAINING_LOG entries for runs 002, 003, 004, and run005_simproc_medianeps with metrics pulled from their respective JSON files, or (b) update the preamble to say "artifacts still on disk under `experiments/run00N_*/` for reproducibility".

---

### W4. Run numbering inconsistency between `TRAINING_LOG.md` and `RESULTS.md`

**Where**: `TRAINING_LOG.md` (Run 005 = `exp001_pbmc_final`) vs `RESULTS.md` running commentary (step 5 = `run005_simproc_medianeps`, step 6 = `exp001_pbmc_final`).

**What**: The two documents disagree on what "Run 005" means. `TRAINING_LOG.md` gives the label "Run 005" to the canonical final run (`exp001_pbmc_final`), while `RESULTS.md`'s running commentary treats the similarity-Procrustes experiment as step 5 and `exp001_pbmc_final` as step 6.

**Severity**: 🟡 Warning — confusing to a reader cross-referencing the two.

**Fix**: Renumber consistently. The most intuitive is to keep `RESULTS.md`'s numbering (5 = simproc attempt, 6 = canonical) and update `TRAINING_LOG.md`'s main entry header to "Run 006".

---

### W5. Single-seed reporting — no statistical uncertainty on metrics

**Where**: All primary metrics in `RESULTS.md` and `paper/results.md`.

**What**: Every number reported is from a single training run at `seed=0`. The paper's comparative claims ("21% better FOSCTTM, 84% better LT, 89% better ARI than SCOT") have no associated confidence interval or seed-variance estimate. For a paper submitted to RECOMB / Cell Systems / Genome Research, reviewers will expect at least 3-seed variance.

**Severity**: 🟡 Warning (expected-for-submission, not for this review gate). The magnitudes are large enough that the comparison will hold across seeds, but this needs to be shown empirically.

**Fix**: Re-run the canonical pipeline at seeds ∈ {0, 1, 2}. Running each takes ~2 minutes, total ~6 minutes. Report mean ± std for every primary metric. This is a 10-minute job that materially strengthens the paper.

---

### W6. `NN_on_IB` baseline wall time is misleading

**Where**: `RESULTS.md` baseline comparison table.

**What**: `NN_on_IB` is reported with a wall time of 0.8 seconds — this is the inference-only time after the IB-VAE was already trained. MOSAIC's 118-second wall time *includes* training. A reader comparing the two might conclude that `NN_on_IB` is 100× faster than MOSAIC, which is false: both depend on the same IB-VAE training and therefore have the same total runtime; `NN_on_IB` just skips the OT step at inference.

**Severity**: 🟡 Warning.

**Fix**: Add a footnote: "NN_on_IB reuses MOSAIC's trained IB-VAE; its wall time excludes the ~95 seconds of IB-VAE training that MOSAIC's wall time includes. Fair total: ~95 s + 0.8 s = ~96 s."

---

### W7. Figure 4 uses hardcoded values rather than reading from JSONs

**Where**: `scripts/make_main_figures.py::figure_4_baselines`.

**What**: The figure hardcodes the metric values in a Python list instead of loading them from `experiments/baselines_pbmc10k_multiome/{baseline_results.json,simple_baseline_results.json}`. If the baseline results are re-run with new seeds or different subsample sizes, the figure will silently become inconsistent with the source-of-truth JSONs.

**Severity**: 🟡 Warning (technical-debt).

**Fix**: Load from JSON in the figure-generation script. Propagates W1 correctness too.

---

### W8. Preprocessing leakage (carried over from Pre-Training review)

**Where**: `src/data/preprocess.py` — `sc.pp.scale`, `sc.tl.pca`, `sc.tl.leiden`, `TruncatedSVD` all fit statistics on the full dataset (including the validation subset).

**What**: This was flagged in the Pre-Training review (W1 there). Does not affect the primary metrics because they use the full dataset for both alignment and evaluation (there is no truly held-out test split). However, the leiden cluster labels that serve as the ground truth for Exp 3's "is this cell in the removed cluster" AUROC computation were learned on the FULL dataset before the cluster was removed. This is a weak form of label leakage: the cluster identity we're trying to detect is partially defined by the cells we're trying to hide.

**Severity**: 🟡 Warning — methodology nit. The AUROC of 0.96 is still a meaningful demonstration because the leave-out happens at the ATAC side while cluster labels come from the RNA side (which still has the cluster), but a rigorous reviewer would ask for the experiment to be repeated with labels computed on the post-removal data only.

**Fix**: Defer to final paper prep. Noted as a limitation in `paper/discussion.md`.

---

## Observations (Nice to Fix)

### O1. Ablation sweep running suggests β=0.001 is better than the reported β=0.01

**Where**: Live output from `scripts/ablation_sweep.py` — `beta_0.001` variant.

**What**: The ablation sweep's `beta_0.001` variant gives (at the time of this review, for the in-progress sweep):
- FOSCTTM 0.1084 (vs base 0.1781)
- LT RNA→ATAC 0.9758 (vs base 0.7249)
- LT ATAC→RNA 0.9302 (vs base 0.6703)
- ARI 0.6913 (vs base 0.6726)

These are substantially BETTER than the reported main result (`exp001_pbmc_final` with β=0.01). The reported main result is therefore *conservative* — a lower β would likely give better numbers. This is a favorable integrity observation (the paper is understating MOSAIC's ceiling), not a concern.

**Recommendation**: For the final paper, either (a) re-run the canonical configuration at β=0.001 and update the headline numbers, or (b) keep β=0.01 and explicitly note in the ablation section that the reported results are conservative.

### O2. Running commentary in `RESULTS.md` references metrics from runs without log entries

**Where**: `RESULTS.md` lines 97–103 (Running Commentary).

**What**: Claims like "Run 002 — ... Cross-modal MSE on the training set dropped to 0.06, but val MSE stayed at 1.2" are verifiable from `experiments/run002_lambda10/results.json` (best_val 1.10 for RNA, 0.97 for ATAC — close to "1.2" but not exactly). The "train MSE 0.06" claim is not in the JSON but IS in the saved `train_log_{rna,atac}.json` epoch histories. So the claims are defensible but would be easier to audit if the numbers were in a structured JSON rather than in prose.

**Severity**: Nice-to-fix.

### O3. Small inconsistencies in Brain ARI (0.367)

**Where**: `paper/results.md` line 14; `RESULTS.md` line 15.

**What**: Brain ARI 0.3665 rounds to **0.367** in paper (3-decimal) but **0.3665** in RESULTS (4-decimal). Both are correct; the different rounding conventions between documents are cosmetic.

### O4. Ablation sweep `base` variant != main experiment

**Where**: `experiments/ablation_pbmc10k_multiome_base/` — uses `epochs=150, ot_subsample=3000` vs main experiment's `epochs=200, ot_subsample=5000`.

**What**: The ablation sweep's "base" configuration is NOT identical to `exp001_pbmc_final`. It trains for 150 epochs instead of 200 and uses a 3000-cell OT subsample instead of 5000. As a consequence, the `base` variant's FOSCTTM of 0.1781 is slightly BETTER than the main result's 0.1947. This means the ablation sweep's results are internally consistent (all variants at 150 epochs) but the `base` variant shouldn't be cited as "the main result's config". It's an *approximation* of the main config at smaller scale.

**Severity**: Nice-to-fix. The ablation's "base" is close enough to the main config that its deltas with the other variants are meaningful. But the paper's ablation section should note that the sweep was run at 150 epochs.

---

## Code-to-Result Traceability Audit (§6.1 detail)

I walked through RESULTS.md line by line. For each quantitative claim, I identified the source JSON and verified the number.

### Exp 1 full-dataset metrics (RESULTS.md §"Exp 1 — Alignment quality")

| Claim | Source | Value in source | Match |
|---|---|---|---|
| PBMC FOSCTTM 0.1947 | `exp001_pbmc_final/results.json` `.metrics.foscttm.foscttm_mean` | 0.19469161058698542 | ✓ |
| PBMC A→B 0.2250 | `.metrics.foscttm.foscttm_a_to_b` | 0.22495509975043854 | ✓ |
| PBMC B→A 0.1644 | `.metrics.foscttm.foscttm_b_to_a` | 0.16442812142353233 | ✓ |
| PBMC LT RNA→ATAC 0.6849 | `.metrics.label_transfer_rna_to_atac` | 0.6848624259046271 | ✓ |
| PBMC LT ATAC→RNA 0.7386 | `.metrics.label_transfer_atac_to_rna` | 0.738564982747943 | ✓ |
| PBMC ARI 0.6802 | `.metrics.joint_clustering_ari` | 0.6801550553548548 | ✓ |
| PBMC Mean H 0.7729 | `.alignment.entropy_mean` | 0.7729314258418157 | ✓ |
| PBMC H std 0.1138 | `.alignment.entropy_std` | 0.11382731382655603 | ✓ |
| PBMC Spearman ρ −0.6481 | `.metrics.entropy_error_corr.spearman_rho` | −0.6481192722287707 | ✓ |
| PBMC wall 117.6s | `.wall_time_sec` | 117.56026196479797 | ✓ |
| Brain FOSCTTM 0.1270 | `exp001_brain_final/results.json` | 0.12700920273046654 | ✓ |
| Brain LT A→B 0.8718 | `.metrics.label_transfer_rna_to_atac` | 0.8717722357095564 | ✓ |
| Brain LT B→A 0.7484 | `.metrics.label_transfer_atac_to_rna` | 0.7483999117192672 | ✓ |
| Brain ARI 0.3665 | `.metrics.joint_clustering_ari` | 0.36650098887905014 | ✓ |
| Brain Mean H 0.7423 | `.alignment.entropy_mean` | 0.7423360051723007 | ✓ |

### Exp 1 baseline comparison (3000-cell subsample)

| Claim | Source | Value in source | Match |
|---|---|---|---|
| MOSAIC FOSCTTM 0.1941 | `simple_baseline_results.json::nn_on_ib.metrics.foscttm_mean` (same config as MOSAIC on subsample) | 0.1940869178615094 | ✓ |
| MOSAIC LT 0.664 | `.nn_on_ib.metrics.label_transfer_rna_to_atac` | 0.664 | ✓ |
| SCOT FOSCTTM 0.2449 | `baseline_results.json::scot.metrics.foscttm.foscttm_mean` | 0.24489268645103923 | ✓ |
| SCOT LT A→B 0.3610 | `.scot.metrics.label_transfer_rna_to_atac` | 0.3606666666666667 | ⚠ (W1) |
| SCOT LT B→A 0.3860 | `.scot.metrics.label_transfer_atac_to_rna` | 0.3863333333333333 | ⚠ (W1) |
| SCOT ARI 0.3452 | `.scot.metrics.joint_clustering_ari` | 0.3452012709978835 | ✓ |
| RawOT FOSCTTM 0.3283 | `simple_baseline_results.json::raw_ot.metrics.foscttm_mean` | 0.3282678114927198 | ✓ |
| RawOT LT A→B 0.1317 | `.raw_ot.metrics.label_transfer_rna_to_atac` | 0.13166666666666665 | ✓ |
| RawOT ARI 0.0931 | `.raw_ot.metrics.joint_clustering_ari` | 0.09305441161925457 | ✓ |

### Exp 2 cluster-resolved entropy

| Claim | Source | Value in source | Match |
|---|---|---|---|
| PBMC argmax cluster acc 98.76% | `exp001_pbmc_final/cluster_entropy_analysis.json::cluster_level_entropy.argmax_cluster_accuracy` | 0.9876 | ✓ |
| PBMC mean H_cluster 0.153 | `.cluster_level_entropy.mean` | 0.15303407609462738 | ✓ |
| PBMC AUROC 0.894 | `.cluster_level_entropy.auroc_entropy_vs_wrong_cluster` | 0.8936914514169247 | ✓ |
| PBMC n_wrong 62 | `.cluster_level_entropy.n_wrong_cluster_cells` | 62 | ✓ |
| Brain argmax cluster acc 96.15% | `exp001_brain_final/cluster_entropy_analysis.json::cluster_level_entropy.argmax_cluster_accuracy` | 0.9615 | ✓ |
| Brain mean H_cluster 0.262 | `.cluster_level_entropy.mean` | 0.2619847059249878 | ✓ |
| Brain AUROC 0.861 | `.cluster_level_entropy.auroc_entropy_vs_wrong_cluster` | 0.8609991152892869 | ✓ |
| Brain n_wrong 154 | `.cluster_level_entropy.n_wrong_cluster_cells` | 154 | ✓ |

### Exp 3 missing type detection

| Claim | Source | Value in source | Match |
|---|---|---|---|
| PBMC mean AUROC 0.960 | `exp001_pbmc_final/exp003_missing_type.json::mean_auroc` | 0.9595604072449682 | ✓ |
| PBMC median AUROC 0.987 | `.median_auroc` | 0.9871704831110291 | ✓ |
| PBMC AUROC range 0.67–1.00 | `.min_auroc`, `.max_auroc` | 0.6707391321995373, 0.9996918980600988 | ✓ |
| Per-cluster table in paper/results.md (all 12 rows) | `.per_cluster[*].target_cluster, .auroc_cluster_entropy, .mean_entropy_target, .mean_entropy_other` | all match | ✓ |
| Brain mean AUROC 0.959 | `exp001_brain_final/exp003_missing_type.json::mean_auroc` | 0.9594185037127164 | ✓ |
| Brain median AUROC 0.964 | `.median_auroc` | 0.9636572080156647 | ✓ |
| Brain range 0.915–0.997 | `.min_auroc`, `.max_auroc` | 0.9146595088975047, 0.9973130377823427 | ✓ |

**Conclusion**: Every number in RESULTS.md and the paper sections that I could check traces to a JSON file on disk with the exact (or correctly-rounded) value. Only the SCOT LT rounding (W1) is off by 0.0003.

---

## Figure Audit

| Figure | Source data | Accurate? |
|---|---|---|
| fig1_aligned_latent_{pbmc,brain} | `experiments/exp001_*/z_rna.npy, z_atac_aligned.npy`, `pbmc/brain_rna.h5ad::obs.cell_type` | ✓ — uses PCA projection of the actual trained latents |
| fig2_entropy_comparison_{pbmc,brain} | `experiments/exp001_*/alignment_plan_subsample.npy, alignment_entropy_subsample.npy, z_rna.npy, z_atac_aligned.npy` | ✓ — recomputes H_cluster inline from the plan and cluster labels |
| fig3_missing_type_auroc | `experiments/exp001_*/exp003_missing_type.json` | ✓ — reads directly |
| fig4_baselines_pbmc | HARDCODED values | ⚠ W7 — should load from JSON; also inherits W1 |

---

## TRAINING_LOG.md Cross-Check

- **Run 001** entry is accurate vs the `run001_failed_lambda1/results.json` values (FOSCTTM 0.4608, LT 0.147/0.136, ARI 0.390 — exact matches in RESULTS.md commentary and TRAINING_LOG).
- **Run 005** entry is accurate vs `exp001_pbmc_final/results.json` — all metrics match to 4 decimals.
- **Runs 002, 003, 004** — not logged. Need W3 fix.
- **Timestamps are plausible**: PBMC 117.6s wall clock for 200 epochs × 2 modalities + OT + metrics is reasonable on RTX 3080; Brain 49.8s for 4531-cell version is plausible (about 2.4× faster than PBMC because of ~2.5× fewer cells). These are consistent with the epoch-level timings I can see streaming from the in-progress ablation sweep (~0.3s/epoch on the same hardware).
- **No fabricated timing**: every `wall_time_sec` I checked came from a real `time.time()` delta in the training script.

---

## Checklist Summary

| Category | Items Checked | Passed | Warned | Failed |
|---|---:|---:|---:|---:|
| Result authenticity | 32 | 30 | 2 (W1) | 0 |
| Code-to-result traceability | 32 | 32 | 0 | 0 |
| Statistical validity | 4 | 1 | 3 (W5, W6, W2) | 0 |
| TRAINING_LOG consistency | 6 | 3 | 3 (W3, W4, O2) | 0 |
| Figure accuracy | 7 | 6 | 1 (W7) | 0 |
| Paper claims match results | 12 | 11 | 1 (W2) | 0 |
| Data leakage (re-check) | 5 | 4 | 1 (W8, pre-existing) | 0 |

**Totals**: 98 checked, 87 passed, 11 warned, 0 failed.

---

## Verdict

🟡 **PASS WITH CONCERNS — paper drafting may proceed.**

Zero fabrication, zero critical data leakage, every number in RESULTS.md traces to a JSON file on disk. The warnings are cosmetic rounding (W1, W2), documentation gaps in TRAINING_LOG.md (W3, W4, W6, O2), single-seed reporting (W5 — significant but known), a figure using hardcoded values (W7), and a pre-existing methodology nit from the Pre-Training review (W8).

Of the warnings, **W5 (single-seed reporting)** is the only one that a RECOMB/Cell Systems reviewer would likely demand be fixed before publication. It's a 10-minute fix: re-run the canonical pipeline at seeds 0, 1, 2 and update the tables to show mean ± std. Everything else is cosmetic or deferred-to-final-prep.

**Recommended next steps before proceeding to the stretch experiments and paper finalization**:
1. (15 min) Re-run seeds 1, 2 on PBMC and Brain; update RESULTS.md with mean ± std.
2. (5 min) Fix W1 SCOT rounding in RESULTS.md and fig4_baselines_pbmc.
3. (5 min) Fix W2 (89% vs 87% on ARI).
4. (10 min) Fix W3 / W4 TRAINING_LOG.md — add entries for runs 002/003/004 and renumber.
5. (5 min) Fix W6 NN_on_IB wall time footnote.

Total time for all warning fixes: ~40 minutes. None are blockers for writing Methods, Results, Discussion, or the convergence theorem appendix — those can proceed in parallel.

**Post-Results Review Gate: ✅ PASSED (with concerns noted).**
