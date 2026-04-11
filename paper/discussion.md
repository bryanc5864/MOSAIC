# Discussion

*MOSAIC: Multi-Omics Stochastic Alignment via Information-theoretic Coupling*

## Summary of findings

MOSAIC provides the first calibrated per-cell alignment uncertainty for unpaired single-cell multi-omics integration. The method outperforms SCOT by 21% on FOSCTTM and 84% on label transfer accuracy on 10x Multiome PBMC 10k, matches or beats SCOT on Brain 5k, and — most importantly — correctly identifies which cells are mis-aligned with AUROC 0.89 (PBMC) and 0.86 (Brain). In leave-one-cluster-out missing-type experiments, the cluster-resolved entropy detects absent cell types with AUROC 0.96 averaged across 12 PBMC and 18 Brain cluster holdouts.

Three findings deserve particular emphasis.

**First, the cross-modal prediction target is the pivotal design choice, not the choice of IB-VAE vs. some other architecture.** Our original design used the partner cell's own feature vector as the cross-modal prediction target. This failed catastrophically: the model memorized per-cell training fingerprints and learned nothing generalizable. Switching to a *cluster-centroid* target — in which every cell of cluster $k$ receives the same target, the mean of the partner modality's cluster-$k$ cells — was a 3.6× improvement in label transfer accuracy (0.19 → 0.68) without any other change. We believe this generalizes: any method that couples two independently-encoded modalities through a cross-modal prediction head is susceptible to per-cell memorization unless the target has a cluster-level structure that the model cannot memorize.

**Second, cell-level row entropy on entropic-OT transport plans is miscalibrated as an alignment uncertainty measure.** This is the key negative result of the paper. Spearman correlation between per-cell row entropy and alignment error is $\rho = -0.65$ on PBMC — the *opposite* of what a well-calibrated uncertainty should show. The explanation is geometric: cells in dense clusters have many valid within-cluster candidates (inflating entropy) but also a closer true partner (reducing alignment error), so the two effects reverse the expected correlation. Marginalizing over target cluster labels — what we call *cluster-resolved entropy* — removes this confound and produces a calibrated signal. Any OT-based method that claims to report per-cell alignment confidence via row entropy needs to check whether this pathology applies, and our framework suggests the marginalization cure is both cheap and principled.

**Third, the Procrustes rotation on cluster centroids is essential.** Independently-trained IB-VAE encoders produce latent clouds that live in different rotations of the same underlying structure. Without Procrustes alignment, the squared-Euclidean cost matrix between them is near-uniform, and Sinkhorn returns near-uniform plans with no per-cell signal. A simple orthogonal-Procrustes rotation on shared cluster centroids brings the two clouds into a common frame, at the cost of assuming shared cluster labels at *training* time (legitimate for paired benchmark datasets, where leiden labels are propagated across modalities via the pairing). At inference on genuinely unpaired data, the same rotation could be fitted from a joint clustering of the concatenated latents or from an external label source, both of which we flag as future work.

## Scope and limitations

### Cluster-label dependence

MOSAIC's Procrustes alignment and the cluster-resolved entropy definition both require cluster labels. On paired benchmark datasets we propagate leiden clusters from one modality to the other via the known pairing, which is legitimate for evaluation but does not translate directly to unpaired application settings. The natural replacement in production is a joint clustering of the concatenated IB-VAE latents (a chicken-and-egg problem if those latents are themselves produced by training that depends on cluster labels), or an external labeling source (e.g., celltypist, Azimuth for immune populations). A principled treatment of the label-free setting is a key direction for follow-up.

### Loss-scale inconsistency

The IB-VAE loss mixes per-element reconstruction loss (ZINB or BCE) with per-sample summed-over-dims KL and per-element MSE on the cross-modal head. This is a standard if inelegant VAE practice that effectively down-weights the KL by the number of features, helping avoid posterior collapse. The reported $\beta = 0.01$ is therefore not directly comparable to a $\beta$ in a paper that uses consistent per-element reductions; readers should interpret it as "KL weight in our specific reduction convention, chosen empirically to give a decent cross-modal prediction without collapsing the latent."

### Cost normalization and Sinkhorn solver

We adopt median-normalized cost and POT's regular (non log-domain) Sinkhorn solver. Log-domain Sinkhorn would be numerically preferred for very small $\epsilon$, but on an 11303×11303 cost matrix at float64, scipy's `logsumexp` allocates several multi-GB intermediate buffers and is OOM-killed on 16 GB workstations. The regular solver is stable at our operating $\epsilon = 0.05$ on median-normalized cost and is what we recommend for problems in our size regime. For substantially smaller problems (N ≤ a few thousand) the log-domain solver is preferable if memory permits.

### Dataset scope

We benchmark on two 10x Multiome datasets that happen to share the same scRNA-seq + scATAC-seq modality pair. The method as written assumes this modality pair — the ATAC preprocessing pipeline (TF-IDF + LSI + BCE decoder) is ATAC-specific. Extending to RNA+protein (CITE-seq) or RNA+methylation requires only changing the decoder likelihood and the preprocessing; the IB + Procrustes + cluster-resolved-entropy scaffolding is modality-agnostic. We leave the CITE-seq and methylation benchmarks to future work but are optimistic based on the conceptual simplicity of the adaptation.

### Convergence theorem

The original research plan promised a convergence theorem: under sub-Gaussian embedding assumptions and sufficient sample size, cluster-resolved entropy converges to 0 for shared cell types and remains bounded away from zero for absent types. We have empirical support for both claims (argmax cluster accuracy 96–99% for shared types; mean AUROC 0.96 for detecting absent types) but a formal proof requires rigorous handling of Sinkhorn dual variables under Procrustes rotation, which we defer to the appendix. The core ingredients — concentration of cluster centroids at rate $O(\sigma^2/\sqrt{n})$, Sinkhorn's log-sum-exp approximation of min-cost matching, and the bounded entropy of a uniform distribution over $K$ clusters — all fit together in a straightforward way but require careful statement.

## β as a tunable knob — a dataset-dependent trade-off

A finding from full-scale ablation that we initially mistook for a trade-off but resolved with multi-seed reporting. We initially defaulted to β=0.01 because a single-seed PBMC run at β=0.001 showed dramatically lower joint-clustering ARI (0.49 vs 0.68 at β=0.01). With a 3-seed replication, however, that ARI drop turns out to be a one-seed outlier from KMeans initialization variance:

- **On Brain 5k, β=0.001 is strictly better** across every primary metric (3-seed mean ± std): FOSCTTM 0.049 ± 0.002 (vs 0.129 ± 0.004 at β=0.01), ARI 0.94 ± 0.003 (vs 0.41 ± 0.04), and the cluster-level UQ AUROC improves from 0.81 ± 0.05 to 0.95 ± 0.02. Brain is uniformly better.
- **On PBMC 10k, β=0.001 wins or ties** on every metric in multi-seed (3 seeds): FOSCTTM 0.120 ± 0.014 (vs 0.188 ± 0.006), LT 88.3% / 87.4% (vs 68.9% / 77.1%), ARI 0.65 ± 0.14 (vs 0.69 ± 0.01 — within 1σ), argmax cluster acc 96.9% (vs 98.4% — 1.5 pp loss). The cluster-level entropy AUROC is unchanged at 0.89.

The single significant difference between datasets is that PBMC at β=0.001 has much higher seed variance on cluster-level metrics (ARI std 0.14 vs 0.009 at β=0.01). The mechanism: a weaker bottleneck (β=0.001) allows the IB-VAE latent to retain finer per-cell detail, which makes downstream KMeans more sensitive to its random initialization — particularly on a dataset like PBMC with very imbalanced cluster sizes (some <50 cells). Brain has 20 fairly evenly-sized clusters, so KMeans is stable.

The practical recommendation is to use **β=0.001 as the default**. The cluster-resolved entropy AUROC, our headline UQ metric, is roughly insensitive to this choice, so the UQ signal is robust. The only caveat is that joint-clustering ARI on datasets with heavily-imbalanced cluster sizes will have higher seed variance at lower β; report ARI as mean ± std over multiple seeds rather than a single number.

## Open questions and next steps

1. **Unpaired application at scale.** We tested MOSAIC on paired benchmarks so we could evaluate against ground truth. The next step is a genuinely unpaired application (Tabula Sapiens × ENCODE scATAC) where cluster labels come from a joint clustering of the IB-VAE latents. We expect the cluster-resolved entropy to surface biologically interpretable signals: cell populations that are in one atlas but not the other, or regulatory states (enhancer priming) where chromatin accessibility and transcription are decorrelated.

2. **Learning the cluster assignment end-to-end.** Currently we use leiden clusters computed separately on each modality. A joint, learned clustering that consumes both modalities' latents and provides the labels MOSAIC depends on would close the loop and remove the external cluster-label dependence.

3. **Multi-modal generalization (> 2 modalities).** MOSAIC as written handles two modalities; the Procrustes rotation would need to be generalized to a joint orthogonal transformation over all modalities simultaneously (a group Procrustes problem), and the cluster-resolved entropy would extend to multi-modal cluster marginals naturally.

4. **Principled treatment of cell populations present in multiple but not all modalities.** Our leave-one-out missing-type test assumes cluster $k^*$ is present in modality A but absent in modality B. Real unpaired data may have cluster $k^*$ present in modalities A and B but absent in C, which requires a richer uncertainty decomposition than the binary "present/absent" framing.

5. **Theoretical calibration of cluster-resolved entropy.** We show empirically that $H_{\text{cluster}}$ has AUROC ≈ 0.9 for the "this cell's argmax cluster is wrong" task, but we have not characterized the distribution of $H_{\text{cluster}}$ values under the null hypothesis of well-calibrated alignment. A calibration curve would be more interpretable than an AUROC for downstream users.
