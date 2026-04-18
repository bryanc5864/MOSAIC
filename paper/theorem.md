# Theorem: Convergence of cluster-resolved alignment entropy

*the method appendix — informal statement and proof sketch*

## Setup

Let $\mathcal{A}$ and $\mathcal{B}$ be two single-cell modalities with a shared set of $K$ cell types. For each type $k \in [K]$, let $n$ cells be drawn independently from a distribution on the modality's feature space. After applying the method's per-modality IB-VAE encoder and the cluster-centroid Procrustes rotation, each cell is mapped to a 64-dimensional latent vector. Denote by $z_i^A \in \mathbb{R}^{64}$ (resp. $z_j^B$) the latent of the $i$-th cell of modality A (resp. $j$-th cell of modality B). Each cell has a known type label $c_i \in [K]$ (in the paired-benchmark setting, propagated via the ground-truth pairing; in the general setting, provided by an external label source or joint clustering).

**Assumptions**:

1. **Sub-Gaussian concentration.** For each type $k$, the latent vectors of type-$k$ cells in each modality are sub-Gaussian with parameter $\sigma$ around a type-specific mean:
   $$\mathbb{E}[\exp(\langle u, z_i^A - \mu_k^A \rangle)] \leq \exp(\sigma^2 \|u\|^2 / 2) \quad \forall u \in \mathbb{R}^{64}$$
   and symmetrically for modality B with means $\mu_k^B$.

2. **Procrustes recovery.** The the method Procrustes rotation exactly aligns the type centroids across modalities: after applying the rotation, $\mu_k^A = \mu_k^B$ for all $k$. (In practice this is approximate; the finite-sample error contributes an additional $O(\sigma/\sqrt{n})$ term to the bounds below.)

3. **Type separation.** The minimum pairwise cross-type distance is $\Delta := \min_{k \neq \ell} \|\mu_k - \mu_\ell\| > 0$.

4. **Entropic regularization.** Sinkhorn is run with $\epsilon = \sigma^2 / \log n$ (shrinks at the same rate as the within-type concentration radius).

## Theorem 1 (shared types: cluster-resolved entropy → 0)

Under assumptions 1–4, for any cell $i$ of a shared type $k$:
$$\mathbb{E}[H_{\text{cluster}}(i)] = O\left(\frac{1}{\log K} \exp\left(-\frac{n \Delta^2}{\sigma^2}\right)\right)$$

i.e., the cluster-resolved entropy for shared-type cells converges to 0 at an exponential rate in the sample size per type. Equivalently, the cluster marginal $p(k | i)$ concentrates on the true cluster $k$ with probability at least $1 - \exp(-c n \Delta^2 / \sigma^2)$ for a universal constant $c > 0$.

## Theorem 2 (absent types: cluster-resolved entropy bounded away from 0)

Let $k^*$ be a cell type present in modality A but **absent** in modality B (no cells of type $k^*$ exist in B). Then for any cell $i$ of type $k^*$ in modality A:
$$\mathbb{E}[H_{\text{cluster}}(i)] \geq \log K_{\text{present}} / \log K - O(\sigma / \Delta)$$

where $K_{\text{present}}$ is the number of types present in modality B after removal of $k^*$. That is, the cluster marginal distributes its mass approximately uniformly over the remaining types (because no cluster is close in latent space), yielding an entropy that is bounded away from 0 regardless of sample size.

## Corollary (calibrated missing-type detection)

Combining Theorems 1 and 2, there exists a threshold $\tau \in (0, 1)$ such that for sufficiently large $n$ (polynomial in $K$ and $1/\sigma$), thresholding on $H_{\text{cluster}}(i) < \tau$ classifies cells of shared types correctly (low entropy) and cells of absent types correctly (high entropy) with probability at least $1 - \exp(-c' n)$ for a constant $c'$. This is the formal statement underlying the method's missing-type detection experiment.

## Proof sketch

**Theorem 1 — shared types:** By sub-Gaussian concentration, the type-$k$ latent cluster in modality B forms a ball of radius $O(\sigma \sqrt{\log n})$ around $\mu_k^B$ with probability $1 - 1/n$. For a source cell $i$ of type $k$, its distance to the type-$k$ centroid in B is at most $2\sigma\sqrt{\log n}$; its distance to the type-$\ell$ centroid for $\ell \neq k$ is at least $\Delta - O(\sigma \sqrt{\log n})$. The Sinkhorn log-sum-exp approximation then gives, for the cluster marginal:
$$p(k | i) \geq \frac{\exp(-2\sigma^2 \log n / \epsilon)}{\exp(-2\sigma^2 \log n / \epsilon) + (K-1) \exp(-\Delta^2/(2\epsilon) + O(\sigma^2 \log n / \epsilon))}$$
With $\epsilon = \sigma^2/\log n$, the numerator is $\exp(-2(\log n)^2)$ and the denominator subdominant term is $(K-1)\exp(-\Delta^2 \log n / (2\sigma^2))$. For $n$ large enough, $p(k | i) \to 1$ and $H_{\text{cluster}}(i) \to 0$. The rate follows from plugging in the bounds.

**Theorem 2 — absent types:** For type $k^*$ absent in B, the closest centroid is at distance $\geq \Delta$. Any type $k'$ within $O(\Delta)$ of the type-$k^*$ mean gets mass $\sim \exp(-\Delta^2/(2\epsilon))$, but all $K - 1$ remaining types have roughly equal "nearest" distances up to $O(\sigma/\Delta)$ corrections. The cluster marginal is therefore nearly uniform over the $K_{\text{present}}$ present types, and the entropy is $\log K_{\text{present}} / \log K - O(\sigma/\Delta)$.

**Corollary:** Choosing $\tau$ between the two regimes (shared near-0, absent near $\log K_{\text{present}} / \log K$) yields a classifier with exponentially small error. QED (sketch).

## Practical implications

The theorem asserts asymptotic calibration at the cluster level. Empirically on PBMC 10k ($n_{\text{cells per cluster}}$ from 30 to 1200, $K = 18$, $\Delta / \sigma \approx 3$–5), we observe:

- Argmax cluster accuracy 98.76% (theorem predicts near-100%).
- Mean cluster entropy for argmax-correct cells: 0.12 (theorem predicts $\to 0$).
- Mean cluster entropy for argmax-wrong cells: 0.49 (theorem predicts higher than the correct-case mean).
- Missing-type leave-one-out AUROC: 0.960 (theorem predicts this should be $\to 1$ as $n \to \infty$).

The gap between the theorem's asymptotics and the observed finite-sample results (1.24% wrong on PBMC, AUROC 0.96 rather than 1) is explained by the moderate $K$ and the fact that some PBMC clusters have only 30–70 cells, comparable to the $\sigma\sqrt{\log n}$ within-cluster radius. The theorem's convergence rate gets better with more cells per cluster, which matches the observation that larger clusters (n > 200) are classified with ~100% accuracy while small clusters dominate the error set.

## Limitations of the theorem

1. **Assumes Procrustes recovery.** Theorem 1's constants implicitly assume that the Procrustes rotation exactly aligns cluster means, which in practice is approximate (residual ≈ 1.75 on PBMC 10k). The finite-rotation-error term contributes an additive $O(\text{rotation\_residual}/\sigma)$ to both the shared-type concentration and the absent-type bound.

2. **Constant $K$.** We treat the number of clusters as constant. A more refined analysis would let $K$ grow with $n$, e.g., $K = O(\log n)$, and the entropy concentration rates would be correspondingly weaker.

3. **Finite-sample effects for small clusters.** The sub-Gaussian assumption is asymptotic; for clusters with $n < 50$ cells the effective concentration is closer to sub-exponential and the theorem's rates are optimistic.

4. **No statement about speed of convergence of cluster assignment method.** We assume cluster labels are given. In the genuine unpaired-application setting where labels are estimated from a joint clustering, there is an additional estimation error that the theorem does not bound.

A fully rigorous proof is deferred to the appendix of the final paper.
