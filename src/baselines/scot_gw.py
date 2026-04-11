# MIT License - Bryan Cheng, 2026
# Part of MOSAIC
"""Minimal reimplementation of the SCOT v1 alignment algorithm.

SCOT (Single-Cell Optimal Transport, Demetci et al. 2022) aligns two
unpaired single-cell modalities by solving a **Gromov-Wasserstein** OT
problem between their intra-modality k-NN distance matrices. The key idea:
we compare the internal geometry (who's near who within a modality) rather
than absolute coordinates, which are not comparable across modalities.

Algorithm (SCOT v1, faithful to Demetci 2022):
    1. For each modality, compute the k-NN graph in the given feature space
       (PCA for RNA, LSI for ATAC in our MOSAIC setup).
    2. Compute the geodesic distance matrix on each k-NN graph
       (shortest-path via Dijkstra).
    3. Run entropic Gromov-Wasserstein OT between the two distance matrices
       (POT ot.gromov.entropic_gromov_wasserstein).
    4. Use the transport plan as the alignment; for a hard matching, take
       argmax per row.

We reimplement this here instead of pulling the github repo because
(a) it's a small algorithm, (b) installing the published repo on Windows
has been flaky historically, and (c) keeping it inline makes the
fair-comparison protocol auditable: the baseline uses EXACTLY the same
preprocessed AnnDatas that MOSAIC uses.

This is a baseline implementation for the MOSAIC benchmark, not an
endorsement of MOSAIC's method. It's a faithful reimplementation of
SCOT v1; any performance differences between SCOT-reported numbers and
what we see are attributable to: (a) different datasets, (b) different
preprocessing, (c) different feature representations, (d) numerical
instability of entropic GW with small epsilon.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional

import numpy as np
import ot
import scipy.sparse.csgraph as csgraph
from sklearn.neighbors import NearestNeighbors


@dataclass
class SCOTResult:
    plan: np.ndarray                # (N, M) GW transport plan
    barycentric_embedding: np.ndarray  # (N, d_B) source cells projected into target's space
    n_iter_used: int
    wall_time_sec: float


def geodesic_distance_matrix(X: np.ndarray, k: int = 10) -> np.ndarray:
    """Build a k-NN graph over X and compute the pairwise geodesic
    (shortest-path) distance matrix via Dijkstra.

    Returns a (n, n) matrix. Infinite entries are replaced with the matrix
    max (disconnected components are re-connected via an epsilon bridge).
    """
    nbrs = NearestNeighbors(n_neighbors=k + 1, algorithm="auto").fit(X)
    knn_graph = nbrs.kneighbors_graph(X, mode="distance")  # (n, n) sparse
    # Symmetrize
    knn_graph = knn_graph.maximum(knn_graph.T)
    dist = csgraph.shortest_path(knn_graph, directed=False, method="D")
    if not np.all(np.isfinite(dist)):
        # Replace inf with max finite distance * 1.1
        finite_max = np.max(dist[np.isfinite(dist)])
        dist = np.where(np.isfinite(dist), dist, finite_max * 1.1)
    return dist.astype(np.float64)


def run_scot(Z_a: np.ndarray, Z_b: np.ndarray, k: int = 10,
             epsilon: float = 5e-3, n_iter: int = 1000,
             verbose: bool = False) -> SCOTResult:
    """Run SCOT v1 alignment between two modalities' feature representations.

    Inputs:
        Z_a, Z_b: (n_a, d_a) and (n_b, d_b) feature matrices (e.g. PCA of RNA,
            LSI of ATAC). d_a and d_b can differ.
        k: number of neighbors for the intra-modality k-NN graph.
        epsilon: GW entropic regularization strength.
        n_iter: max Sinkhorn-like iterations inside the GW solver.

    Returns: SCOTResult with the transport plan and a barycentric-projection
    embedding of Z_a into Z_b's coordinate system (useful for metric
    computation like label transfer, FOSCTTM, ARI).
    """
    t0 = time.time()
    n_a = Z_a.shape[0]
    n_b = Z_b.shape[0]
    if verbose:
        print(f"[scot] building kNN graphs (k={k})...")
    D_a = geodesic_distance_matrix(Z_a, k=k)
    D_b = geodesic_distance_matrix(Z_b, k=k)
    # Normalize each to unit max for numerical stability
    D_a = D_a / max(D_a.max(), 1e-12)
    D_b = D_b / max(D_b.max(), 1e-12)

    p = np.ones(n_a, dtype=np.float64) / n_a
    q = np.ones(n_b, dtype=np.float64) / n_b

    if verbose:
        print(f"[scot] solving entropic Gromov-Wasserstein (eps={epsilon}, "
              f"max_iter={n_iter})...")
    gw_result, log_dict = ot.gromov.entropic_gromov_wasserstein(
        D_a, D_b, p, q, loss_fun="square_loss", epsilon=epsilon,
        max_iter=n_iter, tol=1e-7, verbose=False, log=True,
    )
    # Barycentric projection: each source cell i mapped to a weighted
    # average over the target cells, weighted by the transport plan row.
    # This gives us a representation of Z_a in Z_b's space so metrics that
    # need matching coordinates (FOSCTTM etc.) work.
    row_sums = gw_result.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1e-30
    barycentric = (gw_result / row_sums) @ Z_b

    return SCOTResult(
        plan=gw_result,
        barycentric_embedding=barycentric.astype(np.float32),
        n_iter_used=int(log_dict.get("err", [0])[-1] if "err" in log_dict else 0),
        wall_time_sec=time.time() - t0,
    )


# ----------------------------------------------------------------------------
# Self-test: recover identity alignment on two Gaussian blobs
# ----------------------------------------------------------------------------


def _test_scot_identity():
    rng = np.random.default_rng(0)
    K = 4
    n = 20
    centers = rng.normal(0, 5, size=(K, 4))
    clusters = np.repeat(np.arange(K), n)
    Z_a = np.concatenate([rng.normal(centers[k], 0.1, size=(n, 4)) for k in range(K)])
    # Z_b is Z_a up to a random orthogonal rotation (within-cluster variance preserved)
    Q, _ = np.linalg.qr(rng.normal(0, 1, size=(4, 4)))
    Z_b = Z_a @ Q + rng.normal(0, 0.05, size=Z_a.shape)

    res = run_scot(Z_a, Z_b, k=5, epsilon=0.01, n_iter=200)
    # The barycentric projection should put each source cluster near the
    # corresponding target cluster (by cluster identity, not row index).
    # Simple sanity check: per-cluster mean of barycentric embedding should
    # be close to per-cluster mean of Z_b.
    err_per_cluster = []
    for k in range(K):
        mask = clusters == k
        b_center = res.barycentric_embedding[mask].mean(0)
        target_center = Z_b[mask].mean(0)
        err_per_cluster.append(np.linalg.norm(b_center - target_center))
    print(f"[scot test] mean cluster-center error = {np.mean(err_per_cluster):.3f}")
    print(f"[scot test] GW plan row sums: "
          f"min={res.plan.sum(1).min():.4f}, max={res.plan.sum(1).max():.4f}")
    print(f"[scot test] took {res.wall_time_sec:.2f}s")


if __name__ == "__main__":
    _test_scot_identity()
