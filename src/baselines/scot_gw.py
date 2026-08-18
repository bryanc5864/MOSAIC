# MIT License
# Part of MOSAIC
"""Minimal reimplementation of SCOT v1 (Demetci et al. 2022).

Builds a k-NN graph per modality, takes geodesic distances on each, then solves
entropic Gromov-Wasserstein between the two distance matrices. The point of GW
is that it compares internal geometry rather than coordinates, which don't mean
the same thing across modalities.

Reimplemented rather than installed: it's a short algorithm, the published repo
was flaky on Windows, and inlining it means the baseline provably reads the same
preprocessed AnnDatas we do.

Numbers here won't match the SCOT paper's — different datasets, preprocessing,
feature representations, and entropic GW is unstable at small epsilon.
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
    """k-NN graph then Dijkstra shortest paths.

    Disconnected components come back as inf, which we bridge with 1.1x the max
    finite distance.
    """
    nbrs = NearestNeighbors(n_neighbors=k + 1, algorithm="auto").fit(X)
    knn_graph = nbrs.kneighbors_graph(X, mode="distance")  # (n, n) sparse
    knn_graph = knn_graph.maximum(knn_graph.T)
    dist = csgraph.shortest_path(knn_graph, directed=False, method="D")
    if not np.all(np.isfinite(dist)):
        finite_max = np.max(dist[np.isfinite(dist)])
        dist = np.where(np.isfinite(dist), dist, finite_max * 1.1)
    return dist.astype(np.float64)


def run_scot(Z_a: np.ndarray, Z_b: np.ndarray, k: int = 10,
             epsilon: float = 5e-3, n_iter: int = 1000,
             verbose: bool = False) -> SCOTResult:
    """SCOT v1 between two feature matrices; d_a and d_b may differ.

    Returns the transport plan plus a barycentric projection of Z_a into Z_b's
    coordinates, which is what the metrics need.
    """
    t0 = time.time()
    n_a = Z_a.shape[0]
    n_b = Z_b.shape[0]
    if verbose:
        print(f"[scot] building kNN graphs (k={k})...")
    D_a = geodesic_distance_matrix(Z_a, k=k)
    D_b = geodesic_distance_matrix(Z_b, k=k)
    # unit max for numerical stability
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
    # each source cell becomes a plan-weighted average of target cells
    row_sums = gw_result.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1e-30
    barycentric = (gw_result / row_sums) @ Z_b

    return SCOTResult(
        plan=gw_result,
        barycentric_embedding=barycentric.astype(np.float32),
        n_iter_used=int(log_dict.get("err", [0])[-1] if "err" in log_dict else 0),
        wall_time_sec=time.time() - t0,
    )


def _test_scot_identity():
    rng = np.random.default_rng(0)
    K = 4
    n = 20
    centers = rng.normal(0, 5, size=(K, 4))
    clusters = np.repeat(np.arange(K), n)
    Z_a = np.concatenate([rng.normal(centers[k], 0.1, size=(n, 4)) for k in range(K)])
    # Z_b is Z_a under a random rotation
    Q, _ = np.linalg.qr(rng.normal(0, 1, size=(4, 4)))
    Z_b = Z_a @ Q + rng.normal(0, 0.05, size=Z_a.shape)

    res = run_scot(Z_a, Z_b, k=5, epsilon=0.01, n_iter=200)
    # per-cluster means should land near each other
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
