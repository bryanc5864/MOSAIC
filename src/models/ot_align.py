# MIT License
# Part of MOSAIC
"""Entropic OT alignment and the per-cell alignment entropy.

Given latents Z_A (N x d) and Z_B (M x d): squared-euclidean cost, rescale it,
solve P* = argmin <C, P> + epsilon * KL(P || a b^T) with uniform marginals, then
per source cell take H_i = -sum_j p_ij log p_ij normalized by log M so it sits
in [0, 1].

Dense POT sinkhorn only. Everything here fits in memory at the sizes we run.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import ot                     # POT
import torch


def pairwise_sqeuclidean(Z_a: np.ndarray, Z_b: np.ndarray) -> np.ndarray:
    """Pairwise squared-Euclidean distance, N x M."""
    Z_a = np.asarray(Z_a, dtype=np.float64)
    Z_b = np.asarray(Z_b, dtype=np.float64)
    a_sq = (Z_a ** 2).sum(axis=1, keepdims=True)         # (N, 1)
    b_sq = (Z_b ** 2).sum(axis=1, keepdims=True).T       # (1, M)
    C = a_sq + b_sq - 2.0 * Z_a @ Z_b.T
    np.maximum(C, 0.0, out=C)                             # kill negative roundoff
    return C


def normalize_cost(C: np.ndarray, method: str = "median") -> np.ndarray:
    """Rescale the cost matrix so epsilon means the same thing across datasets.

    'max' is the original behaviour and is outlier-sensitive; run-004 showed it
    inflating the effective epsilon ~30x off a handful of extreme cell pairs.
    'median' is the default now, 'mean' sits in between.
    """
    if method == "max":
        scale = float(C.max())
    elif method == "mean":
        scale = float(C.mean())
    elif method == "median":
        # np.median on a flattened 11k^2 array is too slow, so sample
        flat = C.ravel()
        if flat.size > 1_000_000:
            rng = np.random.default_rng(0)
            sample = rng.choice(flat, size=1_000_000, replace=False)
            scale = float(np.median(sample))
        else:
            scale = float(np.median(flat))
    else:
        raise ValueError(f"unknown normalize_cost method: {method!r}")
    if scale <= 0:
        return C
    return C / scale


@dataclass
class AlignmentResult:
    plan: np.ndarray                  # (N, M) transport plan (marginals a, b)
    entropy: np.ndarray               # (N,) normalized per-row entropy
    top_match: np.ndarray             # (N,) argmax of each row -> index in modality B
    top_match_prob: np.ndarray        # (N,) prob assigned to the argmax
    cost_scale: float                 # scaling applied to C (for reproducibility)
    epsilon: float
    n_iters: int
    converged: bool
    info: dict = field(default_factory=dict)


def sinkhorn_align(Z_a: np.ndarray, Z_b: np.ndarray,
                   epsilon: float = 0.5,
                   n_iter: int = 200,
                   stop_threshold: float = 1e-6,
                   normalize: str = "median",
                   sinkhorn_method: str = "sinkhorn",
                   a: np.ndarray | None = None,
                   b: np.ndarray | None = None) -> AlignmentResult:
    """Entropic OT between two embeddings; returns the plan and per-row entropy.

    normalize: how to rescale the cost before Sinkhorn — 'median' (default),
    'mean', 'max' (legacy, outlier-sensitive), or 'none' (epsilon in raw units).

    sinkhorn_method: plain 'sinkhorn' is the default and the only one that
    survives 11K x 11K on a 16 GB box — sinkhorn_log's scipy.logsumexp blows up
    the RAM and gets OOM-killed with no traceback. 'sinkhorn_log' is stable at
    tiny epsilon but small-problems-only; 'sinkhorn_stabilized' is the middle.
    """
    N, d = Z_a.shape
    M, d2 = Z_b.shape
    assert d == d2, f"embedding dims differ: {d} vs {d2}"

    if a is None:
        a = np.ones(N, dtype=np.float64) / N
    if b is None:
        b = np.ones(M, dtype=np.float64) / M

    C = pairwise_sqeuclidean(Z_a, Z_b)
    if normalize == "none":
        C_norm = C
        scale = 1.0
    else:
        C_norm = normalize_cost(C, method=normalize)
        # only used for logging
        raw_median = float(np.median(C))
        norm_median = float(np.median(C_norm))
        scale = raw_median / max(norm_median, 1e-30)

    P = ot.sinkhorn(a, b, C_norm, reg=epsilon, numItermax=n_iter,
                    stopThr=stop_threshold, method=sinkhorn_method,
                    verbose=False, log=False)
    P = np.asarray(P, dtype=np.float64)

    # rows to probabilities
    row_sums = P.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1e-30
    P_row = P / row_sums
    # log M normalization puts H in [0, 1]
    log_M = float(np.log(max(M, 2)))
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(P_row > 0, np.log(P_row), 0.0)
    H = -(P_row * logp).sum(axis=1) / log_M

    top_match = P_row.argmax(axis=1)
    top_match_prob = P_row[np.arange(N), top_match]

    return AlignmentResult(
        plan=P,
        entropy=H,
        top_match=top_match,
        top_match_prob=top_match_prob,
        cost_scale=scale,
        epsilon=epsilon,
        n_iters=n_iter,
        converged=True,  # POT raises on non-convergence in log mode
        info={"cost_max_raw": float(C.max()), "cost_mean_raw": float(C.mean())},
    )


def _test_identity_alignment():
    """same blobs in the same places -> near-identity plan"""
    rng = np.random.default_rng(0)
    K = 4
    n_per = 10
    centers = np.array([[0, 0], [5, 0], [0, 5], [5, 5]], dtype=np.float64)
    Za = np.concatenate([rng.normal(c, 0.1, size=(n_per, 2)) for c in centers])
    Zb = np.concatenate([rng.normal(c, 0.1, size=(n_per, 2)) for c in centers])
    # shuffle so the answer isn't the identity permutation
    perm = rng.permutation(Za.shape[0])
    Zb = Zb[perm]
    res = sinkhorn_align(Za, Zb, epsilon=0.01)
    cluster_ids_a = np.repeat(np.arange(K), n_per)
    cluster_ids_b = cluster_ids_a[perm]
    predicted_cluster = cluster_ids_b[res.top_match]
    acc = (predicted_cluster == cluster_ids_a).mean()
    print(f"[cluster-match] accuracy = {acc:.3f}")
    assert acc > 0.9, f"argmax-matching accuracy too low: {acc:.3f}"
    print(f"[entropy stats] mean = {res.entropy.mean():.3f}, "
          f"min = {res.entropy.min():.3f}, max = {res.entropy.max():.3f}")


def _test_high_entropy_on_collapsed_blob():
    """identical targets -> uniform rows -> entropy near 1"""
    rng = np.random.default_rng(1)
    Za = rng.normal(0, 1, size=(40, 2))
    Zb = np.zeros((40, 2))
    res = sinkhorn_align(Za, Zb, epsilon=0.1)
    print(f"[collapsed] mean entropy = {res.entropy.mean():.3f} (expect near 1)")
    assert res.entropy.mean() > 0.9, f"collapsed target should yield high entropy, got {res.entropy.mean():.3f}"


def _test_low_entropy_on_well_separated():
    """well-separated blobs with unique partners -> low entropy"""
    rng = np.random.default_rng(2)
    K = 20
    centers = rng.normal(0, 10, size=(K, 8))
    Za = centers + rng.normal(0, 0.01, size=(K, 8))
    Zb = centers + rng.normal(0, 0.01, size=(K, 8))
    Zb = Zb[rng.permutation(K)]
    res = sinkhorn_align(Za, Zb, epsilon=0.005)
    print(f"[well-sep] mean entropy = {res.entropy.mean():.3f} (expect near 0)")
    assert res.entropy.mean() < 0.3, f"well-separated should yield low entropy, got {res.entropy.mean():.3f}"


if __name__ == "__main__":
    _test_identity_alignment()
    _test_high_entropy_on_collapsed_blob()
    _test_low_entropy_on_well_separated()
    print("[ok] ot alignment smoke tests passed")
