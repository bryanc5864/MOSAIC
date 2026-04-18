# MIT License
# Part of MOSAIC
"""Evaluation metrics for single-cell multi-omics alignment.

All metrics here are pure numpy / sklearn functions that take the outputs of
an alignment method and produce a single scalar. The same functions are used
for MOSAIC and for every baseline, so the comparison is apples-to-apples.

Implemented:
  - fosсttm(Z_a, Z_b, pair_idx_a, pair_idx_b)
        Fraction Of Samples Closer Than True Match.
        For each cell in A, compute the fraction of cells in B that are
        closer (in Euclidean distance) to it than its true partner. Average
        over all cells. Range [0, 0.5]; 0.5 is random, 0 is perfect.
  - label_transfer_accuracy(Z_a, Z_b, labels_a, labels_b, k)
        Train a kNN classifier in Z_a's space labeled by labels_a, then
        classify each cell in Z_b by its k nearest neighbors in Z_a and
        compare to labels_b. Range [0, 1].
  - joint_clustering_ari(Z_a, Z_b, labels_a, labels_b, n_clusters)
        Concatenate Z_a and Z_b, cluster jointly (KMeans), compute ARI
        against the union of labels. Range [-1, 1].
  - entropy_error_corr(entropy, Z_a, Z_b, pair_idx_a, pair_idx_b)
        Spearman correlation between per-cell alignment entropy and the
        Euclidean latent distance to the true partner. Positive = calibrated
        (higher entropy <-> higher error).
  - missing_type_auroc(entropy, is_missing)
        Given per-cell entropy and a boolean mask of cells whose true type
        has been removed from the target modality, compute AUROC for
        detecting missing-type cells by thresholding entropy.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, roc_auc_score
from sklearn.neighbors import KNeighborsClassifier


# ----------------------------------------------------------------------------
# FOSCTTM
# ----------------------------------------------------------------------------


def foscttm(Z_a: np.ndarray, Z_b: np.ndarray,
            pair_idx_a: np.ndarray, pair_idx_b: np.ndarray) -> dict:
    """Fraction Of Samples Closer Than True Match.

    Arguments:
        Z_a, Z_b: aligned embeddings for the two modalities, (N, d) and (M, d).
        pair_idx_a, pair_idx_b: integer keys identifying paired cells across
            modalities. For each i in A with pair_idx_a[i] = k, the partner in
            B is the cell whose pair_idx_b equals k.

    Returns dict with:
        - fosсttm_a_to_b: averaged over cells in A
        - fosсttm_b_to_a: averaged over cells in B
        - fosсttm_mean: mean of the two directions
        - n_paired: number of cells that actually have a partner in the
          other modality (a cell type removed from one side contributes 0).
    """
    # Build a lookup: pair_idx_b -> row index in Z_b
    b_lookup = {int(k): i for i, k in enumerate(pair_idx_b)}
    paired_a_rows = []
    paired_b_rows = []
    for i, k in enumerate(pair_idx_a):
        j = b_lookup.get(int(k))
        if j is not None:
            paired_a_rows.append(i)
            paired_b_rows.append(j)
    paired_a_rows = np.asarray(paired_a_rows)
    paired_b_rows = np.asarray(paired_b_rows)
    n = paired_a_rows.size
    if n == 0:
        return {"foscttm_a_to_b": np.nan, "foscttm_b_to_a": np.nan,
                "foscttm_mean": np.nan, "n_paired": 0}

    # Per-cell distances from Z_a[i] to all Z_b (and its true partner at j).
    # Vectorized computation, O(N*M). For M ~ 50K, this is ~2 GB at float32 for N=50K.
    # For our datasets this is fine; we'd chunk for the stretch cross-atlas run.
    Z_a_p = Z_a[paired_a_rows]
    Z_b_p = Z_b  # compare against all Z_b
    d = _sqdist_rows_vs_all(Z_a_p, Z_b_p)       # (n, M)
    true_idx_in_b = paired_b_rows                # index within Z_b of true partner
    true_dist = d[np.arange(n), true_idx_in_b]   # (n,)
    # Count how many cells in Z_b are strictly closer than the true partner,
    # excluding the true partner itself.
    closer = (d < true_dist[:, None]).sum(axis=1)
    # Exclude the true partner (distance is 0 or positive; strict < handles it)
    frac = closer / max(Z_b.shape[0] - 1, 1)
    f_ab = float(frac.mean())

    # Symmetric direction: B -> A
    Z_b_q = Z_b[paired_b_rows]
    d2 = _sqdist_rows_vs_all(Z_b_q, Z_a)
    true_dist_ba = d2[np.arange(n), paired_a_rows]
    closer_ba = (d2 < true_dist_ba[:, None]).sum(axis=1)
    frac_ba = closer_ba / max(Z_a.shape[0] - 1, 1)
    f_ba = float(frac_ba.mean())

    return {
        "foscttm_a_to_b": f_ab,
        "foscttm_b_to_a": f_ba,
        "foscttm_mean": 0.5 * (f_ab + f_ba),
        "n_paired": int(n),
    }


def _sqdist_rows_vs_all(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Squared-Euclidean distance between each row of A and each row of B."""
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    a_sq = (A ** 2).sum(axis=1, keepdims=True)
    b_sq = (B ** 2).sum(axis=1, keepdims=True).T
    D = a_sq + b_sq - 2.0 * A @ B.T
    np.maximum(D, 0.0, out=D)
    return D


# ----------------------------------------------------------------------------
# Label transfer accuracy
# ----------------------------------------------------------------------------


def label_transfer_accuracy(Z_src: np.ndarray, labels_src: np.ndarray,
                             Z_tgt: np.ndarray, labels_tgt: np.ndarray,
                             k: int = 15) -> float:
    """Train kNN on (Z_src, labels_src), predict for Z_tgt, compare to labels_tgt.

    This emulates "transfer cell type labels from RNA to ATAC via the shared
    aligned latent space". Higher is better.
    """
    clf = KNeighborsClassifier(n_neighbors=k, metric="euclidean", n_jobs=-1)
    clf.fit(Z_src, labels_src)
    pred = clf.predict(Z_tgt)
    return float((pred == labels_tgt).mean())


# ----------------------------------------------------------------------------
# Joint-clustering ARI
# ----------------------------------------------------------------------------


def joint_clustering_ari(Z_a: np.ndarray, Z_b: np.ndarray,
                          labels_a: np.ndarray, labels_b: np.ndarray,
                          n_clusters: int | None = None,
                          seed: int = 0) -> float:
    """Cluster the concatenation of Z_a and Z_b with KMeans, compute ARI
    against the concatenated ground-truth labels.

    If n_clusters is None, use the number of distinct labels across both.
    """
    Z = np.concatenate([Z_a, Z_b], axis=0)
    lbl = np.concatenate([np.asarray(labels_a), np.asarray(labels_b)])
    if n_clusters is None:
        n_clusters = len(np.unique(lbl))
    km = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    pred = km.fit_predict(Z)
    return float(adjusted_rand_score(lbl, pred))


# ----------------------------------------------------------------------------
# Entropy calibration
# ----------------------------------------------------------------------------


def entropy_error_corr(entropy: np.ndarray, Z_a: np.ndarray, Z_b: np.ndarray,
                        pair_idx_a: np.ndarray, pair_idx_b: np.ndarray) -> dict:
    """Spearman correlation between per-cell MOSAIC entropy and the
    Euclidean latent distance to the true partner in the other modality.

    Returns dict with spearman rho and p-value.
    """
    b_lookup = {int(k): i for i, k in enumerate(pair_idx_b)}
    errs = []
    ents = []
    for i, k in enumerate(pair_idx_a):
        j = b_lookup.get(int(k))
        if j is None:
            continue
        err = float(np.sqrt(((Z_a[i] - Z_b[j]) ** 2).sum()))
        errs.append(err)
        ents.append(float(entropy[i]))
    if len(errs) < 3:
        return {"spearman_rho": np.nan, "spearman_p": np.nan, "n": len(errs)}
    rho, p = spearmanr(ents, errs)
    return {"spearman_rho": float(rho), "spearman_p": float(p), "n": len(errs)}


# ----------------------------------------------------------------------------
# Missing-type detection
# ----------------------------------------------------------------------------


def missing_type_auroc(entropy: np.ndarray, is_missing: np.ndarray) -> float:
    """AUROC for using per-cell entropy to detect cells whose true type
    has been removed from the other modality.
    """
    y = np.asarray(is_missing).astype(int)
    if y.min() == y.max():
        return float("nan")
    return float(roc_auc_score(y, entropy))


# ----------------------------------------------------------------------------
# Tests on toy data
# ----------------------------------------------------------------------------


def _test_foscttm_identity():
    """Identity alignment: Z_a == Z_b -> foscttm == 0."""
    rng = np.random.default_rng(0)
    Z = rng.normal(0, 1, size=(50, 8))
    pair_idx = np.arange(50)
    out = foscttm(Z, Z, pair_idx, pair_idx)
    print(f"[foscttm identity] mean={out['foscttm_mean']:.4f}")
    assert out["foscttm_mean"] < 1e-6, f"identity alignment should give foscttm 0, got {out['foscttm_mean']}"


def _test_foscttm_random():
    """Random alignment: expected foscttm near 0.5."""
    rng = np.random.default_rng(1)
    Z_a = rng.normal(0, 1, size=(200, 8))
    Z_b = rng.normal(0, 1, size=(200, 8))
    pair_idx = np.arange(200)
    out = foscttm(Z_a, Z_b, pair_idx, pair_idx)
    print(f"[foscttm random] mean={out['foscttm_mean']:.4f}")
    assert 0.4 < out["foscttm_mean"] < 0.6, f"random alignment should be near 0.5, got {out['foscttm_mean']}"


def _test_label_transfer():
    rng = np.random.default_rng(2)
    # 3 clusters, identical embeddings across "modalities"
    centers = rng.normal(0, 5, size=(3, 4))
    labels = np.repeat([0, 1, 2], 20)
    Z = np.concatenate([
        centers[l] + rng.normal(0, 0.1, size=(1, 4)) for l in labels
    ])
    acc = label_transfer_accuracy(Z, labels, Z, labels, k=5)
    print(f"[label transfer identical] acc={acc:.3f}")
    assert acc > 0.95


def _test_ari_identity():
    rng = np.random.default_rng(3)
    centers = rng.normal(0, 10, size=(4, 6))
    labels = np.repeat(np.arange(4), 15)
    Z = np.concatenate([
        centers[l] + rng.normal(0, 0.1, size=(1, 6)) for l in labels
    ])
    ari = joint_clustering_ari(Z, Z, labels, labels, n_clusters=4, seed=0)
    print(f"[ari identity] {ari:.3f}")
    assert ari > 0.95


def _test_entropy_corr_toy():
    rng = np.random.default_rng(4)
    # Construct errors and entropies with known positive correlation
    n = 50
    err = rng.uniform(0, 1, size=n)
    noise = rng.normal(0, 0.1, size=n)
    ent = err + noise
    Z_a = rng.normal(0, 1, size=(n, 4))
    Z_b = Z_a.copy() + err[:, None] * rng.normal(0, 1, size=(n, 4))  # noise scales with err
    pair_idx = np.arange(n)
    res = entropy_error_corr(ent, Z_a, Z_b, pair_idx, pair_idx)
    print(f"[entropy-error] rho={res['spearman_rho']:.3f}")
    assert res["spearman_rho"] > 0.3


def _test_missing_type_auroc():
    rng = np.random.default_rng(5)
    ent = np.concatenate([rng.normal(0.2, 0.05, 30), rng.normal(0.8, 0.05, 30)])
    missing = np.concatenate([np.zeros(30), np.ones(30)]).astype(bool)
    auroc = missing_type_auroc(ent, missing)
    print(f"[missing-type auroc] {auroc:.3f}")
    assert auroc > 0.95


if __name__ == "__main__":
    _test_foscttm_identity()
    _test_foscttm_random()
    _test_label_transfer()
    _test_ari_identity()
    _test_entropy_corr_toy()
    _test_missing_type_auroc()
    print("[ok] evaluation metrics unit tests passed")
