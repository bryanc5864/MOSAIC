# MIT License
# Part of MOSAIC
"""Alignment metrics, shared by the method and every baseline.

Pure numpy/sklearn so the comparison stays apples-to-apples:
  foscttm                 fraction of samples closer than the true match, [0, 0.5]
  label_transfer_accuracy kNN label transfer across the aligned latent, [0, 1]
  joint_clustering_ari    KMeans on the concatenation vs ground truth, [-1, 1]
  entropy_error_corr      spearman of entropy against distance to the true partner;
                          positive means calibrated
  missing_type_auroc      auroc for spotting cells whose type was removed
"""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, roc_auc_score
from sklearn.neighbors import KNeighborsClassifier


def foscttm(Z_a: np.ndarray, Z_b: np.ndarray,
            pair_idx_a: np.ndarray, pair_idx_b: np.ndarray) -> dict:
    """Fraction of samples closer than the true match, both directions.

    pair_idx_a / pair_idx_b are the keys that identify partners: cell i in A
    pairs with the cell in B whose pair_idx_b equals pair_idx_a[i]. Cells with
    no partner (type removed from one side) just drop out of n_paired.
    """
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

    # dense O(N*M); fine at our sizes, would need chunking past ~50K
    Z_a_p = Z_a[paired_a_rows]
    Z_b_p = Z_b
    d = _sqdist_rows_vs_all(Z_a_p, Z_b_p)       # (n, M)
    true_idx_in_b = paired_b_rows
    true_dist = d[np.arange(n), true_idx_in_b]
    # strict < already excludes the true partner
    closer = (d < true_dist[:, None]).sum(axis=1)
    frac = closer / max(Z_b.shape[0] - 1, 1)
    f_ab = float(frac.mean())

    # other direction
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
    """squared euclidean, every row of A against every row of B"""
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    a_sq = (A ** 2).sum(axis=1, keepdims=True)
    b_sq = (B ** 2).sum(axis=1, keepdims=True).T
    D = a_sq + b_sq - 2.0 * A @ B.T
    np.maximum(D, 0.0, out=D)
    return D


def label_transfer_accuracy(Z_src: np.ndarray, labels_src: np.ndarray,
                             Z_tgt: np.ndarray, labels_tgt: np.ndarray,
                             k: int = 15) -> float:
    """kNN fit on the source latent, scored against the target labels."""
    clf = KNeighborsClassifier(n_neighbors=k, metric="euclidean", n_jobs=-1)
    clf.fit(Z_src, labels_src)
    pred = clf.predict(Z_tgt)
    return float((pred == labels_tgt).mean())


def joint_clustering_ari(Z_a: np.ndarray, Z_b: np.ndarray,
                          labels_a: np.ndarray, labels_b: np.ndarray,
                          n_clusters: int | None = None,
                          seed: int = 0) -> float:
    """KMeans on the concatenation, ARI against the concatenated labels.

    n_clusters defaults to the number of distinct labels across both.
    """
    Z = np.concatenate([Z_a, Z_b], axis=0)
    lbl = np.concatenate([np.asarray(labels_a), np.asarray(labels_b)])
    if n_clusters is None:
        n_clusters = len(np.unique(lbl))
    km = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    pred = km.fit_predict(Z)
    return float(adjusted_rand_score(lbl, pred))


def entropy_error_corr(entropy: np.ndarray, Z_a: np.ndarray, Z_b: np.ndarray,
                        pair_idx_a: np.ndarray, pair_idx_b: np.ndarray) -> dict:
    """spearman of per-cell entropy against latent distance to the true partner"""
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


def missing_type_auroc(entropy: np.ndarray, is_missing: np.ndarray) -> float:
    """auroc for entropy detecting cells whose type was dropped from the target"""
    y = np.asarray(is_missing).astype(int)
    if y.min() == y.max():
        return float("nan")
    return float(roc_auc_score(y, entropy))


def _test_foscttm_identity():
    """Z_a == Z_b -> foscttm 0"""
    rng = np.random.default_rng(0)
    Z = rng.normal(0, 1, size=(50, 8))
    pair_idx = np.arange(50)
    out = foscttm(Z, Z, pair_idx, pair_idx)
    print(f"[foscttm identity] mean={out['foscttm_mean']:.4f}")
    assert out["foscttm_mean"] < 1e-6, f"identity alignment should give foscttm 0, got {out['foscttm_mean']}"


def _test_foscttm_random():
    """random alignment sits near 0.5"""
    rng = np.random.default_rng(1)
    Z_a = rng.normal(0, 1, size=(200, 8))
    Z_b = rng.normal(0, 1, size=(200, 8))
    pair_idx = np.arange(200)
    out = foscttm(Z_a, Z_b, pair_idx, pair_idx)
    print(f"[foscttm random] mean={out['foscttm_mean']:.4f}")
    assert 0.4 < out["foscttm_mean"] < 0.6, f"random alignment should be near 0.5, got {out['foscttm_mean']}"


def _test_label_transfer():
    rng = np.random.default_rng(2)
    # 3 clusters, same embedding on both sides
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
    # errors and entropies correlated by construction
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
    print("[ok] metric tests passed")
