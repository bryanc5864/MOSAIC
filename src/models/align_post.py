# MIT License
# Part of MOSAIC
"""Post-hoc latent alignment between the two independently-trained encoders.

Nothing ties the two IB-VAE latent spaces together during training, so the two
clouds usually sit in different regions of the 64-d space. The OT cost matrix
then comes out near-uniform and Sinkhorn returns a max-entropy plan, which kills
the per-cell uncertainty signal.

Fix is an orthogonal Procrustes rotation fit on cluster centroids rather than
individual cells, so within-cluster geometry survives. On the paired benchmarks
the cluster ids are the leiden labels preprocess.py copies across modalities; a
genuinely unpaired run would need joint clustering or external labels.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class ProcrustesAlignment:
    rotation: np.ndarray         # (d, d) orthogonal rotation matrix
    scale: float                 # isotropic scale factor (1.0 = pure rotation)
    src_mean: np.ndarray         # (d,) source centroid mean (subtracted before rotation)
    tgt_mean: np.ndarray         # (d,) target centroid mean (added after rotation)
    n_clusters: int
    residual: float              # mean ||scale * src_centroid @ R - tgt_centroid||


def fit_orthogonal_procrustes(Z_src: np.ndarray, Z_tgt: np.ndarray,
                               labels_src: np.ndarray, labels_tgt: np.ndarray,
                               include_scale: bool = False,
                               ) -> ProcrustesAlignment:
    """Align src latent centroids onto tgt centroids.

    Both modalities have to share a cluster label set, which holds on the paired
    benchmarks because preprocess.py propagates leiden from RNA to ATAC.

    include_scale=True gives the closed-form similarity Procrustes (Schoenemann
    1966, Sibson 1978). It defaults to False: on PBMC 10k the isotropic scale
    over-compresses the ATAC cloud, which nails the centroid residual (0.49 vs
    1.75) but wrecks within-cluster geometry and label transfer (0.19/0.55 vs
    0.68/0.74, run004 -> run005). Rotation alone is worse on centroids and
    better on everything that matters.
    """
    labels_src = np.asarray(labels_src).astype(str)
    labels_tgt = np.asarray(labels_tgt).astype(str)
    common = sorted(set(labels_src) & set(labels_tgt))
    if len(common) < 2:
        raise ValueError("need at least 2 shared clusters for Procrustes alignment")

    src_centroids = np.stack([Z_src[labels_src == c].mean(axis=0) for c in common])
    tgt_centroids = np.stack([Z_tgt[labels_tgt == c].mean(axis=0) for c in common])

    src_mean = src_centroids.mean(axis=0)
    tgt_mean = tgt_centroids.mean(axis=0)
    A = src_centroids - src_mean       # (K, d)
    B = tgt_centroids - tgt_mean       # (K, d)

    # closed form: svd of A^T B = U S V^T, then R = U V^T
    M = A.T @ B                        # (d, d)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    R = U @ Vt

    if include_scale:
        # minimizes ||c A R - B||_F^2
        denom = float((A ** 2).sum())
        scale = float(S.sum() / denom) if denom > 0 else 1.0
    else:
        scale = 1.0

    residual = float(np.linalg.norm(scale * A @ R - B) / max(len(common), 1))

    return ProcrustesAlignment(
        rotation=R.astype(np.float32),
        scale=scale,
        src_mean=src_mean.astype(np.float32),
        tgt_mean=tgt_mean.astype(np.float32),
        n_clusters=len(common),
        residual=residual,
    )


def apply_alignment(Z: np.ndarray, alignment: ProcrustesAlignment) -> np.ndarray:
    """Apply a fitted Procrustes alignment to a (n, d) embedding."""
    Z_centered = np.asarray(Z, dtype=np.float32) - alignment.src_mean
    return alignment.scale * (Z_centered @ alignment.rotation) + alignment.tgt_mean


def _test_recover_rotation():
    """known rotation in, inverse rotation out"""
    rng = np.random.default_rng(0)
    d = 8
    K = 6
    n_per = 30
    centroids = rng.normal(0, 5, size=(K, d))
    A = rng.normal(0, 1, size=(d, d))
    Q, _ = np.linalg.qr(A)
    tgt_data = []
    src_data = []
    labels = []
    for k in range(K):
        c_tgt = centroids[k]
        c_src = c_tgt @ Q       # rotated
        tgt_data.append(rng.normal(c_tgt, 0.1, size=(n_per, d)))
        src_data.append(rng.normal(c_src, 0.1, size=(n_per, d)))
        labels.extend([f"c{k}"] * n_per)
    Z_tgt = np.concatenate(tgt_data)
    Z_src = np.concatenate(src_data)
    labels = np.array(labels)

    align = fit_orthogonal_procrustes(Z_src, Z_tgt, labels, labels)
    Z_aligned = apply_alignment(Z_src, align)

    err = np.mean([
        np.linalg.norm(Z_aligned[labels == c].mean(0) - Z_tgt[labels == c].mean(0))
        for c in [f"c{k}" for k in range(K)]
    ])
    print(f"[procrustes] residual={align.residual:.4f}, scale={align.scale:.4f}, "
          f"cluster centroid distance after alignment={err:.4f}")
    assert err < 0.1, f"Procrustes failed to align centroids: {err}"
    print("[ok] Procrustes alignment recovers known rotation")

    # now src is also scaled by 0.5
    src_data_scaled = []
    for k in range(K):
        c_src = (centroids[k] @ Q) * 0.5
        src_data_scaled.append(rng.normal(c_src, 0.05, size=(n_per, d)))
    Z_src_scaled = np.concatenate(src_data_scaled)

    align2 = fit_orthogonal_procrustes(Z_src_scaled, Z_tgt, labels, labels, include_scale=True)
    Z_aligned2 = apply_alignment(Z_src_scaled, align2)
    err2 = np.mean([
        np.linalg.norm(Z_aligned2[labels == c].mean(0) - Z_tgt[labels == c].mean(0))
        for c in [f"c{k}" for k in range(K)]
    ])
    print(f"[scaled procrustes] recovered scale={align2.scale:.4f} (true=2.0), centroid err={err2:.4f}")
    assert abs(align2.scale - 2.0) < 0.1, f"scale recovery failed: {align2.scale}"
    assert err2 < 0.1, f"scaled Procrustes centroid error: {err2}"
    print("[ok] similarity Procrustes recovers known scale + rotation")


if __name__ == "__main__":
    _test_recover_rotation()
