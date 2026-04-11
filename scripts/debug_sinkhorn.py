#!/usr/bin/env python3
# MIT License - Bryan Cheng, 2026
"""Debug script to test Sinkhorn behavior on MOSAIC embeddings."""
import sys
import time
sys.path.insert(0, '.')
import numpy as np
import ot

print('loading run004 embeddings...', flush=True)
Zr = np.load('experiments/run004_procrustes/z_rna.npy')
Za = np.load('experiments/run004_procrustes/z_atac_aligned.npy')
print(f'shapes {Zr.shape} {Za.shape}', flush=True)

from src.models.ot_align import pairwise_sqeuclidean, normalize_cost

print('computing cost...', flush=True)
t0 = time.time()
C = pairwise_sqeuclidean(Zr, Za)
print(f'raw cost: max={C.max():.1f} mean={C.mean():.1f} median={np.median(C):.1f} took {time.time()-t0:.1f}s', flush=True)

N = Zr.shape[0]
a = np.ones(N) / N
b = np.ones(N) / N

Cm = normalize_cost(C, method='median')
print(f'median-norm cost: max={Cm.max():.2f} mean={Cm.mean():.2f} median={np.median(Cm):.2f}', flush=True)

# Try at several eps, using regular sinkhorn (not log) with small numItermax
for eps in [0.05, 0.1, 0.5, 1.0, 2.0]:
    print(f'\n=== eps={eps} median-norm, sinkhorn ===', flush=True)
    t0 = time.time()
    try:
        P = ot.sinkhorn(a, b, Cm, reg=eps, numItermax=100, stopThr=1e-6,
                        method='sinkhorn', verbose=False)
        dt = time.time() - t0
        rs = P.sum(axis=1, keepdims=True); rs[rs == 0] = 1e-30
        Pn = P / rs
        H = -(Pn * np.log(Pn + 1e-30)).sum(axis=1) / np.log(N)
        top1 = Pn.max(axis=1)
        print(f'  ok time={dt:.1f}s H_mean={H.mean():.4f} top1_mean={top1.mean():.4f}', flush=True)
        print(f'  H quartiles: {np.percentile(H, [5, 25, 50, 75, 95]).round(4).tolist()}', flush=True)
    except Exception as e:
        print(f'  FAIL after {time.time()-t0:.1f}s: {type(e).__name__}: {e}', flush=True)

# Also try sinkhorn_log at a moderate eps
print('\n=== eps=1.0 median-norm sinkhorn_log (50 iter cap) ===', flush=True)
t0 = time.time()
try:
    P = ot.sinkhorn(a, b, Cm, reg=1.0, numItermax=50, stopThr=1e-6,
                    method='sinkhorn_log', verbose=False)
    dt = time.time() - t0
    rs = P.sum(axis=1, keepdims=True); rs[rs == 0] = 1e-30
    Pn = P / rs
    H = -(Pn * np.log(Pn + 1e-30)).sum(axis=1) / np.log(N)
    print(f'  ok time={dt:.1f}s H_mean={H.mean():.4f}', flush=True)
except Exception as e:
    print(f'  FAIL after {time.time()-t0:.1f}s: {type(e).__name__}: {e}', flush=True)

print('\n[done]', flush=True)
