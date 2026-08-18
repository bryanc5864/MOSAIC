# MIT License
# Part of MOSAIC
"""Per-modality AnnData dataloader.

The invariant that matters: modality A's iterator never touches modality B.
Each modality gets its own Dataset, batches are drawn with independent RNGs,
and nothing here reads pair_idx — only the evaluation code does.

Each item is x (encoder input), recon (decoder target: raw counts for RNA,
binary for ATAC), and y_cross (the partner modality's cluster centroid). All
three are precomputed on the AnnData, so this class doesn't care which
modality it's holding.
"""

from __future__ import annotations

from typing import Literal

import anndata as ad
import numpy as np
import scipy.sparse as sp
import torch
from torch.utils.data import DataLoader, Dataset


class ModalityDataset(Dataset):
    """One modality's cells. recon_layer is "counts" for RNA, "binary" for ATAC."""

    def __init__(self, adata: ad.AnnData, recon_layer: str,
                 indices: np.ndarray | None = None):
        self.adata = adata
        self.recon_layer = recon_layer
        if indices is None:
            self.indices = np.arange(adata.n_obs)
        else:
            self.indices = np.asarray(indices, dtype=np.int64)
        # dense up front: ~450 MB at 11K x 10K, but no per-row densify later
        self._X = self._dense_tensor(adata.X)
        self._recon = self._dense_tensor(adata.layers[recon_layer])
        self._y_cross = torch.from_numpy(np.asarray(adata.obsm["y_cross"], dtype=np.float32))
        self._n_obs = adata.n_obs

    @staticmethod
    def _dense_tensor(X) -> torch.Tensor:
        if sp.issparse(X):
            X = X.toarray()
        return torch.from_numpy(np.asarray(X, dtype=np.float32))

    def __len__(self) -> int:
        return len(self.indices)

    def get_batch(self, batch_indices: np.ndarray) -> dict[str, torch.Tensor]:
        """pull a whole batch at once; this is what the training loop uses"""
        rows = torch.from_numpy(np.asarray(batch_indices, dtype=np.int64))
        abs_rows = torch.from_numpy(self.indices[batch_indices].astype(np.int64))
        return {
            "x": self._X.index_select(0, abs_rows),
            "recon": self._recon.index_select(0, abs_rows),
            "y_cross": self._y_cross.index_select(0, abs_rows),
        }

    def __getitem__(self, idx: int) -> dict[str, torch.Tensor]:
        # only here for DataLoader compat, get_batch is much faster
        row = int(self.indices[idx])
        return {
            "x": self._X[row],
            "recon": self._recon[row],
            "y_cross": self._y_cross[row],
        }


def make_split_indices(n: int, val_frac: float = 0.1, seed: int = 0
                       ) -> tuple[np.ndarray, np.ndarray]:
    """random (train_idx, val_idx) over n cells.

    No test split on purpose: final evaluation matches over all cells at OT
    time, so val is only here to drive early stopping.
    """
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    n_val = max(int(n * val_frac), 1)
    return perm[n_val:], perm[:n_val]


def make_loader(ds: ModalityDataset, batch_size: int, shuffle: bool,
                seed: int = 0, num_workers: int = 0, pin_memory: bool = True
                ) -> DataLoader:
    """DataLoader with a seeded generator when shuffling"""
    if shuffle:
        gen = torch.Generator()
        gen.manual_seed(seed)
    else:
        gen = None
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle,
                      num_workers=num_workers, pin_memory=pin_memory,
                      generator=gen, drop_last=False)
