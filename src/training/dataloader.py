# MIT License
# Part of MOSAIC
"""Per-modality AnnData dataloader.

Design invariant (enforced by the pair-leakage test in src/data/validate.py):
  The training-time iterator for modality A has NO access to modality B.
  Each modality is wrapped in its own Dataset; batches are sampled
  independently with independent RNG seeds; and nothing in this module reads
  the `pair_idx` field (only evaluation code does).

The `ModalityDataset` yields per-cell:
  - x         : the modality's MLP-input feature vector (scaled log-norm for
                RNA; log1p TF-IDF for ATAC).
  - recon_tgt : the reconstruction target for the decoder (raw counts for
                RNA, binarized for ATAC).
  - y_cross   : the pre-computed low-dim summary of the OTHER modality, used
                as the cross-modal prediction target.

The three are pre-computed and stored on the AnnData object, so the dataset
does not need to know which modality is which — it just reads the fields.
"""

from __future__ import annotations

from typing import Literal

import anndata as ad
import numpy as np
import scipy.sparse as sp
import torch
from torch.utils.data import DataLoader, Dataset


class ModalityDataset(Dataset):
    """Single-modality Dataset over an AnnData.

    Arguments:
        adata:      a single modality AnnData. Must contain .X (features),
                    a recon target layer (specified by `recon_layer`), and
                    .obsm['y_cross'] (the cross-modal prediction target).
        recon_layer: which layer to use as the reconstruction target. For
                    RNA this is "counts"; for ATAC this is "binary".
        indices:    optional int array selecting a subset of cells
                    (e.g. for train/val/test splits).
    """

    def __init__(self, adata: ad.AnnData, recon_layer: str,
                 indices: np.ndarray | None = None):
        self.adata = adata
        self.recon_layer = recon_layer
        if indices is None:
            self.indices = np.arange(adata.n_obs)
        else:
            self.indices = np.asarray(indices, dtype=np.int64)
        # Preload everything as dense float32 torch tensors at init time.
        # For 11K cells x 10K peaks this is ~450 MB — large but manageable
        # and eliminates per-row sparse-to-dense conversion overhead.
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
        """Vectorized batch access — pull many rows in one shot.

        Preferred path used by the training loop (see train_ibvae.py).
        """
        rows = torch.from_numpy(np.asarray(batch_indices, dtype=np.int64))
        abs_rows = torch.from_numpy(self.indices[batch_indices].astype(np.int64))
        return {
            "x": self._X.index_select(0, abs_rows),
            "recon": self._recon.index_select(0, abs_rows),
            "y_cross": self._y_cross.index_select(0, abs_rows),
        }

    def __getitem__(self, idx: int) -> dict[str, torch.Tensor]:
        # Kept for DataLoader compatibility; the training loop uses get_batch
        # which is much faster.
        row = int(self.indices[idx])
        return {
            "x": self._X[row],
            "recon": self._recon[row],
            "y_cross": self._y_cross[row],
        }


def make_split_indices(n: int, val_frac: float = 0.1, seed: int = 0
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Return (train_idx, val_idx) from a random split of n cells.

    Note: this does NOT produce a test split. The full dataset is used for
    OT-time matching during final evaluation; val exists only for early
    stopping of the IB-VAE training loop.
    """
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    n_val = max(int(n * val_frac), 1)
    return perm[n_val:], perm[:n_val]


def make_loader(ds: ModalityDataset, batch_size: int, shuffle: bool,
                seed: int = 0, num_workers: int = 0, pin_memory: bool = True
                ) -> DataLoader:
    """Wrap a ModalityDataset in a DataLoader with a deterministic generator
    when shuffling is enabled.
    """
    if shuffle:
        gen = torch.Generator()
        gen.manual_seed(seed)
    else:
        gen = None
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle,
                      num_workers=num_workers, pin_memory=pin_memory,
                      generator=gen, drop_last=False)
