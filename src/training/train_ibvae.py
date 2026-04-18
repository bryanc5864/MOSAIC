# MIT License
# Part of MOSAIC
"""Train a single-modality IB-VAE end to end.

Inputs:
  - processed AnnData for the modality
  - hyperparameters (from config yaml or programmatic kwargs)

Produces:
  - trained model state_dict (saved under experiments/<exp>/ckpt_<modality>.pt)
  - per-epoch loss log (saved to experiments/<exp>/train_log_<modality>.json)
  - final latent embedding for ALL cells (saved as experiments/<exp>/z_<modality>.npy)

Training protocol (RESEARCH_PLAN.md section 3.4):
  - AdamW optimizer, lr=1e-3, weight_decay=1e-4
  - Cosine schedule with 10-epoch linear warmup
  - Batch size 512
  - KL beta ramps linearly from 0 to beta over first 30 epochs
  - Early stopping on validation cross-modal prediction MSE, patience=15
  - All seeds set

Deterministic except for cuDNN kernel selection, which is disabled via
torch.use_deterministic_algorithms where possible.
"""

from __future__ import annotations

import json
import math
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import torch
from torch.optim import AdamW

from src.models.ib_vae import IBVAE_ATAC, IBVAE_RNA
from src.training.dataloader import ModalityDataset, make_split_indices


# ----------------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------------


@dataclass
class TrainConfig:
    # Required
    modality: str                    # "rna" or "atac"
    processed_path: str              # path to the h5ad file
    out_dir: str                     # where to save outputs

    # Optimization
    batch_size: int = 512
    lr: float = 1e-3
    weight_decay: float = 1e-4
    n_epochs: int = 200
    warmup_epochs: int = 10
    beta_warmup_epochs: int = 30
    beta: float = 0.01
    lambda_pred: float = 1.0
    patience: int = 15
    val_frac: float = 0.1

    # Model
    latent_dim: int = 64
    hidden_enc: tuple[int, ...] = (512, 256)
    hidden_dec: tuple[int, ...] = (256, 512)
    dropout: float = 0.1

    # Reproducibility
    seed: int = 0
    device: str = "cuda"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _cosine_lr(step: int, total_steps: int, warmup_steps: int, base_lr: float) -> float:
    """Linear warmup then cosine decay to 0."""
    if step < warmup_steps:
        return base_lr * (step + 1) / max(warmup_steps, 1)
    progress = (step - warmup_steps) / max(total_steps - warmup_steps, 1)
    return 0.5 * base_lr * (1 + math.cos(math.pi * progress))


def _beta_schedule(epoch: int, warmup: int, target_beta: float) -> float:
    """Linear warmup of KL weight from 0 to target over `warmup` epochs."""
    if epoch >= warmup:
        return target_beta
    return target_beta * (epoch + 1) / max(warmup, 1)


def _set_seed(seed: int) -> None:
    """Set all RNG seeds and request deterministic cuDNN kernels.

    Note: cuDNN deterministic mode disables algorithm autotuning, which
    slightly slows convolution-heavy models. MOSAIC uses MLPs only, so the
    cost is negligible. We request `warn_only=True` from
    use_deterministic_algorithms so torch will warn (not error) on any op
    without a deterministic kernel.
    """
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    try:
        torch.use_deterministic_algorithms(True, warn_only=True)
    except Exception:
        pass  # torch versions before 1.8 don't support this


# ----------------------------------------------------------------------------
# Main training loop
# ----------------------------------------------------------------------------


def train(cfg: TrainConfig) -> dict:
    _set_seed(cfg.seed)

    device = torch.device(cfg.device if torch.cuda.is_available() else "cpu")
    if cfg.device == "cuda" and not torch.cuda.is_available():
        print("[warn] CUDA requested but not available; falling back to CPU")

    adata = ad.read_h5ad(cfg.processed_path)
    n_vars = adata.n_vars
    n_cells = adata.n_obs
    cross_dim = int(adata.obsm["y_cross"].shape[1])

    train_idx, val_idx = make_split_indices(n_cells, val_frac=cfg.val_frac, seed=cfg.seed)
    recon_layer = "counts" if cfg.modality == "rna" else "binary"
    # One shared dataset (all cells) — the train/val/full "datasets" only
    # differ in which indices the training loop iterates over. This avoids
    # re-allocating the dense tensors three times.
    full_ds = ModalityDataset(adata, recon_layer=recon_layer, indices=None)
    # Move the full feature and target tensors onto the device once. For
    # PBMC 10k this is ~50 MB (RNA) or ~450 MB (ATAC); well within the
    # 10 GB budget on the RTX 3080. Skips CPU->GPU transfer per batch.
    X_dev = full_ds._X.to(device, non_blocking=True)
    recon_dev = full_ds._recon.to(device, non_blocking=True)
    yc_dev = full_ds._y_cross.to(device, non_blocking=True)
    n_train = len(train_idx)
    n_val = len(val_idx)

    if cfg.modality == "rna":
        model = IBVAE_RNA(n_vars=n_vars, cross_dim=cross_dim, latent_dim=cfg.latent_dim,
                          hidden_enc=cfg.hidden_enc, hidden_dec=cfg.hidden_dec,
                          dropout=cfg.dropout)
    elif cfg.modality == "atac":
        model = IBVAE_ATAC(n_vars=n_vars, cross_dim=cross_dim, latent_dim=cfg.latent_dim,
                           hidden_enc=cfg.hidden_enc, hidden_dec=cfg.hidden_dec,
                           dropout=cfg.dropout)
    else:
        raise ValueError(f"modality must be 'rna' or 'atac', got {cfg.modality!r}")
    model = model.to(device)

    opt = AdamW(model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)

    steps_per_epoch = max(1, (n_train + cfg.batch_size - 1) // cfg.batch_size)
    total_steps = steps_per_epoch * cfg.n_epochs
    warmup_steps = steps_per_epoch * cfg.warmup_epochs

    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    epoch_log: list[dict] = []
    best_val = float("inf")
    best_epoch = -1
    patience_counter = 0
    best_state_dict = None

    t0 = time.time()
    step = 0
    gen = torch.Generator()
    gen.manual_seed(cfg.seed)

    def _gather_on_device(idx_in_subset: np.ndarray, subset_indices: np.ndarray):
        abs_rows = torch.from_numpy(subset_indices[idx_in_subset].astype(np.int64)).to(device)
        x = X_dev.index_select(0, abs_rows)
        recon = recon_dev.index_select(0, abs_rows)
        yc = yc_dev.index_select(0, abs_rows)
        return x, recon, yc

    for epoch in range(cfg.n_epochs):
        beta = _beta_schedule(epoch, cfg.beta_warmup_epochs, cfg.beta)

        # ----- train -----
        model.train()
        sum_tot = sum_rec = sum_kl = sum_pred = 0.0
        n_batches = 0
        perm = torch.randperm(n_train, generator=gen).numpy()
        for start in range(0, n_train, cfg.batch_size):
            end = min(start + cfg.batch_size, n_train)
            idx_in_subset = perm[start:end]
            x, recon, yc = _gather_on_device(idx_in_subset, train_idx)
            lr_t = _cosine_lr(step, total_steps, warmup_steps, cfg.lr)
            for g in opt.param_groups:
                g["lr"] = lr_t
            opt.zero_grad()
            out = model(x, recon, yc, beta=beta, lambda_pred=cfg.lambda_pred)
            out.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            opt.step()
            sum_tot += out.total.item()
            sum_rec += out.recon.item()
            sum_kl += out.kl.item()
            sum_pred += out.pred.item()
            n_batches += 1
            step += 1
        train_metrics = {
            "train_total": sum_tot / max(n_batches, 1),
            "train_recon": sum_rec / max(n_batches, 1),
            "train_kl": sum_kl / max(n_batches, 1),
            "train_pred": sum_pred / max(n_batches, 1),
        }

        # ----- val -----
        model.eval()
        val_pred = 0.0
        val_tot = 0.0
        n_val_batches = 0
        with torch.no_grad():
            for start in range(0, n_val, cfg.batch_size):
                end = min(start + cfg.batch_size, n_val)
                idx_in_subset = np.arange(start, end)
                x, recon, yc = _gather_on_device(idx_in_subset, val_idx)
                out = model(x, recon, yc, beta=beta, lambda_pred=cfg.lambda_pred)
                val_pred += out.pred.item()
                val_tot += out.total.item()
                n_val_batches += 1
        val_pred /= max(n_val_batches, 1)
        val_tot /= max(n_val_batches, 1)

        epoch_log.append({
            "epoch": epoch,
            "beta": beta,
            "lr": lr_t,
            **train_metrics,
            "val_pred": val_pred,
            "val_total": val_tot,
        })

        if val_pred < best_val - 1e-5:
            best_val = val_pred
            best_epoch = epoch
            patience_counter = 0
            best_state_dict = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        else:
            patience_counter += 1

        if (epoch + 1) % 10 == 0 or epoch == 0 or epoch == cfg.n_epochs - 1:
            elapsed = time.time() - t0
            print(f"  epoch {epoch+1:3d}/{cfg.n_epochs}  "
                  f"train_total={train_metrics['train_total']:.4f}  "
                  f"train_recon={train_metrics['train_recon']:.4f}  "
                  f"kl={train_metrics['train_kl']:.4f}  "
                  f"pred={train_metrics['train_pred']:.4f}  "
                  f"val_pred={val_pred:.4f}  "
                  f"beta={beta:.4f}  "
                  f"[{elapsed:.1f}s]")

        if patience_counter >= cfg.patience:
            print(f"  [early-stop] patience exhausted at epoch {epoch+1}, "
                  f"best val_pred={best_val:.4f} @ epoch {best_epoch+1}")
            break

    # Restore best checkpoint
    if best_state_dict is not None:
        model.load_state_dict(best_state_dict)
    else:
        best_epoch = epoch  # no early stop, use final

    # ----- full-data embedding -----
    model.eval()
    embeddings = np.zeros((n_cells, cfg.latent_dim), dtype=np.float32)
    with torch.no_grad():
        for start in range(0, n_cells, cfg.batch_size):
            end = min(start + cfg.batch_size, n_cells)
            x = X_dev[start:end]
            mu, _ = model.encoder(x)
            embeddings[start:end] = mu.detach().cpu().numpy()

    # ----- save -----
    ckpt_path = out_dir / f"ckpt_{cfg.modality}.pt"
    emb_path = out_dir / f"z_{cfg.modality}.npy"
    log_path = out_dir / f"train_log_{cfg.modality}.json"

    torch.save({"state_dict": model.state_dict(), "config": asdict(cfg),
                "best_epoch": best_epoch, "best_val": best_val}, ckpt_path)
    np.save(emb_path, embeddings)
    with log_path.open("w") as f:
        json.dump({"config": asdict(cfg),
                   "epochs": epoch_log,
                   "best_val": best_val,
                   "best_epoch": best_epoch,
                   "wall_time_sec": time.time() - t0}, f, indent=2, default=str)

    return {
        "ckpt": str(ckpt_path),
        "embeddings": str(emb_path),
        "log": str(log_path),
        "best_val": float(best_val),
        "best_epoch": int(best_epoch),
        "wall_time_sec": float(time.time() - t0),
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--modality", required=True, choices=["rna", "atac"])
    p.add_argument("--processed", required=True, help="path to processed h5ad")
    p.add_argument("--out", required=True, help="output dir")
    p.add_argument("--epochs", type=int, default=200)
    p.add_argument("--batch-size", type=int, default=512)
    p.add_argument("--beta", type=float, default=0.01)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--device", default="cuda")
    args = p.parse_args()

    cfg = TrainConfig(
        modality=args.modality,
        processed_path=args.processed,
        out_dir=args.out,
        n_epochs=args.epochs,
        batch_size=args.batch_size,
        beta=args.beta,
        seed=args.seed,
        device=args.device,
    )
    result = train(cfg)
    print(json.dumps(result, indent=2))
