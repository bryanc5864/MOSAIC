# MIT License
# Part of MOSAIC
"""Information-Bottleneck VAE for single-cell multi-omics.

Implements the IB-VAE described in RESEARCH_PLAN.md section 3.1.

Per-modality architecture:
  input x (n_vars)
    -> Encoder MLP: n_vars -> 512 -> 256 -> (mu, log_var)  (latent dim 64)
    -> reparameterize z = mu + eps * sigma
    -> Decoder MLP: 64 -> 256 -> 512 -> n_vars (modality-specific heads)
    -> Cross-modal prediction head: 64 -> cross_dim

Loss (per modality):
  L = L_recon + lambda_pred * ||y_cross_pred - y_cross||^2
      + beta * KL(q(z|x) || N(0, I))

where L_recon is:
  - ZINB negative log-likelihood for RNA
  - Binary cross-entropy on binarized peaks for ATAC

The encoder is trained on the preprocessed scaled/log-normalized features
(clean numeric input for MLP), while the decoder's recon loss operates on
the RAW COUNTS layer (for RNA) or the BINARY layer (for ATAC).
"""

from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
import torch.nn.functional as F

from src.models.zinb import log_zinb


# ----------------------------------------------------------------------------
# Building blocks
# ----------------------------------------------------------------------------


class MLP(nn.Module):
    """Simple MLP: Linear -> LayerNorm -> GELU -> Dropout, stacked."""

    def __init__(self, dims: list[int], dropout: float = 0.0):
        super().__init__()
        layers: list[nn.Module] = []
        for i in range(len(dims) - 1):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            # No norm/activation after the final linear layer by default.
            if i < len(dims) - 2:
                layers.append(nn.LayerNorm(dims[i + 1]))
                layers.append(nn.GELU())
                if dropout > 0:
                    layers.append(nn.Dropout(dropout))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class IBEncoder(nn.Module):
    """MLP encoder producing (mu, log_var) for a diagonal Gaussian posterior."""

    def __init__(self, n_vars: int, hidden: tuple[int, ...] = (512, 256),
                 latent_dim: int = 64, dropout: float = 0.1):
        super().__init__()
        self.trunk = MLP([n_vars, *hidden], dropout=dropout)
        self.mu_head = nn.Linear(hidden[-1], latent_dim)
        self.logvar_head = nn.Linear(hidden[-1], latent_dim)
        # Initialize logvar head small so initial KL is modest.
        nn.init.zeros_(self.logvar_head.weight)
        nn.init.constant_(self.logvar_head.bias, -2.0)  # sigma^2 ~ 0.135 initially

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        # Add LayerNorm + GELU on the trunk output for the heads
        h = F.gelu(F.layer_norm(self.trunk(x), (self.trunk.net[-1].out_features,)))
        mu = self.mu_head(h)
        logvar = self.logvar_head(h).clamp(min=-10.0, max=5.0)  # numerical safety
        return mu, logvar


def reparameterize(mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    std = torch.exp(0.5 * logvar)
    eps = torch.randn_like(std)
    return mu + eps * std


def kl_standard_normal(mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    """Mean KL( N(mu, sigma^2) || N(0, 1) ) across the batch, summed over dims."""
    # KL per element = 0.5 * (mu^2 + sigma^2 - 1 - log sigma^2)
    kl = 0.5 * (mu.pow(2) + logvar.exp() - 1.0 - logvar)
    return kl.sum(dim=1).mean()


# ----------------------------------------------------------------------------
# Modality-specific decoders
# ----------------------------------------------------------------------------


class ZINBDecoder(nn.Module):
    """ZINB decoder for RNA. Produces (rate, theta, pi_logits) given z.

    rate and pi_logits come from the decoder MLP; theta is a learnable
    per-gene parameter shared across cells (standard scVI practice).
    """

    def __init__(self, latent_dim: int, n_vars: int, hidden: tuple[int, ...] = (256, 512),
                 dropout: float = 0.1):
        super().__init__()
        self.n_vars = n_vars
        self.trunk = MLP([latent_dim, *hidden], dropout=dropout)
        h_out = hidden[-1]
        self.rate_head = nn.Linear(h_out, n_vars)        # -> softmax to proportions
        self.pi_head = nn.Linear(h_out, n_vars)          # zero-inflation logits
        # Per-gene dispersion, parameterized in log space for positivity.
        self.log_theta = nn.Parameter(torch.zeros(n_vars))

    def forward(self, z: torch.Tensor, lib_size: torch.Tensor
               ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        h = F.gelu(F.layer_norm(self.trunk(z), (self.trunk.net[-1].out_features,)))
        # Rate as fraction of library, then multiplied by observed library size
        # so mu stays in the scale of raw counts.
        rate_prop = F.softmax(self.rate_head(h), dim=-1)      # (B, G), sums to 1 per cell
        mu = rate_prop * lib_size.unsqueeze(1) + 1e-4         # small epsilon for stability
        theta = torch.exp(self.log_theta).unsqueeze(0).expand_as(mu) + 1e-4
        pi_logits = self.pi_head(h)
        return mu, theta, pi_logits


class BernoulliDecoder(nn.Module):
    """Bernoulli (BCE) decoder for ATAC. Produces peak-open logits given z."""

    def __init__(self, latent_dim: int, n_vars: int, hidden: tuple[int, ...] = (256, 512),
                 dropout: float = 0.1):
        super().__init__()
        self.n_vars = n_vars
        self.trunk = MLP([latent_dim, *hidden], dropout=dropout)
        h_out = hidden[-1]
        self.logit_head = nn.Linear(h_out, n_vars)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        h = F.gelu(F.layer_norm(self.trunk(z), (self.trunk.net[-1].out_features,)))
        return self.logit_head(h)  # raw logits; loss uses BCEWithLogits


# ----------------------------------------------------------------------------
# Full IB-VAE per modality
# ----------------------------------------------------------------------------


@dataclass
class IBVAELossBreakdown:
    total: torch.Tensor
    recon: torch.Tensor
    kl: torch.Tensor
    pred: torch.Tensor


class IBVAE_RNA(nn.Module):
    """IB-VAE for an RNA modality with ZINB reconstruction and cross-modal head."""

    def __init__(self, n_vars: int, cross_dim: int, latent_dim: int = 64,
                 hidden_enc: tuple[int, ...] = (512, 256),
                 hidden_dec: tuple[int, ...] = (256, 512),
                 dropout: float = 0.1):
        super().__init__()
        self.latent_dim = latent_dim
        self.encoder = IBEncoder(n_vars, hidden_enc, latent_dim, dropout)
        self.decoder = ZINBDecoder(latent_dim, n_vars, hidden_dec, dropout)
        self.cross_head = nn.Linear(latent_dim, cross_dim)

    @torch.no_grad()
    def embed(self, x: torch.Tensor) -> torch.Tensor:
        mu, _ = self.encoder(x)
        return mu

    def forward(self, x: torch.Tensor, raw_counts: torch.Tensor,
                y_cross: torch.Tensor, beta: float, lambda_pred: float
               ) -> IBVAELossBreakdown:
        mu, logvar = self.encoder(x)
        z = reparameterize(mu, logvar)
        lib_size = raw_counts.sum(dim=1).clamp(min=1.0)
        mu_rate, theta, pi_logits = self.decoder(z, lib_size)
        # Reconstruction loss (negative mean log-likelihood, per element).
        log_px = log_zinb(raw_counts, mu_rate, theta, pi_logits)
        recon = -log_px.mean()
        kl = kl_standard_normal(mu, logvar)
        y_pred = self.cross_head(z)
        pred = F.mse_loss(y_pred, y_cross)
        total = recon + beta * kl + lambda_pred * pred
        return IBVAELossBreakdown(total=total, recon=recon.detach(),
                                  kl=kl.detach(), pred=pred.detach())


class IBVAE_ATAC(nn.Module):
    """IB-VAE for an ATAC modality with Bernoulli (BCE) reconstruction."""

    def __init__(self, n_vars: int, cross_dim: int, latent_dim: int = 64,
                 hidden_enc: tuple[int, ...] = (512, 256),
                 hidden_dec: tuple[int, ...] = (256, 512),
                 dropout: float = 0.1):
        super().__init__()
        self.latent_dim = latent_dim
        self.encoder = IBEncoder(n_vars, hidden_enc, latent_dim, dropout)
        self.decoder = BernoulliDecoder(latent_dim, n_vars, hidden_dec, dropout)
        self.cross_head = nn.Linear(latent_dim, cross_dim)

    @torch.no_grad()
    def embed(self, x: torch.Tensor) -> torch.Tensor:
        mu, _ = self.encoder(x)
        return mu

    def forward(self, x: torch.Tensor, binary: torch.Tensor,
                y_cross: torch.Tensor, beta: float, lambda_pred: float
               ) -> IBVAELossBreakdown:
        mu, logvar = self.encoder(x)
        z = reparameterize(mu, logvar)
        logits = self.decoder(z)
        recon = F.binary_cross_entropy_with_logits(logits, binary, reduction="mean")
        kl = kl_standard_normal(mu, logvar)
        y_pred = self.cross_head(z)
        pred = F.mse_loss(y_pred, y_cross)
        total = recon + beta * kl + lambda_pred * pred
        return IBVAELossBreakdown(total=total, recon=recon.detach(),
                                  kl=kl.detach(), pred=pred.detach())


# ----------------------------------------------------------------------------
# Small-scale synthetic test
# ----------------------------------------------------------------------------


def _synthetic_test():
    """Train for a few steps on synthetic data and assert loss decreases."""
    import numpy as np
    torch.manual_seed(0)
    np.random.seed(0)

    B, G = 256, 100
    # Synthetic RNA counts: each "cell type" has a different Poisson rate
    cell_types = np.random.randint(0, 5, size=B)
    rate_per_type = np.random.rand(5, G).astype(np.float32) * 5 + 0.1
    counts = np.random.poisson(rate_per_type[cell_types]).astype(np.float32)
    x_norm = (np.log1p(counts) - np.log1p(counts).mean(0)) / (np.log1p(counts).std(0) + 1e-6)
    y_cross = np.random.rand(B, 10).astype(np.float32)

    x_norm_t = torch.from_numpy(x_norm)
    counts_t = torch.from_numpy(counts)
    y_cross_t = torch.from_numpy(y_cross)

    model = IBVAE_RNA(n_vars=G, cross_dim=10, latent_dim=8, hidden_enc=(64, 32), hidden_dec=(32, 64))
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)

    losses = []
    for step in range(60):
        opt.zero_grad()
        out = model(x_norm_t, counts_t, y_cross_t, beta=0.01, lambda_pred=1.0)
        out.total.backward()
        opt.step()
        losses.append(out.total.item())
    print(f"[rna] loss {losses[0]:.3f} -> {losses[-1]:.3f}")
    assert losses[-1] < losses[0] - 0.1, "IBVAE_RNA failed to decrease loss on synthetic data"

    # ATAC synthetic: Bernoulli peaks
    P = 80
    peak_prob = np.random.rand(5, P).astype(np.float32) * 0.5
    binary = (np.random.rand(B, P) < peak_prob[cell_types]).astype(np.float32)
    x_atac = (np.log1p(binary * 10) - np.log1p(binary * 10).mean(0)) / (np.log1p(binary * 10).std(0) + 1e-6)
    x_atac = np.nan_to_num(x_atac, nan=0.0, posinf=0.0, neginf=0.0)
    binary_t = torch.from_numpy(binary)
    x_atac_t = torch.from_numpy(x_atac.astype(np.float32))

    model_atac = IBVAE_ATAC(n_vars=P, cross_dim=10, latent_dim=8, hidden_enc=(64, 32), hidden_dec=(32, 64))
    opt_atac = torch.optim.Adam(model_atac.parameters(), lr=1e-3)

    losses = []
    for step in range(60):
        opt_atac.zero_grad()
        out = model_atac(x_atac_t, binary_t, y_cross_t, beta=0.01, lambda_pred=1.0)
        out.total.backward()
        opt_atac.step()
        losses.append(out.total.item())
    print(f"[atac] loss {losses[0]:.3f} -> {losses[-1]:.3f}")
    assert losses[-1] < losses[0] - 0.05, "IBVAE_ATAC failed to decrease loss on synthetic data"

    print("[ok] IB-VAE synthetic train test passed")


if __name__ == "__main__":
    _synthetic_test()
