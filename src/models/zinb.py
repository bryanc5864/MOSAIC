# MIT License
# Part of MOSAIC
"""Zero-inflated negative binomial log-likelihood.

Inlined rather than pulling in scvi-tools, which keeps the env stable. Same
parameterization as scVI (Lopez et al. 2018): mean mu >= 0, inverse-dispersion
theta > 0, zero-inflation pi in [0, 1]. Everything stays in log space.
"""

from __future__ import annotations

import torch
import torch.nn.functional as F


def log_nb(x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor,
           eps: float = 1e-8) -> torch.Tensor:
    """Per-element log NB pmf. x, mu, theta broadcast to (B, G)."""
    log_theta_mu_eps = torch.log(theta + mu + eps)
    log_nb = (
        torch.lgamma(x + theta)
        - torch.lgamma(x + 1)
        - torch.lgamma(theta)
        + theta * (torch.log(theta + eps) - log_theta_mu_eps)
        + x * (torch.log(mu + eps) - log_theta_mu_eps)
    )
    return log_nb


def log_zinb(x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor,
             pi_logits: torch.Tensor, eps: float = 1e-8) -> torch.Tensor:
    """Per-element ZINB log-prob. pi_logits are pre-sigmoid."""
    log_pi = F.logsigmoid(pi_logits)
    log_1mpi = F.logsigmoid(-pi_logits)

    log_nb_zero = theta * (torch.log(theta + eps) - torch.log(theta + mu + eps))
    log_zero = torch.logaddexp(log_pi, log_1mpi + log_nb_zero)

    log_nonzero = log_1mpi + log_nb(x, mu, theta, eps)

    zero_mask = (x < eps).float()
    return zero_mask * log_zero + (1.0 - zero_mask) * log_nonzero


def zinb_loss(x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor,
              pi_logits: torch.Tensor) -> torch.Tensor:
    """scalar mean negative log-likelihood"""
    return -log_zinb(x, mu, theta, pi_logits).mean()


if __name__ == "__main__":
    # loss at the true dispersion should beat a wrong one
    torch.manual_seed(0)
    B, G = 512, 20
    true_mu = torch.full((B, G), 10.0)
    true_theta = torch.full((G,), 5.0)
    # gamma-poisson sampling
    lam = torch.distributions.Gamma(true_theta, true_theta / true_mu).sample()
    x = torch.poisson(lam)
    pi_logits = torch.full((B, G), -10.0)  # basically no zero inflation
    l_true = zinb_loss(x, true_mu, true_theta.expand(B, G), pi_logits)
    l_wrong = zinb_loss(x, true_mu, torch.full((B, G), 0.5), pi_logits)
    print(f"loss(true theta=5.0)   = {l_true.item():.4f}")
    print(f"loss(wrong theta=0.5)  = {l_wrong.item():.4f}")
    assert l_wrong > l_true, "wrong dispersion should give higher loss"
    print("[ok] zinb sanity check passed")
