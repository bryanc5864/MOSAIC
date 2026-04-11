# MIT License — Bryan Cheng, 2026
# Part of MOSAIC
"""Zero-Inflated Negative Binomial (ZINB) log-likelihood.

Implemented inline (no scvi-tools dependency) to keep the environment stable.
Formulation follows the standard single-cell parameterization used in scVI
(Lopez et al. 2018):

    rate  μ  ≥ 0       (Poisson-like mean)
    dispersion  θ  > 0  (negative binomial inverse-dispersion; larger θ = less overdispersion)
    zero-inflation  π ∈ [0, 1]

ZINB likelihood at count x:
    if x == 0: P(x) = π + (1-π) · NB(0; μ, θ)
    if x > 0 : P(x) = (1-π) · NB(x; μ, θ)

where NB(x; μ, θ) = Γ(x+θ)/(Γ(x+1)Γ(θ)) · (θ/(μ+θ))^θ · (μ/(μ+θ))^x.

We work entirely in log-space for stability, using `torch.lgamma` and
`torch.logsumexp` for the zero case.
"""

from __future__ import annotations

import torch
import torch.nn.functional as F


def log_nb(x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor,
           eps: float = 1e-8) -> torch.Tensor:
    """Log of the negative binomial probability mass (per element).

    Shapes: x, mu, theta all broadcastable to (B, G).
    """
    # log Γ(x+θ) - log Γ(x+1) - log Γ(θ) + θ log(θ/(μ+θ)) + x log(μ/(μ+θ))
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
    """Log-probability of a ZINB observation, per element.

    `pi_logits` are the pre-sigmoid logits for the zero-inflation probability.
    """
    # Case x == 0: log[π + (1-π) NB(0; μ, θ)]
    log_pi = F.logsigmoid(pi_logits)       # log π
    log_1mpi = F.logsigmoid(-pi_logits)    # log (1-π)

    # NB log-prob at 0: θ · [log θ - log(μ+θ)]
    log_nb_zero = theta * (torch.log(theta + eps) - torch.log(theta + mu + eps))
    # log[π + (1-π) NB(0)] = logsumexp(log π, log(1-π) + log NB(0))
    log_zero = torch.logaddexp(log_pi, log_1mpi + log_nb_zero)

    # Case x > 0:  log(1-π) + log NB(x)
    log_nonzero = log_1mpi + log_nb(x, mu, theta, eps)

    zero_mask = (x < eps).float()
    return zero_mask * log_zero + (1.0 - zero_mask) * log_nonzero


def zinb_loss(x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor,
              pi_logits: torch.Tensor) -> torch.Tensor:
    """Mean-negative-log-likelihood over elements and batch. Returns a scalar."""
    return -log_zinb(x, mu, theta, pi_logits).mean()


if __name__ == "__main__":
    # Sanity check: synthetic NB-sampled data should give a loss that decreases
    # as θ approaches the true dispersion.
    torch.manual_seed(0)
    B, G = 512, 20
    true_mu = torch.full((B, G), 10.0)
    true_theta = torch.full((G,), 5.0)
    # Sample NB: gamma-poisson
    lam = torch.distributions.Gamma(true_theta, true_theta / true_mu).sample()
    x = torch.poisson(lam)
    pi_logits = torch.full((B, G), -10.0)  # virtually no zero inflation
    # At the true params, loss should be roughly the NB entropy.
    l_true = zinb_loss(x, true_mu, true_theta.expand(B, G), pi_logits)
    # At wrong theta, loss should be higher.
    l_wrong = zinb_loss(x, true_mu, torch.full((B, G), 0.5), pi_logits)
    print(f"loss(true theta=5.0)   = {l_true.item():.4f}")
    print(f"loss(wrong theta=0.5)  = {l_wrong.item():.4f}")
    assert l_wrong > l_true, "wrong dispersion should give higher loss"
    print("[ok] ZINB log-likelihood sanity check passed")
