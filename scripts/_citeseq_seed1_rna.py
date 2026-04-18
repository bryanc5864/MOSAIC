#!/usr/bin/env python3
"""Run RNA IB-VAE training for CITE-seq seed 1."""
import sys
sys.path.insert(0, 'C:/Users/Maozer/projects/MOSAIC')
from src.training.train_ibvae import TrainConfig, train
from src.utils.paths import PROCESSED_DIR
import json, os

os.makedirs('experiments/exp001_citeseq_seed1', exist_ok=True)
cfg = TrainConfig(
    modality='rna',
    processed_path=str(PROCESSED_DIR / 'citeseq_pbmc_rna.h5ad'),
    out_dir='experiments/exp001_citeseq_seed1',
    n_epochs=200, batch_size=512, beta=0.001,
    lambda_pred=5.0, patience=999, seed=1, device='cuda',
)
result = train(cfg)
print(f"RNA done, best_val={result['best_val']:.6f}", flush=True)
