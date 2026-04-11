# MIT License — Bryan Cheng, 2026
# Part of MOSAIC
"""Centralized path configuration for MOSAIC.

All file paths (raw data, processed data, experiments, figures) are resolved
relative to the project root so that scripts can be run from any cwd.
"""

from __future__ import annotations

from pathlib import Path


# Project root = parent of parent of this file (src/utils/paths.py -> src/utils -> src -> root)
ROOT = Path(__file__).resolve().parents[2]

# Data
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"

# Experiments, figures, logs
EXPERIMENTS_DIR = ROOT / "experiments"
FIGURES_DIR = ROOT / "figures"
LOGS_DIR = ROOT / "logs"

# Paper
PAPER_DIR = ROOT / "paper"


def ensure_dirs() -> None:
    """Create all standard data directories if they don't exist."""
    for d in (RAW_DIR, PROCESSED_DIR, EXPERIMENTS_DIR, FIGURES_DIR, LOGS_DIR):
        d.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    ensure_dirs()
    print(f"ROOT          : {ROOT}")
    print(f"RAW_DIR       : {RAW_DIR}")
    print(f"PROCESSED_DIR : {PROCESSED_DIR}")
    print(f"EXPERIMENTS   : {EXPERIMENTS_DIR}")
