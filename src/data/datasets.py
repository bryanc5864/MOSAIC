# MIT License
# Part of MOSAIC
"""Dataset registry: URLs, sizes, and download helpers for the four paired
multi-omics datasets used in MOSAIC benchmarking, plus the stretch Tabula
Sapiens / ENCODE pair.

Every dataset record is a dataclass so that experiment configs can reference
them by ID and the download/preprocess code can look up URLs centrally.

References to each dataset's public source are in RESEARCH_PLAN.md §5.1.
"""

from __future__ import annotations

import hashlib
import os
import urllib.request
from dataclasses import dataclass, field

# 10x Genomics' CDN rejects requests without a User-Agent. Use a browser-like
# header for all downloads. This is a download convenience, not a scraping trick —
# the datasets are publicly hosted.
_USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/124.0.0.0 Safari/537.36"
)
from pathlib import Path
from typing import Optional

from src.utils.paths import RAW_DIR


@dataclass
class DatasetSpec:
    """Metadata for a single dataset to be downloaded and preprocessed."""

    dataset_id: str                 # short id, used in filenames
    description: str                # human-readable one-liner
    url: str                        # direct-download URL
    filename: str                   # local filename under RAW_DIR/<dataset_id>/
    modalities: tuple[str, ...]     # e.g. ("rna", "atac")
    approx_cells: int               # rough cell count
    approx_size_mb: float           # rough on-disk size after download
    md5: Optional[str] = None       # integrity hash if known
    notes: str = ""

    @property
    def raw_dir(self) -> Path:
        return RAW_DIR / self.dataset_id

    @property
    def raw_path(self) -> Path:
        return self.raw_dir / self.filename


# --- Dataset registry --------------------------------------------------------
#
# Note on URLs: these are the canonical public sources as of 2026-04. The 10x
# multiome datasets use the "filtered_feature_bc_matrix.h5" file, which is a
# single HDF5 containing BOTH modalities (Gene Expression + Peaks). scanpy's
# sc.read_10x_h5(..., gex_only=False) splits them by feature_type.
# -----------------------------------------------------------------------------

DATASETS: dict[str, DatasetSpec] = {
    # Primary debugging target: small, well-annotated, both modalities in one file.
    "pbmc10k_multiome": DatasetSpec(
        dataset_id="pbmc10k_multiome",
        description="10x Multiome PBMC granulocyte-sorted 10k (scRNA + scATAC, paired)",
        url="https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
        filename="pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
        modalities=("rna", "atac"),
        approx_cells=11_898,
        approx_size_mb=150.0,
        notes="10x public demo dataset. Contains Gene Expression + Peaks feature types.",
    ),
    # Second paired dataset: brain tissue, similar size.
    "brain3k_multiome": DatasetSpec(
        dataset_id="brain3k_multiome",
        description="10x Multiome E18 mouse brain 5k (scRNA + scATAC, paired)",
        url="https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5",
        filename="e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5",
        modalities=("rna", "atac"),
        approx_cells=4_900,
        approx_size_mb=130.0,
        notes="10x public demo; used as second paired benchmark.",
    ),
    # Third paired dataset: SHARE-seq skin. The original paper's processed
    # count matrices are on GEO (GSE140203); the 'SHARE-seq hair follicle'
    # mouse skin dataset is commonly used for multi-omics benchmarking.
    # We will use a preprocessed version from the scGLUE benchmark if the GEO
    # download is too large for the 10 GB budget; see download_all().
    "shareseq_skin": DatasetSpec(
        dataset_id="shareseq_skin",
        description="SHARE-seq mouse skin (scRNA + scATAC, paired) — Ma et al. 2020",
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140203/suppl/GSE140203_RAW.tar",
        filename="GSE140203_RAW.tar",
        modalities=("rna", "atac"),
        approx_cells=34_774,
        approx_size_mb=3_500.0,
        notes="Large. Subsample to ~15K cells during preprocessing.",
    ),
    # Fourth paired dataset: CITE-seq PBMC (RNA + surface protein).
    # 10x Genomics "10k PBMC TotalSeq-B" is a public CITE-seq demo dataset
    # with RNA + Antibody Capture (ADT) in a single h5.
    "citeseq_pbmc": DatasetSpec(
        dataset_id="citeseq_pbmc",
        description="10x CITE-seq PBMC 10k Protein v3 (scRNA + ADT, paired)",
        url="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5",
        filename="pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5",
        modalities=("rna", "protein"),
        approx_cells=7_865,
        approx_size_mb=70.0,
        notes="10x v3 public demo. 'Gene Expression' + 'Antibody Capture' feature types in one h5.",
    ),
    # Tabula Sapiens lung subset (for cross-atlas application Exp 6).
    "tabula_sapiens_lung": DatasetSpec(
        dataset_id="tabula_sapiens_lung",
        description="Tabula Sapiens lung subset (scRNA, unpaired)",
        url="https://datasets.cellxgene.cziscience.com/cbb0ec38-3c0d-4cdf-b76b-ba7c3b29a0fc.h5ad",
        filename="tabula_sapiens_lung.h5ad",
        modalities=("rna",),
        approx_cells=30_000,
        approx_size_mb=400.0,
        notes="Lung tissue subset from Tabula Sapiens (cellxgene). For cross-atlas Exp 6.",
    ),
    # 10x multiome lung (as an ATAC source for the cross-atlas experiment).
    # We use the 10x 2k lung multiome demo dataset.
    "lung_multiome": DatasetSpec(
        dataset_id="lung_multiome",
        description="10x Multiome Human Lung 2k (scRNA + scATAC, paired — use ATAC side as ENCODE proxy)",
        url="https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.h5",
        filename="human_brain_3k_filtered_feature_bc_matrix.h5",
        modalities=("rna", "atac"),
        approx_cells=3_000,
        approx_size_mb=100.0,
        notes="10x human brain 3k as ATAC partner for cross-atlas experiment with Tabula Sapiens.",
    ),
}


# ----------------------------------------------------------------------------
# Download helpers
# ----------------------------------------------------------------------------


def _md5_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            buf = f.read(chunk)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()


def download_one(spec: DatasetSpec, force: bool = False) -> Path:
    """Download a single dataset to its raw_path. Idempotent."""
    spec.raw_dir.mkdir(parents=True, exist_ok=True)
    path = spec.raw_path
    if path.exists() and not force:
        size_mb = path.stat().st_size / 1e6
        print(f"[skip] {spec.dataset_id}: already present ({size_mb:.1f} MB) at {path}")
        return path
    print(f"[download] {spec.dataset_id}: {spec.url}")
    print(f"           -> {path}")
    tmp = path.with_suffix(path.suffix + ".part")
    try:
        req = urllib.request.Request(spec.url, headers={"User-Agent": _USER_AGENT})
        with urllib.request.urlopen(req) as resp, tmp.open("wb") as out:
            total = int(resp.headers.get("Content-Length", 0))
            read = 0
            chunk = 1 << 20  # 1 MB
            while True:
                buf = resp.read(chunk)
                if not buf:
                    break
                out.write(buf)
                read += len(buf)
                if total:
                    pct = 100.0 * read / total
                    print(f"\r           {read/1e6:.1f} / {total/1e6:.1f} MB  ({pct:5.1f}%)",
                          end="", flush=True)
            print()
        tmp.replace(path)
    except Exception as e:
        if tmp.exists():
            tmp.unlink()
        raise RuntimeError(f"Download failed for {spec.dataset_id}: {e}") from e
    size_mb = path.stat().st_size / 1e6
    print(f"[ok]  {spec.dataset_id}: {size_mb:.1f} MB")
    if spec.md5 is not None:
        got = _md5_of(path)
        if got != spec.md5:
            raise RuntimeError(f"MD5 mismatch for {spec.dataset_id}: expected {spec.md5}, got {got}")
        print(f"[ok]  {spec.dataset_id}: md5 verified")
    return path


def download_all(which: Optional[list[str]] = None, force: bool = False) -> None:
    """Download one or more datasets. If `which` is None, downloads everything."""
    if which is None:
        which = list(DATASETS.keys())
    for name in which:
        if name not in DATASETS:
            raise KeyError(f"unknown dataset id: {name}")
        try:
            download_one(DATASETS[name], force=force)
        except Exception as e:
            print(f"[fail] {name}: {e}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Download MOSAIC datasets")
    p.add_argument("--which", nargs="*", default=None,
                   help="Subset of dataset IDs to download. Default: all.")
    p.add_argument("--force", action="store_true", help="Re-download even if file exists.")
    p.add_argument("--list", action="store_true", help="List available datasets and exit.")
    args = p.parse_args()

    if args.list:
        for k, v in DATASETS.items():
            print(f"{k:25s}  ~{v.approx_cells:>7} cells  ~{v.approx_size_mb:>6.0f} MB  {v.description}")
    else:
        download_all(which=args.which, force=args.force)
