# MIT License — Bryan Cheng, 2026
# Part of MOSAIC
"""Preprocessing pipelines for paired single-cell multi-omics datasets.

Pipeline (see RESEARCH_PLAN.md §5.2):
  1. Load raw h5 (10x Multiome format) or h5ad.
  2. Split into modality-specific AnnData objects.
  3. QC filter cells jointly (so pairing is preserved).
  4. RNA: normalize_total → log1p → top 2000 HVGs.
  5. ATAC: binarize → TF-IDF → LSI to 50 dims, drop first component.
  6. Compute cross-modal prediction targets:
       - y_atac_summary = top-50 LSI of ATAC (used as target when RNA is encoder input)
       - y_rna_summary  = top-50 PCA of RNA (used as target when ATAC is encoder input)
     Both are stored under .obsm so the training loop can read them without
     depending on the other modality's full matrix.
  7. Save processed AnnData to data/processed/<dataset_id>_{rna,atac}.h5ad.

This module is dataset-agnostic for 10x Multiome h5 format. Other dataset
formats (SHARE-seq tar, CITE-seq h5ad) get their own loaders.
"""

from __future__ import annotations

import gc
import hashlib
import warnings
from dataclasses import asdict
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD

from src.data.datasets import DATASETS, DatasetSpec
from src.utils.paths import PROCESSED_DIR

# Silence noisy scanpy warnings during routine preprocessing.
warnings.filterwarnings("ignore", category=UserWarning, module="anndata")
warnings.filterwarnings("ignore", category=FutureWarning, module="scanpy")


# ----------------------------------------------------------------------------
# Config — hyperparameters set by RESEARCH_PLAN.md §5.2.
# ----------------------------------------------------------------------------

RNA_MIN_GENES = 200              # cells must express at least this many genes
RNA_MIN_CELLS_PER_GENE = 3       # genes expressed in at least this many cells
RNA_N_HVG = 2000                 # number of highly variable genes to keep
RNA_NORM_TARGET = 1e4            # total-count normalization target
RNA_PCA_DIM = 50                 # RNA PCA dim used for clustering + cross-modal target

ATAC_MIN_PEAKS = 1000            # cells must have at least this many accessible peaks
ATAC_MIN_CELLS_PER_PEAK = 10     # peaks accessible in at least this many cells
ATAC_N_VAR_PEAKS = 10_000        # top variable peaks to keep (by prevalence variance)
# Note: reduced from 50K (plan default) to 10K to fit the single-GPU compute
# budget — 10K peaks still captures the top variable chromatin landscape while
# making the decoder's output layer tractable on an RTX 3080. Documented as a
# compute-driven deviation in RESEARCH_PLAN.md section 9.
ATAC_LSI_DIM = 50                # total LSI dims kept (first one dropped downstream)

LEIDEN_RES = 0.8                 # leiden resolution for cluster labeling

SEED = 0


# ----------------------------------------------------------------------------
# TF-IDF + LSI for ATAC  (replaces episcanpy, kept inline for stability)
# ----------------------------------------------------------------------------


def tfidf(X: sp.spmatrix) -> sp.spmatrix:
    """TF-IDF transform for sparse peak matrix (cells × peaks).

    Uses the standard Signac/ArchR formulation:
        TF_ij = x_ij / sum_j x_ij       (normalize per cell)
        IDF_j = log(1 + N / sum_i x_ij) (inverse peak frequency, smoothed)
        out_ij = TF_ij * IDF_j

    Binarizes input first (presence/absence of accessibility).
    """
    X = (X > 0).astype(np.float32)           # binarize
    X = sp.csr_matrix(X)
    n_cells = X.shape[0]
    # TF per cell
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0            # avoid div-by-zero
    tf_inv = 1.0 / row_sums
    # Scale rows: row i becomes x_i / sum_j x_ij
    X_tf = X.multiply(tf_inv[:, None]).tocsr()
    # IDF per peak
    peak_sums = np.asarray(X.sum(axis=0)).ravel()
    idf = np.log(1.0 + (n_cells / np.maximum(peak_sums, 1.0)))
    # Scale columns
    X_tfidf = X_tf.multiply(idf[None, :]).tocsr()
    return X_tfidf


def lsi(X_tfidf: sp.spmatrix, n_components: int, seed: int = SEED) -> np.ndarray:
    """Truncated SVD on TF-IDF matrix, producing LSI components.

    Returns a dense (n_cells, n_components) array. Caller should drop the
    first component — this is standard ATAC practice because component 1 is
    correlated with read depth.
    """
    svd = TruncatedSVD(n_components=n_components, random_state=seed, algorithm="arpack")
    Z = svd.fit_transform(X_tfidf)
    return Z.astype(np.float32)


# ----------------------------------------------------------------------------
# Main preprocessing entry points
# ----------------------------------------------------------------------------


def _md5_of_array(X) -> str:
    """MD5 of a sparse or dense array's bytes. For integrity logging."""
    if sp.issparse(X):
        data = np.ascontiguousarray(X.data).tobytes() + \
               np.ascontiguousarray(X.indices).tobytes() + \
               np.ascontiguousarray(X.indptr).tobytes()
    else:
        data = np.ascontiguousarray(X).tobytes()
    return hashlib.md5(data).hexdigest()


def load_10x_multiome_h5(path: Path) -> tuple[ad.AnnData, ad.AnnData]:
    """Load a 10x Multiome filtered_feature_bc_matrix.h5 and split into
    (rna, atac) AnnData objects with matching row order (same cells).
    """
    full = sc.read_10x_h5(path, gex_only=False)
    full.var_names_make_unique()
    ft = full.var["feature_types"].astype(str).values
    rna_mask = ft == "Gene Expression"
    atac_mask = ft == "Peaks"
    if rna_mask.sum() == 0 or atac_mask.sum() == 0:
        raise ValueError(
            f"Expected both 'Gene Expression' and 'Peaks' feature types; "
            f"got {np.unique(ft)}"
        )
    rna = full[:, rna_mask].copy()
    atac = full[:, atac_mask].copy()
    # Sanity: matching cell order
    assert np.array_equal(rna.obs_names.values, atac.obs_names.values), \
        "RNA and ATAC cell orders must match after split"
    # Record original pair index (each cell's row position) — preserved even
    # after later subsetting, used as the ground-truth pairing for evaluation.
    rna.obs["pair_idx"] = np.arange(rna.n_obs)
    atac.obs["pair_idx"] = np.arange(atac.n_obs)
    return rna, atac


def qc_filter_joint(rna: ad.AnnData, atac: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    """Filter cells that fail either modality's QC, then subset both to the
    intersection so pairing is preserved."""
    n0 = rna.n_obs
    # Per-cell counts
    rna_n_genes = np.asarray((rna.X > 0).sum(axis=1)).ravel()
    atac_n_peaks = np.asarray((atac.X > 0).sum(axis=1)).ravel()
    keep = (rna_n_genes >= RNA_MIN_GENES) & (atac_n_peaks >= ATAC_MIN_PEAKS)
    print(f"  [qc] cells: {n0} -> {keep.sum()}  "
          f"(rna_min_genes>={RNA_MIN_GENES}, atac_min_peaks>={ATAC_MIN_PEAKS})")
    rna_f = rna[keep].copy()
    atac_f = atac[keep].copy()
    # Gene/peak filters
    n_genes_before = rna_f.n_vars
    sc.pp.filter_genes(rna_f, min_cells=RNA_MIN_CELLS_PER_GENE)
    n_peaks_before = atac_f.n_vars
    # ATAC: filter peaks by prevalence (min cells)
    peak_prevalence = np.asarray((atac_f.X > 0).sum(axis=0)).ravel()
    atac_f = atac_f[:, peak_prevalence >= ATAC_MIN_CELLS_PER_PEAK].copy()
    print(f"  [qc] genes: {n_genes_before} -> {rna_f.n_vars}, "
          f"peaks: {n_peaks_before} -> {atac_f.n_vars}")
    return rna_f, atac_f


def preprocess_rna(rna: ad.AnnData) -> ad.AnnData:
    """Standard scRNA-seq pipeline: normalize → log1p → HVG → PCA."""
    # Keep raw counts on a layer for the ZINB decoder
    rna.layers["counts"] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=RNA_NORM_TARGET)
    sc.pp.log1p(rna)
    # HVG
    sc.pp.highly_variable_genes(rna, n_top_genes=RNA_N_HVG, flavor="seurat_v3",
                                 layer="counts", subset=True)
    # PCA for clustering and cross-modal target
    sc.pp.scale(rna, max_value=10)
    sc.tl.pca(rna, n_comps=RNA_PCA_DIM, random_state=SEED)
    # Leiden clustering — will provide cell-type labels
    sc.pp.neighbors(rna, n_neighbors=15, n_pcs=RNA_PCA_DIM, random_state=SEED)
    sc.tl.leiden(rna, resolution=LEIDEN_RES, random_state=SEED, flavor="igraph", n_iterations=2)
    # UMAP for later visualization
    sc.tl.umap(rna, random_state=SEED)
    return rna


def preprocess_atac(atac: ad.AnnData) -> ad.AnnData:
    """ATAC pipeline: select variable peaks → binarize → TF-IDF → LSI."""
    # Keep raw binarized counts on a layer (for BCE decoder)
    atac.layers["binary"] = (atac.X > 0).astype(np.float32).tocsr() \
        if sp.issparse(atac.X) else (atac.X > 0).astype(np.float32)
    # Top variable peaks by prevalence variance
    prev = np.asarray((atac.X > 0).sum(axis=0)).ravel().astype(np.float32)
    # Use prevalence-based variance proxy: prev * (1 - prev/n_cells) = Bernoulli var
    n_cells = atac.n_obs
    p_hat = prev / n_cells
    var_proxy = p_hat * (1 - p_hat)
    n_keep = min(ATAC_N_VAR_PEAKS, atac.n_vars)
    top_peaks = np.argsort(-var_proxy)[:n_keep]
    atac = atac[:, np.sort(top_peaks)].copy()
    print(f"  [atac] kept top {atac.n_vars} variable peaks")
    # TF-IDF + LSI
    X_tfidf = tfidf(atac.X)
    Z_lsi = lsi(X_tfidf, n_components=ATAC_LSI_DIM + 1, seed=SEED)
    # Drop the first component (correlated with depth); keep the next ATAC_LSI_DIM
    Z_lsi_final = Z_lsi[:, 1:1 + ATAC_LSI_DIM].astype(np.float32)
    atac.obsm["X_lsi"] = Z_lsi_final
    # For training, also store the log(1 + TF-IDF) matrix as an input feature layer.
    # Use .X so the model can read `adata.X` directly.
    # We keep the TF-IDF-transformed sparse matrix (log-transformed for stability).
    atac.layers["tfidf"] = X_tfidf.tocsr()
    atac.X = sp.csr_matrix(np.log1p(X_tfidf.toarray())).astype(np.float32)
    # UMAP on LSI for visualization (skip first LSI component)
    atac.obsm["X_pca"] = Z_lsi_final  # alias so scanpy's neighbors/umap work
    sc.pp.neighbors(atac, n_neighbors=15, use_rep="X_lsi", random_state=SEED)
    sc.tl.umap(atac, random_state=SEED)
    return atac


def add_cross_modal_targets(rna: ad.AnnData, atac: ad.AnnData) -> None:
    """Store cross-modal prediction targets in each modality's .obsm.

    DESIGN NOTE (revised after run-002 overfitting on per-cell targets):
    We use *cluster-centroid* targets rather than per-cell partner features.
    For each cell in modality A, y_cross is the mean of the partner modality's
    low-dim summary across all cells assigned to the same leiden cluster as A.

    Why cluster centroids:
      - All cells in cluster k share the same target -> the prediction head
        learns cluster identity (a generalizing signal) rather than per-cell
        memorization. With per-cell PCA/LSI targets, the decoder essentially
        memorized 10K unique fingerprints, train_pred -> 0 while val_pred
        increased -- a textbook overfitting failure.
      - Generalizes to new cells: at inference time, a held-out cell of
        cluster k will have the same target as training cells of cluster k,
        provided the cluster assignment is stable.
      - Still cross-modal: knowing the RNA cell's cluster lets you predict
        the average ATAC profile of that cluster.

    Cluster source: leiden on the RNA modality, propagated to ATAC via the
    paired ground truth (legitimate at training time on paired benchmark
    datasets; would need to be replaced by a clustering of one modality
    alone, or by an external label source, in the unpaired-application
    setting).

    Both centroid sets are standardized per-dimension to unit variance over
    the cells (NOT over the centroids), so the MSE loss is on the same scale
    as before.
    """
    assert rna.n_obs == atac.n_obs, "Must be called on paired, QC-matched AnnDatas"
    # Propagate leiden cluster labels (from RNA) to ATAC
    atac.obs["leiden"] = rna.obs["leiden"].values
    atac.obs["cell_type"] = rna.obs["leiden"].values
    rna.obs["cell_type"] = rna.obs["leiden"].values

    # Standardize per-cell modality summaries
    rna_pca = np.asarray(rna.obsm["X_pca"], dtype=np.float32)
    atac_lsi = np.asarray(atac.obsm["X_lsi"], dtype=np.float32)
    rna_pca_std = (rna_pca - rna_pca.mean(0)) / (rna_pca.std(0) + 1e-6)
    atac_lsi_std = (atac_lsi - atac_lsi.mean(0)) / (atac_lsi.std(0) + 1e-6)

    # Per-cluster centroids of the standardized representations
    clusters = rna.obs["leiden"].astype(str).values
    unique_clusters = np.unique(clusters)
    rna_centroid = {}
    atac_centroid = {}
    for c in unique_clusters:
        mask = clusters == c
        rna_centroid[c] = rna_pca_std[mask].mean(axis=0)
        atac_centroid[c] = atac_lsi_std[mask].mean(axis=0)

    # Per-cell y_cross targets: each cell receives its cluster's centroid
    # in the OTHER modality's space.
    rna_y_cross = np.stack([atac_centroid[c] for c in clusters]).astype(np.float32)
    atac_y_cross = np.stack([rna_centroid[c] for c in clusters]).astype(np.float32)
    rna.obsm["y_cross"] = rna_y_cross
    atac.obsm["y_cross"] = atac_y_cross
    # Also keep the per-cell standardized representations for downstream
    # analyses (e.g., entropy correlation diagnostics).
    rna.obsm["X_pca_std"] = rna_pca_std
    atac.obsm["X_lsi_std"] = atac_lsi_std


def save_processed(dataset_id: str, rna: ad.AnnData, atac: ad.AnnData) -> dict:
    """Save both AnnDatas and return a dict of integrity metadata."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    rna_path = PROCESSED_DIR / f"{dataset_id}_rna.h5ad"
    atac_path = PROCESSED_DIR / f"{dataset_id}_atac.h5ad"
    # anndata does not always serialize sparse tfidf layer cleanly; drop it.
    if "tfidf" in atac.layers:
        del atac.layers["tfidf"]
    rna.write_h5ad(rna_path, compression="gzip")
    atac.write_h5ad(atac_path, compression="gzip")
    meta = {
        "dataset_id": dataset_id,
        "rna_path": str(rna_path),
        "atac_path": str(atac_path),
        "rna_shape": list(rna.shape),
        "atac_shape": list(atac.shape),
        "rna_md5": _md5_of_array(rna.X),
        "atac_md5": _md5_of_array(atac.X),
        "n_clusters": int(rna.obs["leiden"].nunique()),
        "cluster_sizes": rna.obs["leiden"].value_counts().to_dict(),
    }
    return meta


def load_10x_citeseq_h5(path: Path) -> tuple[ad.AnnData, ad.AnnData]:
    """Load a 10x CITE-seq filtered_feature_bc_matrix.h5 and split into
    (rna, protein) AnnData objects. Feature types are 'Gene Expression'
    and 'Antibody Capture'.
    """
    full = sc.read_10x_h5(path, gex_only=False)
    full.var_names_make_unique()
    ft = full.var["feature_types"].astype(str).values
    rna_mask = ft == "Gene Expression"
    prot_mask = ft == "Antibody Capture"
    if rna_mask.sum() == 0 or prot_mask.sum() == 0:
        raise ValueError(
            f"Expected 'Gene Expression' and 'Antibody Capture'; got {np.unique(ft)}"
        )
    rna = full[:, rna_mask].copy()
    prot = full[:, prot_mask].copy()
    assert np.array_equal(rna.obs_names.values, prot.obs_names.values), \
        "RNA and protein cell orders must match"
    rna.obs["pair_idx"] = np.arange(rna.n_obs)
    prot.obs["pair_idx"] = np.arange(prot.n_obs)
    return rna, prot


def preprocess_protein(prot: ad.AnnData, drop_isotype: bool = True) -> ad.AnnData:
    """Preprocess CITE-seq ADT counts.

    Protein counts are dense (every cell expresses all measured proteins at
    some level), so we use CLR (centered log-ratio) normalization rather
    than TF-IDF. Isotype controls (negative-control antibodies used to
    measure background binding) are optionally dropped since they don't
    contribute cell-type signal.

    Output:
      - prot.layers['counts']: raw integer counts (for reference)
      - prot.X: CLR-normalized log-ratios, suitable as encoder input
      - prot.obsm['X_pca']: 16-dim PCA for cross-modal target (since only
        ~14 real proteins, we use a smaller PCA dim)
    """
    # Drop isotype controls (standard CITE-seq practice)
    if drop_isotype:
        is_control = np.array(["control" in str(n).lower() or "IgG" in str(n)
                                for n in prot.var_names])
        print(f"  [prot] dropping {is_control.sum()} isotype controls; "
              f"keeping {(~is_control).sum()} protein markers")
        prot = prot[:, ~is_control].copy()

    prot.layers["counts"] = prot.X.copy()
    # CLR normalization: log(x / gm(x)) per cell, where gm = geometric mean
    X = np.asarray(prot.X.toarray() if sp.issparse(prot.X) else prot.X,
                   dtype=np.float64)
    # Add pseudocount of 1
    X_p = X + 1.0
    log_X = np.log(X_p)
    geo_mean = log_X.mean(axis=1, keepdims=True)
    clr = log_X - geo_mean
    prot.X = clr.astype(np.float32)
    # PCA for cross-modal target and joint clustering
    n_comps = min(16, prot.n_vars - 1)
    sc.tl.pca(prot, n_comps=n_comps, random_state=SEED)
    # Add an X_lsi alias so the cross-modal-target code paths (which read
    # atac.obsm['X_lsi']) work without special-casing protein.
    prot.obsm["X_lsi"] = prot.obsm["X_pca"]
    # Neighbors + UMAP for sanity check
    sc.pp.neighbors(prot, n_neighbors=15, use_rep="X_pca", random_state=SEED)
    sc.tl.umap(prot, random_state=SEED)
    return prot


def run_pipeline_citeseq(dataset_id: str) -> dict:
    """End-to-end preprocessing for a 10x CITE-seq dataset."""
    spec = DATASETS[dataset_id]
    print(f"\n=== Preprocessing {dataset_id} ===")
    print(f"  source: {spec.raw_path}")
    rna, prot = load_10x_citeseq_h5(spec.raw_path)
    print(f"  loaded: rna {rna.shape}, protein {prot.shape}")

    # Joint QC: drop cells with few RNA reads (no protein threshold — every
    # cell expresses some level of every measured protein)
    n0 = rna.n_obs
    rna_n_genes = np.asarray((rna.X > 0).sum(axis=1)).ravel()
    keep = rna_n_genes >= RNA_MIN_GENES
    print(f"  [qc] cells: {n0} -> {keep.sum()}")
    rna = rna[keep].copy()
    prot = prot[keep].copy()
    # Gene filter
    sc.pp.filter_genes(rna, min_cells=RNA_MIN_CELLS_PER_GENE)

    # RNA side: same as multiome
    print("  [rna] running scanpy pipeline")
    rna = preprocess_rna(rna)
    # Protein side: CLR
    print("  [prot] running CLR + PCA")
    prot = preprocess_protein(prot)

    # Cross-modal targets: reuse the cluster-centroid mechanism from
    # add_cross_modal_targets, which reads rna.obsm['X_pca'] and
    # atac.obsm['X_lsi']. We aliased prot.obsm['X_lsi'] = prot.obsm['X_pca']
    # so the same function works unchanged.
    add_cross_modal_targets(rna, prot)
    # Save as dataset_id + rna/atac — use the "atac" suffix for protein so
    # downstream training code can load with recon_layer='binary'. We need
    # to also set prot.layers['binary'] to a reasonable binary signal for
    # the BCE decoder. For CITE-seq we binarize each protein around its
    # per-cell-median expression as a simple "expressed high vs low" indicator.
    prot_dense = np.asarray(prot.X)
    median_per_cell = np.median(prot_dense, axis=1, keepdims=True)
    prot.layers["binary"] = (prot_dense > median_per_cell).astype(np.float32)

    meta = save_processed(dataset_id, rna, prot)
    print(f"  [save] {meta['rna_path']}")
    print(f"  [save] {meta['atac_path']}")
    print(f"  [save] {meta['n_clusters']} leiden clusters")
    gc.collect()
    return meta


def run_pipeline_10x_multiome(dataset_id: str) -> dict:
    """End-to-end preprocessing for a 10x Multiome dataset."""
    spec = DATASETS[dataset_id]
    print(f"\n=== Preprocessing {dataset_id} ===")
    print(f"  source: {spec.raw_path}")
    rna, atac = load_10x_multiome_h5(spec.raw_path)
    print(f"  loaded: rna {rna.shape}, atac {atac.shape}")
    rna, atac = qc_filter_joint(rna, atac)
    print("  [rna] running scanpy pipeline")
    rna = preprocess_rna(rna)
    print("  [atac] running TF-IDF + LSI")
    atac = preprocess_atac(atac)
    add_cross_modal_targets(rna, atac)
    meta = save_processed(dataset_id, rna, atac)
    print(f"  [save] {meta['rna_path']}")
    print(f"  [save] {meta['atac_path']}")
    print(f"  [save] {meta['n_clusters']} leiden clusters")
    gc.collect()
    return meta


if __name__ == "__main__":
    import argparse
    import json

    p = argparse.ArgumentParser(description="MOSAIC preprocessing")
    p.add_argument("dataset", help="Dataset id (e.g. pbmc10k_multiome)")
    args = p.parse_args()

    if args.dataset not in DATASETS:
        raise SystemExit(f"unknown dataset: {args.dataset}")

    # Dispatch by modality pair
    spec = DATASETS[args.dataset]
    if "protein" in spec.modalities:
        meta = run_pipeline_citeseq(args.dataset)
    else:
        meta = run_pipeline_10x_multiome(args.dataset)
    # Write metadata
    meta_path = PROCESSED_DIR / f"{args.dataset}_meta.json"
    with meta_path.open("w") as f:
        json.dump(meta, f, indent=2, default=str)
    print(f"\n[done] metadata saved to {meta_path}")
